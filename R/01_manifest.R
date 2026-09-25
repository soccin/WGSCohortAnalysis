# 01_manifest.R
#
# Manifests and the cohort table. A Tempo manifest is a CSV with columns
# TID, NID, ProjNo, PATH, Sig where PATH is the per-pair somatic.final.maf
# (SNV manifest) or final.bedpe (SV manifest) and Sig is the file md5. Two
# manifests are needed because SNV and SV output can come from different
# Tempo runs (the NoSV / DownSV split).

manifest_cols <- c("TID", "NID", "ProjNo", "PATH", "Sig")

#' Read a Tempo manifest CSV
#'
#' @param file Manifest path.
#' @return Tibble with the five manifest columns, `ProjNo` as character.
read_manifest <- function(file) {
  read_csv(file, show_col_types = FALSE, col_types = cols(.default = "c")) |>
    require_cols(manifest_cols, what = basename(file)) |>
    select(all_of(manifest_cols))
}

#' Apply the cohort rules to a manifest
#'
#' Rules are applied in order: drop tumors matching `exclude_tid_regex`,
#' keep only `include_projects` (when given), then drop duplicate rows by
#' `Sig` so a pair reachable by two paths is read once.
#'
#' @param manifest Tibble from `read_manifest()`.
#' @param exclude_tid_regex Regex on `TID`; `NULL` or `""` disables.
#' @param include_projects Character vector of `ProjNo`; `NULL` keeps all.
#' @return Filtered manifest.
apply_cohort_rules <- function(manifest, exclude_tid_regex = NULL,
                               include_projects = NULL) {
  out <- manifest
  if (!is.null(exclude_tid_regex) && nzchar(exclude_tid_regex)) {
    out <- filter(out, !str_detect(TID, exclude_tid_regex))
  }
  if (!is.null(include_projects) && length(include_projects) > 0) {
    out <- filter(out, ProjNo %in% as.character(include_projects))
  }
  distinct(out, Sig, .keep_all = TRUE)
}

#' Check a manifest for duplicate tumors and missing files
#'
#' @param manifest Tibble from `read_manifest()` after cohort rules.
#' @param what Label for messages.
#' @param check_files Verify every `PATH` exists.
#' @return The manifest, invisibly. Errors on duplicated `TID`.
check_manifest <- function(manifest, what = "manifest", check_files = TRUE) {
  dups <- manifest |> count(TID) |> filter(n > 1) |> pull(TID)
  if (length(dups) > 0) {
    wca_abort("{what}: duplicated TID after dedup by Sig: {str_c(dups, collapse = ', ')}")
  }
  if (check_files) {
    missing <- manifest$PATH[!fs::file_exists(manifest$PATH)]
    if (length(missing) > 0) {
      wca_warn("{what}: {length(missing)} PATH entries do not exist (first: {missing[[1]]})")
    }
  }
  invisible(manifest)
}

#' Tumors listed more than once in a manifest
#'
#' A tumor appears twice for one of two reasons: the same file reached by
#' two paths (`exact`, one md5, safe to drop) or two different files
#' (`conflict`, which needs a decision about which run is current).
#' `apply_cohort_rules()` collapses the first kind by `Sig` when the cohort
#' is built, so an exact duplicate is untidy rather than wrong; a conflict
#' is resolved nowhere and reaches the analysis as a duplicated tumor.
#'
#' @param manifest Tibble from `read_manifest()`.
#' @return Tibble with one row per manifest row whose `TID` appears more
#'   than once: `TID`, `dup` (`exact` or `conflict`), `n` (rows for that
#'   tumor) and the remaining manifest columns. Zero rows when the manifest
#'   lists every tumor once.
manifest_duplicates <- function(manifest) {
  require_cols(manifest, manifest_cols, "manifest")
  manifest |>
    select(all_of(manifest_cols)) |>
    add_count(TID, name = "n") |>
    filter(n > 1) |>
    group_by(TID) |>
    mutate(dup = if_else(n_distinct(Sig) == 1, "exact", "conflict")) |>
    ungroup() |>
    arrange(TID) |>
    select(TID, dup, n, NID, ProjNo, PATH, Sig)
}

#' Drop manifest rows that name a file already listed for the same tumor
#'
#' One row of each redundant group is kept and the rest are returned for the
#' record. Rows that share a `TID` but have different md5 are all kept:
#' choosing between two runs of the same tumor is not a decision this
#' function can make, so `manifest_duplicates()` reports them as conflicts
#' and a person picks one.
#'
#' A copy under `out/` wins by default. `out/<run>/somatic/` is where Tempo
#' writes, and a second path for the same md5 is a copy of it, so the run
#' directory is the location to record: it is the one a rerun refreshes and
#' the one every other row of the manifest already points at. `prefer`
#' outranks `keep`, and nothing here checks that the winning file still
#' exists; the caller warns about missing paths afterwards.
#'
#' @param manifest Tibble from `read_manifest()`.
#' @param keep Which row of a redundant group to keep once `prefer` has been
#'   applied: `first` in manifest order, or the `newest` or `oldest` file by
#'   modification time. A row whose file is missing loses to one whose file
#'   exists.
#' @param prefer Regex on `PATH`. A row that matches beats one that does
#'   not, whatever `keep` says. `""` or `NULL` disables the preference and
#'   leaves the choice to `keep` alone.
#' @return List: `manifest` (the kept rows in their original order) and
#'   `removed` (the dropped rows).
dedup_manifest <- function(manifest, keep = c("first", "newest", "oldest"),
                           prefer = "/out/") {
  keep <- match.arg(keep)
  require_cols(manifest, manifest_cols, "manifest")
  rows <- manifest |> select(all_of(manifest_cols)) |> mutate(.row = row_number())
  ranked <- rows |>
    mutate(.prefer = if (is.null(prefer) || !nzchar(prefer)) FALSE else str_detect(PATH, prefer))
  if (keep == "first") {
    ranked <- arrange(ranked, TID, desc(.prefer))
  } else {
    ranked <- mutate(ranked, .mtime = fs::file_info(PATH)$modification_time)
    ranked <- if (keep == "newest") {
      arrange(ranked, TID, desc(.prefer), desc(.mtime))
    } else {
      arrange(ranked, TID, desc(.prefer), .mtime)
    }
  }
  kept <- ranked |> distinct(TID, Sig, .keep_all = TRUE) |> pull(.row)
  list(
    manifest = rows |> filter(.row %in% kept) |> select(all_of(manifest_cols)),
    removed = rows |> filter(!.row %in% kept) |> select(all_of(manifest_cols))
  )
}

#' Read the projects workbook and return the project numbers to include
#'
#' @param file An xlsx with at least `Project` and `Ucode` columns.
#' @param include_ucode Keep projects whose `Ucode` is in this set; `NULL`
#'   keeps all.
#' @return Tibble of the selected project rows.
read_projects_file <- function(file, include_ucode = NULL) {
  proj <- readxl::read_excel(file) |>
    require_cols(c("Project", "Ucode"), what = basename(file)) |>
    mutate(Project = as.character(Project))
  if (!is.null(include_ucode) && length(include_ucode) > 0) {
    proj <- filter(proj, Ucode %in% include_ucode)
  }
  proj
}

#' Resolve the cohort's project list from the three cohort keys
#'
#' The cohort is the union of `include_projects` and the workbook rows whose
#' `Ucode` is in `include_ucode`; an empty union applies no project filter.
#' Three settings would silently select far more than intended, so each
#' stops with a message saying what to change (`docs/METHODS.md`, Cohort
#' selection):
#'
#' * `include_ucode` without a workbook, where it would be ignored;
#' * a workbook with `include_ucode` empty, which selects every workbook
#'   project whatever `include_projects` says;
#' * an `include_ucode` value that matches no workbook row, such as a typo.
#'
#' @param include_projects Character vector of `ProjNo`, or `NULL`.
#' @param include_ucode Character vector of `Ucode` values, or `NULL`.
#' @param workbook Every row of the projects workbook, from
#'   `read_projects_file(file)` with no `include_ucode`; `NULL` when
#'   `cohort.projects_file` is not set.
#' @return List: `projects`, the sorted `ProjNo` values to keep
#'   (`character(0)` means every project in the manifests), and `workbook`,
#'   the selected workbook rows or `NULL`.
select_cohort_projects <- function(include_projects = NULL, include_ucode = NULL,
                                   workbook = NULL) {
  include_projects <- as.character(include_projects)
  include_ucode <- as.character(include_ucode)
  see <- "The rule and a table of cases: docs/METHODS.md, Cohort selection."

  if (is.null(workbook) && length(include_ucode) > 0) {
    wca_stop_user(
      "cohort.include_ucode is set but cohort.projects_file is null",
      detail = c(
        str_glue("include_ucode: {str_c(include_ucode, collapse = ', ')}"),
        "Without a workbook the codes are ignored, and the cohort would be every project in the manifests.",
        see
      ),
      fix = c(
        "set cohort.projects_file to the projects workbook",
        "remove include_ucode and list the projects in cohort.include_projects"
      )
    )
  }

  if (!is.null(workbook) && length(include_ucode) == 0) {
    wca_stop_user(
      "cohort.projects_file is set but cohort.include_ucode is empty",
      detail = c(
        str_glue("An empty include_ucode selects every project in the workbook ({nrow(workbook)} rows)."),
        if (length(include_projects) > 0)
          "include_projects is added to that selection; it cannot narrow it.",
        see
      ),
      fix = c(
        "set cohort.include_ucode to the cohort's code",
        "set cohort.projects_file: null and list the projects in cohort.include_projects"
      )
    )
  }

  if (!is.null(workbook)) {
    known <- sort(unique(as.character(na.omit(workbook$Ucode))))
    unmatched <- setdiff(include_ucode, known)
    if (length(unmatched) > 0) {
      wca_stop_user(
        "cohort.include_ucode has codes that match no row of the projects workbook",
        detail = c(
          str_glue("not in the workbook: {str_c(unmatched, collapse = ', ')}"),
          "codes in the workbook:",
          str_c("    ", wca_truncate_list(known)),
          see
        ),
        fix = c(
          "correct the code in cohort.include_ucode",
          "add the cohort's projects to the workbook under that code"
        )
      )
    }
    workbook <- filter(workbook, as.character(Ucode) %in% include_ucode)
  }

  projects <- sort(unique(c(include_projects, workbook$Project)))
  list(projects = projects, workbook = workbook)
}

#' Build the cohort table from the SNV and SV manifests
#'
#' One row per tumor. `snv_dir` and `sv_dir` are the per-pair Tempo
#' directories (two levels above PATH) resolved independently, because the
#' two manifests can point at different runs. `has_facets` is true when the
#' FACETS gene-level file is found under `snv_dir`.
#'
#' @param snv_manifest,sv_manifest Tibbles from `read_manifest()` after
#'   `apply_cohort_rules()`. Either may be `NULL`.
#' @param require_both When `TRUE` the two tumor sets must be identical;
#'   otherwise the union is kept and the asymmetry is reported.
#' @return Cohort tibble: Sample, NID, ProjNo, snv_file, sv_file, snv_dir,
#'   sv_dir, has_snv, has_sv, has_facets, Sig_snv, Sig_sv.
build_cohort <- function(snv_manifest = NULL, sv_manifest = NULL, require_both = TRUE) {
  if (is.null(snv_manifest) && is.null(sv_manifest)) {
    wca_abort("build_cohort(): at least one manifest is required")
  }
  as_side <- function(m, side) {
    if (is.null(m)) {
      return(tibble(Sample = character(), NID = character(), ProjNo = character(),
        file = character(), Sig = character()))
    }
    check_manifest(m, what = str_glue("{side} manifest"))
    m |> transmute(Sample = TID, NID, ProjNo, file = PATH, Sig)
  }
  snv <- as_side(snv_manifest, "SNV") |>
    rename(snv_file = file, Sig_snv = Sig, NID_snv = NID, ProjNo_snv = ProjNo)
  sv <- as_side(sv_manifest, "SV") |>
    rename(sv_file = file, Sig_sv = Sig, NID_sv = NID, ProjNo_sv = ProjNo)

  only_snv <- setdiff(snv$Sample, sv$Sample)
  only_sv <- setdiff(sv$Sample, snv$Sample)
  if (!is.null(snv_manifest) && !is.null(sv_manifest) &&
      (length(only_snv) > 0 || length(only_sv) > 0)) {
    detail <- c(
      if (length(only_snv) > 0) c(
        str_glue("{length(only_snv)} tumor(s) in the SNV manifest with no SV row:"),
        str_c("    ", wca_truncate_list(only_snv))),
      if (length(only_sv) > 0) c(
        str_glue("{length(only_sv)} tumor(s) in the SV manifest with no SNV row:"),
        str_c("    ", wca_truncate_list(only_sv)))
    )
    if (require_both) {
      wca_stop_user(
        "SNV and SV manifests do not cover the same tumors",
        detail = detail,
        fix = c(
          "add the missing rows to the manifest named under `manifests:` in 00.PARAMS.yml",
          "drop those tumors with `cohort: exclude_tid_regex:` in 00.PARAMS.yml",
          "set `cohort: require_both: false` to keep them with one assay only"
        )
      )
    }
    wca_msg("  NOTE: SNV and SV manifests do not cover the same tumors")
    wca_msg_lines(detail, indent = "    ")
  }

  full_join(snv, sv, by = "Sample") |>
    mutate(
      NID = coalesce(NID_snv, NID_sv),
      ProjNo = coalesce(ProjNo_snv, ProjNo_sv),
      snv_dir = if_else(is.na(snv_file), NA_character_, as.character(fs::path_dir(fs::path_dir(snv_file)))),
      sv_dir = if_else(is.na(sv_file), NA_character_, as.character(fs::path_dir(fs::path_dir(sv_file)))),
      has_snv = !is.na(snv_file),
      has_sv = !is.na(sv_file),
      has_facets = map_lgl(snv_dir, \(d) !is.na(d) && !is.na(tempo_pair_files(d)$facets_gene))
    ) |>
    select(Sample, NID, ProjNo, snv_file, sv_file, snv_dir, sv_dir,
      has_snv, has_sv, has_facets, Sig_snv, Sig_sv) |>
    arrange(ProjNo, Sample)
}

#' Scan Tempo output roots for per-pair MAF and bedpe files
#'
#' Walks each root for `*.somatic.final.maf` and `*.final.bedpe` under a
#' `somatic/<PAIR>/` directory and returns manifest-shaped tables with the
#' file md5 as `Sig`. Use this to build manifests for a new project or the
#' test fixtures; the KEJ manifests are maintained separately.
#'
#' @param roots One or more directories to scan recursively.
#' @param proj_no Project number to record; defaults to the `out/<X>` name
#'   with any `Proj_` prefix removed, so `out/Proj_14887` and `out/14887`
#'   both give `14887`.
#' @return List with `snv` and `sv` manifest tibbles.
scan_tempo_outputs <- function(roots, proj_no = NULL) {
  files <- fs::dir_ls(roots, recurse = TRUE, type = "file",
    regexp = "/somatic/[^/]+/(combined_mutations/[^/]+\\.somatic\\.final\\.maf|combined_svs/[^/]+\\.final\\.bedpe)$")
  if (length(files) == 0) wca_abort("scan_tempo_outputs(): no Tempo files under {str_c(roots, collapse = ', ')}")
  tbl <- tibble(PATH = as.character(files)) |>
    mutate(
      pair = fs::path_file(fs::path_dir(fs::path_dir(PATH))),
      kind = if_else(str_ends(PATH, "\\.maf"), "snv", "sv"),
      ProjNo = proj_no %||% str_remove(fs::path_file(fs::path_dir(fs::path_dir(fs::path_dir(fs::path_dir(PATH))))), "^Proj_"),
      Sig = unname(tools::md5sum(PATH))
    ) |>
    bind_cols(split_pair_name(fs::path_file(fs::path_dir(fs::path_dir(files))))) |>
    transmute(kind, TID = Sample, NID, ProjNo, PATH, Sig)
  list(
    snv = tbl |> filter(kind == "snv") |> select(-kind),
    sv = tbl |> filter(kind == "sv") |> select(-kind)
  )
}

#' Write a pair of manifests as `tempoSNVManifest_<stamp>.csv` and
#' `tempoSVManifest_<stamp>.csv`
#'
#' @param manifests List with `snv` and `sv` manifest tibbles.
#' @param dir Output directory (created if needed).
#' @param stamp Date stamp for the file names.
#' @param overwrite Replace existing files; otherwise abort if either exists.
#' @return Named character vector of the two paths, invisibly.
write_manifests <- function(manifests, dir, stamp = date_stamp(), overwrite = FALSE) {
  fs::dir_create(dir)
  paths <- c(
    snv = fs::path(dir, str_glue("tempoSNVManifest_{stamp}.csv")),
    sv = fs::path(dir, str_glue("tempoSVManifest_{stamp}.csv"))
  )
  if (!overwrite && any(fs::file_exists(paths))) {
    wca_abort("write_manifests(): output exists (use overwrite or another stamp): ",
      "{str_c(paths[fs::file_exists(paths)], collapse = ', ')}")
  }
  write_csv(select(manifests$snv, all_of(manifest_cols)), paths[["snv"]])
  write_csv(select(manifests$sv, all_of(manifest_cols)), paths[["sv"]])
  invisible(paths)
}

#' Write the manifests returned by `scan_tempo_outputs()`
#'
#' @param scan List from `scan_tempo_outputs()`.
#' @param dir Output directory.
#' @param stamp Date stamp for the file names.
#' @return Named character vector of the two paths, invisibly.
manifest_from_scan <- function(scan, dir, stamp = date_stamp()) {
  write_manifests(scan, dir, stamp, overwrite = TRUE)
}

#' Merge newly scanned pairs into an existing manifest
#'
#' Keyed on `TID` (the tumor id is the sample identity throughout the
#' toolkit). A tumor found in `new` that is absent from `old` is added; one
#' that is present is updated in place when its normal, project, path or
#' file md5 differs and left untouched when nothing differs. Rows of `old`
#' whose tumor is not in `new` are kept as they are, duplicates included.
#' An updated tumor replaces every old row for that `TID`, so a pair that
#' used to be listed twice (two paths, one md5) collapses to one row.
#'
#' Two rows in `new` for the same `TID` with different md5 are an error,
#' because the scan cannot decide which run is current; put one root per
#' run in the roots file instead.
#'
#' Duplicate rows carried in from `old` are reported in `duplicates` and are
#' dropped only when `dedup` is set, because an old manifest is the record of
#' what earlier runs used and this function does not quietly rewrite it.
#'
#' @param old Tibble from `read_manifest()`.
#' @param new Manifest-shaped tibble, usually one of the tables from
#'   `scan_tempo_outputs()`.
#' @param dedup Drop rows that name a file already listed for the same tumor
#'   (see `dedup_manifest()`). Tumors listed with different files are never
#'   dropped whatever this is set to.
#' @param keep,prefer Which row of a redundant group `dedup` keeps; see
#'   `dedup_manifest()`.
#' @return List: `manifest` (old order, updated rows in place, added rows
#'   appended by ProjNo and TID), `changes` (one row per `TID` in either
#'   table with `status` added / updated / unchanged / kept, the old and new
#'   NID, ProjNo, PATH and Sig, and `detail` naming the changed fields, plus
#'   one row per dropped duplicate with `status` removed) and `duplicates`
#'   (`manifest_duplicates()` of the manifest as returned, so what is left
#'   for a person to resolve).
update_manifest <- function(old, new, dedup = FALSE, keep = c("first", "newest", "oldest"),
                            prefer = "/out/") {
  keep <- match.arg(keep)
  require_cols(old, manifest_cols, "old manifest")
  require_cols(new, manifest_cols, "new manifest")
  conflicts <- new |> distinct(TID, Sig) |> count(TID) |> filter(n > 1) |> pull(TID)
  if (length(conflicts) > 0) {
    wca_abort("update_manifest(): {length(conflicts)} TID found with different files in the scan: ",
      "{str_c(conflicts, collapse = ', ')}")
  }
  new <- distinct(new, TID, .keep_all = TRUE) |> select(all_of(manifest_cols))
  old <- select(old, all_of(manifest_cols))
  old_idx <- old |> mutate(.row = row_number())
  first_old <- old_idx |> group_by(TID) |> slice_head(n = 1) |> ungroup()

  cmp <- new |>
    left_join(first_old, by = "TID", suffix = c("_new", "_old")) |>
    mutate(
      chg_NID = !is.na(.row) & NID_new != NID_old,
      chg_ProjNo = !is.na(.row) & ProjNo_new != ProjNo_old,
      chg_PATH = !is.na(.row) & PATH_new != PATH_old,
      chg_content = !is.na(.row) & Sig_new != Sig_old,
      status = case_when(
        is.na(.row) ~ "added",
        chg_NID | chg_ProjNo | chg_PATH | chg_content ~ "updated",
        TRUE ~ "unchanged"
      ),
      detail = purrr::pmap_chr(list(chg_NID, chg_ProjNo, chg_PATH, chg_content),
        \(a, b, c, d) str_c(c("NID", "ProjNo", "PATH", "content")[c(a, b, c, d)], collapse = ";"))
    )
  kept <- old_idx |>
    filter(!TID %in% new$TID) |>
    distinct(TID, .keep_all = TRUE) |>
    transmute(TID, NID_old = NID, ProjNo_old = ProjNo, PATH_old = PATH, Sig_old = Sig,
      status = "kept", detail = "")
  changes <- bind_rows(cmp, kept) |>
    mutate(status = factor(status, levels = c("added", "updated", "unchanged", "kept"))) |>
    arrange(status, ProjNo_new, ProjNo_old, TID) |>
    mutate(status = as.character(status)) |>
    select(TID, status, detail, NID_old, NID_new, ProjNo_old, ProjNo_new,
      PATH_old, PATH_new, Sig_old, Sig_new)

  updated <- changes |> filter(status == "updated") |> pull(TID)
  added <- changes |> filter(status == "added") |> pull(TID)
  replaced <- old_idx |>
    filter(TID %in% updated) |>
    distinct(TID, .keep_all = TRUE) |>
    select(TID, .row) |>
    inner_join(new, by = "TID")
  manifest <- bind_rows(filter(old_idx, !TID %in% updated), replaced) |>
    arrange(.row) |>
    select(all_of(manifest_cols)) |>
    bind_rows(new |> filter(TID %in% added) |> arrange(ProjNo, TID))

  if (dedup) {
    dropped <- dedup_manifest(manifest, keep = keep, prefer = prefer)
    manifest <- dropped$manifest
    changes <- bind_rows(changes, dropped$removed |>
      transmute(TID, status = "removed", detail = "duplicate", NID_old = NID,
        ProjNo_old = ProjNo, PATH_old = PATH, Sig_old = Sig))
  }
  list(manifest = manifest, changes = changes,
    duplicates = manifest_duplicates(manifest))
}

#' Markdown report for one or more `update_manifest()` results
#'
#' @param updates Named list (e.g. `snv`, `sv`) of `update_manifest()` results.
#' @param inputs Named character vector of the input manifest paths.
#' @param roots Character vector of scanned roots.
#' @param outputs Named character vector of the written manifest paths.
#' @return Character vector of markdown lines.
manifest_update_report <- function(updates, inputs = character(), roots = character(),
                                   outputs = character()) {
  counts <- purrr::imap(updates, \(u, nm) {
    u$changes |> count(status) |> mutate(manifest = nm)
  }) |>
    bind_rows() |>
    tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0)
  for (col in c("added", "updated", "unchanged", "kept")) if (!col %in% names(counts)) counts[[col]] <- 0L
  counts <- counts |> mutate(total = added + updated + unchanged + kept) |>
    select(manifest, added, updated, unchanged, kept, total)
  md_table <- function(df) {
    df <- mutate(df, across(everything(), as.character))
    c(str_c("| ", str_c(names(df), collapse = " | "), " |"),
      str_c("|", str_c(rep("---", ncol(df)), collapse = "|"), "|"),
      purrr::pmap_chr(df, \(...) str_c("| ", str_c(c(...), collapse = " | "), " |")))
  }
  lines <- c(
    str_glue("# Manifest update {date_stamp()}"), "",
    "## Inputs", "",
    str_c("- ", names(inputs), ": `", inputs, "`"),
    str_c("- root: `", roots, "`"), "",
    "## Outputs", "",
    str_c("- ", names(outputs), ": `", outputs, "`"), "",
    "## Summary", "",
    "`added` = tumor not in the old manifest; `updated` = tumor present with a",
    "different normal, project, path or file md5; `unchanged` = scanned and",
    "identical; `kept` = old row not touched by this scan.", "",
    md_table(counts), ""
  )
  for (nm in names(updates)) {
    ch <- updates[[nm]]$changes
    lines <- c(lines, str_glue("## {nm}: added"), "")
    a <- ch |> filter(status == "added") |> select(TID, NID = NID_new, ProjNo = ProjNo_new, PATH = PATH_new)
    lines <- c(lines, if (nrow(a) == 0) "none" else md_table(a), "")
    lines <- c(lines, str_glue("## {nm}: updated"), "")
    u <- ch |> filter(status == "updated") |>
      select(TID, detail, NID_old, NID_new, ProjNo_old, ProjNo_new, PATH_old, PATH_new)
    lines <- c(lines, if (nrow(u) == 0) "none" else md_table(u), "")

    lines <- c(lines, str_glue("## {nm}: duplicates"), "")
    gone <- ch |> filter(status == "removed") |> select(TID, PATH = PATH_old, Sig = Sig_old)
    left <- updates[[nm]]$duplicates %||% tibble()
    if (nrow(gone) == 0 && nrow(left) == 0) {
      lines <- c(lines, "none: every tumor is listed once", "")
    } else {
      if (nrow(gone) > 0) {
        lines <- c(lines,
          str_glue("Dropped {nrow(gone)} row(s) naming a file already listed for the same tumor."),
          "", md_table(gone), "")
      }
      if (nrow(left) > 0) {
        lines <- c(lines,
          str_glue("{dplyr::n_distinct(left$TID)} tumor(s) still listed more than once. ",
            "`exact` is one file reached by two paths and can be dropped with --dedup; ",
            "`conflict` is different files for one tumor and needs a decision."),
          "", md_table(select(left, TID, dup, n, PATH, Sig)), "")
      }
    }
  }
  if (all(c("snv", "sv") %in% names(updates))) {
    touched <- function(nm) updates[[nm]]$changes |> filter(status %in% c("added", "updated")) |> pull(TID)
    only_snv <- setdiff(touched("snv"), touched("sv"))
    only_sv <- setdiff(touched("sv"), touched("snv"))
    lines <- c(lines, "## SNV / SV asymmetry", "",
      if (length(only_snv) + length(only_sv) == 0) "none: every added or updated tumor changed in both manifests" else c(
        if (length(only_snv) > 0) str_glue("- changed in SNV only: {str_c(only_snv, collapse = ', ')}"),
        if (length(only_sv) > 0) str_glue("- changed in SV only: {str_c(only_sv, collapse = ', ')}")), "")
  }
  lines
}

#' Changed rows of one or more `update_manifest()` results as a table
#'
#' The flat companion to `manifest_update_report()`: one row per tumor this
#' scan added or updated, in every manifest, with the old and new values
#' side by side for a review pass in a spreadsheet.
#'
#' @param updates Named list (e.g. `snv`, `sv`) of `update_manifest()` results.
#' @param statuses Statuses to keep, in the order they are reported. Add
#'   `unchanged` and `kept` to account for every tumor in the output
#'   manifest. `removed` rows are duplicates dropped by `dedup`, and are
#'   reported by default so that nothing leaves the manifest unseen.
#' @return Tibble: `manifest` followed by the `changes` columns of
#'   `update_manifest()`.
manifest_update_changes <- function(updates, statuses = c("added", "updated", "removed")) {
  purrr::imap(updates, \(u, nm) mutate(u$changes, manifest = nm, .before = 1)) |>
    bind_rows() |>
    filter(status %in% statuses) |>
    mutate(
      manifest = factor(manifest, levels = names(updates)),
      status = factor(status, levels = statuses)
    ) |>
    arrange(manifest, status, TID) |>
    mutate(across(c(manifest, status), as.character))
}
