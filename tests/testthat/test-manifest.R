test_that("scan_tempo_outputs finds the three fixture pairs", {
  sc <- scan_tempo_outputs(fixture_root())
  expect_equal(nrow(sc$snv), 3)
  expect_equal(nrow(sc$sv), 3)
  expect_setequal(sc$snv$TID, c("CL01", "CL02", "CL03"))
  expect_equal(sc$snv$NID, str_c(sc$snv$TID, "N"))
  expect_true(all(nchar(sc$snv$Sig) == 32))
  expect_true(all(fs::file_exists(sc$sv$PATH)))
})

test_that("manifest round-trips through manifest_from_scan and read_manifest", {
  sc <- scan_tempo_outputs(fixture_root())
  dir <- withr::local_tempdir()
  paths <- manifest_from_scan(sc, dir, stamp = "000000")
  m <- read_manifest(paths[["snv"]])
  expect_equal(names(m), c("TID", "NID", "ProjNo", "PATH", "Sig"))
  expect_equal(m$TID, sc$snv$TID)
})

test_that("cohort rules dedup by Sig, exclude by regex and project", {
  m <- tibble(
    TID = c("A", "A", "B_CL", "C"), NID = c("An", "An", "Bn", "Cn"),
    ProjNo = c("1", "1", "1", "2"), PATH = c("p1", "p1b", "p2", "p3"),
    Sig = c("s1", "s1", "s2", "s3")
  )
  expect_equal(apply_cohort_rules(m)$TID, c("A", "B_CL", "C"))
  expect_equal(apply_cohort_rules(m, exclude_tid_regex = "_CL")$TID, c("A", "C"))
  expect_equal(apply_cohort_rules(m, include_projects = "2")$TID, "C")
})

test_that("select_cohort_projects takes the union of the two ways", {
  wb <- tibble(Project = c("10_A", "10_B", "20_A", "30_A"),
               Ucode = c("X", "X", "Y", NA))
  expect_equal(select_cohort_projects()$projects, character(0))
  expect_null(select_cohort_projects()$workbook)
  expect_equal(select_cohort_projects(include_projects = c("20_A", "10_A"))$projects,
    c("10_A", "20_A"))
  by_code <- select_cohort_projects(include_ucode = "X", workbook = wb)
  expect_equal(by_code$projects, c("10_A", "10_B"))
  expect_equal(by_code$workbook$Project, c("10_A", "10_B"))
  expect_equal(select_cohort_projects("30_A", "X", wb)$projects, c("10_A", "10_B", "30_A"))
})

test_that("select_cohort_projects stops on settings that fail open", {
  wb <- tibble(Project = c("10_A", "10_B", "20_A"), Ucode = c("X", "X", "Y"))

  err <- expect_error(select_cohort_projects(include_ucode = "X"),
    class = "wca_user_error")
  expect_match(conditionMessage(err), "projects_file is null", fixed = TRUE)

  err <- expect_error(select_cohort_projects(include_projects = "10_A", workbook = wb),
    class = "wca_user_error")
  expect_match(conditionMessage(err), "include_ucode is empty", fixed = TRUE)
  expect_match(conditionMessage(err), "cannot narrow it", fixed = TRUE)
  expect_error(select_cohort_projects(workbook = wb), class = "wca_user_error")

  err <- expect_error(select_cohort_projects(include_ucode = c("X", "Z"), workbook = wb),
    class = "wca_user_error")
  expect_match(conditionMessage(err), "not in the workbook: Z", fixed = TRUE)
})

test_that("update_manifest adds, updates in place, keeps and reports", {
  old <- tibble(
    TID = c("A", "B", "C", "C", "D"), NID = c("An", "Bn", "Cn", "Cn", "Dn"),
    ProjNo = c("1", "1", "2", "2", "2"), PATH = c("pA", "pB", "pC1", "pC2", "pD"),
    Sig = c("sA", "sB", "sC", "sC", "sD")
  )
  new <- tibble(
    TID = c("A", "B", "C", "E"), NID = c("An", "Bn2", "Cn", "En"),
    ProjNo = c("1", "1", "2", "3"), PATH = c("pA", "pB", "pC3", "pE"),
    Sig = c("sA", "sB", "sC9", "sE")
  )
  u <- update_manifest(old, new)
  m <- u$manifest
  expect_equal(m$TID, c("A", "B", "C", "D", "E"))
  expect_equal(m$PATH, c("pA", "pB", "pC3", "pD", "pE"))
  expect_equal(m$NID[m$TID == "B"], "Bn2")
  expect_equal(names(m), c("TID", "NID", "ProjNo", "PATH", "Sig"))
  st <- set_names(u$changes$status, u$changes$TID)
  expect_equal(st[c("A", "B", "C", "D", "E")],
    c(A = "unchanged", B = "updated", C = "updated", D = "kept", E = "added"))
  det <- set_names(u$changes$detail, u$changes$TID)
  expect_equal(unname(det[["B"]]), "NID")
  expect_equal(unname(det[["C"]]), "PATH;content")
  expect_equal(nrow(u$changes), 5)

  # an unchanged tumor keeps its duplicate old rows
  same <- update_manifest(old, filter(new, TID == "A"))
  expect_equal(same$manifest, old)

  # conflicting files for one tumor in the scan are an error
  bad <- bind_rows(new, tibble(TID = "E", NID = "En", ProjNo = "3", PATH = "pE2", Sig = "sE2"))
  expect_error(update_manifest(old, bad), "different files")
  # the same file reached twice is fine
  twice <- bind_rows(new, tibble(TID = "E", NID = "En", ProjNo = "3", PATH = "pE2", Sig = "sE"))
  expect_equal(nrow(update_manifest(old, twice)$manifest), 5)

  rep <- manifest_update_report(list(snv = u, sv = u), inputs = c(snv = "a", sv = "b"),
    roots = "r", outputs = c(snv = "x", sv = "y"))
  expect_true(any(str_detect(rep, "^\\| snv \\| 1 \\| 2 \\| 1 \\| 1 \\| 5 \\|$")))
  expect_true(any(str_detect(rep, "none: every added")))
})

test_that("manifest_update_changes tabulates the added and updated tumors", {
  old <- tibble(
    TID = c("A", "B", "D"), NID = c("An", "Bn", "Dn"), ProjNo = c("1", "1", "2"),
    PATH = c("pA", "pB", "pD"), Sig = c("sA", "sB", "sD")
  )
  new <- tibble(
    TID = c("A", "B", "E"), NID = c("An", "Bn2", "En"), ProjNo = c("1", "1", "3"),
    PATH = c("pA", "pB", "pE"), Sig = c("sA", "sB", "sE")
  )
  u <- update_manifest(old, new)
  tab <- manifest_update_changes(list(snv = u, sv = u))

  # added and updated only, added first, manifests in the order given
  expect_equal(tab$manifest, rep(c("snv", "sv"), each = 2))
  expect_equal(tab$status, rep(c("added", "updated"), 2))
  expect_equal(tab$TID, rep(c("E", "B"), 2))
  expect_equal(names(tab), c("manifest", names(u$changes)))

  # an added tumor has no old side, an updated one carries both
  added <- tab |> filter(status == "added") |> head(1)
  expect_true(is.na(added$PATH_old))
  expect_equal(added$PATH_new, "pE")
  upd <- tab |> filter(status == "updated") |> head(1)
  expect_equal(c(upd$NID_old, upd$NID_new), c("Bn", "Bn2"))
  expect_equal(upd$detail, "NID")

  # unchanged and kept are reachable but off by default
  wide <- manifest_update_changes(list(snv = u), statuses = c("added", "updated", "unchanged", "kept"))
  expect_setequal(wide$TID, c("A", "B", "D", "E"))
  expect_equal(nrow(manifest_update_changes(list(snv = u), statuses = "kept")), 1)
})

test_that("manifest_duplicates separates the same file twice from two different files", {
  m <- tibble(
    TID = c("A", "A", "B", "B", "C"), NID = c("An", "An", "Bn", "Bn", "Cn"),
    ProjNo = c("1", "1", "1", "1", "2"),
    PATH = c("pA1", "pA2", "pB1", "pB2", "pC"),
    Sig = c("sA", "sA", "sB", "sB9", "sC")
  )
  d <- manifest_duplicates(m)
  expect_equal(names(d), c("TID", "dup", "n", "NID", "ProjNo", "PATH", "Sig"))
  expect_equal(d$TID, c("A", "A", "B", "B"))
  expect_equal(unique(d$dup[d$TID == "A"]), "exact")
  expect_equal(unique(d$dup[d$TID == "B"]), "conflict")
  expect_equal(unique(d$n), 2L)
  expect_equal(nrow(manifest_duplicates(filter(m, TID == "C"))), 0)
  expect_equal(nrow(manifest_duplicates(m[0, ])), 0)

  # only the redundant row goes; the two different files for B both stay
  dd <- dedup_manifest(m)
  expect_equal(dd$manifest$PATH, c("pA1", "pB1", "pB2", "pC"))
  expect_equal(dd$removed$PATH, "pA2")
  expect_equal(manifest_duplicates(dd$manifest)$TID, c("B", "B"))
})

test_that("dedup_manifest keeps the newest or oldest file when asked", {
  dir <- withr::local_tempdir()
  older <- fs::path(dir, "older.maf")
  newer <- fs::path(dir, "newer.maf")
  writeLines("x", older)
  writeLines("x", newer)
  Sys.setFileTime(older, Sys.time() - 3600)
  m <- tibble(TID = "A", NID = "An", ProjNo = "1",
    PATH = as.character(c(older, newer)), Sig = "sA")
  expect_equal(dedup_manifest(m)$manifest$PATH, as.character(older))
  expect_equal(dedup_manifest(m, keep = "newest")$manifest$PATH, as.character(newer))
  expect_equal(dedup_manifest(m, keep = "oldest")$manifest$PATH, as.character(older))

  # a row whose file is gone loses to one whose file is there
  gone <- mutate(m, PATH = as.character(c(fs::path(dir, "gone.maf"), older)))
  expect_equal(dedup_manifest(gone, keep = "newest")$manifest$PATH, as.character(older))
})

test_that("dedup_manifest keeps the copy under out/ whatever the tie-break says", {
  dir <- withr::local_tempdir()
  in_out <- fs::path(dir, "out", "s.maf")
  copy <- fs::path(dir, "results", "r_002", "s.maf")
  fs::dir_create(fs::path_dir(c(in_out, copy)))
  writeLines("x", in_out)
  writeLines("x", copy)
  Sys.setFileTime(in_out, Sys.time() - 3600)
  # the out/ copy is both the second row and the older file, so it loses
  # every tie-break and can only win on the path rule
  m <- tibble(TID = "A", NID = "An", ProjNo = "1",
    PATH = as.character(c(copy, in_out)), Sig = "sA")
  expect_equal(dedup_manifest(m)$manifest$PATH, as.character(in_out))
  expect_equal(dedup_manifest(m, keep = "newest")$manifest$PATH, as.character(in_out))
  expect_equal(dedup_manifest(m)$removed$PATH, as.character(copy))

  # the rule can be turned off or pointed somewhere else
  expect_equal(dedup_manifest(m, prefer = "")$manifest$PATH, as.character(copy))
  expect_equal(dedup_manifest(m, keep = "newest", prefer = "")$manifest$PATH, as.character(copy))
  expect_equal(dedup_manifest(m, prefer = "r_002")$manifest$PATH, as.character(copy))

  # neither side matching falls back to the tie-break
  expect_equal(dedup_manifest(mutate(m, PATH = c("pA1", "pA2")))$manifest$PATH, "pA1")
})

test_that("update_manifest reports duplicates and drops them only when asked", {
  old <- tibble(
    TID = c("A", "A", "B"), NID = c("An", "An", "Bn"), ProjNo = "1",
    PATH = c("pA1", "pA2", "pB"), Sig = c("sA", "sA", "sB")
  )
  new <- tibble(TID = "B", NID = "Bn", ProjNo = "1", PATH = "pB", Sig = "sB")

  u <- update_manifest(old, new)
  expect_equal(nrow(u$manifest), 3)
  expect_equal(u$duplicates$TID, c("A", "A"))
  expect_equal(unique(u$duplicates$dup), "exact")
  expect_false("removed" %in% u$changes$status)

  ded <- update_manifest(old, new, dedup = TRUE)
  expect_equal(ded$manifest$PATH, c("pA1", "pB"))
  expect_equal(nrow(ded$duplicates), 0)
  drop <- ded$changes |> filter(status == "removed")
  expect_equal(drop$PATH_old, "pA2")
  expect_equal(drop$detail, "duplicate")
  expect_true(is.na(drop$PATH_new))

  # the dropped row is reported by default, in the table and in the report
  expect_equal(manifest_update_changes(list(snv = ded))$status, "removed")
  expect_true(any(str_detect(manifest_update_report(list(snv = ded)), "Dropped 1 row")))
  expect_true(any(str_detect(manifest_update_report(list(snv = u)), "still listed more than once")))
  expect_true(any(str_detect(manifest_update_report(list(snv = update_manifest(new, new))),
    "none: every tumor is listed once")))
})

test_that("wcaUpdateManifests.R flags duplicates and drops them with --dedup", {
  home <- wca_home()
  sc <- scan_tempo_outputs(fixture_root())
  dir <- withr::local_tempdir()
  copy <- sc$snv |> filter(TID == "CL01") |>
    mutate(PATH = str_replace(PATH, "/out/", "/results/r_002/tempo/"))
  write_manifests(list(snv = bind_rows(sc$snv, copy), sv = sc$sv), dir, stamp = "old")
  roots_file <- fs::path(dir, "roots.txt")
  writeLines(as.character(fixture_root()), roots_file)
  run <- function(stamp, ...) {
    suppressWarnings(system2("Rscript", c(fs::path(home, "bin", "wcaUpdateManifests.R"),
      fs::path(dir, "tempoSNVManifest_old.csv"), fs::path(dir, "tempoSVManifest_old.csv"),
      roots_file, dir, str_c("--stamp=", stamp), ...), stdout = TRUE, stderr = TRUE))
  }

  flagged <- run("flag")
  expect_equal(attr(flagged, "status") %||% 0, 0, info = str_c(flagged, collapse = "\n"))
  expect_true(any(str_detect(flagged, "rerun with --dedup")))
  expect_equal(nrow(read_manifest(fs::path(dir, "tempoSNVManifest_flag.csv"))), 4)

  cut <- run("cut", "--dedup")
  expect_equal(attr(cut, "status") %||% 0, 0, info = str_c(cut, collapse = "\n"))
  expect_true(any(str_detect(cut, "dropped 1 duplicate row")))
  expect_false(any(str_detect(cut, "rerun with --dedup")))
  snv <- read_manifest(fs::path(dir, "tempoSNVManifest_cut.csv"))
  expect_equal(snv$TID, c("CL01", "CL02", "CL03"))
  chg <- readr::read_csv(fs::path(dir, "manifestUpdate_cut.csv"), show_col_types = FALSE)
  expect_equal(chg |> filter(status == "removed") |> pull(PATH_old), copy$PATH)
  expect_true(any(str_detect(readr::read_lines(fs::path(dir, "manifestUpdate_cut.md")),
    "Dropped 1 row")))
  # the copy under out/ is the row that stays
  expect_equal(snv$PATH[snv$TID == "CL01"], sc$snv$PATH[sc$snv$TID == "CL01"])

  # --prefer points the rule somewhere else
  other <- run("other", "--dedup", "--prefer=r_002")
  expect_equal(attr(other, "status") %||% 0, 0, info = str_c(other, collapse = "\n"))
  moved <- read_manifest(fs::path(dir, "tempoSNVManifest_other.csv"))
  expect_equal(moved$PATH[moved$TID == "CL01"], copy$PATH)
})

test_that("scan_tempo_outputs strips a Proj_ prefix from the run directory", {
  root <- withr::local_tempdir()
  pair <- fs::path(root, "Proj_99", "out", "Proj_99", "somatic", "T1__N1")
  fs::dir_create(fs::path(pair, "combined_mutations"))
  fs::dir_create(fs::path(pair, "combined_svs"))
  writeLines("x", fs::path(pair, "combined_mutations", "T1__N1.somatic.final.maf"))
  writeLines("y", fs::path(pair, "combined_svs", "T1__N1.final.bedpe"))
  sc <- scan_tempo_outputs(root)
  expect_equal(sc$snv$ProjNo, "99")
  expect_equal(sc$sv$ProjNo, "99")
  expect_equal(scan_tempo_outputs(root, proj_no = "Proj_99")$snv$ProjNo, "Proj_99")
})

test_that("wcaUpdateManifests.R adds, updates and writes a report", {
  home <- wca_home()
  sc <- scan_tempo_outputs(fixture_root())
  dir <- withr::local_tempdir()
  old_snv <- sc$snv |> filter(TID != "CL03") |> mutate(PATH = if_else(TID == "CL02", "/old/CL02.maf", PATH),
    Sig = if_else(TID == "CL02", "0", Sig))
  old_sv <- sc$sv |> filter(TID != "CL03")
  write_manifests(list(snv = old_snv, sv = old_sv), dir, stamp = "old")
  roots_file <- fs::path(dir, "roots.txt")
  writeLines(c("# fixture root", "", as.character(fixture_root())), roots_file)
  out <- system2("Rscript", c(fs::path(home, "bin", "wcaUpdateManifests.R"),
    fs::path(dir, "tempoSNVManifest_old.csv"), fs::path(dir, "tempoSVManifest_old.csv"),
    roots_file, dir, "--stamp=new"), stdout = TRUE, stderr = TRUE)
  expect_equal(attr(out, "status") %||% 0, 0, info = str_c(out, collapse = "\n"))
  snv <- read_manifest(fs::path(dir, "tempoSNVManifest_new.csv"))
  sv <- read_manifest(fs::path(dir, "tempoSVManifest_new.csv"))
  expect_equal(snv$TID, c("CL01", "CL02", "CL03"))
  expect_equal(snv$PATH[snv$TID == "CL02"], sc$snv$PATH[sc$snv$TID == "CL02"])
  expect_equal(sv$TID, c("CL01", "CL02", "CL03"))
  rep <- readr::read_lines(fs::path(dir, "manifestUpdate_new.md"))
  expect_true(any(str_detect(rep, "^\\| snv \\| 1 \\| 1 \\| 1 \\| 0 \\| 3 \\|$")))
  expect_true(any(str_detect(rep, "^\\| sv \\| 1 \\| 0 \\| 2 \\| 0 \\| 3 \\|$")))
  expect_true(any(str_detect(rep, "changed in SNV only: CL02")))
  chg <- readr::read_csv(fs::path(dir, "manifestUpdate_new.csv"), show_col_types = FALSE)
  expect_equal(names(chg)[1:3], c("manifest", "TID", "status"))
  snv_chg <- chg |> filter(manifest == "snv") |> arrange(TID)
  expect_equal(snv_chg$TID, c("CL02", "CL03"))
  expect_equal(snv_chg$status, c("updated", "added"))
  expect_equal(chg |> filter(manifest == "sv") |> pull(TID), "CL03")
  expect_equal(chg$PATH_new[chg$manifest == "snv" & chg$TID == "CL02"],
    sc$snv$PATH[sc$snv$TID == "CL02"])
  # refuses to overwrite without --force
  out2 <- suppressWarnings(system2("Rscript", c(fs::path(home, "bin", "wcaUpdateManifests.R"),
    fs::path(dir, "tempoSNVManifest_old.csv"), fs::path(dir, "tempoSVManifest_old.csv"),
    roots_file, dir, "--stamp=new"), stdout = TRUE, stderr = TRUE))
  expect_true((attr(out2, "status") %||% 0) != 0)
})

test_that("check_manifest rejects duplicated TIDs", {
  m <- tibble(TID = c("A", "A"), NID = "n", ProjNo = "1", PATH = c("x", "y"), Sig = c("1", "2"))
  expect_error(check_manifest(m, check_files = FALSE), "duplicated TID")
})

test_that("build_cohort resolves dirs and flags, and enforces symmetry", {
  coh <- fixture_cohort()
  expect_equal(nrow(coh), 3)
  expect_true(all(coh$has_snv & coh$has_sv & coh$has_facets))
  expect_equal(fs::path_file(coh$snv_dir), str_c(coh$Sample, "__", coh$NID))
  expect_equal(coh$snv_dir, coh$sv_dir)

  sc <- scan_tempo_outputs(fixture_root())
  dropped <- sc$sv$TID[[1]]
  err <- expect_error(build_cohort(sc$snv, sc$sv[-1, ], require_both = TRUE),
    class = "wca_user_error")
  expect_match(conditionMessage(err), "do not cover the same tumors")
  expect_match(conditionMessage(err), dropped, fixed = TRUE)
  expect_match(conditionMessage(err), "require_both", fixed = TRUE)
  asym <- build_cohort(sc$snv, sc$sv[-1, ], require_both = FALSE)
  expect_equal(sum(!asym$has_sv), 1)
  only_snv <- build_cohort(sc$snv, NULL)
  expect_true(all(!only_snv$has_sv))
})
