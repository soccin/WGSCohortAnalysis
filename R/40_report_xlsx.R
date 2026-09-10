# 40_report_xlsx.R
#
# xlsx report writer: one sheet per table, PCT columns formatted as percent,
# header row bold, auto column widths, and a DataDictionary sheet built from
# R/41_dictionaries.R plus the denominators attached to each table.

#' Write a list of tables to an xlsx workbook
#'
#' Any column whose name matches `pct_regex` gets the Excel number format
#' `0.00%`, so a stored 0.28 renders as 28.00%. A `DataDictionary` sheet is
#' appended describing every column of every sheet (unknown columns get an
#' empty description) and listing the denominator of each table that
#' carries one (see `recurrence_table()`).
#'
#' @param tables Named list of data frames; names become sheet names.
#' @param file Output path; overwritten.
#' @param dictionary Tibble with Sheet, Column, Description; see
#'   `wca_dictionary()`. `NULL` skips the sheet.
#' @param pct_regex Regex selecting percent columns.
#' @param notes Optional named character vector of free-text notes written
#'   at the top of the DataDictionary sheet (e.g. cohort size, cutoffs).
#' @return The openxlsx workbook, invisibly.
write_report_xlsx <- function(tables, file, dictionary = wca_dictionary(), pct_regex = "PCT",
                              notes = NULL) {
  if (is.null(names(tables)) || any(names(tables) == "")) wca_abort("write_report_xlsx(): tables must be named")
  tables <- purrr::map(tables, \(x) as.data.frame(x, stringsAsFactors = FALSE))
  wb <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(textDecoration = "bold", halign = "left", wrapText = TRUE)
  pct_style <- openxlsx::createStyle(numFmt = "0.00%")
  add_sheet <- function(name, d) {
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, d, headerStyle = header_style)
    openxlsx::freezePane(wb, name, firstRow = TRUE)
    pct_cols <- which(str_detect(names(d), pct_regex))
    if (length(pct_cols) > 0 && nrow(d) > 0) {
      openxlsx::addStyle(wb, name, pct_style, rows = seq_len(nrow(d)) + 1, cols = pct_cols,
        gridExpand = TRUE, stack = TRUE)
    }
    widths <- pmin(pmax(nchar(names(d)) + 2, 8), 40)
    openxlsx::setColWidths(wb, name, cols = seq_along(d), widths = widths)
  }
  purrr::iwalk(tables, \(d, name) add_sheet(name, d))
  if (!is.null(dictionary)) {
    add_sheet("DataDictionary", build_dictionary_sheet(tables, dictionary, notes))
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  invisible(wb)
}

#' Assemble the DataDictionary sheet for a set of tables
#'
#' @inheritParams write_report_xlsx
#' @return Data frame: Sheet, Column, Description.
build_dictionary_sheet <- function(tables, dictionary, notes = NULL) {
  cols <- purrr::imap(tables, \(d, sheet) tibble(Sheet = sheet, Column = names(d))) |>
    bind_rows()
  dict <- dictionary |> distinct(Sheet, Column, .keep_all = TRUE)
  generic <- dict |> filter(Sheet == "*") |> select(Column, Description_any = Description)
  described <- cols |>
    left_join(dict |> filter(Sheet != "*"), by = c("Sheet", "Column")) |>
    left_join(generic, by = "Column") |>
    mutate(Description = coalesce(Description, Description_any, "")) |>
    select(Sheet, Column, Description)
  denoms <- purrr::imap(tables, \(d, sheet) {
    den <- attr(d, "denominator")
    if (is.null(den)) NULL else tibble(Sheet = sheet, Column = "(denominator)",
      Description = str_glue("PCT columns are n / {den} samples"))
  }) |> bind_rows()
  note_rows <- if (is.null(notes)) NULL else
    tibble(Sheet = "(notes)", Column = names(notes), Description = unname(notes))
  bind_rows(note_rows, denoms, described) |> as.data.frame()
}
