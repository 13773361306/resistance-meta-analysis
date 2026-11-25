# ================================================================
# SINGLE-PROPORTION META-ANALYSIS: BINOMIAL GLMM (LOGIT), ML
# Original-criteria-preserving GLMM; separate from PFT-REML primary analysis
# Input: one Excel workbook; one antibiotic/outcome per worksheet
# Tested logically against meta 8.x interface; RUN IN R to verify.
# ================================================================


# 1. USER SETTINGS -- RELATIVE PATHS
# Run this script from the repository root directory.
# Example: Rscript code/GLMM_ML.R data/E.coli.xlsx

args <- commandArgs(trailingOnly = TRUE)

file_path <- if (length(args)) args[1] else "data/E.coli.xlsx"

min_total <- 10L

if (!file.exists(file_path)) {
  stop("Input file not found: ", file_path)
}

species <- tools::file_path_sans_ext(basename(file_path))

output_dir <- file.path("results", "GLMM_ML", species)

# 2. REQUIRED PACKAGES
needed <- c("meta", "metafor", "readxl", "openxlsx")
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install missing packages: install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(output_dir)) stop("Could not create output folder: ", output_dir)
plot_dir <- file.path(output_dir, "Forest_Plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(plot_dir)) stop("Could not create forest plot folder: ", plot_dir)

# 3. UTILITIES
first_number <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0L || !is.finite(v[1])) return(NA_real_)
  v[1]
}

inv_logit <- function(x) {
  v <- first_number(x)
  if (is.na(v)) return(NA_real_)
  stats::plogis(v)
}

p_display <- function(x) {
  p <- first_number(x)
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

find_column <- function(nms, candidates) {
  idx <- which(tolower(trimws(nms)) %in% tolower(candidates))
  if (!length(idx)) return(NA_character_)
  nms[idx[1]]
}

numeric_strict <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  z <- trimws(as.character(x))
  z[z == ""] <- NA_character_
  suppressWarnings(as.numeric(z))
}

make_excel_name <- function(prefix, index) {
  # Unique by index; always within Excel's 31-character limit
  paste0(prefix, sprintf("%03d", index))
}

make_summary <- function(sheet, original_n, analyzed_n, audit,
                         total_n = NA_real_, total_events = NA_real_,
                         pooled = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                         pi_lo = NA_real_, pi_hi = NA_real_,
                         tau2 = NA_real_, i2 = NA_real_,
                         q = NA_real_, df_q = NA_real_, p_q = NA_character_,
                         status = "Not analyzed", note = "") {
  counts <- table(factor(audit$Exclusion_Reason,
                         levels = c("Included", "Missing/non-numeric", "Invalid count", "Total < 10", "Missing study ID")))
  data.frame(
    Sheet_Name = sheet, Model = "Binomial GLMM (logit), ML",
    Studies_Original = original_n, Studies_Analyzed = analyzed_n,
    Studies_Removed_NA = as.integer(counts["Missing/non-numeric"]),
    Studies_Removed_Invalid = as.integer(counts["Invalid count"]),
    Studies_Removed_Small = as.integer(counts["Total < 10"]),
    Studies_Removed_Missing_ID = as.integer(counts["Missing study ID"]),
    Total_Removed = sum(audit$Exclusion_Reason != "Included"),
    Total_Participants = total_n, Total_Events = total_events,
    Overall_Proportion = round(pooled, 6), Overall_Percent = round(pooled * 100, 3),
    CI_95_Lower = round(ci_lo, 6), CI_95_Upper = round(ci_hi, 6),
    CI_95_Lower_Percent = round(ci_lo * 100, 3),
    CI_95_Upper_Percent = round(ci_hi * 100, 3),
    PI_95_Lower = round(pi_lo, 6), PI_95_Upper = round(pi_hi, 6),
    PI_95_Lower_Percent = round(pi_lo * 100, 3),
    PI_95_Upper_Percent = round(pi_hi * 100, 3),
    Tau2 = round(tau2, 6), I2 = round(i2, 3),
    Q_statistic = round(q, 4), df = df_q,
    P_value_heterogeneity = p_q, Status = status, Notes = note,
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

analyze_sheet <- function(sheet, index) {
  cat("\n--- ", sheet, " ---\n", sep = "")
  raw <- as.data.frame(readxl::read_excel(file_path, sheet = sheet),
                       stringsAsFactors = FALSE)
  names(raw) <- trimws(names(raw))
  events_col <- find_column(names(raw), c("Events", "Event", "Resistant", "Resistance", "R"))
  total_col  <- find_column(names(raw), c("Total", "N", "Tested", "Number tested", "Sample size"))
  study_col  <- find_column(names(raw), c("study", "Study ID", "Author", "Reference", "studlab"))
  audit <- data.frame(Row_in_Source = seq_len(nrow(raw)),
                      Exclusion_Reason = rep("Included", nrow(raw)),
                      stringsAsFactors = FALSE)
  if (anyNA(c(events_col, total_col, study_col))) {
    audit$Exclusion_Reason <- "Missing/non-numeric"
    return(list(summary = make_summary(sheet, nrow(raw), 0L, audit,
                    status = "Skipped: missing columns",
                    note = paste("Available columns:", paste(names(raw), collapse = ", "))),
                audit = cbind(audit, raw), data = raw, details = NULL, meta_object = NULL))
  }

  # Preserve source columns; standardized columns are added separately.
  dat <- raw
  dat$study <- as.character(raw[[study_col]])
  dat$Events <- numeric_strict(raw[[events_col]])
  dat$Total <- numeric_strict(raw[[total_col]])
  # EXACTLY match the user's original GLMM script's inclusion criteria.
  # Do not add integer-only or non-missing-study-ID exclusions.
  missing_count <- is.na(dat$Events) | is.na(dat$Total)
  invalid_count <- !missing_count & (
    dat$Total <= 0 | dat$Events < 0 | dat$Events > dat$Total)
  small_count <- !missing_count & !invalid_count & dat$Total < min_total
  audit$Exclusion_Reason[missing_count] <- "Missing/non-numeric"
  audit$Exclusion_Reason[invalid_count] <- "Invalid count"
  audit$Exclusion_Reason[small_count] <- "Total < 10"
  audit$Study <- dat$study
  audit$Events <- dat$Events
  audit$Total <- dat$Total
  clean <- dat[audit$Exclusion_Reason == "Included", , drop = FALSE]
  k <- nrow(clean)
  n_sum <- sum(clean$Total)
  event_sum <- sum(clean$Events)
  base <- function(status, note = "", analyzed = 0L) {
    make_summary(sheet, nrow(raw), analyzed, audit, n_sum, event_sum,
                 status = status, note = note)
  }
  if (k == 0L) return(list(summary = base("Skipped: no valid rows"),
                           audit = audit, data = clean, details = NULL, meta_object = NULL))
  # As in the original code, attempt the GLMM fit even for k=1 or
  # all-zero/all-event sheets; report actual fit errors or unavailable PI.
  # No artificial 0-100% interval is ever inserted.

  # Fit failures do not prevent processing subsequent worksheets.
  warnings_seen <- character()
  fit <- tryCatch(withCallingHandlers(
    meta::metaprop(event = Events, n = Total, studlab = study,
                   data = clean, sm = "PLOGIT", method = "GLMM",
                   method.tau = "ML", common = FALSE, random = TRUE,
                   overall = TRUE, prediction = TRUE, level.predict = 0.95,
                   backtransf = TRUE, warn = TRUE),
    warning = function(w) {
      warnings_seen <<- unique(c(warnings_seen, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(summary = base("Error: model fit", conditionMessage(fit), analyzed = k),
                audit = audit, data = clean, details = NULL, meta_object = NULL))
  }

  pooled <- inv_logit(fit$TE.random)
  ci_lo <- inv_logit(fit$lower.random)
  ci_hi <- inv_logit(fit$upper.random)
  pi_lo <- inv_logit(fit$lower.predict)
  pi_hi <- inv_logit(fit$upper.predict)
  tau2 <- first_number(fit$tau2)
  i2 <- first_number(fit$I2) * 100
  q <- first_number(fit$Q)
  df_q <- first_number(fit$df.Q)
  p_q <- p_display(fit$pval.Q)
  note <- paste(warnings_seen, collapse = " | ")
  status <- "Success"
  if (!all(is.finite(c(pooled, ci_lo, ci_hi)))) {
    status <- "Review: nonfinite pooled estimate/CI"
  }
  if (!all(is.finite(c(pi_lo, pi_hi)))) {
    note <- paste(c(note, "Prediction interval unavailable"), collapse = " | ")
  }
  if (is.finite(pi_lo) && is.finite(pi_hi) &&
      (pi_lo <= 0.000001 || pi_hi >= 0.999999)) {
    note <- paste(c(note, "PI touches boundary; interpret cautiously"), collapse = " | ")
  }
  summary <- make_summary(sheet, nrow(raw), k, audit, n_sum, event_sum,
                          pooled, ci_lo, ci_hi, pi_lo, pi_hi,
                          tau2, i2, q, df_q, p_q, status, note)
  details <- data.frame(Study = fit$studlab, Events = fit$event, Total = fit$n,
                        Raw_Proportion = round(fit$event / fit$n, 6),
                        Raw_Percent = round(100 * fit$event / fit$n, 3),
                        stringsAsFactors = FALSE)

  # Plotting is independent of model fitting; a graphics error does not discard estimates.
  safe_name <- gsub("[^[:alnum:]_-]+", "_", sheet)
  plot_file <- file.path(plot_dir, sprintf("%03d_%s_GLMM_forest.pdf", index,
                                            substr(safe_name, 1, 55)))
  # Match the PFT-REML forest layout: numeric study 95% CIs on the left,
  # no right-hand columns, and an explicit numeric 95% PI label.
  pi_available <- is.finite(pi_lo) && is.finite(pi_hi) && k >= 2L
  pi_label <- if (pi_available) {
    sprintf("95%% PI: %.2f%% to %.2f%%", 100 * pi_lo, 100 * pi_hi)
  } else "95% PI: not estimable"
  plot_error <- tryCatch({
    grDevices::pdf(plot_file, width = 14,
                   height = min(70, max(9, 0.35 * k + 8)), family = "Helvetica")
    tryCatch({
      meta::forest(fit,
        layout = "Revman",
        leftcols = c("studlab", "event", "n", "effect.ci"),
        leftlabs = c("Study", "Resistant", "Total", "Proportion [95% CI]"),
        rightcols = FALSE,
        common = FALSE, random = TRUE,
        overall = TRUE, overall.hetstat = TRUE,
        prediction = pi_available,
        text.predict = if (pi_available) pi_label else NULL,
        backtransf = TRUE,
        xlab = "Proportion", xlim = c(0, 1), smlab = "Proportion",
        digits = 2, fontsize = 9.5,
        col.diamond.random = "#669999",
        col.diamond.lines.random = "#669999",
        col.square = "#336699",
        col.square.lines = "#336699",
        col.predict = "#660099")
    }, finally = grDevices::dev.off())
    NULL
  }, error = function(e) conditionMessage(e))
  summary$PDF_Status <- if (is.null(plot_error)) "Saved" else paste("Error:", plot_error)
  summary$PDF_File <- plot_file
  if (!is.null(plot_error)) {
    summary$Notes <- paste(c(summary$Notes, paste("Plot error:", plot_error)), collapse = " | ")
    cat("  Plot error: ", plot_error, "\n", sep = "")
  }
  cat(sprintf("  k=%d; pooled=%.3f%%; 95%% CI=[%.3f, %.3f]%%\n",
              k, 100 * pooled, 100 * ci_lo, 100 * ci_hi))
  list(summary = summary, audit = audit, data = clean, details = details, meta_object = fit)
}

# 4. RUN ALL WORKSHEETS
sheet_names <- readxl::excel_sheets(file_path)
all_results <- lapply(seq_along(sheet_names), function(i) {
  tryCatch(analyze_sheet(sheet_names[i], i), error = function(e) {
    message("Unexpected error in ", sheet_names[i], ": ", conditionMessage(e))
    audit <- data.frame(Exclusion_Reason = character())
    list(summary = make_summary(sheet_names[i], NA_integer_, 0L, audit,
                                status = "Error: worksheet", note = conditionMessage(e)),
         audit = data.frame(Error = conditionMessage(e)), data = NULL, details = NULL, meta_object = NULL)
  })
})
names(all_results) <- sheet_names

# Each worksheet now retains its OWN fitted model. Never rely on an old
# global `meta.glmm` object from a previous analysis.
model_check <- data.frame(
  Sheet = names(all_results),
  Studies = vapply(all_results, function(x) {
    if (is.null(x$meta_object)) return(NA_integer_)
    as.integer(x$meta_object$k)
  }, integer(1)),
  Events = vapply(all_results, function(x) {
    if (is.null(x$meta_object)) return(NA_real_)
    sum(x$meta_object$event)
  }, numeric(1)),
  Total = vapply(all_results, function(x) {
    if (is.null(x$meta_object)) return(NA_real_)
    sum(x$meta_object$n)
  }, numeric(1)),
  Status = vapply(all_results, function(x) as.character(x$summary$Status[1]), ""),
  stringsAsFactors = FALSE
)
cat("\n=== WORKSHEET / MODEL CHECK ===\n")
print(model_check, row.names = FALSE)

# Optional: explicitly redraw ONE worksheet, using its saved model.
# Example after running the script: plot_glmm_sheet("CR-EC")
plot_glmm_sheet <- function(sheet, filename = NULL) {
  if (!sheet %in% names(all_results)) {
    stop("Worksheet not found: ", sheet,
         ". Available: ", paste(names(all_results), collapse = ", "))
  }
  fit <- all_results[[sheet]]$meta_object
  if (is.null(fit)) stop("No successful GLMM model for sheet: ", sheet,
                         ". Check model_check and Summary_All_Sheets.")
  if (is.null(filename)) {
    safe_sheet <- gsub("[^[:alnum:]_-]+", "_", sheet)
    filename <- file.path(plot_dir, paste0(safe_sheet, "_GLMM_selected.pdf"))
  }
  lo <- inv_logit(fit$lower.predict)
  hi <- inv_logit(fit$upper.predict)
  pi_ok <- is.finite(lo) && is.finite(hi) && fit$k >= 2
  pi_text <- if (pi_ok) sprintf("95%% PI: %.2f%% to %.2f%%", 100 * lo, 100 * hi) else NULL
  grDevices::pdf(filename, width = 14, height = min(70, max(9, 0.35 * fit$k + 8)),
                 family = "Helvetica")
  on.exit(grDevices::dev.off(), add = TRUE)
  meta::forest(fit, layout = "Revman",
               leftcols = c("studlab", "event", "n", "effect.ci"),
               leftlabs = c("Study", "Resistant", "Total", "Proportion [95% CI]"),
               rightcols = FALSE, common = FALSE, random = TRUE,
               overall = TRUE, overall.hetstat = TRUE,
               prediction = pi_ok, text.predict = pi_text,
               backtransf = TRUE, fontsize = 9.5, digits = 2,
               xlab = "Proportion", xlim = c(0, 1))
  message("Selected worksheet: ", sheet, "; k=", fit$k,
          "; events=", sum(fit$event), "; total=", sum(fit$n),
          "; PDF: ", filename)
  invisible(filename)
}

# 5. EXCEL EXPORT: summary, audit trail, included rows, study details
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Summary_All_Sheets")
# Harmonize columns: successful fits include PDF_Status/PDF_File,
# while skipped/failed fits do not. Base rbind otherwise errors.
summaries <- lapply(all_results, `[[`, "summary")
all_summary_columns <- unique(unlist(lapply(summaries, names)))
summaries <- lapply(summaries, function(x) {
  for (nm in setdiff(all_summary_columns, names(x))) x[[nm]] <- NA_character_
  x[, all_summary_columns, drop = FALSE]
})
summary_table <- do.call(rbind, summaries)
openxlsx::writeData(wb, "Summary_All_Sheets", summary_table)
openxlsx::freezePane(wb, "Summary_All_Sheets", firstRow = TRUE)

for (i in seq_along(all_results)) {
  item <- all_results[[i]]
  for (kind in c("audit", "data", "details")) {
    value <- item[[kind]]
    if (is.null(value)) next
    prefix <- switch(kind, audit = "Audit_", data = "Data_", details = "Det_")
    nm <- make_excel_name(prefix, i)
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, value)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
  }
}
openxlsx::addWorksheet(wb, "Sheet_Index")
openxlsx::writeData(wb, "Sheet_Index", data.frame(
  Source_Sheet = sheet_names,
  Audit_Sheet = vapply(seq_along(sheet_names), function(i) make_excel_name("Audit_", i), ""),
  Data_Sheet = vapply(seq_along(sheet_names), function(i) make_excel_name("Data_", i), ""),
  Detail_Sheet = vapply(seq_along(sheet_names), function(i) {
    if (is.null(all_results[[i]]$details)) NA_character_ else make_excel_name("Det_", i)
  }, ""), stringsAsFactors = FALSE))

openxlsx::addWorksheet(wb, "Methods_Readme")
methods <- data.frame(
  Item = c("Analysis", "Distribution", "Link", "Between-study variance",
           "Primary comparator", "CI", "PI", "Minimum tested sample",
           "Invalid count rules", "Boundary-only datasets", "R version",
           "meta version", "metafor version", "Input file", "Output folder"),
  Setting = c("Random-effects single-proportion GLMM (meta::metaprop)",
              "Binomial", "Logit", "ML", "PFT-REML (separate primary analysis)",
              "95%, inverse-logit of fitted logit-scale bounds",
              "95%, inverse-logit of fitted logit-scale bounds; not reported when unavailable",
              paste0("Total >= ", min_total),
              "Original rules: NA Events/Total; Total<10; Total<=0; Events<0; Events>Total",
              "GLMM attempted as in original code; failed fits or unavailable PI reported without fabrication",
              R.version.string, as.character(utils::packageVersion("meta")),
              as.character(utils::packageVersion("metafor")),
              normalizePath(file_path), normalizePath(output_dir)),
  stringsAsFactors = FALSE)
openxlsx::writeData(wb, "Methods_Readme", methods)

out_xlsx <- file.path(output_dir, "All_Sheets_Meta_Analysis_GLMM_ML_with_PI.xlsx")
openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat("\nCompleted. Excel: ", out_xlsx, "\nForest plots: ", plot_dir, "\n", sep = "")
cat("Success: ", sum(summary_table$Status == "Success"),
    "; review/failed/skipped: ", sum(summary_table$Status != "Success"), "\n", sep = "")

# 6. FINISH / QUICK VERIFICATION
cat("\nTo redraw a specific worksheet, run e.g.: plot_glmm_sheet(\"CR-EC\")\n")
cat("Check model_check and the PDF_Status column in Summary_All_Sheets.\n")
