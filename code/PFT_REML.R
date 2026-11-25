# ================================================================
# Single-proportion meta-analysis: PFT + random effects + REML
# Compatible target: R 4.4.2, meta 8.2.0
# Input: Excel workbook, one antimicrobial/outcome per sheet.
# Required columns (case-insensitive): study, Events, Total
# Output: forest PDFs + audit-ready Excel workbook.
# ================================================================

required_packages <- c("meta", "readxl", "openxlsx")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install missing packages: install.packages(c(",
       paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
}


# --------------------- USER SETTINGS ----------------------------
# Run this script from the repository root directory.
# Example: Rscript code/PFT_REML.R data/E.coli.xlsx

args <- commandArgs(trailingOnly = TRUE)

file_path <- if (length(args)) args[1] else "data/E.coli.xlsx"

min_total <- 10L

if (!file.exists(file_path)) {
  stop("Input file not found: ", file_path)
}

species <- tools::file_path_sans_ext(basename(file_path))

out_dir <- file.path("results", "PFT_REML", species)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(out_dir)) {
  stop("Could not create output folder: ", out_dir)
}

# ---------------------- HELPERS ---------------------------------
num_parse <- function(x) {
  # Do not silently convert nonnumeric strings to valid numbers.
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

# The PFT back-transformation requires a sample-size parameter.
# We explicitly use the harmonic mean of the included study sizes,
# for the pooled estimate, CI endpoints, and PI endpoints alike.
# This must be reported in the Methods and kept consistent with plots.
pft_backtransform <- function(x, study_n) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return(NA_real_)
  study_n <- study_n[is.finite(study_n) & study_n > 0]
  if (!length(study_n)) return(NA_real_)
  harmonic_n <- 1 / mean(1 / study_n)
  ans <- as.numeric(meta::backtransf(x, sm = "PFT", n = harmonic_n))
  if (!is.finite(ans)) return(NA_real_)
  min(1, max(0, ans))
}

fmt <- function(x, digits = 4) if (is.finite(x)) round(x, digits) else NA_real_

safe_name <- function(x, index, prefix = "") {
  x <- gsub("[\\[\\]:*?/\\\\]", "_", x)
  x <- gsub("[^[:alnum:]_. -]", "_", x)
  x <- trimws(x)
  # Unique by design; Excel worksheet name <=31 chars.
  substr(sprintf("%02d_%s%s", index, prefix, x), 1, 31)
}

safe_filename <- function(x, index) {
  x <- gsub("[^[:alnum:]_. -]", "_", x)
  substr(sprintf("%02d_%s", index, x), 1, 120)
}

make_summary <- function(sheet, n_original, clean, audit, model = NULL,
                         status = "Success", plot_status = NA_character_) {
  excluded <- audit$Reason != "Included"
  n_na <- sum(grepl("Missing/non-numeric", audit$Reason, fixed = TRUE))
  n_small <- sum(audit$Reason == "Total < 10")
  n_invalid <- sum(audit$Reason %in% c("Invalid Events/Total", "Non-integer count"))

  if (is.null(model)) {
    return(data.frame(
      Sheet_Name = sheet, Studies_Original = n_original,
      Studies_Analyzed = NA_integer_, Studies_Eligible = nrow(clean),
      Removed_Missing = n_na, Removed_Small = n_small,
      Removed_Invalid = n_invalid, Total_Removed = sum(excluded),
      Total_Isolates_Tested = sum(clean$Total), Total_Events = sum(clean$Events),
      Overall_Proportion = NA_real_, Overall_Percent = NA_real_,
      CI_95_Lower_Percent = NA_real_, CI_95_Upper_Percent = NA_real_,
      PI_95_Lower_Percent = NA_real_, PI_95_Upper_Percent = NA_real_,
      Tau2 = NA_real_, I2_Percent = NA_real_, Q = NA_real_, Q_df = NA_real_,
      P_value_heterogeneity = NA_character_,
      Status = status, Forest_Plot_Status = plot_status,
      stringsAsFactors = FALSE))
  }

  bt <- function(x) pft_backtransform(x, model$n)
  pooled <- bt(model$TE.random)
  ci_lo <- bt(model$lower.random)
  ci_hi <- bt(model$upper.random)
  # No meaningful random-effects PI when k < 2.
  pi_lo <- if (model$k >= 2) bt(model$lower.predict) else NA_real_
  pi_hi <- if (model$k >= 2) bt(model$upper.predict) else NA_real_
  data.frame(
    Sheet_Name = sheet, Studies_Original = n_original,
    Studies_Analyzed = model$k, Studies_Eligible = nrow(clean),
    Removed_Missing = n_na, Removed_Small = n_small,
    Removed_Invalid = n_invalid, Total_Removed = sum(excluded),
    Total_Isolates_Tested = sum(model$n), Total_Events = sum(model$event),
    Overall_Proportion = fmt(pooled), Overall_Percent = fmt(pooled * 100, 2),
    CI_95_Lower_Percent = fmt(ci_lo * 100, 2),
    CI_95_Upper_Percent = fmt(ci_hi * 100, 2),
    PI_95_Lower_Percent = fmt(pi_lo * 100, 2),
    PI_95_Upper_Percent = fmt(pi_hi * 100, 2),
    Tau2 = fmt(model$tau2, 6),
    I2_Percent = fmt(model$I2 * 100, 2),
    Q = fmt(model$Q, 4), Q_df = model$df.Q,
    P_value_heterogeneity = if (is.finite(model$pval.Q))
      format.pval(model$pval.Q, eps = 0.001, digits = 3) else NA_character_,
    Status = status, Forest_Plot_Status = plot_status,
    stringsAsFactors = FALSE)
}

sheet_names <- readxl::excel_sheets(file_path)
all_results <- vector("list", length(sheet_names))
names(all_results) <- sheet_names

for (i in seq_along(sheet_names)) {
  sheet_name <- sheet_names[i]
  message(sprintf("[%d/%d] %s", i, length(sheet_names), sheet_name))
  raw <- tryCatch(as.data.frame(readxl::read_excel(file_path, sheet = sheet_name)),
                  error = function(e) e)
  if (inherits(raw, "error")) {
    all_results[[i]] <- list(summary = data.frame(Sheet_Name = sheet_name,
      Status = paste("Read error:", conditionMessage(raw))),
      audit = NULL, clean = NULL, details = NULL)
    next
  }
  names(raw) <- trimws(names(raw))
  idx <- match(tolower(c("study", "events", "total")), tolower(names(raw)))
  if (anyNA(idx)) {
    all_results[[i]] <- list(summary = data.frame(Sheet_Name = sheet_name,
      Status = "Skipped: missing study/Events/Total columns"),
      audit = NULL, clean = NULL, details = NULL)
    next
  }
  names(raw)[idx] <- c("study", "Events", "Total")
  raw$.Original_Row <- seq_len(nrow(raw)) + 1L  # Excel row (header is row 1)
  raw$.Original_Events <- as.character(raw$Events)
  raw$.Original_Total <- as.character(raw$Total)
  raw$Events <- num_parse(raw$Events)
  raw$Total <- num_parse(raw$Total)

  # Mutually exclusive exclusion reasons, applied in priority order.
  reason <- rep("Included", nrow(raw))
  reason[is.na(raw$Events) | is.na(raw$Total) |
         is.na(raw$study) | !nzchar(trimws(as.character(raw$study)))] <-
    "Missing/non-numeric Events, Total, or study"
  valid <- reason == "Included"
  reason[valid & (!is.finite(raw$Events) | !is.finite(raw$Total) |
                  raw$Events < 0 | raw$Total <= 0 |
                  raw$Events > raw$Total)] <- "Invalid Events/Total"
  valid <- reason == "Included"
  reason[valid & (raw$Events != floor(raw$Events) |
                  raw$Total != floor(raw$Total))] <- "Non-integer count"
  valid <- reason == "Included"
  reason[valid & raw$Total < min_total] <- "Total < 10"

  audit <- data.frame(Excel_Row = raw$.Original_Row,
                      Study = as.character(raw$study),
                      Original_Events = raw$.Original_Events,
                      Original_Total = raw$.Original_Total,
                      Reason = reason, stringsAsFactors = FALSE)
  clean <- raw[reason == "Included", , drop = FALSE]
  clean <- clean[, !grepl("^\\.Original_", names(clean)), drop = FALSE]
  clean$.Original_Row <- NULL
  if (!nrow(clean)) {
    all_results[[i]] <- list(
      summary = make_summary(sheet_name, nrow(raw), clean, audit,
                             status = "Skipped: no eligible data"),
      audit = audit, clean = clean, details = NULL)
    next
  }

  model <- tryCatch(meta::metaprop(
    event = Events, n = Total, studlab = study, data = clean,
    sm = "PFT", method.tau = "REML", common = FALSE,
    random = TRUE, prediction = nrow(clean) >= 2,
    level.predict = 0.95, backtransf = TRUE,
    method.random.ci = "classic"
  ), error = function(e) e)
  if (inherits(model, "error")) {
    all_results[[i]] <- list(
      summary = make_summary(sheet_name, nrow(raw), clean, audit,
                             status = paste("Model error:", conditionMessage(model))),
      audit = audit, clean = clean, details = NULL)
    next
  }

  details <- data.frame(
    Study = model$studlab, Events = model$event, Total = model$n,
    Raw_Proportion = round(model$event / model$n, 4),
    Raw_Percent = round(100 * model$event / model$n, 2),
    Weight_Random = if (!is.null(model$w.random))
      round(model$w.random, 4) else NA_real_,
    stringsAsFactors = FALSE)

  # Plot errors are independent of the statistical results.
  pdf_path <- file.path(out_dir, paste0(safe_filename(sheet_name, i), "_forest.pdf"))
  plot_status <- tryCatch({
    grDevices::pdf(pdf_path, width = 14,
                   height = max(10, 0.30 * model$k + 6),
                   family = "serif", onefile = TRUE)
    tryCatch({
      meta::forest(model,
        layout = "Revman",
        leftcols = c("studlab", "event", "n", "effect.ci"),
        leftlabs = c("Study", "Resistant", "Total", "Proportion [95% CI]"),
        rightcols = FALSE,
        common = FALSE, random = TRUE,
        prediction = model$k >= 2,
        backtransf = TRUE,
        overall.hetstat = TRUE,
        xlab = "Proportion", digits = 2,
        col.diamond.random = "#669999",
        col.diamond.lines.random = "#669999",
        col.square = "#336699",
        col.square.lines = "#336699")
    }, finally = grDevices::dev.off())
    "Saved"
  }, error = function(e) paste("Plot error:", conditionMessage(e)))

  all_results[[i]] <- list(
    summary = make_summary(sheet_name, nrow(raw), clean, audit,
                           model = model, plot_status = plot_status),
    audit = audit, clean = clean, details = details)
  message("  k = ", model$k, "; forest: ", plot_status)
}

# ---------------------- EXCEL EXPORT ----------------------------
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Summary_All_Sheets")
# bind_rows without dplyr: align varying status-only summaries.
summaries <- lapply(all_results, `[[`, "summary")
all_cols <- unique(unlist(lapply(summaries, names)))
summaries <- lapply(summaries, function(x) {
  for (nm in setdiff(all_cols, names(x))) x[[nm]] <- NA
  x[all_cols]
})
summary_df <- do.call(rbind, summaries)
openxlsx::writeData(wb, "Summary_All_Sheets", summary_df)

for (i in seq_along(all_results)) {
  r <- all_results[[i]]
  if (!is.null(r$details)) {
    nm <- safe_name(sheet_names[i], i, "Det_")
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, r$details)
  }
  if (!is.null(r$clean)) {
    nm <- safe_name(sheet_names[i], i, "Data_")
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, r$clean)
  }
  if (!is.null(r$audit)) {
    nm <- safe_name(sheet_names[i], i, "Audit_")
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, r$audit)
  }
}

openxlsx::addWorksheet(wb, "Methods_Readme")
methods <- data.frame(
  Item = c("Input file", "Model", "Transformation", "Tau2 estimator",
           "Random-effects CI", "Prediction interval", "PFT inverse",
           "Eligibility", "Exclusion audit", "R version", "meta version",
           "Note on PI", "Note on units"),
  Setting = c(basename(file_path), "Single proportion, random effects",
              "Freeman-Tukey double arcsine (PFT)", "REML",
              "Classic normal approximation", "95% (k >= 2)",
              "Harmonic mean of included study sample sizes",
              paste0("Total >= ", min_total, "; integer Events/Total; 0 <= Events <= Total"),
              "One exclusive exclusion reason per input row",
              R.version.string, as.character(utils::packageVersion("meta")),
              "PI not reported for k < 2; interpret cautiously when k is small",
              "Total_Isolates_Tested sums tested isolates, not unique patients"),
  stringsAsFactors = FALSE)
openxlsx::writeData(wb, "Methods_Readme", methods)

result_path <- file.path(out_dir, "All_Sheets_PFT_REML_Results_with_PI.xlsx")
openxlsx::saveWorkbook(wb, result_path, overwrite = TRUE)
message("Finished. Results: ", result_path)
message("Success: ", sum(summary_df$Status == "Success", na.rm = TRUE),
        "; other: ", sum(summary_df$Status != "Success", na.rm = TRUE))
