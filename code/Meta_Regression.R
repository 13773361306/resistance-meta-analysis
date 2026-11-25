# ==============================================================
# FINAL: Stata 18-compatible single-proportion meta-regression
# Revised: centered title and ONE-LINE caption with k, beta/OR 95% CIs, p, R²
# Input sheets: study | Year | Events | Total
# Reference: Stata 'meta set logit_es logit_se; meta regress year'
# Continuity correction ONLY for Events==0 or Events==Total.
# Model: metafor::rma.uni(..., method='REML', test='z')
# ==============================================================

# ---------- 0. Settings: relative paths -----------------------
# Run from the repository root directory.
# Example: Rscript code/Meta_Regression.R data/meta_regression.xlsx

args <- commandArgs(trailingOnly = TRUE)

file_path <- if (length(args)) args[1] else "data/meta_regression.xlsx"

dataset_name <- tools::file_path_sans_ext(basename(file_path))

output_dir <- file.path("results", "Meta_Regression", dataset_name)
min_total <- 1L          # Stata PDF did NOT filter Total >= 10
min_studies <- 3L        # technical minimum; set 10L for conservative policy
plot_bubbles <- TRUE
# Set TRUE only if the source Excel has 'Publication year' instead of 'Year'
# and you INTENTIONALLY want to regress on publication year.
allow_publication_year <- FALSE

# ---------- 1. Dependencies -----------------------------------
packages <- c("readxl", "metafor", "ggplot2", "openxlsx")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install missing packages first: install.packages(c(",
       paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
}
if (!file.exists(file_path)) stop("Excel file not found: ", file_path)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
plot_dir <- file.path(output_dir, "Bubble_Plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---------- 2. Helpers ----------------------------------------
clean_numeric <- function(x) {
  z <- trimws(gsub("\u3000", "", as.character(x), fixed = TRUE))
  z[z %in% c("", "NA", "N/A", "NULL", "-", ".")] <- NA_character_
  suppressWarnings(as.numeric(z))
}

find_column <- function(nms, choices) {
  pos <- match(tolower(choices), tolower(trimws(nms)))
  pos <- pos[!is.na(pos)]
  if (!length(pos)) return(NA_character_)
  nms[pos[1]]
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]", "_", x)
  substr(x, 1, 90)
}

fmt <- function(x, digits = 6) {
  if (!length(x) || !is.finite(x[1])) return(NA_real_)
  round(as.numeric(x[1]), digits)
}

# Stata PDF reference values, only used for comparison; NOT imposed on models
stata_reference <- data.frame(
  Sheet = c("3GCR-EC", "CR-EC", "3GCR-Klebsiella", "CR-Klebsiella"),
  Stata_k = c(81L, 77L, 55L, 53L),
  Stata_Beta = c(0.1355, 0.140872, 0.0638106, 0.1368842),
  Stata_SE = c(0.0353681, 0.0374705, 0.0361385, 0.0367681),
  Stata_OR = c(1.145109, 1.151277, 1.065890, 1.146695),
  Stata_Beta_CI_Lower = c(0.0661799, 0.0674311, -0.0070196, 0.06482),
  Stata_Beta_CI_Upper = c(0.2048201, 0.2143129, 0.1346407, 0.2089483),
  Stata_OR_CI_Lower = c(1.068419, 1.069757, 0.993005, 1.066967),
  Stata_OR_CI_Upper = c(1.227304, 1.239010, 1.144126, 1.232381),
  Stata_Tau2 = c(2.443, 1.882, 1.679, 1.313),
  stringsAsFactors = FALSE
)

# ---------- 3. Batch analysis ---------------------------------
sheets <- readxl::excel_sheets(file_path)
summary_list <- list()
detail_list <- list()
audit_list <- list()
model_list <- list()

for (sheet in sheets) {
  cat("\n--- ", sheet, " ---\n", sep = "")
  raw <- tryCatch(readxl::read_excel(file_path, sheet = sheet),
                  error = function(e) e)
  if (inherits(raw, "error")) {
    summary_list[[sheet]] <- data.frame(Sheet = sheet, Status = "Read error",
                                         Note = conditionMessage(raw))
    next
  }
  names(raw) <- trimws(names(raw))
  sc <- find_column(names(raw), c("study"))
  yc <- find_column(names(raw), c("Year"))
  # Case-insensitive matching makes 'year' equivalent to 'Year'.
  if (is.na(yc) && allow_publication_year)
    yc <- find_column(names(raw), c("Publication year"))
  ec <- find_column(names(raw), c("Events"))
  nc <- find_column(names(raw), c("Total"))
  if (anyNA(c(sc, yc, ec, nc))) {
    note <- paste("Required columns missing; found:", paste(names(raw), collapse = ", "),
                  "; set allow_publication_year=TRUE only if appropriate")
    summary_list[[sheet]] <- data.frame(Sheet = sheet, Status = "Missing columns", Note = note)
    cat(note, "\n")
    next
  }

  dat <- data.frame(
    Row_Excel = seq_len(nrow(raw)) + 1L,
    study = trimws(as.character(raw[[sc]])),
    Year = clean_numeric(raw[[yc]]),
    Events = clean_numeric(raw[[ec]]),
    Total = clean_numeric(raw[[nc]]),
    stringsAsFactors = FALSE
  )
  dat$study[is.na(dat$study) | dat$study == ""] <- NA_character_

  # Stata's drop if missing(events,total,year) is the primary rule.
  # Other checks guard against invalid binomial counts and unusable labels.
  dat$Reason <- "Included"
  miss <- !is.finite(dat$Year) | !is.finite(dat$Events) | !is.finite(dat$Total)
  dat$Reason[miss] <- "Missing/non-numeric Year, Events or Total"
  bad <- !miss & (dat$Total <= 0 | dat$Events < 0 |
                  dat$Events > dat$Total |
                  dat$Events != floor(dat$Events) |
                  dat$Total != floor(dat$Total))
  dat$Reason[bad] <- "Invalid binomial counts"
  small <- !miss & !bad & dat$Total < min_total
  dat$Reason[small] <- paste0("Total < ", min_total)
  no_id <- dat$Reason == "Included" & is.na(dat$study)
  dat$Reason[no_id] <- "Missing study ID"

  keep <- dat[dat$Reason == "Included", , drop = FALSE]
  audit_list[[sheet]] <- data.frame(Sheet = sheet, dat, stringsAsFactors = FALSE)
  n_original <- nrow(dat)
  n_clean <- nrow(keep)
  cat("Original=", n_original, "; included=", n_clean, "; removed=", n_original-n_clean, "\n", sep="")

  if (n_clean > 0) {
    # EXACT Stata PDF formulas: boundary studies only
    boundary <- keep$Events == 0 | keep$Events == keep$Total
    keep$Events_adj <- keep$Events + ifelse(boundary, 0.5, 0)
    keep$Total_adj <- keep$Total + ifelse(boundary, 1, 0)
    keep$yi <- log(keep$Events_adj / (keep$Total_adj - keep$Events_adj))
    keep$sei <- sqrt(1 / keep$Events_adj + 1 / (keep$Total_adj - keep$Events_adj))
    keep$vi <- keep$sei^2
    keep$Raw_Percent <- 100 * keep$Events / keep$Total
    keep$Corrected_Boundary <- boundary
    detail_list[[sheet]] <- data.frame(Sheet = sheet, keep, stringsAsFactors = FALSE)
  }

  base <- data.frame(Sheet = sheet, Year_Source_Column = yc,
                     Original_k = n_original, Included_k = n_clean,
                     Excluded_k = n_original - n_clean,
                     Events_Sum = if (n_clean) sum(keep$Events) else NA_real_,
                     Total_Sum = if (n_clean) sum(keep$Total) else NA_real_,
                     Boundary_Corrected_k = if (n_clean) sum(keep$Corrected_Boundary) else NA_integer_,
                     stringsAsFactors = FALSE)

  if (n_clean < min_studies || length(unique(keep$Year)) < 2L) {
    base$Status <- "Skipped"
    base$Note <- if (n_clean < min_studies) "Too few studies" else "No Year variation"
    summary_list[[sheet]] <- base
    next
  }

  # Uncentered Year reproduces Stata's intercept as well as slope.
  fit <- tryCatch(
    metafor::rma.uni(yi = keep$yi, vi = keep$vi,
                     mods = ~ Year, data = keep,
                     method = "REML", test = "z"),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    base$Status <- "Fit error"
    base$Note <- conditionMessage(fit)
    summary_list[[sheet]] <- base
    cat("Fit error: ", base$Note, "\n", sep="")
    next
  }
  model_list[[sheet]] <- fit
  b <- as.numeric(coef(fit))
  beta <- b[2]
  ci_lb <- as.numeric(fit$ci.lb[2])
  ci_ub <- as.numeric(fit$ci.ub[2])

  null_fit <- tryCatch(metafor::rma.uni(yi = keep$yi, vi = keep$vi,
                                       method = "REML", test = "z"),
                       error = function(e) NULL)
  pseudo_r2 <- if (!is.null(null_fit) && is.finite(null_fit$tau2) && null_fit$tau2 > 0)
    100 * (null_fit$tau2 - fit$tau2) / null_fit$tau2 else NA_real_

  base$Status <- "Success"
  base$Note <- ""
  base$Intercept <- fmt(b[1])
  base$Beta_Year <- fmt(beta)
  base$Beta_SE <- fmt(fit$se[2])
  base$Beta_CI_Lower <- fmt(ci_lb)
  base$Beta_CI_Upper <- fmt(ci_ub)
  base$OR_Year <- fmt(exp(beta))
  base$OR_CI_Lower <- fmt(exp(ci_lb))
  base$OR_CI_Upper <- fmt(exp(ci_ub))
  base$P_Year <- as.numeric(fit$pval[2])
  base$Wald_Chi2 <- fmt(fit$QM)
  base$P_Wald <- as.numeric(fit$QMp)
  base$Tau2_Residual <- fmt(fit$tau2)
  base$I2_Residual_Percent <- fmt(fit$I2, 3)
  base$Q_Residual <- fmt(fit$QE, 3)
  base$Q_Residual_df <- fit$k - fit$p
  base$Pseudo_R2_Percent <- fmt(pseudo_r2, 3)
  summary_list[[sheet]] <- base
  cat(sprintf("k=%d; beta/year=%.6f (95%% CI %.6f to %.6f); OR/year=%.6f (95%% CI %.6f to %.6f); p=%.6g; tau2=%.4f\n",
              n_clean, beta, ci_lb, ci_ub, exp(beta), exp(ci_lb), exp(ci_ub), fit$pval[2], fit$tau2))

  if (plot_bubbles) {
    grid_year <- seq(min(keep$Year), max(keep$Year), length.out = 120)
    pred <- tryCatch(predict(fit, newmods = matrix(grid_year, ncol = 1)),
                     error = function(e) NULL)
    if (!is.null(pred)) {
      trend <- data.frame(Year = grid_year,
                          Predicted = 100 * plogis(pred$pred),
                          Lower = 100 * plogis(pred$ci.lb),
                          Upper = 100 * plogis(pred$ci.ub))
      # Approximate random-effects inverse-variance weights for visualization.
      keep$Weight <- 1 / (keep$vi + fit$tau2)
      # Single-line figure caption: k, beta CI, OR CI, p, and R-squared.
      p_label <- if (is.finite(fit$pval[2]) && fit$pval[2] < 0.001) {
        "p < 0.001"
      } else {
        sprintf("p = %.3f", fit$pval[2])
      }
      r2_label <- if (is.finite(pseudo_r2)) sprintf("%.1f%%", pseudo_r2) else "NA"
      one_line_caption <- sprintf(
        "k = %d; beta/year = %.4f (95%% CI %.4f to %.4f); OR/year = %.3f (95%% CI %.3f to %.3f); %s; R² = %s",
        n_clean, beta, ci_lb, ci_ub, exp(beta), exp(ci_lb), exp(ci_ub), p_label, r2_label
      )
      gp <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data = trend,
          ggplot2::aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.18) +
        ggplot2::geom_line(data = trend,
          ggplot2::aes(x = Year, y = Predicted), linewidth = 0.9) +
        ggplot2::geom_point(data = keep,
          ggplot2::aes(x = Year, y = Raw_Percent, size = Weight), alpha = 0.65) +
        ggplot2::scale_size_continuous(range = c(1.8, 7), name = "Approx. RE weight") +
        ggplot2::scale_y_continuous(limits = c(0, 100),
          labels = function(x) paste0(x, "%")) +
        ggplot2::labs(
          title = paste0("Meta-regression: ", sheet),
          x = "Midpoint year of isolate collection",
          y = "Resistance Rate (%)",
          caption = one_line_caption
        ) +
        ggplot2::theme_minimal(base_size = 17) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 24, margin = ggplot2::margin(b = 18)),
          axis.title = ggplot2::element_text(size = 19, face = "bold"),
          axis.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16),
          legend.text = ggplot2::element_text(size = 15),
          plot.caption = ggplot2::element_text(size = 14, hjust = 0.5,
            margin = ggplot2::margin(t = 18)),
          plot.caption.position = "plot",
          plot.margin = ggplot2::margin(22, 45, 25, 45)
        )
      tryCatch(ggplot2::ggsave(
        filename = file.path(plot_dir, paste0(safe_name(sheet), "_meta_regression.pdf")),
        plot = gp, width = 20, height = 11, units = "in",
        device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"),
        error = function(e) warning("PDF failed for ", sheet, ": ", conditionMessage(e)))
    }
  }
}

# ---------- 4. Assemble results safely -------------------------
bind_rows_base <- function(xs) {
  if (!length(xs)) return(data.frame())
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(d) {
    for (col in setdiff(cols, names(d))) d[[col]] <- NA
    d[, cols, drop = FALSE]
  })
  rownames_out <- do.call(rbind, xs)
  rownames(rownames_out) <- NULL
  rownames_out
}
summary_df <- bind_rows_base(summary_list)
detail_df <- bind_rows_base(detail_list)
audit_df <- bind_rows_base(audit_list)

# ---------- 5. Stata benchmark ---------------------------------
benchmark <- merge(stata_reference, summary_df, by = "Sheet", all.x = TRUE)
for (nm in c("Included_k", "Beta_Year", "Beta_SE", "OR_Year", "Tau2_Residual",
             "Beta_CI_Lower", "Beta_CI_Upper", "OR_CI_Lower", "OR_CI_Upper")) {
  if (!nm %in% names(benchmark)) benchmark[[nm]] <- NA_real_
}
benchmark$Delta_k <- benchmark$Included_k - benchmark$Stata_k
benchmark$Delta_Beta <- benchmark$Beta_Year - benchmark$Stata_Beta
benchmark$Delta_SE <- benchmark$Beta_SE - benchmark$Stata_SE
benchmark$Delta_OR <- benchmark$OR_Year - benchmark$Stata_OR
benchmark$Delta_Beta_CI_Lower <- benchmark$Beta_CI_Lower - benchmark$Stata_Beta_CI_Lower
benchmark$Delta_Beta_CI_Upper <- benchmark$Beta_CI_Upper - benchmark$Stata_Beta_CI_Upper
benchmark$Delta_OR_CI_Lower <- benchmark$OR_CI_Lower - benchmark$Stata_OR_CI_Lower
benchmark$Delta_OR_CI_Upper <- benchmark$OR_CI_Upper - benchmark$Stata_OR_CI_Upper
benchmark$Delta_Tau2 <- benchmark$Tau2_Residual - benchmark$Stata_Tau2
benchmark$Note_Compare <- ifelse(is.na(benchmark$Included_k),
  "Sheet not analyzed; check sheet name and Year column",
  ifelse(benchmark$Delta_k != 0,
         "Study count differs: compare source Excel and audit first",
         "Study count matches; inspect effect sizes and model settings"))

# ---------- 6. Excel output ------------------------------------
wb <- openxlsx::createWorkbook()
write_sheet <- function(name, dat) {
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, dat)
  if (ncol(dat) > 0) openxlsx::freezePane(wb, name, firstRow = TRUE)
}
write_sheet("Summary", summary_df)
write_sheet("Stata_Comparison", benchmark)
write_sheet("Study_Details", detail_df)
write_sheet("Data_Audit", audit_df)
write_sheet("Stata_Reference", stata_reference)
write_sheet("Methods", data.frame(
  Item = c("Input", "Year definition", "Data exclusions", "Continuity correction",
           "Effect size", "Standard error", "Model", "Estimator", "Coefficient test",
           "Default minimum Total", "Minimum studies", "R2", "Input file",
           "R version", "metafor version"),
  Value = c("study, Year, Events, Total", "Sampling year; decimal midpoints preserved",
            "Missing/non-numeric Year/Events/Total, invalid binomial counts, missing study ID",
            "Events+0.5, Total+1 ONLY when Events=0 or Events=Total",
            "log(Events_adj/(Total_adj-Events_adj))",
            "sqrt(1/Events_adj+1/(Total_adj-Events_adj))",
            "Random-effects meta-regression, yi ~ Year", "REML", "Wald z",
            as.character(min_total), as.character(min_studies),
            "Pseudo-R2 from null vs full REML tau2; may differ from Stata's R-squared",
            file_path, R.version.string, as.character(utils::packageVersion("metafor"))),
  stringsAsFactors = FALSE))
excel_file <- file.path(output_dir, "Meta_Regression_STATA_FINAL.xlsx")
openxlsx::saveWorkbook(wb, excel_file, overwrite = TRUE)

cat("\nDONE. Excel: ", excel_file, "\n", sep = "")
cat("PDF folder: ", plot_dir, "\n", sep = "")
cat("Successful models: ", sum(summary_df$Status == "Success", na.rm = TRUE), "\n", sep = "")
cat("Stata comparison:\n")
print(benchmark[, c("Sheet", "Stata_k", "Included_k", "Stata_Beta", "Beta_Year",
                    "Beta_CI_Lower", "Beta_CI_Upper", "OR_Year", "OR_CI_Lower", "OR_CI_Upper", "Delta_Beta")],
      row.names = FALSE)
