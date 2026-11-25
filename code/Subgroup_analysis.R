library(meta)
library(readxl)
library(openxlsx)
library(dplyr)

# 1. File paths and settings
file_path <- "D:/Desktop/E.coli.xlsx"
output_dir <- "D:/Desktop/meta/CP_Final"
pdf_font <- "Times New Roman"
minimum_total <- 10
subgroup_vars <- c(
  "AST method", "Etiology", "Income Level",
  "Year Group", "Region", "AST Guideline"
)

if (!file.exists(file_path)) stop("Input file not found: ", file_path)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if (!capabilities("cairo")) stop("Cairo PDF is unavailable.")
if (as.character(packageVersion("meta")) != "8.2.0")
  warning("This forest patch was validated by the user with meta 8.2.0.")

# Preserve the original forest implementation, skipping ONLY study CI backtransformation.
src <- deparse(getFromNamespace("forest.meta", "meta"), width.cutoff = 500)
idx_lower <- grep("lowTE <- backtransf\\(lowTE, sm, npft,", src)
idx_upper <- grep("uppTE <- backtransf\\(uppTE, sm, npft,", src)
stopifnot(length(idx_lower) == 1L, length(idx_upper) == 1L)
src[idx_lower] <- "lowTE <- lowTE"
src[idx_upper] <- "uppTE <- uppTE"
forest_exact <- eval(parse(text = paste(src, collapse = "\n")), envir = asNamespace("meta"))

make_cp_copy <- function(m) {
  stopifnot(inherits(m, "metaprop"), m$sm == "PFT",
            length(m$event) == length(m$n), length(m$n) == length(m$TE))
  cp <- t(vapply(seq_along(m$n), function(i) {
    unname(stats::binom.test(m$event[i], m$n[i], conf.level = 0.95)$conf.int)
  }, numeric(2)))
  stopifnot(all(cp[,1] <= m$event/m$n), all(cp[,2] >= m$event/m$n))
  m$lower <- cp[,1]
  m$upper <- cp[,2]
  m
}

# Match the native PFT inverse using harmonic sample size.
raw_pft <- function(x, n) {
  if (!length(x) || is.na(x[1])) return(NA_real_)
  meta:::backtransf(x[1], sm = "PFT", n = length(n)/sum(1/n))
}

sheet_names <- excel_sheets(file_path)
cat("Available worksheets:", paste(sheet_names, collapse = ", "), "\n")

# Exact Clopper-Pearson estimate for a single study
exact_binomial <- function(events, total, conf.level = 0.95) {
  z <- stats::binom.test(
    x = as.integer(events),
    n = as.integer(total),
    conf.level = conf.level
  )
  c(
    proportion = events / total,
    lower = unname(z$conf.int[1]),
    upper = unname(z$conf.int[2])
  )
}


# Safely combine nested result tables
combine_tables <- function(x, element) {
  out <- unlist(lapply(x, function(s) {
    lapply(s, function(y) y[[element]])
  }), recursive = FALSE)
  out <- out[vapply(out, function(z) is.data.frame(z) && nrow(z) > 0,
                    logical(1))]
  if (!length(out)) return(data.frame())
  dplyr::bind_rows(out)
}

all_results <- list()

# 2. Analyze each worksheet
for (sheet_name in sheet_names) {

  cat("\n========================================\n")
  cat("Analyzing:", sheet_name, "\n")
  cat("========================================\n")

  dat <- read_excel(file_path, sheet = sheet_name)
  names(dat) <- trimws(names(dat))
  names(dat)[tolower(names(dat)) == "ast guideline"] <- "AST Guideline"

  required <- c("Events", "Total", "study")
  if (!all(required %in% names(dat))) {
    cat("Skipped: missing required columns:",
        paste(setdiff(required, names(dat)), collapse = ", "), "\n")
    next
  }

  dat <- dat %>%
    mutate(
      Events = suppressWarnings(as.numeric(Events)),
      Total = suppressWarnings(as.numeric(Total)),
      study = as.character(study)
    ) %>%
    filter(
      !is.na(Events), !is.na(Total), !is.na(study), study != "",
      Total >= minimum_total, Events >= 0, Events <= Total,
      Events == floor(Events), Total == floor(Total)
    )

  if (!nrow(dat)) {
    cat("Skipped: no eligible studies after data cleaning.\n")
    next
  }

  cat("Eligible studies:", nrow(dat), "\n")
  sheet_results <- list()

  for (sg_var in subgroup_vars) {

    cat("\n--- Subgroup variable:", sg_var, "---\n")

    if (!sg_var %in% names(dat)) {
      cat("Skipped: subgroup variable not found.\n")
      next
    }

    dat_sg <- dat %>%
      filter(!is.na(.data[[sg_var]]),
             trimws(as.character(.data[[sg_var]])) != "") %>%
      mutate("{sg_var}" := as.character(.data[[sg_var]]))

    if (!nrow(dat_sg)) {
      cat("Skipped: no valid subgroup data.\n")
      next
    }

    # Fit PFT + REML for multi-study subgroups and the overall estimate.
    # method.ci = "CP" also requests exact CIs for individual study rows.
    meta.sg <- tryCatch(
      metaprop(
        event = Events,
        n = Total,
        data = dat_sg,
        studlab = study,
        sm = "PFT",
        subgroup = dat_sg[[sg_var]],
        overall = TRUE,
        common = FALSE,
        random = TRUE,
        method.tau = "REML",
        method.ci = "CP",
        prediction = TRUE
      ),
      error = function(e) {
        cat("Meta-analysis model error:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (is.null(meta.sg)) next

    # 3. Build subgroup results.
    # Do not rely on meta.sg$bylevs / meta.sg$k.w positional slots.
    # Each subgroup is extracted directly from the source data.
    detail_list <- list()
    subgroup_levels <- unique(as.character(dat_sg[[sg_var]]))

    # Return the first scalar value, or NA if the field is absent/empty.
    get_scalar <- function(x, default = NA_real_) {
      if (is.null(x) || length(x) == 0 || is.na(x[1])) {
        return(default)
      }
      as.numeric(x[1])
    }

    for (i in seq_along(subgroup_levels)) {

      lev <- subgroup_levels[i]
      one_level <- dat_sg[
        as.character(dat_sg[[sg_var]]) == lev, ,
        drop = FALSE
      ]

      k_level <- nrow(one_level)
      ev <- sum(one_level$Events, na.rm = TRUE)
      tot <- sum(one_level$Total, na.rm = TRUE)

      if (k_level == 1) {

        exact <- exact_binomial(
          events = one_level$Events[1],
          total = one_level$Total[1]
        )

        detail_list[[i]] <- data.frame(
          Sheet_Name = as.character(sheet_name),
          Subgroup_Var = as.character(sg_var),
          Subgroup_Level = as.character(lev),
          Studies = 1L,
          Events = ev,
          Total = tot,
          Proportion = unname(exact["proportion"]),
          CI_lower = unname(exact["lower"]),
          CI_upper = unname(exact["upper"]),
          I2 = NA_real_,
          PI_lower = NA_real_,
          PI_upper = NA_real_,
          Estimation_Method =
            "Exact binomial (Clopper-Pearson); no pooling",
          Status = "Descriptive estimate; one study only",
          stringsAsFactors = FALSE
        )

        cat(lev, ": k = 1; exact binomial 95% CI calculated.\n")

      } else {

        # Extract the same subgroup model used by the forest plot.
        j <- match(lev, as.character(meta.sg$subgroup.levels))
        stopifnot(!is.na(j))
        n_h <- meta.sg$n.harmonic.mean.w[j]
        if (!length(n_h) || !is.finite(n_h))
          n_h <- k_level / sum(1 / one_level$Total)
        bt <- function(field) {
          z <- meta.sg[[field]][j]
          if (!length(z) || is.na(z)) return(NA_real_)
          meta:::backtransf(z, sm = "PFT", n = n_h)
        }
        detail_list[[i]] <- data.frame(
          Sheet_Name = sheet_name, Subgroup_Var = sg_var,
          Subgroup_Level = lev, Studies = k_level, Events = ev, Total = tot,
          Proportion = bt("TE.random.w"), CI_lower = bt("lower.random.w"),
          CI_upper = bt("upper.random.w"), I2 = get_scalar(meta.sg$I2.w[j]),
          PI_lower = bt("lower.predict.w"), PI_upper = bt("upper.predict.w"),
          Estimation_Method = "PFT random-effects meta-analysis (REML)",
          Status = "Pooled estimate", stringsAsFactors = FALSE
        )
      }
    }

    subgroup_details <- dplyr::bind_rows(detail_list)

    # Overall pooled estimate (all eligible studies in this worksheet
    # and subgroup variable). The original combined meta model is retained.
    overall_row <- data.frame(
      Sheet_Name = as.character(sheet_name),
      Subgroup_Var = as.character(sg_var),
      Subgroup_Level = "Overall",
      Studies = get_scalar(meta.sg$k),
      Events = sum(meta.sg$event, na.rm = TRUE),
      Total = sum(meta.sg$n, na.rm = TRUE),
      Proportion = raw_pft(meta.sg$TE.random, meta.sg$n),
      CI_lower = raw_pft(meta.sg$lower.random, meta.sg$n),
      CI_upper = raw_pft(meta.sg$upper.random, meta.sg$n),
      I2 = get_scalar(meta.sg$I2),
      PI_lower = raw_pft(meta.sg$lower.predict, meta.sg$n),
      PI_upper = raw_pft(meta.sg$upper.predict, meta.sg$n),
      Estimation_Method = "PFT random-effects meta-analysis (REML)",
      Status = "Overall pooled estimate",
      stringsAsFactors = FALSE
    )

    subgroup_details <- dplyr::bind_rows(
      subgroup_details,
      overall_row
    )

    # 4. Subgroup-difference tests
    subgroup_test <- data.frame(
      Sheet_Name = sheet_name,
      Subgroup_Var = sg_var,
      Test = c(
        "Test for subgroup differences (Common effect)",
        "Test for subgroup differences (Random effects)"
      ),
      Chi_squared = c(meta.sg$Q.b.common, meta.sg$Q.b.random),
      df = c(meta.sg$df.Q.b.common, meta.sg$df.Q.b.random),
      P_value = c(meta.sg$pval.Q.b.common, meta.sg$pval.Q.b.random),
      Status = "Success",
      stringsAsFactors = FALSE
    )

    # 5. Forest plot
    # Retain the original meta::forest() RevMan layout and styling.
    # method.ci = "CP" requests Clopper-Pearson exact CIs for study rows.
    # The k = 1 subgroup is descriptive only; its pooled subgroup summary
    # is not interpreted as a random-effects estimate.
    safe_sheet <- gsub("[\\\\/:*?\"<>|]", "_", sheet_name)
    safe_var <- gsub("[\\\\/:*?\"<>|]", "_", sg_var)
    pdf_file <- file.path(
      output_dir, paste0(safe_sheet, "_", safe_var, "_forest.pdf")
    )
    pdf_open <- FALSE
    plot_success <- FALSE

    tryCatch({
      grDevices::cairo_pdf(
        filename = pdf_file,
        height = max(10, 0.4 * meta.sg$k + 6),
        width = 18,
        family = pdf_font
      )
      pdf_open <- TRUE

      meta.cp <- make_cp_copy(meta.sg)
      forest_exact(
        meta.cp,
        layout = "Revman",
        leftcols = c("studlab", "event", "n", "effect.ci"),
        rightcols = FALSE,
        family = pdf_font,
        fontsize = 10,
        col.vertical.line = "transparent",
        col.diamond.random = "#669999",
        col.diamond.lines.random = "#669999",
        col.square = "#336699",
        col.square.lines = "#336699",
        col.study = "#006666",
        col.inside = "black",
        col.random = "#660099",
        col.subgroup = "#003399",
        overall = TRUE,
        overall.hetstat = TRUE,
        common = FALSE,
        random = TRUE,
        prediction = TRUE,
        backtransf = TRUE,
        weight.study = "same",
        squaresize = 0.8,
        height = 0.6,
        plotwidth = "12cm",
        colgap.forest.left = "2cm",
        colgap.left = "2cm",
        xlab = "",
        xlim = c(0, 1),
        smlab = "Proportion",
        digits = 2,
        digits.pval = 3,
        digits.pval.Q = 3
      )

      grDevices::dev.off()
      pdf_open <- FALSE
      plot_success <- TRUE
      cat("Forest plot saved:", pdf_file, "\n")

    }, error = function(e) {
      if (pdf_open && grDevices::dev.cur() > 1) {
        try(grDevices::dev.off(), silent = TRUE)
      }
      if (file.exists(pdf_file)) unlink(pdf_file)
      cat("Forest plot error:", conditionMessage(e), "\n")
    })

    sheet_results[[sg_var]] <- list(
      meta_object = meta.sg,
      subgroup_test = subgroup_test,
      subgroup_details = subgroup_details,
      cleaned_data = dat_sg,
      plot_success = plot_success,
      study_ci = data.frame(Sheet_Name = sheet_name, Subgroup_Var = sg_var,
        Study = meta.sg$studlab, Subgroup = meta.sg$subgroup, Events = meta.sg$event,
        Total = meta.sg$n, Proportion = meta.sg$event/meta.sg$n,
        CI_lower = make_cp_copy(meta.sg)$lower, CI_upper = make_cp_copy(meta.sg)$upper),
      analysis_success = TRUE
    )

    cat("Subgroup analysis completed:", sg_var, "; PDF success:", plot_success, "\n")
  }

  all_results[[sheet_name]] <- sheet_results
}

# 6. Export all results
if (!length(all_results)) {
  stop("No worksheet could be analyzed. Check the input file and column names.")
}

all_tests <- combine_tables(all_results, "subgroup_test")
all_details <- combine_tables(all_results, "subgroup_details")

wb <- createWorkbook()
all_study_ci <- combine_tables(all_results, "study_ci")
addWorksheet(wb, "Study_Exact_CI")
writeData(wb, "Study_Exact_CI", all_study_ci)
plot_status <- bind_rows(lapply(names(all_results), function(sh) {
  bind_rows(lapply(names(all_results[[sh]]), function(sg) {
    data.frame(Sheet_Name = sh, Subgroup_Var = sg,
               PDF_success = all_results[[sh]][[sg]]$plot_success)
  }))
}))
addWorksheet(wb, "PDF_Status")
writeData(wb, "PDF_Status", plot_status)

addWorksheet(wb, "Subgroup_Tests")
if (nrow(all_tests)) writeData(wb, "Subgroup_Tests", all_tests)

addWorksheet(wb, "Subgroup_Results")
if (nrow(all_details)) writeData(wb, "Subgroup_Results", all_details)

# Export one cleaned-data sheet per source worksheet.
used_names <- c("Subgroup_Tests", "Subgroup_Results")

for (sheet_name in names(all_results)) {
  sheet_res <- all_results[[sheet_name]]
  clean_data <- NULL

  for (sg_var in names(sheet_res)) {
    if (!is.null(sheet_res[[sg_var]]$cleaned_data)) {
      clean_data <- sheet_res[[sg_var]]$cleaned_data
      break
    }
  }
  if (is.null(clean_data)) next

  base <- paste0("Clean_", substr(
    gsub("[\\\\/:*?\\[\\]]", "_", sheet_name), 1, 24
  ))
  clean_name <- base
  suffix <- 1
  while (clean_name %in% used_names) {
    suffix <- suffix + 1
    clean_name <- paste0(substr(base, 1, 27), "_", suffix)
  }
  clean_name <- substr(clean_name, 1, 31)

  addWorksheet(wb, clean_name)
  writeData(wb, clean_name, clean_data)
  used_names <- c(used_names, clean_name)
}

# Methods documentation
addWorksheet(wb, "Methods_Readme")
methods_df <- data.frame(
  Item = c(
    "Outcome",
    "Single-study estimate",
    "Single-study 95% CI",
    "Multi-study summary measure",
    "Random-effects model",
    "Tau-squared estimator",
    "Subgroup variables",
    "Minimum Total filter",
    "Continuity correction for exact single-study estimates",
    "R version",
    "meta package version"
  ),
  Setting = c(
    "Single proportion",
    "Events / Total",
    "Clopper-Pearson exact binomial method for every individual study",
    "Freeman-Tukey double-arcsine transformation (PFT) when k >= 2",
    "Random effects",
    "REML",
    paste(subgroup_vars, collapse = "; "),
    paste0("Total >= ", minimum_total),
    "None",
    R.version.string,
    as.character(packageVersion("meta"))
  ),
  stringsAsFactors = FALSE
)
writeData(wb, "Methods_Readme", methods_df)

results_file <- file.path(
  output_dir, "All_Sheets_Multi_Subgroup_Analysis.xlsx"
)
saveWorkbook(wb, results_file, overwrite = TRUE)

cat("\nAnalysis completed.\n")
cat("Excel results:", results_file, "\n")
cat("Forest plots:", output_dir, "\n")
print(plot_status)
if (any(!plot_status$PDF_success)) warning("Some PDF plots failed; see PDF_Status.")

# Print significant random-effects subgroup tests
if (nrow(all_tests) > 0) {
  sig <- all_tests %>%
    filter(
      grepl("Random", Test),
      !is.na(P_value),
      P_value < 0.05
    )

  cat("\nSignificant subgroup-difference tests (P < 0.05):\n")
  if (nrow(sig) > 0) {
    print(sig[, c("Sheet_Name", "Subgroup_Var", "P_value")])
  } else {
    cat("No statistically significant subgroup differences.\n")
  }
}

# Print all k = 1 subgroup estimates for verification
if (nrow(all_details) > 0) {
  single_results <- all_details %>%
    filter(Studies == 1, Subgroup_Level != "Overall") %>%
    select(
      Sheet_Name, Subgroup_Var, Subgroup_Level,
      Events, Total, Proportion, CI_lower, CI_upper,
      I2, Estimation_Method
    )

  cat("\nSingle-study subgroups (exact binomial CIs):\n")
  if (nrow(single_results) > 0) {
    print(single_results)
  } else {
    cat("No single-study subgroups found.\n")
  }
}
