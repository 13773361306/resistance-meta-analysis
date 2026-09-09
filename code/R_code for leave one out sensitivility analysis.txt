# ============================================================
# SINGLE-PROPORTION META-ANALYSIS - LEAVE-ONE-OUT SENSITIVITY
# PFT + Random-effects + REML
# ============================================================

library(meta)
library(readxl)
library(openxlsx)
library(dplyr)
library(purrr)

# ============================================================
# 1. CONFIGURATION
# ============================================================

config <- list(
  file_path <- file.path("path", "to", "your", "data.xlsx"),
  required_cols = c("Events", "Total", "study"),
  results_dir = "D:/Desktop/meta/",
  min_studies = 5
)

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

# ------------------------------------------------------------
# 2.1 PFT back-transformation
# Uses harmonic mean of sample sizes for back-transformation
# ------------------------------------------------------------

backtransform_pft <- function(x, n_vector) {

  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  n_vector <- as.numeric(n_vector)
  n_vector <- n_vector[is.finite(n_vector) & n_vector > 0]

  if (length(n_vector) == 0) {
    return(NA_real_)
  }

  harmonic_n <- 1 / mean(1 / n_vector)

  result <- tryCatch(
    {
      as.numeric(
        meta::backtransf(
          x[1],
          sm = "PFT",
          n = harmonic_n
        )
      )
    },
    error = function(e) NA_real_
  )

  return(result)
}


# ------------------------------------------------------------
# 2.2 Safely extract a single numeric value
# ------------------------------------------------------------

first_numeric <- function(x) {

  if (is.null(x) || length(x) == 0) {
    return(NA_real_)
  }

  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  return(x[1])
}


# ============================================================
# 3. DATA READING AND CLEANING
# ============================================================

read_and_clean_data <- function(file_path, sheet_name, required_cols) {

  tryCatch({

    data <- read_excel(
      file_path,
      sheet = sheet_name
    )

    # Check required columns
    missing_cols <- setdiff(required_cols, names(data))

    if (length(missing_cols) > 0) {

      warning(
        sprintf(
          "Sheet %s missing columns: %s",
          sheet_name,
          paste(missing_cols, collapse = ", ")
        )
      )

      return(NULL)
    }

    cleaned <- data %>%
      mutate(
        Events = suppressWarnings(as.numeric(Events)),
        Total = suppressWarnings(as.numeric(Total)),
        study = as.character(study)
      ) %>%
      filter(
        complete.cases(Events, Total, study),
        study != "",
        Total >= 10,
        Events >= 0,
        Events <= Total
      )

    if (nrow(cleaned) == 0) {

      warning(
        sprintf(
          "Sheet %s no valid data after cleaning",
          sheet_name
        )
      )

      return(NULL)
    }

    return(cleaned)

  }, error = function(e) {

    warning(
      sprintf(
        "Error reading Sheet %s: %s",
        sheet_name,
        e$message
      )
    )

    return(NULL)
  })
}


# ============================================================
# 4. MAIN META-ANALYSIS
# ============================================================

perform_main_meta <- function(data) {

  if (nrow(data) < 2) {
    return(NULL)
  }

  meta_result <- tryCatch({

    metaprop(
      event = Events,
      n = Total,
      studlab = study,
      data = data,

      sm = "PFT",

      common = FALSE,
      random = TRUE,

      method.tau = "REML",

      prediction = TRUE,
      level.predict = 0.95,

      backtransf = TRUE
    )

  }, error = function(e) {

    warning(
      sprintf(
        "Main meta-analysis failed: %s",
        e$message
      )
    )

    return(NULL)
  })

  return(meta_result)
}


# ============================================================
# 5. EXTRACT MAIN META RESULTS
# ============================================================

extract_main_results <- function(meta_result) {

  if (is.null(meta_result)) {
    return(NULL)
  }

  pooled <- backtransform_pft(
    meta_result$TE.random,
    meta_result$n
  )

  lower <- backtransform_pft(
    meta_result$lower.random,
    meta_result$n
  )

  upper <- backtransform_pft(
    meta_result$upper.random,
    meta_result$n
  )

  pi_lower <- backtransform_pft(
    meta_result$lower.predict,
    meta_result$n
  )

  pi_upper <- backtransform_pft(
    meta_result$upper.predict,
    meta_result$n
  )

  data.frame(
    Proportion = pooled,
    CI_Lower = lower,
    CI_Upper = upper,
    PI_Lower = pi_lower,
    PI_Upper = pi_upper,
    I2 = first_numeric(meta_result$I2),
    Tau2 = first_numeric(meta_result$tau2),
    Q = first_numeric(meta_result$Q),
    df_Q = first_numeric(meta_result$df.Q),
    P_Heterogeneity = first_numeric(meta_result$pval.Q)
  )
}


# ============================================================
# 6. LEAVE-ONE-OUT SENSITIVITY ANALYSIS
# ============================================================

perform_leave_one_out <- function(data, main_result) {

  if (nrow(data) < 3) {
    return(NULL)
  }

  # Original results
  original_prop <- backtransform_pft(
    main_result$TE.random,
    main_result$n
  )

  original_i2 <- first_numeric(main_result$I2)
  original_tau2 <- first_numeric(main_result$tau2)

  sensitivity_results <- map_dfr(
    seq_len(nrow(data)),
    function(i) {

      data_subset <- data[-i, , drop = FALSE]

      if (nrow(data_subset) < 2) {
        return(NULL)
      }

      meta_subset <- tryCatch({

        metaprop(
          event = Events,
          n = Total,
          studlab = study,
          data = data_subset,

          sm = "PFT",

          common = FALSE,
          random = TRUE,

          method.tau = "REML",

          prediction = TRUE,
          level.predict = 0.95,

          backtransf = TRUE
        )

      }, error = function(e) {

        warning(
          sprintf(
            "Meta-analysis failed after excluding %s: %s",
            data$study[i],
            e$message
          )
        )

        return(NULL)
      })

      if (is.null(meta_subset)) {
        return(NULL)
      }

      subset_prop <- backtransform_pft(
        meta_subset$TE.random,
        meta_subset$n
      )

      subset_lower <- backtransform_pft(
        meta_subset$lower.random,
        meta_subset$n
      )

      subset_upper <- backtransform_pft(
        meta_subset$upper.random,
        meta_subset$n
      )

      subset_pi_lower <- backtransform_pft(
        meta_subset$lower.predict,
        meta_subset$n
      )

      subset_pi_upper <- backtransform_pft(
        meta_subset$upper.predict,
        meta_subset$n
      )

      subset_i2 <- first_numeric(meta_subset$I2)
      subset_tau2 <- first_numeric(meta_subset$tau2)

      change_prop <- subset_prop - original_prop
      change_pp <- change_prop * 100

      data.frame(
        Excluded_Study = data$study[i],

        Proportion = subset_prop,
        CI_Lower = subset_lower,
        CI_Upper = subset_upper,

        PI_Lower = subset_pi_lower,
        PI_Upper = subset_pi_upper,

        I2 = subset_i2,
        Tau2 = subset_tau2,

        Change_Proportion = change_prop,
        Change_PP = change_pp,
        Absolute_Change_PP = abs(change_pp),

        I2_Change = subset_i2 - original_i2,
        Tau2_Change = subset_tau2 - original_tau2
      )
    }
  )

  return(sensitivity_results)
}


# ============================================================
# 7. SUMMARIZE LEAVE-ONE-OUT RESULTS
# ============================================================

summarize_leave_one_out <- function(sensitivity_results) {

  if (
    is.null(sensitivity_results) ||
    nrow(sensitivity_results) == 0
  ) {

    return(
      list(
        loo_min = NA_real_,
        loo_max = NA_real_,
        max_abs_change_pp = NA_real_,
        largest_change_study = NA_character_,
        largest_change_direction = NA_real_,
        min_i2 = NA_real_,
        min_tau2 = NA_real_,
        study_min_i2 = NA_character_,
        study_min_tau2 = NA_character_
      )
    )
  }

  valid_change <- sensitivity_results %>%
    filter(is.finite(Absolute_Change_PP))

  if (nrow(valid_change) == 0) {

    largest_change_study <- NA_character_
    largest_change_direction <- NA_real_
    max_abs_change_pp <- NA_real_

  } else {

    idx_change <- which.max(
      valid_change$Absolute_Change_PP
    )

    largest_change_study <-
      valid_change$Excluded_Study[idx_change]

    largest_change_direction <-
      valid_change$Change_PP[idx_change]

    max_abs_change_pp <-
      valid_change$Absolute_Change_PP[idx_change]
  }


  # Minimum I2
  valid_i2 <- sensitivity_results %>%
    filter(is.finite(I2))

  if (nrow(valid_i2) > 0) {

    idx_i2 <- which.min(valid_i2$I2)

    min_i2 <- valid_i2$I2[idx_i2]
    study_min_i2 <-
      valid_i2$Excluded_Study[idx_i2]

  } else {

    min_i2 <- NA_real_
    study_min_i2 <- NA_character_
  }


  # Minimum Tau2
  valid_tau2 <- sensitivity_results %>%
    filter(is.finite(Tau2))

  if (nrow(valid_tau2) > 0) {

    idx_tau2 <- which.min(valid_tau2$Tau2)

    min_tau2 <- valid_tau2$Tau2[idx_tau2]
    study_min_tau2 <-
      valid_tau2$Excluded_Study[idx_tau2]

  } else {

    min_tau2 <- NA_real_
    study_min_tau2 <- NA_character_
  }


  return(
    list(
      loo_min = min(
        sensitivity_results$Proportion,
        na.rm = TRUE
      ),

      loo_max = max(
        sensitivity_results$Proportion,
        na.rm = TRUE
      ),

      max_abs_change_pp = max_abs_change_pp,

      largest_change_study =
        largest_change_study,

      largest_change_direction =
        largest_change_direction,

      min_i2 = min_i2,

      min_tau2 = min_tau2,

      study_min_i2 = study_min_i2,

      study_min_tau2 = study_min_tau2
    )
  )
}


# ============================================================
# 8. LEAVE-ONE-OUT PLOT
# ============================================================

plot_leave_one_out <- function(
    sensitivity_results,
    sheet_name,
    results_dir) {

  if (
    is.null(sensitivity_results) ||
    nrow(sensitivity_results) == 0
  ) {

    return(NULL)
  }

  plot_data <- sensitivity_results %>%
    filter(is.finite(Change_PP)) %>%
    arrange(Change_PP)

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  safe_file_name <- gsub(
    "[\\\\/:*?\"<>|]",
    "_",
    sheet_name
  )

  pdf_file <- file.path(
    results_dir,
    paste0(
      safe_file_name,
      "_Leave_One_Out.pdf"
    )
  )

  pdf(
    pdf_file,
    width = 10,
    height = max(
      6,
      0.32 * nrow(plot_data) + 2
    )
  )

  par(
    mar = c(5, 12, 4, 2)
  )

  max_abs <- max(
    abs(plot_data$Change_PP),
    na.rm = TRUE
  )

  if (!is.finite(max_abs) || max_abs == 0) {
    max_abs <- 1
  }

  x_lim <- max_abs * 1.15

  plot(
    x = plot_data$Change_PP,
    y = seq_len(nrow(plot_data)),

    xlim = c(-x_lim, x_lim),

    ylim = c(
      0.5,
      nrow(plot_data) + 0.5
    ),

    xlab =
      "Change in pooled resistance proportion (percentage points)",

    ylab = "",

    yaxt = "n",

    pch = 19,

    main =
      paste(
        "Leave-One-Out Sensitivity Analysis:",
        sheet_name
      )
  )

  axis(
    side = 2,
    at = seq_len(nrow(plot_data)),
    labels = plot_data$Excluded_Study,
    las = 2,
    cex.axis = 0.65
  )

  abline(
    v = 0,
    lty = 2,
    lwd = 1.5
  )

  dev.off()

  return(pdf_file)
}


# ============================================================
# 9. SAVE EXCEL RESULTS
# ============================================================

save_results <- function(all_results, config) {

  # ----------------------------------------------------------
  # 9.1 Summary
  # ----------------------------------------------------------

  summary_data <- map_dfr(
    names(all_results),
    function(sheet_name) {

      result <- all_results[[sheet_name]]

      main <- result$main_summary
      loo <- result$loo_summary

      data.frame(

        Sheet = sheet_name,

        Studies = nrow(result$data),

        Sensitivity_Analysis_Performed =
          !is.null(result$sensitivity_results),

        Original_Proportion =
          round(main$Proportion, 4),

        Original_Percent =
          round(main$Proportion * 100, 2),

        CI_95_Lower =
          round(main$CI_Lower, 4),

        CI_95_Upper =
          round(main$CI_Upper, 4),

        CI_95_Lower_Percent =
          round(main$CI_Lower * 100, 2),

        CI_95_Upper_Percent =
          round(main$CI_Upper * 100, 2),

        PI_95_Lower =
          round(main$PI_Lower, 4),

        PI_95_Upper =
          round(main$PI_Upper, 4),

        PI_95_Lower_Percent =
          round(main$PI_Lower * 100, 2),

        PI_95_Upper_Percent =
          round(main$PI_Upper * 100, 2),

        I2 =
          round(main$I2, 2),

        Tau2 =
          round(main$Tau2, 4),

        LOO_Min_Percent =
          round(loo$loo_min * 100, 2),

        LOO_Max_Percent =
          round(loo$loo_max * 100, 2),

        Maximum_Absolute_Change_PP =
          round(
            loo$max_abs_change_pp,
            2
          ),

        Study_Producing_Largest_Change =
          loo$largest_change_study,

        Change_After_Exclusion_PP =
          round(
            loo$largest_change_direction,
            2
          ),

        Minimum_I2_After_Exclusion =
          round(
            loo$min_i2,
            2
          ),

        Study_Producing_Minimum_I2 =
          loo$study_min_i2,

        Minimum_Tau2_After_Exclusion =
          round(
            loo$min_tau2,
            4
          ),

        Study_Producing_Minimum_Tau2 =
          loo$study_min_tau2
      )
    }
  )


  # ----------------------------------------------------------
  # 9.2 Detailed results
  # ----------------------------------------------------------

  detailed_data <- map_dfr(
    names(all_results),
    function(sheet_name) {

      result <- all_results[[sheet_name]]

      if (
        !is.null(result$sensitivity_results) &&
        nrow(result$sensitivity_results) > 0
      ) {

        result$sensitivity_results %>%
          mutate(
            Sheet = sheet_name,
            .before = 1
          )
      }
    }
  )


  # ----------------------------------------------------------
  # 9.3 Workbook
  # ----------------------------------------------------------

  wb <- createWorkbook()


  # Summary
  addWorksheet(
    wb,
    "Summary"
  )

  writeData(
    wb,
    "Summary",
    summary_data
  )

  freezePane(
    wb,
    "Summary",
    firstRow = TRUE
  )

  setColWidths(
    wb,
    "Summary",
    cols = 1:ncol(summary_data),
    widths = "auto"
  )


  # Detailed
  if (
    !is.null(detailed_data) &&
    nrow(detailed_data) > 0
  ) {

    addWorksheet(
      wb,
      "Detailed_Results"
    )

    writeData(
      wb,
      "Detailed_Results",
      detailed_data
    )

    freezePane(
      wb,
      "Detailed_Results",
      firstRow = TRUE
    )

    setColWidths(
      wb,
      "Detailed_Results",
      cols = 1:ncol(detailed_data),
      widths = "auto"
    )
  }


  # ----------------------------------------------------------
  # 9.4 Methods / Software information
  # ----------------------------------------------------------

  methods_info <- data.frame(

    Item = c(
      "Analysis",
      "Effect measure",
      "Meta-analysis model",
      "Between-study variance estimator",
      "Prediction interval",
      "Sensitivity analysis",
      "Minimum studies for leave-one-out",
      "R version",
      "meta package version"
    ),

    Value = c(
      "Single-proportion meta-analysis",
      "Freeman-Tukey double-arcsine transformation (PFT)",
      "Random-effects model",
      "Restricted maximum likelihood (REML)",
      "95% prediction interval",
      "Leave-one-out analysis by sequentially excluding each study",
      config$min_studies,
      R.version.string,
      as.character(packageVersion("meta"))
    )
  )

  addWorksheet(
    wb,
    "Methods_Readme"
  )

  writeData(
    wb,
    "Methods_Readme",
    methods_info
  )

  setColWidths(
    wb,
    "Methods_Readme",
    cols = 1:2,
    widths = "auto"
  )


  # ----------------------------------------------------------
  # 9.5 Per-sheet output
  # ----------------------------------------------------------

  for (sheet_name in names(all_results)) {

    result <- all_results[[sheet_name]]

    safe_sheet_name <- gsub(
      "[:\\\\/?*\\[\\]]",
      "_",
      sheet_name
    )

    safe_sheet_name <- substr(
      paste0(
        "LOO_",
        safe_sheet_name
      ),
      1,
      31
    )

    # Avoid duplicate Excel sheet names
    original_safe_name <- safe_sheet_name
    counter <- 1

    while (
      safe_sheet_name %in% names(wb)
    ) {

      suffix <- paste0("_", counter)

      safe_sheet_name <- paste0(
        substr(
          original_safe_name,
          1,
          31 - nchar(suffix)
        ),
        suffix
      )

      counter <- counter + 1
    }

    addWorksheet(
      wb,
      safe_sheet_name
    )

    main <- result$main_summary
    loo <- result$loo_summary

    info_df <- data.frame(

      Item = c(
        "Outcome / Sheet",
        "Number of studies",
        "Original pooled resistance",
        "95% CI",
        "95% prediction interval",
        "I2",
        "Tau2",
        "Leave-one-out minimum",
        "Leave-one-out maximum",
        "Maximum absolute change (pp)",
        "Study producing largest change",
        "Change after exclusion (pp)"
      ),

      Value = c(
        sheet_name,

        nrow(result$data),

        sprintf(
          "%.2f%%",
          main$Proportion * 100
        ),

        sprintf(
          "%.2f%% - %.2f%%",
          main$CI_Lower * 100,
          main$CI_Upper * 100
        ),

        ifelse(
          is.finite(main$PI_Lower) &
            is.finite(main$PI_Upper),

          sprintf(
            "%.2f%% - %.2f%%",
            main$PI_Lower * 100,
            main$PI_Upper * 100
          ),

          "NA"
        ),

        sprintf(
          "%.2f%%",
          main$I2
        ),

        sprintf(
          "%.4f",
          main$Tau2
        ),

        ifelse(
          is.finite(loo$loo_min),

          sprintf(
            "%.2f%%",
            loo$loo_min * 100
          ),

          "NA"
        ),

        ifelse(
          is.finite(loo$loo_max),

          sprintf(
            "%.2f%%",
            loo$loo_max * 100
          ),

          "NA"
        ),

        ifelse(
          is.finite(
            loo$max_abs_change_pp
          ),

          sprintf(
            "%.2f",
            loo$max_abs_change_pp
          ),

          "NA"
        ),

        ifelse(
          is.na(
            loo$largest_change_study
          ),

          "NA",

          loo$largest_change_study
        ),

        ifelse(
          is.finite(
            loo$largest_change_direction
          ),

          sprintf(
            "%.2f",
            loo$largest_change_direction
          ),

          "NA"
        )
      )
    )

    writeData(
      wb,
      safe_sheet_name,
      info_df,
      startRow = 1
    )

    if (
      !is.null(result$sensitivity_results) &&
      nrow(result$sensitivity_results) > 0
    ) {

      writeData(
        wb,
        safe_sheet_name,
        result$sensitivity_results,
        startRow = 16
      )
    }

    setColWidths(
      wb,
      safe_sheet_name,
      cols = 1:15,
      widths = "auto"
    )
  }


  # ----------------------------------------------------------
  # 9.6 Save Excel
  # ----------------------------------------------------------

  timestamp <- format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )

  output_file <- file.path(
    config$results_dir,
    paste0(
      "Leave_One_Out_Sensitivity_Final_",
      timestamp,
      ".xlsx"
    )
  )

  saveWorkbook(
    wb,
    output_file,
    overwrite = TRUE
  )

  return(output_file)
}


# ============================================================
# 10. MAIN PROGRAM
# ============================================================

run_sensitivity_analysis <- function(config) {

  cat("\n")
  cat("============================================================\n")
  cat("   Leave-One-Out Sensitivity Analysis: PFT + REML\n")
  cat("============================================================\n\n")

  sheet_names <- excel_sheets(
    config$file_path
  )

  cat(
    sprintf(
      "Found %d sheets.\n\n",
      length(sheet_names)
    )
  )

  all_results <- list()


  for (sheet_name in sheet_names) {

    cat(
      sprintf(
        ">>> Processing: %s\n",
        sheet_name
      )
    )


    # --------------------------------------------------------
    # Data
    # --------------------------------------------------------

    data <- read_and_clean_data(
      config$file_path,
      sheet_name,
      config$required_cols
    )

    if (is.null(data)) {
      next
    }

    cat(
      sprintf(
        "    Valid studies: %d\n",
        nrow(data)
      )
    )


    # --------------------------------------------------------
    # Main meta-analysis
    # --------------------------------------------------------

    main_result <- perform_main_meta(
      data
    )

    if (is.null(main_result)) {

      cat(
        "    Main meta-analysis failed. Skipping.\n\n"
      )

      next
    }


    main_summary <- extract_main_results(
      main_result
    )


    cat(
      sprintf(
        paste0(
          "    Pooled resistance: %.2f%% ",
          "(95%% CI %.2f%% - %.2f%%)\n"
        ),

        main_summary$Proportion * 100,

        main_summary$CI_Lower * 100,

        main_summary$CI_Upper * 100
      )
    )


    cat(
      sprintf(
        "    I² = %.1f%%; Tau² = %.4f\n",
        main_summary$I2,
        main_summary$Tau2
      )
    )


    # --------------------------------------------------------
    # Leave-one-out
    # --------------------------------------------------------

    if (
      nrow(data) >=
      config$min_studies
    ) {

      sensitivity_results <-
        perform_leave_one_out(
          data,
          main_result
        )

      loo_summary <-
        summarize_leave_one_out(
          sensitivity_results
        )

      cat(
        sprintf(
          paste0(
            "    Leave-one-out range: ",
            "%.2f%% - %.2f%%\n"
          ),

          loo_summary$loo_min * 100,

          loo_summary$loo_max * 100
        )
      )

      cat(
        sprintf(
          paste0(
            "    Maximum absolute change: ",
            "%.2f percentage points\n"
          ),

          loo_summary$max_abs_change_pp
        )
      )

      cat(
        sprintf(
          paste0(
            "    Study producing largest change: ",
            "%s\n"
          ),

          loo_summary$largest_change_study
        )
      )


      # Plot
      plot_leave_one_out(
        sensitivity_results,
        sheet_name,
        config$results_dir
      )

    } else {

      sensitivity_results <- NULL

      loo_summary <-
        summarize_leave_one_out(
          NULL
        )

      cat(
        sprintf(
          paste0(
            "    Leave-one-out not performed: ",
            "studies %d < %d\n"
          ),

          nrow(data),

          config$min_studies
        )
      )
    }


    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    all_results[[sheet_name]] <- list(

      data = data,

      main_result =
        main_result,

      main_summary =
        main_summary,

      sensitivity_results =
        sensitivity_results,

      loo_summary =
        loo_summary
    )

    cat("\n")
  }


  # ----------------------------------------------------------
  # Output
  # ----------------------------------------------------------

  if (length(all_results) == 0) {

    cat(
      "No sheets were successfully analyzed.\n"
    )

    return(NULL)
  }


  output_file <-
    save_results(
      all_results,
      config
    )


  cat("\n")
  cat("============================================================\n")
  cat("Analysis complete\n")
  cat("============================================================\n")

  cat(
    sprintf(
      "Results file: %s\n",
      output_file
    )
  )

  cat(
    sprintf(
      "Sheets successfully processed: %d\n",
      length(all_results)
    )
  )

  invisible(
    all_results
  )
}


# ============================================================
# 11. EXECUTE
# ============================================================

results <- run_sensitivity_analysis(
  config
)