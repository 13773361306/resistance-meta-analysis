# AMR Meta-Analysis Pipeline

A comprehensive R pipeline for single‑proportion meta‑analysis of antimicrobial resistance (AMR) data, with support for multiple statistical models, subgroup analyses, meta‑regression, publication bias diagnostics, and sensitivity analyses.

## Features

- **Meta‑analysis models**
  - Freeman‑Tukey double arcsine transformation (PFT) with random‑effects REML
- **Subgroup analysis**
  - Any categorical variable (AST guideline, region, income level, year group, etc.)
  - Test for subgroup differences (random‑effects Q‑test)
  - Small‑subgroup warnings (k < 5 studies)
- **Meta‑regression**
  - Temporal trend analysis using `metafor` (REML, Z‑test)
  - Logit transformation with 0.5 continuity correction
  - Reports beta coefficient and odds ratio per year
- **Publication bias assessment**
  - Egger's test
  - Funnel plots (PDF and high‑resolution TIFF)
  - Summary diagnostic plots
- **Sensitivity analyses**
  - Leave‑one‑out analysis
  - Exclusion of high‑risk‑of‑bias studies
  - Alternative proportion model (GLMM)
- **Publication‑ready output**
  - Forest plots with 95% prediction intervals
  - Funnel plots
  - Excel workbooks with multiple sheets
  - Full reproducibility metadata (R version, package versions)

## Repository Structure

```
.
├── code/
│   ├── R_forest_PFT.R                    # Forest plots (PFT + REML)
│   ├── R_forest_GLMM.R                   # Forest plots (GLMM)
│   ├── R_subgroup_analysis.R             # Subgroup meta-analysis
│   ├── R_meta_regression.R               # Temporal meta-regression
│   ├── R_leave_one_out.R                 # Leave-one-out sensitivity
│   ├── R_funnel_egger.R                  # Funnel plots + Egger's test
│   ├── STATA_meta_regression.do          # STATA: meta-regression
│   └── STATA_subgroup_analysis.do        # STATA: subgroup analysis
├── data/
│   └── example_data.xlsx                 # Example dataset (replace with your own)
├── README.md
└── LICENSE
```

## Requirements

- R ≥ 4.4.2
- Required R packages:
  - `meta` (≥ 8.2.0)
  - `metafor` (≥ 4.8.0)
  - `readxl`, `openxlsx`, `writexl`, `dplyr`, `purrr`, `ggplot2`, `grid`

Install all dependencies with:

```r
install.packages(c("meta", "metafor", "readxl", "openxlsx", "writexl",
                   "dplyr", "purrr", "ggplot2", "grid"))
```

## Quick Start

1. **Prepare your data**
   - Save your Excel file with at least these columns: `Events`, `Total`, `study`
   - Optionally add columns for subgroup analysis (`Region`, `AST Guideline`, `Year`, etc.)

2. **Run a basic meta‑analysis (PFT)**
   ```r
   source("R/main_meta_analysis_PFT.R")
   # Configure file_path in the script, then run
   ```

3. **Run a subgroup analysis**
   ```r
   source("R/subgroup_analysis.R")
   # Change `subgroup_var` in the CONFIGURATION section to your variable
   ```

4. **Check publication bias**
   ```r
   source("R/egger_test_funnel.R")
   ```

5. **All outputs are saved to `outputs/`**:
   - PDF forest plots
   - TIFF funnel plots (300 dpi)
   - Excel workbooks with summary tables and detailed results
   - Methods_Readme sheet with full reproducibility metadata

## Configuration

Each script includes a `CONFIGURATION` section at the top. Simply modify the following parameters:

```r
config <- list(
  file_path   = "path/to/your/data.xlsx",
  results_dir = "path/to/output/folder",
  subgroup_var = "AST Guideline",   # change to "Region", "Year", etc.
  min_total   = 10                  # minimum sample size per study
)
```

## Key Methodological Settings

| Parameter | Setting |
|-----------|---------|
| Effect measure | PFT (Freeman‑Tukey) or logit (GLMM) |
| Between‑study variance | REML (PFT) / ML (GLMM) |
| Confidence intervals | Hartung‑Knapp or classic normal approximation |
| Subgroup tau² | Independent across subgroups (`tau.common = FALSE`) |
| Continuity correction | 0.5 (when using logit) |

## Output Files

- **Excel summaries**: `FINAL_*_Subgroup_Meta_Analysis.xlsx` containing:
  - `Subgroup_Results`: pooled proportions by subgroup
  - `Subgroup_Tests`: test for subgroup differences
  - `Methods_Readme`: full methodological record
  - `Clean_*`: cleaned study‑level data
- **Forest plots**: `*_Subgroup_Forest.pdf` (RevMan style, 19pt Times font)
- **Funnel plots**: `*_PFT_funnel.pdf` and `*.tiff`

## Citation

If you use this pipeline in your research, please cite this repository:

> *AMR Meta‑Analysis Pipeline. https://github.com/yourusername/amr-meta-analysis-r*

## License

MIT License — see LICENSE file for details.