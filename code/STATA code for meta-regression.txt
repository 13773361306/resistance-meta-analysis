* ============================================================
* SINGLE-PROPORTION META-REGRESSION (STATA)
* ============================================================

* 1. Data cleaning (keep studies with events=total)
drop if missing(events, total, year)

* 2. 0.5 correction (consistent with R)
gen events_adj = events + 0.5
gen total_adj = total + 1

* 3. Calculate logit-transformed effect sizes
gen logit_es = log(events_adj / (total_adj - events_adj))
gen logit_se = sqrt(1/events_adj + 1/(total_adj - events_adj))

* 4. Declare meta dataset
meta set logit_es logit_se, studylabel(study)

* 5. Run meta-regression
meta regress year

* 6. Calculate and display OR (odds ratio per 1-year increase)
lincom year, eform