*============================================================*
* Subgroup meta-analysis of antimicrobial resistance
* Stata official meta suite
* Freeman-Tukey transformation + Random effects + REML
*============================================================*

* 1. Declare the meta-analysis of proportions
* events = resistant isolates
* total  = total isolates
meta esize events total, esize(ftukeyprop) studylabel(study)

* The output should indicate:
* Effect size: Freeman-Tukey's p
* Model: Random effects
* Method: REML


*============================================================*
* Subgroup analysis by AST guideline
*============================================================*

* Numerical results
meta summarize, subgroup(astguideline) proportion

* Forest plot
meta forestplot, subgroup(astguideline) proportion


*============================================================*
* Other subgroup analyses: replace the variable in subgroup()
*============================================================*

* WHO Region
meta summarize, subgroup(region) proportion
meta forestplot, subgroup(region) proportion

* Income level
meta summarize, subgroup(incomelevel) proportion
meta forestplot, subgroup(incomelevel) proportion

* Year group
meta summarize, subgroup(yeargroup) proportion
meta forestplot, subgroup(yeargroup) proportion

* Etiology
meta summarize, subgroup(etiology) proportion
meta forestplot, subgroup(etiology) proportion

* AST method
meta summarize, subgroup(astmethod) proportion
meta forestplot, subgroup(astmethod) proportion