# Replication Code for NEISS Legal Intervention Injuries Study 

This repository contains replication code and data for the analysis presented in **"Racial/Ethnic Inequalities for Legal Intervention Injuries Treated in U.S. Emergency Departments: United States (2004-2021)"**. The study investigates trends in legal intervention injury rates by race and ethnicity using data from the National Electronic Injury Surveillance Sytem (NEISS), incorporating race imputation and uncertainty propagation using parametric bootstrapping. 

## Overview 
The central script in this repository, `bootstrap_parametric.R`, performs the full analysis pipeline: 
1. **Imputation of race/ethnicity** <br>
   Uses a multinomial logistic model to impute race/ethnicity for injury cases with missing data. 
2. **Parametric Bootstrapping** <br>
   Propagates uncertainty from both the imputation model and the NEISS survey design using simulation-based inference. 
3. **Trend Estimation** <br>
   Computes annual injury rates by race/ethnicity and estimates temporal trends, accounting for combined sources of uncertainty.

## Repository Contents 
- `bootstrap_parametric.R` <br>
  Self-contained R script that executes all analysis steps described above.
- `neiss_legal_race_imputed.csv`<br>
  Incident-level legal intervention injury data including model-derived predicted probabilities for race/ethnicity.
- `postcensus_pop.csv`<br>
  Population estimates by race/ethnicity and year (2004-2021) obtained from the Postcensus data.

## Notes 
This repository does not include the full incident-level NEISS dataset. The original data can be downloaded from the ICPSR repository (https://www.icpsr.umich.edu/web/ICPSR/series/198). We additionally created unique hospital identifiers based on the table of date changes to each Primary Samplint Unit in NEISS, available from the Consumer Product Safety Commission upon request. 

## Citation 
If you use this code or data, please cite: https://www.medrxiv.org/content/10.1101/2025.05.09.25327323v1 
