# flu_cog_HRS

This repo contains all the code necessary to replicate data construction and analyses in the paper: 

Wu Y, Ikesu R, Rojas-Saunero LP, Seamans MJ, Elser H, Pescador Jimenez M, Nianogo R, Gross AL, Ritz BR, Seamans MJ, Elser H, Mayeda ER. Sustained Influenza Vaccination and Memory Decline Among Middle-Aged and Older Adults in the United States: A Target Trial Emulation Study.

## Table of contents
1. [Project Descriptions](#proj_desc)
2. [Scripts Descriptions](#script_desc)

## Project Descriptions <a name="proj_desc"></a>
While prior research suggests influenza vaccination may decrease dementia risk, the existing evidence remains limited with prior studies focusing on one-time-point influenza vaccination, short follow-up periods, and uncontrolled confounding. We emulated a target trial in the U.S. Health and Retirement Study (2000-2020) to estimate effects of sustained influenza vaccination (receiving ≥1 vaccine every 4 years over 20 years) compared to the natural course (vaccination as it occurred in the population without intervention) on 22-year memory change. Using targeted maximum likelihood estimation adjusting for baseline and time-varying confounders, sustained vaccination was associated with slower memory decline.

* Dataset: Health and Retirement Study (HRS), a nationally representative longitudinal study of US adults aged 50 and older. We used data from waves 5 to 16 (2000-2022).
* Exposure: Sustained influenza vaccination.
* Outcome: Composite memory scores, combining direct and proxy memory assessments, from 2000 to 2022 
* Method: Targeted maximum likelihood estimation (TMLE) with time-varying exposures and confounders.


## Scripts Descriptions <a name="script_desc"></a>

We conducted three analyses: the main analysis and two sensitivity analyses. The main analysis included participants aged 52-80 years, the first sensitivity analysis included participants aged 52-90 years, and the second sensitivity analysis included only one participant per household (the participant with the smallest person id number from HRS) aged 52-80 years. The following table summarizes the manuscript analyses and their corresponding script labels.

*R scripts for data construction, TMLE, cluster computing submission and generating tables and figures were in the `scripts` folder. The R scripts are named in the order they should be run. Functions created by the authors were sourced where needed.*

* **1.dataset_construction.R**: Cleans the data and creates the necessary variable for the analysis.
* **2. multiple_imputation.R**: Performs multiple imputation to handle missing data using multiple imputation (predictive mean matching).
* **3.0_data_prep_cluster_SL.R**: Creates the final analytic dataset and other objects after multiple imputation for main and sensitivity analysis.
* **3.1_cluster_TMLE_superlearner.R**: Run on UCLA Hoffman cluster. Performs TMLE with SuperLearner to get effect estimates. Scenarios include sustained and natural course (as observed).
* **3.1.2_parallel_model_mispecification_TMLE_superlearner_nc_noipcw.R**: Run locally. Model specification diagnosis: comparing natural course without applying censoring weights with the observed memory trajectory.
* **3.2_results_TMLE_superlearner_2000.R**: Combines results from all imputations and calculates estimates and confidence intervals. Create figures and tables for modelling results.
* **3.3_RD_RR_calculation_2000.R**: Calculate RD and RR based on classified memory impairment and the potential outcome estimates from TMLE.
* **ager5le80_tmle_superlearner_submission.sh**: Cluster submission script for main analysis
* **Scripts 4.1-4.3**: follow the same structure as scripts 3.1, 3.2, 3.3 expect that the analysis was run on a restricted dataset with participants age < 90.
* **ager5l90_tmle_superlearner_submission.sh**: Cluster submission script for sensitivity analysis 1 (age < 90).
* **Scripts 5.1-5.3**: follow the same structure as scripts 3.1, 3.2, 3.3 expect that the analysis was run on a restricted dataset with only one person (with the smallest person id number from HRS) from the same household.
* **sens_agele80minpn_tmle_superlearner_submission.sh**: Cluster submission script for sensitivity analysis 2 (one person per household).
* **9_descriptive_tables_figures_2000.R**: Creates all descriptive tables and figures for the manuscript.
* **functions.R**: Contains all the custom helper functions used in data construction.
