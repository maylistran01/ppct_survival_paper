## Overview

This repository contains the code supporting the paper:

**Prediction-Powered Inference in Clinical Trials with Survival Outcomes**  
*Under review*

Recruiting patients for clinical trials in rare neurodegenerative diseases is challenging due to ethical concerns surrounding placebo enrollment and difficulties in achieving sufficient statistical power. Prediction-Powered Inference for Clinical Trials (PPCT) provides a principled way to improve efficiency by leveraging external predictive models.

In this work, we extend PPCT to **time-to-event outcomes** by focusing on the **Restricted Mean Survival Time (RMST)**, a survival summary measure that:
- Remains within a linear framework
- Allows direct comparison of average survival via differences in means

We compute individual-level RMST pseudo-values and use external models trained on observational data to predict counterfactual RMST pseudo-values. These predictions are then incorporated into the PPCT framework to strengthen statistical inference.

The method is evaluated through simulation studies and applied to the Trophos ALS clinical trial, where we show that 49 out of 511 patients could be removed while maintaining equivalent statistical power.

---

## Repository Structure

```text
.
├── application_trophos/
│   ├── 0_extract_clean_data.ipynb
│   ├── 1_leaspy_ready.ipynb
│   ├── 2_cohorts_baseline_obs.ipynb
│   ├── 3_predict_rmst_all.ipynb
│   ├── 4_PPCT_computation_all.ipynb
│   ├── 5_estimators_comparison.ipynb
│   ├── estimators.py
|   ├── outputs/
│   └── utils_process.py
│
├── simulations/
│   ├── 0_simulate_ct_df_v2.R
│   ├── 1_simulate_RMST_preds.ipynb
│   ├── 2_simulations_analysis_ADEMP.ipynb
│   └── sim_outputs_1000/
│
└── README.Rmd
```

---

## Details on the content of this repository


Simulation experiments are located in the `simulations/` directory and follow an ADEMP-style workflow.
- `0_simulate_ct_df_v2.R`
    - Generates simulated clinical trial datasets under different conditions.
    - Simulation outputs will be stored in: `simulations/sim_ct_df/`
- `1_simulate_RMST_preds.ipynb`
    - Simulation outputs will be stored in: `simulations/sim_outputs_1000/`
- `2_simulations_analysis_ADEMP.ipynb`
Analyzes simulation outputs and produces:
    - Tables 1, 2 and 3
    - Figures 1 and 2


All analyses related to the real-data application are located in the `application_trophos/` directory.
- `0_extract_clean_data.ipynb`
    - **Important notes:**
        - Access to all ALS cohort datasets is required
        - Access requests must follow the Data Availability Statement in the paper
        - Access to the private [leaspy_converters](https://github.com/aramis-lab/leaspy_converters) repository is required
        - Some dependencies require specific versions (notably, pandas must be downgraded)
- `1_leaspy_ready.ipynb`
- `2_cohorts_baseline_obs.ipynb`
    - Computes baseline descriptive statistics for the ALS cohorts.
    - Corresponds to Table 4 in the paper.
- `3_predict_rmst_all.ipynb`
- `4_PPCT_computation_all.ipynb`
    - Implements the prediction-powered inference procedure using RMST pseudo-values.
    - Corresponds to Table 5 in the paper.
- `5_estimators_comparison.ipynb`
    - Compares standard estimators with prediction-powered estimators in terms of efficiency and variance reduction.
    - Corresponds to Table 6 in the paper.



