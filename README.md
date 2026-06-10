# pg-ap-analytics
 
This repository contains code used to generate the results in the following preprint:
* Alireza Mahboub-Ahari, Joe Hilton, Maria Rodrigues, John McDermott, Richard Body, William G Newman, and Katherine Payne (2026) &quot;Implementing a CYP2C19-guided approach for prescribing dual antiplatelet therapy in acute coronary syndrome for patients undergoing percutaneous coronary intervention: a cost-effectiveness analysis.&quot; <i>Research Square</i> [https://doi.org/10.21203/rs.3.rs-9050351/v1]

The results in our preprint can be replicated as follows:
* The base-case deterministic results for STEMI and UA/NSTEMI in Table 2 are calculated in `deteministic_STEMI_analysis.R` and `deterministic_NSTEMI_analysis.R` respectively, and tabulated in `format_outputs.R`.
* The results on number of clinical events each year for 1,000 patients for each DAPT strategy in STEMI and UA/NSTEMI patients in Table 3 are calculated in `deteministic_STEMI_analysis.R` and `deterministic_NSTEMI_analysis.R` respectively, and tabulated in `format_outputs.R`.
* The probabilistic sensitivity analyses for CYP2C19-guided treatment versus current DAPT in STEMI and UA/NSTEMI in Table 4 are conducted in `STEMI_psa.R` and `NSTEMI_psa.R` and tabulated in `format_outputs.R`. These analyses are the source of the credible intervals for base-case results throughout the manuscript and supplementary appendices. The plots based on these results presented in supplementary Figures S5.1 and S5.2 are generated in `plot_psa_results.R`. These scripts are also the sources of the undiscounted results presented in supplementary Table S6.1.
* The deterministic sensitivity analyses covered in Supplementary Appendix 4 are conducted in `stemi_one_way_sensitivity.R` and `nstemi_one_way_sensitivity.R`. The tornado plots in this appendix are generated in `plot_dsa_results.R`.
* The scenario analyses covered in Supplementary Appendix 7 are conducted in `stemi_scenario_psa.R` and `nstemi_scenario_psa.R` and tabulated in `format_outputs.R`.
* Results quoted in the manuscript text on cohort-level incremental differences in expected numbers of events are tabulated in `generate_additional_outputs.R`.
