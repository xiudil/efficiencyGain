# efficiencyGain
estimating the efficiency gain from covariate-adjusted analysis

the R code here implements the methods in the manuscript "Estimating the Efficiency Gain of Covariate-Adjusted Analyses in Future Clinical Trials Using External Data"

"ordinal_analytical.R" implements the analytical approach to estimate the relative efficiency when the outcome in the future trial is ordinal. "helper_fns.R" contains useful functions in fitting proportional odds working model and generating data.

"ordinal_bootstrap.R" implements the double bootstrap procedure for the relative efficiency of working-model-based adjusted estimator.

"survival_analytical.R" implements the analytical approach to estimate the relative efficiency when the outcome in the future trial is time-to-event.

The folder "covid_analysis" contains the dataset, as well as code to analyze this dataset. In particular, "estimate_fns.R" contains all relevant functions in the previously mentioned R scripts for easy import. "Covid_Analysis.R" contains the preprocessing we have conducted to prepare the dataset and code to estimate the efficiency gain. "Covid_Analysis_MonteCarlo.R" implements the simulation based approach, and estimates the treatment effect using different estimators on a single simulated dataset. To estimate the relative efficiency, 10,000 datasets were simulated.
