##############################################################################
# This script runs a Bayesian Meta analysis in Stan
# and conduct a prior sensitivity analysis 
# https://statmodeling.stat.columbia.edu/2022/02/28/answering-some-questions-about-meta-analysis-using-ivermectin-as-an-example/
# On priors for meta analysis:
# https://statmodeling.stat.columbia.edu/2022/02/28/priors-for-meta-analysis/
# https://mc-stan.org/docs/2_29/stan-users-guide/meta-analysis.html
# Stan documentation:
# https://mc-stan.org/docs/2_29/stan-users-guide-2_29.pdf

# This script runs the same code as the markdown file, but just as a script
# that prints out the results
##############################################################################
#library("rstan")
library("cmdstanr")
library("bayesplot")
library("ggplot2")
library("extraDistr")
library("dplyr")
library("tidyr")
library("kableExtra")

user_dir <- "d:/asus_documents/ku_leuven/courses/meta_analysis/project/bayesian-meta-analysis/src"
setwd(user_dir)
# Functions and model definitions are in here:
source("meta_analysis_stan_functions.R")

# Check out the priors
plot_priors()

# Run each model and print out tables/plots of estimates
run_model(flat_path, flat_string, "uniform")
run_model(ig_path, ig_string, "inverse-gamma")
run_model(hc_path, hc_string, "half-cauchy")
run_model(ht_path, ht_string, "half-t")
    

