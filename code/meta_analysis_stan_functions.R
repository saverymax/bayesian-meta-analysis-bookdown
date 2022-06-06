##############################################################################
# Script containing functions and model definitions used for meta-analysis
##############################################################################

plot_priors <- function(){
    p <- ggplot(data.frame(x=c(0, 5)), aes(x)) + 
        stat_function(fun=dinvgamma, args=list(alpha=1, beta=.015), 
                      aes(colour="IG"), lwd=1.5)+
        stat_function(fun=dhcauchy, args=list(sigma=0.3),
                      aes(colour="HC"),lwd=1.5)+
        stat_function(fun=dht, args=list(nu=01, sigma=0.2),
                      aes(colour="HT"),lwd=1.5)+
        scale_color_manual("Dist", values=c("IG"="#6C97F6", "HC"="#F1A722", "HT"="#94EF7C")) +
        theme_bw()
    print(p)
}

# Flat prior by default
flat_path <- "meta_analysis_model_flat_prior.stan"
flat_string <- "// Stan model for meta analysis
    data {
        int J;
        vector[J] est;
        vector[J] se;
    }
    parameters {
        real mu;
        real<lower=0, upper = 10> tau;
        vector<offset=mu, multiplier=tau>[J] theta;
    }
    model {
        est ~ normal(theta, se);
        theta ~ normal(mu, tau);
        mu ~ normal(0, 1);
    }
    generated quantities {
        real theta_new = normal_rng(mu, tau);
    }
"

# inverse-gamma prior
# From Stan documentation:
# https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-scale-parameters-in-hierarchical-models
# https://discourse.mc-stan.org/t/prior-recommendation-for-scale-parameters-in-hierarchical-models-too-strong/2927/26
ig_path <- "meta_analysis_model_inverse-gamma_prior.stan"
ig_string <- " //Stan model for meta analysis with different priors
    data {
        int J;
        vector[J] est;
        vector[J] se;
    }
    parameters {
        real mu;
        real<lower=0> tau;
        vector<offset=mu, multiplier=tau>[J] theta;
    }
    model {
        est ~ normal(theta, se);
        theta ~ normal(mu, tau);
        mu ~ normal(0, 1);
        tau ~ inv_gamma(1, 0.15);
    }
    generated quantities {
        real theta_new = normal_rng(mu, tau);
    } "

# Half-cauchy prior
# How to fit a half-cauchy in stan?
# https://distribution-explorer.github.io/continuous/halfcauchy.html
hc_path <- "meta_analysis_model_half-cauchy_prior.stan"
hc_string <- "//Stan model for meta analysis with different priors
    data {
        int J;
        vector[J] est;
        vector[J] se;
    }
    parameters {
        real mu;
        real<lower=0> tau;
        vector<offset=mu, multiplier=tau>[J] theta;
    }
    model {
        est ~ normal(theta, se);
        theta ~ normal(mu, tau);
        mu ~ normal(0, 1.0);
        tau ~ cauchy(0, 0.3);
    }
    generated quantities {
        real theta_new = normal_rng(mu, tau);
    } "

# Half-t prior
ht_path <- "meta_analysis_model_half-t_prior.stan"
ht_string <- "//Stan model for meta analysis with different priors
    data {
        int J;
        vector[J] est;
        vector[J] se;
    }
    parameters {
        real mu;
        real<lower=0> tau;
        vector<offset=mu, multiplier=tau>[J] theta;
    }
    model {
        est ~ normal(theta, se);
        theta ~ normal(mu, tau);
        mu ~ normal(0, 1.0);
        tau ~ student_t(10, 0, 0.2);
    }
    generated quantities {
        real theta_new = normal_rng(mu, tau);
    } "

# Main function for running all modelling, with provided stan model
run_model <- function(model_path, model_string, model_desc){
    
    # Debugging lines
    #model_path=flat_path
    #model_string=flat_string
    
    data <- read.table("gelman_meta_analysis_data.txt", header=TRUE)
    write(model_string, model_path)
    stan_data <- list(est=data$est, se=data$se, J=nrow(data))
    
    # rstan usage
    # Should be nearly compatible with all the code below, except for cmdstanr block
    #stanc(model_path)
    # The choice of iterations and thinning is a memory issue for my computer when generating .RMD
    #fit <- stan(file = model_path, data = stan_data, warmup = 500, iter = 5000, chains = 3)
    #fit_summary <- summary(fit)
    #fit_cols <- as.data.frame(round(fit_summary$summary[1:14,c(1, 3, 4, 8)],4))
    #posterior <- rstan::extract(fit, inc_warmup = FALSE, permuted = FALSE)
    
    # cmdstanr usage 
    # Will only diverge from rstan in printing estimates
    model <- cmdstan_model(model_path)
    fit <- model$sample(data=stan_data, seed=13, chains=3, iter_sampling=1000, iter_warmup=100)
    # Use draws to get the samples from cmdstanr object
    
    return(fit)
}

print_estimates <- function(model_desc, fit){
    # Print posterior estimates and RR ratio 
    fit_summary <- fit$summary()
    fit_cols <- fit_summary[2:15,c(1, 2, 4, 6, 7)] %>% mutate(across(2:5, round, 4))
    rr <- fit_summary[c(2, 4:15), c(1,2,4,6,7)] %>% mutate(across(2:5, exp)) %>% mutate(across(2:5, round, 4))
    caption <- paste("Posterior estimates, with", model_desc, "prior on tau") 
    label <- paste("post-est-", model_desc, sep="")
    #label <- paste(model_desc)
    print(kbl(fit_cols, booktabs = T, escape=T, caption=caption, label=label)) 
    caption <- paste("Relative Risk posterior estimates, with", model_desc, "prior on tau") 
    label <- paste("exp-post-est-", model_desc, sep="")
    print(kbl(rr, booktabs = T, escape=T, caption=caption, label=label)) 
    
}

plot_posteriors <- function(model_desc, fit){
    posterior <- fit$draws()
    plot_title <- ggtitle(paste("Posterior distributions,", model_desc, "prior on tau"), "with means and 90% interval")
    p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", regex_pars = c("theta", "mu")) + plot_title
    print(p_post)
}

plot_tau_posterior <- function(model_desc, fit){
    posterior <- fit$draws()
    plot_title <- ggtitle(paste("Posterior distribution of tau,", model_desc, "prior"), "with mean and 90% interval")
    p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", pars = c("tau")) + plot_title
    print(p_post)
}

plot_convergence <- function(model_desc, fit){
    posterior <- fit$draws()
    # Use bayesplot to look at convergence.
    color_scheme_set("mix-blue-pink")
    p_trace <- mcmc_trace(posterior,  pars = c("theta_new", "mu", "tau"),
               facet_args = list(nrow = 2, labeller = label_parsed))
    print(p_trace + facet_text(size = 15))
}
    
plot_post_intervals <- function(model_desc, fit){
    posterior <- fit$draws()
    # Posterior intervals
    #interval_title <- ggtitle(paste("Posterior credible intervals,", model_desc, "prior on tau"))
    p_ints <- mcmc_intervals(exp(posterior), regex_pars = c("theta", "mu"))
    print(p_ints)
}

