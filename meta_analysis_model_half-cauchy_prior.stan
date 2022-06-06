//Stan model for meta analysis with different priors
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
    } 
