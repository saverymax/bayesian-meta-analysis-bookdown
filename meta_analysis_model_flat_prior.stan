// Stan model for meta analysis
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

