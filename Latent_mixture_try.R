library('gtools')
library('splines')
# This script explores the construction of the linear predictor η in a Negative Binomial (NB) model.
# All model components are developed and interpreted on the link scale.

#set.seed(123)
x <- seq(1, 9, length.out = 200) # the windows, since they are between 1 and 9
spline <- TRUE
importance <- c('same') # choose between 'same', 'prenatal', 'postnatal'

# Mixture components ------------------------------------------------------

# Pre e Post natal components have a peak: 
# Pre_peak could reasonably around the 16PWC (previous analyses)
# Post_peak could be everywhere, I couldn't find any useful indication 
peak_prenatal <- rnorm(n = 1, mean = 3, sd = 0.2) # 3 is the value of the window for the 16PWC, it is generated from a hiherarchical prior
peak_postnatal <- rnorm(n = 1, mean = 6.5, sd = 0.6) # 6.5 is a random postnatal expression, since it seems there's no key peak-timing for postnatal component, thus its sd might be larger than the prenatal_sd

# Only pay attention on the sd: if it get's too large the hyperprior is noninformative and could give convergence issues

# Three components: a linear component, and hit-one and hit-two components
# The reasoning is an implementation of the two-hits hypothesis on SCZ with a source of 'error' in the biological theory (the linear component)
# The f_growth assumes that some genes do not strictly follow the two-hits logic
if(spline){
  f_growth <- ns(x, df = 4) %*% c(1,1,1,1) + 50 # spline growth, I used the spline to interpret it as an error
} else{
  f_growth <- 0.2 * x + 50 # linear growth for semplicity
}
f_prenatal <- 3 * exp(-((x - peak_prenatal)^2) / 0.6) + 50 # prenatal peak (lo posso mettere sulla 16PWC visto che si è dimostrata importante)
f_postnatal <- 3 * exp(-((x - peak_postnatal)^2) / 0.6) + 50 # postnatal peak
# il + 50 deriva dal fatto che l'espressione media nei precedenti modelli è ~50 indipendentemente dalla normalizzazione, centro la prior lì


# Weights for the components -----------------------------------------------
importance <- match.arg( # Troubleshooting
  importance,
  choices = c("same", "prenatal", "postnatal")
)

# Weights from a Dirichlet distribution (aka multinomial model)
# This function determines the importance of each component
if(importance == 'same'){
  alpha_unif <- c(50, 50, 50) # same weight each component, I put 50 to have less volatile weights
  w_singolo <- rdirichlet(n = 1, alpha = alpha_unif)
} else if (importance == 'prenatal'){
  alpha_prenatal_priority <- c(10,50,10) # more weight to the prenatal component
  w_singolo <- rdirichlet(n = 1, alpha = alpha_prenatal_priority)
} else if (importance == 'postnatal'){
  alpha_postnatal_priority <- c(10,10,50) # more weight to the postnatal component
  w_singolo <- rdirichlet(n = 1, alpha = alpha_postnatal_priority)
}

# eta building & plots ----------------------------------------------------

# Creating the eta for the NB
eta <- w_singolo[1]*f_growth +
  w_singolo[2]*f_prenatal +
  w_singolo[3]*f_postnatal

# plotting components and final trajectory in the eta space
plot(x, eta, type = "l", lwd = 3,
     main = "Latent mixture trajectory - prior", ylab = 'E[Y|X]', xlab = 'Window', ylim = c(48,55))

lines(x, f_growth, col = "blue", lwd = 2) # linear component
lines(x, f_prenatal, col = "red", lwd = 2) # prenatal component
lines(x, f_postnatal, col = "green", lwd = 2) # postnatal component
abline(v = 5, col = "grey", lty = 2) # birth

# The weights and the components will be our priors
