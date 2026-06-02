library('gtools')
library('splines')
# This script explores the construction of the linear predictor η in a Negative Binomial (NB) model.
# All model components are developed and interpreted on the link scale.

#set.seed(123)
x <- seq(1, 9, length.out = 200) # the windows, since they are between 1 and 9
spline <- TRUE
importance <- c('same') # choose between 'same', 'prenatal', 'postnatal'

# Mixture components ------------------------------------------------------

# Three components: a linear component, and hit-one and hit-two components
# The reasoning is an implementation of the two-hits hypothesis on SCZ with a source of 'error' in the biological theory (the linear component)
# The f_growth assumes that some genes do not strictly follow the two-hits logic
if(spline){
  f_growth <- ns(x, df = 4) %*% c(1,1,1,1) # spline growth, I used the spline to interpret it as an error
} else{
  f_growth <- 0.2 * x # linear growth for semplicity
}
f_prenatal <- 3 * exp(-((x - 3)^2) / 0.6) # prenatal peak
f_postnatal <- -2 * exp(-((x - 6.5)^2) / 0.8) # postnatal peak

# Weights for the compoents -----------------------------------------------
importance <- match.arg( # Troubleshooting
  importance,
  choices = c("same", "prenatal", "postnatal")
)
# Weights from a Dirichlet distribution (aka multinomial model)
# This function determines the importance of each component
if(importance == 'same'){
  alpha_unif <- c(1, 1, 1) # same weight each component
  w_singolo <- rdirichlet(n = 1, alpha = alpha_unif)
} else if (importance == 'prenatal'){
  alpha_prenatal_priority <- c(1,10,1) # more weight to the prenatal component
  w_singolo <- rdirichlet(n = 1, alpha = alpha_prenatal_priority)
} else if (importance == 'postnatal'){
  alpha_postnatal_priority <- c(1,1,10) # more weight to the postnatal component
  w_singolo <- rdirichlet(n = 1, alpha = alpha_postnatal_priority)
}

# eta building & plots ----------------------------------------------------

# Creating the eta for the NB
eta <- w_singolo[1]*f_growth +
  w_singolo[2]*f_prenatal +
  w_singolo[3]*f_postnatal

# plotting components and final trajectory in the eta space
plot(x, eta, type = "l", lwd = 3,
     main = "Latent mixture trajectory", ylab = 'eta of the NB', xlab = 'Window', ylim = c(-2,3))

lines(x, f_growth, col = "blue", lwd = 2) # linear component
lines(x, f_prenatal, col = "red", lwd = 2) # prenatal component
lines(x, f_postnatal, col = "green", lwd = 2) # postnatal component
abline(v = 5, col = "grey", lty = 2) # birth

# The weights and the components will be our priors
