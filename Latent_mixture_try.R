library('gtools')
library('splines')
library('ggplot2')
library('tidyr')
library('dplyr')

# This script explores the construction of the linear predictor η in a Negative Binomial (NB) model.
# All model components are developed and interpreted on the link scale.

#set.seed(123)
x <- seq(1, 9, length.out = 200) # the windows, since they are between 1 and 9
spline <- TRUE # the 'error growth' will be a spline if arg == TRUE

# Setting a dynamic term to simulate the 3 models I would like to compare
# This term acts on the alpha vector of the Dirichlet distribution
importance <- c('same') # choose between 'same', 'prenatal', 'postnatal'
importance <- match.arg( # Troubleshooting
  importance,
  choices = c("same", "prenatal", "postnatal")
)

# Mixture components ------------------------------------------------------

# Pre e Post natal components different peaks: 
# Pre_peak could reasonably be around the 16PWC (knowledge from previous papers)
# Post_peak could be everywhere after birth, I couldn't find any useful indication in the literature
peak_prenatal <- rnorm(n = 1, mean = 3, sd = 0.2) # 3 is the value of the window for the 16PWC, it is generated from a hiherarchical prior
peak_postnatal <- rnorm(n = 1, mean = 6.5, sd = 0.6) # 6.5 is a random postnatal expression, since it seems there's no key peak-timing for postnatal component, thus its sd might be larger than the prenatal_sd

# Only pay attention on the sd: if it get's too large the hyperprior is noninformative and could give convergence issues

# Three components: an 'error' component, and hit-one and hit-two components
# The reasoning is an implementation of the two-hits hypothesis on SCZ with a source of 'error' in the biological theory (the linear/spline component)
# Thus, the f_error assumes that some genes do not strictly follow the two-hits logic
if(spline){
  f_error <- ns(x, df = 4) %*% c(1,1,1,1) # + 50 # spline growth, I used the spline to interpret it as an error
} else{
  f_error <- 0.2 * x # + 50 # linear growth for semplicity
}
f_prenatal <- 3 * exp(-((x - peak_prenatal)^2) / 0.6) #+ 50 # prenatal peak centered on the 16PWC
f_postnatal <- 3 * exp(-((x - peak_postnatal)^2) / 0.6)# + 50 # postnatal peak with hiher varaibility then the corresponding prenatal parameter
# il + 50 deriva dal fatto che l'espressione media nei precedenti modelli è ~50 indipendentemente dalla normalizzazione, centro la prior lì
# L'ho tolta perchè viene catturata dall'intercetta, a cui posso mettere una prior informativa su 50 se voglio
intercept <- rnorm(1, mean = 50, sd = 0.5) # Informative prior on the general intercept of the model

# Weights for the components -----------------------------------------------

# Weights from a Dirichlet distribution (aka multinomial model)
# This function determines the importance of each component and their variability 
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
eta <- intercept + 
  w_singolo[1]*f_error +
  w_singolo[2]*f_prenatal +
  w_singolo[3]*f_postnatal

# plotting components and final trajectory in the eta space
df_plot <- data.frame( # Putting all into one df
  x = x,
  mixture = eta,
  growth = f_error + intercept,
  prenatal = f_prenatal + intercept,
  postnatal = f_postnatal + intercept
)

df_long <- df_plot %>% # df into long format
  pivot_longer(cols = -x,
               names_to = "component",
               values_to = "value")

ggplot(df_long, aes(x = x, y = value)) +
  geom_line(aes(color = component), linewidth = 1.2) +
  geom_vline(aes(xintercept = 5, linetype = "birth"),
             color = "grey40") +
  scale_color_manual(
    values = c(
      mixture = "#1A1A1A",   # quasi nero ma più soft
      growth = "#7AA6D9",    # blu soft (classico scientifico)
      prenatal = "#F28C8C",  # rosso desaturato
      postnatal = "#8FD3CF"  # teal soft (più moderno del verde puro)
    ),
  labels = c(
      mixture = "mixture",
      growth = "error",
      prenatal = "prenatal",
      postnatal = "postnatal"
   ))+
  scale_linetype_manual(
    values = c(birth = "dashed"),
    labels = c(birth = "birth")
  ) +
  labs(
    title = paste0("Latent mixture trajectory - prior with ", importance, " importance"),
    x = "Window",
    y = "E[Y|X]",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )
