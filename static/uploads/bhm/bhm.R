# BIOS 525: Bayesian Hierarchical Model
# Jacob Englert
# 31 October 2022

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(R2jags)
library(lme4)

# Gibbs Sampler for Normal Mean and Variance ------------------------------

# Data (outcome) vector
y <- c(4.2, 6.6, 5.1, 2.0, 2.8, 3.2, 4.7, 4.1, 7.3, 0.8)
n <- length(y)    # Sample size

# Parameters for Y
mu <- 0           # Initial values for mu
sigma2 <- 1       # Initial values for sigma2

# Priors for mu
mu0 <- 0	        # Prior mean for mu
tau2 <- 5^2       # Prior variance for mu

# Priors for sigma2
alpha <- 10       # Prior alpha value for sigma2
beta <- 50        # Prior beta value for sigma2

# Number of MCMC samples
n.iter <- 5000

# Allocate storage for posterior samples
mu.save <- numeric(n.iter)
sigma2.save <- numeric(n.iter)

# Run Gibbs Sampler!
set.seed(2187)
for(i in 1:n.iter){
  
  # Step 1) Update mu
  mu.V <- 1 / (n/sigma2 + 1/tau2)
  mu.M <- mu.V * (n/sigma2*mean(y) + mu0 / tau2)
  mu <- rnorm(1, mean = mu.M, sd = sqrt(mu.V))
  
  # Step 2) Update sigma2
  sigma2 <- 1 / rgamma(1, shape = n/2 + alpha, rate = sum((y - mu)^2)/2 + beta)
  
  # Store updates
  mu.save[i] <- mu
  sigma2.save[i] <- sigma2
}

post <- data.frame(mu = mu.save, sigma2 = sigma2.save)
  
# Trace Plots
post %>%
  mutate(Iteration = 1:n.iter) %>%
  pivot_longer(cols = -Iteration, names_to = 'Type') %>%
  ggplot(aes(x = Iteration, y = value)) +
  geom_line() +
  facet_wrap(~Type, nrow = 2, scales = 'free') +
  theme_bw()

# Joint Posterior Distribution
joint <- ggplot(mapping = aes(x = mu, y = sigma2)) +
  geom_point(data = post, alpha = 0.2) +
  theme_bw() +
  labs(x = expression('Population Mean'~mu),
       y = expression('Population Variance'~sigma^2))
joint.arrow <- joint +
  geom_path(data = post[1:100,], color = 1:100, arrow = arrow())

joint.hist <- ggExtra::ggMarginal(joint, type = 'histogram')
joint.arrow.hist <- ggExtra::ggMarginal(joint.arrow, type = 'histogram')


# Bayesian Linear Regression ----------------------------------------------

# Simulate linear regression data
set.seed(1031)
n <- 100
X <- runif(n, 0, 5)
Y <- 2 + 2*X + rnorm(n)

data.frame(X, Y) %>% ggplot(aes(x = X, y = Y)) + geom_point() + theme_bw()

# Specify JAGS model
jags.SLR <- function(){
  
  # Data Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(mu[i], sigma2.inv)
    mu[i] = beta0 + beta1 * X[i]
  }
  
  # Prior distribution
  beta0 ~ dnorm(0, 1e-06)
  beta1 ~ dnorm(0, 1e-06)
  sigma2.inv ~ dgamma(0.0001, 0.0001)
  sigma2 = 1 / sigma2.inv
}

# Fit JAGS model
SLR.fit <- jags(data = list("Y", "X", "n"), model.file = jags.SLR,
                parameters.to.save = c("beta0", "beta1", "sigma2"), DIC = FALSE,
                n.chain = 1, n.iter = 20000, n.burnin = 10000, n.thin = 1,
                jags.seed = 2187)
SLR.mcmc <- as.mcmc(SLR.fit) # MCMC object (helpful for quick summaries)
plot(SLR.mcmc)
summary(SLR.mcmc)
summary(lm(Y ~ X))



# Toenail Example ---------------------------------------------------------
toe <- read_table("~/OneDrive - Emory University/Documents/Courses/BIOS 560R/Data/toenail.txt")

# Visualize Data
toe %>%
  group_by(treat, time) %>%
  mutate(mean.response = mean(response)) %>%
  ggplot(aes(x = time, group = id)) +
  geom_line(aes(y = response), alpha = 0.1, color = 'blue') +
  geom_line(aes(y = mean.response, color = 'Time-specific Mean')) +
  geom_point(aes(y = mean.response), color = 'red') +
  facet_wrap(~paste('Treat =', treat), nrow = 1) +
  scale_x_continuous(breaks = unique(toe$time)) +
  scale_color_manual(values = c('Time-specific Mean' = 'red')) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = 'Week', y = 'Response', color = '')

# Random Intercept Model
toe.jags1 <- function(){

  # Data Likelihood
  for(i in 1:n.obs){
    response[i] ~ dnorm(mu[i], prec1)
    mu[i] = alpha[id[i]] + beta1*treat[i] + beta2*time[i]+ beta3*time[i]*treat[i]
  }

  for(j in 1:n.id){ alpha[j] ~ dnorm(alpha0, prec2) }

  # Priors
  alpha0 ~ dnorm(0, 1E-10)
  beta1 ~ dnorm(0, 1E-10)
  beta2 ~ dnorm(0, 1E-10)
  beta3 ~ dnorm(0, 1E-10)
  prec1 ~ dgamma(0.0001, 0.0001)
  prec2 ~ dgamma(0.0001, 0.0001)
  sigma2 = 1/prec1
  tau2 = 1/prec2
}

# JAGS data input
toe$id <- as.numeric(as.factor(toe$id)) # Recode ID (so as to not skip any IDs)
response <- toe$response
treat <- toe$treat
time <- toe$time
id <- toe$id
n.id <- max(toe$id)
n.obs <- nrow(toe)

toe.fit1 <- jags(data = list("response", "treat", "time", "id", "n.id", "n.obs"),
                 model.file = toe.jags1, DIC = FALSE, jags.seed = 1031,
                 parameters.to.save = c("alpha0","beta1","beta2","beta3","sigma2","tau2","alpha[1:2]"),
                 n.chain = 1, n.iter = 20000, n.burnin = 10000, n.thin = 4)
toe.mcmc0 <- as.mcmc(toe.fit1)
toe.post0 <- data.frame(as.matrix(toe.mcmc0))

# Create versions without random intercepts for cleaner output
toe.mcmc1 <- as.mcmc(toe.post0[,c("alpha0","beta1","beta2","beta3","sigma2","tau2")])
toe.post1 <- data.frame(as.matrix(toe.mcmc1))

# Data frame of posterior samples
head(toe.post1)

# Assess convergence of posterior samples
plot(toe.mcmc1)

# Summarize posterior samples of parameters
summary(toe.mcmc1)

# Compare to frequentist model
toe.lmer <- lmer(response ~ (1|id) + treat + time + treat*time, data = toe)
summary(toe.lmer)

# Visualize pairwise joint posteriors
GGally::ggpairs(toe.post1, upper = list(continuous = "density"),
                lower = list(continuous = "points"),
                mapping = aes(alpha = 0.001)) + 
  theme_bw()

# 95% CrI of exp(beta1) - beta2
new.samp1 <- exp(toe.post1$beta1) - toe.post1$beta2
mean(new.samp1)
quantile(new.samp1, c(0.025, 0.975))

data.frame(Post = new.samp1) %>%
  ggplot(aes(x = Post)) +
  geom_histogram(fill = 'white', color = 'black') +
  theme_bw() +
  labs(title = 'Posterior Distribution of exp(beta1) - beta2')

# 95% CrI for ICC
toe.icc <- toe.post1$tau2 / (toe.post1$tau2 + toe.post1$sigma2)
mean(toe.icc)
quantile(toe.icc, c(0.025, 0.975))

data.frame(Post = toe.icc) %>%
  ggplot(aes(x = Post)) +
  geom_histogram(fill = 'white', color = 'black') +
  theme_bw() +
  labs(title = 'Posterior Distribution of ICC')

# Pr(ICC > 0.5)
mean(toe.icc > 0.5)

data.frame(Post = toe.icc) %>%
  ggplot(aes(x = Post)) +
  geom_histogram(fill = 'white', color = 'black') +
  geom_vline(xintercept = 0.5, color = 'red') +
  theme_bw() +
  labs(title = 'Posterior Distribution of ICC')

# Compare random intercepts for subjects 1 and 2
mean(toe.post0$alpha.1. > toe.post0$alpha.2.)
data.frame(post = c(toe.post0$alpha.1., toe.post0$alpha.2.),
           id = rep(1:2, each = nrow(toe.post0))) %>%
  ggplot(aes(x = post, fill = paste0('ID = ', id))) +
  geom_histogram(alpha = 0.3, position = "identity") +
  theme_bw() +
  theme(legend.position = c(.1, .8)) +
  labs(fill = 'ID')


# In- and Out-of-Sample Predictions
toe.jags2 <- function(){
  
  # Data Likelihood
  for(i in 1:n.obs){
    response[i] ~ dnorm(mu[i], prec1)
    mu[i] = alpha[id[i]] + beta1*treat[i] + beta2*time[i]+ beta3*time[i]*treat[i]
  }
  
  for(j in 1:n.id){ alpha[j] ~ dnorm(alpha0, prec2) }
  
  # In-sample prediction
  for(k in 1:n.id){
    pred.mu_IN[k] = alpha[k] + beta1 + beta2*15 + beta3*15
    pred_IN[k] ~ dnorm(pred.mu_IN[k], prec1)
  }

  # Out-sample prediction
  pred.alpha ~ dnorm(alpha0, prec2)
  pred.mu_OUT = pred.alpha + beta1 + beta2*15 + beta3*15
  pred_OUT ~ dnorm(pred.mu_OUT, prec1)
  
  # Priors
  alpha0 ~ dnorm(0, 1E-10)
  beta1 ~ dnorm(0, 1E-10)
  beta2 ~ dnorm(0, 1E-10)
  beta3 ~ dnorm(0, 1E-10)
  prec1 ~ dgamma(0.0001, 0.0001)
  prec2 ~ dgamma(0.0001, 0.0001)
  sigma2 = 1/prec1
  tau2 = 1/prec2
}

toe.fit2 <- jags(data = list("response", "treat", "time", "id", "n.id", "n.obs"),
                 model.file = toe.jags2, DIC = FALSE, jags.seed = 1031,
                 parameters.to.save = c("pred_IN[1:3]", "pred_OUT"),
                 n.chain = 1, n.iter = 20000, n.burnin = 10000, n.thin = 4)
toe.mcmc2 <- as.mcmc(toe.fit2)

summary(toe.mcmc2)


# Random Intercept and Slope Model
toe.jags3 <- function(){

  # Data Likelihood
  for(i in 1:n.obs){
    response[i] ~ dnorm(mu[i], prec1)
    mu[i] = alpha[ id[i],1] + beta1*treat[i] + alpha[id[i],2]*time[i]+ beta2*time[i]*treat[i]
  }

  for (j in 1:n.id){ alpha[j,1:2] ~ dmnorm(alpha_mu, Omega) }

  # Priors
  alpha_mu[1] ~ dnorm(0, 1E-10)
  alpha_mu[2] ~ dnorm(0, 1E-10)
  beta1 ~ dnorm(0, 1E-10)
  beta2 ~ dnorm(0, 1E-10)
  prec1 ~ dgamma(0.0001, 0.0001)
  Omega[1:2, 1:2] ~ dwish(R, 2)
  sigma = 1 / prec1
  BigSigma = inverse(Omega)
}

# Get estimate for R
toe.lmer2 <- lmer(response ~ (1 + time|id) + treat + time + treat*time, data = toe)
summary(toe.lmer2)
R <- 2 * matrix(c(7.37, -0.39, -0.39, 0.23), ncol = 2)

toe.fit3 <- jags(data = list("response", "treat", "time", "id", "n.obs", "n.id", "R"), 
                 model.file = toe.jags3, DIC = FALSE, jags.seed = 1031,
                 parameters.to.save = c("alpha_mu", "beta1", "beta2", "sigma", "BigSigma"),
                 n.chain = 1, n.iter = 20000, n.burnin = 10000, n.thin = 4)
toe.mcmc3 <- as.mcmc(toe.fit3)
toe.post3 <- data.frame(as.matrix(toe.mcmc3))

plot(toe.mcmc3)
summary(toe.mcmc3)

post.corr <- toe.post3$BigSigma.2.1. / sqrt(toe.post3$BigSigma.1.1. * toe.post3$BigSigma.2.2.)
mean(post.corr)
quantile(post.corr, c(0.025, .975))




# Crossover Trial Example -------------------------------------------------

cross <- read_table("~/OneDrive - Emory University/Documents/Courses/BIOS 525/Data/2by2.txt")

# JAGS Model
cross.jags <- function(){
  
  # Data Likelihood
  for(i in 1:n.obs){
    outcome[i] ~ dbern(p[i])
    logit(p[i]) = beta0 + theta[ID[i]] + beta1*trt[i] + beta2*period[i] + beta3*trt[i]*period[i]
  }
  
  for(s in 1:n.ID){theta[s] ~ dnorm(0, tau2.inv)}
  
  # Priors
  beta0 ~ dnorm(0, 1e-10)
  beta1 ~ dnorm(0, 1e-10)
  beta2 ~ dnorm(0, 1e-10)
  beta3 ~ dnorm(0, 1e-10)
  tau2.inv ~ dgamma(0.0001, 0.0001)
  tau2 = 1 / tau2.inv
}

# JAGS data input
outcome <- cross$outcome
ID <- cross$ID
trt <- cross$trt
period <- cross$period
n.obs <- nrow(cross)
n.ID <- max(ID)

cross.fit <- jags(data = list ("outcome", "ID", "trt", "period", "n.obs", "n.ID"), 
                  model.file = cross.jags, DIC = FALSE, jags.seed = 1031,
                  parameters.to.save = c("beta0", "beta1", "beta2", "beta3", "tau2"),
                  n.chain = 1, n.iter = 40000, n.burnin = 20000, n.thin = 4)
cross.mcmc <- as.mcmc(cross.fit)
cross.post <- data.frame(as.matrix(cross.mcmc))
summary(cross.mcmc)
plot(cross.mcmc)

# Frequentist model
cross.fit0 <- glmer(outcome ~ trt*period + (1|ID), nAGQ = 100,
                    family = binomial(link = 'logit'),
                    data = cross)
summary(cross.fit0)

# Lung Cancer Example -----------------------------------------------------

lung <- read_csv("~/OneDrive - Emory University/Documents/Courses/BIOS 560R/Data/Cancer.csv")

# JAGS Model
lung.jags <- function(){
  
  # Data Likelihood
  for(i in 1:n.obs){
    death[i] ~ dpois(lambda[i]*pop[i])
    log(lambda[i]) = beta0 + theta[cty[i]] + beta1*sex[i] + beta2*race[i] + beta3*sex[i]*race[i] + beta4*year[i]
  }
  
  for(s in 1:n.cty){theta[s] ~ dnorm(0, tau2.inv)}
  
  # Priors
  beta0 ~ dnorm(0, 1e-10)
  beta1 ~ dnorm(0, 1e-10)
  beta2 ~ dnorm(0, 1e-10)
  beta3 ~ dnorm(0, 1e-10)
  beta4 ~ dnorm(0, 1e-10)
  tau2.inv ~ dgamma(0.0001, 0.0001)
  tau2 = 1 / tau2.inv
}

# JAGS data input
death <- lung$death
pop <- lung$pop
sex <- lung$sex - 1
race <- lung$race - 1
year <- lung$year
cty <- lung$county
n.obs <- nrow(lung)
n.cty <- max(lung$county)

lung.fit <- jags(data = list ("death", "pop", "sex", "race", "year", "cty", "n.obs", "n.cty"), 
                 model.file = lung.jags, DIC = FALSE, jags.seed = 1031,
                 parameters.to.save = c("beta0", "beta1", "beta2", "beta3", "beta4", "tau2"),
                 n.chain = 1, n.iter = 40000, n.burnin = 20000, n.thin = 10)
lung.mcmc <- as.mcmc(lung.fit)
lung.post <- data.frame(as.matrix(lung.mcmc))
summary(lung.mcmc)
plot(lung.mcmc)

# Frequentist model
lung.fit0 <- glmer(death ~ (1|county) + factor(sex) + factor(race) +
                          offset(log(pop))  + factor(sex)*factor(race) + year,
                   nAGQ = 100, family = poisson(link = "log"),
                   data = lung)
summary(lung.fit0)
