# BIOS 509
# Spring 2023
# Poisson Regression - MLS Example

# Load Data ---------------------------------------------------------------
library(tidyverse)

# Calculate GC and eMP
mls <- read_csv('mls.csv') %>%
  filter(Min > 500) %>%
  filter(Pos != 'GK') %>%
  mutate(GC = Gls + Ast,
         eMP = Min / 90,
         Pos = factor(Pos, levels = c('DF','MF','FW')))

# EDA ---------------------------------------------------------------------

# Boxplots
mls %>%
  pivot_longer(cols = c('Gls','Ast','GC','Wage','Min','MP'), names_to = 'Variable', values_to = 'Value') %>%
  ggplot(aes(x = Value)) +
  geom_histogram(fill = 'grey', color = 'black') +
  theme_bw() +
  facet_wrap(~Variable, scales = 'free')
# ggsave('Figures/histograms.png', width = 5.5, height = 4.5)

# Log transform wage to make more normal
mls$lWage <- log(mls$Wage)
mls %>%
  pivot_longer(cols = c('Wage','lWage'), names_to = 'Variable', values_to = 'Value') %>%
  ggplot(aes(x = Value)) +
  geom_histogram(fill = 'grey', color = 'black') +
  theme_bw() +
  facet_wrap(~Variable, scales = 'free')
# ggsave('Figures/wagehist.png', width = 5.5, height = 4.5)

# Frequency of positions
table(mls$Pos)

# Compare MP to eMP
mls %>%
  ggplot(aes(x = MP, y = eMP)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1) +
  theme_bw() +
  labs(x = 'Matches Played',
       y = 'Effective Matches Played')
# ggsave('Figures/eMPvMP.png', width = 5.5, height = 4.5)

# Overall scatter plot
mls %>%
  ggplot(aes(x = lWage, y = GC)) +
  geom_point(alpha = 0.4) +
  # geom_smooth(method = 'lm') +
  theme_bw() +
  labs(x = 'Log Annual Wage ($100Ks)',
       y = 'Goal Contributions')
# ggsave('Figures/scatter.png', width = 5.5, height = 4.5)



# Poisson Regression ------------------------------------------------------

# Transform wage again - this allows for easier interpretation (for a 10% 
# increase in wage ...). This is not required for poisson reg, just helpful 
# since we log-transformed a PREDICTOR
# mls$lWage <- mls$lWage / log(1.1)

# Model P1: Wage
p1 <- glm(GC ~ offset(log(eMP)) + lWage, family = poisson, data = mls)
summary(p1)

# Estimates on the log and response scale
round(cbind(Est = coef(p1), confint.default(p1)), 3)
round(cbind(Est = exp(coef(p1)), exp(confint.default(p1))), 3)


# Scatter plot by Position
mls %>%
  ggplot(aes(x = lWage, y = GC, color = Pos)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_bw() +
  labs(x = 'Log Annual Wage ($100Ks)')

# Looks like there are different trends across positions. Maybe its not fair to 
# expect everyone to score the same amount. Let's add position to the model.

# Model P2: Wage + Position 
p2 <- glm(GC ~ offset(log(eMP)) + lWage + Pos, family = poisson, data = mls)
summary(p2)

# In fact, maybe more money != more goals for all positions.
# Model P3: Wage + Position + Wage*Position
p3 <- glm(GC ~ offset(log(eMP)) + lWage * Pos, family = poisson, data = mls)
summary(p3)

# We can compare the deviance across models.
# More terms = lower deviance (similar to why R^2 always increases in OLS)
anova(p3)

# Let check to see if the interaction terms are necessary
deviance(p2) - deviance(p3)               # Method 1: difference in deviance
as.numeric(-2*(logLik(p2) - logLik(p3)))  # Method 2: LRT (matches)
anova(p2, p3, test = 'Chisq')             # Shortcut
# It looks like the interaction terms are necessary.

# Parameter estimates for the interaction model (on the rate ratio scale)
round(cbind(Est = exp(coef(p3)), exp(confint.default(p3))), 3)

# This is great, but does the interaction model fit well? Use GOF tests
D <- p3$deviance
pchisq(D, p3$df.residual, lower.tail = FALSE)

X2 <- sum(resid(p3, type = "pearson")^2)
pchisq(X2, p3$df.residual, lower.tail = FALSE)

# Still very significant evidence that the model does not fit well.

# Maybe the mean model isn't the problem, but the variance model is. 


# Checking for Overdispersion ---------------------------------------------

# If either of the following are much greater than 1 (rule of thumb is > 1.5),
# there is evidence of overdispersion
phi.hat.X2 <- X2 / p3$df.residual
phi.hat.D <- D / p3$df.residual
phi.hat.X2
phi.hat.D

# Both greater than 1.5

# Cut the data into groups and plot the group means vs variances
mls %>%
  group_by(Pos, WageGroup = cut_number(Wage, 5)) %>%
  summarise(mean = mean(GC), var = var(GC)) %>%
  ggplot(aes(x = mean, y = var)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, color = 'red', se = FALSE) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, color = 'blue', se = FALSE) +
  xlim(c(0, 12)) +
  theme_bw() +
  labs(x = 'Mean Number of GC',
       y = 'Variance of Number of GC')
# ggsave('Figures/overdispersion.png', width = 5.5, height = 4.5)

# The dotted line is the reference for mean = variance. Clearly that is not the
# case. Whether the relationship is still a straight line is up for debate.



# Quasi-Poisson Regression ------------------------------------------------

# If we assume the variance is proportional to the mean, quasi-poisson is good

# Try both the interaction and non-interaction model. Sometimes the significance
# of terms can change due to increased SE, so it will be good to check.
qp2 <- glm(GC ~ offset(log(eMP)) + lWage + Pos, family = quasipoisson, data = mls)
qp3 <- glm(GC ~ offset(log(eMP)) + lWage * Pos, family = quasipoisson, data = mls)
summary(qp2)
summary(qp3)

# It looks like the interaction is still significant.
# The parameter estimates are the same as the poisson models.
# The SEs are larger (and we also are using T statistics now)

# How did R calculate these SEs?
sqrt(diag(vcov(qp3)))             # SEs from quasipoisson fit
sqrt(diag(vcov(p3)) * phi.hat.X2) # SEs from poisson fit scaled by our estimate of phi

# They match! We could use D to estimate the overdispersion as well, but the 
# default in R is the Pearson X2 version. If you go back and look at the model
# summary you will see the "Dispersion parameter" matches.
summary(qp3)$dispersion

# Parameter estimates for the interaction quasi-poisson model
round(cbind(Est = coef(qp3), confint.default(qp3)), 3)
round(cbind(Est = exp(coef(qp3)), exp(confint.default(qp3))), 3)

# We can test for the interaction using an F test for nested QP models
anova(qp2, qp3, test = "F") # interaction is still significant in the qp model

# Calculate this test by hand:
(deviance(qp2) - deviance(qp3) ) / 2 / phi.hat.X2




# Negative Binomial Regression --------------------------------------------

# What is we don't believe the dispersion is proportional to the mean? The
# negative binomial is a good alternative as it is essentially a Poisson model
# with an extra parameter modeling the variability.
nb2 <- MASS::glm.nb(GC ~ offset(log(eMP)) + lWage + Pos, data = mls)
nb3 <- MASS::glm.nb(GC ~ offset(log(eMP)) + lWage * Pos, data = mls)
summary(nb2)
summary(nb3)

# The estimates are similar to the other models, but the SEs are slightly
# lower than the Quasi-Poisson. Interaction is only borderline significant now.
round(cbind(Est = exp(coef(nb3)), exp(confint.default(nb3))), 3)

# Since the negative binomial is likelihood-based, we can conduct a LRT.
as.numeric(-2*(logLik(nb2) - logLik(nb3)))  # Method 1 (LRT)
(nb2$deviance - nb3$deviance)               # Method 2 (deviance)
# Warning! The drop-in-deviance test is no longer equivalent to the LRT.

anova(nb2, nb3, test = 'Chisq') # Shortcut to LRT

# Parameter estimates on the rate ratio scale
round(cbind(Est = exp(coef(nb2)), exp(confint.default(nb2))), 3)

# We can also conduct a LRT between a NB model and the matching poisson model. 
# This test has 1 df since the only difference is the dispersion parameter. 
# This test is one-sided since the dispersion parameter cannot be negative.
X2 <- as.numeric(-2*(logLik(p3) - logLik(nb3)))
pchisq(X2, df = 1, lower.tail = FALSE) / 2

# Significant evidence of overdispersion in the data



# Model Comparison --------------------------------------------------------

# We can compare AIC among candidate likelihood-based models
# Let's also throw linear models in
lm2 <- glm(GC ~ lWage + Pos, data = mls)
lm3 <- glm(GC ~ lWage * Pos, data = mls)

list(lm2, lm3, p2, p3, nb2, nb3) |>
  lapply(\(m) cbind(Family = m$family[[1]], Model = deparse1(m$call[[2]]), AIC = round(AIC(m),2))) |>
  do.call(what = rbind)

# Both NB models look good (low AIC). We will go with the simpler one.


# We might also consider plotting the model fits to see how "nice they look"

# Make predictions on test dataset
wagelist <- seq(min(mls$lWage), max(mls$lWage),.1)
test <- data.frame(lWage = rep(wagelist, times = 3), 
                   Pos = rep(c('DF','MF','FW'), each = length(wagelist)), 
                   eMP = 1)
lm.preds <- predict(lm3, test, se.fit = TRUE)
p.preds <- predict(p3, test, se.fit = TRUE)
qp.preds <- predict(qp3, test, se.fit = TRUE)
nb.preds <- predict(nb2, test, se.fit = TRUE)

preds <- bind_rows(test, test, test, test) %>%
  mutate(Model = rep(c('LM','Pois','QPois','NB'), each = nrow(test)),
         Pred = c(lm.preds$fit, p.preds$fit, qp.preds$fit, nb.preds$fit),
         SE = c(lm.preds$fit, p.preds$se.fit, qp.preds$se.fit, nb.preds$se.fit)) %>%
  mutate(l95 = Pred - 1.96*SE, u95 = Pred + 1.96*SE,
         Pred = ifelse(Model == 'LM', Pred, exp(Pred)),
         l95 = ifelse(Model == 'LM', l95, exp(l95)),
         u95 = ifelse(Model == 'LM', u95, exp(u95)))


# # Plot all predictions
# ggplot(NULL, aes(x = lWage)) +
#   geom_point(data = mls, aes(y = GC/eMP), alpha = 0.3) +
#   geom_line(data = preds, aes(y = Pred, color = Model), lty = 1) +
#   geom_line(data = preds, aes(y = l95, color = Model), lty = 2) +
#   geom_line(data = preds, aes(y = u95, color = Model), lty = 2) +
#   facet_wrap(~Pos, scales = 'free_y') +
#   theme_bw() +
#   theme(legend.position = 'top') +
#   labs(x = 'Log Annual Wage ($100,000s)',
#        y = 'Goal Contributions per 90')

# Remove the linear model fit
preds <- filter(preds, Model != 'LM')
ggplot(NULL, aes(x = lWage)) +
  geom_point(data = mls, aes(y = GC/eMP), alpha = 0.3) +
  geom_line(data = preds, aes(y = Pred, color = Model), lty = 1) +
  geom_line(data = preds, aes(y = l95, color = Model), lty = 2) +
  geom_line(data = preds, aes(y = u95, color = Model), lty = 2) +
  ylim(c(0, NA)) +
  facet_wrap(~Pos, scales = 'free_y') +
  theme_bw() +
  theme(legend.position = 'top') +
  labs(x = 'Log Annual Wage ($100,000s)',
       y = 'Goal Contributions per 90')
# ggsave('Figures/fits.png', width = 6.5, height = 4.5)


# Diagnostics -------------------------------------------------------------

# Store deviance residuals and fitted values
diag <- data.frame(Model = rep(c('LM','Pois','QPois','NB'), each = nrow(mls)),
                   #resid = c(resid(lm3), resid(p3), resid(qp3), resid(nb3)),
                   resid = c(rstandard(lm3), rstandard(p3), rstandard(qp3), rstandard(nb2)),
                   fitted = c(fitted(lm3), fitted(p3), fitted(qp3), fitted(nb2)))

# Examine linearity on the link scale
diag %>%
  ggplot(aes(y = resid, x = fitted)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  facet_wrap(~Model) +
  theme_bw() +
  labs(x = 'Model Prediction',
       y = 'Standardized Deviance Residual')
      # title = 'Checking Linearity on the Link Function Scale')
# ggsave('Figures/linear_link.png', width = 6.5, height = 4.5)

# All but the normal model seem ok

# Examine normality of deviance residuals
diag %>%
  ggplot(aes(sample = resid)) +
  geom_qq(alpha = 0.3) + 
  geom_qq_line(color = 'red') +
  facet_wrap(~Model, scales = 'free') +
  theme_bw() +
  labs(x = 'Normal Theoretical Quantile',
       y = 'Standardized Deviance Residual Quantile')
      # title = 'QQ-Plot of Deviance Residuals')
# ggsave('Figures/qq.png', width = 6.5, height = 4.5)

# all but the normal model seem ok. NB maybe a little better (see top right)

# Examine indepdendence assumption
diag %>%
  group_by(Model) %>%
  mutate(obsid = row_number()) %>%
  ggplot(aes(y = resid, x = obsid)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  facet_wrap(~Model, scales = 'free') +
  theme_bw() +
  labs(x = 'Observation Index in Sample',
       y = 'Standardized Deviance Residual')
     #  title = 'Checking Independence of Deviance Residuals')
# ggsave('Figures/indep.png', width = 6.5, height = 4.5)

# all seem ok

colSums(influence.measures(nb3)$is.inf)
