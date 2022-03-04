library(tidyverse)
library(broom)    # tidy
library(emmeans)  # emmeans, contrast
library(car)      # leveneTest
library(ARTool)   # Aligned-rank transform: art, artlm, art.con
library(WRS2)     # Robust methods: t2way, med2way, mcp2a
library(geepack)  # Generalized Estimating Equations

#===============================================================================
# Data from experiment 2 of Jansen & Hornbæk (2018) 
#      paper: https://doi.org/10.1145/3173574.3173588
#      data: https://github.com/yvonnejansen/posture
#
#   * Independent variable: Posture conditions (expansive vs. constrictive)
#   * Dependent variable: BART measure of risk-taking. It is the percentage change 
#                         of the average number of pumps before the balloon explodes
#                         between the the first 10 balloons (where the participant 
#                         is still calibrating how the system behaves) vs. the last 
#                         10 balloons (where the participants can exercise their
#                         risk-taking)
#                          
#   * Covariate: BIS (Barrett impulsiveness scale) as two categories (low vs. high) 
#                BIS was suggested by the literature to be correlated with BART measures
#
#   * Experimental design: Between-subjects (For each participant, 30 trials are aggregated to 1 BART measure)
#
# Although the paper mentions two covariates BIS and discomfort, 
# only the BIS covariate was presented in the paper and the dataset.
# The discomfort data and its analysis code are unavailable in the public material.
# Chat emailed the author a request (04.03.22).

posture_data <- 
  read_csv(
    "posture/exp2.csv",
    col_types = cols_only(
     participant = col_character(),
     condition = col_character(),
     BIS = col_character(),
     change = col_double())) %>% 
  mutate(
    condition = factor(condition),
    BIS = factor(BIS))

#-------------------------------------------------------------------------------
# visualize

posture_data %>% 
  ggplot(aes(x = change)) +
  geom_density() +
  facet_grid(condition ~ BIS)

#===============================================================================
# Analysis 1: Assuming normality (which is doubtful)
m <- lm(change ~ condition * BIS, 
        data = posture_data)

# model coefficients
conf_m_normal <- tidy(m, conf.int = TRUE)
conf_m_normal

# ANOVA
anova(m)

#-------------------------------------------------------------------------------
# estimates of the each contributing factor (Figure 11, but with a frequentist method)

emmeans(m, ~ 1) %>% 
  summary() %>% 
  ggplot(aes(x = "intercept", y = emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange() +
  coord_flip()

emmeans(m, ~ condition) %>% 
  contrast("eff", by = NULL) %>%
  summary(infer = TRUE) %>% 
  ggplot(aes(x = contrast, y = estimate, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange() +
  coord_flip()

emmeans(m, ~ BIS) %>% 
  contrast("eff", by = NULL) %>% 
  summary(infer = TRUE) %>% 
  ggplot(aes(x = contrast, y = estimate, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange() +
  coord_flip()

emmeans(m, ~ condition * BIS) %>%
  contrast("eff", by = NULL) %>% 
  summary(infer = TRUE, adjust = "bonferroni") %>% 
  ggplot(aes(x = contrast, y = estimate, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange() +
  coord_flip()


#===============================================================================
# Analysis 2: Log-transform (still doesn't solve normality problem)

# figure out the offset that makes all data positive to allow log() ==> 38
abs(min(posture_data$change)) + 1
# WARNING 1: The large offset is problematic for hypothesis testing 
#            (See: [Feng et al. (2014)](http://dx.doi.org/10.3969/j.issn.1002-0829.2014.02.009))

m_log <- lm(log(change + 38) ~ condition * BIS, 
        data = posture_data)

# assumption tests
m_log %>% 
  residuals() %>% 
  shapiro.test()
leveneTest(m_log, center = "median")

# WARNING 2: The log-transformed data still doesn't meet the normality assumption

# With these two warnings, we gave up the log-transformation


#===============================================================================
# Analysis 3: Aligned-rank transform

posture_art <- art(change ~ condition * BIS, 
            data = posture_data)
posture_art

# WARNING: F values of ANOVAs on aligned response indicates that ART is inappropriate

#===============================================================================
# Analysis 4: Robust two-way ANOVA

# trimmed mean test
t2way(change ~ condition * BIS, data = posture_data)
mcp2atm(change ~ condition * BIS, data = posture_data)

# median test
med2way(change ~ condition * BIS, data = posture_data, est = "median")
mcp2a(change ~ condition * BIS, data = posture_data, est = "median")

# M-estimator (doesn't work because of thesingular variance-covariance matrix)
# pbad2way(change ~ condition * BIS, data = posture_data, est = "mom", pro.dis = TRUE)
# mcp2a(change ~ condition * BIS, data = posture_data, est = "mom")


# Results: Two robust analyses (trimmed means, and median test) runs. 
#          The estimates of the interaction effect from these results 
#          indicates that a small interaction effect, but the uncertainty in the 
#          interval estimates makes this interaction effect not statistically significant.
#          When comparing the results of the robust method with those from the
#          model that falsely-assume normality (m, conf_m_normal), the robust method
#          is slightly more powerful in revealing the interaction effect. 



#===============================================================================
# Analysis 5: GEE

m_gee <- geeglm(change ~ condition * BIS, 
                id = participant,
                data = posture_data,
                family = gaussian())

tidy(m_gee, conf.int = TRUE)

anova(m_gee)


# NOTE: The `type = "response"` parameter applies back-transformation after all calculations

emmeans(m_gee, ~ 1, type = "response") %>% 
  summary() %>% 
  ggplot(aes(x = "intercept", y = response, ymin = asymp.LCL, ymax = asymp.UCL)) +
  geom_pointrange() +
  coord_flip()

emmeans(m_gee, ~ condition) %>% 
  contrast("eff", by = NULL) %>% 
  summary(infer = TRUE, type = "response")


# WARNING: Below, regrid() back-transforms the data before calculating contrasts.
#          This method yields symmetric confidence intervals---which I cannot think
#          of a sound justification. (The original data are not.)

emmeans(m_gee, ~ condition) %>% 
  regrid() %>% 
  contrast("eff", by = NULL) %>%
  summary(infer = TRUE) %>% 
  ggplot(aes(x = contrast, y = estimate, ymin = asymp.LCL, ymax = asymp.UCL)) +
  geom_pointrange() +
  coord_flip()

emmeans(m_gee, ~ BIS) %>% 
  regrid() %>% 
  contrast("eff", by = NULL) %>% 
  summary(infer = TRUE) %>% 
  ggplot(aes(x = contrast, y = estimate, ymin = asymp.LCL, ymax = asymp.UCL)) +
  geom_pointrange() +
  coord_flip()

emmeans(m_gee, ~ condition * BIS) %>%
  regrid() %>%
  contrast("eff", by = NULL) %>% 
  summary(infer = TRUE, adjust = "bonferroni") %>% 
  ggplot(aes(x = contrast, y = estimate, ymin = asymp.LCL, ymax = asymp.UCL)) +
  geom_pointrange() +
  coord_flip()
