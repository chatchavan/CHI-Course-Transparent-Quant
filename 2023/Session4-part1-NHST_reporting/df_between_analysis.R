# Scenario: You developed an electrical muscle stimulation (EMS) wristband that 
#           is supposed to help  the users type faster. To evaluate its effectiveness,
#           you conducted an experiment to compare three conditions: 
#             * off: no muscle stimulation
#             * constant: stimulate the wrist with a constant electrical intensity
#             * adaptive: adapt the intensity based on the real-time movement of the wrists
# 

library(tidyverse)
import::from(broom.mixed, tidy)  # install.packages("broom.mixed")


df_between <- read_csv("df_between.csv")
m_between <- lm(wpm ~ cond, data = df_between)

summary(m_between)                  # model summary with p-values
tidy(m_between, conf.int = TRUE)    # estimation of coefficients