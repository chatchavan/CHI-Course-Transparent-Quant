# Answers to Activity 2.3

- The density curves show both parameter and prediction uncertainty because 
there are many prediction distributions implied by the posterior. Each curve is a 
t-distribution with different mu (parameter uncertainty) and all curves combined
show the uncertainty in the outcome (prediction uncertainty).
- An alternative is to visualize a single posterior predictive distribution, 
which is integrated over parameter uncertainty.


```{r}
df1 %>%
  ggplot(aes(x = effectiveness, y = condition)) + 
  stat_interval(
    aes(x = .prediction), 
    data = preds_t_bayes) +
  stat_pointinterval(
    aes(x = .epred), 
    data = means_t_bayes,
    .width = c(.66, .95),
    position = position_nudge(y = -0.2)) +
  geom_point(position = jitterer) +
  # END EXERCISE
  scale_x_continuous(breaks = 1:9)+
  scale_color_brewer("Predictive interval") + 
  labs(title = "Bayesian lm(), estimated effectiveness by condition",
       subtitle = "Credible interval widths: 0.66 and 0.95") +
  NULL
```
