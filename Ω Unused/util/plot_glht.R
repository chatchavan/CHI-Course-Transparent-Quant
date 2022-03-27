plot_glht <- function(glht_confint_tidy, unit = NA, ylim = NA, hide_y_title = FALSE, adjustment_method = NA, reference_value = 0) {
  label_unit <- ""
  adj_text <- "adjusted with the single-step method"
  
  if(!missing(unit))
    label_unit <- sprintf("in %s ", unit)
  
  if(!missing(adjustment_method))
    adj_text <- adjustment_method
  
  y_label <- sprintf("Estimate of the difference %s \n with 95%% CI (%s)", label_unit, adj_text)
  
  glht_confint_tidy %>% 
    ggplot(aes(x = contrast, y = estimate, ymin = conf.low, ymax = conf.high)) +
    geom_hline(yintercept = reference_value, color = "red") +
    geom_pointrange() +
    expand_limits(y = reference_value) +
    coord_flip() +
    ylab(y_label)
}
