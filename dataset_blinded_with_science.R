library(tidyverse)

dataset_blinded_with_science <- function() {
  read_csv("blinded.csv",
           col_types = cols(
             experiment = col_integer(),
             condition = col_character(),
             effectiveness = col_integer())) %>%
    
    # select only experiment 1 (out of 4)
    filter(experiment == "1") %>%             
    mutate(
      
      condition = factor(condition),
      
      # adding a synthetic participant ID to clearly communicate that this data is from a between-subjects design
      participant_id = str_c("P", str_pad(row_number(), 2, pad = "0"))) %>% 
    
    select(participant_id, condition, effectiveness)
  }
