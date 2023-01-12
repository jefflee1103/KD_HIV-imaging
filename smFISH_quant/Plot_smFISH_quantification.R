# * * * * * Plot smFISH quantification results

library(tidyverse)
library(ggbeeswarm)

raw_csvs_files <- list.files(
    "./smFISH_quant/smFISH_quantification_output/final_outpus/",
    pattern = "*.csv",
    full.names = TRUE
    )

raw_csvs <- map_dfr(raw_csvs_files, read_csv)

annotated_df <- raw_csvs %>%
    mutate(time = case_when(
        str_detect(file_name, "24hpi") ~ "24hpi",
        str_detect(file_name, "48hpi") ~ "48hpi",
        str_detect(file_name, "72hpi") ~ "72hpi"
    )) %>%
    mutate(condition = case_when(
        str_detect(file_name, "WT_HIV") ~ "WT HIV",
        str_detect(file_name, "WT_no_virus") ~ "WT Mock",
        str_detect(file_name, "MKRN1KO") ~ "MKRN1KO HIV"
    )) %>%
    mutate(condition = fct_relevel(condition, c("WT Mock", "WT HIV", "MKRN1KO HIV")))

annotated_df %>%
    ggplot(aes(x = condition, y = total_RNAs + 1, colour = condition)) +
    geom_quasirandom(alpha = 0.5, stroke = 0) +
    geom_boxplot(width = 0.1) +
    labs(x = "Condition", y = "Intracellular Viral RNA (+1 log)", colour = "Condition") + 
    scale_y_log10() +
    theme_minimal()

