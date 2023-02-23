library(tidyverse)

## Read in data
data_raw <- list.files(
    "./gag_line_profile/data/Line_quantification_Gag", 
    pattern = "*.csv",
    full.names = TRUE
    ) %>%
    set_names(~basename(.x) %>% str_remove(".csv")) %>%
    map(~read_csv(.x) %>% setNames(c("distance", "gray_value")))

## Normalise intensity and length
data_normalised <- data_raw %>%
    map(~{
        max_intensity <- max(.x$gray_value)
        min_dist <- min(.x$distance)
        max_dist <- max(.x$distance)
        .x %>% 
            mutate(norm_intensity = gray_value/max_intensity) %>%
            mutate(norm_distance = 100 * (distance - min_dist) / (max_dist - min_dist))
    })

## Discretise distance
steps <- 0.25
discrete_range = seq(0, 100, steps)

data_discretised <- data_normalised %>%
    map(~{
        .x$discrete_distance <- cut(.x$norm_distance, breaks = discrete_range, labels = FALSE, include.lowest = TRUE)
        .x %>%
            mutate(discrete_distance = steps * (discrete_distance - 1)) %>%
        group_by(discrete_distance) %>%
        summarise(
            norm_intensity_by_discrete_distance = median(norm_intensity),
            gray_value_by_discrete_distance = median(gray_value)
            ) 
    })

## Merge to a single dataframe
df_discretised <- data_discretised %>%
    map_dfr(~.x, .id = "file_name") %>%
    mutate(time = case_when(
        str_detect(file_name, "24hpi") ~ "24hpi",
        str_detect(file_name, "48hpi") ~ "48hpi",
        str_detect(file_name, "72hpi") ~ "72hpi"
    )) %>%
    mutate(channel = case_when(
        str_detect(file_name, "rna") ~ "RNA",
        str_detect(file_name, "gag") ~ "Gag",
    )) %>%
    mutate(condition = case_when(
        str_detect(file_name, "WT_HIV") ~ "WT HIV",
        str_detect(file_name, "WT_mock") ~ "WT Mock",
        str_detect(file_name, "MKRN1") ~ "MKRN1KO HIV"
    )) %>%
    mutate(condition = fct_relevel(
        condition,
        c("WT Mock", "WT HIV", "MKRN1KO HIV")
    )) %>%
    mutate(repeat_code = str_sub(file_name, 1, 6)) %>%
    mutate(repeats = case_when(
        repeat_code == "221123" ~ "Rep1",
        repeat_code == "221218" ~ "Rep2",
        repeat_code == "221221" ~ "Rep3"
    ))

## Adjust replicate variability
rep_adj_df <- df_discretised %>%
    group_by(channel, time, condition, repeats) %>%
    summarise(
        mean_intensity = mean(gray_value_by_discrete_distance),
        median_intensity = median(gray_value_by_discrete_distance)
    ) %>%
    ungroup() %>%
    group_by(channel, time, condition) %>%
    mutate(group_mean = mean(mean_intensity)) %>%
    mutate(adj_factor = group_mean / mean_intensity) %>%
    ungroup()

df_discretised_adjusted <- df_discretised %>%
    left_join(rep_adj_df, by = c("time", "condition", "channel", "repeats")) %>%
    mutate(gray_value_by_discrete_distance_adj = gray_value_by_discrete_distance * adj_factor)

## Plot
df_discretised_adjusted <- read_csv("./gag_line_profile/df_discretised_adjusted.csv")
breaks_to_plot <- c(TRUE, rep(FALSE, 19))

### For gag protein
df_discretised_adjusted %>%
    filter(channel == "Gag") %>%
    filter(condition != "WT Mock") %>%
    filter(!is.na(discrete_distance)) %>%
    ggplot(aes(x = as.factor(discrete_distance), y = gray_value_by_discrete_distance_adj, colour = condition, group = condition)) +
    geom_line(stat = "summary", alpha = 0.2, size = 0.3) + 
    geom_pointrange(stat = "summary", size = 0.1, alpha = 0.3) +
    facet_wrap(~time, nrow = 3) +
    labs(
        x = "Cell cross-secion distance (AU)",
        y = "Normalised Gag fluorescence intensity (AU)",
        colour = "Condition"
        ) +
    scale_colour_manual(values = c("dodgerblue", "coral")) +
    scale_x_discrete(breaks = function(x){x[breaks_to_plot]}) +
    theme_classic()

ggsave("./gag_line_profile/gag_line_profile_plot.pdf", height = 8, width = 10)

### For HIV RNA
df_discretised_adjusted %>%
    filter(channel == "RNA") %>%
    filter(condition != "WT Mock") %>%
    filter(!is.na(discrete_distance)) %>%
    ggplot(aes(x = as.factor(discrete_distance), y = gray_value_by_discrete_distance_adj, colour = condition, group = condition)) +
    geom_line(stat = "summary", alpha = 0.2, size = 0.3) + 
    geom_pointrange(stat = "summary", size = 0.1, alpha = 0.3) +
    facet_wrap(~time, nrow = 3) +
    labs(
        x = "Cell cross-secion distance (AU)",
        y = "Normalised Gag fluorescence intensity (AU)",
        colour = "Condition"
        ) +
    scale_colour_manual(values = c("dodgerblue", "coral")) +
    scale_x_discrete(breaks = function(x){x[breaks_to_plot]}) +
    theme_classic()

ggsave("./gag_line_profile/rna_line_profile_plot.pdf", height = 8, width = 10)

## Save df output
write_csv(df_discretised_adjusted, "./gag_line_profile/df_discretised_adjusted.csv")

## Statistics
### Testing ratio of Area under curve between inner and outer regions
library(rstatix)

df_ratio <- df_discretised_adjusted %>%
    filter(!is.na(discrete_distance)) %>%
    # mutate(localisation = if_else(discrete_distance > 25 | discrete_distance < 75, "inner", "outer")) %>%
    group_by(channel, condition, time, discrete_distance, repeats) %>%
    summarise(mean_norm_intensity = mean(norm_intensity_by_discrete_distance)) %>%
    ungroup() %>%
    mutate(localisation = if_else(discrete_distance > 25 & discrete_distance < 75, "inner", "outer")) %>%
    group_by(channel, condition, time, repeats, localisation) %>%
    summarise(sum_aoc = sum(mean_norm_intensity)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(channel, condition, time, repeats), names_from = localisation, values_from = sum_aoc) %>%
    group_by(channel, condition, time, repeats) %>%
    summarise(ratio = outer / inner)

df_ratio %>%
    ggplot(aes(x = condition, y = ratio)) +
    geom_point() +
    theme_classic() +
    facet_wrap(~channel + time)

df_ratio %>%
    filter(condition != "WT Mock") %>%
    group_by(channel, time) %>%
    t_test(ratio ~ condition)
