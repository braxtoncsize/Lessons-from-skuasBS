##############################################################################
# SKUA CASE STUDY -- CONSIDERATIONS FOR DISEASE/MOVEMENT ECOLOGY STUDY #
# Falkland Islands/Islas Malvinas Bayesian Multilevel Regression Analyses #
# Last updated: 22Mar26 BS #
# R version 4.5.2 (2025-10-31) # 

##############################################################################

library(brms)
library(bayesplot)
library(shinystan)
library(ggplot2)
library(readr)
library(dplyr)
library(tidytable)

FLK23_MasterSheet <- read_csv("Mastersheets/merged_3day_windows_multi_UD.csv")
master <- FLK23_MasterSheet

## add capture type
master$breeding_status <- NA
allids <- unique(master$ID)

onnest <- c("912", "916", "920", "927", "928", "929", "933", "934", "935", "941", "952", "954", "958")
offnest <- setdiff(allids, onnest)

for(i in 1:nrow(master)) {
  if(master$ID[i] %in% onnest) {
    master$breeding_status[i] <- "On-nest"
  }
  else if (master$ID[i] == 915) {
    master$breeding_status[i] <- "Off-nest?"
  }
  else {
    master$breeding_status[i] <- "Off-nest"
  }
}


## add island, add as column
new <- c("930", "932", "933", "934", "935", "941", "942", "944", "956", "959")
saunders <- c("911", "912", "914", "918", "927", "928", "929")
bleaker <- setdiff(allids, union(new, saunders))

master <- master %>%
  mutate(
    island = case_when(
      ID %in% new ~ "New Island",
      ID %in% saunders ~ "Saunders",
      TRUE ~ "Bleaker"
    ))  

## make factors for analyses
master$breeding_status <- as.factor(master$breeding_status)
master$island <- as.factor(master$island)
master$fluaca <- factor(master$fluaca, levels = c("Negative", "Positive"))

# Ensure date is in Date format
master$date <- as.Date(master$window_start)

# Find first date in dataset
season_start <- min(master$date)
master <- master %>% select(-center_date)

# Compute day of season
master$day_of_season <- as.numeric(master$date - season_start) + 1
# Day 1 = first observation date, Day 2 = second day, etc.

# Create log-transformed response variable
master$log_distance <- log(master$total_distance_km)
master$log_area_50 <- log(master$area_50_overall_km2)
master$log_area_95 <- log(master$area_95_overall_km2)
master$log_area_state_1_50 <- log(master$area_50_state_1_km2)
master$log_area_state_1_95 <- log(master$area_95_state_1_km2)
master$log_area_state_2_50 <- log(master$area_50_state_2_km2)
master$log_area_state_2_95 <- log(master$area_95_state_2_km2)
master$log_area_state_3_50 <- log(master$area_50_state_3_km2)
master$log_area_state_3_95 <- log(master$area_95_state_3_km2)

############################################################
# log_distance model
bayes_distance <- brm(bf(log_distance ~ island * breeding_status + day_of_season + (1 | ID),
                         sigma ~ island * breeding_status),
                      data = master,
                      family = gaussian(),
                      prior = c(
                        prior(normal(0, 2), class = b),
                        prior(normal(0, 2), class = b, dpar = sigma)
                      ),
                      chains = 4, iter = 2000, cores = 4, seed = 123,
                      save_pars = save_pars(all = TRUE))

launch_shinystan(bayes_distance) 
  # no diverging iterations
  # traceplot shows good chain mixing
  # no parameters with Monte Carlo standard error > 10% of posterior SD
  # no parameters with Rhat above 1.1

# log_area_50 model
bayes_area_50 <- brm(bf(log_area_50 ~ island * breeding_status + day_of_season + (1 | ID),
                        sigma ~ island * breeding_status),
                     data = master,
                     family = gaussian(),
                     prior = c(
                       prior(normal(0, 2), class = b),
                       prior(normal(0, 2), class = b, dpar = sigma)),
                     chains = 4, iter = 2000, cores = 4, seed = 123,
                     save_pars = save_pars(all = TRUE))

launch_shinystan(bayes_area_50) 

# log_area_95 model
bayes_area_95 <- brm(bf(log_area_95 ~ island * breeding_status + day_of_season + (1 | ID),
                        sigma ~ island * breeding_status),
                     data = master,
                     family = gaussian(),
                     prior = c(
                       prior(normal(0, 2), class = b),
                       prior(normal(0, 2), class = b, dpar = sigma)),
                     chains = 4, iter = 2000, cores = 4, seed = 123,
                     save_pars = save_pars(all = TRUE))

launch_shinystan(bayes_area_95) 

# log_area_state_1_50 model
bayes_area_state_1_50 <- brm(bf(log_area_state_1_50 ~ island * breeding_status + day_of_season + (1 | ID),
                                sigma ~ island * breeding_status),
                             data = master,
                             family = gaussian(),
                             prior = c(
                               prior(normal(0, 2), class = b),
                               prior(normal(0, 2), class = b, dpar = sigma)),
                             chains = 4, iter = 2000, cores = 4, seed = 123,
                             save_pars = save_pars(all = TRUE))

launch_shinystan(bayes_area_state_1_50) 

# log_area_state_1_95 model
bayes_area_state_1_95 <- brm(bf(log_area_state_1_95 ~ island * breeding_status + day_of_season + (1 | ID),
                                sigma ~ island * breeding_status),
                             data = master,
                             family = gaussian(),
                             prior = c(
                               prior(normal(0, 2), class = b),
                               prior(normal(0, 2), class = b, dpar = sigma)),
                             chains = 4, iter = 2000, cores = 4, seed = 123,
                             save_pars = save_pars(all = TRUE))

launch_shinystan(bayes_area_state_1_95) 

# log_area_state_2_50 model
bayes_area_state_2_50 <- brm(bf(log_area_state_2_50 ~ island * breeding_status + day_of_season + (1 | ID),
                                sigma ~ island * breeding_status),
                             data = master,
                             family = gaussian(),
                             prior = c(
                               prior(normal(0, 2), class = b),
                               prior(normal(0, 2), class = b, dpar = sigma)),
                             chains = 4, iter = 2000, cores = 4, seed = 123,
                             save_pars = save_pars(all = TRUE))

launch_shinystan(bayes_area_state_2_50) 
  # on parameter with effective sample size less than 10% of total samples size, z_2[1,8], z_2[1,9]

# log_area_state_2_95 model
bayes_area_state_2_95 <- brm(bf(log_area_state_2_95 ~ island * breeding_status + day_of_season + (1 | ID),
                                sigma ~ island * breeding_status),
                             data = master,
                             family = gaussian(),
                             prior = c(
                               prior(normal(0, 2), class = b),
                               prior(normal(0, 2), class = b, dpar = sigma)),
                             chains = 4, iter = 2000, cores = 4, seed = 123,
                             save_pars = save_pars(all = TRUE))

launch_shinystan(bayes_area_state_2_95) 

# log_area_state_3_50 model
bayes_area_state_3_50 <- brm(bf(log_area_state_3_50 ~ island * breeding_status + day_of_season + (1 | ID),
                                sigma ~ island * breeding_status),
                             data = master,
                             family = gaussian(),
                             prior = c(
                               prior(normal(0, 2), class = b),
                               prior(normal(0, 2), class = b, dpar = sigma)),
                             chains = 4, iter = 2000, cores = 4, seed = 123,
                             save_pars = save_pars(all = TRUE))

launch_shinystan(bayes_area_state_3_50) 

# log_area_state_3_95 model
bayes_area_state_3_95 <- brm(bf(log_area_state_3_95 ~ island * breeding_status + day_of_season + (1 | ID),
                                sigma ~ island * breeding_status),
                             data = master,
                             family = gaussian(),
                             prior = c(
                               prior(normal(0, 2), class = b),
                               prior(normal(0, 2), class = b, dpar = sigma)),
                             chains = 4, iter = 2000, cores = 4, seed = 123,
                             save_pars = save_pars(all = TRUE))

launch_shinystan(bayes_area_state_3_95) 

######################
# figures...
library(tidyverse)
library(viridis)
library(brms)

# clean term labels 
clean_term <- function(term) {
  case_when(
    str_detect(term, "^sigma_Intercept$") ~ "Intercept",
    term == "Intercept" ~ "Intercept",
    str_detect(term, "^sigma_") ~ {
      clean_sigma <- str_remove(term, "^sigma_")
      case_when(
        clean_sigma == "islandNewIsland" ~ "New Island",
        clean_sigma == "islandSaunders" ~ "Saunders", 
        clean_sigma == "breeding_statusOnMnest" ~ "On-nest",
        clean_sigma == "islandNewIsland:breeding_statusOnMnest" ~ "New Island × On-nest",
        clean_sigma == "islandSaunders:breeding_statusOnMnest" ~ "Saunders × On-nest",
        TRUE ~ str_replace_all(clean_sigma, c("_" = " ", ":" = " × "))
      )
    },
    term == "islandNewIsland" ~ "New Island",
    term == "islandSaunders" ~ "Saunders", 
    term == "breeding_statusOnMnest" ~ "On-nest",
    term == "day_of_season" ~ "Day of Season",
    term == "islandNewIsland:breeding_statusOnMnest" ~ "New Island × On-nest",
    term == "islandSaunders:breeding_statusOnMnest" ~ "Saunders × On-nest",
    TRUE ~ str_replace_all(term, c("_" = " ", ":" = " × "))
  )
}

term_order <- c(
  "Saunders × On-nest",
  "New Island × On-nest", 
  "On-nest",
  "Saunders",
  "New Island",
  "Day of Season",    # Add this line
  "Intercept"
)

# clean model names for plotting
clean_model_name <- function(name) {
  name %>%
    str_replace_all("_", " ") %>%
    str_replace("state 1", "Resting") %>%
    str_replace("state 2", "Foraging") %>%
    str_replace("state 3", "Commuting") %>%
    str_replace("log area ", "Log ") %>%
    str_replace("log distance", "Log Distance") %>%
    str_replace(" 50", " 50% UD Area") %>%
    str_replace(" 95", " 95% UD Area") %>%
    str_to_title() %>%
    str_replace_all("Ud", "UD")              # Ensure UD stays uppercase
}

models <- list(
  "log_area_50"          = bayes_area_50,
  "log_area_95"          = bayes_area_95,
  "log_area_state_1_50"  = bayes_area_state_1_50,
  "log_area_state_1_95"  = bayes_area_state_1_95,
  "log_area_state_2_50"  = bayes_area_state_2_50,
  "log_area_state_2_95"  = bayes_area_state_2_95,
  "log_area_state_3_50"  = bayes_area_state_3_50,
  "log_area_state_3_95"  = bayes_area_state_3_95,
  "log_distance"         = bayes_distance
)

# extract & clean all fixed effects
effects_df <- map2_df(names(models), models, function(model_name, model) {
  fe <- as.data.frame(fixef(model, robust = TRUE))
  fe %>%
    rownames_to_column("term") %>%
    mutate(
      model      = clean_model_name(model_name),
      submodel   = if_else(str_detect(term, "^sigma_"), "Sigma sub-model", "Mean sub-model"),
      term_clean = clean_term(term)
    )
})

# Add this after the effects_df creation and before the grouping variables
effects_df <- effects_df %>%
  mutate(
    term_clean = factor(term_clean, levels = term_order),
    # Create ordered factor for submodel to control dodge order (Mean first)
    submodel_ordered = factor(submodel, levels = c("Mean sub-model", "Sigma sub-model"))
  )

effects_df <- effects_df %>%
  mutate(term_clean = factor(term_clean, levels = term_order))

# add grouping variables for facet layout
effects_df <- effects_df %>%
  mutate(
    group = case_when(
      str_detect(model, "50% UD Area|95% UD Area") & 
        !str_detect(model, "Resting|Foraging|Commuting") ~ "Whole",
      str_detect(model, "Resting") ~ "Resting",
      str_detect(model, "Foraging") ~ "Foraging",
      str_detect(model, "Commuting") ~ "Commuting",
      str_detect(model, "Distance") ~ "Distance"
    ),
    version = case_when(
      str_detect(model, "50% UD Area") ~ "50% UD Area",
      str_detect(model, "95% UD Area") ~ "95% UD Area",
      str_detect(model, "Distance") ~ "Distance"
    )
  )

# split datasets
ud_df   <- effects_df %>% filter(group != "Distance")
dist_df <- effects_df %>% filter(group == "Distance")

# Filter out Day of Season for main plots
ud_df_no_day <- ud_df %>% filter(term_clean != "Day of Season")
dist_df_no_day <- dist_df %>% filter(term_clean != "Day of Season")

# Create Day of Season only datasets
ud_df_day <- ud_df %>% filter(term_clean == "Day of Season")
dist_df_day <- dist_df %>% filter(term_clean == "Day of Season")

# Main plots without Day of Season
p_ud <- ggplot(ud_df_no_day,
               aes(y = term_clean, x = Estimate, xmin = Q2.5, xmax = Q97.5, color = submodel_ordered)) +
  geom_vline(xintercept = 0, color = "darkgrey", linewidth = 0.8, alpha = 0.8) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_errorbarh(height = 0.3, position = position_dodge(width = 0.6)) +
  facet_grid(rows = vars(group), cols = vars(version), scales = "free_y") +
  scale_color_manual(values = c("Mean sub-model" = "#30123B", "Sigma sub-model" = "#1AE4B6")) +
  theme_bw() +
  theme(
    strip.text      = element_text(face = "bold", size = 10),
    legend.position = "top",
    axis.text.y     = element_text(size = 9),
    panel.spacing.x = unit(1, "lines")
  ) +
  labs(
    x = "Effect Size",
    y = NULL,
    color = "Sub-model"
  )

p_dist <- ggplot(dist_df_no_day,
                 aes(y = term_clean, x = Estimate, xmin = Q2.5, xmax = Q97.5, color = submodel_ordered)) +
  geom_vline(xintercept = 0, color = "darkgrey", linewidth = 0.8, alpha = 0.8) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_errorbarh(height = 0.3, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("Mean sub-model" = "#30123B", "Sigma sub-model" = "#1AE4B6")) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.y     = element_text(size = 9)
  ) +
  labs(
    x = "Effect Size",
    y = NULL,
    color = "Sub-model"
  )

# Combined Day of Season plot (UD + Distance models, mean model only)
day_df_combined <- bind_rows(ud_df_day, dist_df_day) %>%
  filter(submodel == "Mean sub-model")  # Only keep mean sub-model

p_day_combined <- ggplot(day_df_combined,
                         aes(y = term_clean, x = Estimate, xmin = Q2.5, xmax = Q97.5, color = model)) +
  geom_vline(xintercept = 0, color = "darkgrey", linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 2, position = position_dodge(width = 0.4)) +
  geom_errorbarh(height = 0.3, position = position_dodge(width = 0.4)) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y     = element_text(size = 9),
    legend.text     = element_text(size = 8)
  ) +
  labs(
    x = "Effect Size",
    y = NULL,
    color = "Model",
    title = "Day of Season Effects"
  )

# Save all plots
ggsave("fixed_effects_ud_models.png", p_ud, width = 8, height = 6, dpi = 600)
ggsave("fixed_effects_distance.png", p_dist, width = 6, height = 5, dpi = 600)
ggsave("day_of_season_effects.png", p_day_combined, width = 8, height = 6, dpi = 600)

#### ridge plots
library(ggridges)

# Extract posterior samples function (same as before)
extract_posterior_samples <- function(model_name, model) {
  posterior_samples <- as_draws_df(model)
  
  # Get fixed effect columns
  fe_cols <- posterior_samples %>%
    select(starts_with("b_")) %>%
    names()
  
  posterior_samples %>%
    select(all_of(fe_cols), .draw) %>%
    pivot_longer(cols = -.draw, names_to = "term", values_to = "value") %>%
    mutate(
      model = clean_model_name(model_name),
      term = str_remove(term, "^b_"),
      submodel = if_else(str_detect(term, "^sigma_"), "Sigma sub-model", "Mean sub-model"),
      term_clean = clean_term(term)
    )
}

# Extract all posterior samples
posterior_df <- map2_df(names(models), models, extract_posterior_samples)

# Add the same grouping and ordering as before
posterior_df <- posterior_df %>%
  mutate(
    term_clean = factor(term_clean, levels = term_order),
    submodel_ordered = factor(submodel, levels = c("Mean sub-model", "Sigma sub-model")),
    group = case_when(
      str_detect(model, "50% UD Area|95% UD Area") & 
        !str_detect(model, "Resting|Foraging|Commuting") ~ "Whole",
      str_detect(model, "Resting") ~ "Resting",
      str_detect(model, "Foraging") ~ "Foraging",
      str_detect(model, "Commuting") ~ "Commuting",
      str_detect(model, "Distance") ~ "Distance"
    ),
    version = case_when(
      str_detect(model, "50% UD Area") ~ "50% UD Area",
      str_detect(model, "95% UD Area") ~ "95% UD Area",
      str_detect(model, "Distance") ~ "Distance"
    )
  )

# Filter datasets (excluding Day of Season from main plots)
ud_posterior <- posterior_df %>% filter(group != "Distance", term_clean != "Day of Season")
dist_posterior <- posterior_df %>% filter(group == "Distance", term_clean != "Day of Season")

# Day of Season data (mean sub-model only)
day_posterior <- posterior_df %>% 
  filter(term_clean == "Day of Season", submodel == "Mean sub-model")

# Ridge plot for UD models
p_ud_ridges <- ggplot(ud_posterior, 
                      aes(x = value, y = term_clean, 
                          fill = submodel_ordered, color = submodel_ordered)) +
  geom_density_ridges(alpha = 0.7, scale = 0.9, size = 0.5, 
                      position = position_points_jitter(width = 0.05, height = 0)) +
  geom_vline(xintercept = 0, color = "darkgrey", linewidth = 0.8, alpha = 0.8) +
  facet_grid(rows = vars(group), cols = vars(version), scales = "free") +
  scale_fill_manual(values = c("Mean sub-model" = "#30123B", "Sigma sub-model" = "#1AE4B6")) +
  scale_color_manual(values = c("Mean sub-model" = "#30123B", "Sigma sub-model" = "#1AE4B6")) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "top",
    axis.text.y = element_text(size = 9),
    panel.spacing.x = unit(1, "lines")
  ) +
  labs(
    x = "Effect Size",
    y = NULL,
    fill = "Sub-model",
    color = "Sub-model"
  )

# Ridge plot for Distance model
p_dist_ridges <- ggplot(dist_posterior, 
                        aes(x = value, y = term_clean, 
                            fill = submodel_ordered, color = submodel_ordered)) +
  geom_density_ridges(alpha = 0.7, scale = 0.9, size = 0.5) +
  geom_vline(xintercept = 0, color = "darkgrey", linewidth = 0.8, alpha = 0.8) +
  scale_fill_manual(values = c("Mean sub-model" = "#30123B", "Sigma sub-model" = "#1AE4B6")) +
  scale_color_manual(values = c("Mean sub-model" = "#30123B", "Sigma sub-model" = "#1AE4B6")) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 9)
  ) +
  labs(
    x = "Effect Size",
    y = NULL,
    fill = "Sub-model",
    color = "Sub-model"
  )

# Ridge plot for Day of Season effects
p_day_ridges <- ggplot(day_posterior, 
                       aes(x = value, y = model, fill = model)) +
  geom_density_ridges(alpha = 0.8, scale = 0.9, size = 0.5) +
  geom_vline(xintercept = 0, color = "darkgrey", linewidth = 0.8, alpha = 0.8) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  theme_bw() +
  theme(
    legend.position = "none",  # Too many models for a clean legend
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    x = "Effect Size",
    y = NULL
  )

# Save all ridge plots
ggsave("fixed_effects_ud_ridges.png", p_ud_ridges, width = 6, height = 6, dpi = 600)
ggsave("fixed_effects_distance_ridges.png", p_dist_ridges, width = 7, height = 5, dpi = 600)
ggsave("day_of_season_ridges.png", p_day_ridges, width = 6, height = 5, dpi = 600)

