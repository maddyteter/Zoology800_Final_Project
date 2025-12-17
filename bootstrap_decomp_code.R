#Bootstrapping decomp data
library(dplyr)
library(minpack.lm)
library(ggplot2)
library(purrr)

#Set working directory 
setwd("C:/Users/maddy/OneDrive/Grad Classes/Zoology 800/Zoology800_Final_Project")


#PSD25 Decomposition values
decomp <- read.csv("PSD25_decomposition.csv", header = TRUE)

# Read in site data to group by other variables in site dataframe
site_data <- read.csv("PSD25_site_data.csv")

#Create column for remaining tea only - no bag weight included (0.21g)
decomp <- decomp %>%
  mutate(final_tea_dw_g = final_tb_dw_g - empty_bag_weight_g)

# calculate se and mean for tea_decomp
se <- function(x) sd(x, na.rm = TRUE) / sqrt(length((x)))

#Join in site data to decomp df
decomp <- decomp %>%
  left_join(site_data %>% select(site_id, treatment, 
                                 distance_from_scar_m, distance_from_edge_m), 
            by = "site_id")

#create column for mass fraction/proportion at each time point
decomp <- decomp %>%
  mutate(tea_dw_proportion = final_tea_dw_g/initial_tea_dw_g)

#Find decomposable fraction using t7 mass (where the curve starts to 
#level out in all plots). This fraction represents 
#the portion of tea that is labile and easily decomposable

decomp_a_value <- decomp %>%
  filter(time_interval_wks == "t7") %>%
  mutate(a = (1 - tea_dw_proportion))

#Create average a value for each tea type to use in equation
decomp_a_value <- decomp_a_value %>%
  group_by(tea_type) %>%
  summarise(a_average = mean(a, na.rm = TRUE), .groups = "drop")

#join average a values back into decomp df
decomp <- decomp %>%
  left_join(decomp_a_value %>% select(tea_type, a_average), 
            by = c("tea_type"))

#Visualize distribution of a values by tea type 
decomp_a_value %>%
  ggplot(aes(tea_type, a_average)) +
  geom_boxplot() +
  theme_classic()
#A values are very tight by tea type, not a lot of variation. For now use a single
#a value in formulas per tea type (might change at 6 month point)

#Filter for green tea data only 
green_model_data <- decomp %>%
  filter(tea_type == "G")

#Remove single row with tea bag note that says "do not use" due to tea 
#loss during cleaning
green_model_data <- subset(green_model_data, sample_ID != "PSD25_S1.2_t7_B_G")


#Create function to fit k to green tea
fit_k_once <- function(df) {
  df <- df %>% filter(!is.na(tea_dw_proportion))
  if (nrow(df) < 3) return(NA_real_) #requires at least 3 mass loss values to calculate k 
  a_val <- unique(df$a_average) #pulling labile portion (a) for green tea 
  if (length(a_val) != 1) return(NA_real_)
  fit <- tryCatch( #fit k using mass loss values at each timepoint plus a to account for labile proportion
    nlsLM(
      tea_dw_proportion ~ a_val * exp(-k * total_deployment_time_days) + (1 - a_val),
      data = df,
      start = list(k = 0.1),
      lower = 0,
      control = nls.lm.control(maxiter = 200)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_) #return na rather than crashing if fit can't be calculated
  coef(fit)[["k"]]
}

#Run bootstrap for various green tea proportion dw grouped by site*site_position 
#to get a range of kg values
set.seed(123)
n_boot <- 1000

green_k_boot_refit <- green_model_data %>%
  group_by(site, site_position) %>%
  group_modify(~ {
    df <- .x
    n <- nrow(df)
    tibble(
      k_boot = replicate(n_boot, {
        df_boot <- df[sample(seq_len(n), replace = TRUE), ]
        fit_k_once(df_boot)
      })
    )
  }) %>%
  ungroup()


#group kg values by site*site position and find kg mean, sd, and CI's
green_k_boot_summary <- green_k_boot_refit %>%
  group_by(site, site_position) %>%
  summarise(
    k_mean = mean(k_boot, na.rm = TRUE),
    k_sd   = sd(k_boot, na.rm = TRUE),
    k_lwr  = quantile(k_boot, 0.025, na.rm = TRUE),
    k_upr  = quantile(k_boot, 0.975, na.rm = TRUE),
    n_fail = sum(is.na(k_boot)),
    .groups = "drop"
  )
green_k_boot_summary

#Join in site data
green_k_boot_summary <- green_k_boot_summary %>%
  left_join(site_data %>% select(site, site_position, treatment, distance_from_edge_m,
                                 distance_from_scar_m),
            by = c("site", "site_position"))

green_k_boot_refit <- green_k_boot_refit %>%
  left_join(site_data %>% select(site, site_position, treatment, distance_from_edge_m,
                                 distance_from_scar_m),
            by = c("site", "site_position"))

#Make treatment and distance factors to compare across groups
green_k_boot_summary <- green_k_boot_summary %>%
  mutate(treatment = as.factor(treatment),
         distance_from_edge_m = as.factor(distance_from_edge_m),
         distance_from_scar_m = as.factor(distance_from_scar_m))

green_k_boot_refit <- green_k_boot_refit %>%
  mutate(treatment = as.factor(treatment),
         distance_from_edge_m = as.factor(distance_from_edge_m),
         distance_from_scar_m = as.factor(distance_from_scar_m))

#Plot green k value CI's and mean across all sites and positions
ggplot(green_k_boot_summary,
       aes(x = site_position, y = k_mean,
           ymin = k_lwr, ymax = k_upr, color = treatment)) +
  geom_pointrange() +
  facet_wrap(~ site) +
  theme_classic() +
  labs(
    x = "Site position",
    y = "Decay rate (k)",
    title = "Green Tea Bootstrap-refit k with 95% CI"
  )

#Create df to compare edge vs interior for kg (using bootstrapped k values, 
#not summary stats)
green_k_boot_edge <- green_k_boot_refit %>%
  filter(distance_from_scar_m == c("0", "10"))


delta_kg_edge <- green_k_boot_edge %>%
  group_by(site, treatment) %>%  # only group by site
  group_modify(~ {
    vals_edge     <- .x %>% filter(distance_from_edge_m == "0") %>% pull(k_boot)
    vals_interior <- .x %>% filter(distance_from_edge_m == "30") %>% pull(k_boot)
    
    delta <- outer(vals_edge, vals_interior, `-`) %>% as.vector()
    
    tibble(
      median_delta = median(delta, na.rm = TRUE),
      lwr = quantile(delta, 0.025, na.rm = TRUE),
      upr = quantile(delta, 0.975, na.rm = TRUE),
      prob_edge_higher = mean(delta > 0, na.rm = TRUE)
    )
  }) %>%
  ungroup()
delta_kg_edge


#create summary dataset to view k value mean and distributions by group
green_k_boot_edge_summary <- green_k_boot_summary %>%
  filter(distance_from_scar_m == c("0", "10"))

ggplot(green_k_boot_edge_summary,
       aes(x = distance_from_edge_m, y = k_mean,
           ymin = k_lwr, ymax = k_upr, color = treatment)) +
  geom_pointrange(position = position_dodge(width = 1.5)) +
  facet_wrap(~ site) +
  scale_color_manual(
    values = c(scarred = "black", vegetated = "seagreen3")) +
  theme_classic() +
  labs(
    x = "Distance from Edge (m)",
    y = "Green Tea Decay rate (kg)",
    title = "Kg vs Distance from Meadow Edge by Treatment (Bootstrap w 95% CI)",
    color = "Treatment"
  )

#Create df to compare green tea scar vs distance for kg (using bootstrapped k values, 
#not summary stats)
green_k_boot_scar <- green_k_boot_refit %>%
  filter(distance_from_edge_m == "30")

#Calculate the change in k, grouped by site and distance from scar 
delta_kg_scar <- green_k_boot_scar %>%
  filter(distance_from_scar_m != 0) %>% #filter out 0m (scar) since this is the values we're comparing kr to
  group_by(site, distance_from_scar_m) %>% #give output for each site*distance grouping
  group_modify(~ {
    scar_vals <- green_k_boot_scar %>%
      filter(site == .y$site,
             distance_from_scar_m == 0) %>% 
      pull(k_boot) #pull scar kr values from every site - pulled from kr bootstrap values
    tibble(
      delta_kg = scar_vals - .x$k_boot #difference between scar kr bootstrap values and values at all other sites and distances
    )
  }) %>%
  ungroup()

#create summary table of kr values comparing distance from scar
delta_kg_scar_summary <- delta_kg_scar %>%
  group_by(site, distance_from_scar_m) %>%
  summarise(
    median_delta = median(delta_kg, na.rm = TRUE),
    lwr = quantile(delta_kg, 0.025, na.rm = TRUE),
    upr = quantile(delta_kg, 0.975, na.rm = TRUE),
    prob_scar_higher = mean(delta_kg > 0, na.rm = TRUE),
    .groups = "drop"
  )

delta_kg_scar_summary
#Very inconsistent - sometimes the scar is higher, sometimes the vegetation at 
#various distances is higher. No clear trend in distance or treatment across sites

#create summary df to view scar vs distances adjacent to scar 
green_k_boot_scar_summary <- green_k_boot_summary %>%
  filter(distance_from_edge_m == "30")


ggplot(green_k_boot_scar_summary,
       aes(x = factor(distance_from_scar_m,
                      levels = c(0, 0.25, 1, 10)),
           y = k_mean,
           ymin = k_lwr, ymax = k_upr,
           color = treatment)) +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  scale_color_manual(
    values = c(scarred = "black", vegetated = "seagreen3")
  ) +
  facet_wrap(~ site) +
  theme_classic() +
  labs(
    x = "Distance from Scar (m)",
    y = "Green Tea Decay rate (kg)",
    title = "Kg vs Distance from Scar by Treatment (Bootstrap w 95% CI)"
  )



####Repeat bootstrapping process with red tea
red_model_data <- decomp %>%
  filter(tea_type == "R")

#Remove rows that have holes/tears causing abnormal weights (s2.3 had multiple bags
#that was causing a range of k values that were extreme outliers compared to all other
#red tea)
red_model_data <- red_model_data %>%
  filter(!sample_ID %in% c("PSD25_S2.3_t7_C_R", "PSD25_S2.3_t2_A_R",
                           "PSD25_S2.3_t4_B_R"))

#Create function to fit k to green tea
fit_k_once <- function(df) {
  df <- df %>% filter(!is.na(tea_dw_proportion))
  if (nrow(df) < 3) return(NA_real_) #requires at least 3 mass loss values at diff time points
  a_val <- unique(df$a_average) #pulling labile portion (a) for red tea 
  if (length(a_val) != 1) return(NA_real_)
  fit <- tryCatch( #fit k using mass loss values at each timepoint plus a to account for labile proportion
    nlsLM(
      tea_dw_proportion ~ a_val * exp(-k * total_deployment_time_days) + (1 - a_val),
      data = df,
      start = list(k = 0.1),
      lower = 0,
      control = nls.lm.control(maxiter = 200)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  coef(fit)[["k"]]
}

#Run bootstrap for every k grouped by site*site_position to get a range of k values
set.seed(123)
n_boot <- 1000

red_k_boot_refit <- red_model_data %>%
  group_by(site, site_position) %>%
  group_modify(~ {
    df <- .x
    n <- nrow(df)
    tibble(
      k_boot = replicate(n_boot, {
        df_boot <- df[sample(seq_len(n), replace = TRUE), ]
        fit_k_once(df_boot)
      })
    )
  }) %>%
  ungroup()


#group k values by site*site position and find kg mean, sd, and CI's
red_k_boot_summary <- red_k_boot_refit %>%
  group_by(site, site_position) %>%
  summarise(
    k_mean = mean(k_boot, na.rm = TRUE),
    k_sd   = sd(k_boot, na.rm = TRUE),
    k_lwr  = quantile(k_boot, 0.025, na.rm = TRUE),
    k_upr  = quantile(k_boot, 0.975, na.rm = TRUE),
    n_fail = sum(is.na(k_boot)),
    .groups = "drop"
  )

#Join in site data
red_k_boot_summary <- red_k_boot_summary %>%
  left_join(site_data %>% select(site, site_position, treatment, distance_from_edge_m,
                                 distance_from_scar_m),
            by = c("site", "site_position"))

red_k_boot_refit <- red_k_boot_refit %>%
  left_join(site_data %>% select(site, site_position, treatment, distance_from_edge_m,
                                 distance_from_scar_m),
            by = c("site", "site_position"))

#Make treatment and distances factors (groups of interest for comparison)
red_k_boot_summary <- red_k_boot_summary %>%
  mutate(treatment = as.factor(treatment),
         distance_from_edge_m = as.factor(distance_from_edge_m),
         distance_from_scar_m = as.factor(distance_from_scar_m))

red_k_boot_refit <- red_k_boot_refit %>%
  mutate(treatment = as.factor(treatment),
         distance_from_edge_m = as.factor(distance_from_edge_m),
         distance_from_scar_m = as.factor(distance_from_scar_m))

#Plot red k value CI's and mean for all sites and positions (summary)
ggplot(red_k_boot_summary,
       aes(x = site_position, y = k_mean,
           ymin = k_lwr, ymax = k_upr, color = treatment)) +
  geom_pointrange() +
  facet_wrap(~ site) +
  theme_classic() +
  labs(
    x = "Site position",
    y = "Decay rate (k)",
    title = "Red Tea Bootstrap-refit k with 95% CI"
  )

#Create df to compare stats for kr between edge and interior meadow
red_k_boot_edge <- red_k_boot_refit %>%
  filter(distance_from_scar_m == c("0", "10"))

delta_kr_edge <- red_k_boot_edge %>%
  group_by(site, treatment) %>%  # only group by site
  group_modify(~ {
    vals_edge     <- .x %>% filter(distance_from_edge_m == "0") %>% pull(k_boot)
    vals_interior <- .x %>% filter(distance_from_edge_m == "30") %>% pull(k_boot)
    
    delta <- outer(vals_edge, vals_interior, `-`) %>% as.vector()
    
    tibble(
      median_delta = median(delta, na.rm = TRUE),
      lwr = quantile(delta, 0.025, na.rm = TRUE),
      upr = quantile(delta, 0.975, na.rm = TRUE),
      prob_edge_higher = mean(delta > 0, na.rm = TRUE)
    )
  }) %>%
  ungroup()
delta_kr_edge


#Create df to visualize differences in kr between edge and interior meadow
red_k_boot_edge_summary <- red_k_boot_summary %>%
  filter(distance_from_scar_m == c("0", "10"))

ggplot(red_k_boot_edge_summary,
       aes(x = distance_from_edge_m, y = k_mean,
           ymin = k_lwr, ymax = k_upr, color = treatment)) +
  geom_pointrange(position = position_dodge(width = 1.5)) +
  scale_color_manual(
    values = c(scarred = "black", vegetated = "firebrick") 
  ) +
  facet_wrap(~ site) +
  theme_classic() +
  labs(
    x = "Distance from Edge (m)",
    y = "Red Tea Decay rate (kr)",
    title = "Kr vs Distance from Meadow Edge by Treatment (Bootstrap w 95% CI)"
  )

#create df to statistically compare kr between scar and meadow 
red_k_boot_scar <- red_k_boot_refit %>%
  filter(distance_from_edge_m == "30")

#Calculate the change in k, grouped by site and distance from scar 
delta_kr_scar <- red_k_boot_scar %>%
  filter(distance_from_scar_m != 0) %>% #filter out 0m (scar) since this is the values we're comparing kr to
  group_by(site, distance_from_scar_m) %>% #give output for each site*distance grouping
  group_modify(~ {
    scar_vals <- red_k_boot_scar %>%
      filter(site == .y$site,
             distance_from_scar_m == 0) %>% 
      pull(k_boot) #pull scar kr values from every site - pulled from kr bootstrap values
    tibble(
      delta_kr = scar_vals - .x$k_boot #difference between scar kr bootstrap values and values at all other sites and distances
    )
  }) %>%
  ungroup()

#create summary table of kr values comparing distance from scar
delta_kr_scar_summary <- delta_kr_scar %>%
  group_by(site, distance_from_scar_m) %>%
  summarise(
    median_delta = median(delta_kr, na.rm = TRUE),
    lwr = quantile(delta_kr, 0.025, na.rm = TRUE),
    upr = quantile(delta_kr, 0.975, na.rm = TRUE),
    prob_scar_higher = mean(delta_kr > 0, na.rm = TRUE),
    .groups = "drop"
  )

delta_kr_scar_summary
#Median delta shows the effect size, prop_scar_higher is the probability of all bootstraped kr's for that pairing (site * distance) will be higher in scar 

#Create df to visualize differences between scar vs meadow for kr
red_k_boot_scar_summary <- red_k_boot_summary %>%
  filter(distance_from_edge_m == "30")

ggplot(red_k_boot_scar_summary,
       aes(x = factor(distance_from_scar_m,
                      levels = c(0, 0.25, 1, 10)), y = k_mean,
           ymin = k_lwr, ymax = k_upr, color = treatment)) +
  geom_pointrange(position = position_dodge(width = 1.5)) +
  scale_color_manual(
    values = c(scarred = "black", vegetated = "firebrick") 
  ) +
  facet_wrap(~ site) +
  theme_classic() +
  labs(
    x = "Distance from Scar (m)",
    y = "Red Tea Decay rate (Kr)",
    title = "Kr vs Distance from Scar by Treatment (Bootstrap w 95% CI)",
    color = "Treatment"
  )









