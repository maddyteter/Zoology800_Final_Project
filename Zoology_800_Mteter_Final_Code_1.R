#####Zoology 800 Final Project Code - True k values grouped by site PLUS 
#environmental variable comparison 


##### Part 1: Calculate a single, true k (no bootstrapping) for each tea type at 
#each unique site position (site*site_position = 18 per tea type). 
#Use these k values to compare against environmental data

library(tidyr)
library(broom)
library(lmerTest)
library(car)

a_green <- 0.7455474  
a_red <- 0.2930870

green_model_data_2 <- decomp %>%
  filter(tea_type == "G") %>%
  filter(sample_ID != "PSD25_S1.2_t7_B_G")


# Fit the exponential decay model by site_position only
green_k_fits_2 <- green_model_data_2 %>%
  group_by(site, site_position) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ tryCatch(
      nls(
        formula = tea_dw_proportion ~ a_green * exp(-k * total_deployment_time_days) +
          (1 - a_green),
        data = .x,
        start = list(k = 0.02),
        control = nls.control(maxiter = 200, warnOnly = TRUE)
      ),
      error = function(e) NULL
    )),
    tidied = map(fit, ~ {
      if (is.null(.x)) {
        tibble(term = "k", estimate = NA_real_)
      } else {
        tryCatch(tidy(.x), error = function(e) tibble(term = "k", estimate = NA_real_))
      }
    }),
    augmented = map2(fit, data, ~ {
      if (is.null(.x)) return(NULL)
      tryCatch(augment(.x, data = .y), error = function(e) NULL)
    }),
    rsq = map_dbl(augmented, ~ if (!is.null(.x)) {
      y <- .x$tea_dw_proportion
      y_hat <- .x$.fitted
      1 - sum((y - y_hat)^2) / sum((y - mean(y))^2)
    } else NA_real_)
  ) %>%
  unnest(tidied) %>%
  filter(term == "k") %>%
  rename(k = estimate)

# Create prediction grid for site_position-based curves
predict_time_green_2 <- expand_grid(
  green_model_data_2 %>% distinct(site, site_position),
  total_deployment_time_days = seq(
    min(green_model_data_2$total_deployment_time_days, na.rm = TRUE),
    max(green_model_data_2$total_deployment_time_days, na.rm = TRUE),
    length.out = 100
  )
)

# Add k values and generate predictions
predict_decay_green_2 <- predict_time_green_2 %>%
  left_join(
    green_k_fits_2 %>%
      select(site, site_position, k),
    by = c("site", "site_position")
  ) %>%
  mutate(
    pred_frac_green_2 = a_green * exp(-k * total_deployment_time_days) + (1 - a_green),
    pred_2 = pred_frac_green_2 * 3
  )

# Merge k values into the original dataset
green_model_data_2 <- green_model_data_2 %>%
  left_join(
    predict_decay_green_2 %>% distinct(site, site_position, k),
    by = c("site", "site_position"))


# Plot observed vs predicted for green tea decomposition by site_position
ggplot(green_model_data_2, aes(x = total_deployment_time_days)) +
  geom_point(aes(y = tea_dw_proportion), alpha = 0.4) +
  geom_line(
    data = predict_decay_green_2,
    aes(x = total_deployment_time_days, y = pred_frac_green_2),
    color = "seagreen4"
  ) +
  facet_wrap(~ site + site_position) +
  labs(
    x = "Time (days)",
    y = "Proportion Mass Remaining",
    title = "Observed vs Predicted Green Tea Mass Proportion (Grouped by Site Position)"
  ) +
  theme_minimal()


#### Plot mean k values for distance from scar by treatment
kg_summary_treat_2 <- green_model_data_2 %>%
  filter(distance_from_edge_m == "30") %>%
  group_by(distance_from_scar_m, treatment) %>%
  summarise(
    k_mean = mean(k, na.rm = TRUE),
    k_sd = sd(k, na.rm = TRUE),
    n = n()
  )
kg_summary_treat_2

#Create a time sequence for the curve to fit to 
time_seq_2 <- seq(
  min(green_model_data_2$total_deployment_time_days, na.rm = TRUE),
  max(green_model_data_2$total_deployment_time_days, na.rm = TRUE),
  length.out = 100
)

#Create predicted curve
predict_decay_treat_green_2 <- kg_summary_treat_2 %>%
  rowwise() %>%
  mutate(
    pred_df = list(
      tibble(
        total_deployment_time_days = time_seq_2,
        tea_dw_proportion = a_green * exp(-k_mean * time_seq_2) + (1 - a_green)
      )
    )
  ) %>%
  unnest(pred_df)

green_model_data_2 <- green_model_data_2 %>%
  mutate(distance_from_scar_m = factor(distance_from_scar_m,
                                       levels = c("0", "0.25", "1", "10")))
predict_decay_treat_green_2 <- predict_decay_treat_green_2 %>%
  mutate(distance_from_scar_m = factor(distance_from_scar_m,
                                       levels = c("0", "0.25", "1", "10")))


#Plot green tea data decomp curve with sites grouped together, separate curve by
#distance from scar
ggplot() +
  geom_point(data = green_model_data_2,
             aes(x = total_deployment_time_days, 
                 y = tea_dw_proportion, 
                 color = distance_from_scar_m),
             alpha = 0.4) +
  geom_line(data = predict_decay_treat_green_2,
            aes(x = total_deployment_time_days,
                y = tea_dw_proportion,
                color = distance_from_scar_m,
                linetype = treatment),
            linewidth = 0.8) +
  scale_color_manual(values = c("0" = "brown", "0.25" = "forestgreen", 
                                "1" = "skyblue", "10" = "orange1")) +
  scale_linetype_manual(values = c("scarred" = "dotted", "vegetated" = "solid")) +
  theme_classic() +
  labs(
    title = "Mean Green Tea Decomposition Curves by Distance From Scar and Treatment",
    x = "Time (days)",
    y = "Proportion Mass Remaining",
    color = "Distance from Scar (m)",
    linetype = "Treatment"
  )






#### Repeat red tea by site position only 
# Filter for green tea and remove “do not use” row
red_model_data_2 <- decomp %>%
  filter(tea_type == "R") %>%
  filter(!sample_ID %in% c("PSD25_S2.3_t7_C_R", "PSD25_S2.3_t2_A_R",
                           "PSD25_S2.3_t4_B_R"))


# Fit the exponential decay model by site_position only
red_k_fits_2 <- red_model_data_2 %>%
  group_by(site, site_position) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ tryCatch(
      nls(
        formula = tea_dw_proportion ~ a_red * exp(-k * total_deployment_time_days) +
          (1 - a_red),
        data = .x,
        start = list(k = 0.02),
        control = nls.control(maxiter = 200, warnOnly = TRUE)
      ),
      error = function(e) NULL
    )),
    tidied = map(fit, ~ {
      if (is.null(.x)) {
        tibble(term = "k", estimate = NA_real_)
      } else {
        tryCatch(tidy(.x), error = function(e) tibble(term = "k", estimate = NA_real_))
      }
    }),
    augmented = map2(fit, data, ~ {
      if (is.null(.x)) return(NULL)
      tryCatch(augment(.x, data = .y), error = function(e) NULL)
    }),
    rsq = map_dbl(augmented, ~ if (!is.null(.x)) {
      y <- .x$tea_dw_proportion
      y_hat <- .x$.fitted
      1 - sum((y - y_hat)^2) / sum((y - mean(y))^2)
    } else NA_real_)
  ) %>%
  unnest(tidied) %>%
  filter(term == "k") %>%
  rename(k = estimate)

# Create prediction grid for site_position-based curves
predict_time_red_2 <- expand_grid(
  red_model_data_2 %>% distinct(site, site_position),
  total_deployment_time_days = seq(
    min(red_model_data_2$total_deployment_time_days, na.rm = TRUE),
    max(red_model_data_2$total_deployment_time_days, na.rm = TRUE),
    length.out = 100
  )
)

# Add k values and generate predictions
predict_decay_red_2 <- predict_time_red_2 %>%
  left_join(
    red_k_fits_2 %>%
      select(site, site_position, k),
    by = c("site", "site_position")
  ) %>%
  mutate(
    pred_frac_red_2 = a_red * exp(-k * total_deployment_time_days) + (1 - a_red),
    pred_2 = pred_frac_red_2 * 3
  )

# Merge k values into the original dataset
red_model_data_2 <- red_model_data_2 %>%
  left_join(
    predict_decay_red_2 %>% distinct(site, site_position, k),
    by = c("site", "site_position")
  )

# Plot observed vs predicted for green tea decomposition by site_position
ggplot(red_model_data_2, aes(x = total_deployment_time_days)) +
  geom_point(aes(y = tea_dw_proportion), alpha = 0.4) +
  geom_line(
    data = predict_decay_red_2,
    aes(x = total_deployment_time_days, y = pred_frac_red_2),
    color = "firebrick3"
  ) +
  facet_wrap(~ site + site_position) +
  labs(
    x = "Time (days)",
    y = "Proportion Mass Remaining",
    title = "Observed vs Predicted Red Tea Mass Proportion (Grouped by Site Position)"
  ) +
  theme_minimal()



#Red tea model with distance from scar * treatment type

kr_summary_treat_2 <- red_model_data_2 %>%
  filter(distance_from_edge_m == "30") %>%
  group_by(distance_from_scar_m, treatment) %>%
  summarise(
    k_mean = mean(k, na.rm = TRUE),
    k_sd = sd(k, na.rm = TRUE),
    n = n()
  )
kr_summary_treat_2

#Create a time sequence for the curve to fit to 
time_seq_2 <- seq(
  min(red_model_data_2$total_deployment_time_days, na.rm = TRUE),
  max(red_model_data_2$total_deployment_time_days, na.rm = TRUE),
  length.out = 100
)

#Create predicted curve
predict_decay_treat_red_2 <- kr_summary_treat_2 %>%
  rowwise() %>%
  mutate(
    pred_df = list(
      tibble(
        total_deployment_time_days = time_seq_2,
        tea_dw_proportion = a_red * exp(-k_mean * time_seq_2) + (1 - a_red)
      )
    )
  ) %>%
  unnest(pred_df)

red_model_data_2 <- red_model_data_2 %>%
  mutate(distance_from_scar_m = factor(distance_from_scar_m,
                                       levels = c("0", "0.25", "1", "10")))
predict_decay_treat_red_2 <- predict_decay_treat_red_2 %>%
  mutate(distance_from_scar_m = factor(distance_from_scar_m,
                                       levels = c("0", "0.25", "1", "10")))


#Plot red tea data decomp curve with sites grouped together, separate curve by
#distance from scar
ggplot() +
  geom_point(data = red_model_data_2,
             aes(x = total_deployment_time_days, 
                 y = tea_dw_proportion, 
                 color = distance_from_scar_m),
             alpha = 0.4) +
  geom_line(data = predict_decay_treat_red_2,
            aes(x = total_deployment_time_days,
                y = tea_dw_proportion,
                color = distance_from_scar_m,
                linetype = treatment),
            linewidth = 0.8) +
  scale_color_manual(values = c("0" = "brown", "0.25" = "forestgreen", 
                                "1" = "skyblue", "10" = "orange1")) +
  scale_linetype_manual(values = c("scarred" = "dotted", "vegetated" = "solid")) +
  theme_classic() +
  labs(
    title = "Mean Red Tea Decomposition Curves by Distance From Scar and Treatment",
    x = "Time (days)",
    y = "Proportion Mass Remaining",
    color = "Distance from Scar (m)",
    linetype = "Treatment"
  )






################# Part Two - Load Habitat Characteristics Data 
read.csv("PSD25_aboveground_biomass.csv")

#create agb df
agb <- read.csv("PSD25_aboveground_biomass.csv")
View(agb)

#create column to get final DW 
agb <- agb %>%
  mutate(final_agb_dw = (dry_mass_g - dish_mass_g)*100)


#create standard error function
se <- function(x) sd(x, na.rm = TRUE) / sqrt(length((x)))

#calculate average agb dw 
agb_average <- agb%>%
  group_by(site_id, site, site_position, species, time_interval_wks) %>%
  summarise(
    mean_agb_species_dw = mean(final_agb_dw, na.rm = TRUE),
    se_agb_species_dw = se(final_agb_dw),
    .groups = "drop"
  )

#merge columns distance from edge, distance from scar, and treatment from site
# data df to agb_average

agb <- agb %>%
  left_join(site_data %>% select(site_id, treatment, 
                                 distance_from_edge_m, distance_from_scar_m), by = "site_id")

#merge columns mean_agb_species_dw and se_agb_species_dw into agb df
agb <- agb %>%
  group_by(species, time_interval_wks) %>%
  left_join(
    agb_average %>% 
      select(site_id, species, time_interval_wks, mean_agb_species_dw, se_agb_species_dw),
    by = c("site_id", "species", "time_interval_wks")  
  )

agb <- agb %>%
  group_by(site_id, site, site_position, time_interval_wks) %>%
  mutate(mean_agb_total_dw = mean(final_agb_dw),
         se_agb_total_dw = se(final_agb_dw),
         .groups = "drop"
  )

#Filter for only thalassia 
agb_tt_t7 <- agb %>%
  filter(time_interval_wks == "t7", species == "Tt") %>%
  select(site_id, site, site_position, mean_agb_species_dw, se_agb_species_dw)
#Collapse down into one row per site id

agb_tt_t7 <- agb_tt_t7 %>%
  group_by(site_id) %>%
  mutate(mean_agb_species_dw = mean(mean_agb_species_dw, na.rm = TRUE)) %>%
  slice(1) %>%              # keep only one representative row per site_id
  ungroup()

#######Belowground Biomass

#create bgb df
bgb <- read.csv("PSD25_belowground_biomass.csv")

#add column to get final bgb DW (without dish)
bgb <- bgb%>%
  mutate(final_bgb_dw = dry_weight_g - dish_weight_g)

#Merge site data to bgb df
bgb <- bgb %>%
  left_join(site_data %>% select(site_id, treatment, 
                                 distance_from_edge_m, 
                                 distance_from_scar_m), by = "site_id")

######Sediment Bulk Density
sedbd <- read.csv("PSD25_sediment_bulk_density.csv")

View(sedbd)

#create column for final dry mass (g)
sedbd <- sedbd %>%
  mutate(final_bd_gml = (dry_mass_g - wp_mass_g)/sample_vol_ml)

#Calculate mean and se for sedbd (mg/L)
sedbd_ave <- sedbd %>%
  group_by(site_id) %>%
  summarise(
    mean_sedbd_gml = mean(final_bd_gml, na.rm = TRUE),
    se_sedbd_gml = se(final_bd_gml))


#Merge treatment, dist from edge, and dist from scar in site_data to sedbd

sedbd <- sedbd %>%
  left_join(site_data %>% select(site_id, treatment, distance_from_edge_m,
                                 distance_from_scar_m), by = "site_id")

#Merge mean_sedbd_mgl and se_sedbd_mgl to sedbd df

sedbd <- sedbd %>%
  left_join(sedbd_ave %>% select(site_id, mean_sedbd_gml, 
                                 se_sedbd_gml), by = "site_id")

######Sediment Grain Size 
sedgs <- read.csv("PSD25_sediment_grain_size.csv")
View(sedgs)

#Create column for final dw (g)
sedgs <- sedgs %>%
  mutate(final_sedgs_dw_g = dry_mass_g - wp_mass_g)

#Merge site data into sediment grain size df

sedgs <- sedgs %>%
  left_join(site_data %>% select(site_id, treatment, 
                                 distance_from_scar_m, distance_from_edge_m), 
            by = "site_id")


#Join mean_sedbd column from sedbd df to sedgs df

sedgs <- sedgs %>%
  left_join(sedbd_ave %>% select(site_id, mean_sedbd_gml, 
                                 se_sedbd_gml), by = "site_id")

#Alter dataset from long to wide to allow calculation of fine particle gs dw (g)


#Pivot to wide format — one row per sample_id/site_id with gravel type dw listed
#as multiple columns

#Pivot to wide format — one row per sample_id
sedgs_wide <- sedgs %>%
  pivot_wider(
    id_cols = c(sample_id, mean_sedbd_gml),  # keep sample-specific columns
    names_from = grain_type,
    values_from = final_sedgs_dw_g,
    names_prefix = "",
    values_fn = list(final_sedgs_dw_g = mean)  # in case of duplicates
  )

#Calculate fine where it's NA 
#((mean_sedbd_gml*60(sedgs sample size)) - gravel_dw_g - sand_dw_g = fine_dw_g
#for each sample id/ grain size sample collected
sedgs_wide <- sedgs_wide %>%
  mutate(
    fine = ifelse(
      is.na(fine),
      (mean_sedbd_gml * 60) - gravel - sand,
      fine
    )
  )

# Step 3: Join back to original metadata (optional)
# This brings back all other columns that were dropped in pivot_wider
sedgs_updated <- sedgs %>%
  pivot_wider(names_from = grain_type, values_from = final_sedgs_dw_g) %>%
  left_join(sedgs_wide, by = c("sample_id", "mean_sedbd_gml")) %>%
  mutate(
    gravel = coalesce(gravel.y, gravel.x),
    sand = coalesce(sand.y, sand.x),
    fine = coalesce(fine.y, fine.x)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"))

# Step 4: Pivot back to long format (optional)
sedgs_final <- sedgs_updated %>%
  pivot_longer(
    cols = c(gravel, sand, fine),
    names_to = "grain_type",
    values_to = "final_sedgs_dw_g"
  )

#Create column for proportion_sedgs_dw for each grain type in each sample

sedgs_final <- sedgs_final %>%
  mutate(proportion_sedgs_dw_g = final_sedgs_dw_g/(mean_sedbd_gml*60))

###### Sediment Organic Matter
sed_om <- read.csv("PSD25_sediment_om.csv")
View(sed_om)

#Calculate mass (g) of total sample before ashing
sed_om <- sed_om %>%
  mutate(total_sample_mass_g = dry_mass_crucible_sample - crucible_mass)

#Calculate mass(g) of organic sample after ashing (what was remaining in crucible)
sed_om <- sed_om %>%
  mutate(organic_sample_mass_g = dry_mass_crucible_sample - ashed_mass_crucible_sample) 

#Calculate mass (g) of inorganic sample after ashing 
sed_om <- sed_om %>%
  mutate(inorganic_sample_mass_g = total_sample_mass_g - organic_sample_mass_g)

#Calculate average sed om (g) for each site_id
sed_om_average <- sed_om %>%
  group_by(site_id, site, site_position) %>%
  mutate(mean_sed_om = mean(organic_sample_mass_g, na.rm = TRUE),
         se_sed_om = se(organic_sample_mass_g))

#Calculate proportion sed om 
sed_om_average <- sed_om_average %>%
  group_by(site_id, site, site_position) %>%
  mutate(sed_om_proportion = organic_sample_mass_g/total_sample_mass_g)

#Calculate average proporiton sed om
sed_om_average <- sed_om_average %>%
  group_by(site_id, site, site_position) %>%
  mutate(mean_sed_om_proportion = mean(sed_om_proportion, na.rm = TRUE),
         se_sed_om = se(sed_om_proportion))


#Join site data 
sed_om_average <- sed_om_average %>%
  left_join(site_data %>% select(site_id, treatment, 
                                 distance_from_edge_m, distance_from_scar_m),
            by = "site_id")

######Shoot Density 
shoot_density <- read.csv("PSD25_shoot_density.csv")
View(shoot_density)

#Format wide dataframe to long dataframe so species and count become separate columns

#Reformat csv from wide to long 
long_shoot_density <- shoot_density %>%
  pivot_longer(
    cols = c("Tt_count", "Tt_seed", "Sf_count", "Hw_count", "Hal_count", 
             "Udo_count", "Pen_count"), 
    names_to = "species",
    values_to = "shoot_count"
  )

#Rename species column values 
long_shoot_density <- long_shoot_density %>%
  mutate(species = recode(species, "Tt_count" = "Tt",
                          "Sf_count" = "Sf",
                          "Hw_count" = "Hw",
                          "Tt_seed_count" = "Tt Seed",
                          "Hal_count" = "Hal",
                          "Udo_count" = "Udo",
                          "Pen_count" = "Pen"))

#Remove drift algae and sponge count columns
long_shoot_density <- long_shoot_density %>%
  select(-Drift_count, -Sponge_count)


#Calculate average shoot count by species
shootdensity_average <- long_shoot_density %>%
  group_by(site_id, site, site_position, species) %>%
  mutate(
    mean_shoot_count = mean(shoot_count, na.rm = TRUE),
    se_shoot_count = se(shoot_count))


#Calculate average total shoot density (not grouped by species)
shootdensity_average <- shootdensity_average %>%
  group_by(site_id, site, site_position) %>%
  mutate(total_mean_shoot_density = mean(shoot_count/quadrat_area_m2),
         se_shoot_density = se(shoot_count/quadrat_area_m2))

#Calculate average shoot density by species type

shootdensity_average <- shootdensity_average %>%
  group_by(site_id, site, site_position, species) %>%
  mutate(species_mean_shoot_density = mean(shoot_count/quadrat_area_m2),
         species_se_shoot_count = se(shoot_count/quadrat_area_m2))


#Join site data
shootdensity_average <- shootdensity_average %>%
  left_join(site_data %>% select(site_id, treatment, 
                                 distance_from_edge_m, distance_from_scar_m),
            by = "site_id")


###### Thalassia Blade Morphometry
tt_morphometry <- read.csv("PSD25_thalassia_morphometry.csv")

#Calculate average blade width (cm) 

tt_morphometry <- tt_morphometry %>%
  group_by(site_id, site, site_position, time_interval_wks) %>%
  mutate(
    mean_blade_width_cm = mean(blade_width_cm, na.rm = TRUE),
    se_blade_width_cm = se(blade_width_cm))

#Calculate average blade length (cm)

tt_morphometry <- tt_morphometry %>%
  group_by(site_id, site, site_position, time_interval_wks) %>%
  mutate(
    mean_blade_length_cm = mean(blade_length_cm, na.rm = TRUE),
    se_blade_length_cm = se(blade_length_cm))

#Join site data
tt_morphometry <- tt_morphometry %>%
  left_join(site_data %>% select(site_id, treatment, 
                                 distance_from_edge_m, distance_from_scar_m), by = "site_id")

tt_morphometry_t7 <- tt_morphometry %>%
  filter(time_interval_wks == "t7")








############# Part Three - Load Environmental Variable Characteristics 
#####AGB
# collapse AGB data so there’s one row per site_id
kg_tt_agb <- agb_tt_t7 %>%
  group_by(site_id) %>%
  summarise(mean_agb_species_dw = mean(mean_agb_species_dw, na.rm = TRUE))

#Join to green_explanatory_data
green_explanatory_data <- green_model_data %>%
  left_join(kg_tt_agb, by = "site_id")

#####BGB

kg_bgb <- bgb %>%
  group_by(site_id) %>%
  summarise(final_bgb_dw = mean(final_bgb_dw, na.rm = TRUE))

green_explanatory_data <- green_explanatory_data %>%
  left_join(kg_bgb %>% select(site_id, final_bgb_dw),
            by = "site_id")

####SEDBD
kg_bd <- sedbd %>%
  group_by(site_id) %>%
  summarise(mean_sedbd_gml = mean(mean_sedbd_gml, na.rm = TRUE))

#Join to green_explanatory_data
green_explanatory_data <- green_explanatory_data %>%
  left_join(kg_bd, by = "site_id")


###SED OM
kg_sed_om <- sed_om_average %>%
  group_by(site_id) %>%
  summarise(mean_sed_om_proportion = mean(mean_sed_om_proportion, na.rm = TRUE))

green_explanatory_data <- green_explanatory_data %>%
  left_join(kg_sed_om, by = "site_id")

###SEDGS FINE

kg_sedgs_fine <- sedgs_final %>%
  filter(grain_type == "fine") %>%
  group_by(site_id) %>%
  rename(proportion_sedgs_dw_g_fine = proportion_sedgs_dw_g) %>%
  summarise(proportion_sedgs_dw_g_fine = mean(proportion_sedgs_dw_g_fine, na.rm = TRUE))

#Join to green_explanatory_data
green_explanatory_data <- green_explanatory_data %>%
  left_join(kg_sedgs_fine, by = "site_id")

###SEDGS SAND

kg_sedgs_sand <- sedgs_final %>%
  filter(grain_type == "sand") %>%
  group_by(site_id) %>%
  rename(proportion_sedgs_dw_g_sand = proportion_sedgs_dw_g) %>%
  summarise(proportion_sedgs_dw_g_sand = mean(proportion_sedgs_dw_g_sand, na.rm = TRUE))

#Join to green_explanatory_data
green_explanatory_data <- green_explanatory_data %>%
  left_join(kg_sedgs_sand, by = "site_id")


###SEDGS GRAVEL

kg_sedgs_gravel <- sedgs_final %>%
  filter(grain_type == "gravel") %>%
  group_by(site_id) %>%
  rename(proportion_sedgs_dw_g_gravel = proportion_sedgs_dw_g) %>%
  summarise(proportion_sedgs_dw_g_gravel = mean(proportion_sedgs_dw_g_gravel, na.rm = TRUE))

#Join to green_explanatory_data
green_explanatory_data <- green_explanatory_data %>%
  left_join(kg_sedgs_gravel, by = "site_id")

###SHOOT DENSITY
kg_shoot_dens <- shootdensity_average %>%
  group_by(site_id) %>%
  summarise(total_mean_shoot_density = mean(total_mean_shoot_density, na.rm = TRUE))

green_explanatory_data <- green_explanatory_data %>%
  left_join(kg_shoot_dens, by = "site_id")

###TT BLADE LENGTH
kg_blade_length <- tt_morphometry_t7 %>%
  group_by(site_id) %>%
  summarise(mean_blade_length_cm = mean(mean_blade_length_cm, na.rm = TRUE))

green_explanatory_data <- green_explanatory_data %>%
  left_join(kg_blade_length, by = "site_id")


#RED TEA - Merge explanatory variables into one dataframe 
#####AGB
# collapse AGB data so there’s one row per site_id
kr_tt_agb <- agb_tt_t7 %>%
  group_by(site_id) %>%
  summarise(mean_agb_species_dw = mean(mean_agb_species_dw, na.rm = TRUE))

#Join to green_explanatory_data
red_explanatory_data <- red_model_data %>%
  left_join(kr_tt_agb, by = "site_id")

#####BGB

kr_bgb <- bgb %>%
  group_by(site_id) %>%
  summarise(final_bgb_dw = mean(final_bgb_dw, na.rm = TRUE))

red_explanatory_data <- red_explanatory_data %>%
  left_join(kr_bgb %>% select(site_id, final_bgb_dw),
            by = "site_id")

####SEDBD

kr_bd <- sedbd %>%
  group_by(site_id) %>%
  summarise(mean_sedbd_gml = mean(mean_sedbd_gml, na.rm = TRUE))

#Join to green_explanatory_data
red_explanatory_data <- red_explanatory_data %>%
  left_join(kr_bd, by = "site_id")


###SED OM

kr_sed_om <- sed_om_average %>%
  group_by(site_id) %>%
  summarise(mean_sed_om_proportion = mean(mean_sed_om_proportion, na.rm = TRUE))

red_explanatory_data <- red_explanatory_data %>%
  left_join(kr_sed_om, by = "site_id")


###SEDGS FINE

kr_sedgs_fine <- sedgs_final %>%
  filter(grain_type == "fine") %>%
  group_by(site_id) %>%
  rename(proportion_sedgs_dw_g_fine = proportion_sedgs_dw_g) %>%
  summarise(proportion_sedgs_dw_g_fine = mean(proportion_sedgs_dw_g_fine, na.rm = TRUE))

#Join to green_explanatory_data
red_explanatory_data <- red_explanatory_data %>%
  left_join(kr_sedgs_fine, by = "site_id")

###SEDGS SAND

kr_sedgs_sand <- sedgs_final %>%
  filter(grain_type == "sand") %>%
  group_by(site_id) %>%
  rename(proportion_sedgs_dw_g_sand = proportion_sedgs_dw_g) %>%
  summarise(proportion_sedgs_dw_g_sand = mean(proportion_sedgs_dw_g_sand, na.rm = TRUE))

#Join to green_explanatory_data
red_explanatory_data <- red_explanatory_data %>%
  left_join(kr_sedgs_sand, by = "site_id")


###SEDGSGRAVEL

kr_sedgs_gravel <- sedgs_final %>%
  filter(grain_type == "gravel") %>%
  group_by(site_id) %>%
  rename(proportion_sedgs_dw_g_gravel = proportion_sedgs_dw_g) %>%
  summarise(proportion_sedgs_dw_g_gravel = mean(proportion_sedgs_dw_g_gravel, na.rm = TRUE))

#Join to green_explanatory_data
red_explanatory_data <- red_explanatory_data %>%
  left_join(kr_sedgs_gravel, by = "site_id")

###SHOOT DENSITY
kr_shoot_dens <- shootdensity_average %>%
  group_by(site_id) %>%
  summarise(total_mean_shoot_density = mean(total_mean_shoot_density, na.rm = TRUE))

red_explanatory_data <- red_explanatory_data %>%
  left_join(kr_shoot_dens, by = "site_id")

###TT BLADE LENGTH
kr_blade_length <- tt_morphometry_t7 %>%
  group_by(site_id) %>%
  summarise(mean_blade_length_cm = mean(mean_blade_length_cm, na.rm = TRUE))

red_explanatory_data <- red_explanatory_data %>%
  left_join(kr_blade_length, by = "site_id")

#collapse green and red explanatory data so it's only 1 row per k value 


#combine red and green explanatory data

combined_explanatory_data <- bind_rows(green_explanatory_data, red_explanatory_data)




################ Part 5 -Compare Kg and Kr against environmental variables 
#to see which ones have the most influence


####Green Tea Data 
green_explanatory_data <- green_explanatory_data %>%
  left_join(green_k_fits_2 %>% select(site, site_position, k),
            by = c("site", "site_position"))

green_explanatory_data <- green_explanatory_data %>%
  group_by(site_id) %>% 
  summarize(
    across(
      everything(),
      ~ if (is.numeric(.x)) mean(.x, na.rm = TRUE) else first(.x)),
    .groups = "drop")


#Site as random effect 
kg_enviro_site_random_results <- lmer(
  k ~ proportion_sedgs_dw_g_fine + 
    mean_sed_om_proportion + (1 | site), data = green_explanatory_data
)
summary(kg_enviro_site_random_results)

vif(kg_enviro_site_random_results)

par(mfrow = c(2, 2))
hist(residuals(kg_enviro_site_random_results), main = "A")
kg_enviro_site_random_results_residuals <- residuals(kg_enviro_site_random_results, type = 'pearson')
kg_enviro_site_random_results_fitted <- fitted.values(kg_enviro_site_random_results)
qqnorm(kg_enviro_site_random_results_residuals, main = "B")
qqline(kg_enviro_site_random_results_residuals)
plot(
  kg_enviro_site_random_results_fitted,
  kg_enviro_site_random_results_residuals,
  main = "C: Fitted vs. Residuals",
  xlab = "Fitted values",
  ylab = "Pearson residuals"
)
abline(h = 0, col = "red", lwd = 2)
shapiro.test(residuals(kg_enviro_site_random_results))

#Repeat with site as fixed effect
kg_enviro_site_fixed_results <- lm(
  k ~ proportion_sedgs_dw_g_fine + 
    mean_sed_om_proportion + site, data = green_explanatory_data
)
summary(kg_enviro_site_fixed_results)

vif(kg_enviro_site_fixed_results)

#Test assumptions
par(mfrow = c(2, 2))
hist(residuals(kg_enviro_site_fixed_results), main = "A")
kg_enviro_site_fixed_results_residuals <- residuals(kg_enviro_site_fixed_results, type = 'pearson')
kg_enviro_site_fixed_results_fitted <- fitted.values(kg_enviro_site_fixed_results)
qqnorm(kg_enviro_site_fixed_results_residuals, main = "B")
qqline(kg_enviro_site_fixed_results_residuals)
plot(
  kg_enviro_site_fixed_results_fitted,
  kg_enviro_site_fixed_results_residuals,
  main = "C: Fitted vs. Residuals",
  xlab = "Fitted values",
  ylab = "Pearson residuals"
)
abline(h = 0, col = "red", lwd = 2)
shapiro.test(residuals(kg_enviro_site_fixed_results))
#normal residuals (small sample size but qq plot and residuals look pretty good
#for not having many kg values - looks normal enought)


AIC(kg_enviro_site_fixed_results) #site as fixed effect is better model 
AIC(kg_enviro_site_random_results)



#####Test Red Tea Data 
red_explanatory_data <- red_explanatory_data %>%
  left_join(red_k_fits_2 %>% select(site, site_position, k),
            by = c("site", "site_position"))

red_explanatory_data <- red_explanatory_data %>%
  group_by(site_id) %>% 
  summarize(
    across(
      everything(),
      ~ if (is.numeric(.x)) mean(.x, na.rm = TRUE) else first(.x)),
    .groups = "drop")

kr_enviro_site_random_results <- lmer(
  k ~ final_bgb_dw +  
     (1 | site), data = red_explanatory_data
)
summary(kr_enviro_site_random_results)
vif(kr_enviro_site_random_results)

#Repeat with site as fixed variable 
kr_enviro_site_fixed_results <- lm(
  k ~ final_bgb_dw + 
    site, data = red_explanatory_data
)
summary(kr_enviro_site_fixed_results)
vif(kr_enviro_site_fixed_results)

AIC(kr_enviro_site_fixed_results) #better AIC score
AIC(kr_enviro_site_random_results)

#Check assumptions
par(mfrow = c(2, 2))
hist(residuals(kr_enviro_site_fixed_results), main = "A")
kr_enviro_site_fixed_results_residuals <- residuals(kr_enviro_site_fixed_results, type = 'pearson')
kr_enviro_site_fixed_results_fitted <- fitted.values(kr_enviro_site_fixed_results)
qqnorm(kr_enviro_site_fixed_results_residuals, main = "B")
qqline(kr_enviro_site_fixed_results_residuals)
plot(
  kr_enviro_site_fixed_results_fitted,
  kr_enviro_site_fixed_results_residuals,
  main = "C: Fitted vs. Residuals",
  xlab = "Fitted values",
  ylab = "Pearson residuals"
)
abline(h = 0, col = "red", lwd = 2)
shapiro.test(residuals(kr_enviro_site_fixed_results))


#Test variables to see how they differ by treatment or distance from scar
fine_distance <- lm(proportion_sedgs_dw_g_fine ~ mean_agb_species_dw,
                    data = green_explanatory_data)
summary(fine_distance)
#fine sediments has to be pieced a part more spatially - but we do know that agb
#is very correlated with fine sed grain size through particle trapping -->
#more agb in vegetation = more fine sed in vegetation

bgb_distance <- lm(final_bgb_dw ~ distance_from_scar_m + site,
                   data = red_explanatory_data)
summary(bgb_distance)

#bgb associated with distance from scar - could explain decrease in kr between 
#scar and vegetated positions
