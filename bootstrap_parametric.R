library(surveybootstrap)
library(mgcv)
library(tidyr)
library('dplyr')
library(readr)
library(brms)
library(ggplot2)
library(scales)
library(tibble)
library(survey)

#set working directory
setwd()

###############################################################################
######################## PREP DATA AND IMPUTE RACE ############################
###############################################################################

#import neiss and population data

neiss_legal <- read_csv("neiss_legal_race_imputed.csv")
neiss_legal <- neiss_legal %>% 
  mutate(
    hid = as.factor(hid),
    male = as.factor(male),
    STRATUM = as.factor(STRATUM)
  )

population <- read_csv("postcensus_pop.csv")
population <- population %>%
  rename(
    total_pop = total,
    white = nhw, 
    black = b, 
    latinx = hisp
  ) %>%
  mutate(other = total_pop - (white + black + latinx)) %>%
  pivot_longer(
    cols = c(latinx, white, black, other), 
    names_to = "race",  
    values_to = "pop"  
  ) 

######### Predict race with multilevel, multinomial logistic model ############
# (No need to specify complete case, as mgcv will drop rows w/ any missingness)

# Note: For a multinomial model, GAM wants a list of predictors
# for all but the reference level (raceeth = 0; NH white) of the outcome.
# Only the first formula should specify the outcome variable name.
# The formula will otherwise be identical for all 3 Ks

# Modify this to change the imputation formula
# Manually specifying number of spline knots to potentially increase speed
race_imp_formula <- "s(age, k=5) + male + STRATUM + s(injcount_monthly, k=5) + s(pct_assault, k=5) + s(pct_legal, k=5) + s(mean_age, k=5) + s(pct_male, k=5) + s(hid, bs = 're')"

# Combine formulas into a list or mcgv
race_formulae <- list(
  as.formula(paste("raceeth ~", race_imp_formula)),
  as.formula(paste("~", race_imp_formula)),
  as.formula(paste("~", race_imp_formula))
)

race_imp_model <- gam(race_formulae, family = multinom(K=3), 
                      data = neiss_legal, method = "GCV.Cp")

predicted_race <- predict(race_imp_model,
                          newdata = neiss_legal,
                          type = "response")

# Convert the predicted_probs matrix to a data frame
predicted_race_df <- as.data.frame(predicted_race)

# Name the columns based on your categories (e.g., Category 1, Category 2, etc.)
colnames(predicted_race_df) <- c("pred_white", "pred_black", "pred_latinx", "pred_other")

# Combine with the original dataframe
neiss_legal <- cbind(neiss_legal, predicted_race_df)

#replace predicted probabilities for observations with known race 
neiss_legal <- neiss_legal %>%
  mutate(
    pred_white = case_when(
      !is.na(raceeth) & raceeth == 0 ~ 1,
      !is.na(raceeth) & raceeth != 0 ~ 0,
      TRUE ~ pred_white
    ),
    pred_black = case_when(
      !is.na(raceeth) & raceeth == 1 ~ 1,
      !is.na(raceeth) & raceeth != 1 ~ 0,
      TRUE ~ pred_black
    ),
    pred_latinx = case_when(
      !is.na(raceeth) & raceeth == 2 ~ 1,
      !is.na(raceeth) & raceeth != 2 ~ 0,
      TRUE ~ pred_latinx
    ),
    pred_other = case_when(
      !is.na(raceeth) & raceeth == 3 ~ 1,
      !is.na(raceeth) & raceeth != 3 ~ 0,
      TRUE ~ pred_other
    )
  )

#save legal injury data with imputed race 
write.csv(neiss_legal, "neiss_legal_race_imputed.csv", row.names = FALSE)

###############################################################################
######################## PARAMETRIC BOOTSTRAPPING  ############################
###############################################################################

#import population data and race-imputed neiss data 

population <- read_csv("postcensus_pop.csv")
population <- population %>%
  rename(
    total_pop = total,
    white = nhw, 
    black = b, 
    latinx = hisp
  ) %>%
  mutate(other = total_pop - (white + black + latinx)) %>%
  pivot_longer(
    cols = c(latinx, white, black, other), 
    names_to = "race",  
    values_to = "pop"  
  ) %>%
  select(race, year, pop, total_pop) 

population <- population %>%
  group_by(year) %>%
  bind_rows(
    population %>%
      group_by(year) %>%
      summarise(race = "all", pop = unique(total_pop), .groups = "drop")  # Create "all" category
  ) %>%
  arrange(year, race) %>%
  select(year, race, pop)

neiss_legal<-read.csv("neiss_legal_race_imputed.csv")

neiss_legal <- neiss_legal %>%
  filter(
    !is.na(pred_white) & 
      !is.na(pred_black) & 
      !is.na(pred_latinx) & 
      !is.na(pred_other)
  )

#define survey design object 
neiss_design <- svydesign(
  id = ~PSU,            # Primary Sampling Unit
  strata = ~STRATUM,     # Strata
  weights = ~WEIGHT,    # Survey weights
  data = neiss_legal,
  nest = TRUE           # Use TRUE if PSU is nested within strata
)

set.seed(123)

#### define function to calculate yearly injury counts by race/ethnicity #####

calculate_count_se <- function(survey_design) {
  totals_by_year <- svyby(
    ~pred_white + pred_black + pred_latinx + pred_other,  # Variables to sum
    ~year,          # Group by year
    survey_design,  # Survey design object
    svytotal        # Function to apply
  )
  
  # Convert output into a tidy format
  results <- totals_by_year %>%
    pivot_longer(
      cols = starts_with("pred_"),   # Pivot race columns
      names_to = "race",  
      values_to = "total_count"
    ) %>%
    mutate(
      race = case_when(
        race == "pred_white"  ~ "white",
        race == "pred_black"  ~ "black",
        race == "pred_latinx" ~ "latinx", 
        race == "pred_other"  ~ "other"
      ),
      SE = pivot_longer(
        totals_by_year %>%
          select(starts_with("se.pred_")), 
        cols = everything(),
        values_to = "SE"
      )$SE
    )
  
  results<-results %>%
    select(year, race, total_count, SE) %>%
    rename(SE_survey=SE) %>%
    mutate(total_count = round(total_count))
  return(results)
}

neiss_byraceyear<-calculate_count_se(neiss_design)

####parametric bootstrapping to account for uncertainty from race imputation####
#set bootstrap sample number 
n_bootstrap = 10000

simulate_race_counts <- function(data) {
  set.seed(123)
  # Assign each row a race based on probabilities
  data<- data %>%
    rowwise() %>%
    mutate(
      race = sample(
        c("white", "black", "latinx", "other"), 
        size = 1, 
        prob = c(pred_white, pred_black, pred_latinx, pred_other), 
        replace = TRUE
      )
    ) %>%
    ungroup()
  
  # Step 2: Calculate weighted totals by year and race
  total_counts <- data %>%
    group_by(year, race) %>%
    summarize(total_count = round(sum(WEIGHT, na.rm = TRUE)), .groups = "drop") 
  
  return(total_counts)
}

race_simulation <- data.frame()

# Loop through bootstrap iterations
for (i in 1:n_bootstrap) {
  # Apply the simulate_race_counts function
  sim_result <- simulate_race_counts(neiss_legal) 
  
  # Add a column indicating the simulation number (for reference)
  sim_result$sim_id <- i
  
  race_simulation <- bind_rows(race_simulation, sim_result)
}

se_impute <- race_simulation %>%
  group_by(year, race) %>%
  summarize(
    SE_impute = sd(total_count, na.rm = TRUE),  
    .groups = "drop"  
  )

neiss_counts_se<-neiss_byraceyear %>% 
  left_join(se_impute, by=c("race", "year")) %>% 
  mutate(SE_total = sqrt(SE_survey^2 + SE_impute^2))

#get SE for total population 

  totals_by_year <- svyby(
    ~legal,  # Variables to sum
    ~year,          # Group by year
    neiss_design,  # Survey design object
    svytotal        # Function to apply
  )
  
  # Convert output into a tidy format
  results_all <- totals_by_year %>%
    pivot_longer(
      cols = starts_with("legal"),   # Pivot race columns
      names_to = "race",  
      values_to = "total_count"
    ) %>%
    mutate(
      race = "all",
    ) %>%
    select(year, race, total_count, se) %>%
    rename(SE_total=se) %>%
    mutate(total_count = round(total_count))
  
#merge population data to calculate injury rate by race and year 
  merged_raceyear <- neiss_counts_se %>% 
    select(year, race, total_count, SE_total) %>% 
    rbind(results_all) %>% 
    merge(population, by = c("year", "race")) %>% 
    mutate(injury_rate_100k = total_count/pop*100000, 
           SE_rate_100k = SE_total/pop*100000)

###############second parametric bootstrapping for modelling trend##############
#define parameters for gamma distribution based on the counts and SE 
merged_gamma <- merged_raceyear %>% 
  mutate(alpha = (total_count^2) / (SE_total^2),
         beta = total_count / (SE_total^2))

#generate parametric bootstrap counts 
simulate_gamma_counts <- function(data, num_simulations = n_bootstrap) {
  set.seed(123)
  simulations <- data.frame()
  # Loop through each row (each combination of race and year)
  for (i in 1:nrow(data)) {
    # Extract the alpha and beta for the current row
    alpha <- data$alpha[i]
    beta <- data$beta[i]
    
    # Generate random values from the Gamma distribution
    simulated_values <- rgamma(num_simulations, shape = alpha, rate = beta)
    
    # Create a data frame with the simulation results for the current observation
    sim_df <- data.frame(
      year = data$year[i],
      race = data$race[i],
      simulation_id = 1:num_simulations,
      simulated_value = simulated_values
    )
    
    # Append to the simulations data frame
    simulations <- bind_rows(simulations, sim_df)
  }
  
  return(simulations)
}

boot_gamma<-simulate_gamma_counts(merged_gamma) 
boot_gamma<- boot_gamma %>% 
  rename(total_count=simulated_value)

simulate_lnormal_counts <- function(data, num_simulations = n_bootstrap) {
  set.seed(123)
  simulations <- data.frame()
  # Loop through each row (each combination of race and year)
  for (i in 1:nrow(data)) {
    # Extract the alpha and beta for the current row
    mu <- data$mu[i]
    sigma <- data$sigma[i]
    
    # Generate random values from the Gamma distribution
    simulated_values <- rnorm(num_simulations, mean = mu, sd = sigma)
    
    # Create a data frame with the simulation results for the current observation
    sim_df <- data.frame(
      year = data$year[i],
      race = data$race[i],
      simulation_id = 1:num_simulations,
      simulated_value = simulated_values
    )
    
    # Append to the simulations data frame
    simulations <- bind_rows(simulations, sim_df)
  }
  
  return(simulations)
}

boot_lnormal<-simulate_lnormal_counts(merged_lnormal) 
boot_lnormal<- boot_lnormal %>% 
  mutate(simulated_value = exp(simulated_value)) %>% 
  rename(total_count=simulated_value)

###############################################################################
############################ TREND ESTIMATION #################################
###############################################################################

################estimate trends for race/ethnicity injury rate#################
fit_time_trend <- function(data_merged) {
  data_merged <- data_merged %>%
    filter(race != "all") 
  data_merged$race <- as.factor(data_merged$race)
  
  #spline
  race_trend_model <- gam(total_count ~ s(year) + race + s(year, by = race),
                          offset = log(pop), family = quasipoisson(),
                          data = data_merged, method = "REML")
  
  #linear trend
  race_linear_model <- gam(total_count ~ year + race + year:race,
                          offset = log(pop), family = quasipoisson(),
                          data = data_merged)
  
  #generate predicted estimates for group-specific time trend (injured per thousand)
  predicted_trend_100k <- predict(race_trend_model,
                                      newdata = data_merged, 
                                      type = "response") * 100000
  
  predicted_linear_100k <- predict(race_linear_model,
                                  newdata = data_merged, 
                                  type = "response") * 100000
  
  # Convert the predicted_trend matrix to a data frame
  predicted_trend_df <- as.data.frame(predicted_trend_100k)
  predicted_linear_df <- as.data.frame(predicted_linear_100k)
  
  # Combine with the original dataframe
  data_merged <- cbind(data_merged, predicted_trend_df)
  data_merged <- cbind(data_merged, predicted_linear_df)
  
  return(data_merged)
}

#predict trend for each simulation of the bootstrapped data
boot_trends <- boot_gamma %>%
  filter(race != "all") %>%
  merge(population, by = c("year", "race")) %>%
  group_by(simulation_id) %>%  # Group by simulation_id
  do({
    # Apply the fit_time_trend function to each subset of data
    fit_time_trend(.)  # This applies the function to the current simulation group
  })

#calculate confidence interval for the smooth trend 
bootstrap_CI <- boot_trends %>%
  group_by(race, year) %>%
  summarise(
    mean_trend = mean(predicted_trend_100k),
    lower_trend = quantile(predicted_trend_100k, 0.025),
    upper_trend = quantile(predicted_trend_100k, 0.975), 
  )

#merge with trend and injury rate from original sample 
yearly_stats_original <- fit_time_trend(merged_raceyear)
bootstrap_CI <- bootstrap_CI %>%
  left_join(yearly_stats_original %>% select(race, year, count_original = total_count, trend_original = predicted_trend_100k, rate_original = injury_rate_100k), by = c("race", "year"))


#estimate trends for all population 
boot_trends_all<- boot_gamma %>%
  merge(population, by = c("year", "race")) %>%
  mutate(injury_rate_100k=total_count/pop*100000)  %>%
  filter(race=="all") %>%
  group_by(simulation_id) %>%
  do({
    trend_model <- gam(total_count ~ s(year),
                            offset = log(pop), family = quasipoisson(),
                            data = ., method = "REML")
    linear_model <- gam(total_count ~ year,
                             offset = log(pop), family = quasipoisson(),
                             data = .)
    predicted_trend_100k <- predict(trend_model,
                                     newdata = ., 
                                     type = "response") * 100000
    predicted_linear_100k <- predict(linear_model,
                                     newdata = ., 
                                     type = "response") * 100000
    predicted_trend_df <- as.data.frame(predicted_trend_100k)
    predicted_linear_df <- as.data.frame(predicted_linear_100k)
    data_merged <- cbind(., predicted_trend_df)
    data_merged <- cbind(data_merged, predicted_linear_df)
    data_merged
  }) %>%
  ungroup()

bootstrap_CI_all <- boot_trends_all %>%
  group_by(year) %>%
  summarise(
    mean_trend = mean(predicted_trend_100k),
    lower_trend = quantile(predicted_trend_100k, 0.025),
    upper_trend = quantile(predicted_trend_100k, 0.975), 
  ) %>%
  ungroup() %>%
  mutate(race="all") 
  
#merge with trend and injury rate from original sample 
yearly_stats_all <- merged_raceyear %>%
  filter(race=="all")
linear_model <- gam(total_count ~ year,
                   offset = log(pop), family = quasipoisson(),
                   data = yearly_stats_all, method = "REML")
trend_model <- gam(total_count ~ s(year),
                   offset = log(pop), family = quasipoisson(),
                   data = yearly_stats_all, method = "REML")
predicted_trend_100k <- predict(trend_model,
                                  newdata = yearly_stats_all, 
                                  type = "response") * 100000
predicted_linear_100k <- predict(linear_model,
                                newdata = yearly_stats_all, 
                                type = "response") * 100000
predicted_linear_df <- as.data.frame(predicted_linear_100k)
predicted_trend_df <- as.data.frame(predicted_trend_100k)
yearly_stats_all<-cbind(yearly_stats_all, predicted_trend_df)
yearly_stats_all<-cbind(yearly_stats_all, predicted_linear_df)

    
bootstrap_CI_all <- bootstrap_CI_all %>%
  left_join(yearly_stats_original %>% select(race, year, count_original = total_count, trend_original = predicted_trend_100k, rate_original = injury_rate_100k), by = c("race", "year")) 

bootstrap_CI_all<-rbind(bootstrap_CI, bootstrap_CI_all)


###################plot smooth trends with confidence intervals#################
#plot estimated trend with CI
bootstrap_CI <- bootstrap_CI_all  %>%
  mutate(race = case_when(
    race == "all"  ~ "All",
    race == "black"  ~ "Black",
    race == "latinx" ~ "Hispanic", 
    race == "white" ~ "White"
  ))
plot_trend<-ggplot(bootstrap_CI %>% filter(race != "other"), aes(x = year, y = trend_original, color = race)) +
  geom_line(size = 1.2) + # Plot the trend line
  geom_point(aes(y = rate_original, shape = "Injured per 100,000"), size = 1.5) + # Plot points for injury rate
  geom_ribbon(aes(ymin = lower_trend, ymax = upper_trend, fill = race), 
              alpha = 0.2, color = NA) + # Plot confidence intervals as shaded area
  labs(
    title = "Yearly Estimated Legal Injury Trend by Race/Ethnicity",
    x = "Year",
    y = "Estimated trend (per 100,000)",
    color = "Race/Ethnicity",
    fill = "Race/Ethnicity",
    shape = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") # Move the legend to the bottom

file_name <- paste0("yearly_trend_parametric_", n_bootstrap, ".png")
ggsave(filename = file_name, plot = plot_trend, width = 8, height = 6, dpi = 300)


####################calculate injury rate ratio and CI#########################
bootstrap_RR <- boot_trends %>%
  mutate(injury_rate_100k = total_count/pop*100000) %>%
  group_by(race, simulation_id) %>%
  summarise(
    mean_rate = mean(injury_rate_100k)
  )

bootstrap_RR <- bootstrap_RR %>%
  pivot_wider(names_from = race, values_from = mean_rate) %>%
  mutate(black_white = black/white, 
         latinx_white = latinx/white) %>%
  summarise(
    mean_rr_bw = mean(black_white),
    lower_rr_bw = quantile(black_white, 0.025),
    upper_rr_bw = quantile(black_white, 0.975), 
    mean_rr_lw = mean(latinx_white),
    lower_rr_lw = quantile(latinx_white, 0.025),
    upper_rr_lw = quantile(latinx_white, 0.975), 
  )

# merge with injury counts and trend from original sample 

rr_original<- yearly_stats_original %>%
  group_by(race) %>%
  summarise(
    mean_rate = mean(injury_rate_100k)
  ) %>%
  pivot_wider(names_from = race, values_from = mean_rate) %>%
  mutate(black_white_original = black/white, 
         latinx_white_original = latinx/white) %>%
  select(black_white_original, latinx_white_original) 

bootstrap_RR<-cbind(bootstrap_RR, rr_original)
file_name <- paste0("RR_parametric_", n_bootstrap, ".csv")
write.csv(bootstrap_RR, file_name, row.names = FALSE)

#calculate CI for rate (ratio) change based on estimated linear trend

boot_linear <- boot_trends %>%
  rbind(boot_trends_all) %>%
  group_by(year, simulation_id) %>%
  summarise(
    Black_Trend = sum(predicted_linear_100k[race == "black"], na.rm = TRUE),
    White_Trend = sum(predicted_linear_100k[race == "white"], na.rm = TRUE),
    Hispanic_Trend = sum(predicted_linear_100k[race == "latinx"], na.rm = TRUE),
    All_Trend = sum(predicted_linear_100k[race == "all"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    black_white_rr = Black_Trend / White_Trend,
    hispanic_white_rr = Hispanic_Trend / White_Trend
  ) %>%
  select(year, simulation_id, black_white_rr, hispanic_white_rr, All_Trend) 

boot_linear_CI <- boot_linear %>% 
  filter(year %in% c(2004, 2021)) %>%
  pivot_wider(
    names_from = year, 
    values_from = c(black_white_rr, hispanic_white_rr, All_Trend)
  ) %>%
  mutate(
    bw_rr_diff = black_white_rr_2021-black_white_rr_2004,
    hw_rr_diff = hispanic_white_rr_2021-hispanic_white_rr_2004,
    all_trend_diff = All_Trend_2021-All_Trend_2004,
    bw_rr_diff_pct = black_white_rr_2021/black_white_rr_2004 -1,
    hw_rr_diff_pct = hispanic_white_rr_2021/hispanic_white_rr_2004 -1,
    all_trend_diff_pct = All_Trend_2021/All_Trend_2004 -1
  ) %>%
  summarise(
    mean_diff_bw = mean(bw_rr_diff),
    lower_diff_bw = quantile(bw_rr_diff, 0.025),
    upper_diff_bw = quantile(bw_rr_diff, 0.975), 
    mean_diff_hw = mean(hw_rr_diff),
    lower_diff_hw = quantile(hw_rr_diff, 0.025),
    upper_diff_hw = quantile(hw_rr_diff, 0.975), 
    mean_diff_all = mean(all_trend_diff),
    lower_diff_all = quantile(all_trend_diff, 0.025),
    upper_diff_all = quantile(all_trend_diff, 0.975),
    mean_diff_pct_bw = mean(bw_rr_diff_pct),
    lower_diff_pct_bw = quantile(bw_rr_diff_pct, 0.025),
    upper_diff_pct_bw = quantile(bw_rr_diff_pct, 0.975), 
    mean_diff_pct_hw = mean(hw_rr_diff_pct),
    lower_diff_pct_hw = quantile(hw_rr_diff_pct, 0.025),
    upper_diff_pct_hw = quantile(hw_rr_diff_pct, 0.975), 
    mean_diff_pct_all = mean(all_trend_diff_pct),
    lower_diff_pct_all = quantile(all_trend_diff_pct, 0.025),
    upper_diff_pct_all = quantile(all_trend_diff_pct, 0.975)
  )

#merge with linear trend RR from original sample 
yearly_stats_original <- fit_time_trend(merged_raceyear)
linear_original<-yearly_stats_original %>%
  rbind(yearly_stats_all) %>%
  group_by(year) %>%
  summarise(
    Black_Trend = sum(predicted_linear_100k[race == "black"], na.rm = TRUE),
    White_Trend = sum(predicted_linear_100k[race == "white"], na.rm = TRUE),
    Hispanic_Trend = sum(predicted_linear_100k[race == "latinx"], na.rm = TRUE),
    All_Trend = sum(predicted_linear_100k[race == "all"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    black_white_rr = Black_Trend / White_Trend,
    hispanic_white_rr = Hispanic_Trend / White_Trend
  ) %>%
  select(year, black_white_rr, hispanic_white_rr, All_Trend) %>%
  filter(year %in% c(2004, 2021)) %>%
  pivot_wider(
    names_from = year, 
    values_from = c(black_white_rr, hispanic_white_rr, All_Trend)
  ) %>%
  mutate(
    bw_rr_diff_pct_original = black_white_rr_2021/black_white_rr_2004 -1,
    hw_rr_diff_pct_original = hispanic_white_rr_2021/hispanic_white_rr_2004 -1,
    all_trend_diff_pct_original = All_Trend_2021/All_Trend_2004 -1,
    bw_rr_diff_original = black_white_rr_2021-black_white_rr_2004,
    hw_rr_diff_original = hispanic_white_rr_2021-hispanic_white_rr_2004,
    all_trend_diff_original = All_Trend_2021-All_Trend_2004 
  ) %>%
  select(bw_rr_diff_original, hw_rr_diff_original, all_trend_diff_original, bw_rr_diff_pct_original, hw_rr_diff_pct_original, all_trend_diff_pct_original) 

boot_linear_RR <- cbind(boot_linear_CI, linear_original)

file_name <- paste0("RR_linearchange_parametric_", n_bootstrap, ".csv")
write.csv(boot_linear_RR, file_name, row.names = FALSE)

