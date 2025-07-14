library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(viridis)
library(cowplot)
library(ggpubr)
library(utils)
library(tidyverse)
library(MetBrewer)
library(ggh4x)
library(data.table)

# Set working directory. 
setwd("/Users/sambents/Desktop/NIH/processed")

###############################################################################################
# color palettes
re_pal <- met.brewer(name="Archambault", n=6)
scenario_pal <- c("gray50","#0A79AA", "#CE4A4E","#6F3D75")

###############################################################################################
#########################################################################################################
# Prerequisite to load data. 
# local path to up-to-date research data repo

# Updated data Feb 2025 - now includes UT and NIH
data_path <- "../processed/data_mar_2025"

# Date used to tag the round (file name - date round id)
round_id_date <- "2024-06-25" #Phase 1 
round_id_date <- "2024-07-16" #Phase 2

# Name of the ensemble(s) to exclude from the analysis
ens_to_excl <- c("Ensemble", "Ensemble_LOP")

# Target on which to run the analysis (cumulative and incidence version)
target_cum <- "cum death"
target_inc <- "inc death"

# Number of sample use to generate the ensemble
n_sample <- 100

# Max horizon
max_horizon <- 20

obs_data_path <-
  paste0("https://raw.githubusercontent.com/midas-network/",
         "covid19-smh-research/main/target-data/time-series.csv")

# Day - 1 to start observed data
start_obs <- as.Date("2020-11-14")

###########################################################################################################
# Load data 
# Connection to processed data
dc <- arrow::open_dataset(paste0(data_path, "/model-processed/"),
                          partitioning = c("round_id", "model_id", "target",
                                           "location"))
# Create re-code vector for model id
mod_encode <- c(dir(paste0(data_path, "/model-output/")),
                "Ensemble", "Ensemble_LOP", "Ensemble_LOP_untrimmed") %>%
  setNames(c(LETTERS[1:length(dir(paste0(data_path, "/model-output/")))],
             "Ensemble", "Ensemble_LOP", "Ensemble_LOP_untrimmed"), .)


# Data frame of cumulative deaths by race/ethnicity at final time point 
df <- dplyr::filter(dc, output_type == "sample", race_ethnicity != "overall",
                    target == target_cum, round_id == round_id_date,
                    horizon == max_horizon) %>%
  dplyr::collect() %>%
  # dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1,
  #               model_id = mod_encode[model_id]) %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1,
                model_id = model_id) %>%
  dplyr::mutate(overall = sum(value),
                .by = c("model_id", "location", "output_type_id", "scenario_id",
                        "target")) %>%
  dplyr::mutate(model_ratio = value / overall)
head(df)
print(unique(df$model_id))

# Generate enesemble by taking 100 simulations from each grouping. 
ensemble <-
  dplyr::slice_sample(df, n = n_sample, replace = FALSE,
                      by = c("model_id", "location", "race_ethnicity",
                             "scenario_id", "target")) %>%
  dplyr::mutate(model_id = "Ensemble")

# Generative quantiles on incident data 
df_quant <- dplyr::filter(dc, output_type == "quantile",
                          !model_id %in% ens_to_excl,
                          round_id == round_id_date,
                          target == target_inc) %>%
  dplyr::collect() %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1,
                model_id = mod_encode[model_id])

##############################################################################################
# Format observed data
# Time series by loaction and race ethnicity of incident deaths 

obs_death_data <- read.csv(obs_data_path) 
head(obs_death_data)

obs_death_data <- read.csv(obs_data_path) %>%
  dplyr::mutate(obs = observation + min_suppressed,
                date = as.Date(date),
                location = as.numeric(location)) %>%
  dplyr::filter(date > "2020-09-01") %>% # was start_obs
  dplyr::filter(target == "inc death") %>%
  mutate(time_value = date) %>%
  dplyr::select(location, race_ethnicity, time_value, obs)
head(obs_death_data)

gold_standard_data_ts <-
  rbind(obs_death_data, dplyr::summarise(obs_death_data, obs = sum(obs),
                                         .by = c("time_value", "location")) %>%
          dplyr::mutate(race_ethnicity = "overall"))
head(gold_standard_data_ts)

# Cumulative ratios of deaths by race/ethnicity by location
gold_standard_data <- obs_death_data %>%
  dplyr::summarise(obs = sum(obs), .by = c("race_ethnicity", "location")) %>%
  dplyr::mutate(tot = sum(obs), .by = location) %>%
  dplyr::mutate(obs_ratio = obs / tot) %>%
  dplyr::select(-obs, -tot)
head(gold_standard_data)

#############################################################################################
# Figure 1: Introduction ot observed data                                                   #
#############################################################################################

pop_table <- tibble(
  race_ethnicity = c("Overall","Asian", "Black", "Latino", "White", "Other", "Overall","Asian", "Black", "White", "Other"),
  location = c("California","California", "California", "California", "California", "California","North Carolina", "North Carolina", "North Carolina", "North Carolina", "North Carolina"),
  pop = c(39346023,5743983, 2142371, 15380929, 14365145, 1713595,
          10698973, 341052, 2155650, 6497519, 1704752
  )
) %>% as.data.frame()

# Census population data 
## CA
Pop_T_CA <- 39346023
Pop_Asian_CA <- 5743983
Pop_Black_CA <- 2142371
Pop_Latino_CA <- 15380929
Pop_White_CA <- 14365145
Pop_Other_CA <- 1713595

## NC
Pop_T_NC <- 10698973
Pop_Asian_NC <- 341052
Pop_Black_NC <- 2155650
Pop_White_NC <- 6497519
Pop_Other_NC <- 1704752

figure1 = gold_standard_data_ts %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "asian", "Asian"))%>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "black", "Black")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "latino", "Latino")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "other", "Other")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "white", "White")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "overall", "Overall")) %>%
  mutate(location = replace(location, location == 6, "California")) %>%
  mutate(location = replace(location, location == 37, "North Carolina")) %>%
  left_join(pop_table) %>%
  mutate(race_ethnicity = fct_relevel(race_ethnicity, "Overall")) 

fig1a = ggplot(data = figure1, aes(x = time_value, y = 100000*obs/pop, col = race_ethnicity)) +
  geom_line(lwd =1.5) +
  geom_vline(xintercept = as.Date("2020-11-15"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-04-01"), lty = "dashed") +
  facet_wrap(vars(location), scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = re_pal) +
  guides(col=guide_legend(title="Race/ethnicity", nrow = 1))  +
  ylab("Weekly deaths per 100k") +
  xlab("Date") + 
  theme_bw() +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16)) +
  #ggtitle("Observed data") +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(colour = "black", size = 12, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  # guides(col = guide_legend(nrow = 1)) +
  ggtitle("b")
fig1a

# Serology over time
serology_plot <- read_csv("serology_data_complete.csv") %>%
  mutate(location = if_else(location == 6, "California", "North Carolina")) %>%
  filter(median_donation > as.Date("2020-09-13", "%Y-%m-%d")) %>%
  group_by(month, year, location) %>%
  mutate(median_donation = median(median_donation)) %>%
  ungroup() %>%
  group_by(median_donation, location, race_ethnicity) %>%
  summarise(value = mean(value), upper = mean(upper), lower = mean(lower)) %>%
  ungroup() %>%
  mutate(race_ethnicity = str_to_sentence(race_ethnicity)) %>%
  mutate(race_ethnicity = fct_relevel(race_ethnicity, "Overall")) %>%
  filter(median_donation < as.Date("2021-05-13", "%Y-%m-%d")) %>%
  ggplot(aes(x = median_donation, color = race_ethnicity)) + 
  geom_line(aes(y = value), size = 1.5) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = race_ethnicity), alpha = 0.08, color = NA) +
  facet_grid(~location) + 
  #facet_grid(vars(race_ethnicity, location, ncol = 2)) +
  scale_color_manual(values = re_pal) +
  scale_fill_manual(values = re_pal) +
  # scale_color_manual("State",values = state_colors) + 
  #  scale_fill_manual("State",values = state_colors) + 
  theme_bw() + 
  geom_vline(xintercept = as.Date("2020-11-15"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-04-01"), lty = "dashed") +
  ggtitle("a") +
  xlab("Date") + 
  guides(fill=guide_legend(title="Race/ethnicity"))  +
  guides(col=guide_legend(title="Race/ethnicity"))  +
  ylab("Proportion seropositive") + 
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "none",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(colour = "black", size = 12, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.spacing = unit(1, "cm")) 
#  aspect.ratio=1)
serology_plot

# calculate IFR: 
# Relative severity of under 65 compared to over 65: 0.046
# IFR = IFR_under65 * (% under 65) + IFR_over65 * (% over 65)
# IFR_under65 = rel_sev * IFR_over65

## Population percentages by age group
age_pops_table <- tibble(
  race_ethnicity = c("Asian", "Black", "White", "Other", "Latino", 
                     "Asian", "Black", "White", "Other", "Latino",
                     "Asian", "Black", "White", "Other",
                     "Asian", "Black", "White", "Other"),
  age = c(rep("65+", 5), rep("< 65", 5), rep("65+", 4), rep("< 65", 4)),
  location = c(rep("California", 10), rep("North Carolina", 8)),
  perc = c(0.1671,0.1317,0.2075, 0.1054,0.0804,
           0.8329,0.8683,0.7925, 0.8946, 0.9196,
           0.0917, 0.1307, 0.1906, 0.0808,
           0.9083,0.8693, 0.8094,0.9192
  )
) %>% as.data.frame() %>% 
  mutate(race_ethnicity = tolower(race_ethnicity))

## Estimating change in serology proportion between november 2020 to april 2021
## using conservative n = 1500 for CA, n = 500 for NC
serology_file <- read_csv("serology_data_complete.csv") %>% 
  filter((month == 11 & year == 2020) | (month == 4 & year == 2021) ) %>%
  mutate(location = if_else(location == 6, "California", "North Carolina")) %>%
  group_by(location, race_ethnicity, year) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = "year", values_from = "value") %>%
  mutate(mean_prop = (`2021` + `2020`/2),
         denom = if_else(location == "California",2/1500,2/500),
         SE_sqrd = mean_prop*(1-mean_prop)*denom,
         SE = sqrt(SE_sqrd),
         diff = `2021`-`2020`,
         upper = diff + SE,
         lower = diff-SE
  ) %>%
  dplyr::select(location, race_ethnicity, upper, lower, diff) %>%
  mutate(race_ethnicity = str_to_sentence(race_ethnicity)) %>%
  left_join(pop_table, by = c("race_ethnicity", "location")) %>%
  mutate(infs = diff*pop, infs_lower = lower*pop, infs_upper = upper*pop)

## Estimating overall IFR between 11/2020 and 4/2021
get_IFR <- gold_standard_data_ts %>%
  filter(time_value >= "2020-11-15") %>%
  group_by(race_ethnicity, location) %>%
  summarise(obs = sum(obs, na.rm = T)) %>%
  mutate(location = if_else(location == 6, "California", "North Carolina")) %>%
  mutate(race_ethnicity = str_to_sentence(race_ethnicity)) %>%
  left_join(serology_file, by = c("race_ethnicity", "location")) %>%
  group_by(race_ethnicity, location) %>%
  summarise(IFR_value = 100*obs/infs, IFR_upper = 100*obs/infs_lower, IFR_lower = 100*obs/infs_upper) %>%
  ungroup() %>%
  filter(race_ethnicity != "Overall") %>%
  left_join(age_pops_table %>% mutate(race_ethnicity = str_to_sentence(race_ethnicity)), by = c("location", "race_ethnicity")) %>%
  pivot_wider(values_from = perc, names_from = "age") %>%
  mutate(adjust_val = 1/(0.046*`< 65` + `65+`),
         IFR_over65_value = IFR_value*adjust_val,
         IFR_over65_lower = IFR_lower*adjust_val,
         IFR_over65_upper = IFR_upper*adjust_val,
         
         IFR_under65_value = 0.046*IFR_over65_value,
         IFR_under65_lower = 0.046*IFR_over65_lower,
         IFR_under65_upper = 0.046*IFR_over65_upper,
         
         IFR_overall_value = IFR_value,
         IFR_overall_lower = IFR_lower,
         IFR_overall_upper = IFR_upper
  ) %>%
  pivot_longer(9:17) %>%
  separate(name, into = c("Metric", "Age", "Estimate"), "_") %>%
  dplyr::select(race_ethnicity, location, Age, Estimate, value) %>%
  mutate(Age = case_when(
    Age == "over65" ~ "65+",
    Age == "under65" ~ "< 65",
    Age == "overall" ~ "Overall"
  ))

race_IFR = get_IFR %>%
  filter(Estimate == "value", Age == "65+")

fig1b <- get_IFR %>%
  pivot_wider(names_from = Estimate, values_from = value) %>%
  mutate(Age = fct_relevel(Age, c("Overall"))) %>%
  ggplot(aes(x = location, fill = race_ethnicity)) + 
  geom_col(aes(y = value, color = race_ethnicity), position = position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.9),
                width = .4) + 
  scale_color_manual("Race/ethnicity",values = re_pal[2:6]) +
  scale_fill_manual("Race/ethnicity",values = re_pal[2:6]) +
  facet_wrap(~Age, scales = "free") + 
  xlab("State") + 
  ylab("Infection fatality rate (%)") + 
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "none",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(colour = "black", size = 12, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.spacing = unit(0.5, "cm")) +
  ggtitle("c")
fig1b

# Compile Figure 1 
plot_grid(serology_plot,fig1a, fig1b, ncol = 1, rel_heights = c(.7, .9))
# 1100 x 1100 dimensions


#############################################################################################
# Supplementary File:  IFR comparison across models                                         #
#############################################################################################
dc2 <- arrow::open_dataset(paste0(data_path, "/model-processed/"),
                           partitioning = c("round_id", "model_id", "target",
                                            "location"))
# Create re-code vector for model id
mod_encode <- c(dir(paste0(data_path, "/model-output/")),
                "Ensemble", "Ensemble_LOP", "Ensemble_LOP_untrimmed") %>%
  setNames(c(LETTERS[1:length(dir(paste0(data_path, "/model-output/")))],
             "Ensemble", "Ensemble_LOP", "Ensemble_LOP_untrimmed"), .)


# Data frame of cumulative deaths by race/ethnicity at final time point 
df2 <- dplyr::filter(dc2, output_type == "sample", race_ethnicity != "overall",
                     round_id == round_id_date,
                     horizon == max_horizon) %>%
  dplyr::collect() %>%
  # dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1,
  #               model_id = mod_encode[model_id]) %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1) %>%
  dplyr::mutate(overall = sum(value),
                .by = c("model_id", "location", "output_type_id", "scenario_id",
                        "target")) %>%
  dplyr::mutate(model_ratio = value / overall)
head(df)
print(unique(df$model_id))

model_IFRs <- df2 %>% 
  filter(scenario_id == "A-2020-11-15") %>%
  dplyr::select(race_ethnicity, output_type_id, value, model_id, target, location) %>%
  filter(target %in% c("cum inf", "cum death")) %>%
  pivot_wider(names_from = "target", values_from = "value") %>%
  mutate(IFR = 100*`cum death`/`cum inf`) %>%
  group_by(location, race_ethnicity, model_id) %>%
  summarise(
    IFR_value = mean(IFR),
    IFR_lower = quantile(IFR, 0.025),
    IFR_upper = quantile(IFR, 0.975)
  ) %>%
  ungroup() %>%
  mutate(location = if_else(location == 6, "California", "North Carolina")) %>%
  left_join(age_pops_table, by = c("race_ethnicity", "location")) %>%
  pivot_wider(values_from = perc, names_from = "age") %>%
  mutate(adjust_val = 1/(0.046*`< 65` + `65+`),
         IFR_over65_value = IFR_value*adjust_val,
         IFR_over65_lower = IFR_lower*adjust_val,
         IFR_over65_upper = IFR_upper*adjust_val,
         
         IFR_under65_value = 0.046*IFR_over65_value,
         IFR_under65_lower = 0.046*IFR_over65_lower,
         IFR_under65_upper = 0.046*IFR_over65_upper,
         
         IFR_overall_value = IFR_value,
         IFR_overall_lower = IFR_lower,
         IFR_overall_upper = IFR_upper
  ) %>%
  dplyr::select(-IFR_value, -IFR_lower, -IFR_upper, -`65+`, -`< 65`) %>%
  pivot_longer(5:13) %>%
  separate(name, into = c("Metric", "Age", "Estimate"), "_") %>%
  dplyr::select(model_id, race_ethnicity, location, Age, Estimate, value) %>%
  mutate(Age = case_when(
    Age == "over65" ~ "65+",
    Age == "under65" ~ "< 65",
    Age == "overall" ~ "Overall"
  ))

plot_model_IFR <- model_IFRs %>%
  pivot_wider(names_from = Estimate, values_from = value) %>%
  mutate(Age = fct_relevel(Age, c("Overall"))) %>%
  mutate(model_id = replace(model_id, model_id == "cumt-seivrcm", "A"))  %>%
  mutate(model_id = replace(model_id, model_id == "JHU_UNC-flepiMoP", "B"))  %>%
  mutate(model_id = replace(model_id, model_id == "MOBS_NEU-COVACS_SEIR", "C"))  %>%
  mutate(model_id = replace(model_id, model_id == "NIH_UIUC-RIFTcov", "D"))  %>%
  mutate(model_id = replace(model_id, model_id == "UTA-ImmunoSEIRS", "E"))  %>%
  mutate(model_id = replace(model_id, model_id == "UVA-EpiHiper", "F"))  %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "asian", "Asian"))%>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "black", "Black")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "latino", "Latino")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "other", "Other")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "white", "White")) %>%
  ggplot(aes(x = location, fill = race_ethnicity)) + 
  geom_col(position = position_dodge(width = 0.9), aes(y = value, color = race_ethnicity)) + 
  geom_errorbar(position = position_dodge(width = 0.9), aes(ymin = lower, ymax = upper), width = .4) + 
  scale_color_manual("Race/ethnicity",values = re_pal[2:6]) +
  scale_fill_manual("Race/ethnicity",values = re_pal[2:6]) +
  ggh4x::facet_grid2(Age~model_id, scales = "free_y") + 
  xlab("State") + 
  ylab("Infection fatality rate (%)") + 
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(colour = "black", size = 12, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.spacing = unit(0.5, "cm"), 
        aspect.ratio=1)
plot_model_IFR



#############################################################################################
# Figure 2: Phase 2 projections                                                             #
#############################################################################################

################## California 
figure2_plot_ca = ggplot (data = scenario_a_plot_p1 %>% filter(location == "California")) + 
  annotate("rect", xmin = as.Date("2020-09-01"), xmax = as.Date("2021-04-01"), ymin = -Inf, ymax = Inf, alpha = .6, fill = "gray89") +
  geom_ribbon(aes(x = time_value, ymin = low, ymax = up, fill = race_ethnicity), alpha = 0.4) +
  geom_ribbon(aes(x = time_value, ymin = low, ymax = up, fill = race_ethnicity), alpha = 0.4) +
  geom_ribbon(aes(x = time_value, ymin = q25, ymax = q75, fill = race_ethnicity), alpha = 0.6) +
  geom_ribbon(aes(x = time_value, ymin = q25, ymax = q75, fill = race_ethnicity), alpha = 0.6) +
  geom_line(aes(x = time_value, y = obs, col = "Observed"), col = "black", lwd = .9, lty = "dashed") +
  geom_line(data = scenario_a_plot_p1 %>% filter(location == "California"), aes(x = time_value, y = obs, col = "Observed"), col = "black", lwd = .9, lty = "dashed") +
  facet_wrap(vars(race_ethnicity), ncol = 3, scales = "free_y") + 
  theme_bw()  +
  guides(col=guide_legend(title="Race/ethnicity")) +
  guides(fill=guide_legend(title="Race/ethnicity")) +
  labs( x = "Date", y = "Incident deaths")  +
  scale_color_manual(values = re_pal[2:6]) +
  scale_fill_manual(values = re_pal[2:6])  +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        strip.text = element_text(colour = "black", size = 10, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  ggtitle("a")
figure2b_plot_ca 

################## North Carolina
figure2_plot_nc = ggplot (data = scenario_a_plot_p1 %>% filter(location == "North Carolina")) + 
  annotate("rect", xmin = as.Date("2020-09-01"), xmax = as.Date("2021-04-01"), ymin = -Inf, ymax = Inf, alpha = .6, fill = "gray89") +
  geom_ribbon(aes(x = time_value, ymin = low, ymax = up, fill = race_ethnicity), alpha = 0.4) +
  geom_ribbon(aes(x = time_value, ymin = low, ymax = up, fill = race_ethnicity), alpha = 0.4) +
  geom_ribbon(aes(x = time_value, ymin = q25, ymax = q75, fill = race_ethnicity), alpha = 0.6) +
  geom_ribbon(aes(x = time_value, ymin = q25, ymax = q75, fill = race_ethnicity), alpha = 0.6) +
  geom_line(aes(x = time_value, y = obs, col = "Observed"), col = "black", lwd = .9, lty = "dashed") +
  geom_line(data = scenario_a_plot_p1 %>% filter(location == "North Carolina"), aes(x = time_value, y = obs, col = "Observed"), col = "black", lwd = .9, lty = "dashed") +
  facet_wrap(vars(race_ethnicity), ncol  = 3, scales = "free_y") + 
  theme_bw()  +
  guides(col=guide_legend(title="Race/ethnicity")) +
  guides(fill=guide_legend(title="Race/ethnicity")) +
  labs( x = "Date", y = "Incident deaths")  +
  scale_color_manual(values = re_pal[2:6]) +
  scale_fill_manual(values = re_pal[2:6])  +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        strip.text = element_text(colour = "black", size = 10, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  ggtitle("b")
figure2_plot_nc

# Join together
plot_grid(figure2_plot_ca, figure2_plot_nc, ncol = 1, rel_heights = c(.25, .26))

##############################################################################################
# Supplementary File:  Projections for Phase 1                                               #
##############################################################################################
round_id_date <- "2024-06-25" #Phase 1 

df_quantiles_a_p1 <- dplyr::filter(dc, output_type == "quantile", race_ethnicity != "overall",
                                   target ==target_inc, round_id == round_id_date) %>%
  dplyr::collect() %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1,
                model_id = mod_encode[model_id]) %>%
  dplyr::mutate(overall = sum(value),
                .by = c("model_id", "location", "output_type_id", "scenario_id",
                        "target"))  %>%
  mutate(time_value = as.Date(time_value)) %>%
  filter(scenario_id == "A-2020-05-01") # phase 1 
print(unique(df_quantiles_a_p1$model_id))

median_quantile = df_quantiles_a_p1 %>% filter(output_type_id == .5) %>% mutate(mid = value)
upper_quantile = df_quantiles_a_p1 %>% filter(output_type_id == .975) %>% mutate(up = value)
lower_quantile = df_quantiles_a_p1 %>% filter(output_type_id == .025) %>% mutate(low = value)
q75_quantile = df_quantiles_a_p1 %>% filter(output_type_id == .75) %>% mutate(q75 = value)
q25_quantile = df_quantiles_a_p1 %>% filter(output_type_id == .25) %>% mutate(q25 = value)

join = c("origin_date", "scenario_id", 
         "race_ethnicity", "horizon", "output_type",
         "location", "model_id", "time_value")
quantile_scenario_a <- median_quantile %>%
  left_join(upper_quantile, by = join) %>%
  left_join(lower_quantile, by = join) %>%
  left_join(q75_quantile, by = join) %>%
  left_join(q25_quantile , by = join)

# Join with observed data.
scenario_a_plot_p1 = left_join(gold_standard_data_ts, quantile_scenario_a %>% filter(model_id == "Ensemble") , by = c("location", "race_ethnicity", 
                                                                                                                      "time_value") ) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "asian", "Asian"))%>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "black", "Black")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "latino", "Latino")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "other", "Other")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "white", "White")) %>%
  filter(race_ethnicity != "overall") %>%
  mutate(location = replace(location, location == 6, "California")) %>%
  mutate(location = replace(location, location == 37, "North Carolina"))

phase1 = ggplot (data = scenario_a_plot_p1) + 
  annotate("rect", xmin = as.Date("2020-09-01"), xmax = as.Date("2020-11-15"), ymin = -Inf, ymax = Inf, alpha = .6, fill = "gray89") +
  geom_ribbon(aes(x = time_value, ymin = low, ymax = up, fill = race_ethnicity), alpha = 0.4) +
  geom_ribbon(aes(x = time_value, ymin = low, ymax = up, fill = race_ethnicity), alpha = 0.4) +
  geom_ribbon(aes(x = time_value, ymin = q25, ymax = q75, fill = race_ethnicity), alpha = 0.6) +
  geom_ribbon(aes(x = time_value, ymin = q25, ymax = q75, fill = race_ethnicity), alpha = 0.6) +
  geom_line(aes(x = time_value, y = obs, col = "Observed"), col = "black", lwd = .9, lty = "dashed") +
  geom_line(data = scenario_a_plot_p1 %>% filter(time_value > "2020-11-15"), aes(x = time_value, y = obs, col = "Observed"), col = "black", lwd = .9) +
  facet_wrap(location ~ race_ethnicity, nrow = 3, scales = "free_y") + 
  theme_bw()  +
  guides(col=guide_legend(title="Race/ethnicity")) +
  guides(fill=guide_legend(title="Race/ethnicity")) +
  labs( x = "Date", y = "Incident deaths")  +
  scale_color_manual(values = re_pal[2:6]) +
  scale_fill_manual(values = re_pal[2:6])  +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        strip.text = element_text(colour = "black", size = 10, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) 
phase1

##############################################################################################
# Supplementary File:  Individual models (only for phase 2)                                 #
##############################################################################################
round_id_date <- "2024-07-16" #Phase 2
df_quantiles_a <- dplyr::filter(dc, output_type == "quantile", race_ethnicity != "overall",
                                target ==target_inc, round_id == round_id_date) %>%
  dplyr::collect() %>%
  #  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1,
  #              model_id = mod_encode[model_id]) %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1) %>%
  dplyr::mutate(overall = sum(value),
                .by = c("model_id", "location", "output_type_id", "scenario_id",
                        "target"))  %>%
  mutate(time_value = as.Date(time_value)) %>%
  filter(scenario_id == "A-2020-11-15")  # phase 2 
print(unique(df_quantiles_a$model_id))

median_quantile = df_quantiles_a %>% filter(output_type_id == .5) %>% mutate(mid = value)
upper_quantile = df_quantiles_a %>% filter(output_type_id == .975) %>% mutate(up = value)
lower_quantile = df_quantiles_a %>% filter(output_type_id == .025) %>% mutate(low = value)
q75_quantile = df_quantiles_a %>% filter(output_type_id == .75) %>% mutate(q75 = value)
q25_quantile = df_quantiles_a %>% filter(output_type_id == .25) %>% mutate(q25 = value)

join = c("origin_date", "scenario_id", 
         "race_ethnicity", "horizon", "output_type",
         "location", "model_id", "time_value")
quantile_scenario_a <- median_quantile %>%
  left_join(upper_quantile, by = join) %>%
  left_join(lower_quantile, by = join) %>%
  left_join(q75_quantile, by = join) %>%
  left_join(q25_quantile , by = join) %>%
  mutate(model_id = replace(model_id, model_id == "cumt-seivrcm", "A"))  %>%
  mutate(model_id = replace(model_id, model_id == "JHU_UNC-flepiMoP", "B"))  %>%
  mutate(model_id = replace(model_id, model_id == "MOBS_NEU-COVACS_SEIR", "C"))  %>%
  mutate(model_id = replace(model_id, model_id == "NIH_UIUC-RIFTcov", "D"))  %>%
  mutate(model_id = replace(model_id, model_id == "UTA-ImmunoSEIRS", "E"))  %>%
  mutate(model_id = replace(model_id, model_id == "UVA-EpiHiper", "F"))  %>%
  mutate(model_id = replace(model_id, model_id == "Ensemble_LOP_untrimmed", "Ensemble linear opinion pool"))  

model_plots = list()
models = print(unique(quantile_scenario_a$model_id ))

for(i in models) {
  
  # Join with observed data.
  scenario_a_plot_p1 = left_join(gold_standard_data_ts, quantile_scenario_a %>% filter(model_id == i) , by = c("location", "race_ethnicity", 
                                                                                                               "time_value") ) %>%
    mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "asian", "Asian"))%>%
    mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "black", "Black")) %>%
    mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "latino", "Latino")) %>%
    mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "other", "Other")) %>%
    mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "white", "White")) %>%
    # filter(model_id == "Ensemble") %>%
    filter(race_ethnicity != "overall") %>%
    mutate(location = replace(location, location == 6, "California")) %>%
    mutate(location = replace(location, location == 37, "North Carolina"))
  
  re_pal <- met.brewer(name="Archambault", n=6)
  
  figure2a_plot = ggplot (data = scenario_a_plot_p1) + 
    annotate("rect", xmin = as.Date("2020-09-01"), xmax = as.Date("2021-04-01"), ymin = -Inf, ymax = Inf, alpha = .6, fill = "gray89") +
    geom_ribbon(aes(x = time_value, ymin = low, ymax = up, fill = race_ethnicity), alpha = 0.4) +
    geom_ribbon(aes(x = time_value, ymin = low, ymax = up, fill = race_ethnicity), alpha = 0.4) +
    geom_ribbon(aes(x = time_value, ymin = q25, ymax = q75, fill = race_ethnicity), alpha = 0.6) +
    geom_ribbon(aes(x = time_value, ymin = q25, ymax = q75, fill = race_ethnicity), alpha = 0.6) +
    geom_line(aes(x = time_value, y = obs, col = "Observed"), col = "black", lwd = .9, lty = "dashed") +
    geom_line(data = scenario_a_plot_p1, aes(x = time_value, y = obs, col = "Observed"), col = "black", lwd = .9, lty = "dashed") +
    facet_wrap(location ~ race_ethnicity, nrow = 2, scales = "free_y") + 
    theme_bw()  +
    ggtitle(paste("Model:", i)) +
    guides(col=guide_legend(title="Race/ethnicity")) +
    guides(fill=guide_legend(title="Race/ethnicity")) +
    labs( x = "Date", y = "Incident deaths")  +
    scale_color_manual(values = re_pal[2:6]) +
    scale_fill_manual(values = re_pal[2:6])  +
    theme(axis.text = element_text(size = 7, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "black", size = 12.5),
          legend.position = "none",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 8, color = "black"),
          legend.title = element_text(size = 10, color = "black"),
          strip.text = element_text(colour = "black", size = 10, hjust = 0),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill=NA)) 
  
  model_plots[[i]] = figure2a_plot
  
} 

# Plot all
combined_plot <- ggarrange(plotlist = model_plots[1:8], ncol = 1, nrow = length(models))
combined_plot

# Plot individual 
m1 = model_plots[7]
m1

###########################################################################################
# Figure 3a: Differences in cumulative number of deaths between scenarios                 #                             
###########################################################################################
## Take the Vincent average (pull quantiles and average values across them)
vincent <- dplyr::filter(dc, output_type == "quantile", race_ethnicity != "overall",
                         target == target_cum, round_id == round_id_date,
                         horizon == max_horizon) %>%
  dplyr::collect() %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1)  %>%
  filter(model_id != "Ensemble_LOP" & model_id != "Ensemble" & model_id != "Ensemble_LOP_untrimmed") %>%
  filter(output_type_id %in% c(.025, .50, .975)) %>%
  group_by(scenario_id, race_ethnicity, location, output_type_id) %>%
  summarize(across(value, ~mean(.x))) 

## Apply this to Overall population
setlow_v = vincent %>% filter(output_type_id == .025) %>%
  group_by(location, scenario_id) %>%
  mutate(overall_025 = sum(value)) %>%
  dplyr::select(location, overall_025, scenario_id) %>%
  ungroup() %>%
  distinct(location, scenario_id, overall_025)

setmid_v = vincent %>% filter(output_type_id == .5) %>%
  group_by(location, scenario_id) %>%
  mutate(overall_50 = sum(value)) %>%
  dplyr::select(location, overall_50, scenario_id)%>%
  ungroup() %>%
  distinct(location, scenario_id, overall_50)

sethigh_v = vincent %>% filter(output_type_id == .975) %>%
  group_by(location, scenario_id) %>%
  mutate(overall_975 = sum(value)) %>%
  dplyr::select(location, overall_975, scenario_id)%>%
  ungroup() %>%
  distinct(location, scenario_id, overall_975)

full_dat_v_overall = left_join(setlow_v , setmid_v, by = c( "location" , "scenario_id")) %>%
  left_join(sethigh_v, by = c("location", "scenario_id" ) ) %>%
  # mutate(scenario_id = replace(scenario_id, scenario_id == "A-2020-11-15", "Inequity-driven transmission and severity")) %>%
  #  mutate(scenario_id = replace(scenario_id, scenario_id == "B-2020-11-15", "Inequity-mitigated transmission")) %>%
  #  mutate(scenario_id = replace(scenario_id, scenario_id == "C-2020-11-15", "Inequity-mitigated severity")) %>%
  #  mutate(scenario_id = replace(scenario_id, scenario_id == "D-2020-11-15", "Inequity-mitigated transmission and severity"))  %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "A-2020-11-15", "A")) %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "B-2020-11-15", "B")) %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "C-2020-11-15", "C")) %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "D-2020-11-15", "D"))  %>%
  mutate(location = replace(location, location == "6", "California")) %>%
  mutate(location = replace(location, location == "37", "North Carolina"))
head(full_dat_v_overall)

vincent_overall = ggplot(data =full_dat_v_overall) +
  facet_wrap(vars(location), scales = "free", nrow = 1)+
  geom_col(data =full_dat_v_overall, aes(x = scenario_id , y = overall_50, fill = scenario_id),  position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(x = scenario_id, ymin = overall_025, ymax = overall_975), 
                position = position_dodge(width = 1.4), width = .2)  +
  theme_bw() +
  scale_fill_manual(values = scenario_pal[1:4]) +
  theme(
    axis.text.x=element_blank()) + 
  ylab("Cumulative deaths") +
  guides(fill=guide_legend(title="Scenario"), position = "bottom")  +
  theme(legend.position = "bottom")  +
  theme(axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 9, color = "black"),
        #legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "black", size = 12, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))   +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("a")
vincent_overall

###########################################################################################
# Supplementary File: Individual model cumualtive deaths using Vincet method.             #              
###########################################################################################

# model A: "cumt-seivrcm"
# model B: "JHU_UNC-flepiMoP"
# model C: "MOBS_NEU-COVACS_SEIR"
# model D: "NIH_UIUC-RIFTcov"
# model E: "UTA-ImmunoSEIRS"
# model F: "UVA-EpiHiper"

#################################################################################
#### Vincent 
vincent <- dplyr::filter(dc, output_type == "quantile", race_ethnicity != "overall",
                         target == target_cum, round_id == round_id_date,
                         horizon == max_horizon) %>%
  dplyr::collect() %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1)  %>%
  filter(model_id != "Ensemble_LOP" & model_id != "Ensemble" & model_id != "Ensemble_LOP_untrimmed") %>%
  filter(model_id == "cumt-seivrcm") %>% # Choose model of interest 
  filter(output_type_id %in% c(.025, .50, .975)) %>%
  group_by(scenario_id, race_ethnicity, location, output_type_id) %>%
  summarize(across(value, ~mean(.x))) 

setlow_v = vincent %>% filter(output_type_id == .025) %>%
  mutate(quant_025 = value) %>%
  dplyr::select(race_ethnicity, location, quant_025, scenario_id) %>%
  ungroup()
head(setlow_v)
setmid_v = vincent %>% filter(output_type_id == .5) %>%
  mutate(quant_50 = value) %>%
  dplyr::select(race_ethnicity, location, quant_50, scenario_id)%>%
  ungroup()
sethigh_v = vincent %>% filter(output_type_id == .975) %>%
  mutate(quant_975 = value) %>%
  dplyr::select(race_ethnicity, location, quant_975, scenario_id)%>%
  ungroup()

full_dat_v = left_join(setlow_v , setmid_v, by = c("race_ethnicity", "location" , "scenario_id")) %>%
  left_join(sethigh_v, by = c("race_ethnicity", "location", "scenario_id" ) ) %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "A-2020-11-15", "Inequity-driven transmission and severity")) %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "B-2020-11-15", "Inequity-mitigated transmission")) %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "C-2020-11-15", "Inequity-mitigated severity")) %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "D-2020-11-15", "Inequity-mitigated transmission and severity")) %>% # %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "asian", "Asian"))%>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "black", "Black")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "latino", "Latino")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "other", "Other")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "white", "White"))  

vincent_ca = ggplot(data = full_dat_v %>% filter(location == 6) %>% ungroup()) +
  facet_wrap(vars(race_ethnicity), scales = "free" , nrow = 1)+
  geom_col(data =full_dat_v %>% filter(location == 6) , aes(x = scenario_id , y = quant_50, fill = scenario_id),  position = position_dodge(width = 0.9)) +
  geom_errorbar(data =full_dat_v %>% filter(location == 6) , aes(x = scenario_id, ymin = quant_025, ymax = quant_975), 
                position = position_dodge(width = 1.4), width = .2)  +
  theme_bw() +
  scale_fill_manual(values = scenario_pal[1:4]) +
  theme(
    axis.text.x=element_blank()) + ggtitle("Model A: California") +
  ylab("Cumulative deaths") +
  guides(fill=guide_legend(title="scenario"), position = "bottom")  +
  theme(legend.position = "bottom") +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "none",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        # legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "black", size = 12, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
vincent_ca

vincent_nc = ggplot(data =full_dat_v %>% filter(location == 37)) +
  facet_wrap(vars(race_ethnicity), scales = "free" , nrow = 1)+
  geom_col(data =full_dat_v %>% filter(location == 37), aes(x = scenario_id , y = quant_50, fill = scenario_id),  position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(x = scenario_id, ymin = quant_025, ymax = quant_975), 
                position = position_dodge(width = 1.4), width = .2)  +
  theme_bw() +
  scale_fill_manual(values = scenario_pal[1:4]) +
  theme(
    axis.text.x=element_blank()) + ggtitle("Model A: North Carolina") +
  ylab("Cumulative deaths") +
  guides(fill=guide_legend(title="Scenario"), position = "bottom")  +
  theme(legend.position = "bottom")  +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        #  axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        #   legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "black", size = 12, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))   +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
vincent_nc

plot_grid(vincent_ca, vincent_nc, nrow = 2, rel_heights = c(.85, 1))

###########################################################################################
# Figure 3: Percent disease averted by scenario using paired trajectories.               #              
###########################################################################################
# JHU = paired 
# cumt-seivrcm = paired
# MOBS_NEU-COVACS_SEIR = paired 
# NIH_UIUC-RIFTcov = unpaired 
# UTA-ImmunoSEIRS = paired
# UVA- Epi-Hiper = unpaired

# Pull out scenario A cumulative deaths 
a_cumulative_deaths <- dplyr::filter(dc, output_type == "sample", race_ethnicity != "overall",
                                     target == target_cum, round_id == round_id_date,
                                     horizon == max_horizon) %>%
  dplyr::collect() %>% filter(scenario_id == "A-2020-11-15") %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1 ) %>%
  filter(model_id !=  "UVA-EpiHiper") %>% # unpaired 
  mutate(scenario_A_value = value) %>% dplyr::select(-value) 
print(unique(a_cumulative_deaths$model_id))

# Calculate difference between scenarios B-D compared to A. 
df_cumulative_deaths_comparison <- dplyr::filter(dc, output_type == "sample", race_ethnicity != "overall",
                                                 target == target_cum, round_id == round_id_date,
                                                 horizon == max_horizon) %>%
  dplyr::collect() %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1)  %>%
  filter(model_id != "UVA-EpiHiper") %>%
  left_join(a_cumulative_deaths, by  = c("output_type", "origin_date", "race_ethnicity", "horizon", "output_type_id",
                                         "round_id", "model_id", "target", "location", "time_value"
  )) %>%
  mutate(scenario_difference = scenario_A_value - value) %>%
  filter(scenario_id.x != "A-2020-11-15") %>%
  mutate(relative_difference = scenario_difference/scenario_A_value ) %>% # %>%
  mutate(location = replace(location, location == 6, "California")) %>%
  mutate(location = replace(location, location == 37, "North Carolina")) %>%
  mutate(scenario_id.x = replace(scenario_id.x, scenario_id.x == "B-2020-11-15", "B")) %>%
  mutate(scenario_id.x = replace(scenario_id.x, scenario_id.x == "C-2020-11-15", "C")) %>%
  mutate(scenario_id.x = replace(scenario_id.x, scenario_id.x == "D-2020-11-15", "D")) %>% # %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "asian", "Asian"))%>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "black", "Black")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "latino", "Latino")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "other", "Other")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "white", "White"))  %>%
  dplyr::slice_sample(n = n_sample, replace = FALSE,
                      by = c("model_id", "location", "race_ethnicity",
                             "scenario_id.x", "target")) 

# Calculate the quantiles on each model 
quantiles = c(.025, .5, .975) 
quantiles_by_model = df_cumulative_deaths_comparison %>%
  group_by(model_id, race_ethnicity, location, scenario_id.x) %>%
  dplyr::summarize(
    quant_025 = quantile(relative_difference, probs = quantiles[1]),
    quant_50 = quantile(relative_difference, probs = quantiles[2]),
    quant_975 = quantile(relative_difference, probs = quantiles[3])) %>%
  ungroup() %>% 
  group_by(race_ethnicity, location, scenario_id.x) %>%
  mutate(avg_025 = mean(quant_025),
         avg_50 = mean(quant_50),
         avg_975 = mean(quant_975)) %>%
  distinct(race_ethnicity, location, avg_025, avg_50, avg_975 )

ensemble_diff_ca = ggplot(data = quantiles_by_model %>% filter(location == "California"), 
                          aes(x = race_ethnicity, y = avg_50, fill = scenario_id.x)) +
  facet_wrap(vars(race_ethnicity), scales = "free", nrow =1) +
  theme_bw() +
  geom_errorbar(aes(ymin = avg_025, ymax = avg_975), position=position_dodge(.9), width = .2) +
  xlab("Race/ethnicity") +
  ylab("Percent deaths averted") +
  guides(col=guide_legend(title="scenario"), position = "bottom") +
  #scale_fill_viridis(discrete = TRUE, option = "G", begin = .3, end = .9)   +
  scale_color_manual(values = scenario_pal[2:4]) +
  geom_point(position = position_dodge(width = 0.9), aes(col = scenario_id.x), cex = 3)  + 
  theme(legend.position = "bottom") +
  guides(fill = FALSE) + 
  scale_y_continuous(limits =c(-.1,1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "none",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        #  legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "white", size = 12, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))  
ensemble_diff_ca 

ensemble_diff_nc = ggplot(data = quantiles_by_model %>% filter(location == "North Carolina"), 
                          aes(x = race_ethnicity, y = avg_50, fill = scenario_id.x)) +
  facet_wrap(vars(race_ethnicity), scales = "free", nrow =1) +
  theme_bw() +
  geom_errorbar(aes(ymin = avg_025, ymax = avg_975), position=position_dodge(.9), width = .2) +
  xlab("Race/ethnicity") +
  ylab("Percent deaths averted") +
  guides(col=guide_legend(title="Scenario"), position = "bottom") +
  scale_color_manual(values = scenario_pal[2:4]) +
  geom_point(position = position_dodge(width = 0.9), aes(col = scenario_id.x), cex = 3)  + 
  theme(legend.position = "bottom") +
  guides(fill = FALSE) + 
  scale_y_continuous(limits =c(-.1,1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        #  legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "white", size = 12, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))
ensemble_diff_nc 

bottom_fig3= plot_grid(ensemble_diff_ca , ensemble_diff_nc, nrow = 2,
                       rel_heights = c(1.4, 1.8))
bottom_fig3

# Compile Figure 3
plot_grid(vincent_overall, bottom_fig3, ncol = 1, rel_heights = c(1.4, 1.8), labels = c("", "b"))

###########################################################################################
# Supplementary File: Model-specific infection rates                                      #        
###########################################################################################
dc2 <- arrow::open_dataset(paste0(data_path, "/model-processed/"),
                           partitioning = c("round_id", "model_id", "target",
                                            "location"))
# Create re-code vector for model id
mod_encode <- c(dir(paste0(data_path, "/model-output/")),
                "Ensemble", "Ensemble_LOP", "Ensemble_LOP_untrimmed") %>%
  setNames(c(LETTERS[1:length(dir(paste0(data_path, "/model-output/")))],
             "Ensemble", "Ensemble_LOP", "Ensemble_LOP_untrimmed"), .)

# Data frame of cumulative deaths by race/ethnicity at final time point 
df2 <- dplyr::filter(dc2, output_type == "sample", race_ethnicity != "overall",
                     round_id == round_id_date,
                     horizon == max_horizon) %>%
  dplyr::collect() %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1) %>%
  dplyr::mutate(overall = sum(value),
                .by = c("model_id", "location", "output_type_id", "scenario_id",
                        "target")) %>%
  dplyr::mutate(model_ratio = value / overall)
head(df2)

head(age_pops_table)
model_infections <- df2 %>% 
  dplyr::select(race_ethnicity, output_type_id, value, model_id, target, location, scenario_id) %>%
  filter(target %in% c("cum inf", "cum death")) %>%
  pivot_wider(names_from = "target", values_from = "value") %>% 
  mutate(location = replace(location, location == 6, "California")) %>%
  mutate(location = replace(location, location == 37, "North Carolina")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "asian", "Asian"))%>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "black", "Black")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "latino", "Latino")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "other", "Other")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "white", "White"))  %>%
  left_join(pop_table, by = c("race_ethnicity", "location")) %>%
  group_by(location, race_ethnicity, model_id, scenario_id) %>%
  summarise(
    cumf_value = mean(`cum inf`),
    cumf_lower = quantile(`cum inf`, 0.025),
    cumf_upper = quantile(`cum inf`, 0.975)
  ) %>%
  ungroup() %>% 
  left_join(pop_table, by = c("race_ethnicity", "location")) 
head(model_infections)

plot_cum_inf_models <- model_infections %>%
  mutate(model_id = replace(model_id, model_id == "cumt-seivrcm", "A"))  %>%
  mutate(model_id = replace(model_id, model_id == "JHU_UNC-flepiMoP", "B"))  %>%
  mutate(model_id = replace(model_id, model_id == "MOBS_NEU-COVACS_SEIR", "C"))  %>%
  mutate(model_id = replace(model_id, model_id == "NIH_UIUC-RIFTcov", "D"))  %>%
  mutate(model_id = replace(model_id, model_id == "UTA-ImmunoSEIRS", "E"))  %>%
  mutate(model_id = replace(model_id, model_id == "UVA-EpiHiper", "F"))  %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "A-2020-11-15", "A")) %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "B-2020-11-15", "B")) %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "C-2020-11-15", "C")) %>%
  mutate(scenario_id = replace(scenario_id, scenario_id == "D-2020-11-15", "D")) %>% # %>%
  ggplot(aes(x = scenario_id, fill = race_ethnicity)) + 
  geom_col(position = position_dodge(width = 0.9), aes(y = cumf_value/pop, color = race_ethnicity)) + 
  geom_errorbar(position = position_dodge(width = 0.9), aes(ymin = cumf_lower/pop, ymax = cumf_upper/pop), width = .4) + 
  scale_color_manual("Race/ethnicity",values = re_pal[2:6]) +
  scale_fill_manual("Race/ethnicity",values = re_pal[2:6]) +
  facet_wrap(vars(location, model_id), nrow = 4) +
  xlab("State") + 
  ylab("Infection per populations") + 
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(colour = "black", size = 12, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.spacing = unit(0.5, "cm"), 
        aspect.ratio=1)
plot_cum_inf_models

###########################################################################################
# EVALULATION.                                                                           #        
###########################################################################################

###########################################################################################
# Supplementary File: Coverage                                                            #        
###########################################################################################
#data_path <- "../processed/data_feb_2025"
data_path <- "../processed/data_mar_2025"

#round_id_date <- "2024-06-25" #Phase 1 
round_id_date <- "2024-07-16" #Phase 2

# Name of the ensemble(s) to exclude from the analysis
ens_to_excl <- c("Ensemble", "Ensemble_LOP")

# Target on which to run the analysis (cumulative and incidence version)
target_cum <- "cum death"
target_inc <- "inc death"

# Number of sample use to generate the ensemble
n_sample <- 100

# Max horizon
max_horizon <- 20

# Path to observed data (csv format)
obs_data_path <-
  paste0("https://raw.githubusercontent.com/midas-network/",
         "covid19-smh-research/main/target-data/target_data_phase2.csv")

# Day - 1 to start observed data
start_obs <- as.Date("2020-11-14")

###########################################################################################################
# Load data 
# Connection to processed data
dc <- arrow::open_dataset(paste0(data_path, "/model-processed/"),
                          partitioning = c("round_id", "model_id", "target",
                                           "location"))

# Data frame of cumulative deaths by race/ethnicity at final time point 
df_models <- dplyr::filter(dc, output_type == "quantile", race_ethnicity != "overall",
                           target == target_inc, round_id == round_id_date) %>%
  dplyr::collect() %>%
  dplyr::mutate(time_value = as.Date(origin_date) + horizon * 7 - 1) %>%
  #   model_id = mod_encode[model_id]) %>%
  dplyr::mutate(overall = sum(value),
                .by = c("model_id", "location", "output_type_id", "scenario_id",
                        "target")) %>%
  dplyr::mutate(model_ratio = value / overall) %>%
  mutate(time_value = as.Date(time_value)) %>%
  mutate(model_id = replace(model_id, model_id == "cumt-seivrcm", "A"))  %>%
  mutate(model_id = replace(model_id, model_id == "JHU_UNC-flepiMoP", "B"))  %>%
  mutate(model_id = replace(model_id, model_id == "MOBS_NEU-COVACS_SEIR", "C"))  %>%
  mutate(model_id = replace(model_id, model_id == "NIH_UIUC-RIFTcov", "D"))  %>%
  mutate(model_id = replace(model_id, model_id == "UTA-ImmunoSEIRS", "E"))  %>%
  mutate(model_id = replace(model_id, model_id == "UVA-EpiHiper", "F"))  %>%
  mutate(model_id = replace(model_id, model_id == "Ensemble_LOP_untrimmed", "Ensemble LOP"))  %>%
  filter(model_id != "Ensemble_LOP")
head(df_models)

# Time series by loaction and race ethnicity of incident deaths 
obs_death_data <- read.csv(obs_data_path) %>%
  dplyr::mutate(obs = value + min_suppressed,
                date = as.Date(date),
                location = as.numeric(location)) %>%
  dplyr::filter(date > start_obs, target == target_inc) %>%
  dplyr::select(location, race_ethnicity, time_value = date, obs)
head(obs_death_data)

gold_standard_data_ts <-
  rbind(obs_death_data, dplyr::summarise(obs_death_data, obs = sum(obs),
                                         .by = c("time_value", "location")) %>%
          dplyr::mutate(race_ethnicity = "overall")) %>%
  mutate(time_value = as.Date(time_value)) %>%
  ungroup()
head(gold_standard_data_ts)

# add observations to gold star data 
df_gs = left_join( gold_standard_data_ts, df_models, by = c("time_value", "location", "race_ethnicity")) 

# coverage table 
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98)) # find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile") #cov %>% dplyr::rename(quantile = value) -> cov
cov$quantile = round(cov$quantile, 3) 
cov$quantile = as.numeric(cov$quantile)
setDT(cov)

head(df_gs)
cov_modes <- cov[df_gs %>% filter(scenario_id == "A-2020-11-15"), on = .(quantile = output_type_id), allow.cartesian=TRUE] %>%
  .[quantile != 0.5] %>%
  data.table::dcast(location + race_ethnicity + model_id  + time_value + alpha + obs ~ variable, value.var = "value") %>%
  .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>%
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs = NULL)]
head(cov_modes)

# Plot coverage by model against alpha. 
cov_modes %>%
  mutate(location = replace(location, location == 6, "California")) %>%
  mutate(location = replace(location, location == 37, "North Carolina")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "asian", "Asian"))%>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "black", "Black")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "latino", "Latino")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "other", "Other")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "white", "White")) %>%
  group_by(alpha, race_ethnicity, location, model_id) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  ggplot(aes(x = alpha, y = sum_cov, col = race_ethnicity)) +
  geom_line(linewidth = 1) + 
  theme_bw() +
  scale_color_manual(values = re_pal[2:6]) +
  ylab("Coverage") +
  xlab("Alpha value")+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  #facet_grid(vars(location), vars(model_id)) +
  facet_grid(vars(factor(model_id, levels = c("A", "B", "C", "D", "E", "F", "Ensemble", "Ensemble LOP"))), vars(location)) +
  #  scale_color_viridis(discrete = TRUE, option = "G", begin = 0, end = .9) +
  theme(legend.position = "bottom") +
  guides(color =guide_legend(title="Race/ethnicity")) +
  theme(axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 9, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 9),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 9, color = "black"),
        legend.title = element_text(size = 9, color = "black"),
        #  legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "black", size = 10, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))


###########################################################################################
# Supplementary File: Weighted interval score calculation.                                 #        
###########################################################################################
wis <- function(q,v,o, 
                a = 2*c(0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.5), 
                IS_components = FALSE){
  q <- round(q,3)
  # check a
  a <- sapply(a/2, function(i){ifelse(i %in% q, 2*i, NA)})
  a <- sapply(a/2, function(i){ifelse((1-i) %in% q, 2*i, NA)})
  a <- a[!is.na(a)]
  # define weight s.t. WIS approximates CRPS (see Bracher et al.)
  w <- a # this assumes w0 = 1/2
  w[length(w)] <- w[length(w)]/2 # do not double count median
  # prepare inputs
  o <- o[1]
  q <- round(q, 4)
  a <- round(a, 4)
  # lower and upper interval bounds
  l <- sapply(a, function(i){v[q == i/2]})
  u <- sapply(a, function(i){v[q == (1-i/2)]})
  # IS components
  IS <- list(
    disp = u-l,
    underpred = 2/a*(o - u)*ifelse(o > u, 1, 0),
    overpred = 2/a*(l - o)*ifelse(o<l, 1, 0)
  )
  # weight
  IS <- lapply(IS, function(i){w*i})
  # sum
  if(IS_components){
    return(list(IS_disp = (1/length(q)) * sum(IS$disp), 
                IS_underpred = (1/length(q)) * sum(IS$underpred), 
                IS_overpred = (1/length(q)) * sum(IS$overpred),
                WIS = 1/(length(q)) * do.call(sum,IS)))
  }
  else{
    return(1/(length(q)) * do.call(sum,IS))
  }
}

df_gs = df_gs %>%
  mutate(value = as.numeric(value), obs = as.numeric(obs), output_type_id = as.numeric(output_type_id)) %>%
  filter(race_ethnicity != "overall")
head(df_gs)
setDT(df_gs)
print(unique(df_gs$model_id))

scores = df_gs[, wis(output_type_id, value,obs,IS_components = TRUE),
               by=.( model_id, location, race_ethnicity, time_value)] %>%
  group_by(location, race_ethnicity) %>%
  mutate(avg_model = mean(WIS)) %>%
  ungroup() %>%
  group_by(model_id, location, race_ethnicity) %>%
  mutate(avg_winmodel = mean(WIS)) %>%
  mutate(score = avg_winmodel/avg_model)  %>%
  ungroup() %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "asian", "Asian"))%>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "black", "Black")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "latino", "Latino")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "other", "Other")) %>%
  mutate(race_ethnicity = replace(race_ethnicity, race_ethnicity == "white", "White")) 

# CA
ca_wis = ggplot(data = scores %>%
                  filter(location == 6), aes( factor(model_id, levels = c("A", "B", "C", "D", "E", "F", "Ensemble", "Ensemble LOP")), race_ethnicity)) +
  geom_tile(aes(fill = score), colour = "white") +
  scale_fill_gradient2(
    low = "royalblue1",
    mid = "white",
    high = "sienna2",
    midpoint = 1,
    space = "Lab",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  theme_bw() + ylab("Race/ethnicity") + xlab("Model") +
  ggtitle("California") +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 12, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "black", size = 10, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))


# NC
nc_wis = ggplot(data = scores %>%
                  filter(location == 37), aes( factor(model_id, levels = c("A", "B", "C", "D", "E", "F", "Ensemble", "Ensemble LOP")), race_ethnicity)) +
  geom_tile(aes(fill = score), colour = "white") +
  scale_fill_gradient2(
    low = "royalblue1",
    mid = "white",
    high = "sienna2",
    midpoint = 1,
    space = "Lab",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  theme_bw() + ylab("Race/ethnicity") + xlab("Model") +
  ggtitle("North Carolina") +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 12, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "black", size = 10, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))

plot_grid(ca_wis, nc_wis, ncol = 1)

