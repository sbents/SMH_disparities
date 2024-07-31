setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/Year_2024/equity/post")
library(arrow)
library(dplyr)
library(ggplot2)
library(data.table)
library(viridis)

# Disparities Round 1 Evaluation. 

# 1. Distribution by race/ethnicity analysis
###################################################### 
# CA
m1_6 <- arrow::read_parquet("MOBS_NEU-COVACS_SEIR6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "A")
m2_6 = arrow::read_parquet("NIH_UIUC-RIFTcov6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "B")
m3_6 = arrow::read_parquet("USC-SIkJalpha6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "C")
m4_6 = arrow::read_parquet("UTA-ImmunoSEIRS6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "D")
m5_6 = arrow::read_parquet("cumt-seivrcm6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "E")

# NC
m1_37 <- arrow::read_parquet("MOBS_NEU-COVACS_SEIR37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "A")
m2_37 = arrow::read_parquet("NIH_UIUC-RIFTcov37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "B")
m3_37 = arrow::read_parquet("USC-SIkJalpha37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "C")
m4_37 = arrow::read_parquet("UTA-ImmunoSEIRS37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "D")
m5_37 = arrow::read_parquet("cumt-seivrcm37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "E")


# Foin NC and CA models together and calculate cumulative deaths by model, location, race/ethnicity, and sample during last week. 
df_models_cum_ratio = rbind(m1_6, m2_6, m3_6, m4_6, m5_6, 
                            m1_37, m2_37, m3_37, m4_37, m5_37) %>%
  mutate(time_value = as.Date(origin_date) + horizon*7-1 ) %>%
  filter(output_type == "sample") %>%
  mutate(location = as.numeric(location), time_value = as.Date(time_value),
         output_type_id = as.numeric(output_type_id)) %>%
  filter(race_ethnicity != "overall") %>%
  group_by(model, location, output_type_id, race_ethnicity) %>%
  mutate(cumulative = cumsum(value)) %>%
  filter(time_value == max(time_value)) %>%
  ungroup() %>%
  group_by(model, location, output_type_id) %>%
  mutate(overall = sum(cumulative)) %>%
  ungroup() %>%
  mutate(model_ratio = cumulative/overall)

# Generative ratios for ensemble by selecting 100 samples from each model grouping. 
ensemble = df_models_cum_ratio %>%
  group_by(model, location, race_ethnicity) %>%
  slice_sample(n = 100, replace = FALSE)

# Coverage table 
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98)) # find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile") #cov %>% dplyr::rename(quantile = value) -> cov
cov$quantile = round(cov$quantile, 3) 
cov$quantile = as.numeric(cov$quantile)
setDT(cov)
# Quantile set 
quantile_bins <- print(cov$quantile)
print(quantile_bins)
qtib = tibble(q1=0.55,q2=0.60,q3=0.65,q4=0.7,q5=0.75, q6=0.8,q7= 0.850,
              q8=0.9,q9=0.95,q10= 0.975, q11=0.99,q13= 0.45,  q14=0.400,q15= 0.350,q16= 0.300,q17=0.250,q18= 0.200,
              q19=0.15,q20=0.1,q21= 0.05,q22= 0.025,q23= 0.011)

# Calculate quantiles on ratios by model, location, and race/ethnicity. 
mod_ratio <- df_models_cum_ratio %>%
  group_by(model, location, race_ethnicity) %>%
  summarise(model_cum_ratio =quantile(model_ratio, qtib),
            quantile = quantile_bins)%>%
  ungroup() 
# Calculate quantiles on ensemble. 
ens_ratio = ensemble %>%
  group_by(location, race_ethnicity) %>%
  summarise(model_cum_ratio =quantile(model_ratio, qtib),
            quantile = quantile_bins)%>%
  ungroup()  %>%
  mutate(model ="Ensemble")

# Join models and ensemble quantile datasets. 
full_ratio = rbind(mod_ratio, ens_ratio)

## Load gold standard data. 
gold_standard_data = read.csv("target_data_phase2.csv") %>%
  mutate(eval = value + min_suppressed) %>%
  filter(date > "2020-11-14" ) %>%
  filter(target == "inc death") %>%
  mutate(time_value = as.Date(date), location = as.numeric(location)) %>%
  mutate(obs = eval) %>% 
  dplyr::select(location, race_ethnicity, time_value, obs)
overall = gold_standard_data %>%
  group_by(location, time_value) %>%
  summarize(across(obs, ~sum(.x, na.rm = TRUE))) %>%
  mutate(race_ethnicity = "overall")
gs_cum = rbind(gold_standard_data, overall) %>%
  ungroup() %>%
  mutate(obs = as.numeric(obs), time_value = as.Date(time_value)) %>%
  group_by(race_ethnicity, location) %>%
  mutate(cum = cumsum(obs)) %>%
  filter(time_value == max(time_value)) %>%
  mutate(obs_ratio = ifelse(location == 6, cum/41937 , cum/7850)) %>%
  filter(race_ethnicity != "overall")

# Join model data and gold standard data. 
df_gs_ratio = left_join(full_ratio, gs_cum, by = c("location", "race_ethnicity")) %>%
  dplyr::select(-obs, -cum) 
setDT(df_gs_ratio)

# Plot states.
## NC
ggplot(data = df_gs_ratio %>% filter(location == 37))  + 
  geom_point(aes(x = model_cum_ratio, y = quantile,col = race_ethnicity), cex =2) +
  facet_grid(vars(model), vars(race_ethnicity), scales  = "free") + 
  xlab("Proportion") +
  ylab("Quantile") + theme_bw() +
  ggtitle("North Carolina") +
  geom_vline(aes(xintercept = obs_ratio), lwd = 1.2)  +
  scale_color_viridis(discrete = TRUE, option = "H", begin = .15, end = .85)  +
  labs(col ='Race/ethnicity')

## CA
ggplot(data = df_gs_ratio %>% filter(location == 6))  + 
  geom_point(aes(x = model_cum_ratio, y = quantile,col = race_ethnicity), cex =2) +
  facet_grid(vars(model), vars(race_ethnicity), scales  = "free") + 
  xlab("Proportion") +
  ylab("Quantile") + theme_bw() +
  ggtitle("California") +
  geom_vline(aes(xintercept = obs_ratio), lwd = 1.2)  +
  scale_color_viridis(discrete = TRUE, option = "G", begin = .2, end = .8)  +
  labs(col ='Race/ethnicity') 

# Add in null expectation 
race_ethnicity = c("white", "black", "asian", "latino", "other", "white", "asian", "black", "other")
location = c(6, 6, 6, 6, 6,  37, 37, 37, 37)
prop = c(14365145/39346023, 2142371/39346023, 5743983/39346023, 15380929/39346023, 1713595/39346023,6497519/10698973, 341052/10698973, 2155650/10698973, 1704752/10698973)
tab_race = data.frame(race_ethnicity, location, prop)

null = left_join(df_gs_ratio, tab_race, by = c("race_ethnicity", "location"))

## NC
ggplot(data = null %>% filter(location == 37))  + 
  geom_point(aes(x = model_cum_ratio, y = quantile,col = race_ethnicity), cex =2) +
  facet_grid(vars(model), vars(race_ethnicity), scales  = "free") + 
  xlab("Proportion") +
  ylab("Quantile") + theme_bw() +
  ggtitle("North Carolina") +
  geom_vline(aes(xintercept = obs_ratio), lwd = 1.2)  +
  geom_vline(aes(xintercept = prop), col = "red", alpha = .7, lwd = 1.2, lty = "dashed")  +
  scale_color_viridis(discrete = TRUE, option = "G", begin = .2, end = .8)  +
  labs(col ='Race/ethnicity')

## CA
ggplot(data = null %>% filter(location == 6))  + 
  geom_vline(aes(xintercept = prop), col = "red", alpha = .7, lwd = 1.2, lty = "dashed")  +
  geom_point(aes(x = model_cum_ratio, y = quantile,col = race_ethnicity), cex =2) +
  facet_grid(vars(model), vars(race_ethnicity), scales  = "free") + 
  xlab("Proportion") +
  ylab("Quantile") + theme_bw() +
  ggtitle("California") +
  geom_vline(aes(xintercept = obs_ratio), lwd = 1.2)  +
  scale_color_viridis(discrete = TRUE, option = "G", begin = .2, end = .8)  +
  labs(col ='Race/ethnicity') 


# Check ensemble 95% CI.
ggplot() +
  geom_point(data = df_gs_ratio %>% filter(model == "Ensemble"), aes(x = race_ethnicity, y = obs_ratio), col = "red", cex = 4, shape = 23) +
  geom_point(data = df_gs_ratio %>% filter(model == "Ensemble") %>% filter(quantile == .975 | quantile == .025 | quantile == .5)
             , aes(x = race_ethnicity, y = model_cum_ratio, col = as.character(quantile)), cex = 4) +
  facet_wrap(vars(location)) +
  scale_color_manual(values = c("darkslategray2", "darkslategray3","darkslategray4")) +
  theme_bw()  +
  guides(col=guide_legend(title="Quantile"))

# Calculate coverage 
cov_ratio <- cov[df_gs_ratio, on = .(quantile = quantile), allow.cartesian=TRUE] %>%
  data.table::dcast(location + model +  race_ethnicity  + alpha + obs_ratio ~ variable, value.var = "model_cum_ratio") %>%
  .[, ":=" (cov = ifelse(obs_ratio < upr & obs_ratio > lwr, 1, 0))] %>%
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs_ratio = NULL)]

# Coverage by race/ethnicity and location. 
cov_ratio %>%
  group_by(alpha, location, race_ethnicity) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  ggplot(aes(x = alpha, y = sum_cov), col = "orchid") +
  geom_line(lwd = 3, col = "orchid") +
  theme_bw() +
  ylab("Coverage") +
  xlab("Alpha value")+
  facet_grid(vars(location), vars(race_ethnicity)) +
  ggtitle("Ratio of cumulative race:overall")

# Save weighted interval score calculation. 
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


# Calculate WIS by location and model. 
#scores_cum <-  df_gs_ratio[, wis(quantile,model_cum_ratio,obs_ratio,IS_components = TRUE),
#                           by=.( model, location)] %>%
#  group_by(location) %>%
#  mutate(avg_model = mean(WIS)) %>%
#  mutate(score = WIS/avg_model)#

#ggplot(data = scores_cum, aes( model, as.character(location ))) +
#  geom_tile(aes(fill = score), colour = "white") +
#  scale_fill_gradient2(
#    low = "royalblue1",
#    mid = "white",
#    high = "sienna2",
#    midpoint = 1,
#    space = "Lab",
#    guide = "colourbar",
#    aesthetics = "fill"
#  ) +
#  theme_bw() + ylab("Location") + xlab("Model") +
#  ggtitle("WIS by location and model compared to average")

head(df_gs_ratio)
#Calculate WIS by location and model for race/ethnicity. 
scores_cum <-  df_gs_ratio[, wis(quantile,model_cum_ratio,obs_ratio,IS_components = TRUE),
                           by=.( model, location, race_ethnicity)] %>%
  group_by(location, race_ethnicity) %>%
  mutate(avg_model = mean(WIS)) %>%
  mutate(score = WIS/avg_model)

# CA
ggplot(data = scores_cum %>% filter(location == 6), aes( model, as.character(race_ethnicity))) +
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
  theme_bw() + ylab("Location") + xlab("Model") +
  ggtitle("CA: WIS by r/e and model compared to average")

# NC
ggplot(data = scores_cum %>% filter(location == 37), aes( model, as.character(race_ethnicity))) +
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
  theme_bw() + ylab("Location") + xlab("Model") +
  ggtitle("NC: WIS by r/e and model compared to average")

##############################################################################
##################################################################################
##################################################################################
# 2. Coverage and WIS on deaths across all weeks. 

# CA
df_ensemble6 <- arrow::read_parquet("ELOP_UT6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "Ensemble")
# NC
df_ensemble37 <- arrow::read_parquet("ELOP_UT37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "Ensemble")

# Join NC and CA. 
df_models = rbind(df_ensemble6, df_ensemble37, m1_6, m2_6, m3_6, m4_6, m5_6, 
                  m1_37, m2_37, m3_37, m4_37, m5_37) %>%
  mutate(time_value = as.Date(origin_date) + horizon*7-1 ) %>%
  filter(output_type == "quantile") %>%
  mutate(location = as.numeric(location), time_value = as.Date(time_value),
         output_type_id = as.numeric(output_type_id))

# truth data, value = obs 
gold_standard_data = read.csv("target_data_phase2.csv") %>%
  mutate(eval = value + min_suppressed) %>%
  filter(date > "2020-11-14" ) %>%
  filter(target == "inc death") %>%
  mutate(time_value = as.Date(date), location = as.numeric(location)) %>%
  mutate(obs = eval) %>% 
  dplyr::select(location, race_ethnicity, time_value, obs)
# calculate for overall 
overall = gold_standard_data %>%
  group_by(location, time_value) %>%
  summarize(across(obs, ~sum(.x, na.rm = TRUE))) %>%
  mutate(race_ethnicity = "overall")
gold_standard_data = rbind(gold_standard_data, overall) %>%
  mutate(obs = as.numeric(obs))
table(gold_standard_data$race_ethnicity)


# add observations to gold star data 
df_gs = left_join(df_models, gold_standard_data, by = c("time_value", "location", "race_ethnicity")) 
setDT(df_gs)
table(df_gs$time_value)
table(df_gs$model)

# calcuate MAE for 50% median 
head(df_gs)
mae = df_gs %>%
  ungroup() %>%
  filter(output_type_id == .5) %>%
  mutate(error = abs(value - obs)) %>%
  group_by(model, location, race_ethnicity) %>%
  mutate(sum_error = sum(error)) %>%
  mutate(mae = sum_error/length(unique(time_value))) %>%
  distinct(model, mae, location, race_ethnicity) %>%
  mutate(location = replace(location , location == 6, "California")) %>%
  mutate(location = replace(location , location == 37, "North Carolina"))
head(mae)

ggplot(data = mae, aes(x = race_ethnicity, y = log(mae), fill = model)) +
  geom_col( position = "dodge" )+
  facet_wrap(vars(location), scales = "free_y", ncol = 1) +
  scale_fill_viridis(discrete = TRUE, option = "H", begin = .2, end = 1, alpha = .8) +
  theme_bw() +
  ylab("Log mean absolute error") +
  xlab("Race/ethnicity")

mae_6 = ggplot(data = mae %>% filter(location == "California"), aes(x = race_ethnicity, y = mae, fill = model)) +
  geom_col( position = "dodge" )+
  facet_wrap(vars(race_ethnicity), scales = "free", nrow = 1) +
  scale_fill_viridis(discrete = TRUE, option = "H", begin = .2, end = 1, alpha = .8) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab("Race/ethnicity") +
  ggtitle("Calfornia") +
  theme(legend.position="none")

mae_37 = ggplot(data = mae %>% filter(location == "North Carolina"), aes(x = race_ethnicity, y = mae, fill = model)) +
  geom_col( position = "dodge" )+
  facet_wrap(vars(race_ethnicity), scales = "free", nrow = 1) +
  scale_fill_viridis(discrete = TRUE, option = "H", begin = .2, end = 1, alpha = .8) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab("Race/ethnicity")  +
  ggtitle("North Carolina")
mae_plot = plot_grid(mae_6, mae_37, nrow = 2)
mae_plot          

# coverage table 
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98)) # find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile") #cov %>% dplyr::rename(quantile = value) -> cov
cov$quantile = round(cov$quantile, 3) 
cov$quantile = as.numeric(cov$quantile)
setDT(cov)


cov_modes <- cov[df_gs, on = .(quantile = output_type_id), allow.cartesian=TRUE] %>%
  .[quantile != 0.5] %>%
  data.table::dcast(location + race_ethnicity + model  + time_value + alpha + obs ~ variable, value.var = "value") %>%
  .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>%
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs = NULL)]

# Plot coverage by model against alpha. 
cov_modes %>%
  mutate(location = replace(location, location == 6, "California")) %>%
  mutate(location = replace(location, location == 37, "North Carolina")) %>%
  group_by(alpha, race_ethnicity, location, model) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  ggplot(aes(x = alpha, y = sum_cov, col = race_ethnicity)) +
  geom_line(linewidth = 1) + 
  theme_bw() +
  ylab("Coverage") +
  xlab("Alpha value")+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  facet_grid(vars(model), vars(location)) +
  scale_color_viridis(discrete = TRUE, option = "H", begin = .2, end = .9) 

# Plot coverage by race against alpha. 
cov_modes %>%
  mutate(location = replace(location, location == 6, "California")) %>%
  mutate(location = replace(location, location == 37, "North Carolina")) %>%
  group_by(alpha, race_ethnicity, location, model) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  ggplot(aes(x = alpha, y = sum_cov, col = model)) +
  geom_line(linewidth = 1) + 
  theme_bw() +
  ylab("Coverage") +
  xlab("Alpha value")+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  facet_grid(vars(race_ethnicity), vars(location)) +
  scale_color_viridis(discrete = TRUE, option = "H", begin = .2, end = .9) 

# Calculate WIS by race ethnicity and model 
#scores <-  df_gs[, wis(output_type_id,value,obs,IS_components = TRUE),
#                  by=.( model, location, race_ethnicity, time_value)] %>%
#  group_by(location, race_ethnicity) %>%
#  mutate(avg_model = mean(WIS)) %>%
#  mutate(score = WIS/avg_model) %>%
#  mutate(avg_winmodel = mean(WIS), .by = "model") %>%
#  mutate(score = avg_winmodel/avg_model)
head(df_gs)
table(df_gs$time_value)

scores = df_gs[, wis(output_type_id,value,obs,IS_components = TRUE),
               by=.( model, location, race_ethnicity, time_value)] %>%
  group_by(location, race_ethnicity) %>%
  mutate(avg_model = mean(WIS)) %>%
  mutate(score = WIS/avg_model)
head(scores)


scores = df_gs[, wis(output_type_id,value,obs,IS_components = TRUE),
      by=.( model, location, race_ethnicity, time_value)] %>%
  group_by(location, race_ethnicity) %>%
  mutate(avg_model = mean(WIS)) %>%
  ungroup() %>%
  group_by(model, location, race_ethnicity) %>%
  mutate(avg_winmodel = mean(WIS)) %>%
  mutate(score = avg_winmodel/avg_model)


# CA
ggplot(data = scores %>%
         filter(location == 6), aes( model, race_ethnicity)) +
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
  ggtitle("California")


# NC
ggplot(data = scores %>%
         filter(location == 37), aes( model, race_ethnicity)) +
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
  ggtitle("North Carolina: WIS by r/e and model compared to average")



##################### plot overall deaths 


# CA
df_ensemble6 <- arrow::read_parquet("ELOP_UT6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "Ensemble")
m1_6 <- arrow::read_parquet("MOBS_NEU-COVACS_SEIR6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "A")
m2_6 = arrow::read_parquet("NIH_UIUC-RIFTcov6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "B")
m3_6 = arrow::read_parquet("USC-SIkJalpha6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "C")
m4_6 = arrow::read_parquet("UTA-ImmunoSEIRS6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "D")
m5_6 = arrow::read_parquet("cumt-seivrcm6.gz.parquet") %>% mutate(location = 6) %>% mutate(model = "E")

# NC
df_ensemble37 <- arrow::read_parquet("ELOP_UT37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "Ensemble")
m1_37 <- arrow::read_parquet("MOBS_NEU-COVACS_SEIR37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "A")
m2_37 = arrow::read_parquet("NIH_UIUC-RIFTcov37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "B")
m3_37 = arrow::read_parquet("USC-SIkJalpha37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "C")
m4_37 = arrow::read_parquet("UTA-ImmunoSEIRS37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "D")
m5_37 = arrow::read_parquet("cumt-seivrcm37.gz.parquet") %>% mutate(location = 37) %>% mutate(model = "E")

# join NC and CA
df_models = rbind(df_ensemble6, df_ensemble37, m1_6, m2_6, m3_6, m4_6, m5_6, 
                  m1_37, m2_37, m3_37, m4_37, m5_37) %>%
  mutate(time_value = as.Date(origin_date) + horizon*7-1 ) %>%
  filter(output_type == "quantile") %>%
  mutate(location = as.numeric(location), time_value = as.Date(time_value),
         output_type_id = as.numeric(output_type_id))
table(df_models$model)

# truth data, value = obs 
gold_standard_data = read.csv("target_data_phase2.csv") %>%
  mutate(eval = value + min_suppressed) %>%
  filter(date > "2020-11-14" ) %>%
  filter(target == "inc death") %>%
  mutate(time_value = as.Date(date), location = as.numeric(location)) %>%
  mutate(obs = eval) %>% 
  dplyr::select(location, race_ethnicity, time_value, obs)
# calculate for overall 
overall = gold_standard_data %>%
  group_by(location, time_value) %>%
  summarize(across(obs, ~sum(.x, na.rm = TRUE))) %>%
  mutate(race_ethnicity = "overall")
gold_standard_data = rbind(gold_standard_data, overall) %>%
  mutate(obs = as.numeric(obs))
table(gold_standard_data$race_ethnicity)


# add observations to gold star data 
df_gs = left_join(df_models, gold_standard_data, by = c("time_value", "location", "race_ethnicity")) 
setDT(df_gs)
head(df_gs)

### 

mid = df_gs %>% filter(output_type_id == .5) %>% mutate(mid = value)
up = df_gs %>% filter(output_type_id == .975) %>% mutate(up = value)
low = df_gs %>% filter(output_type_id == .025) %>% mutate(low = value)

hey = left_join(mid, up, by = c("origin_date", "scenario_id", 
                                "race_ethnicity", "horizon", "output_type",
                                "location", "model", "time_value", "obs"))
hey1 = left_join(hey, low, by = c("origin_date", "scenario_id", 
                                 "race_ethnicity", "horizon", "output_type",
                                 "location", "model", "time_value", "obs"))

ggplot(hey1 %>% filter(location == 6), aes(x = time_value, y = mid, col = race_ethnicity)) +
  geom_line(lwd = 1.25, alpha = .8) +
  geom_line(aes(x = time_value, y = obs), col = "black") +
  geom_ribbon(aes(ymin = low, ymax = up, fill = race_ethnicity), alpha = 0.2) +
  labs(
       x = "Date",
       y = "Incident daeths") +
  facet_grid(vars(race_ethnicity), vars(model), scales = "free_y") + 
  theme_bw() +
  scale_color_viridis(discrete = TRUE, option = "H", begin = .3, end = .9) +
  scale_fill_viridis(discrete = TRUE, option = "H", begin = .3, end = .9) 

ggplot(hey1 %>% filter(location == 37), aes(x = time_value, y = mid, col = race_ethnicity)) +
  geom_line(lwd = 1.25, alpha = .8) +
  geom_line(aes(x = time_value, y = obs), col = "black") +
  geom_ribbon(aes(ymin = low, ymax = up, fill = race_ethnicity), alpha = 0.2) +
  labs(
    x = "Date",
    y = "Incident daeths") +
  facet_grid(vars(race_ethnicity), vars(model), scales = "free_y") + 
  theme_bw() +
  scale_color_viridis(discrete = TRUE, option = "H", begin = .3, end = .9) +
  scale_fill_viridis(discrete = TRUE, option = "H", begin = .3, end = .9) 


## Density plots #####################################################


density = df_models_cum_ratio  %>% filter(location == 37) %>% 
  left_join(gs_cum, by = c("race_ethnicity", "location"))
head(density)

a = ggplot(data = density %>% filter(model == "A"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(data = density %>% filter(model == "A") %>% group_by(race_ethnicity) %>%
               mutate(mean = mean(model_ratio)),
             aes(xintercept = mean), lty = "dashed", lwd= 1, col = "purple" ) +
  geom_vline(aes(xintercept = obs_ratio), col = "black", lwd= 1 ) +
  facet_wrap(vars(race_ethnicity), scales = "free", nrow = 1) +
  theme_bw() + ggtitle("A")   +
  xlab("Proportion")

b = ggplot(data = density %>% filter(model == "B"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(data = density %>% filter(model == "B") %>% group_by(race_ethnicity) %>%
               mutate(mean = mean(model_ratio)),
             aes(xintercept = mean), lty = "dashed", lwd= 1, col = "purple" ) +
  geom_vline(aes(xintercept = obs_ratio), col = "black",  lwd= 1 ) +
  facet_wrap(vars(race_ethnicity), scales = "free", nrow = 1) +
  theme_bw() + ggtitle("B")  +   xlab("Proportion")
c = ggplot(data = density %>% filter(model == "C"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(data = density %>% filter(model == "C") %>% group_by(race_ethnicity) %>%
               mutate(mean = mean(model_ratio)),
             aes(xintercept = mean), lty = "dashed", lwd= 1, col = "purple" ) +
  geom_vline(aes(xintercept = obs_ratio), col = "black", lwd= 1 ) +
  facet_wrap(vars(race_ethnicity), scales = "free", nrow = 1) +
  theme_bw() + ggtitle("C")  +   xlab("Proportion")

d = ggplot(data = density %>% filter(model == "D"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(data = density %>% filter(model == "D") %>% group_by(race_ethnicity) %>%
               mutate(mean = mean(model_ratio)),
             aes(xintercept = mean), lty = "dashed", lwd= 1, col = "purple" ) +
  geom_vline(aes(xintercept = obs_ratio),col = "black", lwd= 1 ) +
  facet_wrap(vars(race_ethnicity), scales = "free", nrow = 1) +
  theme_bw() + ggtitle("D")   +   xlab("Proportion")


e = ggplot(data = density %>% filter(model == "E"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(data = density %>% filter(model == "E") %>% group_by(race_ethnicity) %>%
               mutate(mean = mean(model_ratio)),
             aes(xintercept = mean), lty = "dashed", lwd= 1, col = "purple" ) +
  geom_vline(aes(xintercept = obs_ratio),col = "black", lwd= 1 ) +
  facet_wrap(vars(race_ethnicity), scales = "free", nrow = 1) +
  theme_bw() + ggtitle("E")  +   xlab("Proportion")


plot_grid(a, b, c, d, e, ncol = 1)

## Second plot for Justin. 

# Add in null expectation 
race_ethnicity = c("white", "black", "asian", "latino", "other", "white", "asian", "black", "other")
location = c(6, 6, 6, 6, 6,  37, 37, 37, 37)
prop = c(14365145/39346023, 2142371/39346023, 5743983/39346023, 15380929/39346023, 1713595/39346023,6497519/10698973, 341052/10698973, 2155650/10698973, 1704752/10698973)
tab_race = data.frame(race_ethnicity, location, prop)

null = left_join(df_models_cum_ratio, tab_race, by = c("race_ethnicity", "location"))

density = null  %>% filter(location == 6) %>% 
  left_join(gs_cum, by = c("race_ethnicity", "location"))

a = ggplot(data = density %>% filter(race_ethnicity == "white"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio), col = "black", lwd= 1 ) +
  facet_wrap(vars(model), scales = "free_y", ncol = 1) +
  theme_bw() + ggtitle("White")   +
  xlab("Proportion") +
  xlim(c(.1, .75))
a

b = ggplot(data = density %>% filter(race_ethnicity == "black"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio), col = "black",  lwd= 1 ) +
  facet_wrap(vars(model), scales = "free_y", ncol = 1) +
  theme_bw() + ggtitle("Black")  +   xlab("Proportion") +
  xlim(c(0, .15))
b

c = ggplot(data = density %>% filter(race_ethnicity == "asian"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio), col = "black", lwd= 1 ) +
  facet_wrap(vars(model), scales = "free_y", ncol = 1) +
  theme_bw() + ggtitle("Asian")  +   xlab("Proportion")+
 # xlim(c(0, 1))
  xlim(c(0, .36))
c

d = ggplot(data = density %>% filter(race_ethnicity == "other"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio),col = "black", lwd= 1 ) +
  facet_wrap(vars(model), scales = "free_y", ncol = 1) +
  theme_bw() + ggtitle("Other")   +   xlab("Proportion")+
  xlim(c(0,  .3))
d

e = ggplot(data = density %>% filter(race_ethnicity == "latino"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio),col = "black", lwd= 1 ) +
  facet_wrap(vars(model), scales = "free_y", ncol = 1) +
  theme_bw() + ggtitle("Latino")  +   xlab("Proportion")+
  xlim(c(.1, .75))
e

plot_grid(a, b, c, d, e, ncol = 5)

# 1400 x 1100


# NC
density = null  %>% filter(location == 37) %>% 
  left_join(gs_cum, by = c("race_ethnicity", "location"))

a = ggplot(data = density %>% filter(race_ethnicity == "white"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio), col = "black", lwd= 1 ) +
  facet_wrap(vars(model), scales = "free_y", ncol = 1) +
  theme_bw() + ggtitle("White")   +
  xlab("Proportion") +
  xlim(c(.5, 1))
a

b = ggplot(data = density %>% filter(race_ethnicity == "black"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio), col = "black",  lwd= 1 ) +
  facet_wrap(vars(model), scales = "free_y", ncol = 1) +
  theme_bw() + ggtitle("Black")  +   xlab("Proportion") +
  xlim(c(0, .4))
b

c = ggplot(data = density %>% filter(race_ethnicity == "asian"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio), col = "black", lwd= 1 ) +
  facet_wrap(vars(model), scales = "free_y", ncol = 1) +
  theme_bw() + ggtitle("Asian")  +   xlab("Proportion")+
  # xlim(c(0, 1))
  xlim(c(0, .1))
c

d = ggplot(data = density %>% filter(race_ethnicity == "other"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio),col = "black", lwd= 1 ) +
  facet_wrap(vars(model), scales = "free_y", ncol = 1) +
  theme_bw() + ggtitle("Other")   +   xlab("Proportion")+
  xlim(c(0,  .18))
d


plot_grid(a, b, c, d, ncol = 4)












####### dont need
a = ggplot(data = density %>% filter(model == "A"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio), col = "black", lwd= 1 ) +
  facet_wrap(vars(race_ethnicity), scales = "free_y", nrow = 1) +
  theme_bw() + ggtitle("A")   +
  xlab("Proportion") +
  xlim(c(0, 1))
a

b = ggplot(data = density %>% filter(model == "B"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio), col = "black",  lwd= 1 ) +
  facet_wrap(vars(race_ethnicity), scales = "free_y", nrow = 1) +
  theme_bw() + ggtitle("B")  +   xlab("Proportion") +
  xlim(c(0, 1))
b

c = ggplot(data = density %>% filter(model == "C"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio), col = "black", lwd= 1 ) +
  facet_wrap(vars(race_ethnicity), scales = "free_y", nrow = 1) +
  theme_bw() + ggtitle("C")  +   xlab("Proportion")+
  xlim(c(0, 1))

d = ggplot(data = density %>% filter(model == "D"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio),col = "black", lwd= 1 ) +
  facet_wrap(vars(race_ethnicity), scales = "free_y", nrow = 1) +
  theme_bw() + ggtitle("D")   +   xlab("Proportion")+
  xlim(c(0,  1))

e = ggplot(data = density %>% filter(model == "E"), aes(x = model_ratio)) +
  geom_density(lwd =1, fill = "lightsteelblue1") +
  geom_vline(aes(xintercept = prop), col = "red", lwd= 1 , lty = "dashed") +
  geom_vline(aes(xintercept = obs_ratio),col = "black", lwd= 1 ) +
  facet_wrap(vars(race_ethnicity), scales = "free_y", nrow = 1) +
  theme_bw() + ggtitle("E")  +   xlab("Proportion")+
  xlim(c(0, 1))


plot_grid(a, b, c, d, e, ncol = 1)