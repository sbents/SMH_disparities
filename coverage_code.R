
setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/Year_2024/equity/post")
library(arrow)
library(dplyr)
library(ggplot2)
library(data.table)

### relative ratio 
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


# join NC and CA
df_models_cum = rbind(m1_6, m2_6, m3_6, m4_6, m5_6, 
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

# coverage table 
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98)) # find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile") #cov %>% dplyr::rename(quantile = value) -> cov
cov$quantile = round(cov$quantile, 3) 
cov$quantile = as.numeric(cov$quantile)
setDT(cov)
quantile_bins <- print(cov$quantile)
print(quantile_bins)
qtib = tibble(q1=0.55,q2=0.60,q3=0.65,q4=0.7,q5=0.75, q6=0.8,q7= 0.850,
              q8=0.9,q9=0.95,q10= 0.975, q11=0.99,q13= 0.45,  q14=0.400,q15= 0.350,q16= 0.300,q17=0.250,q18= 0.200,
              q19=0.15,q20=0.1,q21= 0.05,q22= 0.025,q23= 0.011)

mod_ratio <- df_models_cum %>%
  group_by(model, location, race_ethnicity) %>%
  summarise(model_cum_ratio =quantile(model_ratio, qtib),
            quantile = quantile_bins)%>%
  ungroup() 

ggplot(data = mod_ratio %>% filter(location == 37))  + 
  geom_point(aes(x = model_cum_ratio, y = quantile,col = race_ethnicity)) +
  facet_grid(vars(race_ethnicity), vars(model)) + 
  xlab("Ratio") +
  ylab("Quantile") + theme_bw() +
  ggtitle("North Carolina") 

## Real data 
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
head(gs_cum)

ggplot(data = gs_cum) +
  geom_point(aes(x = time_value, y = obs_ratio, col = race_ethnicity)) +
  facet_wrap(vars(location))


# join them 

df_gs_ratio = left_join(mod_ratio, gs_cum, by = c("location", "race_ethnicity")) %>%
  dplyr::select(-obs, -cum) 
setDT(df_gs_ratio)
head(df_gs_ratio)

ggplot(data = df_gs_ratio %>% filter(location == 37))  + 
  geom_point(aes(x = model_cum_ratio, y = quantile,col = race_ethnicity), cex =2) +
  facet_grid(vars(model), vars(race_ethnicity), scales  = "free_x") + 
  xlab("Ratio") +
  ylab("Quantile") + theme_bw() +
  ggtitle("North Carolina") +
  geom_vline(aes(xintercept = obs_ratio), lwd = 1.2)

ggplot(data = df_gs_ratio %>% filter(location == 6))  + 
  geom_point(aes(x = model_cum_ratio, y = quantile,col = race_ethnicity), cex =2) +
  facet_grid(vars(model), vars(race_ethnicity), scales  = "free_x") + 
  xlab("Ratio") +
  ylab("Quantile") + theme_bw() +
  ggtitle("California") +
  geom_vline(aes(xintercept = obs_ratio), lwd = 1.2)

ggplot() +
  geom_point(data = df_gs_ratio %>% filter(model == "B"), aes(x = race_ethnicity, y = obs_ratio), col = "red", cex = 4, shape = 23) +
  geom_point(data = df_gs_ratio %>% filter(model == "B") %>% filter(quantile == .975 | quantile == .025 | quantile == .5)
             , aes(x = race_ethnicity, y = model_cum_ratio, col = as.character(quantile)), cex = 4) +
  facet_wrap(vars(location)) +
  scale_color_manual(values = c("darkslategray2", "darkslategray3","darkslategray4")) +
  theme_bw()  +
  guides(col=guide_legend(title="Quantile"))

head(df_gs_ratio)
cov_ratio <- cov[df_gs_ratio, on = .(quantile = quantile), allow.cartesian=TRUE] %>%
#  .[quantile != 0.5] %>%
  data.table::dcast(location + model +  race_ethnicity  + alpha + obs_ratio ~ variable, value.var = "model_cum_ratio") %>%
  .[, ":=" (cov = ifelse(obs_ratio < upr & obs_ratio > lwr, 1, 0))] %>%
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs_ratio = NULL)]
head(cov_ratio)


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


### WIS on ratio 
head(df_gs_ratio)
# Calculate WIS by location and model
scores_cum <-  df_gs_ratio[, wis(quantile,model_cum_ratio,obs_ratio,IS_components = TRUE),
                 by=.( model, location)] %>%
  group_by(location) %>%
  mutate(avg_model = mean(WIS)) %>%
  mutate(score = WIS/avg_model)

ggplot(data = scores_cum, aes( model, as.character(location ))) +
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
  ggtitle("WIS by location and model compared to average")


### WIS on ratio NC
head(df_gs_ratio)
# Calculate WIS by location and model
scores_cum <-  df_gs_ratio[, wis(quantile,model_cum_ratio,obs_ratio,IS_components = TRUE),
                           by=.( model, location, race_ethnicity)] %>%
  group_by(location, race_ethnicity) %>%
  mutate(avg_model = mean(WIS)) %>%
  mutate(score = WIS/avg_model)

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

exp(-10)


###########################################################
###########################################################
# All models coverage and WIS on incident deaths 

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


# plot coverage by model against alpha 
cov_modes %>%
  group_by(alpha, race_ethnicity, location, model) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  ggplot(aes(x = alpha, y = sum_cov, col = race_ethnicity)) +
  geom_line(linewidth = 1) + 
  theme_bw() +
  ylab("Coverage") +
  xlab("Alpha value")+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  facet_grid(vars(model), vars(location))
# plot coverage by race against alpha
cov_modes %>%
    filter(race_ethnicity != "overall") %>%
    group_by(alpha, race_ethnicity, location, model) %>%
    summarize(sum_cov=sum(cov)/n()) %>%
    ggplot(aes(x = alpha, y = sum_cov, col = model)) +
    geom_line(linewidth = 1) + 
    theme_bw() +
    ylab("Coverage") +
    xlab("Alpha value")+
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    facet_grid(vars(location), vars(race_ethnicity))

#  WIS function
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
  
  
# Calculate WIS by location and model
scores <-  df_gs[, wis(output_type_id,value,obs,IS_components = TRUE),
                  by=.( model, location)] %>%
  group_by(location) %>%
  mutate(avg_model = mean(WIS)) %>%
  mutate(score = WIS/avg_model)

ggplot(data = scores, aes( model, as.character(location ))) +
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
  ggtitle("WIS by location and model compared to average")

# Calculate WIS by race ethnicity and model 
scores1 <-  df_gs[, wis(output_type_id,value,obs,IS_components = TRUE),
                 by=.( model, location, race_ethnicity)] %>%
  group_by(location, race_ethnicity) %>%
  mutate(avg_model = mean(WIS)) %>%
  mutate(score = WIS/avg_model)

ggplot(data = scores_plot1 %>%
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
 # ggtitle("California: WIS by r/e and model compared to average") 
  ggtitle("North Carolina: WIS by r/e and model compared to average")









##########################################################################################


#################################################################################
## Sung-mok's code 

dc <- arrow::open_dataset(paste0(repo_data, folder_path), partitioning = "model_name", 
                          format = "parquet", schema = schema, 
                          factory_options = list(
                            exclude_invalid_files = TRUE))

df_quantile <- dplyr::filter(dc, type == "quantile") %>%  #& (target == "inc death" |  target == "inc hosp")
  dplyr::collect()
setDT(df_quantile)
df_quantile <- df_quantile %>%
  .[, origin_date := as.IDate(origin_date)]

# truth data (from fluview - check with Lucie)
gs_inc_hosp <- setDT(read.csv(paste0(repo_data, "visualization/data-goldstandard/hospitalization.csv"))) %>% 
  .[time_value > "2023-04-16" & time_value <= "2023-12-16"]
gs_inc_death <-  setDT(read.csv(paste0(repo_data,"visualization/data-goldstandard/fv_death_incidence_num.csv"))) %>%
  .[time_value > "2023-04-16" & time_value <= "2023-12-16"]
gold_standard_data <- rbindlist(l = list( copy(gs_inc_hosp) %>% 
                                            .[, target := "inc hosp"], 
                                          copy(gs_inc_death) %>% 
                                            .[, target := "inc death"], 
                                          copy(gs_inc_hosp) %>% 
                                            .[, value := cumsum(value), 
                                              by = .(geo_value_fullname, fips)] %>% 
                                            .[, target := "cum hosp"], 
                                          copy(gs_inc_death) %>% 
                                            .[, value := cumsum(value), 
                                              by = .(geo_value_fullname, fips)] %>% 
                                            .[, target := "cum death"])) %>%
  .[, time_value := as.IDate(time_value)] %>%
  dplyr::rename(obs = "value")


# add observations to df_quantile
df_quantile <- df_quantile %>% 
  .[, time_value :=  origin_date + horizon*7-1] %>%
  .[gold_standard_data, on = .(time_value, target, location = fips)] %>%
  .[!is.na(origin_date)]

#### CALCULATE COVERAGE --------------------------------------------------------
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98))
# find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile")
#cov %>% dplyr::rename(quantile = value) -> cov
cov$quantile = round(cov$quantile, 3)
setDT(cov)
head(cov)

# merge with proj to assign alpha and lwr/upr to each quantile in proj
# cov <- cov[df_quantile, on = .(quantile = type_id)] %>% 
#   .[quantile != 0.5] %>% 
#   # reshape to make lwr and upr columnns
#   data.table::dcast(origin_date + scenario_id + location + target + horizon + time_value  + 
#                       model_name + alpha + obs ~ variable, value.var = "value")  %>%
#   .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>% 
#   .[, ":=" (upr = NULL,
#             lwr = NULL, 
#             obs = NULL)]

cov <- cov[df_quantile, on = .(quantile = type_id), allow.cartesian=TRUE] %>%
  .[quantile != 0.5] %>%
  data.table::dcast(origin_date + scenario_id + location + target + horizon + time_value +
                      model_name + alpha + obs ~ variable, value.var = "value") %>%
  .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>%
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs = NULL)]


#### PLOT RESULTS --------------------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# COVERAGE ACROSS LOCATIONS AND WEEKS FOR EACH SCENARIO/TARGET
options(repr.plot.width=10,repr.plot.height=8)


cov %>%
  .[, max_time_value := max(time_value), by = .(target)] %>%
  .[!is.na(cov) &
      time_value == max_time_value &
      substr(scenario_id,1,1) %in% c("E", "F") &
      substr(target, 1,3) == "cum" &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, model_name, scenario_id)] %>%
  
  mutate(target=case_when(target==c("cum death")~c("Deaths"), target==c("cum hosp")~c("Hospitalizations")),
         scenario_id=case_when(scenario_id==c("E-2023-04-16")~c("Booster for all \n Low immune escape"), 
                               scenario_id==c("F-2023-04-16")~c("Booster for all \n High immune escape"))) %>%
  mutate(model_name=case_when(model_name==c("Ensemble_LOP")~c("Ensemble"), TRUE~model_name)) %>%
  filter(model_name %in% c('Ensemble', 'JHU_IDD-CovidSP', 'MOBS_NEU-GLEAM_COVID', 'NotreDame-FRED', 
                           'UNCC-hierbin', 'USC-SIkJalpha', 'UTA-ImmunoSEIRS','UVA-adaptive','UVA-EpiHiper')) %>%
  ggplot(aes(x = alpha, y = cov_summ, color = model_name)) +
  geom_abline(linetype = "dashed") +
  geom_line(linewidth = 1) + 
  facet_grid(cols = vars(scenario_id), 
             rows = vars(factor(target, levels = c("Hospitalizations", "Deaths")))) + 
  labs(x = "Expected coverage", y = "Actual coverage") +
  scale_x_continuous(expand = c(0,0), label = percent) +
  scale_y_continuous(expand = c(0,0), label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=11, family="sans",color="black"),
        plot.title = element_text(size=17, family="sans",color="black"),
        strip.text = element_text(size = 13, family="sans",color="black"),
        axis.title.y = element_text(size=15, family="sans",color="black"),
        axis.text.y = element_text(size=10, family="sans",color="black"),
        axis.text.x = element_text(angle=45, hjust=1), 
        legend.title = element_blank(), legend.text = element_text(size=11),
        panel.spacing = unit(2, "lines")) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5)))


#############################################################################



















































## summarize to quantiles 
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98))
quantile_bins <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

#quantile_bins <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
#quantile_bin_names <- as.character(c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99))
# bin labels
#qtib = tibble(q1=0.010,q2=0.025,q3=0.050,q4=0.100,q5=0.150,q6=0.200,q7= 0.250,
#              q8=0.300,q9=0.350,q10= 0.400, q11=0.450,q12= 0.500,q13= 0.550,
#              q14=0.600,q15= 0.650,q16= 0.700,q17=0.750,q18= 0.800,
#              q19=0.850,q20=0.900,q21= 0.950,q22= 0.975,q23= 0.990)
# quantiles
quant_mobs <- df_mobs %>%
  ungroup() %>%
  group_by(location, race_ethnicity, time_value) %>%
  summarise(quant=quantile(value, qtib),
            probs = quantile_bin_names)%>%
  as.data.frame()
head(quant_mobs)


# truth data, value = obs and cumsum 
gold_standard_data = read.csv("target_data_phase2.csv") %>%
  mutate(eval = value + min_suppressed) %>%
  filter(date > "2020-11-14" ) %>%
  filter(target == "inc death") %>%
  mutate(time_value = as.Date(date), location = as.character(location)) %>%
  group_by(location, race_ethnicity) %>%
  mutate(obs = value)
# plot
ggplot(data = gold_standard_data, aes(x  = time_value, y = obs, col = race_ethnicity))+
  geom_line()+
  facet_wrap(vars(location))









####################################################################
## Incidence
####################################################################
library(tidyr)
#A 
A_incident_hosp_combined = stack_hvoi %>% 
  as_tibble()%>%
  group_by(target_end_date)%>%
  summarise(quant=quantile(value, qtib),
            probs = quantile_bin_names)%>%
  as.data.frame()%>%
  tidyr::pivot_wider(values_from = quant,names_from=probs)%>%
  mutate(target_end_date = row_number()) %>%
  pivot_longer(c(`0.01`:`0.99`),names_to = "quantile",values_to = "value") %>%
  mutate(location="US",scenario_id = "A-2022-11-13",
         model_projection_date="2022-11-13",
         type="quantile",scenario_name="highVac_optImm",age_group="0-130",
         value = value*4.33)%>%
  dplyr::select("location","quantile","scenario_id","model_projection_date",
                "target_end_date","type","scenario_name","value","age_group")

A_incident_hosp_combined_final = write.csv(A_incident_hosp_combined, "C:/Users/bentssj/OneDrive - National Institutes of Health/flumodel/A_incident_hosp_combined_final.csv" )




# truth data, value = obs and cumsum 
gold_standard_data = read.csv("target_data_phase2.csv") %>%
  mutate(eval = value + min_suppressed) %>%
  filter(date > "2020-11-14" ) %>%
  filter(target == "inc death") %>%
  mutate(time_value = as.Date(date), location = as.character(location)) %>%
  group_by(location, race_ethnicity) %>%
  mutate(obs = cumsum(value))
# plot
ggplot(data = gold_standard_data, aes(x  = time_value, y = obs, col = race_ethnicity))+
  geom_line()+
  facet_wrap(vars(location))


# add observations to gold star data 
df_mobs = df_mobs %>% left_join(gold_standard_data, by = c("time_value", "location", "race_ethnicity", "target"))

# calculate coverage 


#### CALCULATE COVERAGE --------------------------------------------------------
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98))
# find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile")
#cov %>% dplyr::rename(quantile = value) -> cov
cov$quantile = round(cov$quantile, 3)
setDT(cov)
head(cov)


cov <- cov[df_quantile, on = .(quantile = type_id), allow.cartesian=TRUE] %>%
  .[quantile != 0.5] %>%
  data.table::dcast(origin_date + scenario_id + location + target + horizon + time_value +
                      model_name + alpha + obs ~ variable, value.var = "value") %>%
  .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>%
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs = NULL)]





#### PLOT RESULTS --------------------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# COVERAGE ACROSS LOCATIONS AND WEEKS FOR EACH SCENARIO/TARGET
options(repr.plot.width=10,repr.plot.height=8)

cov %>%
  .[, max_time_value := max(time_value), by = .(target)] %>%
  .[!is.na(cov) &
      time_value == max_time_value &
      substr(scenario_id,1,1) %in% c("E", "F") &
      substr(target, 1,3) == "cum" &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, model_name, scenario_id)] %>%
  
  mutate(target=case_when(target==c("cum death")~c("Deaths"), target==c("cum hosp")~c("Hospitalizations")),
         scenario_id=case_when(scenario_id==c("E-2023-04-16")~c("Booster for all \n Low immune escape"), 
                               scenario_id==c("F-2023-04-16")~c("Booster for all \n High immune escape"))) %>%
  mutate(model_name=case_when(model_name==c("Ensemble_LOP")~c("Ensemble"), TRUE~model_name)) %>%
  filter(model_name %in% c('Ensemble', 'JHU_IDD-CovidSP', 'MOBS_NEU-GLEAM_COVID', 'NotreDame-FRED', 
                           'UNCC-hierbin', 'USC-SIkJalpha', 'UTA-ImmunoSEIRS','UVA-adaptive','UVA-EpiHiper')) %>%
  ggplot(aes(x = alpha, y = cov_summ, color = model_name)) +
  geom_abline(linetype = "dashed") +
  geom_line(linewidth = 1) + 
  facet_grid(cols = vars(scenario_id), 
             rows = vars(factor(target, levels = c("Hospitalizations", "Deaths")))) + 
  labs(x = "Expected coverage", y = "Actual coverage") +
  scale_x_continuous(expand = c(0,0), label = percent) +
  scale_y_continuous(expand = c(0,0), label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=11, family="sans",color="black"),
        plot.title = element_text(size=17, family="sans",color="black"),
        strip.text = element_text(size = 13, family="sans",color="black"),
        axis.title.y = element_text(size=15, family="sans",color="black"),
        axis.text.y = element_text(size=10, family="sans",color="black"),
        axis.text.x = element_text(angle=45, hjust=1), 
        legend.title = element_blank(), legend.text = element_text(size=11),
        panel.spacing = unit(2, "lines")) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5)))

















#################################################################################



# add observations to gold star data 
df_mobs <- df_mobs %>% 
  .[, time_value :=  origin_date + horizon*7-1] %>%
  .[gold_standard_data, on = .(time_value, target, location = fips)] %>%
  .[!is.na(origin_date)]


df_quantile <- dplyr::filter(dc, type == "quantile") %>%  #& (target == "inc death" |  target == "inc hosp")
  dplyr::collect()
setDT(df_quantile)
df_quantile <- df_quantile %>%
  .[, origin_date := as.IDate(origin_date)]

# truth data (from fluview - check with Lucie)
gs_inc_hosp <- setDT(read.csv(paste0(repo_data, "visualization/data-goldstandard/hospitalization.csv"))) %>% 
  .[time_value > "2023-04-16" & time_value <= "2023-12-16"]
gs_inc_death <-  setDT(read.csv(paste0(repo_data,"visualization/data-goldstandard/fv_death_incidence_num.csv"))) %>%
  .[time_value > "2023-04-16" & time_value <= "2023-12-16"]
gold_standard_data <- rbindlist(l = list( copy(gs_inc_hosp) %>% 
                                            .[, target := "inc hosp"], 
                                          copy(gs_inc_death) %>% 
                                            .[, target := "inc death"], 
                                          copy(gs_inc_hosp) %>% 
                                            .[, value := cumsum(value), 
                                              by = .(geo_value_fullname, fips)] %>% 
                                            .[, target := "cum hosp"], 
                                          copy(gs_inc_death) %>% 
                                            .[, value := cumsum(value), 
                                              by = .(geo_value_fullname, fips)] %>% 
                                            .[, target := "cum death"])) %>%
  .[, time_value := as.IDate(time_value)] %>%
  dplyr::rename(obs = "value")


# add observations to df_quantile
df_quantile <- df_quantile %>% 
  .[, time_value :=  origin_date + horizon*7-1] %>%
  .[gold_standard_data, on = .(time_value, target, location = fips)] %>%
  .[!is.na(origin_date)]

#### CALCULATE COVERAGE --------------------------------------------------------
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98))
# find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile")
cov %>% dplyr::rename(quantile = value) -> cov
cov$quantile = round(cov$quantile, 3)
setDT(cov)

# merge with proj to assign alpha and lwr/upr to each quantile in proj
# cov <- cov[df_quantile, on = .(quantile = type_id)] %>% 
#   .[quantile != 0.5] %>% 
#   # reshape to make lwr and upr columnns
#   data.table::dcast(origin_date + scenario_id + location + target + horizon + time_value  + 
#                       model_name + alpha + obs ~ variable, value.var = "value")  %>%
#   .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>% 
#   .[, ":=" (upr = NULL,
#             lwr = NULL, 
#             obs = NULL)]

cov <- cov[df_quantile, on = .(quantile = type_id), allow.cartesian=TRUE] %>%
  .[quantile != 0.5] %>%
  data.table::dcast(origin_date + scenario_id + location + target + horizon + time_value +
                      model_name + alpha + obs ~ variable, value.var = "value") %>%
  .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>%
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs = NULL)]


#### PLOT RESULTS --------------------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# COVERAGE ACROSS LOCATIONS AND WEEKS FOR EACH SCENARIO/TARGET
options(repr.plot.width=10,repr.plot.height=8)

cov %>%
  .[, max_time_value := max(time_value), by = .(target)] %>%
  .[!is.na(cov) &
      time_value == max_time_value &
      substr(scenario_id,1,1) %in% c("E", "F") &
      substr(target, 1,3) == "cum" &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, model_name, scenario_id)] %>%
  
  mutate(target=case_when(target==c("cum death")~c("Deaths"), target==c("cum hosp")~c("Hospitalizations")),
         scenario_id=case_when(scenario_id==c("E-2023-04-16")~c("Booster for all \n Low immune escape"), 
                               scenario_id==c("F-2023-04-16")~c("Booster for all \n High immune escape"))) %>%
  mutate(model_name=case_when(model_name==c("Ensemble_LOP")~c("Ensemble"), TRUE~model_name)) %>%
  filter(model_name %in% c('Ensemble', 'JHU_IDD-CovidSP', 'MOBS_NEU-GLEAM_COVID', 'NotreDame-FRED', 
                           'UNCC-hierbin', 'USC-SIkJalpha', 'UTA-ImmunoSEIRS','UVA-adaptive','UVA-EpiHiper')) %>%
  ggplot(aes(x = alpha, y = cov_summ, color = model_name)) +
  geom_abline(linetype = "dashed") +
  geom_line(linewidth = 1) + 
  facet_grid(cols = vars(scenario_id), 
             rows = vars(factor(target, levels = c("Hospitalizations", "Deaths")))) + 
  labs(x = "Expected coverage", y = "Actual coverage") +
  scale_x_continuous(expand = c(0,0), label = percent) +
  scale_y_continuous(expand = c(0,0), label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=11, family="sans",color="black"),
        plot.title = element_text(size=17, family="sans",color="black"),
        strip.text = element_text(size = 13, family="sans",color="black"),
        axis.title.y = element_text(size=15, family="sans",color="black"),
        axis.text.y = element_text(size=10, family="sans",color="black"),
        axis.text.x = element_text(angle=45, hjust=1), 
        legend.title = element_blank(), legend.text = element_text(size=11),
        panel.spacing = unit(2, "lines")) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5)))



##########################################################################################
# adapt for disparities 

https://github.com/midas-network/covid19-smh-research/tree/main/model-output/MOBS_NEU-COVACS_SEIR

2024-06-25-MOBS_NEU-COVACS_SEIR.parquet
# open a model's data 
dc <- arrow::open_dataset(paste0(2024-06-25-MOBS_NEU-COVACS_SEIR, midas-network/covid19-smh-research/tree/main/model-output), partitioning = "MOBS_NEU-COVACS_SEIR", 
                          format = "parquet", schema = schema, 
                          factory_options = list(
                            exclude_invalid_files = TRUE))

df_quantile <- dplyr::filter(dc, type == "quantile") %>%  #& (target == "inc death" |  target == "inc hosp")
  dplyr::collect()
setDT(df_quantile)
df_quantile <- df_quantile %>%
  .[, origin_date := as.IDate(origin_date)]


gs_inc_death <-  setDT(read.csv(paste0(repo_data,"visualization/data-goldstandard/fv_death_incidence_num.csv"))) %>%
  .[time_value > "2023-04-16" & time_value <= "2023-12-16"]
gold_standard_data <- rbindlist(l = list( copy(gs_inc_hosp) %>% 
                                            .[, target := "inc hosp"], 
                                          copy(gs_inc_death) %>% 
                                            .[, target := "inc death"], 
                                          copy(gs_inc_hosp) %>% 
                                            .[, value := cumsum(value), 
                                              by = .(geo_value_fullname, fips)] %>% 
                                            .[, target := "cum hosp"], 
                                          copy(gs_inc_death) %>% 
                                            .[, value := cumsum(value), 
                                              by = .(geo_value_fullname, fips)] %>% 
                                            .[, target := "cum death"])) %>%
  .[, time_value := as.IDate(time_value)] %>%
  dplyr::rename(obs = "value")


### code that works 


# 2024-06-25-NIH_UIUC-RIFTcov.parquet
# 2024-06-25-MOBS_NEU-COVACS_SEIR.parquet
# model data 
df_mod <- arrow::read_parquet("2024-06-25-NIH_UIUC-RIFTcov.parquet") %>%
  filter(target == "inc death") %>%
  mutate(time_value = as.Date(origin_date) + horizon*7-1 ) %>%
  mutate(location = as.numeric(location), time_value = as.Date(time_value))

# plot
ggplot(data = df_mod, aes(x  = time_value, y = value, col = race_ethnicity))+
  geom_line()+
  facet_wrap(vars(location)) 

# coverage table 
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98)) # find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile") #cov %>% dplyr::rename(quantile = value) -> cov
cov$quantile = round(cov$quantile, 3) 
setDT(cov)
head(cov)


## summarize to quantiles 
#cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98))
#quantile_bins <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

#quantile_bins <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
#quantile_bin_names <- as.character(c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99))
# bin labels
#qtib = tibble(q1=0.010,q2=0.025,q3=0.050,q4=0.100,q5=0.150,q6=0.200,q7= 0.250,
#             q8=0.300,q9=0.350,q10= 0.400, q11=0.450,q12= 0.500,q13= 0.550,              q14=0.600,q15= 0.650,q16= 0.700,q17=0.750,q18= 0.800,
#             q19=0.850,q20=0.900,q21= 0.950,q22= 0.975,q23= 0.990)


quantile_bins <- print(cov$quantile)
print(quantile_bins)
qtib = tibble(q1=0.55,q2=0.60,q3=0.65,q4=0.7,q5=0.75, q6=0.8,q7= 0.850,
              q8=0.9,q9=0.95,q10= 0.975, q11=0.99,q13= 0.45,  q14=0.400,q15= 0.350,q16= 0.300,q17=0.250,q18= 0.200,
              q19=0.15,q20=0.1,q21= 0.05,q22= 0.025,q23= 0.011)


#qtib = tibble(q1=0.010,q2=0.025,q3=0.050,q4=0.100,q5=0.150,q6=0.200,q7= 0.250,
#              q8=0.300,q9=0.350,q10= 0.400, q11=0.450,q13= 0.550,  q14=0.600,q15= 0.650,q16= 0.700,q17=0.750,q18= 0.800,
#              q19=0.850,q20=0.900,q21= 0.950,q22= 0.975,q23= 0.990)
#quantile_bin_names <- c(.01, .025, seq(.05, .45, .05), seq(.55, .95, .05), .975, .99)


#quantile_name = print(cov$alpha)
# quantiles
quant_mod <- df_mod %>%
  group_by(location, race_ethnicity, time_value) %>%
  summarise(value =quantile(value, qtib),
            quantile = quantile_bins)%>%
  ungroup() 

ggplot(data = quant_mod %>%
         filter(location == 6) %>%
         filter(race_ethnicity == "asian"), aes(x = time_value, y = value, col = as.character(quantile))) +
  geom_line() + 
  geom_line(data = gold_standard_data %>% filter(location == 6) %>%
              filter(race_ethnicity == "asian"), aes(x = time_value, y = obs), col = "black", lwd = 3)
head(quant_mod)


# truth data, value = obs 
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

gold_standard_data = rbind(gold_standard_data, overall) %>%
  mutate(obs = as.numeric(obs))

# plot
ggplot(data = gold_standard_data, aes(x  = time_value, y = cumsum(obs), col = race_ethnicity))+
  geom_line()+
  facet_wrap(vars(location))


# add observations to gold star data 
df_gs = left_join(quant_mod, gold_standard_data, by = c("time_value", "location", "race_ethnicity")) %>%
  mutate(quantile = as.numeric(quantile))
setDT(df_gs)
head(df_gs)

#ggplot(data = df_gs %>% filter(race_ethnicity == "overall") %>% filter(quantile == .55)) +
#  geom_line(aes(x  = time_value, y = obs, col = race_ethnicity), col = "black") +
#  geom_line(aes(x = time_value, y = value, col = race_ethnicity))

#### CALCULATE COVERAGE --------------------------------------------------------
# coverage table 
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98)) # find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile") #cov %>% dplyr::rename(quantile = value) -> cov
cov$quantile = round(cov$quantile, 3) 
cov$quantile = as.numeric(cov$quantile)
setDT(cov)



cov2 <- cov[df_gs, on = .(quantile = quantile), allow.cartesian=TRUE] %>%
  #  .[quantile != 0.5] %>%
  data.table::dcast(location + race_ethnicity  + time_value + alpha + obs ~ variable, value.var = "value") %>%
  .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>%
  #  filter(race_ethnicity == "asian") %>%
  # filter(location == 6) %>%
  #  filter(alpha == .95)
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs = NULL)]

cov2 %>%
  #  filter(race_ethnicity == "asian") %>%
  #  filter(location == 6) %>%
  # filter(alpha == .95) %>%
  # .[, max_time_value := max(time_value)] %>% # set max time value 
  #  .[!is.na(cov) &
  #    time_value == max_time_value &
  #   # substr(scenario_id,1,1) %in% c("E", "F") &
  #  substr(target, 1,3) == "cum" &
  # !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  group_by(alpha, race_ethnicity, location) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  #  dplyr::summarize(across(cov, ~sum(.x, na.rm = TRUE)))
  # .[, .(cov_summ = sum(cov)/.N), 
  #   by = .( alpha)] %>%
  
  # mutate(target=case_when(target==c("cum death")~c("Deaths"), target==c("cum hosp")~c("Hospitalizations")),
  #         scenario_id=case_when(scenario_id==c("E-2023-04-16")~c("Booster for all \n Low immune escape"), 
  #                               scenario_id==c("F-2023-04-16")~c("Booster for all \n High immune escape"))) %>%
  #  mutate(model_name=case_when(model_name==c("Ensemble_LOP")~c("Ensemble"), TRUE~model_name)) %>%
  #  filter(model_name %in% c('Ensemble', 'JHU_IDD-CovidSP', 'MOBS_NEU-GLEAM_COVID', 'NotreDame-FRED', 
  #                           'UNCC-hierbin', 'USC-SIkJalpha', 'UTA-ImmunoSEIRS','UVA-adaptive','UVA-EpiHiper')) %>%
  ggplot(aes(x = alpha, y = sum_cov, col = race_ethnicity)) +
  #  geom_abline(linetype = "dashed") +
  geom_line(linewidth = 1) + 
  theme_bw() +
  ylab("Coverage") +
  xlab("Alpha value")+
  facet_wrap(vars(location)) 
# ylab("95% coverage")

## cumulative ratio
#################################################################################
#################################################################################
head(gold_standard_data)
# truth data, value = obs 
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
  mutate(obs_ratio = ifelse(location == 6, cum/41937 , cum/7850))
head(gs_cum)

ggplot(data = gs_cum) +
  geom_point(aes(x = time_value, y = obs_ratio, col = race_ethnicity)) +
  facet_wrap(vars(location))


### observed

df_ratio <- arrow::read_parquet("2024-06-25-MOBS_NEU-COVACS_SEIR.parquet") %>%
  filter(target == "inc death") %>%
  mutate(time_value = as.Date(origin_date) + horizon*7-1 ) %>%
  mutate(location = as.numeric(location), time_value = as.Date(time_value)) %>%
  group_by(location, race_ethnicity, run_grouping) %>%
  mutate(cum = cumsum(value)) %>%
  filter(time_value == max(time_value)) %>%
  ungroup() %>%
  group_by(location, run_grouping) %>%
  mutate(overall_value = max(cum)) %>%
  mutate(sample_ratio = cum/overall_value)


quantile_bins <- print(cov$quantile)
print(quantile_bins)
qtib = tibble(q1=0.55,q2=0.60,q3=0.65,q4=0.7,q5=0.75, q6=0.8,q7= 0.850,
              q8=0.9,q9=0.95,q10= 0.975, q11=0.99,q13= 0.45,  q14=0.400,q15= 0.350,q16= 0.300,q17=0.250,q18= 0.200,
              q19=0.15,q20=0.1,q21= 0.05,q22= 0.025,q23= 0.011)

quant_ratio <- df_ratio %>%
  filter(race_ethnicity != "overall") %>%
  group_by(location, race_ethnicity) %>%
  summarise(value =quantile(sample_ratio, qtib),
            quantile = quantile_bins) %>%
  ungroup() 
table(quant_ratio$value)

ggplot(data = quant_ratio, aes(x =race_ethnicity, y = value, col = as.character(quantile))) +
  geom_point() 

# join them 

df_gs_ratio = left_join(quant_ratio, gs_cum, by = c("location", "race_ethnicity")) %>%
  mutate(quantile = as.numeric(quantile)) %>%
  dplyr::select(-obs, -cum)
setDT(df_gs_ratio)
head(df_gs_ratio)

cov3 <- cov[df_gs_ratio, on = .(quantile = quantile), allow.cartesian=TRUE] %>%
  #  .[quantile != 0.5] %>%
  data.table::dcast(location + race_ethnicity  + alpha + obs_ratio ~ variable, value.var = "value") %>%
  .[, ":=" (cov = ifelse(obs_ratio < upr & obs_ratio > lwr, 1, 0))] %>%
  #  filter(race_ethnicity == "asian") %>%
  # filter(location == 6) %>%
  #  filter(alpha == .95)
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs_ratio = NULL)]

head(cov3)


cov3 %>%
  #  substr(target, 1,3) == "cum" &
  # !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  group_by(alpha, race_ethnicity, location) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  #  dplyr::summarize(across(cov, ~sum(.x, na.rm = TRUE)))
  # .[, .(cov_summ = sum(cov)/.N), 
  #   by = .( alpha)] %>%
  
  # mutate(target=case_when(target==c("cum death")~c("Deaths"), target==c("cum hosp")~c("Hospitalizations")),
  #         scenario_id=case_when(scenario_id==c("E-2023-04-16")~c("Booster for all \n Low immune escape"), 
  #                               scenario_id==c("F-2023-04-16")~c("Booster for all \n High immune escape"))) %>%
  #  mutate(model_name=case_when(model_name==c("Ensemble_LOP")~c("Ensemble"), TRUE~model_name)) %>%
  #  filter(model_name %in% c('Ensemble', 'JHU_IDD-CovidSP', 'MOBS_NEU-GLEAM_COVID', 'NotreDame-FRED', 
  #                           'UNCC-hierbin', 'USC-SIkJalpha', 'UTA-ImmunoSEIRS','UVA-adaptive','UVA-EpiHiper')) %>%
  ggplot(aes(x = alpha, y = sum_cov, col = race_ethnicity)) +
  #  geom_abline(linetype = "dashed") +
  geom_line(linewidth = 1) + 
  theme_bw() +
  ylab("Coverage") +
  xlab("Alpha value")+
  facet_wrap(vars(location))

##### load ensemble 


## cumulative ratio
#################################################################################
#################################################################################
head(gold_standard_data)
# truth data, value = obs 
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
  mutate(obs_ratio = ifelse(location == 6, cum/41937 , cum/7850))
head(gs_cum)

ggplot(data = gs_cum) +
  geom_point(aes(x = time_value, y = obs_ratio, col = race_ethnicity)) +
  facet_wrap(vars(location))


### observed

df_ratio <- arrow::read_parquet("2024-06-25-MOBS_NEU-COVACS_SEIR.parquet") %>%
  filter(target == "inc death") %>%
  mutate(time_value = as.Date(origin_date) + horizon*7-1 ) %>%
  mutate(location = as.numeric(location), time_value = as.Date(time_value)) %>%
  group_by(location, race_ethnicity, run_grouping) %>%
  mutate(cum = cumsum(value)) %>%
  filter(time_value == max(time_value)) %>%
  ungroup() %>%
  group_by(location, run_grouping) %>%
  mutate(overall_value = max(cum)) %>%
  mutate(sample_ratio = cum/overall_value)


quantile_bins <- print(cov$quantile)
print(quantile_bins)
qtib = tibble(q1=0.55,q2=0.60,q3=0.65,q4=0.7,q5=0.75, q6=0.8,q7= 0.850,
              q8=0.9,q9=0.95,q10= 0.975, q11=0.99,q13= 0.45,  q14=0.400,q15= 0.350,q16= 0.300,q17=0.250,q18= 0.200,
              q19=0.15,q20=0.1,q21= 0.05,q22= 0.025,q23= 0.011)

quant_ratio <- df_ratio %>%
  filter(race_ethnicity != "overall") %>%
  group_by(location, race_ethnicity) %>%
  summarise(value =quantile(sample_ratio, qtib),
            quantile = quantile_bins) %>%
  ungroup() 
table(quant_ratio$value)

ggplot(data = quant_ratio, aes(x =race_ethnicity, y = value, col = as.character(quantile))) +
  geom_point() 

# join them 

df_gs_ratio = left_join(quant_ratio, gs_cum, by = c("location", "race_ethnicity")) %>%
  mutate(quantile = as.numeric(quantile)) %>%
  dplyr::select(-obs, -cum)
setDT(df_gs_ratio)
head(df_gs_ratio)

cov3 <- cov[df_gs_ratio, on = .(quantile = quantile), allow.cartesian=TRUE] %>%
  #  .[quantile != 0.5] %>%
  data.table::dcast(location + race_ethnicity  + alpha + obs_ratio ~ variable, value.var = "value") %>%
  .[, ":=" (cov = ifelse(obs_ratio < upr & obs_ratio > lwr, 1, 0))] %>%
  #  filter(race_ethnicity == "asian") %>%
  # filter(location == 6) %>%
  #  filter(alpha == .95)
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs_ratio = NULL)]

head(cov3)


cov3 %>%
  #  substr(target, 1,3) == "cum" &
  # !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  group_by(alpha, race_ethnicity, location) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  #  dplyr::summarize(across(cov, ~sum(.x, na.rm = TRUE)))
  # .[, .(cov_summ = sum(cov)/.N), 
  #   by = .( alpha)] %>%
  
  # mutate(target=case_when(target==c("cum death")~c("Deaths"), target==c("cum hosp")~c("Hospitalizations")),
  #         scenario_id=case_when(scenario_id==c("E-2023-04-16")~c("Booster for all \n Low immune escape"), 
  #                               scenario_id==c("F-2023-04-16")~c("Booster for all \n High immune escape"))) %>%
  #  mutate(model_name=case_when(model_name==c("Ensemble_LOP")~c("Ensemble"), TRUE~model_name)) %>%
  #  filter(model_name %in% c('Ensemble', 'JHU_IDD-CovidSP', 'MOBS_NEU-GLEAM_COVID', 'NotreDame-FRED', 
  #                           'UNCC-hierbin', 'USC-SIkJalpha', 'UTA-ImmunoSEIRS','UVA-adaptive','UVA-EpiHiper')) %>%
  ggplot(aes(x = alpha, y = sum_cov, col = race_ethnicity)) +
  #  geom_abline(linetype = "dashed") +
  geom_line(linewidth = 1) + 
  theme_bw() +
  ylab("Coverage") +
  xlab("Alpha value")+
  facet_wrap(vars(location))



##### load ensemble 
# CA
df_ensemble6 <- arrow::read_parquet("ELOP_UT6.gz.parquet") %>% mutate(location = 6)
# NC
df_ensemble37 <- arrow::read_parquet("ELOP_UT37.gz.parquet") %>% mutate(location = 37)
# join
df_ensemble = rbind(df_ensemble6, df_ensemble37) %>%
  mutate(time_value = as.Date(origin_date) + horizon*7-1 ) %>%
  mutate(location = as.numeric(location), time_value = as.Date(time_value))


# truth data, value = obs 
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

gold_standard_data = rbind(gold_standard_data, overall) %>%
  mutate(obs = as.numeric(obs))


# add observations to gold star data 
head(df_ensemble)
df_gs = left_join(df_ensemble, gold_standard_data, by = c("time_value", "location", "race_ethnicity")) 
setDT(df_gs)
head(df_gs)

# coverage table 
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98)) # find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile") #cov %>% dplyr::rename(quantile = value) -> cov
cov$quantile = round(cov$quantile, 3) 
cov$quantile = as.numeric(cov$quantile)
setDT(cov)

head(cov_e)
cov_e <- cov[df_gs, on = .(quantile = output_type_id), allow.cartesian=TRUE] %>%
  .[quantile != 0.5] %>%
  data.table::dcast(location + race_ethnicity  + time_value + alpha + obs ~ variable, value.var = "value") %>%
  .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>%
  #  filter(race_ethnicity == "asian") %>%
  # filter(location == 6) %>%
  #  filter(alpha == .95)
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs = NULL)]

# by race against alpha 
cov_e %>%
  group_by(alpha, race_ethnicity, location) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  ggplot(aes(x = alpha, y = sum_cov, col = race_ethnicity)) +
  geom_line(linewidth = 1) + 
  theme_bw() +
  ylab("Coverage") +
  xlab("Alpha value")+
  facet_wrap(vars(location)) 



# overall over time  
table(cov_e$alpha)
cov_e %>%
  mutate(alpha = as.character(alpha)) %>%
  # filter(alpha  == .5 | alpha == .95) %>%
  filter(alpha  == .95) %>%
  # filter(race_ethnicity == "overall") %>%
  group_by(alpha, time_value, location, race_ethnicity) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  ggplot(aes(x = time_value, y = sum_cov, col = race_ethnicity)) +
  geom_line(linewidth = 1) + 
  theme_bw() +
  ylab("Coverage") +
  xlab("Time value")+
  facet_wrap(vars(location)) 


#######################################



##### load processed cumulative ensembles 
# CA
dfc_ensemble6 <- arrow::read_parquet("ELOP_UT6_CUM.parquet") %>% mutate(location = 6)
# NC
dfc_ensemble37 <- arrow::read_parquet("ELOP_UT37_CUM.parquet") %>% mutate(location = 37)

# join
dfc_ensemble = rbind(dfc_ensemble6, dfc_ensemble37) %>%
  mutate(time_value = as.Date(origin_date) + horizon*7-1 ) %>%
  mutate(location = as.numeric(location), time_value = as.Date(time_value)) %>%
  filter(time_value == max(time_value)) %>%
  group_by(location, output_type_id) %>%
  mutate(overall_value = max(value)) %>%
  ungroup() %>%
  mutate(sample_ratio = value/overall_value) %>%
  head(dfc_ensemble)
table(dfc_ensemble$race_ethnicity)

ggplot(data = dfc_ensemble %>%
         filter(out_p), aes(x = race_ethnicity, y= sample_ratio, col = output_type_id ))+ 
  geom_point() + 
  facet_wrap(vars(location))




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
  mutate(obs_ratio = ifelse(location == 6, cum/41937 , cum/7850))
head(gs_cum)

ggplot(data = gs_cum) +
  geom_point(aes(x = time_value, y = obs_ratio, col = race_ethnicity)) +
  facet_wrap(vars(location))


# join them 

df_gs_ratio = left_join(dfc_ensemble, gs_cum, by = c("location", "race_ethnicity", "time_value")) %>%
  dplyr::select(-obs, -cum, -value) %>%
  drop_na(race_ethnicity) %>%
  filter(sample_ratio != 1)
setDT(df_gs_ratio)
head(df_gs_ratio)

ggplot(data = df_gs_ratio, aes(x = race_ethnicity, y = sample_ratio, col = output_type_id))+
  geom_boxplot(cex = 1) +
  geom_point(aes(x = race_ethnicity, y = obs_ratio), col = "red", cex = 4)+
  facet_wrap(vars(location)) +
  theme_bw() +
  ylab("Proportion race/ethnicity:overall") +
  xlab("Race_ethnicity") + 
  ggtitle("Ensemble LOP untrimmed")

ggplot() +
  # geom_point(cex = 4) +
  ylab("Proportion race/ethnicity:overall") +
  xlab("Race_ethnicity") + 
  ggtitle("Ensemble LOP untrimmed") +
  #  guides(col=guide_legend(title="Quantile")) +
  #  scale_color_discrete(name = "New Legend Title") +
  geom_point(data = df_gs_ratio, aes(x = race_ethnicity, y = obs_ratio), fill = "red", cex = 4, shape = 23) +
  geom_point(data = df_gs_ratio %>% filter(output_type_id == .975 | output_type_id == .025 | output_type_id == .5)
             , aes(x = race_ethnicity, y = sample_ratio, col = as.character(output_type_id)), cex = 4) +
  facet_wrap(vars(location)) +
  scale_color_manual(values = c("darkslategray2", "darkslategray3","darkslategray4")) +
  theme_bw()  +
  guides(col=guide_legend(title="Quantile"))


cov_ratio <- cov[df_gs_ratio, on = .(quantile = output_type_id), allow.cartesian=TRUE] %>%
  .[quantile != 0.5] %>%
  data.table::dcast(location + race_ethnicity  + alpha + obs_ratio ~ variable, value.var = "sample_ratio") %>%
  .[, ":=" (cov = ifelse(obs_ratio < upr & obs_ratio > lwr, 1, 0))] %>%
  .[, ":=" (upr = NULL,
            lwr = NULL,
            obs_ratio = NULL)]
head(cov_ratio)


cov_ratio %>%
  #  substr(target, 1,3) == "cum" &
  # !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  group_by(alpha, location, race_ethnicity) %>%
  summarize(sum_cov=sum(cov)/n()) %>%
  ggplot(aes(x = alpha, y = sum_cov, col = race_ethnicity)) +
  #  geom_abline(linetype = "dashed") +
  #  geom_line(col = "mediumorchid2", linewidth = 3) + 
  geom_line() +
  theme_bw() +
  ylab("Coverage") +
  xlab("Alpha value")+
  facet_wrap(vars(location)) +
  ggtitle("Ratio of cumulative race:overall")



