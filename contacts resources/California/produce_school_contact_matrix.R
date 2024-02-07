library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(readr)
library(viridis)
library(ggpubr)
library(tidyverse)
devtools::install_github("seroanalytics/serosolver", ref= "multiple_obs_types") 
library(serosolver)


setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/South Africa")
serosolver_wd = setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/South Africa")

#serosolver_wd <- "~/Documents/GitHub/serosolver/"
#devtools::load_all(serosolver_wd)

# prep data ___________________________________________________________
# singular obs

setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/Scenario Hub/Seattle Flu Study/serosolver")
serology_samples = read.csv("serology_samples_all.csv")
unique(serology_samples$Assay)
head(serology_samples)
rsv = serology_samples %>% filter(Assay== "Flu-H3-HongKong")%>%
  filter(logtiter > 0)  %>%
  filter(age < 6)%>%
  filter(age > .25)
summary(rsv$logtiter)
#  filter(age > .5)
 # filter(age < 11)
rsv_2samp = rsv %>% group_by(sample_ID) %>% sample_n(1, 1)


sample_ID = unique(rsv_2samp$sample_ID)
individual = seq(1, length(sample_ID), 1)
in_dat = data.frame(sample_ID, individual)
rsv_in = left_join(rsv_2samp, in_dat, by = "sample_ID")
head(rsv_in)

ggplot(data = rsv, aes(x = age, y = logtiter))+
  geom_point()

rsv_sample = rsv_in %>% 
  mutate(samples = year)%>% 
  mutate(samples = as.numeric(samples)) %>%
  mutate(titre = logtiter) %>%
  mutate(virus = 2020, run = 1, group = 1) %>%
  mutate(DOB = 2015) %>%
  mutate(obs_type = 1) %>%
  mutate(obs_type = as.numeric(obs_type)) %>%
  select(individual, samples, virus, obs_type, titre, run, group, DOB)%>%
  dplyr::select(-sample_ID)

titre_dat = rsv_sample
titre_dat = data.frame(rsv_sample[c(2:9)])
print(head(titre_dat))

ggplot(data = titre_dat) + geom_violin(aes(x=samples,y=titre,fill=obs_type,group=samples)) + facet_wrap(~obs_type)+
  theme_bw()


#_________________________________

# run serosolver model 
run_name <- "byam_under5_seattle_nov1"
main_wd <- "C:/Users/bentssj/OneDrive - National Institutes of Health/South Africa"
chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/")

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

buckets = 1  ## Time resolution
prior_version <- 2 ## Which prior on the infection histories to use? Prior version 2 is generally preferred
n_chains <- 5 ## Number of MCMC chains to run

rerun <- TRUE

setwd(main_wd)
print(paste0("In directory: ", main_wd))
print(paste0("Saving to: ", save_wd))

set.seed(1)

## Run multiple chains in parallel
cl <- makeCluster(n_chains)
registerDoParallel(cl)

## MCMC settings, not super important but can tweak number of iterations
mcmc_pars <- c("save_block"=100,"thin"=10,"thin_hist"=50,
               "iterations"=500000,
               "adaptive_period"=50000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)



n_obs_types <- 1 ## 1 observation type
obs_type_dists <- c(1) # # number of samples per individual 
sampled_viruses <- c(2020)*buckets

table(titre_dat$samples)
table(titre_dat$DOB)


antigenic_map <- data.frame(inf_times= seq(2020,  2022, 1),x_coord=1,y_coord=1)
strain_isolation_times <- antigenic_map$inf_times
n_times <- length(strain_isolation_times)



## Set up parameter table
par_tab <- read.csv("par_tab_base_new.csv",stringsAsFactors=FALSE)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names %in% c("alpha","beta"),c("values")] <- c(1/3,1/3) ## Can also try c(1,1), or something informative. Just the parameters of a beta distribution which acts as the prior on the per-time attack rate.

## Just some setup of the parameter table to get parameters vaguely like you showed me
par_tab$fixed <- 1
par_tab[par_tab$names %in% c("obs_sd","mu_short","wane"),"fixed"] <- 0


par_tab[par_tab$names %in% c("MIN_TITRE"),"values"] <- 0
par_tab[par_tab$names %in% c("MAX_TITRE"),"values"] <-16

par_tab[par_tab$names %in% c("mu"),"values"] <- 0

summary(titre_dat$titre)


## Create some random measurement offsets -- these add a systematic shift to each observation
## These are unknown parameters to be estimated, assuming the offsets are drawn from ~ norm(0, 1)
par_tab_rhos <- data.frame(names="rho",values= 8,fixed=0,steps=0.1,
                           lower_bound= 6.86, upper_bound= 15.406, lower_start= 10,
                           upper_start=12,type=3) # lower bound is first quadrant

par_tab_rhos$values <- rnorm(length(sampled_viruses), 1)
measurement_indices <- data.frame(virus = sampled_viruses, 
                                  obs_type = rep(1:n_obs_types, each=length(sampled_viruses)),
                                  rho_index=1:(length(sampled_viruses)*n_obs_types))

par_tab <- bind_rows(par_tab, par_tab_rhos)

## Extend parameter table for each aditional observation type
par_tab$obs_type <- 1
par_tab_tmp <- par_tab
if(n_obs_types > 1){
  for(i in 2:n_obs_types){
    par_tab_tmp2 <- par_tab_tmp
    antigenic_map_tmp2 <- antigenic_map_tmp
    antigenic_map_tmp2$obs_type <- i 
    par_tab_tmp2$obs_type <- i
    par_tab <- bind_rows(par_tab, par_tab_tmp2 %>% filter(!(names %in% c("alpha","beta"))))
    antigenic_map <- bind_rows(antigenic_map, antigenic_map_tmp2)
  }
}

## Randomize model parameters
par_tab <- generate_start_tab(par_tab)
head(par_tab)
#par_tab[par_tab$names %in% c("rho"),"values"] <- .95

prior_func <- function(par_tab){
  f <- function(pars){
    pr <- dlnorm(pars["wane"],log(.32),0.34,log=TRUE) #1.8
    pr
  }
}


#sampled_viruses = 2016
#n_obs_types <- 1 ##  one observation type
#obs_type_dists <- c(1)





# create posterior function 
f <- create_posterior_func(par_tab,titre_dat,antigenic_map=antigenic_map,
                           strain_isolation_times = strain_isolation_times,
                           version=prior_version,solve_likelihood=TRUE,n_alive=NULL,
                           measurement_indices_by_time=measurement_indices,
                           data_type=obs_type_dists)

## Time runs and use dopar to run multiple chains in parallel
t1 <- Sys.time()
filenames <- paste0(chain_wd, "/",run_name, "_",1:n_chains)



if(rerun){
  res <- foreach(x = filenames, .packages = c('data.table','plyr',"dplyr","tidyverse","serosolver")) %dopar% {
    #devtools::load_all(serosolver_wd)
    index <- 1
    lik <- -Inf
    inf_hist_correct <- 1
    while((!is.finite(lik) || inf_hist_correct > 0) & index < 100){
      start_tab <- generate_start_tab(par_tab)
      start_inf <- setup_infection_histories_total(titre_dat,strain_isolation_times,2,3) #2,3 
      
      inf_hist_correct <- sum(check_inf_hist(titre_dat, strain_isolation_times, start_inf))
      y <- f(start_tab$values, start_inf)
      lik <- sum(y[[1]])
      index <- index + 1
    }
    
    write.csv(start_tab, paste0(x, "_start_tab.csv"))
    write.csv(start_inf, paste0(x, "_start_inf_hist.csv"))
    write.csv(antigenic_map, paste0(x, "_antigenic_map.csv"))
    write.csv(titre_dat, paste0(x, "_titre_dat.csv"))
    
    res <- run_MCMC(start_tab, 
                    titre_dat, 
                    antigenic_map, 
                    start_inf_hist=start_inf,
                    filename=x,
                    CREATE_POSTERIOR_FUNC=create_posterior_func, 
                    CREATE_PRIOR_FUNC = NULL, # was null, CAN BE PRIOR FUNC
                    version=prior_version,
                    mcmc_pars=mcmc_pars,
                    measurement_indices= measurement_indices, ## measurement_indices, ## NULL
                    measurement_random_effects = FALSE, ## TRUE, ## FALSE
                    solve_likelihood=TRUE,
                    data_type=obs_type_dists)
  }
}


run_time_fast <- Sys.time() - t1
run_time_fast


## load in chains 
chains <- load_mcmc_chains(chain_wd,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)
#print(summary(chains))
list_chains = chains$theta_list_chains
head(list_chains)


list_chains1 <- lapply(list_chains, function(x) x[,c("mu_short", "wane", "total_infections", "rho")])
#list_chains1 <- lapply(list_chains, function(x) x[,c("mu", "wane", "total_infections", "rho")])
#list_chains1 <- lapply(list_chains, function(x) x[,c( "mu")])

#check diagnositcs 
print(gelman.diag(as.mcmc.list(list_chains1)))
effectiveSize(as.mcmc.list(list_chains1))
print(summary(as.mcmc.list(list_chains1)))





# plot
mu_lower = c(1.13, .09, .02)
mu_short = c(2.59, 1.38, .8)
mu_upper = c(4.43, 3.20, 2.77)
wane_lower = c(.47, .10, .04)
wane_short = c(.84, .58, .52)
wane_upper = c(.99, .99, .98)
age = c("<5 yo", "5-10 yo", "20+ yo")

seattle_age = data.frame(mu_lower, mu_short, mu_upper, 
                         wane_lower, wane_short, wane_upper, age)



boost = ggplot(data = seattle_age) + 
  geom_col(aes(x = factor(age, level = c("<5 yo", "5-10 yo", "20+ yo")), y = mu_short), fill = c("lightblue2", "lightblue3", "lightblue4") )+
  geom_errorbar(
    aes(x=age, 
        ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "darkgray" )+ 
  theme_light()+
  ylab("Titer boost")+
  xlab("Age group") #+
 # theme(axis.title.x=element_blank(),
 #       axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank())

wane = ggplot(data = seattle_age) + 
  geom_point(aes(x = factor(age, level = c("<5 yo", "5-10 yo", "20+ yo")), y = wane_short),
             col = c("salmon2", "salmon3", "salmon4"), fill = c("salmon2", "salmon3", "salmon4"), cex = 6,
             shape = 23)+
  geom_errorbar(
    aes(x=age, 
        ymin = wane_lower, 
        ymax = wane_upper), lwd = .8, width = .2, col = c("salmon2", "salmon3", "salmon4"))+ 
  theme_light()+
  geom_point(aes(x = factor(age, level = c("<5 yo", "5-10 yo", "20+ yo")), y = wane_short), fill = c("salmon2", "salmon3", "salmon4"), cex = 7,
             shape = 23)+
  ylab("Percent boost lost in one year")+
  xlab("Age group")
  


h_plots = plot_grid(boost, wane, nrow= 1)
h_plots

# by subtyoe


mu_lower = c(1.13, .52, .1, .08)
mu_short = c (2.59, 1.86, 1.81, 1.38 )
mu_upper = c(4.43, 3.88, 3.84, 3.12)
wane_lower = c(.47, .18, .07, .12)
wane_short = c(.84, .7,.58, .64 )
wane_upper = c(.99, .99, .99,.99)
sub = c("A/H3", "A/H1", "B/Vic", "B/Yam")


seattle_sub = data.frame(mu_lower, mu_short, mu_upper, 
                         wane_lower, wane_short, wane_upper, sub)

boost = ggplot(data = seattle_sub) + 
  geom_col(aes(x = factor(sub, level = c("A/H3", "A/H1", "B/Vic", "B/Yam")), y = mu_short), fill = c("lightblue1", "lightblue2", "lightblue3", "lightblue4") )+
  geom_errorbar(
    aes(x=sub, 
        ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "darkgray" )+ 
  theme_light()+
  ylab("Titer boost")+
  xlab("Subtype") 

wane = ggplot(data = seattle_sub) + 
  geom_point(aes(x = factor(sub, level = c("A/H3", "A/H1", "B/Vic", "B/Yam")), y = wane_short),
             col = c("salmon1", "salmon2", "salmon3", "salmon4"), fill = c("salmon1", "salmon2", "salmon3", "salmon4"), cex = 6,
             shape = 23)+
  geom_errorbar(
    aes(x=sub, 
        ymin = wane_lower, 
        ymax = wane_upper), lwd = .8, width = .2, col = c("salmon1", "salmon2", "salmon3", "salmon4"))+ 
  theme_light()+
  geom_point(aes(x = factor(sub, level = c("A/H3", "A/H1", "B/Vic", "B/Yam")), y = wane_short), fill = c("salmon1", "salmon2", "salmon3", "salmon4"), cex = 7,
             shape = 23)+
  ylab("Percent boost lost in one year")+
  xlab("Subtype")

h_plots = plot_grid(boost, wane, nrow= 1)
h_plots




# south africa age plot 

mu_lower = c(1.1, .81, .87)
mu_short = c(1.25, .89, .92)
mu_upper = c(1.37, .97, .97)
wane_lower = c(0.0, .11, .07)
wane_short = c(.33, .15, .11)
wane_upper = c(.96, .2, .14)
age = c("<5 yo", "5-10 yo", "11+ yo")

seattle_age = data.frame(mu_lower, mu_short, mu_upper, 
                         wane_lower, wane_short, wane_upper, age)



boost = ggplot(data = seattle_age) + 
  geom_col(aes(x = factor(age, level = c("<5 yo", "5-10 yo", "11+ yo")), y = mu_short), fill = c("lightblue2", "lightblue3", "lightblue4") )+
  geom_errorbar(
    aes(x=age, 
        ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "darkgray" )+ 
  theme_light()+
  ylab("Titer boost")+
  xlab("Age group") #+
# theme(axis.title.x=element_blank(),
#       axis.text.x=element_blank(),
#      axis.ticks.x=element_blank())

wane = ggplot(data = seattle_age) + 
  geom_point(aes(x = factor(age, level = c("<5 yo", "5-10 yo", "11+ yo")), y = wane_short),
             col = c("salmon2", "salmon3", "salmon4"), fill = c("salmon2", "salmon3", "salmon4"), cex = 6,
             shape = 23)+
  geom_errorbar(
    aes(x=age, 
        ymin = wane_lower, 
        ymax = wane_upper), lwd = .8, width = .2, col = c("salmon2", "salmon3", "salmon4"))+ 
  theme_light()+
  geom_point(aes(x = factor(age, level = c("<5 yo", "5-10 yo", "11+ yo")), y = wane_short), fill = c("salmon2", "salmon3", "salmon4"), cex = 7,
             shape = 23)+
  ylab("Percent boost lost in one year")+
  xlab("Age group")



h_plots = plot_grid(boost, wane, nrow= 1)
h_plots



