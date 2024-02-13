library(dplyr)

#setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/Year_2024/equity/contacts")

#### total number of effective contacts
eff_contacts= 2.79

#### assortativity coefficient
assort_coeff=0.0 # 1 indicates a perfect assortative mixing, while 0 indicates a proportionate mixing

#### list of race/ethnicity
race_list <- c("Asian", "White", "Black", "Other")

#### the number of 5-17 year old children in each race/ethnicity by GEOID
census_tract_nc = read.csv("censustract_data_NC.csv")
head(census_tract_nc)

census_tract_nc %>% group_by(GEOID) %>%
  mutate(contact_Asian = case_when(demographic == "Asian" ~ 
                                     eff_contacts*assort_coeff + (1-assort_coeff)*percent[demographic == "Asian"]*eff_contacts,
                                   TRUE ~ (1-assort_coeff)*percent[demographic == "Asian"]*eff_contacts),
         contact_White = case_when(demographic == "White" ~ 
                                     eff_contacts*assort_coeff + (1-assort_coeff)*percent[demographic == "White"]*eff_contacts,
                                   TRUE ~ (1-assort_coeff)*percent[demographic == "White"]*eff_contacts),
         contact_Black = case_when(demographic == "Black" ~ 
                                     eff_contacts*assort_coeff + (1-assort_coeff)*percent[demographic == "Black"]*eff_contacts,
                                   TRUE ~ (1-assort_coeff)*percent[demographic == "Black"]*eff_contacts),
         contact_Other = case_when(demographic == "Other" ~ 
                                     eff_contacts*assort_coeff + (1-assort_coeff)*percent[demographic == "Other"]*eff_contacts,
                                   TRUE ~ (1-assort_coeff)*percent[demographic == "Other"]*eff_contacts)) %>%
  ungroup() -> census_tract_contacts

#### total number of children 5-17 yo by race/ethnicity in the state
ct_n = census_tract_nc %>%
  group_by(demographic)%>%
  summarize(across(pop_size, ~ sum(as.numeric(.x), na.rm = TRUE)))

contact_race_list <- list()
for(i in 1:length(race_list)){
  community_dem = census_tract_contacts %>%
    filter(demographic == race_list[i]) %>%
    mutate(Asian_k = pop_size*contact_Asian, 
           White_k = pop_size*contact_White,
           Black_k = pop_size*contact_Black,
           Other_k = pop_size*contact_Other) %>%
    summarize(across(Asian_k:Other_k, ~ sum(as.numeric(.x), na.rm = TRUE))) 
  
  contact_race_list[[i]] = community_dem/ct_n$pop_size[ct_n$demographic == race_list[i]]
}

cbind(t(as.matrix(contact_race_list[[1]])), t(as.matrix(contact_race_list[[2]])), 
      t(as.matrix(contact_race_list[[3]])), t(as.matrix(contact_race_list[[4]]))) %>% as.data.frame() -> contact_matrix
colnames(contact_matrix) <- c("Asian", "White", "Black", "Other")
rownames(contact_matrix) <- colnames(contact_matrix)

contact_matrix
write.csv(contact_matrix, "community_contact_matrix_NC_proportionate.csv")

