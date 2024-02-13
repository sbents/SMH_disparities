library(dplyr)

setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/Year_2024/equity/contacts")

#### total number of effective contacts
eff_contacts=11.41

#### assortativity coefficient
assort_coeff=0.0 # 1 indicates a perfect assortative mixing, while 0 indicates a proportionate mixing

#### list of race/ethnicity
race_list <- c("Asian", "White", "Black", "Other")

#### the number of 5-17 year old children in each race/ethnicity by GEOID
publicschool_data = read.csv("publicschool_data_NC_clean.csv")

publicschool_data %>% group_by(School_District) %>%
  mutate(contact_Asian = case_when(demographic == "Asian" ~ 
                                     eff_contacts*assort_coeff + (1-assort_coeff)*percent_final[demographic == "Asian"]*eff_contacts,
                                   TRUE ~ (1-assort_coeff)*percent_final[demographic == "Asian"]*eff_contacts),
         contact_White = case_when(demographic == "White" ~ 
                                     eff_contacts*assort_coeff + (1-assort_coeff)*percent_final[demographic == "White"]*eff_contacts,
                                   TRUE ~ (1-assort_coeff)*percent_final[demographic == "White"]*eff_contacts),
         contact_Black = case_when(demographic == "Black" ~ 
                                     eff_contacts*assort_coeff + (1-assort_coeff)*percent_final[demographic == "Black"]*eff_contacts,
                                   TRUE ~ (1-assort_coeff)*percent_final[demographic == "Black"]*eff_contacts),
         contact_Other = case_when(demographic == "Other" ~ 
                                     eff_contacts*assort_coeff + (1-assort_coeff)*percent_final[demographic == "Other"]*eff_contacts,
                                   TRUE ~ (1-assort_coeff)*percent_final[demographic == "Other"]*eff_contacts)) %>%
  ungroup() -> school_probability

#### total number of children 5-17 yo by race/ethnicity in the state
children_n = school_probability %>%
  group_by(demographic)%>%
  summarize(across(children_5_17, ~ sum(as.numeric(.x), na.rm = TRUE)))

contact_race_list <- list()
for(i in 1:length(race_list)){
  school_dem = school_probability %>%
    filter(demographic == race_list[i]) %>%
    mutate(Asian_k = children_5_17*contact_Asian, 
           White_k = children_5_17*contact_White,
           Black_k = children_5_17*contact_Black,
           Other_k = children_5_17*contact_Other) %>%
    summarize(across(Asian_k:Other_k, ~ sum(as.numeric(.x), na.rm = TRUE))) 
  
  contact_race_list[[i]] = school_dem/children_n$children_5_17[children_n$demographic == race_list[i]]
}

cbind(t(as.matrix(contact_race_list[[1]])), t(as.matrix(contact_race_list[[2]])), 
      t(as.matrix(contact_race_list[[3]])), t(as.matrix(contact_race_list[[4]]))) %>% as.data.frame() -> contact_matrix
colnames(contact_matrix) <- c("Asian", "White", "Black", "Other")
rownames(contact_matrix) <- colnames(contact_matrix)

contact_matrix
write.csv(contact_matrix, "school_contact_matrix_NC_proportionate.csv")


