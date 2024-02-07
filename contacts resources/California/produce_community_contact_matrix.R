
library(dplyr)

# load in census tract data 
census_tract_ca = read.csv("censustract_data.csv")

#### total number of effective contacts
eff_contacts=2.79

#### assortativity coefficient
assort_coeff=0.0 # 1 indicates a perfect assortative mixing, while 0 indicates a proportionate mixing

#### list of race/ethnicity
race_list <- c("Asian", "White", "Black", "Latino", "Other")

census_tract_ca %>% group_by(GEOID) %>%
  mutate(contact_Asian = eff_contacts*assort_coeff + (1-assort_coeff)*percent_final[demographic == "Asian"]*eff_contacts, 
         contact_White = (1-assort_coeff)*percent_final[demographic == "White"]*eff_contacts,
         contact_Black = (1-assort_coeff)*percent_final[demographic == "Black"]*eff_contacts, 
         contact_Latino = (1-assort_coeff)*percent_final[demographic == "Latino"]*eff_contacts,
         contact_Other = (1-assort_coeff)*percent_final[demographic == "Other"]*eff_contacts) %>% ungroup() %>%
  filter(demographic==race_list[1]) -> Asian_contacts

census_tract_ca %>% group_by(GEOID) %>%
  mutate(contact_Asian = (1-assort_coeff)*percent_final[demographic == "Asian"]*eff_contacts, 
         contact_White = eff_contacts*assort_coeff + (1-assort_coeff)*percent_final[demographic == "White"]*eff_contacts,
         contact_Black = (1-assort_coeff)*percent_final[demographic == "Black"]*eff_contacts, 
         contact_Latino = (1-assort_coeff)*percent_final[demographic == "Latino"]*eff_contacts,
         contact_Other = (1-assort_coeff)*percent_final[demographic == "Other"]*eff_contacts) %>% ungroup() %>%
  filter(demographic==race_list[2]) -> White_contacts

census_tract_ca %>% group_by(GEOID) %>%
  mutate(contact_Asian = (1-assort_coeff)*percent_final[demographic == "Asian"]*eff_contacts, 
         contact_White = (1-assort_coeff)*percent_final[demographic == "White"]*eff_contacts,
         contact_Black = eff_contacts*assort_coeff + (1-assort_coeff)*percent_final[demographic == "Black"]*eff_contacts, 
         contact_Latino = (1-assort_coeff)*percent_final[demographic == "Latino"]*eff_contacts,
         contact_Other = (1-assort_coeff)*percent_final[demographic == "Other"]*eff_contacts) %>% ungroup() %>%
  filter(demographic==race_list[3]) -> Black_contacts

census_tract_ca %>% group_by(GEOID) %>%
  mutate(contact_Asian = (1-assort_coeff)*percent_final[demographic == "Asian"]*eff_contacts, 
         contact_White = (1-assort_coeff)*percent_final[demographic == "White"]*eff_contacts,
         contact_Black = (1-assort_coeff)*percent_final[demographic == "Black"]*eff_contacts, 
         contact_Latino = eff_contacts*assort_coeff + (1-assort_coeff)*percent_final[demographic == "Latino"]*eff_contacts,
         contact_Other = (1-assort_coeff)*percent_final[demographic == "Other"]*eff_contacts) %>% ungroup() %>%
  filter(demographic==race_list[4]) -> Latino_contacts

census_tract_ca %>% group_by(GEOID) %>%
  mutate(contact_Asian = (1-assort_coeff)*percent_final[demographic == "Asian"]*eff_contacts, 
         contact_White = (1-assort_coeff)*percent_final[demographic == "White"]*eff_contacts,
         contact_Black = (1-assort_coeff)*percent_final[demographic == "Black"]*eff_contacts, 
         contact_Latino = (1-assort_coeff)*percent_final[demographic == "Latino"]*eff_contacts,
         contact_Other = eff_contacts*assort_coeff + (1-assort_coeff)*percent_final[demographic == "Other"]*eff_contacts) %>% ungroup() %>%
  filter(demographic==race_list[5]) -> Other_contacts

rbind(Asian_contacts, White_contacts, Black_contacts, Latino_contacts, Other_contacts) -> census_tract_contacts 



#### total number of children 5-17 yo by race/ethnicity in the state
children_n = census_tract_ca %>%
  group_by(demographic)%>%
  summarize(across(pop_size, ~ sum(as.numeric(.x), na.rm = TRUE)))


contact_race_list <- list()
for(i in 1:length(race_list)){
  community_dem = census_tract_contacts %>%
    filter(demographic == race_list[i]) %>%
    mutate(Asian_k = pop_size*contact_Asian, 
           White_k = pop_size*contact_White,
           Black_k = pop_size*contact_Black,
           Latino_k = pop_size*contact_Latino,
           Other_k = pop_size*contact_Other) %>%
    summarize(across(Asian_k:Other_k, ~ sum(as.numeric(.x), na.rm = TRUE))) 
  
  contact_race_list[[i]] = community_dem/children_n$pop_size[children_n$demographic == race_list[i]]
}

cbind(t(as.matrix(contact_race_list[[1]])), t(as.matrix(contact_race_list[[2]])), 
      t(as.matrix(contact_race_list[[3]])), t(as.matrix(contact_race_list[[4]])), 
      t(as.matrix(contact_race_list[[5]]))) %>% as.data.frame() -> contact_matrix
colnames(contact_matrix) <- c("Asian", "White", "Black", "Latino", "Other")
rownames(contact_matrix) <- colnames(contact_matrix)

contact_matrix
#write.csv(contact_matrix, "community_contact_matrix_CA_proportionate.csv")

library(ggplot2); library(viridis); library(reshape2)

contact_matrix_melt <- melt(as.matrix(contact_matrix))

options(repr.plot.width=7,repr.plot.height=6)
ggplot(data = contact_matrix_melt, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  labs(x="\n Race of individual", y="Race of contact\n") +
  ggtitle("Mixing pattern by race (CT in CA)") +
  theme(axis.title = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=13, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15),
        plot.title = element_text(size=18, family="sans",color="black")) +
  scale_fill_viridis("Rate", direction = 1)


