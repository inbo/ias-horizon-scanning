#---------------------------------------------------------------------
#-----------Add reason why species were excluded from shortlist-------
#---------------------------------------------------------------------

#This script is used to create the google sheets document "total_gridcells_with_info"

#-------load data-------
initial_CM_list<-read_sheet("1TevRYtBiBUhiPI-lRbJPuhkNAAqBhz3g3kEyKZY42YU", sheet="Data")

climate_matching_EU<-read_sheet("1Gf5-sOxQTznOUAutQFNXCaqyxZnY7PXr56v6mdXpzwE",sheet="Data")%>%
  left_join(.,initial_CM_list, by="speciesKey")%>%
  filter(!ID %in% Irena_plants_exclude)

#Read in species that were excluded during climate matching because they are marine
marine_16<-read_sheet("14DoltaPfvDA_ZMf_x7uDRLZhMLYFG8J-KtfGEdwBi4c", sheet="Data")
marine_7<-read_sheet("1Mkyt4mzp1lr-AEG8A86ERck_pXZihqHeO0tbtiPSyx0", sheet="Data")
marine_cm<-rbind(marine_16, marine_7)%>%
  left_join(.,initial_CM_list, by="speciesKey")

#Read in species that were excluded from climate matching
marine<-read_sheet("1limsIJcHBtrJdKMMKwcK8_usdCOs0_v8u71MW1bC5w4",sheet="Data")

marine<-rbind(marine, marine_cm)


#------------------Create vectors with IDs-----------------

#all keys before climate matching (including marine)
climate_matching_IDs<-read_sheet("1TevRYtBiBUhiPI-lRbJPuhkNAAqBhz3g3kEyKZY42YU", sheet="Data")%>%
  filter(ID!="HS7556")%>%
  pull(ID) #4571

#keys of marine species stored during climate matching
marine_cm_IDs<-marine_cm%>%
  pull(ID)

#all IDs right after climate matching
climate_matching_EU_IDs<-climate_matching_EU%>%
  distinct(ID)%>%
  pull(ID) #3603

#All IDs after filtering out those with 20 or less occurences
climate_matching_final_IDs<-climate_matching_EU%>% #3110 left
  filter(n_totaal>20)%>%
  distinct(ID)%>%
  pull(ID)

#IDs that were removed because they had less than 20 occurrences
CM_removed_20records_IDs<-climate_matching_EU_IDs[!climate_matching_EU_IDs%in%climate_matching_final_IDs]

#All IDs on shortlist
shortlist_IDs<-read_sheet("1Wo-SPMJaa8M3R3PnyOkPbVA5ethnWYLxNVxNyQ4_w-g", sheet="Data")%>%
  pull(ID)#3985

#All IDs that were climatematched but were removed during matching but not on the marine list
climate_matching_removed_IDs<-climate_matching_IDs[!climate_matching_IDs%in%climate_matching_EU_IDs]
climate_matching_removed_IDs<-climate_matching_removed_IDs[!climate_matching_removed_IDs%in%marine_cm_IDs]

#unicellular algae IDs
unicellular_algae_IDs<-marine%>%
  filter(kingdom!="Animalia")%>%
  filter(scientificName %in%c("Alexandrium monilatum (J.F.Howell) Balech, 1995",
                              "Marteilia refringens Grizel, Comps, Bonami, Cousserans, Duthoit & Le Pennec, 1974",
                              "Asteromphalus sarcophagus Wallich, 1860",
                              "Olisthodiscus luteus N.Carter, 1937",
                              "Isochrysis galbana Parke",
                              "Alexandrium leei Balech"))%>%
  pull(ID)

#All IDs that were climate matched but were removed from shortlist (to check afterwards if it adds up)
climate_matching_not_on_shortlist_IDs<-climate_matching_IDs[!climate_matching_IDs%in%shortlist_IDs]#1091
rlang::is_true(length(climate_matching_removed_IDs)+length(unicellular_algae_IDs)+length(CM_removed_20records_IDs)==length(climate_matching_not_on_shortlist_IDs))


#-----------Add info column to final dataset----------------
check<-total_gridcells%>%
  mutate(Info = case_when(total_grid_cells>2~"Excluded: present in >2 grid cells in the EU",
                          ID %in% climate_matching_removed_IDs~ "Excluded: either 0 records that could be climate matched or 0% match with European current and/or future climate", #968
                          ID %in% CM_removed_20records_IDs ~ "Excluded: only 20 or less reliable occurrence records globally",
                          ID %in% unicellular_algae_IDs~ "Excluded: unicellular algae",
                          ID %in% shortlist_IDs~ "Included on shortlist",
                          .default=NA))


#-----------Export final dataset----------------
ss0<-gs4_create("total_gridcells_with_info", sheets="Data")
write_sheet(check, ss=ss0, sheet="Data")

