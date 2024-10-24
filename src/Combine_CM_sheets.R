#Climate matching was performed in 7 blocks. While the first 6 were small enough to be run together and exported together
#the final block was too large and was run separately. The data for groups 1-6 and group 7 were therefore stored separately
#and are combined together into one large dataset in this script

climate_match_EU_16<-read_sheet("1NkYZHCkt-Q5zbSu8IX01xQYziY0G1txGCSHCzFr7cK0", sheet="Data")
climate_match_EU_7<-read_sheet("1TMrMk3OcMOHt8joxoodMxfYGuQlJgwqftH5OIFZtdVs",sheet="Data")

climate_matching_EU<-rbind(climate_match_EU_16, climate_match_EU_7)


climate_match_logs_16<-read_sheet("11_7txbcELAJ3-U5PPFhoZvI_5bSrRpQR09sQN_uZCpY", sheet="Data")
climate_match_logs_7<-read_sheet("1AzhZHybhUN4j3wOvNZTE_uZpXJnXZO08PSuZaBhkEJ0", sheet="Data")%>%
  filter(!is.na(Group))

climate_match_logs<-rbind(climate_match_logs_16, climate_match_logs_7)

marine_16<-read_sheet("14DoltaPfvDA_ZMf_x7uDRLZhMLYFG8J-KtfGEdwBi4c", sheet="Data")
marine_7<-read_sheet("1Mkyt4mzp1lr-AEG8A86ERck_pXZihqHeO0tbtiPSyx0", sheet="Data")

marine_cm<-rbind(marine_16, marine_7)

rm(climate_match_EU_7,climate_match_EU_16,climate_match_logs_16, climate_match_logs_7)
