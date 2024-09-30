#Create empty data frame to store results
climate_matching_EU<-data.frame()

#Create empty logs dataframe to store logs
climate_match_logs<-data.frame(Group=as.integer(""),
                               n_taxa=as.integer(""),
                               n_taxa_downloaded=as.integer(""),
                               n_taxa_after_filtering=as.integer(""),
                               n_taxa_after_removing_marine=as.integer(""),
                               n_taxa_climate_matched=as.integer(""),
                               n_taxa_climate_matched_EU=as.integer("")
)

#create an empty vector for species that have a majority of marine records
mostly_marine_specieskeys <- c()

system.time({
  #Run loop
  for (i in 1:6) {

    #Add group to logs
    climate_match_logs[i,1]<-i

    if(i ==1){
      gbif_download<-"0033829-240906103802322" #9,831,287
    }

    if(i ==2){
      gbif_download<-"0033861-240906103802322" #9,795,168
    }

    if(i ==3){
      gbif_download<-"0033912-240906103802322" #7,673,778
    }

    if(i ==4){
      gbif_download<-"0033954-240906103802322" #4,286,457
    }

    if(i ==5){
      gbif_download<-"0033995-240906103802322" #8,552,687
    }

    if(i ==6){
      gbif_download<-"0034031-240906103802322" #15 600 100
    }

    if(i ==7){
      gbif_download<-"0034087-240906103802322" #25 785 531
    }

#Store download locally: gbif_download<-"0028120-240906103802322"
download_path <-occ_download_get(gbif_download, overwrite = TRUE, path=tempdir()) %>%
  unzip( exdir = file.path(tempdir(),paste0(gbif_download)))

#Import downloaded records
cm_records<-fread(file.path(tempdir(),paste0(gbif_download,"/occurrence.txt")),
                  select=c("speciesKey","scientificName","year","decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", "kingdom","phylum","class","order", "family", "genus", "taxonKey","identificationVerificationStatus")
)

# Remove files from local computer to create space
files<-list.files(file.path(tempdir(),paste0(gbif_download)), pattern = "\\.txt$", full.names=TRUE)
file.remove(files)

#Check how many different species are present
print(paste0("Occurrences for ",length(unique(cm_records$speciesKey))," different species were downloaded.")) #4312

#Add downloaded records to logs
climate_match_logs[i,3]<-length(unique(cm_records$speciesKey))


#-------------------------------------------------
#----STEP 4: Clean data before climate matching----
#-------------------------------------------------
#If there are more than 10 000 000 records, do the cleaning in steps of 10 000 000
if(nrow(cm_records)>10000000){

  cm_store<-data.frame()
  #Count how many blocks we need if we want them to contain around 10 000 000 keys
  nblocks<-ceiling(nrow(cm_records)/10000000)
  for(x in 1:nblocks){

    # Define the start and end row index for each block
    start_row <- (x - 1) * 10000000 + 1
    end_row <- min(x * 10000000, nrow(cm_records))

    # Subset the data for the current block
    cm_block <- cm_records[start_row:end_row, ]

    cm_block<- cm_block %>%
      filter(!identificationVerificationStatus %in% identificationverificationstatus_to_discard)%>%
      cc_cen(buffer=50) %>% # remove points within a buffer of 50m of country centroids
      cc_cap(buffer=50) %>% # remove capitals centroids (buffer 50m)
      cc_inst(buffer=50) %>% # remove zoo and herbaria
      cc_gbif(buffer=50)%>%
      cc_zero()%>%
      select(c("year",
               "speciesKey",
               "decimalLatitude",
               "decimalLongitude",
               "coordinateUncertaintyInMeters"
      )) %>%
      mutate(year_cat = case_when(year <= 1925 ~ "1901-1925",
                                  year <= 1950 ~ "1926-1950",
                                  year <= 1975 ~ "1951-1975",
                                  year <= 2000 ~ "1976-2000",
                                  year > 2000 ~ "2001-2025",
                                  TRUE ~ NA_character_)) %>%
      group_by(year_cat,speciesKey, decimalLatitude,decimalLongitude) %>%
      mutate(coordinateUncertaintyInMeters = ifelse(all(is.na(coordinateUncertaintyInMeters)),
                                                    NA_real_,
                                                    min(coordinateUncertaintyInMeters, na.rm = TRUE)),
             n_obs= 1)%>%
      summarize(n_obs = first(n_obs),
                coordinateUncertaintyInMeters = first(coordinateUncertaintyInMeters)) %>%
      ungroup()

    if(nrow(cm_store)==0){
      cm_store<-cm_block
    }else{
      cm_store<-rbind(cm_store, cm_block)
    }
  rm(cm_block)
  gc()
  }
  cm_records<-cm_store

}else{
  cm_records <- cm_records %>%
    filter(!identificationVerificationStatus %in% identificationverificationstatus_to_discard)%>%
    cc_cen(buffer=50) %>% # remove points within a buffer of 50m of country centroids
    cc_cap(buffer=50) %>% # remove capitals centroids (buffer 50m)
    cc_inst(buffer=50) %>% # remove zoo and herbaria
    cc_gbif(buffer=50)%>%
    cc_zero()%>%
    select(c("year",
             "speciesKey",
             "decimalLatitude",
             "decimalLongitude",
             "coordinateUncertaintyInMeters"
    )) %>%
    mutate(year_cat = dplyr::case_when(year <= 1925 ~ "1901-1925",
                                       year <= 1950 ~ "1926-1950",
                                       year <= 1975 ~ "1951-1975",
                                       year <= 2000 ~ "1976-2000",
                                       year > 2000 ~ "2001-2025",
                                       TRUE ~ NA_character_)) %>%
    dplyr::group_by(year_cat,speciesKey, decimalLatitude,decimalLongitude) %>%
    dplyr::mutate(coordinateUncertaintyInMeters = ifelse(all(is.na(coordinateUncertaintyInMeters)),
                                                         NA_real_,
                                                         min(coordinateUncertaintyInMeters, na.rm = TRUE)),
                  n_obs= 1)%>%
    dplyr::summarize(n_obs = first(n_obs),
                     coordinateUncertaintyInMeters = first(coordinateUncertaintyInMeters)) %>%
    dplyr::ungroup()
}


#Remove datasets that we don't need anymore
if (exists("cm_store")) {
  rm(cm_store)
  print("Removed cm_store")
}


#Convert to sf dataframe
cm_records <- cm_records %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)

#Add number of species left after filtering to logs
climate_match_logs[i,4]<-length(unique(cm_records$speciesKey))


#--------------------------------------------------------------
#----STEP 5: Remove species that have mostly marine records----
#-------------------------------------------------------------

#Remove occurrences in the ocean and remove species for which the majority of records fall in the ocean
cm_records<-cm_records%>%
  st_join(world) %>%
  group_by(speciesKey)%>%
  mutate(n_obs_in_sea=sum(is.na(featurecla) * n_obs, na.rm = TRUE),
         n_obs_totaal=sum(n_obs),
         percent_obs_in_sea=n_obs_in_sea/n_obs_totaal*100)%>%
  select(c("year_cat","speciesKey","decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", "n_obs", "percent_obs_in_sea","featurecla","geometry"))

#Store marine species in separate vector
marine_records<-cm_records%>%
  filter(percent_obs_in_sea>=50)%>%
  pull(speciesKey)%>%
  unique()
mostly_marine_specieskeys <- c(mostly_marine_specieskeys,marine_records)
remove(marine_records)

#Remove marine species and records from climate matching dataframe
cm_records<-cm_records%>%
  filter(!is.na(featurecla), #remove records that fall in the ocean
         percent_obs_in_sea<50)%>% #remove species that have more than 50% of records in the ocean
  select(-c(percent_obs_in_sea,featurecla))

#Store in logs
climate_match_logs[i,5]<-length(unique(cm_records$speciesKey))

#-------------------------------------------------
#----STEP 6: Calculate % climate matching     ----
#-------------------------------------------------

# Define time periods
timeperiodes <- c("1901-1925",
                  "1926-1950",
                  "1951-1975",
                  "1976-2000",
                  "2001-2025")

#Initiate empty dataframe
cm_current <- data.frame()

for(t in timeperiodes){

  print(paste0("Time period: ", t))

  # Determine subset parameters
  start <- as.numeric(substr(t, 0, 4))

  # Subset spatial data
  data_sf_sub <- cm_records%>%
    dplyr::filter(year_cat == t)

  print(paste0("Number of observations: ", nrow(data_sf_sub)))


  if(nrow(data_sf_sub)>0){

    #Load appropriate climate layer
    if(start <= 2000){
      obs_shape <- observed[[t]]
    }else{
      t <- "2001-2025-A1FI"
      obs_shape <- future[[t]]
    }

    #Add GRIDCODE and information of climatic region to each observation, coordinates that don't fall within a climatic region are removed
    data_sf_sub<- st_join(data_sf_sub, obs_shape, join = st_within, left=FALSE)%>%
      select(-ID)%>%
      mutate(GRIDCODE= as.double(GRIDCODE))%>%
      left_join(KG_Rubel_Kotteks_Legend, by = c("GRIDCODE"))

    #Combine these data in a data frame
    if(nrow(cm_current)==0){
      cm_current<-data_sf_sub
    }else{
      cm_current<-rbind(cm_current, data_sf_sub)
    }

  }else {
    warning(paste0("No data was present in the GBIF dataset for ", t))
    next
  }

}

# Cleanup
remove(obs_shape, data_sf_sub, cm_records)

#Calculate percentage climate matching
cm_current<- cm_current %>%
  st_drop_geometry()%>%
  group_by(speciesKey, Classification) %>%
  mutate(n_climate = sum(n_obs, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(speciesKey) %>%
  mutate(n_totaal = sum(n_obs, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(perc_climate = (n_climate/n_totaal)*100) %>%
  distinct(speciesKey, Classification,
           .keep_all = TRUE)  %>%
  select("speciesKey",
         "Classification",
         "Description",
         "n_climate",
         "n_totaal",
         "perc_climate")

#Add number of species that were climate matched to logs
climate_match_logs[i,6]<-length(unique(cm_current$speciesKey))

#----------------------------------------------------------------------------------------
#----STEP 7: Climate matching in EU member states under current and future conditions----
#----------------------------------------------------------------------------------------
# Determine current and nearest future scenarios
scenarios <- c("2001-2025-A1FI",
               "2026-2050-A1FI")

# Create empty future_climate
future_climate <- data.frame() %>%
  mutate(scenario = "",
         KG_GridCode = as.integer(""))

# Calculate KG codes
for (s in scenarios) {

  print(paste0("Time period: ", s))

  #Extract Köppen Geiger layer for the selected future scenario
  shape <- future[[s]]%>%
    st_set_crs(4326)

  #Make sure column gridcode is indicated in capital letters and of class integer
  if (c("gridcode") %in% colnames(shape)) {
    shape <- shape %>%
      rename(GRIDCODE = "gridcode")%>%
      mutate(GRIDCODE = as.integer(GRIDCODE))
  }else{
    shape<- shape%>%
      mutate(GRIDCODE = as.integer(GRIDCODE))
  }


  #Set attributes of shapefile (data columns) to constant throughout the geometry to avoid warning  during st_intersection
  st_agr(shape) = "constant"

  #Intersect the Köppen Geiger layer with the region of interest (EU member states)
  gridcode_intersect<-st_intersection(shape, region_shape)

  #Create an future_climate dataframe holding the gridcodes of the Köppen Geiger zones present in the region of interest and the associated future scenario
  for (g in gridcode_intersect$GRIDCODE) {
    future_climate <- future_climate %>%
      add_row(scenario = s,
              KG_GridCode = g)
  }
}

future_climate <- future_climate%>%
  left_join(y = KG_Rubel_Kotteks_Legend,
            by = c("KG_GridCode" = "GRIDCODE"))%>%
  filter(!is.na(Classification))

#-----Add % climate match data------
cm <- data.frame()

for (b in unique(future_climate$scenario)) {
  future_scenario <- future_climate %>%
    filter(scenario == b)

  cm_int <- cm_current %>%
    filter(Classification %in% future_scenario$Classification
    ) %>%
    mutate(scenario = b)

  if (nrow(cm) == 0) {
    cm <- cm_int
  } else {
    cm <- rbind(cm, cm_int)
  }
}

remove (cm_int)

#----Set thresholds-----

#Number of occurrences used to calculate climate matching
n_limit <- 0
#Minimum percent of climate matching
cm_limit <- 0

#Filter based on thresholds
cm<- cm %>%
  filter(n_totaal >= n_limit,
         perc_climate >= cm_limit)

#Add number of species that were climate matched to logs
climate_match_logs[i,7]<-length(unique(cm$speciesKey))

#Store data in final dataframe
if(nrow(climate_matching_EU)==0){
  climate_matching_EU<-cm
}else{
  climate_matching_EU<-rbind(climate_matching_EU,cm)
}

remove(cm, cm_current, future_climate, gridcode_intersect, shape)
gc()
}

})
#----------------------------------
#----STEP 8: Visualize outcome ----
#----------------------------------

#Add total sum of species keys to logs
addsum<-c("Sum", colSums(climate_match_logs[, 2:6]))
climate_match_logs<-rbind(climate_match_logs,addsum)

length(unique(climate_matching_EU$speciesKey)) #3907

climate_matching_EU%>%
  filter(scenario=="2001-2025-A1FI")%>%
  ggplot( aes(n_totaal)) +
  geom_histogram(binwidth=1, color="black", fill="green") +
  xlab("Number of occurrences used for climate matching") +
  ylab("Number of species")+
  geom_magnify(from = c(xmin =0, xmax = 25, ymin = 0, ymax = 200),
               to = c(xmin =70000, xmax = 100000, ymin = 25, ymax = 175),
               axes="xy")+
  theme_bw()

#----------------------------------
#----STEP 9: add a filter to the data ----
#----------------------------------

climate_matching_final<-climate_matching_EU%>%
  distinct(speciesKey, Classification, .keep_all = TRUE)%>%
  group_by(speciesKey) %>%
  summarise(total_obs_EU_climate = sum(n_climate, na.rm = TRUE))%>%
  filter(total_obs_EU_climate>20)%>%
  pull(speciesKey)


#TO DO:
#3. Check coordinate cleaner steps
#4. Do we want to keep all occurrences at the same coordinate in a time category?
#5. Download DWCA for identificatonverificationstatus

#6. Visualisations (only important at a later stage)
