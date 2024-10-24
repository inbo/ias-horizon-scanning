#----------------------------------------------------------------------------------
#Run climate matching without having to download occurrence data anew
#-----------------------------------------------------------------------------------
#This script has the same functionality as the Climate_matching_exercise script, with the difference that it
#starts from the downloads that were generated in the former script to save time.


#------------------Load packages-----------------
library(dplyr)
library(rgbif)
library(sf)
library(data.table)
library(mapview)
library(here)
library(stringr)
library(CoordinateCleaner)
library(ggplot2)
library(rnaturalearth)
library(ggmagnify)


#-----------------Load data-----------------
# Load climate layers
load("./data/input/climate_matching/future.rda")
load("./data/input/climate_matching/legends.rda")
load("./data/input/climate_matching/observed.rda")

# Import legends
KG_Rubel_Kotteks_Legend <- legends$KG_A1FI
KG_Beck <- legends$KG_Beck

# Load shape of the world
world<-ne_countries(scale=10) %>%
  st_make_valid()

#Load shapefile of EU member states
memberstates<-read.csv2(here("./data/input/GRIIS_checklists.csv"))%>%
  pull(EU_state)
region_shape<-st_read("./data/spatial/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")%>%
  filter(ADMIN%in%c(memberstates,"Czech Republic"))%>%
  st_crop(xmin=-13, ymin=20, xmax=41, ymax= 81)%>%
  select(c(ADMIN, geometry))
st_agr(region_shape) = "constant"


#------------Create empty datasets for loop--------------
#Create empty data frame to store results
climate_matching_EU<-data.frame()

#Create empty logs dataframe to store logs
climate_match_logs<-data.frame(Group=as.integer(""),
                               n_taxa=as.integer(""),
                               n_taxa_downloaded=as.integer(""),
                               n_occurrences_downloaded=as.integer(""),
                               n_taxa_after_filtering=as.integer(""),
                               n_taxa_after_removing_marine=as.integer(""),
                               n_taxa_climate_matched=as.integer(""),
                               n_taxa_climate_matched_EU=as.integer("")
)

#create an empty vector for species that have a majority of marine records
mostly_marine_specieskeys <- c()


#---------------------Run climate matching loop-----------------------
#Create a vector with the GBIF downloadkeys of interest
gbif_downloads<-c("0033829-240906103802322", #n = 9,831,287
                  "0033861-240906103802322", #n = 9,795,168
                  "0033912-240906103802322", #n = 7,673,778
                  "0033954-240906103802322", #n = 4,286,457
                  "0033995-240906103802322", #n = 8,552,687
                  "0034031-240906103802322", #n = 15 600 100
                  "0034087-240906103802322" #n = 25 785 531
)
system.time({
  #Run loop
  for (i in 1:length(gbif_downloads)) {

    #Add group to logs
    climate_match_logs[i,1]<-i

    #Extract downloadkey
    gbif_download<-gbif_downloads[i]

    #Store download locally
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
    print(paste0(nrow(cm_records)," occurrences for ",length(unique(cm_records$speciesKey))," different species were downloaded.")) #4312

    #Add downloaded records to logs
    climate_match_logs[i,3]<-length(unique(cm_records$speciesKey))
    climate_match_logs[i,4]<-nrow(cm_records)

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
          mutate(year_cat = cut(year,
                                breaks = c(1900, 1925, 1950, 1975, 2000, 2025),
                                labels = c("1901-1925", "1926-1950", "1951-1975", "1976-2000", "2001-2025"),
                                right = TRUE)) %>%
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
        mutate(year_cat = cut(year,
                              breaks = c(1900, 1925, 1950, 1975, 2000, 2025),
                              labels = c("1901-1925", "1926-1950", "1951-1975", "1976-2000", "2001-2025"),
                              right = TRUE)) %>%
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
    climate_match_logs[i,5]<-length(unique(cm_records$speciesKey))


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
    climate_match_logs[i,6]<-length(unique(cm_records$speciesKey))

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
    climate_match_logs[i,7]<-length(unique(cm_current$speciesKey))

    #----------------------------------------------------------------------------------------
    #----STEP 7: Climate matching in EU member states under current and future conditions----
    #----------------------------------------------------------------------------------------
    # Determine current and nearest future scenarios
    scenarios <- c("2001-2025-A1FI",
                   "2026-2050-A1FI")

    # Create empty future_climate
    future_climate <- data.frame() %>%
      mutate(scenario = "",
             KG_GridCode = as.integer(""),
             surface_area_climatezone=as.integer(""),
             total_surface_area=as.integer(""))

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

      #Calculate the total area occupied by the region shape in km2
      gridcode_intersect_area<-gridcode_intersect%>%
        st_transform(crs = 3035) %>%
        summarise(area_km2 = sum(as.numeric(st_area(.))) / 1e6) %>%
        pull(area_km2)

      #Create an future_climate dataframe holding the gridcodes of the Köppen Geiger zones present in the region of interest and the associated future scenario
      for (g in unique(gridcode_intersect$GRIDCODE)) {

        #Get the surface area of the selected climate zone
        gridcode_area<-gridcode_intersect%>%
          filter(GRIDCODE==g)%>%
          st_transform(crs = 3035) %>%
          summarise(area_km2 = sum(as.numeric(st_area(.))) / 1e6) %>%
          pull(area_km2)

        future_climate <- future_climate %>%
          add_row(scenario = s,
                  KG_GridCode = g,
                  surface_area_climatezone=gridcode_area,
                  total_surface_area= gridcode_intersect_area)
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
        filter(Classification %in% future_scenario$Classification ) %>%
        mutate(scenario = b)%>%
        left_join(., future_scenario, by=c("Description", "Classification","scenario"))%>%
        select(-c("Group", "Precipitation Type", "Level of Heat", "KG_GridCode"))

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
    climate_match_logs[i,8]<-length(unique(cm$speciesKey))

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

#Remove all files except for climate_matching and marine
rm(list = setdiff(ls(), c("climate_match_logs", "climate_matching_EU", "future", "mostly_marine_specieskeys")))

#Load data that were the result of this workflow
climate_matching_EU<-read_sheet("1Gf5-sOxQTznOUAutQFNXCaqyxZnY7PXr56v6mdXpzwE",sheet="Data")
length(unique(climate_matching_EU$speciesKey)) #3604

climate_match_logs<-read_sheet("18SWIn_TXu64rtFFnSunFioUlmgTk8VC0shqDTHgZxVA", sheet="Data") %>%
  mutate_at(2:8, as.integer)


#----------------------------------
#----STEP 8: Process logs ----
#----------------------------------
#Add total sum of species keys to logs
addsum<-c("Sum", NA, colSums(climate_match_logs[, 3:8]))
climate_match_logs<-rbind(climate_match_logs,addsum)
rm(addsum)


#--------------------------------------------------------------------------
#----STEP 9: Visualize number of occurrences used for climate matching ----
#--------------------------------------------------------------------------
climate_matching_EU%>%
  #filter(n_totaal<500000)%>%
  filter(scenario=="2001-2025-A1FI")%>%
  distinct(speciesKey, n_totaal)%>%
  ggplot( aes(n_totaal)) +
  geom_histogram(binwidth=10, color="black", fill="green") +
  xlab("Number of occurrences used for climate matching") +
  ylab("Number of species")+
  geom_magnify(from = c(xmin =0, xmax = 100, ymin = 0, ymax = 210),
               to = c(xmin =1000000, xmax = 2000000, ymin = 50, ymax = 150),
               axes="xy")+
  theme_bw()


#-----------------------------------------------------------------
#----STEP 10: Process and filter data to obtain final dataset ----
#-----------------------------------------------------------------

#---Load datasets needed to append information---

#List with unknown Thematic groups
thematic_groups<-read_sheet("1DkyHPumpbIeqgOhkDXU8F8hjfoLQSDnScC9Ybc_tq9s", sheet="Data")%>%
  select(c(Thematic_group, ID))%>%
  distinct(ID, .keep_all=TRUE)

#Initial list used to extract taxonkeys for climate matching
initial_CM_list<-read_sheet("1TevRYtBiBUhiPI-lRbJPuhkNAAqBhz3g3kEyKZY42YU", sheet="Data")
length(unique(initial_CM_list$speciesKey)) #4572

#List Irena to exclude plants that may have slipped in via correction sheet
Irena_plants_exclude<-read_sheet("1xnQW5rk50Unt9e2P6I3ULrEdxApZ1yj_nnzbSjgzlE0", sheet = "Data")%>%
  filter(DECISION%in% c("arch","native","excludeHigherTaxa"))%>%
  pull(ID)


#--------------Add information-------------
initial_CM_list<-initial_CM_list%>%
  left_join(thematic_groups, by = "ID", suffix = c("_old", "_new")) %>%
  mutate(Thematic_group = ifelse(Thematic_group_old=="Unknown", Thematic_group_new, Thematic_group_old)) %>%
  select(-Thematic_group_old, -Thematic_group_new)%>%
  mutate(Thematic_group=case_when(Thematic_group=="plants"~"Plants",
                                  Thematic_group=="FW Verts/ Marine"~"FW Verts; Marine",
                                  .default = Thematic_group))# Clean up the extra columns

climate_matching_EU<-climate_matching_EU%>%
  left_join(.,initial_CM_list, by="speciesKey")%>%
  filter(!ID %in% Irena_plants_exclude)


length(unique(climate_matching_EU$speciesKey)) #3603

#---------Prepare final datasets------------
#Remove species that had less than 20 records to climate match with
#taxonkey 1420478 not included first time due to n_totaal 17, second time n_totaal 24: check later
climate_matching_final<-climate_matching_EU%>% #3110 left
  filter(n_totaal>20)
length(unique(climate_matching_final$ID))

#Calculate a weighted surface area and percentage climate match with general EU
climate_matching_weighed<-climate_matching_final%>%
  mutate(weighted_surface=perc_climate*surface_area_climatezone/100)%>%
  group_by(scenario, speciesKey,n_totaal)%>%
  summarize(sum_weighted_surface=sum(weighted_surface),
            percentage_climatematch_EU = sum(perc_climate, na.rm = TRUE))%>%
  ungroup()%>%
  group_by(scenario)%>%
  mutate(max_sum_weighted_surface=max(sum_weighted_surface))%>%
  ungroup()%>%
  group_by(speciesKey, scenario)%>%
  mutate(ranked_weighted_surface=sum_weighted_surface/max_sum_weighted_surface*100)%>%
  ungroup()%>%
  select(-c(sum_weighted_surface, max_sum_weighted_surface))%>%
  left_join(.,initial_CM_list, by="speciesKey")


#----------------------------------------------------
#----STEP 11: Plot to get an overview of the data----
#----------------------------------------------------
# Bin the percentage climate match into ranges of 10% and plot the number of species in each range using ggplot2
plot1<-climate_matching_weighed %>%
  group_by(scenario, speciesKey)%>%
  mutate(EU_climatematch_range = cut(percentage_climatematch_EU,
                                     breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                     labels = c("0<-10", "10<-20","20<-30","30<-40","40<-50","50<-60","60<-70","70<-80","80<-90","90<-100"),
                                     right = TRUE)) %>%
  ungroup()%>%
  group_by(EU_climatematch_range, scenario) %>%
  summarise(n_species = n_distinct(speciesKey))%>%
  ggplot(., aes(x = EU_climatematch_range, y = n_species)) +
  geom_bar(stat = "identity", fill = "#f5f29f", color="black") +
  geom_text(aes(label =n_species), vjust = 2) +
  labs(x=NULL,
       y = "Number of Species") +
  facet_grid(~scenario)+
  theme_bw()
plot1

# Bin the percentage climate match into ranges of 10% per thematic group and plot it
plot2 <- climate_matching_weighed %>%
  group_by(scenario, speciesKey)%>%
  mutate(EU_climatematch_range = cut(percentage_climatematch_EU,
                                     breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                     labels = c("0<-10", "10<-20","20<-30","30<-40","40<-50","50<-60","60<-70","70<-80","80<-90","90<-100"),
                                     right = TRUE)) %>%
  ungroup()%>%
  group_by(EU_climatematch_range, scenario, Thematic_group) %>%
  summarise(n_species = n_distinct(speciesKey))%>%
  ggplot(., aes(x = EU_climatematch_range, y = n_species, fill=Thematic_group)) +
  geom_bar(stat = "identity",  color="black") +
  labs(x=NULL,
       y = "Number of Species") +
  facet_grid(~scenario)+
  scale_fill_brewer(palette = "Set3", direction = -1) +
  theme_bw()

# Bin the percentage climate match into ranges of 10% per occupancy potential range and plot it
plot3 <- climate_matching_weighed %>%
  group_by(scenario, speciesKey)%>%
  mutate(EU_climatematch_range = cut(percentage_climatematch_EU,
                                     breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                     labels = c("0<-10", "10<-20","20<-30","30<-40","40<-50","50<-60","60<-70","70<-80","80<-90","90<-100"),
                                     right = TRUE),
         potential_occ_range = cut(ranked_weighted_surface,
                                   breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                   labels = c("0<-10", "10<-20","20<-30","30<-40","40<-50","50<-60","60<-70","70<-80","80<-90","90<-100"),
                                   right = TRUE)%>%
           as.factor() %>%  # Ensure it's a factor
           fct_rev()) %>%  # Then reverse the factor levels
  ungroup()%>%
  group_by(EU_climatematch_range, scenario, potential_occ_range) %>%
  summarise(n_species = n_distinct(speciesKey))%>%
  ggplot(., aes(x = EU_climatematch_range, y = n_species, fill=potential_occ_range)) +
  geom_bar(stat = "identity",  color="black") +
  labs(
    x = "Percentage climate match with EU climate",
    y = "Number of Species") +
  facet_grid(~scenario)+
  scale_fill_brewer(palette = "Spectral", direction = 1) +
  theme_bw()

library(patchwork)
plot1/plot2/plot3

# Bin the percentage climate match into ranges of 10% and plot the number of species in each range using ggplot2
plot4<-climate_matching_weighed %>%
  group_by(scenario, speciesKey)%>%
  mutate(potential_occ_range = cut(ranked_weighted_surface,
                                   breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                   labels = c("0<-10", "10<-20","20<-30","30<-40","40<-50","50<-60","60<-70","70<-80","80<-90","90<-100"),
                                   right = TRUE)) %>%
  ungroup()%>%
  group_by(potential_occ_range, scenario) %>%
  summarise(n_species = n_distinct(speciesKey))%>%
  ggplot(., aes(x = potential_occ_range, y = n_species)) +
  geom_bar(stat = "identity", fill = "#f5f29f", color="black") +
  geom_text(aes(label =n_species), vjust = -0.5) +
  labs(x=NULL,
       y = "Number of Species") +
  facet_grid(~scenario)+
  theme_bw()
plot4

# Bin the percentage climate match into ranges of 10% per thematic group and plot it
plot5 <- climate_matching_weighed %>%
  group_by(scenario, speciesKey)%>%
  mutate(potential_occ_range = cut(ranked_weighted_surface,
                                   breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                   labels = c("0<-10", "10<-20","20<-30","30<-40","40<-50","50<-60","60<-70","70<-80","80<-90","90<-100"),
                                   right = TRUE)) %>%
  ungroup()%>%
  group_by(potential_occ_range, scenario, Thematic_group) %>%
  summarise(n_species = n_distinct(speciesKey))%>%
  ggplot(., aes(x = potential_occ_range, y = n_species, fill=Thematic_group)) +
  geom_bar(stat = "identity",  color="black") +
  labs(x=NULL,
       y = "Number of Species") +
  facet_grid(~scenario)+
  scale_fill_brewer(palette = "Set3", direction = -1) +
  theme_bw()

# Bin the percentage climate match into ranges of 10% per occupancy potential range and plot it
plot6 <- climate_matching_weighed %>%
  group_by(scenario, speciesKey)%>%
  mutate(EU_climatematch_range = cut(percentage_climatematch_EU,
                                     breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                     labels = c("0<-10", "10<-20","20<-30","30<-40","40<-50","50<-60","60<-70","70<-80","80<-90","90<-100"),
                                     right = TRUE)%>%
           as.factor() %>%  # Ensure it's a factor
           fct_rev(),
         potential_occ_range = cut(ranked_weighted_surface,
                                   breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                   labels = c("0<-10", "10<-20","20<-30","30<-40","40<-50","50<-60","60<-70","70<-80","80<-90","90<-100"),
                                   right = TRUE))%>%  # Then reverse the factor levels
  ungroup()%>%
  group_by(EU_climatematch_range, scenario, potential_occ_range) %>%
  summarise(n_species = n_distinct(speciesKey))%>%
  ggplot(., aes(x = potential_occ_range, y = n_species, fill=EU_climatematch_range)) +
  geom_bar(stat = "identity",  color="black") +
  labs(
    x = "Potential occupancy range (%)",
    y = "Number of Species") +
  facet_grid(~scenario)+
  scale_fill_brewer(palette = "Spectral", direction = 1) +
  theme_bw()

plot4/plot5/plot6

#----------------------------------------------------
#----STEP 12: filter data to get to a shortlist----
#----------------------------------------------------

climate_matching_weighed%>%
  filter(scenario=="2026-2050-A1FI") %>%
  group_by(scenario, speciesKey)%>%
  mutate(potential_occ_range = cut(ranked_weighted_surface,
                                   breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                   labels = c("0<-10", "10<-20","20<-30","30<-40","40<-50","50<-60","60<-70","70<-80","80<-90","90<-100"),
                                   right = TRUE))%>%  # Then reverse the factor levels
  ungroup()%>%
  group_by(Thematic_group, scenario, potential_occ_range) %>%
  summarise(n_species = n_distinct(speciesKey))%>%
  ggplot(., aes(x = potential_occ_range, y = n_species)) +
  geom_bar(stat = "identity",  color="black", fill= "#f5f29f") +
  labs(
    x = "Potential occupancy range (%)",
    y = "Number of Species") +
  facet_wrap(~Thematic_group, scales="free_y",  nrow = 4,ncol = 3)+
  scale_fill_brewer(palette = "Spectral", direction = 1) +
  theme_bw()


#-----------------------------------------
#----STEP 13: Marine species ----
#-----------------------------------------
#Read in corrected thematic groups by Diederik
marine_groups_diederik<-read_sheet("1GeiSwNwojQa8YfaF29ebzDyNRwNJzqKiFCZinz2m_Bk", sheet="Data")%>%
  select(c("Thematic_group", "ID"))

#Read in species that were excluded from climate matching
marine<-read_sheet("1limsIJcHBtrJdKMMKwcK8_usdCOs0_v8u71MW1bC5w4",sheet="Data")

#Read in species that were excluded during climate matching
marine_16<-read_sheet("14DoltaPfvDA_ZMf_x7uDRLZhMLYFG8J-KtfGEdwBi4c", sheet="Data")
marine_7<-read_sheet("1Mkyt4mzp1lr-AEG8A86ERck_pXZihqHeO0tbtiPSyx0", sheet="Data")
marine_cm<-rbind(marine_16, marine_7)

#Bind them together with species that were excluded during climate matching
#Do a left join
marine_cm<-marine_cm%>%
  left_join(.,initial_CM_list, by="speciesKey")

marine<-rbind(marine, marine_cm)

#Remove marine species that are on list Irena (none)
marine<-marine%>%
  filter(!ID %in% Irena_plants_exclude)

#Remove unicellular algae (n = 6)
toremove<-marine%>%
  filter(kingdom!="Animalia")%>%
  filter(scientificName %in%c("Alexandrium monilatum (J.F.Howell) Balech, 1995",
                              "Marteilia refringens Grizel, Comps, Bonami, Cousserans, Duthoit & Le Pennec, 1974",
                              "Asteromphalus sarcophagus Wallich, 1860",
                              "Olisthodiscus luteus N.Carter, 1937",
                              "Isochrysis galbana Parke",
                              "Alexandrium leei Balech"))%>%
  pull(ID)

#This yields 875 marine species
marine<-marine%>%
  filter(!ID%in% toremove)%>%
  select(c(ID,scientificName, speciesKey,kingdom,phylum,class,order,family,genus,Source,Thematic_group,total_grid_cells))

#Add corrections Diederik for Thematic group
marine<-marine%>%
  left_join(marine_groups_diederik, by = "ID", suffix = c("_old", "_new")) %>%
  mutate(Thematic_group = ifelse(is.na(Thematic_group_old), Thematic_group_new, Thematic_group_old)) %>%
  select(-Thematic_group_old, -Thematic_group_new)%>%
  mutate(Thematic_group=case_when(Thematic_group=="plants"~"Plants",
                                  Thematic_group=="FW Verts/ Marine"~"FW Verts; Marine",
                                  Thematic_group=="FW Inverts/ Marine"~"FW Inverts; Marine",
                                  .default = Thematic_group))# Clean up the extra columns


table(marine$Thematic_group)


#----------------------------------------------------
#----STEP 12: create list for Ana ----
#----------------------------------------------------

cm_short<-climate_matching_weighed%>%
  filter(scenario=="2026-2050-A1FI") %>%
  group_by(scenario, speciesKey)%>%
  mutate(stand_potential_occupancy_range = cut(ranked_weighted_surface,
                                               breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                               labels = c(">0-10", ">10-20",">20-30",">30-40",">40-50",">50-60",">60-70",">70-80",">80-90",">90-100"),
                                               right = TRUE))%>%  # Then reverse the factor levels
  ungroup()%>%
  rename(standardised_potential_occupancy="ranked_weighted_surface")%>%
  select(c(ID,scientificName, speciesKey,kingdom,phylum,class,order,family,genus,Source,Thematic_group,total_grid_cells,scenario,percentage_climatematch_EU,standardised_potential_occupancy,stand_potential_occupancy_range))

list_ana<-bind_rows(cm_short, marine)

#----add original name as it was included on longlist----
longlist <- read_excel("data/input/Species_long_list_Ana.xlsx")%>%
  select(c("ID", "Long name"))

#Ends up with 3985 observations of 17 variables
list_ana<-list_ana%>%
  left_join(longlist,by="ID")%>%
  rename(original_name_longlist="Long name")%>%
  relocate(original_name_longlist, .after = ID)



#----------------------------------------------------------------------------------------------
#----STEP 13: create visuals to illustrate importance of standardised_potential_occuppancy
#----------------------------------------------------------------------------------------------

#Example 1: species Characodon lateralis (ID HS1915), climate match EU is 100% but standardised potential occupancy is 4.05%
#Example 2: species Clitarchus hookeri (ID HS2132), both climate match and potential occupancy are 100%



#----------------------------------------------------------------------------------------------
#----STEP 14: Calculate number of species left after applying cutoff
#----------------------------------------------------------------------------------------------

#Read in list of Ana that has corrected thematic groups
list_ana<-read_sheet("1Wo-SPMJaa8M3R3PnyOkPbVA5ethnWYLxNVxNyQ4_w-g", sheet="Data")

vertebrates<-filter(list_ana, Thematic_group_lumped=="Vertebrates")
table(vertebrates$stand_potential_occupancy_range, vertebrates$class)

vertebrates <- vertebrates %>%
  mutate(stand_potential_occupancy_range = factor(stand_potential_occupancy_range, levels = c(">0-10", ">10-20",">20-30",">30-40",">40-50",">50-60",">60-70",">70-80",">80-90",">90-100") ))

#Filter out vertebrates for which we have climate matching data
terrest_vertebrates<-filter(vertebrates, !is.na(stand_potential_occupancy_range))

#Count the number of species per category of standardized potential occupancy range
cat_terrest_vertebrates<-terrest_vertebrates%>%
  group_by(stand_potential_occupancy_range) %>%
  summarise(species_count = n())

#Count the number of species that would be left after applying each cutoff
cumulative_counts <-cat_terrest_vertebrates %>%
  mutate(cumulative_count = rev(cumsum(rev(species_count))))%>%
  mutate(cutoff = sub("-.*", "", stand_potential_occupancy_range))

ggplot(cumulative_counts, aes(x = cutoff, y = cumulative_count)) +
  geom_bar(stat = "identity", fill = "#f5f29f", color= "black") +
  geom_text(aes(label = cumulative_count), vjust = -0.5) +
  labs(title = "Number of species left at different cutoff levels",
       x = "Cutoff",
       y = "Number of species left") +
  theme_bw() +
  scale_x_discrete(limits = cumulative_counts$cutoff)




#Only keep vertebrates with a cutoff of >20%
vertebrates_selection<-filter(vertebrates, standardised_potential_occupancy>20 )

vertebrates_selection<-vertebrates_selection%>%
  mutate(class=case_when(order=="Siluriformes"~"Fishes",
                         order=="Esociformes"~"Fishes",
                         order=="Lepisosteiformes"~"Fishes",
                         order=="Cyprinodontiformes"~"Fishes",
                         order=="Characiformes"~"Fishes",
                         order=="Cypriniformes"~"Fishes",
                         order=="Perciformes"~"Fishes",
                         order=="Acipenseriformes"~"Fishes",
                         order=="Pleuronectiformes"~"Fishes",
                         order=="Atheriniformes"~"Fishes",
                         order=="Osmeriformes"~"Fishes",
                         order=="Salmoniformes"~"Fishes",
                         order=="Gasterosteiformes"~"Fishes",
                         order=="Anguilliformes"~"Fishes",
                         .default=class))

# Step 2: Count the number of species per class
species_counts <- vertebrates_selection %>%
  group_by(class) %>%
  summarise(count = n())

# Step 3: Create the bar plot
ggplot(species_counts, aes(x = class, y = count)) +
  geom_bar(stat = "identity", aes(fill = class), color = "black") +
  geom_text(aes(label = count), vjust = -0.5) +  # Optional: add counts above bars
  labs(title = "Number of Species per class for cutoff level > 20%",
       x = "Class",
       y = "Number of species") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()+
  ylim(0, max(species_counts$count) * 1.1)

#---------------------------------------------------------
#Per class and per percentage climate match EU

# Step 2: Count the number of species per class
species_counts <- vertebrates_selection %>%
  group_by(speciesKey)%>%
  mutate(EU_climatematch_range = cut(percentage_climatematch_EU,
                                     breaks = c(0,10,20,30,40,50,60,70,80,90,100.001),
                                     labels = c("0<-10", "10<-20","20<-30","30<-40","40<-50","50<-60","60<-70","70<-80","80<-90","90<-100"),
                                     right = TRUE))%>%
  group_by(class, EU_climatematch_range) %>%
  summarise(count = n())

# Step 3: Create the bar plot
ggplot(species_counts, aes(x = class, y = count)) +
  geom_bar(stat = "identity", aes(fill = EU_climatematch_range), color = "black") +
  labs(title = "Number of Species per class for cutoff level > 20%",
       x = "Class",
       y = "Number of species") +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  theme_bw()


