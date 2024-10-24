#----------------------------------------------------------------------------------
#Run climate matching for 286 species that occupy 3-5 grid cells in the EU
#-----------------------------------------------------------------------------------
#This script has the same functionality as the Climate_matching_exercise script, with the difference that it
#is run for the species that occupy 3-5 grid cells in the EU

#------------Load packages------------
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
library(googlesheets4)


#----------Load data----------
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

#Get list with species occupying 3-5 grid cells
total_gridcells<-read_sheet("1pyutiXpEk3Y4-t3Nw3tg-SoZRqz5PSL3lEf7GGgL1HY", sheet="Data")
climate_matching<-filter(total_gridcells, total_grid_cells>2, total_grid_cells<6)


#---- Prepare blocks of taxonkeys for download ----
#Indicate basis of record for download
basis_of_record <- paste0(
  "OBSERVATION;",
  "HUMAN_OBSERVATION;",
  "MATERIAL_SAMPLE;",
  "LITERATURE;",
  "PRESERVED_SPECIMEN;",
  "UNKNOWN;",
  "MACHINE_OBSERVATION")

#Extract taxonkeys as vector (length 4584 or 4555)
taxonkeys<-climate_matching%>%
  pull(speciesKey)%>%
  unique()

#Count how many blocks we need if we want them to contain around 150 keys
nblocks<-ceiling(length(taxonkeys)/150)

# There are too many taxonkeys to run occ_search all at once so we need to loop over the blocks
# Note: coordinateuncertainty not included
total_records<-data.frame(record_number=rep(as.numeric(NA), nblocks),
                          start_index=rep(as.numeric(NA), nblocks),
                          end_index=rep(as.numeric(NA), nblocks))
for (i in 1:nblocks) {
  # Calculate start and end index for the current chunk
  start_idx <- (i - 1) * 150 + 1
  end_idx <- min(i * 150, length(taxonkeys))

  # Extract the current chunk of taxonkeys and put in right format
  chunk <- taxonkeys[start_idx:end_idx]%>%
    paste0(.,collapse=";")

  #Count number of records (this will be higher than the real amount because we can't filter on identificationverif... and coordinate uncertainty)
  count_chunk<-occ_count(taxonKey=chunk,
                         occurrenceStatus = "PRESENT",
                         basisOfRecord = basis_of_record,
                         year="1901,2024",#Between the years 1901-2024
                         hasGeospatialIssue = FALSE,
                         hasCoordinate = TRUE
  )

  total_records[i,1]<-count_chunk
  total_records[i,2]<-start_idx
  total_records[i,3]<-end_idx
}

# Total records
print(paste0("The expected number of records that will be downloaded is: ",sum(total_records$record_number)))

#-----------------------------------------------------------------
#NOTE: we'll do climate matching for around 10 000 000 records at a time
#----------------------------------------------------------------

#Put groups of taxonkeys together so they'll yield around 10 000 000 records
total_records <- total_records %>%
  arrange(record_number)

#Initialize variables for summation
target_records <- 10000000
current_sum <- 0
current_group <- data.frame()
groups <- list()

# Step 3: Iterate through the sorted data
for (i in 1:nrow(total_records)) {
  # Check if adding the current record_number would exceed the target
  if (current_sum + total_records$record_number[i] < target_records) {
    current_sum <- current_sum + total_records$record_number[i]
    current_group <- rbind(current_group, total_records[i, ])
  } else {
    # If the current group is not empty, save it
    if (nrow(current_group) > 0) {
      groups[[length(groups) + 1]] <- current_group
    }

    # Reset for the new group
    current_group <- total_records[i, , drop = FALSE]  # Start new group with the current row
    current_sum <- total_records$record_number[i]  # Start sum with the current record_number
  }
}

# If there's any remaining group not yet added
if (nrow(current_group) > 0) {
  groups[[length(groups) + 1]] <- current_group
}


#--------Prepare empty dataframes for loop--------
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


#--------Climate matching in loop--------
system.time({
  for (i in seq_along(groups)) {

    #Add group to logs
    climate_match_logs[i,1]<-i

    #Get dataframe of the right group
    Group <- groups[[i]]

    # Initialize an empty vector to store the extracted taxonkeys
    extracted_taxonkeys <- c()

    #  Loop through each row in the dataframe
    for (r in 1:nrow(Group)) {
      # Extract values from the vector based on start_index and end_index
      values <- taxonkeys[Group$start_index[r]:Group$end_index[r]]

      # Append the extracted values to the vector
      extracted_taxonkeys<- c(extracted_taxonkeys, values)
    }

    #Print info
    print(paste0("Group ",i,": ",length(extracted_taxonkeys)," species"))

    #Add number of taxonkeys to logs
    climate_match_logs[i,2]<-length(extracted_taxonkeys)

    #----------------------------------------
    #----STEP 3: Download occurrence data----
    #----------------------------------------

    #Indicate basis of record for download
    basis_of_record <- c(
      "OBSERVATION",
      "HUMAN_OBSERVATION",
      "MATERIAL_SAMPLE",
      "LITERATURE",
      "PRESERVED_SPECIMEN",
      "UNKNOWN",
      "MACHINE_OBSERVATION")

    # Identification_verification_status to discard
    identificationverificationstatus_to_discard<- c(
      "unverified",
      "unvalidated",
      "not validated",
      "under validation",
      "not able to validate",
      "control could not be conclusive due to insufficient knowledge",
      "Control could not be conclusive due to insufficient knowledge",
      "0",
      "uncertain",
      "unconfirmed",
      "Douteux",
      "Invalide",
      "Non r\u00E9alisable",
      "verification needed" ,
      "Probable",
      "unconfirmed - not reviewed",
      "validation requested")

    #Perform download
    gbif_download <- occ_download(
      pred_in("taxonKey", extracted_taxonkeys),
      pred("hasGeospatialIssue", FALSE), #Remove default geospatial issues
      pred("hasCoordinate", TRUE), # Keep only records with coordinates
      pred("occurrenceStatus", "PRESENT"),
      pred_in("basisOfRecord", basis_of_record),
      pred_gt("year", 1900),
      pred_or(
        pred_lt("coordinateUncertaintyInMeters",50000),
        pred_isnull("coordinateUncertaintyInMeters"))
    )

    #Follow the status of the download
    occ_download_wait(gbif_download)

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
    print(paste0("Occurrences for ",length(unique(cm_records$speciesKey))," different species out of ",length(extracted_taxonkeys)," were downloaded.")) #4312

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
climate_matching_EU<-read_sheet("1gM6j7XS4KEQY8hZPQOJXDiysp3B3IYV3FHzYCcTbJMw",sheet="Data")
length(unique(climate_matching_EU$speciesKey)) #254

climate_match_logs<-read_sheet("1TXEIq6roqzi_0CtfGfYzU7NKWRxSuxiZyN3xjn7qGhk", sheet="Data") %>%
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
  geom_magnify(from = c(xmin =0, xmax = 100, ymin = 0, ymax = 15),
               to = c(xmin =1000000, xmax = 2000000, ymin = 5, ymax = 12),
               axes="xy")+
  theme_bw()


#-----------------------------------------------------------------
#----STEP 10: Process and filter data to obtain final dataset ----
#-----------------------------------------------------------------

#---Load datasets needed to append information---

#Initial list used to extract taxonkeys for climate matching of species occupying 3-5 grid cells
initial_CM_list<-read_sheet("1XfE2LmyOOrgl-cxk0EmzWIG_9DHNtFK1r9PNXuPbKJU", sheet="Data")%>%
  filter(total_grid_cells>2 & total_grid_cells<6)
length(unique(initial_CM_list$speciesKey)) #286


#--------------Add taxonomic information-------------

climate_matching_EU<-climate_matching_EU%>%
  left_join(.,initial_CM_list, by="speciesKey")

length(unique(climate_matching_EU$speciesKey)) #254


#---------Remove thematic group = Marine from list as I forgot to do that beforehand
marine<-climate_matching_EU%>%
  filter(Thematic_group=="Marine")%>%
  select(c(1,10:21))
length(unique(marine$ID)) #5

climate_matching_EU<-climate_matching_EU%>% #249 species that were climate matched
  filter(Thematic_group!="Marine")
length(unique(climate_matching_EU$ID))

#---------Prepare final datasets------------
#Remove species that had less than 20 records to climate match with
climate_matching_final<-climate_matching_EU%>% #246 left, 3 removed
  filter(n_totaal>20)
length(unique(climate_matching_final$ID))

#Calculate a weighted surface area and percentage climate match with general EU
#max: 1951756 in 2001-2025; the max in 2026-2050 = 1843108
#For climate matching of 0-2 cells we used 1951756 as max for 2001-2025 and 1843108 for 2026-2050 so we'll do the same here
climate_matching_weighed<-climate_matching_final%>%
  mutate(weighted_surface=perc_climate*surface_area_climatezone/100)%>%
  group_by(scenario, speciesKey,n_totaal)%>%
  summarize(sum_weighted_surface=sum(weighted_surface),
            percentage_climatematch_EU = sum(perc_climate, na.rm = TRUE))%>%
  ungroup()%>%
  mutate(max_sum_weighted_surface= case_when(scenario=="2001-2025-A1FI"~1951756,
                                             scenario=="2026-2050-A1FI"~1843108))%>% #This is the maximum possible value, equal to 100% in Cfb in each scenario
  mutate(ranked_weighted_surface=sum_weighted_surface/max_sum_weighted_surface*100)%>%
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
           forcats::fct_rev()) %>%  # Then reverse the factor levels
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
           forcats::fct_rev(),
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


#-----------------------------------------
#----STEP 12: Marine species ----
#-----------------------------------------

#Read in species that were excluded during climate matching
marine_cm<-read_sheet("1moGtxLqtoSGhKpWcVs_L96FHudRGWMNIlKsfp3cxpy8", sheet="Data")

#Do a left join to add taxonomic information
marine_cm<-marine_cm%>%
  left_join(.,initial_CM_list, by="speciesKey")

##Bind them together with species that were excluded because they had thematic_group "Marine"
marine<-rbind(marine, marine_cm)

#Remove unicellular algae (n = 3)
toremove<-marine%>%
  filter(kingdom!="Animalia")%>%
  filter(scientificName %in%c("Alexandrium catenella (Whedon & Kofoid) Balech",
                              "Dicroerisma psilonereiella F.J.R.Taylor & Cattell",
                              "Pseudochattonella verruculosa (Y.Hara & Chihara) Hosoi-Tanabe, D.Honda, S.Fukaya, Y.Inagaki & Y.Sako"
                              ))%>%
  pull(ID)

#This yields 34 marine species (32 removed because of 50% occurrences, 5 included but had thematic group marine, 3 removed because unicellular algae)
marine<-marine%>%
  filter(!ID%in% toremove)%>%
  distinct(speciesKey, .keep_all=TRUE)%>%
  select(c(ID,scientificName, speciesKey,kingdom,phylum,class,order,family,genus,Source,Thematic_group,total_grid_cells))


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
longlist <- readxl::read_excel("data/input/Species_long_list_Ana.xlsx")%>%
  select(c("ID", "Long name"))

#Ends up with 3985 observations of 17 variables
list_ana<-list_ana%>%
  left_join(longlist,by="ID")%>%
  rename(original_name_longlist="Long name")%>%
  relocate(original_name_longlist, .after = ID)

list_ana$Info<-"Included on shortlist"
#----------------------------------------------------------
#------------------STEP 13: Add info column--------------------
#----------------------------------------------------------

#----------Create vectors with IDs----------

#all keys before climate matching (including marine)
climate_matching_IDs<-initial_CM_list%>%
  pull(ID) #286


#all IDs right after climate matching
climate_matching_EU_IDs<-climate_matching_EU%>%
  distinct(ID)%>%
  pull(ID) #249

#All IDs after filtering out those with 20 or less occurences
climate_matching_final_IDs<-climate_matching_EU%>% #246 left
  filter(n_totaal>20)%>%
  distinct(ID)%>%
  pull(ID)

#IDs that were removed because they had less than 20 occurrences: 3
CM_removed_20records_IDs<-climate_matching_EU_IDs[!climate_matching_EU_IDs%in%climate_matching_final_IDs]

#All IDs on shortlist
shortlist_IDs<-list_ana%>%
  pull(ID)#

#unicellular algae IDs: 3
unicellular_algae_IDs<-marine_cm%>%
  filter(kingdom!="Animalia")%>%
  filter(scientificName %in%c("Alexandrium catenella (Whedon & Kofoid) Balech",
                              "Dicroerisma psilonereiella F.J.R.Taylor & Cattell",
                              "Pseudochattonella verruculosa (Y.Hara & Chihara) Hosoi-Tanabe, D.Honda, S.Fukaya, Y.Inagaki & Y.Sako"
  ))%>%
  pull(ID)

#All IDs that were climate matched but were removed from shortlist (to check afterwards if it adds up)
climate_matching_not_on_shortlist_IDs<-climate_matching_IDs[!climate_matching_IDs%in%shortlist_IDs]#6
rlang::is_true(length(unicellular_algae_IDs)+length(CM_removed_20records_IDs)==length(climate_matching_not_on_shortlist_IDs))


#-----------Add info column to final dataset----------------
initial_CM_list<-initial_CM_list%>%
  left_join(longlist,by="ID")%>%
  rename(original_name_longlist="Long name")%>%
  relocate(original_name_longlist, .after = ID)

initial_CM_list<-initial_CM_list[,colnames(list_ana[,1:13])]

check<-initial_CM_list%>%
  mutate(Info = case_when(ID %in% CM_removed_20records_IDs ~ "Excluded: only 20 or less reliable occurrence records globally",
                          ID %in% unicellular_algae_IDs~ "Excluded: unicellular algae",
                          ID %in% shortlist_IDs~ "Included on shortlist",
                          .default=NA))%>%
  filter(Info!="Included on shortlist")

list_ana<-bind_rows(list_ana, check)

list_ana<-list_ana%>%
  mutate(Thematic_group=case_when(Thematic_group=="FW Inverts/ Marine"~"FW Inverts; Marine",
                                  .default=Thematic_group))

#Fix thematic groups
#Load dataset with fixed groups
thematic_groups<-read_sheet("1k6pSd9quV7moZ7K9sLjghbQJRihZxPP48_5obcRYS4A", sheet="Data all")%>%
  select(ID,Thematic_group)

list_ana_final<-list_ana%>%
  left_join(thematic_groups, by = "ID", suffix = c("_old", "_new")) %>%
  mutate(Thematic_group = ifelse(!is.na(Thematic_group_new), Thematic_group_new, Thematic_group_old)) %>%
  select(-Thematic_group_old, -Thematic_group_new)%>%
  relocate(Thematic_group, .after = Source)

ss0<-gs4_create("species_occupying_3_5_gridcells", sheets="Data")
write_sheet(list_ana_final, ss=ss0, sheet="Data")
