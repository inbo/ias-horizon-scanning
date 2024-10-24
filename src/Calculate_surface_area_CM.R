#---------------------------------------------------------------------
#-----------Calculate the surface area of each EU climate zone -------
#---------------------------------------------------------------------

#This script was used to create the file surface_area_climatezones.csv under the ./data/input folder

#--------Load packages---------
library(dplyr)
library(sf)

#-------------------Load data---------------------
# Load climate layers
load("./data/input/climate_matching/future.rda")
load("./data/input/climate_matching/legends.rda")

# Import legends
KG_Rubel_Kotteks_Legend <- legends$KG_A1FI

#Load shapefile of EU member states
memberstates<-read.csv2(here("./data/input/GRIIS_checklists.csv"))%>%
  pull(EU_state)

region_shape<-st_read("./data/spatial/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")%>%
  filter(ADMIN%in%c(memberstates,"Czech Republic"))%>%
  st_crop(xmin=-13, ymin=20, xmax=41, ymax= 81)%>%
  select(c(ADMIN, geometry))

st_agr(region_shape) = "constant"


#-----------------Prepare data for loop--------------
# Determine current and nearest future scenarios
scenarios <- c("2001-2025-A1FI",
               "2026-2050-A1FI")

# Create empty future_climate
future_climate <- data.frame() %>%
  mutate(scenario = "",
         KG_GridCode = as.integer(""),
         surface_area_climatezone=as.integer(""),
         total_surface_area=as.integer(""))# Calculate KG codes


#-----------Run loop to extract surface area of EU climate zones under each scenario--------
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

  gridcode_intersect_area<-gridcode_intersect%>%
    st_transform(crs = 3035) %>%
    summarise(area_km2 = sum(as.numeric(st_area(.))) / 1e6) %>%
    pull(area_km2)

  #Create an future_climate dataframe holding the gridcodes of the Köppen Geiger zones present in the region of interest and the associated future scenario
  for (g in unique(gridcode_intersect$GRIDCODE)) {

    #To get the surface area it is important that the shape file is in a projected CRS (not WGS84)
    gridcode_area<-gridcode_intersect%>%
      filter(GRIDCODE==g)%>%
      st_transform(crs = 3035) %>%
      summarise(area_km2 = sum(as.numeric(st_area(.))) / 1e6) %>%
      pull(area_km2)

    #Store gridcode, scenario, and associated occupied area
    future_climate <- future_climate %>%
      add_row(scenario = s,
              KG_GridCode = g,
              surface_area_climatezone=gridcode_area,
              total_surface_area= gridcode_intersect_area)
  }
}

future_climate <- future_climate%>%
  left_join(y = KG_Rubel_Kotteks_Legend,
            by = c("KG_GridCode" = "GRIDCODE"))

write.csv(future_climate, "./data/input/surface_area_climatezones.csv")
