
#####################################################################
######### Janek Bruker, Methods IV, ETH Zurich, FS 2021 #############
#####################################################################


#####################################################################
#####################################################################
################## DATA FOR PREDICTING OSV IN MALI ##################
#####################################################################
#####################################################################



### Load Packages and Data

library(here)
library(readr)
library(tidyverse)
library(sf)
library(lubridate)
library(ncdf4)
library(states)
library(countrycode)



### GEO DATA 
grid = st_read(here("data/priogrid_cellshp/priogrid_cell.shp"))


### PREDICTOR & OUTCOME DATA

# Static PRIO Grid 
prio_stat = read.csv(here("data/priogrid_static.csv"))

# SPEI 
spei = nc_open(here("data/spei01.nc"))

# ACLED events
acled <- read.csv(here("data/2012-01-01-2021-05-10-Burkina_Faso-Chad-Mali-Mauritania-Niger.csv"))

# EPR Groups 
epr2prio = read.csv(here("data/geoepr2priogrid.csv"))
epr = read.csv(here("data/EPR-2019.csv"))


### FOR DATA MANIPULATION

# GID to GW/COW codes
prio_gw = read.csv(here("data/gid-gw.csv"))
data("cowstates")




############## PREPARE SHAPE OF DATA ##############

# Unique GW/COW codes
cowstates = cowstates %>% 
  select(cowcode, cowc, country_name) %>% 
  distinct(.keep_all = TRUE)

# Merge info and add ISO3 format
prio_cntry = left_join(prio_gw, cowstates, by = c("gwno" = "cowcode")) %>% 
  mutate(iso3 = countrycode(cowc, origin = "cowc", destination = "iso3c")) %>% 
  select(gid, country_name, iso3, cowc, gwno)
# Warning due to ambiguous code, but not relevant for countries under study

# Filter Mali cells by ISO3 names
mali_grid = prio_cntry %>% 
  filter(iso3 == "MLI") %>% 
  left_join(grid, by = "gid") %>% 
  select(gid, geometry, ycoord, xcoord)

# Save geometric grid of Mali
mali_grid = st_as_sf(mali_grid)
save(mali_grid, file = here("newdata/mali_grid.RData"))
mali_grid = mali_grid %>% select(-c(ycoord, xcoord))

# Save cell IDs for Mali
mali_gid = data.frame(gid = mali_grid$gid)

# Time dimension: months in 2015/16/17/18 (15 for lags)
time_grid = expand.grid(year = c(2015, 2016, 2017, 2018), month_num = c(1:12),
                      gid = mali_grid$gid) %>% 
  mutate(month = ifelse(month_num < 10, 
                       paste(year, month_num, sep = "_0"),
                       paste(year, month_num, sep = "_"))) %>% 
  select(gid, month, year, month_num) %>% 
  distinct(.keep_all = T)

# Empty geodata in right format
mali_blank = left_join(mali_grid, time_grid, by = "gid") 


############## PRIO STATIC VARIABLES ##############

prio_mali = prio_stat %>% 
  select(gid,
    # Land use
    agri_gc, shrub_gc, herb_gc, barren_gc, 
    # Climate
    rainseas, 
    # Accessibility
    mountains_mean, ttime_mean)


predosv = left_join(mali_blank, prio_mali, by = "gid") %>% 
  # Indicate each month of rain season
  mutate(rain = ifelse(rainseas == month_num |
                       (rainseas+1) == month_num |
                         (rainseas+2) == month_num, 1, 0)) %>% 
  select(-rainseas) %>% 
  # Add monthly dummies
  mutate(january = ifelse(month_num == 1, 1, 0), 
         february = ifelse(month_num == 2, 1, 0),
         march = ifelse(month_num == 3, 1, 0), 
         april = ifelse(month_num == 4, 1, 0), 
         may = ifelse(month_num == 5, 1, 0), 
         june = ifelse(month_num == 6, 1, 0), 
         july = ifelse(month_num == 7, 1, 0), 
         august = ifelse(month_num == 8, 1, 0), 
         september = ifelse(month_num == 9, 1, 0),
         october = ifelse(month_num == 10, 1, 0), 
         november = ifelse(month_num == 11, 1, 0), 
         december = ifelse(month_num == 12, 1, 0))
  


################ ACLED EVENTS ###################

## Violence against civilians (ACLED) in Sahel
viol_mali = acled %>% filter((event_type == "Violence against civilians" |
                               event_type == "Explosions/Remote violence" |
                                event_type == "Battles" |
                                event_type == "Riots") &
                               (year > 2014 & year < 2019) &
                               iso3 == "MLI") %>% 
  mutate(date = dmy(event_date),
         month = ifelse(month(date) < 10, 
                       paste(year(date), month(date), sep = "_0"), 
                       paste(year(date), month(date), sep = "_")),
         violence = ifelse(event_type == "Violence against civilians" | 
                             event_type == "Explosions/Remote violence", 1, 0),
         battle = ifelse(event_type == "Battles", 1, 0),
         riot = ifelse(event_type == "Riots", 1, 0)) %>% 
  select(month, event_date, violence, battle, riot, iso3, event_type, sub_event_type, fatalities, notes, 
         longitude, latitude)

viol_mali = st_as_sf(x = viol_mali, 
                     coords = c("longitude", "latitude"))

viol_mali = st_set_crs(viol_mali, "+proj=longlat +datum=WGS84 +no_defs")


# Locate events in grid
viol_mali = st_join(mali_grid, viol_mali) %>% 
  group_by(gid, month) %>% 
  summarise(violence = ifelse(sum(violence, na.rm = T) > 0, 1, 0),
            battle = ifelse(sum(battle, na.rm = T) > 0, 1, 0),
            riot = ifelse(sum(riot, na.rm = T) > 0, 1, 0)) %>% 
  select(gid, month, violence, battle, riot)



# Merge event information to grid/month dataframe
viol_mali = as.data.frame(viol_mali) %>% select(gid, month, violence, battle, riot)

predosv = left_join(predosv, viol_mali, by = c("gid", "month")) %>% 
  mutate(violence = replace_na(violence, 0),
         battle = replace_na(battle, 0),
         riot = replace_na(riot, 0))




################ SPEI (Drought) Predictor ###################

### Get all information to build dataframe
spei_lon = ncvar_get(spei, "lon")
nlon = dim(spei_lon)

spei_lat = ncvar_get(spei, "lat")
nlat = dim(spei_lat)

spei_time = ncvar_get(spei, "time")
tunits = ncatt_get(spei, "time", "units")
ntime = dim(spei_time)

# Get all SPEI values and attributes
spei_arr = ncvar_get(spei, "spei")
spei_lname = ncatt_get(spei, "spei", "long_name")
spei_units = ncatt_get(spei, "spei", "units")
fillvalue = ncatt_get(spei, "spei", "_FillValue")


### Create SPEI Dataframe for Mali 

# replace fill values with NA
spei_arr[spei_arr==fillvalue$value] <- NA

lonlat <- as.matrix(expand.grid(spei_lon,spei_lat))

spei_vec = as.vector(spei_arr)

# Matrix with all SPEI values of last 48 months (latest date: 12/2018)
spei_mat <- matrix(spei_vec, nrow=nlon*nlat, ncol=ntime)[, (ntime-47):ntime]

spei_df = data.frame(cbind(lonlat, spei_mat))

varnames = sort(unique(time_grid$month))

names(spei_df) = c("lon", "lat", varnames)


### Filter SPEI for Mali 
spei_df = st_as_sf(spei_df, coords = c("lon", "lat"))
spei_df = st_set_crs(spei_df, "+proj=longlat +datum=WGS84 +no_defs")
mali_spei = st_join(mali_grid, spei_df)


mali_spei_long = as.data.frame(mali_spei)
mali_spei_long = mali_spei_long %>% 
  select(gid, contains("2015"), contains("2016"), contains("2017"), contains("2018")) %>% 
  pivot_longer(!gid, names_to = "month", values_to = "spei")


### Merge SPEI Data 
predosv = left_join(predosv, mali_spei_long, by = c("gid", "month"))



################ NEIGHBOURING OBSERVATIONS ###################

# Index of neighbouring grids
neighb_ind = st_intersects(mali_grid, mali_grid)

### Drop diagonal and add gid index 
neighb_gid = list()
nb_gid = data.frame(gid = rep(NA, nrow(mali_grid)), nbs = NA)

for(i in 1:nrow(mali_grid)){
  neighb_ind[[i]] = neighb_ind[[i]][neighb_ind[[i]] != i]
  neighb_gid[[i]] = mali_grid$gid[neighb_ind[[i]]]
  nb_gid$gid[i] = mali_grid$gid[i]
  nb_gid$nbs[i] = paste(neighb_gid[[i]], collapse = ", ")
}

names(neighb_gid) = mali_grid$gid




# Function to detect binary pattern of neighbors of X
nb_sum <- function(X, id, variable){
  n_events = sum(variable[id %in% neighb_gid[[as.character(X)]]])
  binary = ifelse(n_events > 0, 1, 0)
  return(binary)
}

# Function to compute average of neighbors of X
nb_mean <- function(X, id, variable){
  average = mean(variable[id %in% neighb_gid[[as.character(X)]]], na.rm = TRUE)
  return(average)
}


# Function to compute minimum of neighbors of X
nb_min <- function(X, id, variable){
  average = min(variable[id %in% neighb_gid[[as.character(X)]]], na.rm = TRUE)
  return(average)
}


# Loop over months and apply neighbor fun. to every cell 
for(m in unique(predosv$month)){
  
  # Violence in neighboring cells
  predosv$nb_viol[predosv$month == m] = as.numeric(sapply(X = subset(predosv, month == m)$gid, 
                                                               FUN = nb_sum, 
                                                               id = subset(predosv, month == m)$gid,
                                                               variable = subset(predosv, month == m)$violence))
  
  # Battles in neighboring cells
  predosv$nb_battle[predosv$month == m] = as.numeric(sapply(X = subset(predosv, month == m)$gid, 
                                                            FUN = nb_sum, 
                                                            id = subset(predosv, month == m)$gid,
                                                            variable = subset(predosv, month == m)$battle))
  
  # Riots in neighboring cells
  predosv$nb_riot[predosv$month == m] = as.numeric(sapply(X = subset(predosv, month == m)$gid, 
                                                            FUN = nb_sum, 
                                                            id = subset(predosv, month == m)$gid,
                                                            variable = subset(predosv, month == m)$riot))
  
  # Droughts in neighboring cells (minimum)
  # Minimum for SPEI indicates strongest drought potential
  predosv$nb_spei[predosv$month == m] = as.numeric(sapply(X = subset(predosv, month == m)$gid, 
                                                            FUN = nb_mean, 
                                                            id = subset(predosv, month == m)$gid,
                                                            variable = subset(predosv, month == m)$spei))
}





################ EPR Groups ###################

epr2prio = epr2prio %>% select(gid, gwgroupid) %>% 
  distinct(.keep_all = T)
  
epr_mali = epr %>% 
  filter(statename == "Mali" & 
         from == 2013)

epr_mali = left_join(epr_mali, epr2prio, by = "gwgroupid") %>%
  select(gid, group)

epr_mali = left_join(mali_gid, epr_mali, by = "gid") %>% 
  mutate(group = replace(group, group == "Blacks (Mande, Peul, Voltaic etc.)", "blacks"),
         group = replace(group, group == "Tuareg", "tuareg"), 
         group = replace(group, group == "Arabs/Moors", "arabs"), 
         value = 1) %>% 
  pivot_wider(names_from = "group", values_from = "value") %>% 
  mutate(blacks = replace_na(blacks, 0), 
         tuareg = replace_na(tuareg, 0),
         arabs = replace_na(arabs, 0)) %>% 
  select(-`NA`)


predosv = left_join(predosv, epr_mali, by = "gid")




################ MISSING VALUES ###################

nona = predosv %>% drop_na()
check_na = anti_join(as.data.frame(predosv), as.data.frame(nona))

# Replace missing SPEI value by mean of neighbors
predosv$spei[predosv$gid == check_na$gid & 
               predosv$month == check_na$month] = 
  nb_mean(X = as.character(check_na$gid), 
        id = subset(predosv, month == check_na$month)$gid,
        variable = subset(predosv, month == check_na$month)$spei)




################ LAGGING ###################

predosv = predosv %>% 
  group_by(gid) %>% 
  mutate(viol_l1 = lag(violence, n = 1, order_by = month),
         viol_l2 = lag(violence, n = 2, order_by = month),
         viol_l3 = lag(violence, n = 3, order_by = month),
         viol_l4 = lag(violence, n = 4, order_by = month),
         viol_l5 = lag(violence, n = 5, order_by = month),
         viol_l6 = lag(violence, n = 6, order_by = month),
         viol_l7 = lag(violence, n = 7, order_by = month),
         viol_l8 = lag(violence, n = 8, order_by = month),
         viol_l9 = lag(violence, n = 9, order_by = month),
         viol_l10 = lag(violence, n = 10, order_by = month),
         viol_l11 = lag(violence, n = 11, order_by = month),
         viol_l12 = lag(violence, n = 12, order_by = month),
         nb_viol_l1 = lag(nb_viol, n = 1, order_by = month),
         nb_viol_l2 = lag(nb_viol, n = 2, order_by = month),
         nb_viol_l3 = lag(nb_viol, n = 3, order_by = month),
         nb_viol_l4 = lag(nb_viol, n = 4, order_by = month),
         nb_viol_l5 = lag(nb_viol, n = 5, order_by = month),
         nb_viol_l6 = lag(nb_viol, n = 6, order_by = month),
         nb_viol_l7 = lag(nb_viol, n = 7, order_by = month),
         nb_viol_l8 = lag(nb_viol, n = 8, order_by = month),
         nb_viol_l9 = lag(nb_viol, n = 9, order_by = month),
         nb_viol_l10 = lag(nb_viol, n = 10, order_by = month),
         nb_viol_l11 = lag(nb_viol, n = 11, order_by = month),
         nb_viol_l12 = lag(nb_viol, n = 12, order_by = month),
         battle_l1 = lag(battle, n = 1, order_by = month),
         battle_l2 = lag(battle, n = 2, order_by = month),
         battle_l3 = lag(battle, n = 3, order_by = month),
         battle_l4 = lag(battle, n = 4, order_by = month),
         battle_l5 = lag(battle, n = 5, order_by = month),
         battle_l6 = lag(battle, n = 6, order_by = month),
         battle_l7 = lag(battle, n = 7, order_by = month),
         battle_l8 = lag(battle, n = 8, order_by = month),
         battle_l9 = lag(battle, n = 9, order_by = month),
         battle_l10 = lag(battle, n = 10, order_by = month),
         battle_l11 = lag(battle, n = 11, order_by = month),
         battle_l12 = lag(battle, n = 12, order_by = month),
         nb_battle_l1 = lag(nb_battle, n = 1, order_by = month),
         nb_battle_l2 = lag(nb_battle, n = 2, order_by = month),
         nb_battle_l3 = lag(nb_battle, n = 3, order_by = month),
         nb_battle_l4 = lag(nb_battle, n = 4, order_by = month),
         nb_battle_l5 = lag(nb_battle, n = 5, order_by = month),
         nb_battle_l6 = lag(nb_battle, n = 6, order_by = month),
         nb_battle_l7 = lag(nb_battle, n = 7, order_by = month),
         nb_battle_l8 = lag(nb_battle, n = 8, order_by = month),
         nb_battle_l9 = lag(nb_battle, n = 9, order_by = month),
         nb_battle_l10 = lag(nb_battle, n = 10, order_by = month),
         nb_battle_l11 = lag(nb_battle, n = 11, order_by = month),
         nb_battle_l12 = lag(nb_battle, n = 12, order_by = month),
         riot_l1 = lag(riot, n = 1, order_by = month),
         riot_l2 = lag(riot, n = 2, order_by = month),
         riot_l3 = lag(riot, n = 3, order_by = month),
         riot_l4 = lag(riot, n = 4, order_by = month),
         riot_l5 = lag(riot, n = 5, order_by = month),
         riot_l6 = lag(riot, n = 6, order_by = month),
         riot_l7 = lag(riot, n = 7, order_by = month),
         riot_l8 = lag(riot, n = 8, order_by = month),
         riot_l9 = lag(riot, n = 9, order_by = month),
         riot_l10 = lag(riot, n = 10, order_by = month),
         riot_l11 = lag(riot, n = 11, order_by = month),
         riot_l12 = lag(riot, n = 12, order_by = month),
         nb_riot_l1 = lag(nb_riot, n = 1, order_by = month),
         nb_riot_l2 = lag(nb_riot, n = 2, order_by = month),
         nb_riot_l3 = lag(nb_riot, n = 3, order_by = month),
         nb_riot_l4 = lag(nb_riot, n = 4, order_by = month),
         nb_riot_l5 = lag(nb_riot, n = 5, order_by = month),
         nb_riot_l6 = lag(nb_riot, n = 6, order_by = month),
         nb_riot_l7 = lag(nb_riot, n = 7, order_by = month),
         nb_riot_l8 = lag(nb_riot, n = 8, order_by = month),
         nb_riot_l9 = lag(nb_riot, n = 9, order_by = month),
         nb_riot_l10 = lag(nb_riot, n = 10, order_by = month),
         nb_riot_l11 = lag(nb_riot, n = 11, order_by = month),
         nb_riot_l12 = lag(nb_riot, n = 12, order_by = month),
         spei_l1 = lag(spei, n = 1, order_by = month),
         spei_l2 = lag(spei, n = 2, order_by = month),
         spei_l3 = lag(spei, n = 3, order_by = month),
         spei_l4 = lag(spei, n = 4, order_by = month),
         spei_l5 = lag(spei, n = 5, order_by = month),
         spei_l6 = lag(spei, n = 6, order_by = month),
         spei_l7 = lag(spei, n = 7, order_by = month),
         spei_l8 = lag(spei, n = 8, order_by = month),
         spei_l9 = lag(spei, n = 9, order_by = month),
         spei_l10 = lag(spei, n = 10, order_by = month),
         spei_l11 = lag(spei, n = 11, order_by = month),
         spei_l12 = lag(spei, n = 12, order_by = month),
         nb_spei_l1 = lag(nb_spei, n = 1, order_by = month),
         nb_spei_l2 = lag(nb_spei, n = 2, order_by = month),
         nb_spei_l3 = lag(nb_spei, n = 3, order_by = month),
         nb_spei_l4 = lag(nb_spei, n = 4, order_by = month),
         nb_spei_l5 = lag(nb_spei, n = 5, order_by = month),
         nb_spei_l6 = lag(nb_spei, n = 6, order_by = month),
         nb_spei_l7 = lag(nb_spei, n = 7, order_by = month),
         nb_spei_l8 = lag(nb_spei, n = 8, order_by = month),
         nb_spei_l9 = lag(nb_spei, n = 9, order_by = month),
         nb_spei_l10 = lag(nb_spei, n = 10, order_by = month),
         nb_spei_l11 = lag(nb_spei, n = 11, order_by = month),
         nb_spei_l12 = lag(nb_spei, n = 12, order_by = month))


### Drop observations from 2015 (after computing lag)
predosv = predosv %>% 
  filter(year > 2015)


### Assemble lagged variables within window of last m months 
# Sum (accumulation) for events
# Average for SPEI 

predosv = predosv %>%
  group_by(gid, month) %>%
  mutate(viol_m1 = viol_l1,
         viol_m2 = sum(viol_l1, viol_l2),
         viol_m3 = sum(viol_l1, viol_l2, viol_l3),
         viol_m4 = sum(viol_l1, viol_l2, viol_l3, viol_l4),
         viol_m5 = sum(viol_l1, viol_l2, viol_l3, viol_l4, viol_l5),
         viol_m6 = sum(viol_l1, viol_l2, viol_l3, viol_l4, viol_l5, viol_l6),
         viol_m7 = sum(viol_l1, viol_l2, viol_l3, viol_l4, viol_l5, viol_l6, viol_l7),
         viol_m8 = sum(viol_l1, viol_l2, viol_l3, viol_l4, viol_l5, viol_l6, viol_l7, viol_l8),
         viol_m9 = sum(viol_l1, viol_l2, viol_l3, viol_l4, viol_l5, viol_l6, viol_l7, viol_l8, viol_l9),
         viol_m10 = sum(viol_l1, viol_l2, viol_l3, viol_l4, viol_l5, viol_l6, viol_l7, viol_l8, viol_l9, viol_l10),
         viol_m11 = sum(viol_l1, viol_l2, viol_l3, viol_l4, viol_l5, viol_l6, viol_l7, viol_l8, viol_l9, viol_l10, viol_l11),
         viol_m12 = sum(viol_l1, viol_l2, viol_l3, viol_l4, viol_l5, viol_l6, viol_l7, viol_l8, viol_l9, viol_l10, viol_l11, viol_l12),
         nb_viol_m1 = nb_viol_l1,
         nb_viol_m2 = sum(nb_viol_l1, nb_viol_l2),
         nb_viol_m3 = sum(nb_viol_l1, nb_viol_l2, nb_viol_l3),
         nb_viol_m4 = sum(nb_viol_l1, nb_viol_l2, nb_viol_l3, nb_viol_l4),
         nb_viol_m5 = sum(nb_viol_l1, nb_viol_l2, nb_viol_l3, nb_viol_l4, nb_viol_l5),
         nb_viol_m6 = sum(nb_viol_l1, nb_viol_l2, nb_viol_l3, nb_viol_l4, nb_viol_l5, nb_viol_l6),
         nb_viol_m7 = sum(nb_viol_l1, nb_viol_l2, nb_viol_l3, nb_viol_l4, nb_viol_l5, nb_viol_l6, nb_viol_l7),
         nb_viol_m8 = sum(nb_viol_l1, nb_viol_l2, nb_viol_l3, nb_viol_l4, nb_viol_l5, nb_viol_l6, nb_viol_l7, nb_viol_l8),
         nb_viol_m9 = sum(nb_viol_l1, nb_viol_l2, nb_viol_l3, nb_viol_l4, nb_viol_l5, nb_viol_l6, nb_viol_l7, nb_viol_l8, nb_viol_l9),
         nb_viol_m10 = sum(nb_viol_l1, nb_viol_l2, nb_viol_l3, nb_viol_l4, nb_viol_l5, nb_viol_l6, nb_viol_l7, nb_viol_l8, nb_viol_l9, nb_viol_l10),
         nb_viol_m11 = sum(nb_viol_l1, nb_viol_l2, nb_viol_l3, nb_viol_l4, nb_viol_l5, nb_viol_l6, nb_viol_l7, nb_viol_l8, nb_viol_l9, nb_viol_l10, nb_viol_l11),
         nb_viol_m12 = sum(nb_viol_l1, nb_viol_l2, nb_viol_l3, nb_viol_l4, nb_viol_l5, nb_viol_l6, nb_viol_l7, nb_viol_l8, nb_viol_l9, nb_viol_l10, nb_viol_l11, nb_viol_l12),
         battle_m1 = battle_l1,
         battle_m2 = sum(battle_l1, battle_l2),
         battle_m3 = sum(battle_l1, battle_l2, battle_l3),
         battle_m4 = sum(battle_l1, battle_l2, battle_l3, battle_l4),
         battle_m5 = sum(battle_l1, battle_l2, battle_l3, battle_l4, battle_l5),
         battle_m6 = sum(battle_l1, battle_l2, battle_l3, battle_l4, battle_l5, battle_l6),
         battle_m7 = sum(battle_l1, battle_l2, battle_l3, battle_l4, battle_l5, battle_l6, battle_l7),
         battle_m8 = sum(battle_l1, battle_l2, battle_l3, battle_l4, battle_l5, battle_l6, battle_l7, battle_l8),
         battle_m9 = sum(battle_l1, battle_l2, battle_l3, battle_l4, battle_l5, battle_l6, battle_l7, battle_l8, battle_l9),
         battle_m10 = sum(battle_l1, battle_l2, battle_l3, battle_l4, battle_l5, battle_l6, battle_l7, battle_l8, battle_l9, battle_l10),
         battle_m11 = sum(battle_l1, battle_l2, battle_l3, battle_l4, battle_l5, battle_l6, battle_l7, battle_l8, battle_l9, battle_l10, battle_l11),
         battle_m12 = sum(battle_l1, battle_l2, battle_l3, battle_l4, battle_l5, battle_l6, battle_l7, battle_l8, battle_l9, battle_l10, battle_l11, battle_l12),
         nb_battle_m1 = nb_battle_l1,
         nb_battle_m2 = sum(nb_battle_l1, nb_battle_l2),
         nb_battle_m3 = sum(nb_battle_l1, nb_battle_l2, nb_battle_l3),
         nb_battle_m4 = sum(nb_battle_l1, nb_battle_l2, nb_battle_l3, nb_battle_l4),
         nb_battle_m5 = sum(nb_battle_l1, nb_battle_l2, nb_battle_l3, nb_battle_l4, nb_battle_l5),
         nb_battle_m6 = sum(nb_battle_l1, nb_battle_l2, nb_battle_l3, nb_battle_l4, nb_battle_l5, nb_battle_l6),
         nb_battle_m7 = sum(nb_battle_l1, nb_battle_l2, nb_battle_l3, nb_battle_l4, nb_battle_l5, nb_battle_l6, nb_battle_l7),
         nb_battle_m8 = sum(nb_battle_l1, nb_battle_l2, nb_battle_l3, nb_battle_l4, nb_battle_l5, nb_battle_l6, nb_battle_l7, nb_battle_l8),
         nb_battle_m9 = sum(nb_battle_l1, nb_battle_l2, nb_battle_l3, nb_battle_l4, nb_battle_l5, nb_battle_l6, nb_battle_l7, nb_battle_l8, nb_battle_l9),
         nb_battle_m10 = sum(nb_battle_l1, nb_battle_l2, nb_battle_l3, nb_battle_l4, nb_battle_l5, nb_battle_l6, nb_battle_l7, nb_battle_l8, nb_battle_l9, nb_battle_l10),
         nb_battle_m11 = sum(nb_battle_l1, nb_battle_l2, nb_battle_l3, nb_battle_l4, nb_battle_l5, nb_battle_l6, nb_battle_l7, nb_battle_l8, nb_battle_l9, nb_battle_l10, nb_battle_l11),
         nb_battle_m12 = sum(nb_battle_l1, nb_battle_l2, nb_battle_l3, nb_battle_l4, nb_battle_l5, nb_battle_l6, nb_battle_l7, nb_battle_l8, nb_battle_l9, nb_battle_l10, nb_battle_l11, nb_battle_l12),
         riot_m1 = riot_l1,
         riot_m2 = sum(riot_l1, riot_l2),
         riot_m3 = sum(riot_l1, riot_l2, riot_l3),
         riot_m4 = sum(riot_l1, riot_l2, riot_l3, riot_l4),
         riot_m5 = sum(riot_l1, riot_l2, riot_l3, riot_l4, riot_l5),
         riot_m6 = sum(riot_l1, riot_l2, riot_l3, riot_l4, riot_l5, riot_l6),
         riot_m7 = sum(riot_l1, riot_l2, riot_l3, riot_l4, riot_l5, riot_l6, riot_l7),
         riot_m8 = sum(riot_l1, riot_l2, riot_l3, riot_l4, riot_l5, riot_l6, riot_l7, riot_l8),
         riot_m9 = sum(riot_l1, riot_l2, riot_l3, riot_l4, riot_l5, riot_l6, riot_l7, riot_l8, riot_l9),
         riot_m10 = sum(riot_l1, riot_l2, riot_l3, riot_l4, riot_l5, riot_l6, riot_l7, riot_l8, riot_l9, riot_l10),
         riot_m11 = sum(riot_l1, riot_l2, riot_l3, riot_l4, riot_l5, riot_l6, riot_l7, riot_l8, riot_l9, riot_l10, riot_l11),
         riot_m12 = sum(riot_l1, riot_l2, riot_l3, riot_l4, riot_l5, riot_l6, riot_l7, riot_l8, riot_l9, riot_l10, riot_l11, riot_l12),
         nb_riot_m1 = nb_riot_l1,
         nb_riot_m2 = sum(nb_riot_l1, nb_riot_l2),
         nb_riot_m3 = sum(nb_riot_l1, nb_riot_l2, nb_riot_l3),
         nb_riot_m4 = sum(nb_riot_l1, nb_riot_l2, nb_riot_l3, nb_riot_l4),
         nb_riot_m5 = sum(nb_riot_l1, nb_riot_l2, nb_riot_l3, nb_riot_l4, nb_riot_l5),
         nb_riot_m6 = sum(nb_riot_l1, nb_riot_l2, nb_riot_l3, nb_riot_l4, nb_riot_l5, nb_riot_l6),
         nb_riot_m7 = sum(nb_riot_l1, nb_riot_l2, nb_riot_l3, nb_riot_l4, nb_riot_l5, nb_riot_l6, nb_riot_l7),
         nb_riot_m8 = sum(nb_riot_l1, nb_riot_l2, nb_riot_l3, nb_riot_l4, nb_riot_l5, nb_riot_l6, nb_riot_l7, nb_riot_l8),
         nb_riot_m9 = sum(nb_riot_l1, nb_riot_l2, nb_riot_l3, nb_riot_l4, nb_riot_l5, nb_riot_l6, nb_riot_l7, nb_riot_l8, nb_riot_l9),
         nb_riot_m10 = sum(nb_riot_l1, nb_riot_l2, nb_riot_l3, nb_riot_l4, nb_riot_l5, nb_riot_l6, nb_riot_l7, nb_riot_l8, nb_riot_l9, nb_riot_l10),
         nb_riot_m11 = sum(nb_riot_l1, nb_riot_l2, nb_riot_l3, nb_riot_l4, nb_riot_l5, nb_riot_l6, nb_riot_l7, nb_riot_l8, nb_riot_l9, nb_riot_l10, nb_riot_l11),
         nb_riot_m12 = sum(nb_riot_l1, nb_riot_l2, nb_riot_l3, nb_riot_l4, nb_riot_l5, nb_riot_l6, nb_riot_l7, nb_riot_l8, nb_riot_l9, nb_riot_l10, nb_riot_l11, nb_riot_l12),
         spei_m1 = spei_l1,
         spei_m2 = mean(c(spei_l1, spei_l2)),
         spei_m3 = mean(c(spei_l1, spei_l2, spei_l3)),
         spei_m4 = mean(c(spei_l1, spei_l2, spei_l3, spei_l4)),
         spei_m5 = mean(c(spei_l1, spei_l2, spei_l3, spei_l4, spei_l5)),
         spei_m6 = mean(c(spei_l1, spei_l2, spei_l3, spei_l4, spei_l5, spei_l6)),
         spei_m7 = mean(c(spei_l1, spei_l2, spei_l3, spei_l4, spei_l5, spei_l6, spei_l7)),
         spei_m8 = mean(c(spei_l1, spei_l2, spei_l3, spei_l4, spei_l5, spei_l6, spei_l7, spei_l8)),
         spei_m9 = mean(c(spei_l1, spei_l2, spei_l3, spei_l4, spei_l5, spei_l6, spei_l7, spei_l8, spei_l9)),
         spei_m10 = mean(c(spei_l1, spei_l2, spei_l3, spei_l4, spei_l5, spei_l6, spei_l7, spei_l8, spei_l9, spei_l10)),
         spei_m11 = mean(c(spei_l1, spei_l2, spei_l3, spei_l4, spei_l5, spei_l6, spei_l7, spei_l8, spei_l9, spei_l10, spei_l11)),
         spei_m12 = mean(c(spei_l1, spei_l2, spei_l3, spei_l4, spei_l5, spei_l6, spei_l7, spei_l8, spei_l9, spei_l10, spei_l11, spei_l12)),
         nb_spei_m1 = nb_spei_l1,
         nb_spei_m2 = mean(c(nb_spei_l1, nb_spei_l2)),
         nb_spei_m3 = mean(c(nb_spei_l1, nb_spei_l2, nb_spei_l3)),
         nb_spei_m4 = mean(c(nb_spei_l1, nb_spei_l2, nb_spei_l3, nb_spei_l4)),
         nb_spei_m5 = mean(c(nb_spei_l1, nb_spei_l2, nb_spei_l3, nb_spei_l4, nb_spei_l5)),
         nb_spei_m6 = mean(c(nb_spei_l1, nb_spei_l2, nb_spei_l3, nb_spei_l4, nb_spei_l5, nb_spei_l6)),
         nb_spei_m7 = mean(c(nb_spei_l1, nb_spei_l2, nb_spei_l3, nb_spei_l4, nb_spei_l5, nb_spei_l6, nb_spei_l7)),
         nb_spei_m8 = mean(c(nb_spei_l1, nb_spei_l2, nb_spei_l3, nb_spei_l4, nb_spei_l5, nb_spei_l6, nb_spei_l7, nb_spei_l8)),
         nb_spei_m9 = mean(c(nb_spei_l1, nb_spei_l2, nb_spei_l3, nb_spei_l4, nb_spei_l5, nb_spei_l6, nb_spei_l7, nb_spei_l8, nb_spei_l9)),
         nb_spei_m10 = mean(c(nb_spei_l1, nb_spei_l2, nb_spei_l3, nb_spei_l4, nb_spei_l5, nb_spei_l6, nb_spei_l7, nb_spei_l8, nb_spei_l9, nb_spei_l10)),
         nb_spei_m11 = mean(c(nb_spei_l1, nb_spei_l2, nb_spei_l3, nb_spei_l4, nb_spei_l5, nb_spei_l6, nb_spei_l7, nb_spei_l8, nb_spei_l9, nb_spei_l10, nb_spei_l11)),
         nb_spei_m12 = mean(c(nb_spei_l1, nb_spei_l2, nb_spei_l3, nb_spei_l4, nb_spei_l5, nb_spei_l6, nb_spei_l7, nb_spei_l8, nb_spei_l9, nb_spei_l10, nb_spei_l11, nb_spei_l12)))

# Drop lags (after computing time-window variables)
predosv = predosv %>% 
  select(-c(contains("viol_l"), contains("riot_l"), contains("battle_l"), contains("spei_l")))

        

################ KEEP DATA FOR FORECASTING ###################

### Drop non-lagged dynamic variables (non-usable for forecasting), year and month number
predosv = predosv %>% 
  select(-c(riot, battle, spei, nb_viol, nb_battle, nb_riot,
            nb_spei, year, month_num))


################ SAVE DATAFRAME ###################

# Save geocoded dataset separately 
geo_predosv = predosv
save(geo_predosv, file = here("newdata/geo_mali_osv.RData"))

# Save dataset without geometric attributes
predosv = as.data.frame(predosv) %>% 
  select(-geometry)
save(predosv, file = here("newdata/mali_osv.RData"))


rm(list = ls())

