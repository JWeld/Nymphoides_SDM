library(biomod2)
library(tidyverse)
library(readxl)
library(sf)
library(corrplot)
library(usdm) #masks select from dplyr!
library(raster)
library(terra)
library(tidyterra)
library(doParallel)
library(rasterVis)
library(mda)
library(earth)
library(stringr)
library(ggtext)
library(maxnet)
library(xgboost)

#1. CORRELATIONS OF EXPLANATORY VARIABLES  ##########
Lakes_info <-
  st_read("Data/Lakes_info_SWE_proportions_CLC.shp")
colnames(Lakes_info)

#Import explanatory variables
SWE_transformed <- read_delim(
  "Data/SWE_transformed.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

SWE_transformed <-
  dplyr::select(
    SWE_transformed,
    "OBJECTID",
    "ApproxID",
    "X",
    "Y",
    "Shape_leng",
    "Mud",
    "Nb_lakesup",
    "Temp",
    "Depth_mean",
    "gdd",
    "pH",
    #"Alk",
    "AMM",
    "NIT",
    "SULF",
    "NITRO",
    "PopdensSCB",
    "prop_AGRI",
    "prop_OWATE"
  )


# Variable for temperature tolerance boundary
SWE_transformed$Temp_limit <- ifelse(SWE_transformed$Temp >= 5, 1, 0)
#Leave this as numeric for correlation analysis but convert to categorical later
#Drop temp
SWE_transformed$Temp <- NULL

#run hol homes raster.R to get this data on holiday home density (already log transformed +1)
#holiday_home <- read_csv("Data/holiday_home.csv")
#SWE_transformed$HH <- holiday_home$hh

#Use measured calcium to predict values everywhere
slu_mvm_kemi_data <- read_excel("Data/slu_mvm_kemi_data.xlsx",
                                col_types = c("text", "text", "numeric",
                                              "numeric", "text", "text", "text",
                                              "text", "text", "numeric", "numeric",
                                              "numeric", "numeric", "numeric"))

library(gstat)
library(sp)

calc <- slu_mvm_kemi_data %>% select(N,E, Ca_mekv_l) %>% drop_na()
calc <- calc  %>% rename(Y = N, X = E, calcium = Ca_mekv_l)

coordinates(calc) <- ~X + Y
Lakes_info_df <- select(Lakes_info, OBJECTID, X, Y) %>% as.data.frame()
coordinates(Lakes_info_df) <- ~X + Y

# use measured calcium values to predict values at test_sites locations
idw_results <- idw(calcium ~ 1, calc, newdata = Lakes_info_df)

# Extract the results as a data frame
predicted_values <- data.frame(
  OBJECTID = Lakes_info_df$OBJECTID,
  E = coordinates(Lakes_info_df)[, 1],
  N = coordinates(Lakes_info_df)[, 2],
  Ca = idw_results$var1.pred  # var1.pred is the column for predicted values
)

SWE_transformed <- left_join(SWE_transformed, predicted_values, by = "OBJECTID")

#drop NA
SWE_transformed <- drop_na(SWE_transformed)

## CALCULATE CORRELATIONS BETWEEN VARIABLES ####

corrdat <-
  dplyr::select(SWE_transformed,-c(OBJECTID, ApproxID, X, Y)) %>% drop_na %>% cor()
corrplot(corrdat)

#Not suprisingly, some likely have a colinearity problem
#the following choose variables that are not highly correlated to each other
r1 <-dplyr::select(SWE_transformed,-c(OBJECTID, ApproxID, X, Y, N, E)) %>% as.data.frame()
vif(r1)

#Use vifcor function to select from highly correlated varibles

nocorrvar <-
 vifstep(r1, th = 5)#, keep = c('Ca'))  # needs usdm

#nocorrvar <-
#  vifcor(r1, th = 0.7)#, keep = c('gdd'))  # needs usdm

predictors <-
  as.character(nocorrvar@results$Variables)  # the names of the remaining variables
predictors
r2 <- dplyr::select(r1, one_of(predictors))
# [1] "Shape_leng" "Mud"        "Nb_lakesup" "Temp"       "Depth_mean" "pH"         "Alk"        "prop_AGRI"  "prop_BUILT"
# [10] "prop_OWATE" "prop_WETLA" "HH"        
cor(drop_na(r2))
corrplot(cor(drop_na(r2)))
vif(drop_na(r2))
#Take the variables suggested by vifstep further in modelling ########


#give dplyr select function back
select <- dplyr::select

#2. IMPORT SHAPEFILE ##########################################################


# Lakes_info <-
#   st_read("Data/Lakes_info_SWE_proportions_CLC.shp")
# colnames(Lakes_info)

#Take columns found in SWE15.csv, including transformations
OBJECTID <- Lakes_info$OBJECTID
Shape_leng <- log(Lakes_info$Shape_Leng)
ApproxID <- Lakes_info$ApproxID
Mud <- log(Lakes_info$Mud + 1)
Peat <- log(Lakes_info$Peat + 1)
Nb_lakesup <-log(as.numeric(Lakes_info$Nb_lakesup) + 1)
Depth_mean <-log(Lakes_info$Depth_mean +1)
gdd <- log(Lakes_info$gdd)
nT0 <- log(Lakes_info$nT0)
pH <- as.numeric(Lakes_info$pH)
Alk <- log((as.numeric(Lakes_info$Alk)+1))
NITRO <- log(Lakes_info$NITRO + 1)
SULF <- log(Lakes_info$SULF + 1)
NIT <- log(Lakes_info$NIT + 1)
AMM <- log(Lakes_info$AMM + 1)
Temp <- Lakes_info$Temp
Ca <- log(predicted_values$Ca +1) #or just log?
#temp limit based on 16 degree july/ 5 degree annual
Lakes_info$Temp_limit <- ifelse(Lakes_info$Temp >= 5, 1, 0)
Temp_limit <- as.factor(Lakes_info$Temp_limit)
#Temp_limit <- Lakes_info$Temp_limit

PopdensSCB <- log(Lakes_info$PopdensSCB + 1)
prop_AGRI <- log(Lakes_info$prop_AGRI + 1)
prop_BUILT <- log(Lakes_info$prop_BUILT + 1)
prop_OWATE <- log(Lakes_info$prop_OWATE + 1)
prop_WETLA <- log(Lakes_info$prop_WETLA + 1)
prop_FORES <- log(Lakes_info$prop_WETLA + 1)
E <- Lakes_info$X
N <- Lakes_info$Y
X <- Lakes_info$X
Y <- Lakes_info$Y

#Combine transformed data from shapefiles into a dataframe, taking account of VIF earlier

Exp_vars <-
  data.frame(
    OBJECTID,
    Shape_leng,
    ApproxID,
    Mud,
    Peat,
    Nb_lakesup,
    Depth_mean,
    NITRO,
    Temp,
    Temp_limit,
    gdd,
    pH,
    Alk,
    AMM,
    PopdensSCB,
    prop_AGRI,
    prop_BUILT,
    prop_OWATE,
    prop_WETLA,
    Ca,
    E,
    N
  ) %>% as.data.frame() %>% drop_na()


Exp_vars$ApproxID <- as_factor(Exp_vars$ApproxID)
Exp_vars$OBJECTID <- as_factor(Exp_vars$OBJECTID)

#keep the variable identified in vifstep
Exp_vars <- Exp_vars %>% dplyr::select(any_of(predictors), E, N, OBJECTID, ApproxID)

#3. MAKE ENV DATA INTO RASTERS #################################################
#Make rasterstack of transformed environmental variables ##
shape.file <- read_sf("Data/Lakes_info_SWE_proportions_CLC.shp")
st_crs(shape.file)
shape.file <-
  st_transform(shape.file, crs = 3006) #SWEREF99
proj <- crs(shape.file)
#make empty raster to fill with data
r <-
  rast(
    xmin = 279123,
    xmax = 916264,
    ymin = 6140440,
    ymax = 7670710,
    crs = proj
  )
crs(r) <- proj

#make dataframe for each variable and use rasterize function to convert to raster
Exp_vars2 <-
  Exp_vars %>% dplyr::select(OBJECTID, ApproxID, E, N, everything()) %>% drop_na()
Exp_vars2 <- rename(Exp_vars2, X = E, Y = N)
#Exp_vars2 <- left_join(Exp_vars2, select(SWE_chem, -c(E,N)), by = "OBJECTID")

make.rast <- function(x,y){
  dat.x <- Exp_vars2[,c("X","Y", x)]
  x.vec <- vect(dat.x, geom=c("X", "Y"), crs=proj)
  x.rast <- rasterize(x.vec, r, field = x, fun = "last")
  crs(x.rast) <- proj
  plot(x.rast)
  writeRaster(x.rast, y,overwrite=TRUE)
}

 make.rast("Shape_leng","Shape_leng.grd")
 make.rast("Mud","Mud.grd")
# make.rast("Peat","Peat.grd")
 make.rast("Nb_lakesup","Nb_lakesup.grd")
 make.rast("Depth_mean","Depth_mean.grd")
 make.rast("gdd","gdd.grd")
# make.rast("Temp","Temp.grd")
 make.rast("Temp_limit","Temp_limit.grd")
 make.rast("pH","pH.grd")
# make.rast("Alk","Alk.grd")
# make.rast("AMM","AMM.grd")
# make.rast("NITRO","NITRO.grd")
 make.rast("PopdensSCB","PopdensSCB.grd")
 make.rast("prop_AGRI","prop_AGRI.grd")
# make.rast("prop_BUILT","prop_BUILT.grd")
 make.rast("prop_OWATE","prop_OWATE.grd")
 make.rast("Ca","Ca.grd")


rasterlist <-
  c(
    "Shape_leng.grd",
    "Mud.grd",
    "Peat.grd",
    "Nb_lakesup.grd",
    "Depth_mean.grd",
    "NITRO.grd",
    "gdd.grd",
    "Temp.grd",
    "pH.grd",
    "Alk.grd",
    "AMM.grd",
    "PopdensSCB.grd",
    "prop_AGRI.grd",
    "prop_BUILT.grd",
    "prop_OWATE.grd",
    "prop_WETLA.grd",
    "Temp_limit.grd",
    "Ca.grd"
  )


# Extract the names of each raster file without the ".grd" suffix
raster_names <- sapply(rasterlist, function(x) str_remove(basename(x), "\\.grd"))

# Filter Rasterlist to include only those rasters whose names (without ".grd") match the names in `B`
FilteredRasterlist <- rasterlist[raster_names %in% predictors]

envdat.r <- rast(FilteredRasterlist)
crs(envdat.r) <- proj
names(envdat.r)
envdat.sr = stack(envdat.r)
names(envdat.sr) <- FilteredRasterlist
names(envdat.sr)
plot(envdat.sr)


#4. Import Nymphoides peltata data #############################################

#From MVM #####
#slu_mvm_data <- read_csv("Data/slu_mvm_macro_data.csv") #SWEREF99

slu_mvm_data <- read_csv("Data/slu_mvm_data2.csv") #SWEREF99

slu_mvm_data <- slu_mvm_data[, -1]
slu_mvm_data <-
  slu_mvm_data %>% rename(E = "Stationskoordinat E/Y", N = "Stationskoordinat N/X")

Nymphoides_mvm <-
  dplyr::select(slu_mvm_data,
                `MD-MVM Id`,
                N,
                E,
                `Provtagningens startdatum`,
                Taxonnamn) %>%
  rename(ApproxID = `MD-MVM Id`, Date = `Provtagningens startdatum`)

Nymphoides_mvm <- drop_na(Nymphoides_mvm)

Nymphoides.pa <-
  Nymphoides_mvm %>% group_by(ApproxID) %>% mutate(Nymphoides_pa = {
    if (any(Taxonnamn == "Nymphoides peltata"))
      "1"
    else
      "0"
  })

Nymphoides.pa$Nymphoides_pa <- as.numeric(Nymphoides.pa$Nymphoides_pa)
Nymphoides.pa <-
  dplyr::select(Nymphoides.pa, `ApproxID`, N, E, Nymphoides_pa) %>% distinct() %>% ungroup()
Nymphoides.pa$ApproxID <- NULL

#some locations show first absence and later presence - combine to show presence only
Nymphoides.pa <-
  Nymphoides.pa %>% group_by(N) %>% mutate(nymtot = sum(Nymphoides_pa))
Nymphoides.pa$nymtot <-
  sign(Nymphoides.pa$nymtot) #turn count into PA
Nymphoides.pa$Nymphoides_pa <- NULL
Nymphoides.pa <-
  Nymphoides.pa %>% rename("Nymphoides_pa" = "nymtot") %>%
  dplyr::select(N, E, Nymphoides_pa) %>% distinct() %>% ungroup()

#From artdatabanken ####
artportalen <-
  read_excel("Data/Artdata_Nymphoides.xlsx") 
artportalen2 <-
  artportalen %>% dplyr::select(Ost, Nord, Slutdatum) %>% mutate(Nymphoides_pa = 1) %>%
  rename(N = Nord, E = Ost) %>% drop_na() %>% distinct()

artportalen2$Slutdatum <- str_sub(artportalen2$Slutdatum,1,4) %>% as.numeric()
artportalen2 <- filter(artportalen2, Slutdatum >= 2000)
artportalen2$Slutdatum <- NULL

artportcoord <- artportalen2 %>% select(E, N) %>% as.matrix()

# Create sf object with RT90 projection
df_sf <- st_as_sf(artportalen2, coords = c("E", "N"), crs = 3021)  # 3021 is the EPSG code for RT90

# Transform the coordinates to SWEREF99 (EPSG: 3006)
df_sweref99_sf <- st_transform(df_sf, crs = 3006)

# Extract transformed coordinates
df_sweref99_coords <- st_coordinates(df_sweref99_sf)
df_sweref99_coords
df_sweref99_coords <- as.data.frame(df_sweref99_coords)

#copy new coords to artdatabank data
artportalen2$E <- df_sweref99_coords$X
artportalen2$N <- df_sweref99_coords$Y

#Combine SLU and Artdatabank data ####

combi_data <- full_join(Nymphoides.pa, artportalen2) %>% distinct()
ggplot(combi_data, aes(E, N)) +
  theme_bw() +
  geom_point() +
  geom_point(aes(colour = factor(Nymphoides_pa)))

#5. MATCH OBSERVATIONS TO POLYGONS IN SHAPEFILE##################################
shape.file <- read_sf("Data/Lakes_info_SWE_proportions_CLC.shp")
st_crs(shape.file)#check is SWEREF99

#Combined data coordinates

combi_pnts_sf3 <-
  st_as_sf(combi_data,
           coords = c('E', 'N'),
           crs = st_crs(3006))
st_crs(combi_pnts_sf3)
plot(combi_pnts_sf3)

match_combi <- st_join(combi_pnts_sf3,
                         shape.file,
                         join = st_within,
                         left = F)
plot(match_combi)

#Extract geometry to dataframe columns###
match_combi <- match_combi %>%
  mutate(E = unlist(map(match_combi$geometry, 1)),
         N = unlist(map(match_combi$geometry, 2)))

match_combi$geometry <- NULL
match_combi <-
  dplyr::select(match_combi, OBJECTID, E, N, Nymphoides_pa) %>% distinct()
match_combi$OBJECTID <- as.factor(match_combi$OBJECTID)

#many observations are repeats from the same polygon, for example well-visited recreation
#areas near a major population centre (so lots of reports to Artdatabanken).
#Combine to PA per polygon...
combined_Nymphoides <-
  match_combi %>% group_by(OBJECTID) %>% mutate(nymtot = sum(Nymphoides_pa))
combined_Nymphoides$nymtot <-
  sign(combined_Nymphoides$nymtot) #turn count into PA
combined_Nymphoides$Nymphoides_pa <- NULL
combined_Nymphoides <-
  combined_Nymphoides %>% rename("Nymphoides_pa" = "nymtot") %>%
  dplyr::select(OBJECTID, Nymphoides_pa) %>% distinct() %>% ungroup()

#Add the PA data matched to polygon to the transformed data from the shapefiles
#to get only prescence and confirmed abscence
Exp_vars3 <-
  left_join(combined_Nymphoides, Exp_vars, by = "OBJECTID") %>%
  distinct() %>%
  drop_na() %>%
  ungroup()

#all observations
swe_shape <- st_read("Data/LAomraden_2018_region.shp", as_tibble = T)

ggplot() +
  geom_sf(data = swe_shape) +
  coord_sf(label_axes = "----")+
  theme_bw()+
  scale_color_manual(name = 'Nymphoides\npresent', values = c("0" = "darkblue",
                                                              "1"="red")) +
  geom_point(data=combi_data, aes(E, N, colour = factor(Nymphoides_pa)))

#one record per lake polygon
ggplot() +
  geom_sf(data = swe_shape) +
  coord_sf(label_axes = "----")+
  theme_bw()+
  scale_color_manual(name = 'Nymphoides\npresent', values = c("0" = "darkblue",
                                                              "1"="red")) +
  geom_point(data=Exp_vars3, aes(E, N, colour = factor(Nymphoides_pa)))

#all data (including NAs for presence/absence)
Exp_vars_projdat <-
  left_join(Exp_vars, combined_Nymphoides, by = "OBJECTID") %>% distinct()

#Calcium as limiting factor
 calc <- filter(predicted_values, Ca > 0.755) #from Smits 1992

ggplot() +
  geom_sf(data = swe_shape) +
  coord_sf(label_axes = "----")+
  theme_bw()+
  geom_point(data=calc, aes(E, N))

# 6     MODELLING      ########## 
####  LOADING DATA

#Option to keep all data including NAs for Nymphoides presence
#eg. if we want to generate pseudo-absences
data_all_XY <- Exp_vars_projdat

#Do this to keep only presence and absence
data_XY <- Exp_vars3
data_XY <- rename(data_XY, X = E, Y = N)
nadata_XY = na.omit(data_XY)
nadata_XY <- nadata_XY %>% distinct()
nrow(nadata_XY) 


#### FORMATING DATA RASTERS VERSION #### 
#Specify biomod data ###
speName='Nymphoides_pa'
espPA = nadata_XY[,speName]
espXY = dplyr::select(nadata_XY, X, Y)
espEXPL = envdat.sr
#Combine to biomod data object
BiomodData = BIOMOD_FormatingData(resp.var = espPA, 
                                  expl.var = espEXPL,
                                  resp.xy = espXY, 
                                  resp.name = speName)

(BiomodData)

plot(BiomodData)   #green presence red absence


####  MODELING ####

#start cluster
#uses doParallel
cl <- makeCluster(8) #use 8 cores
doParallel::registerDoParallel(cl)

# stopCluster(cl)

 
BiomodModelOut_all = BIOMOD_Modeling (BiomodData, 
                                   # models = c("ANN", "CTA", "FDA", "GAM", "GLM", "MARS", "MAXNET", "RF",
                                    #            "SRE", "XGBOOST"),
                                                 models = c("ANN", "CTA", "FDA", "GAM", "GLM", "MARS", "MAXNET"
                                                 ),#SRE and XGBOOST don't work well with categorical variables like temp limit
                                                # "RF") #"RF" badly overfitting
                                      CV.nb.rep = 20,   
                                      OPT.strategy = 'bigboss', 
                                      CV.perc = 0.8,                   #80 - 20 test train split 
                                      prevalence = 0.5,                       #absence has same weight as presences 
                                      var.import = 50,                        #nb of permutations for variable importance
                                      metric.eval = c('TSS', 'ROC'),     
                                      scale.models = T,                  
                                      CV.do.full.models = TRUE,                 
                                      nb.cpu = 6,
                                      modeling.id = paste(speName, "FirstModeling", sep=""))



####  EVALUATION ####

BiomodModelEval_all = get_evaluations(BiomodModelOut_all)

#plot all TSS and ROC by model
bm_PlotEvalMean(BiomodModelOut_all, metric.eval = c("ROC","TSS"),
                group.by = "algo",
                do.plot = TRUE)

results <- bm_PlotEvalMean(bm.out = BiomodModelOut_all)
mean(results$tab$mean1)
mean(results$tab$mean2)


####  VARIABLE IMPORTANCE  ########################### 

##check variable importance 
VarImportance <- get_variables_importance(BiomodModelOut_all)
VarImportance

bm_PlotVarImpBoxplot(bm.out = BiomodModelOut_all, group.by = c('expl.var', 'expl.var', 'PA'))


get_variables_importance(BiomodModelOut_all) %>% 
  group_by(expl.var) %>% 
  summarise(var.imp = median(var.imp)) %>% 
  arrange(desc(var.imp))

####  ENSEMBLE MODELING #####################################

BiomodEM = BIOMOD_EnsembleModeling(bm.mod = BiomodModelOut_all,
                                   models.chosen = 'all',
                                   em.by = 'all',
                                   metric.select = c('TSS','ROC'),
                                   metric.select.thresh = c(0.6, 0.8),
                                   metric.eval = c('TSS', 'ROC'),
                                   var.import = 50,
                                   EMci.alpha = 0.05,
                                   EMwmean.decay	= 'proportional',
                                   seed.val = 123,
                                   nb.cpu = 6,
                                   
)                                   


BiomodEM

# Get evaluation scores & variables importance
get_evaluations(BiomodEM)
get_variables_importance(BiomodEM)

##check variable importance ensemble
VarImportanceEM <- get_variables_importance(BiomodEM)

(imp.tab <- get_variables_importance(BiomodEM) %>% 
  group_by(expl.var) %>% 
  summarise(var.imp = mean(var.imp)) %>% 
  arrange(desc(var.imp)))

sum(imp.tab$var.imp)

# Represent evaluation scores & variables importance

bm_PlotVarImpBoxplot(bm.out = BiomodEM, group.by = c('expl.var', 'algo', 'merged.by.PA'))

# Represent response curves
mods <- get_built_models(BiomodEM)
bm_PlotResponseCurves(bm.out = BiomodEM, 
                      models.chosen = mods,
                      fixed.var = 'median')


####  PROJECTION (single models) #################################
speName='Nymphoides_pa'

BiomodProjec <- BIOMOD_Projection(bm.mod = BiomodModelOut_all,
                                  proj.name = 'current',
                                  new.env = envdat.sr,
                                  new.env.xy = espXY,
                                  #models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = T)

#plot(BiomodProjec, str.grep = 'RF')

ClampingMask <- rast('Nymphoides.pa/proj_current/proj_current_ClampingMask.tif')
#plot(ClampingMask)


####  ENSEMBLE FORECASTING #######################

# Project ensemble models (from single projections)
BiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = BiomodEM, 
                                           bm.proj = BiomodProjec,
                                           models.chosen = 'all',
                                           metric.binary = c('all'), 
                                           metric.filter = c('all'),
                                           output.format = '.grd') 
#plot(BiomodEMProj)

# Project ensemble models (building single projections)
BiomodEMProj2 <- BIOMOD_EnsembleForecasting(bm.em = BiomodEM,
                                            proj.name = 'CurrentEM',
                                            new.env = envdat.sr,
                                            models.chosen = 'all',
                                            metric.binary = 'all',
                                            metric.filter = 'all',
                                            output.format = '.grd',
                                            nb.cpu = 4)

# plot(BiomodEMProj2)
# plot(BiomodEMProj2, str.grep = 'Nymphoides.pa_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData')
# plot(BiomodEMProj2, str.grep = 'Nymphoides.pa_EMwmeanByROC_mergedAlgo_mergedRun_mergedData')

#### Alternative map ###### 

mygrd_ensemble <- raster('Nymphoides.pa/proj_CurrentEM/proj_CurrentEM_Nymphoides.pa_ensemble.grd')
mygrd_rescaled <- mygrd_ensemble/1000

plot(mygrd_rescaled, useRaster=T, interpolate = T,col = hcl.colors(20, "Viridis", rev = TRUE))

#get borders
swe_shape2 <- st_read("Data//hd324xr2654.shp", as_tibble = T)
swe_shape2 <-
  st_transform(swe_shape2, crs = 3006) #SWEREF99

#Using ggplot #
#convert the raster to points for plotting
map.p <- rasterToPoints(mygrd_rescaled)

#Make the points a dataframe for ggplot
df <- data.frame(map.p)

#Make appropriate column headings
colnames(df) <- c('E', 'N', 'layer')

#without borders
ggplot(data=df, aes(y=N, x=E)) +
  geom_raster(aes(fill=layer), interpolate = T) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  coord_equal()+
  scale_fill_viridis_c(direction = -1, name = "Probability")

#with borders
ggplot() +
  geom_sf(data = swe_shape2, fill= "white") +
  coord_sf(label_axes = "----")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_viridis_c(direction = -1, name = "Probability")+
  geom_raster(data=df, aes(y=N, x=E, fill = layer), interpolate = T)

#with borders, tile function
# ggplot() +
#   geom_sf(data = swe_shape2, fill= "white") +
#   coord_sf(label_axes = "----")+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank())+
#   scale_fill_viridis_c(direction = -1, name = "Probability") +
#   geom_tile(data = df, aes(y=N, x=E, fill = layer))
# geom_raster(data=df, aes(y=N, x=E, fill = layer), interpolate = T)


