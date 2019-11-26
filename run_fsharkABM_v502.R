############################################################
# Shark MPA model version 5.02 by Ross Dwyer, 30 Sep 2019 #
############################################################

# in contrast to v1 this version can calculate exposure also based on movement probabilities
# which is initiated by setting probs.per.ind to 1
# 
# 4.01 vs. 4.02: - fixed location of individuals at BRUVs according to movement profile
#                - included option to decrease resolution and 

# 4.04 vs 4.05 - Added capacity to define write folder name and location
#              - Bring in Vinay's day, week, months, year, month, 

# 4.05 vs 5.01- fixed bug where rewriting species dispersal profiles on grey reef sharks 7/05/2019
#             - added histogram plot of different species dispersal capacity

# 5.01 vs 5.02- edits following R1 @ Current Biology
#             - added site specific data for BRS and GRS


rm(list = ls()) # delete parameters in the workspace

versionfolder <- "shark_mpa_model"

library(stringr)

# load data
paste('Need to specifiy datadir - wherever you place the folder')
#datadir <-  paste0('C:/Users/uqnkruec_local/Dropbox/Papers/Shark MPAs/', version,'/')
#datadir <- #paste0('C:/Users/weizenkeim/Dropbox/Papers/Shark MPAs/', version,'/')
#datadir <- paste0('E:/OneDrive/Australia/SharkRay MPAs/Agent based models/', versionfolder,'/')
#datadir <- paste0('C:/Users/uqnkruec_local/Dropbox/Papers/Shark MPAs/', versionfolder ,'/')
#setwd(datadir)
datadir <- ''

#maxndata <- read.csv(paste0(datadir, "maxndata.csv")) # maxndata2.csv - added Caribbean Species and Removed Location Classifications
#maxddata <- read.csv(paste0(datadir, "maxddata.csv"))

# Read in Vinay's new dispersal dataset which has different timeframes
Dispersal_Timescales <- readRDS("2018-11-01_Dispersal_Timescales.RDS") # 
TagMData <- read.csv("2019-09-18_TagMetadata_installations.csv") # 

maxddata.d <- data.frame(Dispersal_Timescales[[1]]) # Load the daily data
maxddata.w <- data.frame(Dispersal_Timescales[[2]]) # Load the weekly data
maxddata.m <- data.frame(Dispersal_Timescales[[3]]) # Load the monthly data
maxddata.y <- data.frame(Dispersal_Timescales[[4]]) # Load the yearly data

# Standardise the Date field column
names(maxddata.d)[1] <- "Date"
names(maxddata.w)[1] <- "Date"
names(maxddata.m)[1] <- "Date"
names(maxddata.y)[1] <- "Date"
maxddata.d[,1] <- as.character(maxddata.d[,1])
maxddata.w[,1] <- as.character(maxddata.w[,1])
maxddata.m[,1] <- as.character(maxddata.m[,1])
maxddata.y[,1] <- as.character(maxddata.y[,1])
# Add column detailing time period
maxddata.d$Ftime <- 'day'
maxddata.w$Ftime <- 'week'
maxddata.m$Ftime <- 'month'
maxddata.y$Ftime <- 'year'

# Combine data into a dataframe for plotting
maxddata.df <- rbind(maxddata.d,
                     maxddata.w,
                     maxddata.m,
                     maxddata.y)
# Reorder time period factor to day - year
maxddata.df$Ftime <- factor(maxddata.df$Ftime)
maxddata.df$Ftime <- factor(maxddata.df$Ftime ,c("day","month","week","year"))

# Reorder species factor to alphabetical latin names
maxddata.df$species_name <- maxddata.df$species
maxddata.df <- maxddata.df[!is.na(maxddata.df$common_name),] # for some odd reason there are NAs in the dispersal file
maxddata.df$common_name <- factor(maxddata.df$common_name)
maxddata.df$common_name <- factor(maxddata.df$common_name ,c("Whitetip Reef Shark",
                                                             "Blacktip Reef Shark",
                                                             "Grey Reef Shark",
                                                             "Caribbean Reef Shark",
                                                             "Nurse Shark"))

# Plot species histograms of dispersal distances under different temporal periods
#spec.dist.hist <- ggplot(maxddata.df, 
#                         aes(max.dispersal_km, fill = common_name)) + 
#  geom_histogram(binwidth = 5) + 
#  xlab("Maximum dispersal distance (km)")+
#  ylab("Number of observations")+
#  facet_wrap(common_name~Ftime,nrow=5, scales = "free") +
#  theme_classic(base_size = 13) + 
#  labs(linetype = "F")+
#  theme(legend.position="none",
#        axis.text.x = element_text(color="black"),
#        axis.text.y = element_text(color="black"),
#        strip.background = element_blank(),
        #strip.placement = "inside",
#        strip.text=element_text(vjust=-1),
#        panel.background = element_blank()) 

#spec.dist.hist
# Save the plot
#ggsave(paste0(datadir,"/Images/Species dispersal histograms.png"), spec.dist.hist)


####

### This is where you set the temporal dataset and folder for the modelling

# Set the temporal period 
# maxddata <- maxddata.d
# resfoldername <- 'Ross/results/day' # Which folder it saves to

maxddata <- maxddata.w
resfoldername <- 'Ross/results/week' # Which folder it saves to

#maxddata <- maxddata.m
#resfoldername <- 'Ross/results/month' # Which folder it saves to

#maxddata <- maxddata.y
#resfoldername <- 'Ross/results/year' # Which folder it saves to
  
library(dplyr)
maxddata <- left_join(maxddata,TagMData)

maxddata$species_name <- maxddata$species
maxddata <- maxddata[!is.na(maxddata$common_name),] # for some odd reason there are NAs in the dispersal file
maxddata$species <- str_extract(maxddata$species, '[^ ]+$')

#Load the species specific maxn data
maxndata_WRS <- read.csv(paste0(datadir, "maxndata_WRS.csv"))
maxndata_BRS <- read.csv(paste0(datadir, "maxndata_BRS.csv"))
maxndata_GRS <- read.csv(paste0(datadir, "maxndata_GRS.csv"))
maxndata_GRS <- maxndata_GRS[!is.na(maxndata_GRS$location_code),] # for some odd reason there are NAs in the count file for Grey reefs
maxndata_CRS <- read.csv(paste0(datadir, "maxndata_CRS.csv"))
maxndata_NS <- read.csv(paste0(datadir, "maxndata_NS.csv"))

# get the data into list format to enable  all 5 to run in parallel 
# 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
maxndata2 <- list(0)
maxndata2[[1]] <- maxndata_WRS
maxndata2[[2]] <- maxndata_BRS
maxndata2[[3]] <- maxndata_GRS
maxndata2[[6]] <- maxndata_CRS
maxndata2[[7]] <- maxndata_NS

#setwd(paste0(datadir, 'Nils/'))

#source('Ross/fsharkABM_v404.R')
source('Ross/fsharkABM_v502.R')

# One species at a time
meanmaxextent <- 1
probsperind <- 1
replicates <- 1000
ddata <- maxddata



#### Run the model for single species and all sites

# Whitetip
fsharkABM_v502(bruvdat = maxndata2[[1]],
               movedat = ddata,
               abundancecat=2,
               speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
               nreplicates = replicates,
               probs.per.ind = probsperind,
               mean.max.extent = meanmaxextent,
               resolution = 1,
               resfolder=resfoldername)

# Blacktip
fsharkABM_v502(bruvdat = maxndata2[[2]],
               movedat = ddata,
               abundancecat=2,
               speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
               nreplicates = replicates,
               probs.per.ind = probsperind,
               mean.max.extent = meanmaxextent,
               resolution = 1,
               resfolder=resfoldername)

# Grey
fsharkABM_v502(bruvdat = maxndata2[[3]],
               movedat = ddata,
               abundancecat=2,
               speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
               nreplicates = replicates,
               probs.per.ind = probsperind,
               mean.max.extent = meanmaxextent,
               resolution = 1,
               resfolder=resfoldername)

# Caribbean
fsharkABM_v502(bruvdat = maxndata2[[6]],
               movedat = ddata,
               abundancecat=2,
               speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
               nreplicates = replicates,
               probs.per.ind = probsperind,
               mean.max.extent = meanmaxextent,
               resolution = 1,
               resfolder=resfoldername)

# Nurse
fsharkABM_v502(bruvdat = maxndata2[[7]],
               movedat = ddata,
               abundancecat=2,
               speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
               nreplicates = replicates,
               probs.per.ind = probsperind,
               mean.max.extent = meanmaxextent,
               resolution = 1,
               resfolder=resfoldername)
# 
# 
# library(parallel)
# 
# # Use the detectCores() function to find the number of cores in system
# no_cores <- detectCores()
# 
# # Setup cluster
# clust <- makeCluster(no_cores) #This line will take 
# 
# #To ensure function and data are sent to the workers, export the names directly
# clusterExport(clust, list("fsharkABM","maxndata2","maxddata","datadir"))
# 
# #The parallel version of lapply() is parLapply() and needs an additional cluster argument.
# parLapply(clust,c(2), function(x) fsharkABM(bruvdat = maxndata2[[x]],#1,2,3,6,7
#                                                   movedat = maxddata,
#                                                   abundancecat=2,
#                                                   speciesToSimulate = x,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
#                                                   rBRUVcatchment = 400, # radius of BRUV plume catchment area: 200,400,600 
#                                                   nreplicates = 1000,
#                                                   probs.per.ind = 0,
#                                                   years = 25, # lifetime over which mortality is assessed
#                                                   mortRate = c("pMortsInst"), # "pMorts","pMortsInst"
#                                                   save.results = 1, # save all results
#                                                   plot.example = 0, # plot example mpa scenario to validate modelling procedure
#                                                   plot.results = 1, # plot key mortality outcomes and surviving individuals
#                                                   save.plots = 1, # save all plots to png file - not showing then
#                                                   maxd.per.ind = 1, # calculate max travel distance per individual at any day
#                                                   probs.per.ind = 1,
#                                                   nreplicates = 10)) #This output is a vector)
# 
# 

#########################################

##Now run the dispersal data for each installation seperately

unique(ddata$installation) #Brazil   Belize   Ningaloo Rowley   Florida  Heron    FNQ      Scott  

ddata_Nin <- ddata %>%
  filter(installation=='Ningaloo')
ddata_Row <- ddata %>%
  filter(installation=='Rowley')
ddata_Her <- ddata %>%
  filter(installation=='Heron')
ddata_Fnq <- ddata %>%
  filter(installation=='FNQ')
ddata_Sco <- ddata %>%
  filter(installation=='Scott')

ddata_site <- ddata %>%
  #filter(common_name %in% c("Blacktip Reef Shark","Grey Reef Shark")) %>%
  group_by(tag_id,common_name,installation) %>%
  summarise(array_area = mean(array_area_m2),
            length_mm=mean(length_mm),
            nrecords = n(),
            meandisp = mean(max.dispersal_km),
            sddisp = sd(max.dispersal_km),
            mediandisp = median(max.dispersal_km))

ddata_site %>%
  group_by(common_name,installation) %>%
  summarise(array_area = mean(array_area),
            length_mm=mean(length_mm),
            nrecords = n(),
            meand = mean(meandisp),
            sddisp = sd(meandisp),
            mediandisp = mean(mediandisp),
            nsharks = n_distinct(tag_id))

###
## First Grey Reef sharks

# Grey @ Ningaloo
fsharkABM_v502_sites(bruvdat = maxndata2[[3]],
                     movedat = ddata_Nin,
                     abundancecat=2,
                     speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
                     nreplicates = 1000,
                     probs.per.ind = probsperind,
                     mean.max.extent = meanmaxextent,
                     resolution = 1,
                     resfolder=paste0(resfoldername,"/Ningaloo"))


# Grey @ Rowley
fsharkABM_v502_sites(bruvdat = maxndata2[[3]],
                     movedat = ddata_Row,
                     abundancecat=2,
                     speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
                     nreplicates = 1000,
                     probs.per.ind = probsperind,
                     mean.max.extent = meanmaxextent,
                     resolution = 1,
                     resfolder=paste0(resfoldername,"/Rowley"))

# Grey @ Heron
fsharkABM_v502_sites(bruvdat = maxndata2[[3]],
                     movedat = ddata_Her,
                     abundancecat=2,
                     speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
                     nreplicates = 1000,
                     probs.per.ind = probsperind,
                     mean.max.extent = meanmaxextent,
                     resolution = 1,
                     resfolder=paste0(resfoldername,"/Heron"))

# Grey @ FNQ
fsharkABM_v502_sites(bruvdat = maxndata2[[3]],
                     movedat = ddata_Fnq,
                     abundancecat=2,
                     speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
                     nreplicates = 1000,
                     probs.per.ind = probsperind,
                     mean.max.extent = meanmaxextent,
                     resolution = 1,
                     resfolder=paste0(resfoldername,"/FNQ"))

# Grey @ Scott
fsharkABM_v502_sites(bruvdat = maxndata2[[3]],
                     movedat = ddata_Sco,
                     abundancecat=2,
                     speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
                     nreplicates = 1000,
                     probs.per.ind = probsperind,
                     mean.max.extent = meanmaxextent,
                     resolution = 1,
                     resfolder=paste0(resfoldername,"/Scott"))

###
## Second Blacktip Reef sharks
# Blacktip @ Ningaloo
fsharkABM_v502_sites(bruvdat = maxndata2[[2]],
                     movedat = ddata_Nin,
                     abundancecat=2,
                     speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
                     nreplicates = 1000,
                     probs.per.ind = probsperind,
                     mean.max.extent = meanmaxextent,
                     resolution = 1,
                     resfolder=paste0(resfoldername,"/Ningaloo"))


# Blacktip @ Rowley
fsharkABM_v502_sites(bruvdat = maxndata2[[2]],
                     movedat = ddata_Row,
                     abundancecat=2,
                     speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
                     nreplicates = 1000,
                     probs.per.ind = probsperind,
                     mean.max.extent = meanmaxextent,
                     resolution = 1,
                     resfolder=paste0(resfoldername,"/Rowley"))

# Blacktip @ Heron
fsharkABM_v502_sites(bruvdat = maxndata2[[2]],
                     movedat = ddata_Her,
                     abundancecat=2,
                     speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
                     nreplicates = 1000,
                     probs.per.ind = probsperind,
                     mean.max.extent = meanmaxextent,
                     resolution = 1,
                     resfolder=paste0(resfoldername,"/Heron"))

# Blacktip @ FNQ
fsharkABM_v502_sites(bruvdat = maxndata2[[2]],
                     movedat = ddata_Fnq,
                     abundancecat=2,
                     speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
                     nreplicates = 1000,
                     probs.per.ind = probsperind,
                     mean.max.extent = meanmaxextent,
                     resolution = 1,
                     resfolder=paste0(resfoldername,"/FNQ"))

# Blacktip @ Scott
fsharkABM_v502_sites(bruvdat = maxndata2[[2]],
                     movedat = ddata_Sco,
                     abundancecat=2,
                     speciesToSimulate = 1,  #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse shark 
                     nreplicates = 1000,
                     probs.per.ind = probsperind,
                     mean.max.extent = meanmaxextent,
                     resolution = 1,
                     resfolder=paste0(resfoldername,"/Scott"))

