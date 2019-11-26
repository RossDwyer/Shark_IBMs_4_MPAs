############################################################
# Shark MPA model version 5.02 by Ross Dwyer @ Nils Krueck, 30 Sep 2019 #
############################################################

# in contrast to v1 this version can calculate exposure also based on movement probabilities
# which is initiated by setting probs.per.ind to 1
# 
# 4.01 vs. 4.02: - fixed location of individuals at BRUVs according to movement profile
#                - included option to decrease resolution and speed up computation time 
#                - double checking cumulative probability 

#4.05 vs 4.04 - added results folder definitiion foldr
#5.01 vs 4.05 - fix bug where species dispersal probabilities not erased prevvious species dispersal probabilities
#5.02 vs 5.01 - add site specific dispersal data


fsharkABM_v502 <- function(bruvdat,
                      movedat,
                      dispdat_limit=1,
                      abundancecat=2, # How many abundance categories are there? 1,2,3?
                      speciesToSimulate, # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse 
                      rBRUVcatchment=400/resolution, # diameter of BRUV plume catchment area: 200,400,600 
                      years=25, # lifetime over which mortality is assessed
                      save.results = 1, # save all results 
                      plot.example = 0, # plot example mpa scenario to validate modelling procedure
                      plot.results = 1, # plot key mortality outcomes and surviving individuals
                      save.plots = 1, # save all plots to png file - not showing then
                      maxd.per.ind = 1, # calculate max travel distance per individual at any day 
                      probs.per.ind = 1, # The activity Range: 0 = uniform (use max dispersal distance) or 1 = ^ shaped (uses daily distances around fixed point)  
                      nreplicates = 1000,
                      mean.max.extent = 1, # use mean of max daily dispersal instead of max of max
                      resolution = 1,
                      resfolder = 'results') # resolution of modelling environment in m
{ # calculate movement probabilities per individual)
  
  # resolution = 1
  # bruvdat = maxndata_NS
  # movedat = maxddata
  # dispdat_limit=7
  # abundancecat = 2 #1,2,3 (same, highlow, highmediumlow)
  # speciesToSimulate = 7   #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse
  # rBRUVcatchment = 400/resolution # radius of BRUV plume catchment area: 200,400,600
  # years = 25 # lifetime over which mortality is assessed
  # mortRate = c("pMortsInst") # "pMorts","pMortsInst"
  # save.results = 1 # save all results
  # plot.example = 0 # plot example mpa scenario to validate modelling procedure
  # plot.results = 1 # plot key mortality outcomes and surviving individuals
  # save.plots = 1 # save all plots to png file - not showing then
  # maxd.per.ind = 1 # calculate max travel distance per individual at any day
  # probs.per.ind = 0
  # nreplicates = 10 # number of replicate simulations per species, location and MPA size
  # #
  # Modelling parameters
  
  # Removes individuals with fewer than X days of mobements
  dfmaxd <- as.data.frame(table(movedat$tag_id))
  exclude_ids <- dfmaxd$Var1[which(dfmaxd$Freq<dispdat_limit)]
  '%!in%' <- function(x,y)!('%in%'(x,y))
  movedat1 <- movedat[movedat$tag_id %!in% exclude_ids, ]
  
  MPAsizes = c(0,seq(500,2000,500),seq(3000,5000,1000),seq(10000,20000,5000),c(30000,40000,50000,100000))/resolution
  
  pMortsInst = seq(0,.3,.05) # probability of mortality when exposed to fishing
  
  #if(mortRate=="pMorts")
  #  pMorts = 1-exp(-pMortsInst)
  #if(mortRate=="pMortsInst")
  #  pMortsInst = -log(1-pMorts) 
  pMorts = 1-exp(-pMortsInst)
  pMortsInst = -log(1-pMorts) # convert discrete to instantaneous mortality rates (Q for Nils - do I have to do this step?)
  
  # Select abundance category
  if(abundancecat==1) # No category
  {
    bruvdat$location_code <- 2 #Standardises locations
    abundance = c(2) #   2 average abundance
  }
  if(abundancecat==2){ # 2 categories
    abundance = unique(bruvdat$location_code)
    bruvdat$location_code # 1 = low abundance, 2 = average abundance
  }
  if(abundancecat==3){ # 3 categories
    abundance = unique(bruvdat$location_code_high) #  1 low abundance, 2 average abundance, 3 high abundance 
    location_code = bruvdat$location_code_high # Rewrite column with 2 levels as 3 levels
  }
  
  # explore nmax data
  lat.spec.names <- unique(bruvdat$species)
  comm.spec.names <- unique(bruvdat$common_name)
  #speciesToAnalyse = c(as.character(lat.spec.names[speciesToSimulate]),"") # has to include empty cells melanopterus","") #c("obesus","melanopterus","amblyrhynchos","perezi") 
  speciesToAnalyse = c(as.character(lat.spec.names[1])) # has to include empty cells melanopterus","") #c("obesus","melanopterus","amblyrhynchos","perezi") 
  
  region.names <- unique(bruvdat$region)
  location.names <- unique(bruvdat$location_name)
  location.code.names <- c("low", "average", "high") 
  nspecs <- length(speciesToAnalyse) #updated when turning into function
  nlocations <- length(abundance) # coded as 1 (low), 2(med), 3(high) in vector location_code for low, medium, high abundance 
  
  # allocate output and working matrices
  sprobslist <- vector("list", nspecs)
  probslist <- vector("list")
  probsdflist <- vector("list")
  
  ioutputHeader <- c("species", "abundance", "replicate", "mpa_size", "home_range", "protection")
  ioutput = matrix(nrow=0,ncol=length(ioutputHeader)) # output table for individual data 
  statstableHeader <- c("species", "abundance", "mean_maxn", "sd_maxn", "max_hr", "mpa_size", "mean_nsharks", "mean_nfull", "mean_npart","mean_prot", "sd_prot", "median_prot")
  statstable = matrix(nrow=nspecs*nlocations*length(MPAsizes),ncol=length(statstableHeader)) # allocate output table for individual data 
  statsmatsHeader = c("species", "abundance", "mpa_size",as.character(pMortsInst))
  
  meanmeanp_survival = matrix(data=0,nrow=nspecs*nlocations*length(MPAsizes),ncol=length(statsmatsHeader))
  
  minminp_survival = meanmeanp_survival
  maxmaxp_survival = meanmeanp_survival
  meanminp_survival = meanmeanp_survival
  meanmaxp_survival = meanmeanp_survival
  p05meanp_survival = meanmeanp_survival
  p95meanp_survival = meanmeanp_survival
  minmeanp_survival = meanmeanp_survival
  maxmeanp_survival = meanmeanp_survival
  meansdp_survival = meanmeanp_survival
  meanp05p_survival = meanmeanp_survival
  meanp95p_survival = meanmeanp_survival
  
  meansumn_survival = meanmeanp_survival
  mediansumn_survival = meanmeanp_survival
  sdsumn_survival = meanmeanp_survival
  p05sumn_survival = meanmeanp_survival
  p95sumn_survival = meanmeanp_survival
  
  meanpropn_mortality = meanmeanp_survival
  medianpropn_mortality = meanmeanp_survival
  sdpropn_mortality = meanmeanp_survival
  p05propn_mortality = meanmeanp_survival
  p95propn_mortality = meanmeanp_survival
  
  mean_mortality_offset = matrix(nrow=nspecs*nlocations*length(MPAsizes),ncol=length(statsmatsHeader))
  sd_mortality_offset = mean_mortality_offset
  mean_mortalityInst_offset = sd_mortality_offset
  sd_mortalityInst_offset = sd_mortality_offset
  
  # create folders to save results
  if(!dir.exists(paste0(datadir, resfolder))){dir.create(paste0(datadir, resfolder))}
  if(!dir.exists(paste0(datadir,resfolder, "/rBRUVcatchment_",rBRUVcatchment))){dir.create(paste0(datadir,resfolder, "/rBRUVcatchment_",rBRUVcatchment))}
  
  
  ##################
  # START OF MODEL #
  ##################
  
  ## Only have to run this once at the beginning
  
  # count for indexing purposes
  count = 0
  
  # species loop
  #for (s in 1:nspecs){  # Removed species loop to run as part of a function
  s <- 1 # This has been standardised to a 1 as a result of turning the loop into a function which runs independently on each species
  
  # create and assign folders to save results
  if(!dir.exists(paste0(datadir,resfolder, "/rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))){
    dir.create(paste0(datadir,resfolder, "/rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))}
  savedir = paste0(datadir,resfolder, "/rBRUVcatchment_", rBRUVcatchment,"/", speciesToAnalyse[s],"/") 
  
  starttime <- Sys.time() # save time
  
  # select data
  spec <- speciesToSimulate[s] 
  speclocs <- which(bruvdat$species == speciesToAnalyse[s] | bruvdat$species == speciesToAnalyse[nspecs+1])
  #specdata <- bruvdat[speclocs,] # step set up to deal with old rel abundance data for all species
  specdata <- bruvdat # new step set up to deal with new species-specific abundance data for all species
  specname <- comm.spec.names[1]
  speclatname <- speciesToAnalyse[s]
  
  if (maxd.per.ind == 1){
    ids = unique(movedat1$tag_id[movedat1$species == speciesToAnalyse[s]]) # get shark ids
    dists = matrix(nrow=length(ids),ncol=1)
    for (id in 1:length(ids)) {dists[id] = max(movedat1$max.dispersal_km[movedat1$tag_id == ids[id]]) * 1000}
  }else{
    dists = movedat1$max.dispersal_km[movedat1$species == speciesToAnalyse[s]] * 1000 / resolution
  } # travel distances in m
  #hist(dists,breaks=length(dists)/4)
  #length(dists)
  

  if (maxd.per.ind == 1 & probs.per.ind == 1){
    # generate probability of occurrence along the distance spectrum
    for (id in 1:length(ids)) {
      alldists = movedat1$max.dispersal_km[movedat1$tag_id == ids[id]] * 1000 / resolution
      # hist(alldists, breaks = 20, xlab = "Travel distance (m)", ylab = "Frequency")
      kernel = as.matrix(rep(1,1,max(alldists)))
      centreloc = ceiling(max(alldists)/2)
      for (dshark in 1:length(unique(alldists))) { # for each observed travel distance (RD WHY dshark???)
        addval = 1 # set initial number of observed vals to 1
        maxlocs = which(alldists == max(alldists)) # identify all distance locs
        startloc = round(centreloc - alldists[maxlocs[1]]/2) # define kernel start loc
        endloc = floor(startloc + alldists[maxlocs[1]]) # define kernel end loc
        if (length(maxlocs) > 1){addval = length(maxlocs)} # if observed more than once update value
        kernel[startloc:endloc] = kernel[startloc:endloc] + addval # add observations along kernel spectrum
        alldists = alldists[alldists != alldists[maxlocs[1]]] # delete observed distance from list
      }
      probslist[[id]] = kernel/sum(kernel,na.rm=TRUE) # get cumulative probability spectrum by dividing by sum
      
      # generate movement probability dataframes for plotting
      probsdf <- data.frame(ID=as.factor(ids[id]),
                            DISTANCE=1:length(probslist[[id]]),
                            P=probslist[[id]])
      probsdflist[[id]] <- probsdf
      
    } # end of individual loop
    sprobslist[[s]] = probslist
    
    probslist <- NULL #RD Reset these values to fix bug in v404
    
  } # end of probs calculations
  
  # location loop (# coded as 1, 2, 3 in vector location_code for low, medium, high abundance)
  for (lc in 1:nlocations){
    
    # select data (low to high abundance)
    #slocs = sum(which(specdata$location_code == lc | specdata$location_code == nlocations+1))
    
    
    # THis is where to add flexibility to set the different High Medium Low categories
    
    loc = abundance[lc]
    locname = location.code.names[loc]
    #locdata <- specdata[which(specdata$location_code == loc | specdata$location_code == 4),] # | specdata$location_code == nlocations+1),]
    locdata <- specdata[which(specdata$location_code == loc),] # | specdata$location_code == nlocations+1),]
    
    for (mpas in 1:length(MPAsizes)){
      
      # count for indexing 
      count = count+1
      
      # determine mpa size
      mpasize = MPAsizes[mpas]
      
      # allocate data matrix
      mpastats = matrix(nrow=length(MPAsizes),ncol=6)
      
      # prepare modelling environment
      maxdist = max(dists) # maximum travel distance by any individual 
      meanmaxn = mean(locdata$maxn) # mean number of individuals at BRUV
      #medianmaxn = median(locdata$maxn) # mean number of individuals at BRUV
      sdmaxn = sd(locdata$maxn) # std of individuals at BRUV
      minlength = mpasize + 2 * maxdist # minimum extent of simulated coastline
      if (mean.max.extent == 1) {minlength = mpasize + 2 * mean(dists)}
      ndrops = 2*ceiling((minlength/rBRUVcatchment)/2)+1 # ensure sufficent hypothetical BRUV samples - odd number
      coastline = seq(1,(ndrops-1)*rBRUVcatchment+1,1) # index for hypothetical coastline in 1 m resolution
      mpalocs = floor(seq(length(coastline)/2-mpasize/2,length(coastline)/2+mpasize/2,1)) # index of mpa locations
      droplocs = seq(1,length(coastline),rBRUVcatchment) # BRUV locations
      
      # allocate data matrices
      snsharks = matrix(nrow = nreplicates,ncol=1)
      nfull = snsharks
      npart = snsharks
      meanp_survival = matrix(ncol=length(pMorts),nrow=nreplicates)
      sdp_survival = meanp_survival
      p95p_survival = meanp_survival
      p05p_survival = meanp_survival
      sumn_survival = meanp_survival
      propn_mortality = meanp_survival 
      minp_survival = meanp_survival
      maxp_survival = meanp_survival
      repdata = matrix(nrow = nreplicates, ncol= 3)
      pidata = matrix(nrow=0,ncol=2) # allocate matrix for individual protection across all replicates
      
      for (rep in 1:nreplicates){
        
        # hypothtetical sampling - drop by drop
        nsharks = 0 # set to zero before each sampling
        hpidata = matrix(nrow=0,ncol=2) # allocate matrix for individual home range and protection
        ssharklocs = matrix(nrow=0,ncol=2)
        
        for (drop in 1:ndrops){
          
          #overlap = 0 # could keep track of individuals overlapping multiple BRUVS
          droploc = droplocs[drop] # index of BRUV location
          sharksAtBRUV = locdata$maxn[round(runif(1,1,dim(locdata)[1]))] # shark abundance at BRUV
          #nsharks = nsharks + sharksAtBRUV  # keep track of total shark numbers sampled
          
          if (sharksAtBRUV > 0){ # if there are any sharks at BRUV
            
            # determine home ranges
            idata = matrix(nrow=sharksAtBRUV,ncol=2) # empty matrix for saving individual-based data 
            hrids = round(runif(sharksAtBRUV,1,length(dists)))
            idata[,1] = floor(dists[hrids]) # individual home ranges
            
            # determine exposure of each individualioutput
            for (sharkAtBRUV in 1:sharksAtBRUV){
              
              hr = idata[sharkAtBRUV,1] # home range of individual
              #hr = 2*floor(hr/2)+1 # round home range to nearest odd integer
              hrprobs = sprobslist[[s]][[hrids[sharkAtBRUV]]] # relative probability of occurrence along hr 
              hrprobs <- hrprobs[!hrprobs %in% NA]
              
              # randloc = round(runif(1,1,hr)) # determine random location
              #print(length(seq(1,hr,1)))
              #print(length(hrprobs))
              
              hrlocs = seq(1,hr,1) # assign home range vector of cell ids
              randloc = sample(hrlocs,size=1,replace=TRUE,prob=hrprobs) # pick a cell where shark is according to probability 
              
              sharklocs = seq(droploc-randloc+1,droploc+hr-randloc,1) # assign activity space (overlapping BRUV) along modelling environment
              
              if(plot.example == 1 && s == nspecs && lc == nlocations && mpasize == 1000) {ssharklocs = rbind(ssharklocs,c(min(sharklocs),max(sharklocs)))}
              
              # calculate overlap of MPA with max travel distance
              idata[sharkAtBRUV,2] = sum(is.element(sharklocs,mpalocs))/hr # calculate proportion of activity space overlapping MPA 
              
              # recalculate overlap by considering movement probability along the travel distance continuum
              if (probs.per.ind == 1){
                iprotect = idata[sharkAtBRUV,2]
                #idata[sharkAtBRUV,2] = sum(sprobslist[[s]][[hrids[sharkAtBRUV]]][1:round(iprotect*hr)],na.rm = TRUE)
                #idata[sharkAtBRUV,2] = sum(hrprobs[1:round(iprotect*hr)],na.rm = TRUE)
                idata[sharkAtBRUV,2] = sum(hrprobs[hrlocs[is.element(sharklocs,mpalocs)]],na.rm = TRUE)
                #print(c(iprotect,idata[sharkAtBRUV,2]))
              }
              
            } # individual loop
            hpidata = rbind(hpidata,idata) # store individual data at all bruvs in matrix
          } # if sharks sampled
          nsharks = nsharks + sharksAtBRUV  # keep track of total shark numbers sampled
        } # drop loop
        
        #print(c(ndrops,nsharks))
        
        if (nsharks > 0) { # if there were any sharks at any BRUV 
          # save all individual results in big table: c("species", "location_code", "replicate", "maxdist", "mpa_size", "home range", "protection")
          itable = matrix(c(spec,loc,rep,mpasize),ncol=4)
          if (nsharks > 1) {itable = cbind(itable[rep(1:nrow(itable), times = nsharks), ],hpidata)}
          #pidata = rbind(pidata,hpidata[,2])}
          else {itable = cbind(itable,hpidata)}
          ioutput = rbind(ioutput,itable)
          #pidata = rbind(pidata,hpidata[,2])
        }
        
        pidata = rbind(pidata,hpidata)
        
        # save stats per replicate
        
        if (drop == ndrops){
          
          snsharks[rep] = nsharks 
          nfull[rep] = sum(hpidata[,2]==1) # number of fully protected individuals
          npart[rep] = sum(hpidata[,2]>0) # number of partially protected individuals
          repdata[rep,] = c(nsharks,nfull[rep],npart[rep]) 
          
          # quantify chances of survival and mortality of fully and partially protected individuals
          exposure = (1-hpidata[hpidata[,2]>0,2]) # of partially protected individuals
          exposuremat = do.call(cbind, replicate(length(pMorts), exposure, simplify=FALSE)) 
          pmortmat = do.call(rbind, replicate(length(exposure), pMorts, simplify=FALSE))
          pstable = (1 - exposuremat * pmortmat) # probability of survival 
          stables = matrix(rbinom(1:length(pstable)*years,size = 1, prob = pstable),dim(pstable)[1],dim(pstable)[2]) # survival drawn from probability
          for (addyears in 2:years) {
            stables = stables+matrix(rbinom(1:length(pstable)*years,size = 1, prob = pstable),dim(pstable)[1],dim(pstable)[2]) # survival drawn from probability
          }    
          stable = matrix(as.integer(stables==years),dim(pstable)[1],dim(pstable)[2]) 
          
          meanp_survival[rep,] = apply(pstable,2,mean)
          sdp_survival[rep,] = apply(pstable,2,sd)
          p95p_survival[rep,] = apply(pstable,2,quantile,0.95)
          p05p_survival[rep,] = apply(pstable,2,quantile,0.05)
          minp_survival[rep,] = apply(pstable,2,min)
          maxp_survival[rep,] = apply(pstable,2,max)
          sumn_survival[rep,] = apply(stable,2,sum)
          propn_mortality[rep,] = 1-sumn_survival[rep,]/sumn_survival[rep,1] 
        }
        
      } # replicate loop
      
      # store stats for each replicate
      reptable = matrix(c(spec,loc,rep,round(maxdist),mpasize),ncol=5)
      reptable = cbind(reptable[rep(1:nrow(reptable), times = nreplicates), ],repdata)
      
      # calculate stats across replicates
      mnsharks = mean(snsharks)
      mnfull = mean(nfull)
      mnpart = mean(npart)
      mprot = mean(pidata[pidata[,2]>0,2])
      sdprot = sd(pidata[pidata[,2]>0,2])
      medprot = median(pidata[pidata[,2]>0,2])
      mpastats[mpas,] = c(mnsharks,mnfull,mnpart,mprot,sdprot,medprot)
      
      meanmeanp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,mean,na.rm=TRUE))
      minminp_survival[count,] = c(spec,loc,mpasize,apply(minp_survival,2,min,na.rm=TRUE))
      maxmaxp_survival[count,] = c(spec,loc,mpasize,apply(maxp_survival,2,max,na.rm=TRUE))
      meanminp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,min,na.rm=TRUE))
      meanmaxp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,max,na.rm=TRUE))
      p05meanp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,quantile,0.05,na.rm=TRUE))
      p95meanp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,quantile,0.95,na.rm=TRUE))
      minmeanp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,min,na.rm=TRUE))
      maxmeanp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,max,na.rm=TRUE))
      meansdp_survival[count,] = c(spec,loc,mpasize,apply(sdp_survival,2,mean,na.rm=TRUE))
      meanp05p_survival[count,] = c(spec,loc,mpasize,apply(p05p_survival,2,mean,na.rm=TRUE))
      meanp95p_survival[count,] = c(spec,loc,mpasize,apply(p95p_survival,2,mean,na.rm=TRUE))
      
      meansumn_survival[count,] = c(spec,loc,mpasize,apply(sumn_survival,2,mean,na.rm=TRUE))
      mediansumn_survival[count,] = c(spec,loc,mpasize,apply(sumn_survival,2,median,na.rm=TRUE))
      sdsumn_survival[count,] = c(spec,loc,mpasize,apply(sumn_survival,2,sd,na.rm=TRUE))
      p05sumn_survival[count,] = c(spec,loc,mpasize,apply(sumn_survival,2,quantile,0.05,na.rm=TRUE))
      p95sumn_survival[count,] = c(spec,loc,mpasize,apply(sumn_survival,2,quantile,0.95,na.rm=TRUE))
      
      meanpropn_mortality[count,] = c(spec,loc,mpasize,apply(propn_mortality,2,mean,na.rm=TRUE))
      medianpropn_mortality[count,] = c(spec,loc,mpasize,apply(propn_mortality,2,median,na.rm=TRUE))
      sdpropn_mortality[count,] = c(spec,loc,mpasize,apply(propn_mortality,2,sd,na.rm=TRUE))
      p05propn_mortality[count,] = c(spec,loc,mpasize,apply(propn_mortality,2,quantile,0.05,na.rm=TRUE))
      p95propn_mortality[count,] = c(spec,loc,mpasize,apply(propn_mortality,2,quantile,0.95,na.rm=TRUE))
      
      statstable[count,] = c(spec,loc,meanmaxn,sdmaxn,round(maxdist),mpasize,mpastats[mpas,])
      
      moffset = (1-mprot)*pMorts
      sdoffset = (1-mprot)*pMorts
      mean_mortality_offset[count,] = c(spec,loc,mpasize,moffset)
      sd_mortality_offset[count,] = c(spec,loc,mpasize,sdoffset)
      mean_mortalityInst_offset[count,] = c(spec,loc,mpasize,-log(1-moffset))
      sd_mortalityInst_offset[count,] = c(spec,loc,mpasize,-log(1-sdoffset))
      
      # plot results
      if (plot.results == 1 && mpas == length(MPAsizes)) {
        pMortIDs = c(1:(length(pMorts)-1)+4) # fishing risks to plot
        legendtext =  as.character(pMortsInst[pMortIDs-3]) # text for legend
        xmaxlimval = 100000
        mpaloc = which(MPAsizes==xmaxlimval)
        
        # plot survival
        ymaxlim = ceiling(max(meansumn_survival[seq(count-mpas+1,count-mpas+mpaloc,1),pMortIDs[1]]))
        if (save.plots == 1) {
          png(paste0(savedir, speclatname, "_", locname , "_abundance_mean_protected_individuals.png"))
        }
        plot(MPAsizes,meansumn_survival[seq(count-mpas+1,count,1),pMortIDs[1]],type="l",lty=1,xlab="MPA size (m)",ylab="Protected individuals",main=paste0(specname, "; ", locname , " abundance."), xlim=c(0,xmaxlimval),ylim=c(0,ymaxlim)) 
        for (i in 2:length(pMortIDs)){
          lines(MPAsizes,meansumn_survival[seq(count-mpas+1,count,1),pMortIDs[i]],lty=i)   
        } 
        legend("topleft", legend = legendtext, lty = seq(1,length(legendtext),1), title = "Fishing mortality")
        if (save.plots == 1) {dev.off()}
        
        # plot mean protection (proportion of range within MPA)
        sdlow = statstable[seq(count-mpas+1,count,1),10] - statstable[seq(count-mpas+1,count,1),11] #/sqrt(nreplicates) 
        sdhigh = statstable[seq(count-mpas+1,count,1),10] + statstable[seq(count-mpas+1,count,1),11]#/sqrt(nreplicates) 
        sdhigh[sdhigh > 1] = 1 # cap values at maximum
        sdlow[sdlow < 0] = 0 # cap values at minimum
        
        if (save.plots == 1) {
          png(paste0(savedir, speclatname, "_", locname , "_abundance_mean_range_in_mpa.png"))
        }
        plot(MPAsizes,statstable[seq(count-mpas+1,count,1),10],type="l",lty = 1, xlab="MPA size (m)",ylab="Protection in MPA (mean ? std)",main=paste0(specname, "; ", locname , " abundance."), xlim=c(0,xmaxlimval), ylim=c(0,1))
        #lines(MPAsizes,meanp05p_survival[seq(count-mpas+1,count,1),pMortIDs[length(pMortIDs)]],lty=2)
        #lines(MPAsizes,meanp95p_survival[seq(count-mpas+1,count,1),pMortIDs[length(pMortIDs)]],lty=2)
        lines(MPAsizes,sdlow,lty=2)
        lines(MPAsizes,sdhigh,lty=2)
        if (save.plots == 1) {dev.off()}
        
        # plot mortality
        ymaxlim = 1
        if (save.plots == 1) {
          png(paste0(savedir, speclatname, "_", locname , "_abundance_mean_Fmortality.png"))
        }
        plot(MPAsizes,meanpropn_mortality[seq(count-mpas+1,count,1),pMortIDs[1]],type="l",lty = 1, xlab="MPA size (m)",ylab="Lifetime fishing mortality",main=paste0(specname, "; ", locname , " abundance."), xlim=c(0,xmaxlimval), ylim=c(0,ymaxlim))
        for (i in 2:length(pMortIDs)){
          lines(MPAsizes,meanpropn_mortality[seq(count-mpas+1,count,1),pMortIDs[i]],lty=i) 
        }
        legend("topright", legend = legendtext, lty = seq(1,length(legendtext),1), title = "Fishing mortality")
        if (save.plots == 1) {dev.off()}
        
        # plot mean probability of survival
        # if (save.plots == 1) {
        #    png(paste0(savedir, speclatname, "_", locname , "_abundance_mean_survival_probability.png"))
        #  }
        #  plot(MPAsizes,meanmeanp_survival[seq(count-mpas+1,count,1),pMortIDs[1]],type="b",lty = 1, xlab="MPA size (m)",ylab="Probability of survival",main=paste0(specname, "; ", locname , " abundance."), xlim=c(0,xmaxlimval), ylim=c(0,1))
        #  for (i in 2:length(pMortIDs)){
        #    lines(MPAsizes,meanmeanp_survival[seq(count-mpas+1,count,1),pMortIDs[i]],lty=i) 
        #  }
        #  legend("bottomright", legend = legendtext, lty = seq(1,length(legendtext),1), title = "Probability of catch ")
        #  if (save.plots == 1) {dev.off()}
        
        
        # plot mortality offset
        if (save.plots == 1) {
          png(paste0(savedir, speclatname, "_", locname , "_abundance_mean_mortality_offset.png"))
        }
        plot(MPAsizes,mean_mortalityInst_offset[seq(count-mpas+1,count,1),pMortIDs[1]],type="l",lty = 1, xlab="MPA size (m)",ylab="Fishing mortality",main=paste0(specname, "; ", locname , " abundance."), xlim=c(0,xmaxlimval), ylim=c(0,max(pMortsInst)))
        for (i in 2:length(pMortIDs)){
          lines(MPAsizes,mean_mortalityInst_offset[seq(count-mpas+1,count,1),pMortIDs[i]],lty=i) 
        }
        #legend("topright", legend = legendtext, lty = seq(1,length(legendtext),1), title = "Fishing morality ")
        if (save.plots == 1) {dev.off()}
        
      }
      
      # computing time status report
      mpatime <- Sys.time()
      print(paste0("Species ", s, "/" , nspecs, ", abundance ", lc, "/" , nlocations, ": ", "MPA size ", mpas, "/" , length(MPAsizes) , " completed."))
      print(mpatime - starttime)  
      
    } # mpa loop   
    
  } # location loop
  
  #} # species loop
  
  endtime <- Sys.time()
  print(endtime - starttime)
  
  if (save.results == 1){
    
    # convert results to data frames
    ioutput <- as.data.frame(ioutput)
    colnames(ioutput) = ioutputHeader
    
    statstable <- as.data.frame(statstable)
    colnames(statstable) = statstableHeader
    
    mean_survival_probability = as.data.frame(meanmeanp_survival, row.names = FALSE) 
    colnames(mean_survival_probability) = statsmatsHeader
    
    sd_survival_probability = as.data.frame(meansdp_survival, row.names = FALSE) 
    colnames(sd_survival_probability) = statsmatsHeader
    
    mean_protected_individuals = as.data.frame(meansumn_survival, row.names = FALSE)
    colnames(mean_protected_individuals) = statsmatsHeader
    
    median_protected_individuals = as.data.frame(mediansumn_survival, row.names = FALSE)
    colnames(median_protected_individuals) = statsmatsHeader
    
    sd_protected_individuals = as.data.frame(sdsumn_survival, row.names = FALSE)
    colnames(sd_protected_individuals) = statsmatsHeader
    
    mean_Fmortality = as.data.frame(meanpropn_mortality, row.names = FALSE)
    colnames(mean_Fmortality) = statsmatsHeader
    
    median_Fmortality = as.data.frame(medianpropn_mortality, row.names = FALSE)
    colnames(median_Fmortality) = statsmatsHeader
    
    sd_Fmortality = as.data.frame(sdpropn_mortality, row.names = FALSE)
    colnames(sd_Fmortality) = statsmatsHeader
    
    mean_Foffset = as.data.frame(mean_mortalityInst_offset, row.names = FALSE)
    colnames(mean_Foffset) = statsmatsHeader
    
    sd_Foffset = as.data.frame(sd_mortalityInst_offset, row.names = FALSE)
    colnames(sd_Foffset) = statsmatsHeader
    
    other.parameters <-  as.data.frame(t(c(nreplicates,years,maxd.per.ind, probs.per.ind)))
    colnames(other.parameters) <- c("num_replicates", "lifetime_years", "max_distance_per_individual","movement_probability_per_individual")
    
    write.table(other.parameters,paste0(savedir,"parameters.csv"), sep =",", row.names = FALSE)
    write.table(ioutput,paste0(savedir,"results_by_individual.csv"), sep =",", row.names = FALSE)
    write.table(statstable,paste0(savedir,"mean_results.csv"), sep =",", row.names = FALSE)
    write.table(mean_survival_probability,paste0(savedir,"mean_survival_probability.csv"), sep =",", row.names = FALSE)
    write.table(sd_survival_probability,paste0(savedir,"sd_survival_probability.csv"), sep =",", row.names = FALSE)
    write.table(mean_protected_individuals,paste0(savedir,"mean_protected_individuals.csv"), sep =",", row.names = FALSE)
    write.table(sd_protected_individuals,paste0(savedir,"sd_protected_individuals.csv"), sep =",", row.names = FALSE)
    #write.table(median_protected_individuals,paste0(savedir,"median_protected_individuals.csv"), sep =",", row.names = FALSE)
    write.table(mean_Fmortality,paste0(savedir,"mean_Fmortality.csv"), sep =",", row.names = FALSE)
    write.table(sd_Fmortality,paste0(savedir,"sd_Fmortality.csv"), sep =",", row.names = FALSE)
    #write.table(median_Fmortality,paste0(savedir,"median_Fmortality.csv"), sep =",", row.names = FALSE)
    write.table(mean_Foffset,paste0(savedir,"mean_Foffset.csv"), sep =",", row.names = FALSE)
    write.table(sd_Foffset,paste0(savedir,"sd_Foffset.csv"), sep =",", row.names = FALSE)
    
  } # save results
  #return(probsdflist)
  #return(do.call(rbind,probsdflist)) # Brings back list of movement probabilities for each shark
  
}# Function end


#######################################################################


fsharkABM_v502_sites <- function(bruvdat,
                           movedat,
                           dispdat_limit=1,
                           abundancecat=2, # How many abundance categories are there? 1,2,3?
                           speciesToSimulate, # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse 
                           rBRUVcatchment=400/resolution, # diameter of BRUV plume catchment area: 200,400,600 
                           years=25, # lifetime over which mortality is assessed
                           save.results = 1, # save all results 
                           plot.example = 0, # plot example mpa scenario to validate modelling procedure
                           plot.results = 1, # plot key mortality outcomes and surviving individuals
                           save.plots = 1, # save all plots to png file - not showing then
                           maxd.per.ind = 1, # calculate max travel distance per individual at any day 
                           probs.per.ind = 1, # The activity Range: 0 = uniform (use max dispersal distance) or 1 = ^ shaped (uses daily distances around fixed point)  
                           nreplicates = 1000,
                           mean.max.extent = 1, # use mean of max daily dispersal instead of max of max
                           resolution = 1,
                           resfolder = 'results') # resolution of modelling environment in m
{ # calculate movement probabilities per individual)
  
  # resolution = 1
  # bruvdat = maxndata_NS
  # movedat = maxddata
  # dispdat_limit=7
  # abundancecat = 2 #1,2,3 (same, highlow, highmediumlow)
  # speciesToSimulate = 7   #speciesToSimulate = c(1,2,3,6,7) # 1 whitetip, 2 blacktip, 3 grey reef shark, 6 caribbean reef, 7 nurse
  # rBRUVcatchment = 400/resolution # radius of BRUV plume catchment area: 200,400,600
  # years = 25 # lifetime over which mortality is assessed
  # mortRate = c("pMortsInst") # "pMorts","pMortsInst"
  # save.results = 1 # save all results
  # plot.example = 0 # plot example mpa scenario to validate modelling procedure
  # plot.results = 1 # plot key mortality outcomes and surviving individuals
  # save.plots = 1 # save all plots to png file - not showing then
  # maxd.per.ind = 1 # calculate max travel distance per individual at any day
  # probs.per.ind = 0
  # nreplicates = 10 # number of replicate simulations per species, location and MPA size
  # #
  # Modelling parameters
  
  # Removes individuals with fewer than X days of mobements
  dfmaxd <- as.data.frame(table(movedat$tag_id))
  exclude_ids <- dfmaxd$Var1[which(dfmaxd$Freq<dispdat_limit)]
  '%!in%' <- function(x,y)!('%in%'(x,y))
  movedat1 <- movedat[movedat$tag_id %!in% exclude_ids, ]
  
  MPAsizes = c(0,seq(500,2000,500),seq(3000,5000,1000),seq(10000,20000,5000),c(30000,40000,50000,100000))/resolution
  # MPAsizes = c(2000,20000,100000)/resolution
  
  #pMortsInst = seq(0,.1,.3) # probability of mortality when exposed to fishing
  pMortsInst = seq(0,.3,.05) # probability of mortality when exposed to fishing
  
  #if(mortRate=="pMorts")
  #  pMorts = 1-exp(-pMortsInst)
  #if(mortRate=="pMortsInst")
  #  pMortsInst = -log(1-pMorts) 
  pMorts = 1-exp(-pMortsInst)
  pMortsInst = -log(1-pMorts) # convert discrete to instantaneous mortality rates (Q for Nils - do I have to do this step?)
  
  # Select abundance category
  if(abundancecat==1) # No category
  {
    bruvdat$location_code <- 2 #Standardises locations
    abundance = c(2) #   2 average abundance
  }
  if(abundancecat==2){ # 2 categories
    abundance = unique(bruvdat$location_code)
    bruvdat$location_code # 1 = low abundance, 2 = average abundance
  }
  if(abundancecat==3){ # 3 categories
    abundance = unique(bruvdat$location_code_high) #  1 low abundance, 2 average abundance, 3 high abundance 
    location_code = bruvdat$location_code_high # Rewrite column with 2 levels as 3 levels
  }
  
  # explore nmax data
  lat.spec.names <- unique(bruvdat$species)
  comm.spec.names <- unique(bruvdat$common_name)
  #speciesToAnalyse = c(as.character(lat.spec.names[speciesToSimulate]),"") # has to include empty cells melanopterus","") #c("obesus","melanopterus","amblyrhynchos","perezi") 
  speciesToAnalyse = c(as.character(lat.spec.names[1])) # has to include empty cells melanopterus","") #c("obesus","melanopterus","amblyrhynchos","perezi") 
  
  region.names <- unique(bruvdat$region)
  location.names <- unique(bruvdat$location_name)
  location.code.names <- c("low", "average", "high") 
  nspecs <- length(speciesToAnalyse) #updated when turning into function
  nlocations <- length(abundance) # coded as 1 (low), 2(med), 3(high) in vector location_code for low, medium, high abundance 
  
  # allocate output and working matrices
  sprobslist <- vector("list", nspecs)
  probslist <- vector("list")
  probsdflist <- vector("list")
  
  ioutputHeader <- c("species", "abundance", "replicate", "mpa_size", "home_range", "protection")
  ioutput = matrix(nrow=0,ncol=length(ioutputHeader)) # output table for individual data 
  statstableHeader <- c("species", "abundance", "mean_maxn", "sd_maxn", "max_hr", "mpa_size", "mean_nsharks", "mean_nfull", "mean_npart","mean_prot", "sd_prot", "median_prot")
  statstable = matrix(nrow=nspecs*nlocations*length(MPAsizes),ncol=length(statstableHeader)) # allocate output table for individual data 
  statsmatsHeader = c("species", "abundance", "mpa_size",as.character(pMortsInst))
  
  meanmeanp_survival = matrix(data=0,nrow=nspecs*nlocations*length(MPAsizes),ncol=length(statsmatsHeader))
  
  minminp_survival = meanmeanp_survival
  maxmaxp_survival = meanmeanp_survival
  meanminp_survival = meanmeanp_survival
  meanmaxp_survival = meanmeanp_survival
  p05meanp_survival = meanmeanp_survival
  p95meanp_survival = meanmeanp_survival
  minmeanp_survival = meanmeanp_survival
  maxmeanp_survival = meanmeanp_survival
  meansdp_survival = meanmeanp_survival
  meanp05p_survival = meanmeanp_survival
  meanp95p_survival = meanmeanp_survival
  
  meansumn_survival = meanmeanp_survival
  mediansumn_survival = meanmeanp_survival
  sdsumn_survival = meanmeanp_survival
  p05sumn_survival = meanmeanp_survival
  p95sumn_survival = meanmeanp_survival
  
  meanpropn_mortality = meanmeanp_survival
  medianpropn_mortality = meanmeanp_survival
  sdpropn_mortality = meanmeanp_survival
  p05propn_mortality = meanmeanp_survival
  p95propn_mortality = meanmeanp_survival
  
  mean_mortality_offset = matrix(nrow=nspecs*nlocations*length(MPAsizes),ncol=length(statsmatsHeader))
  sd_mortality_offset = mean_mortality_offset
  mean_mortalityInst_offset = sd_mortality_offset
  sd_mortalityInst_offset = sd_mortality_offset
  
  # create folders to save results
  if(!dir.exists(paste0(datadir, resfolder))){dir.create(paste0(datadir, resfolder))}
  if(!dir.exists(paste0(datadir,resfolder, "/rBRUVcatchment_",rBRUVcatchment))){dir.create(paste0(datadir,resfolder, "/rBRUVcatchment_",rBRUVcatchment))}
  
  
  ##################
  # START OF MODEL #
  ##################
  
  ## Only have to run this once at the beginning
  
  # count for indexing purposes
  count = 0
  
  # species loop
  #for (s in 1:nspecs){  # Removed species loop to run as part of a function
  s <- 1 # This has been standardised to a 1 as a result of turning the loop into a function which runs independently on each species
  
  # create and assign folders to save results
  if(!dir.exists(paste0(datadir,resfolder, "/rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))){
    dir.create(paste0(datadir,resfolder, "/rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))}
  savedir = paste0(datadir,resfolder, "/rBRUVcatchment_", rBRUVcatchment,"/", speciesToAnalyse[s],"/") 
  
  starttime <- Sys.time() # save time
  
  # select data
  spec <- speciesToSimulate[s] 
  speclocs <- which(bruvdat$species == speciesToAnalyse[s] | bruvdat$species == speciesToAnalyse[nspecs+1])
  #specdata <- bruvdat[speclocs,] # step set up to deal with old rel abundance data for all species
  specdata <- bruvdat # new step set up to deal with new species-specific abundance data for all species
  specname <- comm.spec.names[1]
  speclatname <- speciesToAnalyse[s]
  
  if (maxd.per.ind == 1){
    ids = unique(movedat1$tag_id[movedat1$species == speciesToAnalyse[s]]) # get shark ids
    dists = matrix(nrow=length(ids),ncol=1)
    for (id in 1:length(ids)) {dists[id] = max(movedat1$max.dispersal_km[movedat1$tag_id == ids[id]]) * 1000}
  }else{
    dists = movedat1$max.dispersal_km[movedat1$species == speciesToAnalyse[s]] * 1000 / resolution
  } # travel distances in m
  #hist(dists,breaks=length(dists)/4)
  #length(dists)
  
  
  if (maxd.per.ind == 1 & probs.per.ind == 1){
    # generate probability of occurrence along the distance spectrum
    for (id in 1:length(ids)) {
      alldists = movedat1$max.dispersal_km[movedat1$tag_id == ids[id]] * 1000 / resolution
      # hist(alldists, breaks = 20, xlab = "Travel distance (m)", ylab = "Frequency")
      kernel = as.matrix(rep(1,1,max(alldists)))
      centreloc = ceiling(max(alldists)/2)
      for (dshark in 1:length(unique(alldists))) { # for each observed travel distance (RD WHY dshark???)
        addval = 1 # set initial number of observed vals to 1
        maxlocs = which(alldists == max(alldists)) # identify all distance locs
        startloc = round(centreloc - alldists[maxlocs[1]]/2) # define kernel start loc
        endloc = floor(startloc + alldists[maxlocs[1]]) # define kernel end loc
        if (length(maxlocs) > 1){addval = length(maxlocs)} # if observed more than once update value
        kernel[startloc:endloc] = kernel[startloc:endloc] + addval # add observations along kernel spectrum
        alldists = alldists[alldists != alldists[maxlocs[1]]] # delete observed distance from list
      }
      probslist[[id]] = kernel/sum(kernel,na.rm=TRUE) # get cumulative probability spectrum by dividing by sum
      
      # generate movement probability dataframes for plotting
      probsdf <- data.frame(ID=as.factor(ids[id]),
                            DISTANCE=1:length(probslist[[id]]),
                            P=probslist[[id]])
      probsdflist[[id]] <- probsdf
      
    } # end of individual loop
    sprobslist[[s]] = probslist
    
    probslist <- NULL #RD Reset these values to fix bug in v404
    
  } # end of probs calculations
  
  # location loop (# coded as 1, 2, 3 in vector location_code for low, medium, high abundance)
  for (lc in 1:nlocations){
    
    # select data (low to high abundance)
    #slocs = sum(which(specdata$location_code == lc | specdata$location_code == nlocations+1))
    
    
    # THis is where to add flexibility to set the different High Medium Low categories
    
    loc = abundance[lc]
    locname = location.code.names[loc]
    #locdata <- specdata[which(specdata$location_code == loc | specdata$location_code == 4),] # | specdata$location_code == nlocations+1),]
    locdata <- specdata[which(specdata$location_code == loc),] # | specdata$location_code == nlocations+1),]
    
    for (mpas in 1:length(MPAsizes)){
      
      # count for indexing 
      count = count+1
      
      # determine mpa size
      mpasize = MPAsizes[mpas]
      
      # allocate data matrix
      mpastats = matrix(nrow=length(MPAsizes),ncol=6)
      
      # prepare modelling environment
      maxdist = max(dists) # maximum travel distance by any individual 
      meanmaxn = mean(locdata$maxn) # mean number of individuals at BRUV
      #medianmaxn = median(locdata$maxn) # mean number of individuals at BRUV
      sdmaxn = sd(locdata$maxn) # std of individuals at BRUV
      minlength = mpasize + 2 * maxdist # minimum extent of simulated coastline
      if (mean.max.extent == 1) {minlength = mpasize + 2 * mean(dists)}
      ndrops = 2*ceiling((minlength/rBRUVcatchment)/2)+1 # ensure sufficent hypothetical BRUV samples - odd number
      coastline = seq(1,(ndrops-1)*rBRUVcatchment+1,1) # index for hypothetical coastline in 1 m resolution
      mpalocs = floor(seq(length(coastline)/2-mpasize/2,length(coastline)/2+mpasize/2,1)) # index of mpa locations
      droplocs = seq(1,length(coastline),rBRUVcatchment) # BRUV locations
      
      # allocate data matrices
      snsharks = matrix(nrow = nreplicates,ncol=1)
      nfull = snsharks
      npart = snsharks
      meanp_survival = matrix(ncol=length(pMorts),nrow=nreplicates)
      sdp_survival = meanp_survival
      p95p_survival = meanp_survival
      p05p_survival = meanp_survival
      sumn_survival = meanp_survival
      propn_mortality = meanp_survival 
      minp_survival = meanp_survival
      maxp_survival = meanp_survival
      repdata = matrix(nrow = nreplicates, ncol= 3)
      pidata = matrix(nrow=0,ncol=2) # allocate matrix for individual protection across all replicates
      
      for (rep in 1:nreplicates){
        
        # hypothtetical sampling - drop by drop
        nsharks = 0 # set to zero before each sampling
        hpidata = matrix(nrow=0,ncol=2) # allocate matrix for individual home range and protection
        ssharklocs = matrix(nrow=0,ncol=2)
        
        for (drop in 1:ndrops){
          
          #overlap = 0 # could keep track of individuals overlapping multiple BRUVS
          droploc = droplocs[drop] # index of BRUV location
          sharksAtBRUV = locdata$maxn[round(runif(1,1,dim(locdata)[1]))] # shark abundance at BRUV
          #nsharks = nsharks + sharksAtBRUV  # keep track of total shark numbers sampled
          
          if (sharksAtBRUV > 0){ # if there are any sharks at BRUV
            
            # determine home ranges
            idata = matrix(nrow=sharksAtBRUV,ncol=2) # empty matrix for saving individual-based data 
            hrids = round(runif(sharksAtBRUV,1,length(dists)))
            idata[,1] = floor(dists[hrids]) # individual home ranges
            
            # determine exposure of each individualioutput
            for (sharkAtBRUV in 1:sharksAtBRUV){
              
              hr = idata[sharkAtBRUV,1] # home range of individual
              #hr = 2*floor(hr/2)+1 # round home range to nearest odd integer
              hrprobs = sprobslist[[s]][[hrids[sharkAtBRUV]]] # relative probability of occurrence along hr 
              hrprobs <- hrprobs[!hrprobs %in% NA]
              
              # randloc = round(runif(1,1,hr)) # determine random location
              #print(length(seq(1,hr,1)))
              #print(length(hrprobs))
              
              hrlocs = seq(1,hr,1) # assign home range vector of cell ids
              randloc = sample(hrlocs,size=1,replace=TRUE,prob=hrprobs) # pick a cell where shark is according to probability 
              
              sharklocs = seq(droploc-randloc+1,droploc+hr-randloc,1) # assign activity space (overlapping BRUV) along modelling environment
              
              if(plot.example == 1 && s == nspecs && lc == nlocations && mpasize == 1000) {ssharklocs = rbind(ssharklocs,c(min(sharklocs),max(sharklocs)))}
              
              # calculate overlap of MPA with max travel distance
              idata[sharkAtBRUV,2] = sum(is.element(sharklocs,mpalocs))/hr # calculate proportion of activity space overlapping MPA 
              
              # recalculate overlap by considering movement probability along the travel distance continuum
              if (probs.per.ind == 1){
                iprotect = idata[sharkAtBRUV,2]
                #idata[sharkAtBRUV,2] = sum(sprobslist[[s]][[hrids[sharkAtBRUV]]][1:round(iprotect*hr)],na.rm = TRUE)
                #idata[sharkAtBRUV,2] = sum(hrprobs[1:round(iprotect*hr)],na.rm = TRUE)
                idata[sharkAtBRUV,2] = sum(hrprobs[hrlocs[is.element(sharklocs,mpalocs)]],na.rm = TRUE)
                #print(c(iprotect,idata[sharkAtBRUV,2]))
              }
              
            } # individual loop
            hpidata = rbind(hpidata,idata) # store individual data at all bruvs in matrix
          } # if sharks sampled
          nsharks = nsharks + sharksAtBRUV  # keep track of total shark numbers sampled
        } # drop loop
        
        #print(c(ndrops,nsharks))
        
        if (nsharks > 0) { # if there were any sharks at any BRUV 
          # save all individual results in big table: c("species", "location_code", "replicate", "maxdist", "mpa_size", "home range", "protection")
          itable = matrix(c(spec,loc,rep,mpasize),ncol=4)
          if (nsharks > 1) {itable = cbind(itable[rep(1:nrow(itable), times = nsharks), ],hpidata)}
          #pidata = rbind(pidata,hpidata[,2])}
          else {itable = cbind(itable,hpidata)}
          ioutput = rbind(ioutput,itable)
          #pidata = rbind(pidata,hpidata[,2])
        }
        
        pidata = rbind(pidata,hpidata)
        
        # save stats per replicate
        
        if (drop == ndrops){
          
          snsharks[rep] = nsharks 
          nfull[rep] = sum(hpidata[,2]==1) # number of fully protected individuals
          npart[rep] = sum(hpidata[,2]>0) # number of partially protected individuals
          repdata[rep,] = c(nsharks,nfull[rep],npart[rep]) 
          
          # quantify chances of survival and mortality of fully and partially protected individuals
          exposure = (1-hpidata[hpidata[,2]>0,2]) # of partially protected individuals
          exposuremat = do.call(cbind, replicate(length(pMorts), exposure, simplify=FALSE)) 
          pmortmat = do.call(rbind, replicate(length(exposure), pMorts, simplify=FALSE))
          pstable = (1 - exposuremat * pmortmat) # probability of survival 
          stables = matrix(rbinom(1:length(pstable)*years,size = 1, prob = pstable),dim(pstable)[1],dim(pstable)[2]) # survival drawn from probability
          for (addyears in 2:years) {
            stables = stables+matrix(rbinom(1:length(pstable)*years,size = 1, prob = pstable),dim(pstable)[1],dim(pstable)[2]) # survival drawn from probability
          }    
          stable = matrix(as.integer(stables==years),dim(pstable)[1],dim(pstable)[2]) 
          
          meanp_survival[rep,] = apply(pstable,2,mean)
          sdp_survival[rep,] = apply(pstable,2,sd)
          p95p_survival[rep,] = apply(pstable,2,quantile,0.95)
          p05p_survival[rep,] = apply(pstable,2,quantile,0.05)
          minp_survival[rep,] = apply(pstable,2,min)
          maxp_survival[rep,] = apply(pstable,2,max)
          sumn_survival[rep,] = apply(stable,2,sum)
          propn_mortality[rep,] = 1-sumn_survival[rep,]/sumn_survival[rep,1] 
        }
        
      } # replicate loop
      
      # store stats for each replicate
      reptable = matrix(c(spec,loc,rep,round(maxdist),mpasize),ncol=5)
      reptable = cbind(reptable[rep(1:nrow(reptable), times = nreplicates), ],repdata)
      
      # calculate stats across replicates
      mnsharks = mean(snsharks)
      mnfull = mean(nfull)
      mnpart = mean(npart)
      mprot = mean(pidata[pidata[,2]>0,2])
      sdprot = sd(pidata[pidata[,2]>0,2])
      medprot = median(pidata[pidata[,2]>0,2])
      mpastats[mpas,] = c(mnsharks,mnfull,mnpart,mprot,sdprot,medprot)
      
      meanmeanp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,mean,na.rm=TRUE))
      minminp_survival[count,] = c(spec,loc,mpasize,apply(minp_survival,2,min,na.rm=TRUE))
      maxmaxp_survival[count,] = c(spec,loc,mpasize,apply(maxp_survival,2,max,na.rm=TRUE))
      meanminp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,min,na.rm=TRUE))
      meanmaxp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,max,na.rm=TRUE))
      p05meanp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,quantile,0.05,na.rm=TRUE))
      p95meanp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,quantile,0.95,na.rm=TRUE))
      minmeanp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,min,na.rm=TRUE))
      maxmeanp_survival[count,] = c(spec,loc,mpasize,apply(meanp_survival,2,max,na.rm=TRUE))
      meansdp_survival[count,] = c(spec,loc,mpasize,apply(sdp_survival,2,mean,na.rm=TRUE))
      meanp05p_survival[count,] = c(spec,loc,mpasize,apply(p05p_survival,2,mean,na.rm=TRUE))
      meanp95p_survival[count,] = c(spec,loc,mpasize,apply(p95p_survival,2,mean,na.rm=TRUE))
      
      meansumn_survival[count,] = c(spec,loc,mpasize,apply(sumn_survival,2,mean,na.rm=TRUE))
      mediansumn_survival[count,] = c(spec,loc,mpasize,apply(sumn_survival,2,median,na.rm=TRUE))
      sdsumn_survival[count,] = c(spec,loc,mpasize,apply(sumn_survival,2,sd,na.rm=TRUE))
      p05sumn_survival[count,] = c(spec,loc,mpasize,apply(sumn_survival,2,quantile,0.05,na.rm=TRUE))
      p95sumn_survival[count,] = c(spec,loc,mpasize,apply(sumn_survival,2,quantile,0.95,na.rm=TRUE))
      
      meanpropn_mortality[count,] = c(spec,loc,mpasize,apply(propn_mortality,2,mean,na.rm=TRUE))
      medianpropn_mortality[count,] = c(spec,loc,mpasize,apply(propn_mortality,2,median,na.rm=TRUE))
      sdpropn_mortality[count,] = c(spec,loc,mpasize,apply(propn_mortality,2,sd,na.rm=TRUE))
      p05propn_mortality[count,] = c(spec,loc,mpasize,apply(propn_mortality,2,quantile,0.05,na.rm=TRUE))
      p95propn_mortality[count,] = c(spec,loc,mpasize,apply(propn_mortality,2,quantile,0.95,na.rm=TRUE))
      
      statstable[count,] = c(spec,loc,meanmaxn,sdmaxn,round(maxdist),mpasize,mpastats[mpas,])
      
      moffset = (1-mprot)*pMorts
      sdoffset = (1-mprot)*pMorts
      mean_mortality_offset[count,] = c(spec,loc,mpasize,moffset)
      sd_mortality_offset[count,] = c(spec,loc,mpasize,sdoffset)
      mean_mortalityInst_offset[count,] = c(spec,loc,mpasize,-log(1-moffset))
      sd_mortalityInst_offset[count,] = c(spec,loc,mpasize,-log(1-sdoffset))
      
      # plot results
      if (plot.results == 1 && mpas == length(MPAsizes)) {
        pMortIDs = c(1:(length(pMorts)-1)+4) # fishing risks to plot
        legendtext =  as.character(pMortsInst[pMortIDs-3]) # text for legend
        xmaxlimval = 100000
        mpaloc = which(MPAsizes==xmaxlimval)

        # plot survival
        ymaxlim = ceiling(max(meansumn_survival[seq(count-mpas+1,count-mpas+mpaloc,1),pMortIDs[1]]))
        if (save.plots == 1) {
          png(paste0(savedir, speclatname, "_", locname , "_abundance_mean_protected_individuals.png"))
        }
        plot(MPAsizes,meansumn_survival[seq(count-mpas+1,count,1),pMortIDs[1]],type="l",lty=1,xlab="MPA size (m)",ylab="Protected individuals",main=paste0(specname, "; ", locname , " abundance."), xlim=c(0,xmaxlimval),ylim=c(0,ymaxlim))
        for (i in 2:length(pMortIDs)){
          lines(MPAsizes,meansumn_survival[seq(count-mpas+1,count,1),pMortIDs[i]],lty=i)
        }
        legend("topleft", legend = legendtext, lty = seq(1,length(legendtext),1), title = "Fishing mortality")
        if (save.plots == 1) {dev.off()}

        # plot mean protection (proportion of range within MPA)
        sdlow = statstable[seq(count-mpas+1,count,1),10] - statstable[seq(count-mpas+1,count,1),11] #/sqrt(nreplicates)
        sdhigh = statstable[seq(count-mpas+1,count,1),10] + statstable[seq(count-mpas+1,count,1),11]#/sqrt(nreplicates)
        sdhigh[sdhigh > 1] = 1 # cap values at maximum
        sdlow[sdlow < 0] = 0 # cap values at minimum

        if (save.plots == 1) {
          png(paste0(savedir, speclatname, "_", locname , "_abundance_mean_range_in_mpa.png"))
        }
        plot(MPAsizes,statstable[seq(count-mpas+1,count,1),10],type="l",lty = 1, xlab="MPA size (m)",ylab="Protection in MPA (mean ? std)",main=paste0(specname, "; ", locname , " abundance."), xlim=c(0,xmaxlimval), ylim=c(0,1))
        #lines(MPAsizes,meanp05p_survival[seq(count-mpas+1,count,1),pMortIDs[length(pMortIDs)]],lty=2)
        #lines(MPAsizes,meanp95p_survival[seq(count-mpas+1,count,1),pMortIDs[length(pMortIDs)]],lty=2)
        lines(MPAsizes,sdlow,lty=2)
        lines(MPAsizes,sdhigh,lty=2)
        if (save.plots == 1) {dev.off()}

        # plot mortality
        ymaxlim = 1
        if (save.plots == 1) {
          png(paste0(savedir, speclatname, "_", locname , "_abundance_mean_Fmortality.png"))
        }
        plot(MPAsizes,meanpropn_mortality[seq(count-mpas+1,count,1),pMortIDs[1]],type="l",lty = 1, xlab="MPA size (m)",ylab="Lifetime fishing mortality",main=paste0(specname, "; ", locname , " abundance."), xlim=c(0,xmaxlimval), ylim=c(0,ymaxlim))
        for (i in 2:length(pMortIDs)){
          lines(MPAsizes,meanpropn_mortality[seq(count-mpas+1,count,1),pMortIDs[i]],lty=i)
        }
        legend("topright", legend = legendtext, lty = seq(1,length(legendtext),1), title = "Fishing mortality")
        if (save.plots == 1) {dev.off()}

        # plot mean probability of survival
        # if (save.plots == 1) {
        #    png(paste0(savedir, speclatname, "_", locname , "_abundance_mean_survival_probability.png"))
        #  }
        #  plot(MPAsizes,meanmeanp_survival[seq(count-mpas+1,count,1),pMortIDs[1]],type="b",lty = 1, xlab="MPA size (m)",ylab="Probability of survival",main=paste0(specname, "; ", locname , " abundance."), xlim=c(0,xmaxlimval), ylim=c(0,1))
        #  for (i in 2:length(pMortIDs)){
        #    lines(MPAsizes,meanmeanp_survival[seq(count-mpas+1,count,1),pMortIDs[i]],lty=i)
        #  }
        #  legend("bottomright", legend = legendtext, lty = seq(1,length(legendtext),1), title = "Probability of catch ")
        #  if (save.plots == 1) {dev.off()}


        # plot mortality offset
        if (save.plots == 1) {
          png(paste0(savedir, speclatname, "_", locname , "_abundance_mean_mortality_offset.png"))
        }
        plot(MPAsizes,mean_mortalityInst_offset[seq(count-mpas+1,count,1),pMortIDs[1]],type="l",lty = 1, xlab="MPA size (m)",ylab="Fishing mortality",main=paste0(specname, "; ", locname , " abundance."), xlim=c(0,xmaxlimval), ylim=c(0,max(pMortsInst)))
        for (i in 2:length(pMortIDs)){
          lines(MPAsizes,mean_mortalityInst_offset[seq(count-mpas+1,count,1),pMortIDs[i]],lty=i)
        }
        #legend("topright", legend = legendtext, lty = seq(1,length(legendtext),1), title = "Fishing morality ")
        if (save.plots == 1) {dev.off()}

      }

      # computing time status report
      mpatime <- Sys.time()
      print(paste0("Species ", s, "/" , nspecs, ", abundance ", lc, "/" , nlocations, ": ", "MPA size ", mpas, "/" , length(MPAsizes) , " completed."))
      print(mpatime - starttime)  
      
    } # mpa loop   
    
  } # location loop
  
  #} # species loop
  
  endtime <- Sys.time()
  print(endtime - starttime)
  
  if (save.results == 1){
    
    # convert results to data frames
    ioutput <- as.data.frame(ioutput)
    colnames(ioutput) = ioutputHeader
    
    statstable <- as.data.frame(statstable)
    colnames(statstable) = statstableHeader
    
    mean_survival_probability = as.data.frame(meanmeanp_survival, row.names = FALSE) 
    colnames(mean_survival_probability) = statsmatsHeader
    
    sd_survival_probability = as.data.frame(meansdp_survival, row.names = FALSE) 
    colnames(sd_survival_probability) = statsmatsHeader
    
    mean_protected_individuals = as.data.frame(meansumn_survival, row.names = FALSE)
    colnames(mean_protected_individuals) = statsmatsHeader
    
    median_protected_individuals = as.data.frame(mediansumn_survival, row.names = FALSE)
    colnames(median_protected_individuals) = statsmatsHeader
    
    sd_protected_individuals = as.data.frame(sdsumn_survival, row.names = FALSE)
    colnames(sd_protected_individuals) = statsmatsHeader
    
    mean_Fmortality = as.data.frame(meanpropn_mortality, row.names = FALSE)
    colnames(mean_Fmortality) = statsmatsHeader
    
    median_Fmortality = as.data.frame(medianpropn_mortality, row.names = FALSE)
    colnames(median_Fmortality) = statsmatsHeader
    
    sd_Fmortality = as.data.frame(sdpropn_mortality, row.names = FALSE)
    colnames(sd_Fmortality) = statsmatsHeader
    
    mean_Foffset = as.data.frame(mean_mortalityInst_offset, row.names = FALSE)
    colnames(mean_Foffset) = statsmatsHeader
    
    sd_Foffset = as.data.frame(sd_mortalityInst_offset, row.names = FALSE)
    colnames(sd_Foffset) = statsmatsHeader
    
    other.parameters <-  as.data.frame(t(c(nreplicates,years,maxd.per.ind, probs.per.ind)))
    colnames(other.parameters) <- c("num_replicates", "lifetime_years", "max_distance_per_individual","movement_probability_per_individual")
    
    write.table(other.parameters,paste0(savedir,"parameters.csv"), sep =",", row.names = FALSE)
    write.table(ioutput,paste0(savedir,"results_by_individual.csv"), sep =",", row.names = FALSE)
    write.table(statstable,paste0(savedir,"mean_results.csv"), sep =",", row.names = FALSE)
    write.table(mean_survival_probability,paste0(savedir,"mean_survival_probability.csv"), sep =",", row.names = FALSE)
    write.table(sd_survival_probability,paste0(savedir,"sd_survival_probability.csv"), sep =",", row.names = FALSE)
    write.table(mean_protected_individuals,paste0(savedir,"mean_protected_individuals.csv"), sep =",", row.names = FALSE)
    write.table(sd_protected_individuals,paste0(savedir,"sd_protected_individuals.csv"), sep =",", row.names = FALSE)
    #write.table(median_protected_individuals,paste0(savedir,"median_protected_individuals.csv"), sep =",", row.names = FALSE)
    write.table(mean_Fmortality,paste0(savedir,"mean_Fmortality.csv"), sep =",", row.names = FALSE)
    write.table(sd_Fmortality,paste0(savedir,"sd_Fmortality.csv"), sep =",", row.names = FALSE)
    #write.table(median_Fmortality,paste0(savedir,"median_Fmortality.csv"), sep =",", row.names = FALSE)
    write.table(mean_Foffset,paste0(savedir,"mean_Foffset.csv"), sep =",", row.names = FALSE)
    write.table(sd_Foffset,paste0(savedir,"sd_Foffset.csv"), sep =",", row.names = FALSE)
    
  } # save results
  #return(probsdflist)
  #return(do.call(rbind,probsdflist)) # Brings back list of movement probabilities for each shark
  
}# Function end


####