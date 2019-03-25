#### code adjusted from datalimited2 package.
##Free CM (2018) datalimited2: More stock assessment methods for data-limited fisheries. R package version 0.1.0.
### CMSY where the final accepted biomass trajectories that fall in the CI for the most likely values for 
### intrinsic rate r are sampled to describe the desired stock status prior

CMSY_new<-function (x, Data, reps = 100) { ## add this to the end
  
  dependencies = "Data@MaxAge, Data@vbK, Data@L50, Data@Cat, Data@Dep,Data@CV_Dep"
  
  ########### Helping functions ################################ 
  
  # ##### sampling the retained biomass trajectories so dep fits the desired prior  ######
  ##### addition by AC ##########
  Weights<-function(Dep,CV_Dep,dep_last,kv.all,rv.all) {
    
    exp_pr<-rbeta(1000,alphaconv(Data@Dep,Data@Dep*Data@CV_Dep),
                  betaconv(Data@Dep,Data@Dep*Data@CV_Dep)) 
    exp_pr[is.na(exp_pr)]<-0
    ##### parameters for the intervals #######
    
    t<-0
    
    ##########################  
    val_D<-list(1:20)
    val_K<-list(1:20)
    val_r<-list(1:20)
    count<-vector()
    prc<-vector()
    for (i in 1:20){
      
      count_exp <-exp_pr>=t & exp_pr<t+0.05
      count_mod <-dep_last>=t & dep_last<t+0.05    ######### dep is the depletion from the model for all values between the predefined boundaries
      prc[i]<-round(sum(count_exp))        ## number of samples that fall within each interval from the desired prior
      count[i]<-sum(count_mod)             ## number of samples that fall within each interval from the retained final trajectories in dep
      #print(paste("percent",prc[i]))
      #print(paste("samples",count[i]))
      vec<-vector(mode = "numeric", length = prc[i])
      
          ### if no density in the prior then 0 samples with me taken from dep
      if (prc[i]==0) {
        val_D[[i]]<-0
        val_K[[i]]<- 0
        val_r[[i]] <-0
      }else if (prc[i]<=count[i]){         ### if samples in the prior<=samples in dep then as many samples as in prior will me taken from dep
        vec<-sample(1:count[i],prc[i])
        val_D[[i]]<-dep_last[count_mod][vec]  
        val_K[[i]]<- kv.all[count_mod][vec]   
        val_r[[i]]<- rv.all[count_mod][vec]
      } else if (prc[i]>count[i]){        ### if samples in the prior>samples in dep then all samples in dep will me retained
        
        val_D[[i]]<-dep_last[count_mod]   # Depletion
        val_K[[i]]<- kv.all[count_mod]   # Carrying capacity
        val_r[[i]]<- rv.all[count_mod]   #Intrinsic growth rate
      }
      
      
      t<-t+0.05
      
    }
    
    ######### retained values################## 
    DW<-unlist(val_D)
    DW<-DW[DW!=0]
    KW<-unlist(val_K)
    KW<-KW[KW!=0]
    rW<-unlist(val_r)
    rW<-rW[rW!=0]
    All_W<-cbind(rW,KW,DW)
    colnames(All_W)<-c("rW","KW","DW")
    return(All_W)
  }
  
  # SchaeferParallelSearch
  ################################################################################
  
  # Monte Carlo filtering with Schaefer Function
  `%dopar%` <- foreach::`%dopar%`
  SchaeferParallelSearch <- function(ni, nyr, sigR, duncert, ct, int.yr, intbio, startbt, ki, i, ri, int.yr.i, nstartbt, yr, end.yr, endbio, npoints, pt){
    ptm <- proc.time()
    # create vectors for viable r, k and bt
    inmemorytable <- vector()
    # parallelised for the points in the r-k space
    inmemorytable <- foreach::foreach(i = 1 : npoints, .combine='rbind', .packages='foreach', .inorder=TRUE) %dopar%{
      nsbt = length(startbt)
      VP   <- FALSE
      for(nj in 1:nsbt) {
        # create empty vector for annual biomasses
        bt <- vector()
        j<-startbt[nj]
        # set initial biomass, including 0.1 process error to stay within bounds
        bt[1]=j*ki[i]*exp(stats::rnorm(1,0, 0.1*sigR))  ## set biomass in first year
        # repeat test of r-k-startbt combination to allow for different random error
        for(re in 1:ni)   {
          #loop through years in catch time series
          for (t in 1:nyr)  {  # for all years in the time series
            xt=stats::rnorm(1,0, sigR) # set new process error for every year
            zlog.sd = sqrt(log(1+(duncert)^2))
            zt=stats::rlnorm(1,meanlog = 0, sdlog = zlog.sd) # model the catch error as a log normal distribution.
            # calculate biomass as function of previous year's biomass plus surplus production minus catch
            bt[t+1] <- ifelse(bt[t]/ki[i] >= 0.25,
                              bt[t]+ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt,
                              bt[t]+(4*bt[t]/ki[i])*ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt) # assuming reduced r at B/k < 0.25
            # if biomass < 0.01 k, discard r-k-startbt combination
            if(bt[t+1] < 0.01*ki[i]){
              break
            } # stop looping through years, go to next upper level
            # intermediate year check
            if((t+1)==int.yr.i && (bt[t+1]>(intbio[2]*ki[i]) || bt[t+1]<(intbio[1]*ki[i]))){
              break
            }
          } # end of loop of years
          # if loop was broken or last biomass falls outside of expected ranges
          # do not store results, go directly to next startbt
          if(t < nyr || bt[yr==end.yr] > (endbio[2]*ki[i]) || bt[yr==end.yr] < (endbio[1]*ki[i])){
            next
          }else{
            #each vector will be finally appended to the others found by the threads - this is done by the .combine='rbind' option
            inmemorytablerow<-c(i,j,ri[i],ki[i],bt[1:(nyr+1)]/ki[i])
            if(length(inmemorytablerow)==(4+nyr+1)){
              if(VP==FALSE){
                inmemorytable <- inmemorytablerow
              }else{
                inmemorytable <- rbind(inmemorytable,inmemorytablerow)
              }
              VP <- TRUE
            }
          }
        } # end of repetition for random error
      } # end of j-loop of initial biomasses
      # instruction necessary to make the foreach loop see the variable:
      if(length(inmemorytable)==0){
        inmemorytable <- vector(length=4+nyr+1)*NA
      }else{
        inmemorytable
      }
    }#end loop on points
    
    #create the output matrix
    mdat        <- matrix(data=NA, nrow = npoints*nstartbt, ncol = 2+nyr+1)
    npointsmem = dim(inmemorytable)[1]
    npointscols = dim(inmemorytable)[2]
    #reconstruction of the processing matrix after the parallel search
    if (npointsmem>0 && npointscols>0){
      for (idxr in 1:npointsmem){
        i = inmemorytable[idxr,1]
        if (!is.na(i)){
          j = inmemorytable[idxr,2]
          mdatindex<-((i-1)*nstartbt)+which(startbt==j)
          mdat[mdatindex,1]           <- inmemorytable[idxr,3]
          mdat[mdatindex,2]           <- inmemorytable[idxr,4]
          mdat[mdatindex,3:(2+nyr+1)] <- inmemorytable[idxr,5:(4+nyr+1)]
          # if(pt==T) points(x=ri[i], y=ki[i], pch=".", cex=4, col="gray")
        }
      }
    }
    ptm<-proc.time()-ptm
    mdat <- stats::na.omit(mdat)
    return(mdat)
  }
  
  # SchaeferMC
  ################################################################################
  
  SchaeferMC <- function(ri, ki, startbio, int.yr, intbio, endbio, sigR, pt, duncert, startbins, ni, yr, nyr, ct, end.yr, verbose){
    # create vector for initial biomasses
    startbt <- seq(from =startbio[1], to=startbio[2], by = (startbio[2]-startbio[1])/startbins)
    nstartbt <- length(startbt)
    npoints <- length(ri)
    # get index of intermediate year
    int.yr.i <- which(yr==int.yr)
    #loop through r-k pairs with parallel search
    mdat<-SchaeferParallelSearch(ni, nyr,sigR,duncert,ct,int.yr,intbio, startbt, ki, i, ri, int.yr.i, nstartbt, yr, end.yr, endbio, npoints,pt)
    if(verbose==T){cat("\n")}
    return(list(mdat))
  }
  
  # Moving average function
  ################################################################################
  
  # Calculate moving average
  ma <- function(x){
    x.1 <- stats::filter(x,rep(1/3,3),sides=1)
    x.1[1] <- x[1]
    x.1[2] <- (x[1]+x[2])/2
    return(x.1)
  }
  
  # Functions for setting priors
  ################################################################################
  
  # Set R prior
  r_prior <- function(r.low, r.hi){
    # initial range of r from input file
    
    start.r <- c(r.low,r.hi)
    return(start.r)
  }
  
  # Set start saturation prior
  startbio_prior <- function(stb.low, stb.hi, start.yr){
    # use initial biomass range from input file if stated
    if(is.na(stb.low)==F & is.na(stb.hi)==F){
      startbio <- c(stb.low,stb.hi)
    }else{
      # if start year < 1960 assume high biomass
      if(start.yr < 1960){
        startbio <- c(0.5,0.9)
      }else{
        # else use medium prior biomass range
        startbio <- c(0.2,0.6)
      }
    }
    return(startbio)
  }
  
  # Set intermediate saturation prior
  intbio_prior <- function(intb.low, intb.hi, int.yr, start.yr, end.yr, startbio, yr, ct){
    # get index of years with lowest and highest catch between start+3 and end-3 years
    min.yr.i <- which.min(ct[4:(length(ct)-3)])+3
    max.yr.i <- which.max(ct[4:(length(ct)-3)])+3
    min.ct <- ct[min.yr.i]
    max.ct <- ct[max.yr.i]
    # use year and biomass range for intermediate biomass from input file
    if(is.na(intb.low)==F & is.na(intb.hi)==F){
      int.yr   <- int.yr
      intbio   <- c(intb.low,intb.hi)
      # if contrast in catch is low, use initial range again in mid-year
    }else if(min(ct)/max(ct) > 0.6) {
      int.yr    <- as.integer(mean(c(start.yr, end.yr)))
      intbio    <- startbio
      # else if year of minimum catch is after max catch then use min catch
    }else if(min.yr.i > max.yr.i) {
      int.yr    <- yr[min.yr.i-1]
      if(startbio[1]>=0.5 &  (int.yr-start.yr) < (end.yr-int.yr) &
         (min.ct/max.ct) > 0.3) intbio <- c(0.2,0.6) else intbio <- c(0.01,0.4)
         # else use max catch
    } else {
      # assume that biomass range in year before maximum catch was high or medium
      int.yr    <- yr[max.yr.i-1]
      intbio    <- if((startbio[1]>=0.5 & (int.yr-start.yr) < (end.yr-int.yr))| # if initial biomass is high, assume same for intermediate
                      # ((min.ct/max.ct < 0.3 & (max.yr.i - min.yr.i) < 25))) c(0.5,0.9) else c(0.2,0.6) }
                      (((max.ct-min.ct)/max.ct)/(max.yr.i-min.yr.i) > 0.04)) c(0.5,0.9) else c(0.2,0.6) } # if incease is steep, assume high, else medium
    out <- list(intbio, int.yr)
    return(out)
  }
  
  # Set K prior
  k_prior <- function(endbio, start.r, ct){
    # initial prior range of k values, assuming min k will be larger than max catch / prior for r
    if(mean(endbio) <= 0.5){
      start.k <- c(max(ct)/start.r[2],4*max(ct)/start.r[1])
    }else{
      start.k <- c(2*max(ct)/start.r[2],10*max(ct)/start.r[2])
    }
    return(start.k)
  }
  
  ########## function for age at maturity from DLMtool ###########
  iVB <- function(t0, K, Linf, L) {
    max(1, ((-log(1 - L/Linf))/K + t0))  # Inverse Von-B
  }
  ########################################################## 
  
  # CMSY catch-only stock assessment model from Froese et al. 2017.
   
  verbose <- T ### 
  ################## rul3 for initial value for r #########
  
 ##### intrinsic growth rate is selected as in SPMSY on DLMtool and not qualitative resielence estimates###
  ### added by AC ###
  
  rule <- rep(4, 3)
  if (Data@vbK[x] > 0.3) {
    rule[1] <- 1
  }else if (Data@vbK[x] < 0.3 & Data@vbK[x] > 0.16) {
    rule[1] <- 2
  }else if (Data@vbK[x] < 0.16 & Data@vbK[x] > 0.05) {
    rule[1] <- 3
  }
  
  AM <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], Data@L50[x])
  if (AM < 1.5) {
    rule[2] <- 1
  }else if (AM < 4.5 & AM > 1.5) {
    rule[2] <- 2
  }else if (AM < 10 & AM > 4.5) {
    rule[2] <- 3
  }
  if (Data@MaxAge < 4) {
    rule[3] <- 1
  }else if (Data@MaxAge < 11 & Data@MaxAge > 3) {
    rule[3] <- 2
  }else if (Data@MaxAge < 31 & Data@MaxAge > 10) {
    rule[3] <- 3
  }
  #### setting values for intrinsic rate r###
  if (mean(rule) < 1.5) 
    r<-c(0.6,1.5)
  
  if (mean(rule) > 1.5 & mean(rule) < 2.5) 
    r<-c(0.2,1)
  
  if (mean(rule) > 2.5 & mean(rule) < 3.5) 
    r<-c(0.05,0.5)
  
  if (mean(rule) > 3.5) 
    r<-c(0.015,0.1) 
#####################################################  
  
  stb.low <- NA
  stb.hi <- NA
  int.yr <- NA 
  intb.low<-NA
  intb.hi <- NA 
  
  #### the desired stock status priors used to define the lower abd upper bounds in the model ##
  Dist<-rbeta(10000,alphaconv(Data@Dep,Data@Dep*Data@CV_Dep),
              betaconv(Data@Dep,Data@Dep*Data@CV_Dep)) 
  Dist[is.na(Dist)]<-0   ## removes nas from non-beta distributions
  endb.low<-quantile(Dist,0.01)
  endb.hi<-quantile(Dist,0.99)
  verbose <- T
  
  year<-Data@Year
  catch<-as.vector(Data@Cat[x,])
  
  # Display 3 digits
  options(digits=3)
  
  # Perform a few error checks
  #  res_vals <- c("High", "Medium", "Low", "Very low")
  #  if(sum(is.na(catch))>0){stop("Error: NA in catch time series. Fill or interpolate.")}
  #  if(is.na(resilience) & (is.na(r.low) | is.na(r.hi))){stop("Error: Either a resilience estimate or a complete user-specified r prior (both r.low and r.hi) must be provided.")}
  #  if(!is.na(resilience) & !(resilience%in%res_vals)){stop("Error: Resilience must be 'High', 'Medium', 'Low', or 'Very low' (case-sensitive).")}
  
  # Set model parameters
  #############################################################
  
  # Setup parallel processing
  # Use 3 chains in JAGS if more than 2 cores are available
  n.cores <- pmin(2, parallel::detectCores())
  n.chains <- ifelse(n.cores > 2,3,2)
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl, cores = n.cores)
  
  # Set model parameters
  FullSchaefer <- F # will automatically change to TRUE if enough abundance data available
  dataUncert <- 0.1  # set observation error as uncertainty in catch - default is SD=0.1
  sigmaR <- 0.1 # overall process error for CMSY; SD=0.1 is the default
  n <- 10*reps #0 # initial number of r-k pairs
  n.new <- n # initialize n.new
  ni <- 3 # iterations for r-k-startbiomass combinations, to test different variability patterns; no improvement seen above 3
  nab <- 5 # default=5; minimum number of years with abundance data to run BSM
  duncert <- dataUncert # global defaults for uncertainty
  sigR <- sigmaR # global defaults for uncertainty
  
  # Setup data
  #############################################################
  
  # Build catch data 
  catchData <- data.frame(yr=year, ct=catch)
  
  # Transform catch data
  # 1. Convert to 1000s tons (or other units)
  # 2. Calculate 3-yr moving average (average of past 3 years)
  ct.raw <- catchData$ct / 1000 
  ct <- ma(ct.raw)
  
  # Identify number of years and start/end years
  yr <- catchData$yr # functions use this quantity
  nyr <- length(yr)
  start.yr <- min(yr)
  end.yr <- max(yr)
  
  # Determine initial ranges for parameters and biomass
  #############################################################
  
  # Set priors
  start.r <- r_prior(r[1], r[2])
  startbio <- startbio_prior(stb.low, stb.hi, start.yr)
  int_params <- intbio_prior(intb.low, intb.hi, int.yr, start.yr, end.yr, startbio, yr, ct)
  intbio <- c(0.01,1)  ### intbio with no limitations, changed from AC
  int.yr <- int_params[[2]] 
  endbio <- c(endb.low,endb.hi)
  start.k <- k_prior(endbio, start.r, ct)
  
  # Record priors into dataframe
  priors <- data.frame(cbind(c("r", "k", "startbio", "intbio", "endbio"), source="default",
                             rbind(start.r, start.k, startbio, intbio, endbio)), year=NA, stringsAsFactors=F)
  colnames(priors) <- c("param", "source", "lo", "hi", "year")
  rownames(priors) <- NULL
  priors$year[priors$param=="intbio"] <- int.yr
  priors$lo <- as.numeric(priors$lo)
  priors$hi <- as.numeric(priors$hi)
  if(!is.na(r[1])){priors$source[priors$param=="r"] <- "expert"}
  if(!is.na(stb.low)){priors$source[priors$param=="startbio"] <- "expert"}
  if(!is.na(intb.low)){priors$source[priors$param=="intbio"] <- "expert"}
  if(!is.na(endb.low)){priors$source[priors$param=="endbio"] <- "expert"}
  
  # Print priors (if desired)
  if(verbose==T){
    cat("startbio=",startbio,ifelse(is.na(stb.low)==T,"default","expert"),
        ", intbio=",int.yr,intbio,ifelse(is.na(intb.low)==T,"default","expert"),
        ", endbio=",endbio,ifelse(is.na(endb.low)==T,"default","expert"),"\n")
  }
  
  # Monte Carlo procedure
  #############################################################
  
  # Initialize other vectors anew for each stock
  current.attempts <- NA
  
  
  
  # Get random set of r and k from log space distribution
  ri1 = exp(stats::runif(n, log(start.r[1]), log(start.r[2])))
  ki1 = exp(stats::runif(n, log(start.k[1]), log(start.k[2])))
  
  # 1 - Call CMSY-SchaeferMC function to preliminary explore the r-k space
  if(verbose==T){cat("First Monte Carlo filtering of r-k space with ",n," points...\n")}
  MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR,
                     pt=T, duncert=dataUncert, startbins=10, ni=ni, yr=yr, nyr=nyr, ct=ct, end.yr=end.yr, verbose=verbose)
  #mdat.all <- rbind(mdat.all,MCA[[1]])
  rv.all   <- MCA[[1]][,1]
  kv.all   <- MCA[[1]][,2]
  btv.all  <-MCA[[1]][,3:(2+nyr+1)]
  dep_last<-btv.all[,length(year)] ### depletion in the last year
  
  # Initialize vectors for viable r, k, bt, and all in a list
  n.viable.b   <- length(rv.all)
  n.viable.pt <- length(unique(rv.all))
  W.all <- vector("list",3) ### saves all the retained trajectories
  W.all[[1]]<-rv.all
  W.all[[2]]<-kv.all
  W.all[[3]]<-dep_last
  
  
  if(verbose==T){cat("Found ",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")}
  
  # 2 - if the lower bound of k is too high, reduce it by half and rerun
  if(length(kv.all[kv.all < 1.1*start.k[1] & rv.all < mean(start.r)]) > 10) {
    if(verbose==T){cat("Reducing lower bound of k, resampling area with",n,"additional points...\n")}
    start.k <- c(0.5*start.k[1],start.k[2])
    ri1 = exp(stats::runif(n, log(start.r[1]), log(start.r[2])))
    ki1 = exp(stats::runif(n, log(start.k[1]), log(start.k[2])))
    MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR,
                       pt=T, duncert=dataUncert, startbins=10, ni=ni, yr=yr, nyr=nyr, ct=ct, end.yr=end.yr, verbose=verbose)
    rv.all   <- MCA[[1]][,1]
    kv.all   <- MCA[[1]][,2]
    btv.all  <-MCA[[1]][,3:(2+nyr+1)]
    dep_last<-btv.all[,length(year)] ### depletion in the last year
    
    # Initialize vectors for viable r, k, bt, and all in a list
    
    W.all[[1]]<-c(W.all[[1]],rv.all)
    W.all[[2]]<-c(W.all[[2]],kv.all)
    W.all[[3]]<-c(W.all[[3]],dep_last)
    n.viable.b   <- length(rv.all)
    n.viable.pt <- length(unique(rv.all))
    
    if(verbose==T){cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")}
  }
  
  # 3 - if few points were found then resample and shrink the log k space
  if (n.viable.b < 1000){
    log.start.k.new  <- log(start.k)
    max.attempts     <- 3
    current.attempts <- 1
    startbins        <- 10
    while (n.viable.b <= 1000 && current.attempts <= max.attempts){
      if(n.viable.pt > 0) {
        log.start.k.new[1] <- mean(c(log(start.k[1]), min(log(kv.all))))
        log.start.k.new[2] <- mean(c(log.start.k.new[2], max(log(kv.all)))) }
      n.new <- n*current.attempts #add more points
      ri1 = exp(stats::runif(n.new, log(start.r[1]), log(start.r[2])))
      ki1 = exp(stats::runif(n.new, log.start.k.new[1], log.start.k.new[2]))
      if(verbose==T){cat("Shrinking k space: repeating Monte Carlo in the interval [",exp(log.start.k.new[1]),",",exp(log.start.k.new[2]),"]\n")}
      if(verbose==T){cat("Attempt ",current.attempts," of ",max.attempts," with ",n.new," additional points...","\n")}
      if(current.attempts==2 & n.viable.b < 50){
        duncert   <- 2*dataUncert
        sigR      <- 2*sigmaR
        startbins <- 20
        if(verbose==T){cat("Doubling startbins, catch and process error, and number of variability patterns \n")}
      }
      MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR,
                         pt=T, duncert=duncert, startbins=startbins, ni=2*ni, yr=yr, nyr=nyr, ct=ct, end.yr=end.yr, verbose=verbose)
      rv.all   <- MCA[[1]][,1]
      kv.all   <- MCA[[1]][,2]
      btv.all  <-MCA[[1]][,3:(2+nyr+1)]
      dep_last<-btv.all[,length(year)] ### depletion in the last year
      
        # Initialize vectors for viable r, k, bt, and all in a list
      
      W.all[[1]]<-c(W.all[[1]],rv.all)
      W.all[[2]]<-c(W.all[[2]],kv.all)
      W.all[[3]]<-c(W.all[[3]],dep_last)
      n.viable.b   <- length(rv.all)
      n.viable.pt <- length(unique(rv.all))
      
      if(verbose==T){cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")}
      current.attempts=current.attempts+1 #increment the number of attempts
    }
    if(n.viable.b < 5) {
      if(verbose==T){cat("Only",n.viable.pt,"viable r-k pairs found, check data and settings \n")}
      next
    }
  }
  
  # 4 - if tip of viable r-k pairs is 'thin', do extra sampling there
  if(length(rv.all[rv.all > 0.9*start.r[2]]) < 5) {
    l.sample.r        <- stats::quantile(rv.all,0.6)
    add.points        <- ifelse(is.na(current.attempts)==T,n,ifelse(current.attempts==2,2*n,ifelse(length(rv.all)>500,3*n,6*n)))
    if(verbose==T){cat("Final sampling in the tip area above r =",l.sample.r,"with",add.points,"additional points...\n")}
    log.start.k.new <- c(log(0.8*min(kv.all)),log(max(kv.all[rv.all > l.sample.r])))
    
    ri1 = exp(stats::runif(add.points, log(l.sample.r), log(start.r[2])))
    ki1 = exp(stats::runif(add.points, log.start.k.new[1], log.start.k.new[2]))
    MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR,
                       pt=T, duncert=duncert, startbins=10, ni=ni, yr=yr, nyr=nyr, ct=ct, end.yr=end.yr, verbose=verbose)
    rv.all   <- MCA[[1]][,1]
    kv.all   <- MCA[[1]][,2]
    btv.all  <-MCA[[1]][,3:(2+nyr+1)]
    dep_last<-btv.all[,length(year)] ### depletion in the last year
      
    W.all[[1]]<-c(W.all[[1]],rv.all)
    W.all[[2]]<-c(W.all[[2]],kv.all)
    W.all[[3]]<-c(W.all[[3]],dep_last)
    n.viable.b   <- length(rv.all)
    n.viable.pt <- length(unique(rv.all))
    
    if(verbose==T){cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")}
  }
  
  ############ weight the results based on expert prior ################
  
  
  
  
  # Extract model results
  #############################################################
  
  # get estimate of most probable r as 75th percentile of mid log.r-classes
  # get unique combinations of r-k
  
  rkB<-matrix(unlist(W.all), ncol = 3, byrow = FALSE)
  
  unique.rk         <- unique(rkB)
  # get remaining viable log.r and log.k
  log.rs           <- log(unique.rk[,1])
  log.ks           <- log(unique.rk[,2])
  # get vectors with numbers of r and mid values in classes
  # determine number of classes as a function of r-width
  r.width         <- (max(unique.rk[,1])-start.r[1])/(start.r[2]-start.r[1])
  classes         <- ifelse(r.width>0.8,100,ifelse(r.width>0.5,50,ifelse(r.width>0.3,25,12)))
  hist.log.r      <- graphics::hist(x=log.rs, breaks=classes, plot=F)
  log.r.counts    <- hist.log.r$counts
  log.r.mids      <- hist.log.r$mids
  # get most probable log.r as 75th percentile of mids with counts > 0
  log.r.est       <- as.numeric(stats::quantile(log.r.mids[which(log.r.counts > 0)],0.75))
  median.log.r    <- as.numeric(stats::quantile(x=log.r.mids[which(log.r.counts > 0)], 0.50))
  lcl.log.r       <- as.numeric(stats::quantile(x=log.r.mids[which(log.r.counts > 0)], 0.5125))
  ucl.log.r       <- as.numeric(stats::quantile(x=log.r.mids[which(log.r.counts > 0)], 0.9875))
  sd.log.r.est    <- (ucl.log.r - log.r.est) / 1.96
  r.est           <- exp(log.r.est)
  lcl.r.est       <- exp(log.r.est-1.96*sd.log.r.est)
  ucl.r.est       <- exp(log.r.est+1.96*sd.log.r.est)
  
  
  # get r-k pairs above median of mids
  rem            <- which(unique.rk[,1] > exp(median.log.r))
  rem.log.r      <- log(unique.rk[,1][rem])
  rem.log.k      <- log(unique.rk[,2][rem])
  # do linear regression of log k ~ log r with slope fixed to -1 (from Schaefer)
  reg            <- stats::lm(rem.log.k ~ 1 + offset(-1*rem.log.r))
  int.reg        <- as.numeric(reg[1])
  sd.reg      <- stats::sd(stats::resid(reg))
  # get estimate of log(k) from y where x = log.r.est
  log.k.est      <- int.reg + (-1) * log.r.est
  # get estimates of ucl of log.k.est from y + SD where x = ucl.log.r
  ucl.log.k     <- int.reg + (-1) * lcl.log.r + sd.reg
  # get estimates of sd.log.k.est from upper confidence limit of log.k.est
  sd.log.k.est   <- (ucl.log.k - log.k.est) / 1.96
  lcl.log.k      <- log.k.est - 1.96*sd.log.k.est
  ucl.log.k      <- log.k.est + 1.96*sd.log.k.est
  k.est       <- exp(log.k.est)
  lcl.k.est   <- exp(lcl.log.k)
  ucl.k.est   <- exp(ucl.log.k)
  
  
  # get MSY from remaining log r-k pairs
  log.MSY.est     <- mean(rem.log.r + rem.log.k - log(4))
  sd.log.MSY.est  <- stats::sd(rem.log.r + rem.log.k - log(4))
  lcl.log.MSY.est <- log.MSY.est - 1.96*sd.log.MSY.est
  ucl.log.MSY.est <- log.MSY.est + 1.96*sd.log.MSY.est
  MSY.est         <- exp(log.MSY.est)
  lcl.MSY.est     <- exp(lcl.log.MSY.est)
  ucl.MSY.est     <- exp(ucl.log.MSY.est)
  
  ########### 

  ### estimation of OFL using only biomass trajectories and r-k pairs within the confidence limits ###
  ####### added by AC #####
   
  btv.all<-rkB[which(rkB[,1] > lcl.r.est & rkB[,1] < ucl.r.est
                     & rkB[,2] > lcl.k.est & rkB[,2] < ucl.k.est),3]
  r.all<-rkB[which(rkB[,1] > lcl.r.est & rkB[,1] < ucl.r.est
                   & rkB[,2] > lcl.k.est & rkB[,2] < ucl.k.est),1]
  K.all<-rkB[which(rkB[,1] > lcl.r.est & rkB[,1] < ucl.r.est
                   & rkB[,2] > lcl.k.est & rkB[,2] < ucl.k.est),2]
  
  ### sampling the retained values to describe the desires stock status prior ###
  All_W <-Weights(Data@Dep,Data@CV_Dep,btv.all,K.all,r.all)
  rv.all   <- All_W[,1]
  kv.all   <- All_W[,2]
  Bt_K     <- All_W[,3]
  
### sampling equal to reps to estimate OFL ##
  mean_K<-kv.all
  sd_K<-density(kv.all)$bw
  K_D<-rnorm(reps,mean_K,sd_K)
  K_D[is.na(K_D)]<-0
  mean_r<-rv.all
  sd_r<-density(rv.all)$bw
  r_D<-rnorm(reps,mean_r,sd_r)
  r_D[is.na(r_D)]<-0
  mean_Bt_K<-Bt_K
  sd_Bt_K<-density(Bt_K)$bw
  Bt_K_D<-rbeta(reps,alphaconv(mean_Bt_K,sd_Bt_K),betaconv(mean_Bt_K,sd_Bt_K))
  Bt_K_D[is.na(Bt_K_D)]<-0
  TAC<-K_D*(r_D/2)*Bt_K_D*1000 ### turning to iriginal scale ###
 ######################################################################### 
  
  
  Rec <- new("Rec") 
  Rec@TAC <- TACfilter(TAC) 
  Rec  
  
} 

class(CMSY_new)<-"MP"

sfExport('CMSY_new') 

