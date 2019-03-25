SPMSY_new<-function (x, Data, reps = 100) 
{
  dependencies = "Data@MaxAge, Data@vbK, Data@L50, Data@Cat, Data@Dep,Data@CV_Dep"
  
  ############## inverse von Bertalanffy ###########
  iVB <- function(t0, K, Linf, L) {
    max(1, ((-log(1 - L/Linf))/K + t0))  
  }
  ########################################################## 
 
  nsamp <- reps*1000  # large enough to make sure there is always enough r-K pairs to sample from 
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
    rsamp <- runif(nsamp, 0.6, 1.5)
  if (mean(rule) > 1.5 & mean(rule) < 2.5) 
    rsamp <- runif(nsamp, 0.2, 1)
  if (mean(rule) > 2.5 & mean(rule) < 3.5) 
    rsamp <- runif(nsamp, 0.05, 0.5)
  if (mean(rule) > 3.5) 
    rsamp <- runif(nsamp, 0.015, 0.1)
  ### setting values for carrying capacity ### 
  Ksamp <- runif(nsamp, mean(Data@Cat[x, ])/rsamp, (10 * mean(Data@Cat[x, 
                                                                       ]))/rsamp)
  
  nyears <- length(Data@Cat[x, ])
  B <- array(NA, dim = c(nsamp, nyears))
  
  ### rule for setting stock status at B0 ###
  if (Data@Cat[x,1] < (0.5 * max(Data@Cat[x, ]))) {
    B[, 1] <- Ksamp * runif(nsamp, 0.5, 0.9)
  }else {
    B[, 1] <- Ksamp * runif(nsamp, 0.3, 0.6)
  }
  
  ### rule for setting stock status the last year ##
  ## setting the LB and UB at the 1-99% of the desired prior
  
  Dist<-rbeta(20000,alphaconv(Data@Dep,Data@Dep*Data@CV_Dep),
              betaconv(Data@Dep,Data@Dep*Data@CV_Dep)) 
  Dist[is.na(Dist)]<-0   # removed nas created in the non-beta distributions
  LB<-quantile(Dist,0.01)
  UB<-quantile(Dist,0.99)
  
  
  for (i in 2:nyears) {
    B[, i] <- B[, i - 1] - Data@Cat[x,i - 1]
    B[, i] <- B[, i] + rsamp * B[, i] * (1 - B[, i]/Ksamp)
  }
  B <- B/rep(Ksamp, nyears) 
  cond <- (B[, nyears] >= LB) & (B[, nyears] <= UB) 
  if (sum(cond) < 1) {            
    B[B[, nyears] >= UB, nyears] <- UB 
    cond <- (B[, nyears] >= LB) & (B[, nyears] <= UB)
  }
  
  
  dep <- B[cond, nyears]
  MSY <- rsamp[cond] * Ksamp[cond]/4
  Kc <- Ksamp[cond]
  rc <- rsamp[cond]
  
  
  
  #### sampling from the retained trajectories in intervals of 0.05 to describe the desired shape of stock depletion 
  ##### addition by AC ##########
  
  prior<-rbeta(20000,alphaconv(Data@Dep,Data@Dep*Data@CV_Dep),
                betaconv(Data@Dep,Data@Dep*Data@CV_Dep)) 
  prior[is.na(prior)]<-0     # removes nas in the non-beta priors
  
   t<-0
  val_D<-list(1:20)
  val_K<-list(1:20)
  val_r<-list(1:20)
  count<-vector()
  prc<-vector()
 
  for (i in 1:20){
    
    count_exp <-prior>=t & prior<t+0.05
    count_mod <-dep>=t & dep<t+0.05    
    prc[i]<-round(sum(count_exp)) ## number of samples that fall within each interval from the desired prior
    count[i]<-sum(count_mod)      ## number of samples that fall within each interval from the retained final trajectories in dep

    vec<-vector(mode = "numeric", length = prc[i])
   
    ### if no density in the prior then 0 samples with me taken from dep
    if (prc[i]==0) {
      val_D[[i]]<-0
      val_K[[i]]<- 0
      val_r[[i]] <-0
    }else if (prc[i]<=count[i]){  ### if samples in the prior<=samples in dep then as many samples as in prior will me taken from dep
      vec<-sample(1:count[i],prc[i])
      val_D[[i]]<-dep[count_mod][vec]  
      val_K[[i]]<- Kc[count_mod][vec]   
      val_r[[i]]<- rc[count_mod][vec]
    } else if (prc[i]>count[i]){  ### if samples in the prior>samples in dep then all samples in dep will me retained
      
      val_D[[i]]<-dep[count_mod]   # Depletion
      val_K[[i]]<- Kc[count_mod]   # Carrying capacity
      val_r[[i]]<- rc[count_mod]   #Intrinsic growth rate
    }
    
    
    t<-t+0.05
    
  }
  
  
  ######### retained values################## 
  DW<-unlist(val_D)
  KcW<-unlist(val_K)
  rcW<-unlist(val_r)
  DW<-sample(DW,reps)
  KcW<-sample(KcW,reps)
  rcW<-sample(rcW,reps)
  ####weighted TAC
  TACW <- KcW * DW * rcW/2
  #############################################
  
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TACW)
  Rec
}

class(SPMSY_new)<-"MP"

sfExport('SPMSY_new')






