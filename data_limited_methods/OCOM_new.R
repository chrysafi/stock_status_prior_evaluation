#### code adjusted from datalimited2 package.
##Free CM (2018) datalimited2: More stock assessment methods for data-limited fisheries. R package version 0.1.0. 
## OCOM where only change is the prior for stock status. Instead of being estimated from BRTs from Zhou et al. 2017
## the desired beta prior is the input

OCOM_new<-function (x,Data,reps){ 
  dependencies="Data@Year,Data@Cat, Data@Mort, Data@FMSY_M"
  ####helping function #################### 
  # Biomass dynamics model
  BDM = function(K, r, S, b, C) {
    nyr = length(C)
    B = vector()
    B[1] = K*b
    for (i in 1:nyr) {
      B[i+1] = max(min(B[i]+r*B[i]*(1-B[i]/K)-C[i], K), 0)
    }
    if (all(B[-nyr]>C) & all(B<=K)) abs(B[nyr]/K-S) else max(K)*10^4
  }
  ##########################################
  
  
  if (sum(is.na(Data@Cat[x,])) > 0) {
    stop("Error: NA in catch time series. Fill or interpolate.")
  }
  nsim = reps*10
  summ = matrix(NA, 5, 4)
  C = Data@Cat[x,]
  yr = Data@Year
  M = Data@Mort[x]
  nyr = length(yr)
  FMSYM = Data@FMSY_M[x]
  
  ## here the original rule for setting r is kept but could also be changed to be the same as in SPMSY and CMSY ##
  r_median = 2 * (FMSYM) * M
  r_sig2 = (2 * FMSYM)^2 * (0.0012 + 0.23)
  r_sig = sqrt(r_sig2)
  r_mu = log(r_median)
  ri = stats::rlnorm(nsim, r_mu, r_sig) 
  #########stock status prior #######
  ### Based on the desired prior not BRT as in the original model ###
  ## changed by AC ###
  
  si = rbeta(nsim, alphaconv(Data@Dep, Data@Dep*Data@CV_Dep), betaconv(Data@Dep, Data@Dep*Data@CV_Dep))  ### not using [x]for the non-beta distributions
  si[is.na(si)]<-0  ## removing nas produced for the non-beta distributions
  rs = cbind(r = ri, s = si)
  k.low = max(C)
  k.up = max(C) * 200
  opt = apply(rs, 1, function(x) {
    stats::optimize(BDM, c(k.low, k.up), r = x[["r"]], S = x[["s"]], #### need this BDM function #####
                    b = 1, C = C)
  })
  ki = sapply(opt, "[[", 1)
  msy = ki * ri/4
  obji = sapply(opt, "[[", 2)
  kr = data.frame(k = ki, r = ri, msy, s = si, obj = obji)
  kr2 = kr
  kr2[kr2$k < 1.01 * k.low | kr2$k > 0.99 * k.up | kr2$obj > 
        0.01, ] = NA
  kr2 = stats::na.omit(kr2)
  summ <- apply(kr2, 2, function(x) stats::quantile(x, c(0.025, 
                                                         0.25, 0.5, 0.75, 0.975)))[, 1:4]
  B = Bmed = vector()
  ntrajs = 1000
  oc0 <- kr
  oc1 = oc0[oc0$obj < 0.01 & oc0$k > 1.01 * min(oc0$k) & oc0$k < 
              0.99 * max(oc0$k), ]
  oc2 = oc1[oc1$k > stats::quantile(oc1$k, 0.25) & oc1$k < 
              stats::quantile(oc1$k, 0.75), ]
  smp = sample(1:nrow(oc2), ntrajs, replace = T)
  k = oc2[smp, 1]
  r = oc2[smp, 2]
  kmed = stats::median(k)
  rmed = stats::median(r)
  B <- matrix(NA, nrow = nyr, ncol = ntrajs, dimnames = list(yr, 
                                                             1:ntrajs))
  for (j in 1:ntrajs) {
    B[1, j] = k[j]
    for (t in 1:(nyr - 1)) {
      B[t + 1, j] = B[t, j] + r[j] * B[t, j] * (1 - B[t, 
                                                      j]/k[j]) - C[t]
    }
  }
  
  b_trajs <- B ### Biomass trajectory for all years
 
  # Saturation time series 
  s_trajs <- apply(b_trajs, 2, function(x) x/x[1]) 
  S_last<-s_trajs[length(yr),] ### depletion in the last year of the time series
  
  TAC<-k*S_last*r/2
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
}

class(OCOM_new)<-"MP"

sfExport('OCOM_new')
