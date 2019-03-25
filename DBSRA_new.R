
DBSRA_new<-function (x, Data, reps = 100) 
{
  dependencies = "Data@Cat, Data@Dep, Data@CV_Dep, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M,Data@BMSY_B0, Data@CV_BMSY_B0, Data@L50"
  
  ##### DBSRA supporting functions #################################
  iVB <- function(t0, K, Linf, L) {
    max(1, ((-log(1 - L/Linf))/K + t0))  # Inverse Von-B
  }
  
  fn <- function(n, BMSY_K) {
    # optimizer to find parameter n according to sampled BMSY/B0 (theta)
    thetapred <- n^(-1/(n - 1))
    (BMSY_K - thetapred)^2
  }
  
  getn <- function(BMSY_K) {
    # wrapper for n finder
    optimize(fn, c(0.01, 6), BMSY_K = BMSY_K)$minimum  #get the optimum
  }
  
  gety <- function(n) (n^(n/(n - 1)))/(n - 1)  # More DBSRA code: get the y parameter for n
  
  prodPTF <- function(depletion, n, MSY) {
    # Pella-Tomlinson production function required for DB-SRA
    y <- (n^(n/(n - 1)))/(n - 1)
    MSY * y * depletion - MSY * y * depletion^n
  }
  
  DBSRAopt <- function(lnK, C_hist, nys, Mdb, FMSY_M, BMSY_K, Bt_K, adelay) {
    # the optimization for B0 given DBSRA assumptions
    Kc <- exp(lnK)
    n <- getn(BMSY_K)
    g <- gety(n)
    FMSY <- FMSY_M * Mdb
    UMSY <- (FMSY/(FMSY + Mdb)) * (1 - exp(-(FMSY + Mdb)))
    MSY <- Kc * BMSY_K * UMSY
    # Bjoin rules from Dick & MacCall 2011
    Bjoin_K <- 0.5
    if (BMSY_K < 0.3)
      Bjoin_K <- 0.5 * BMSY_K
    if (BMSY_K > 0.3 & BMSY_K < 0.5)
      Bjoin_K <- 0.75 * BMSY_K - 0.075
    Bjoin <- Bjoin_K * Kc
    PBjoin <- prodPTF(Bjoin_K, n, MSY)
    cp <- (1 - n) * g * MSY * (Bjoin^(n - 2)) * Kc^-n
    Bc <- rep(NA, nys)
    Bc[1] <- Kc
    obj <- 0
    for (yr in 2:nys) {
      yref <- max(1, yr - adelay)
      if (Bc[yref] > Bjoin | BMSY_K > 0.5) {
        Bc[yr] <- Bc[yr - 1] + g * MSY * (Bc[yref]/Kc) - g * MSY *
          (Bc[yref]/Kc)^n - C_hist[yr - 1]
      } else {
        Bc[yr] <- Bc[yr - 1] + Bc[yref] * ((PBjoin/Bjoin) + cp * (Bc[yref] -
                                                                    Bjoin)) - C_hist[yr - 1]
      }
      if (Bc[yr] < 0)
        obj <- obj + log(-Bc[yr])
      Bc[yr] <- max(1e-06, Bc[yr])
    }
    obj + ((Bc[nys]/Kc) - Bt_K)^2
  } # end of DBSRA optimization function
  ####################################################### 
  C_hist <- Data@Cat[x, ]
  TAC <- rep(NA, reps)
  DBSRAcount <- 1
  if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) 
    return(NA)
  while (DBSRAcount < (reps + 1)) {
    #### small modificitation in the prior for the non-beta distributions by AC #### 
    Bt_K <- rbeta(10000, alphaconv(Data@Dep, Data@Dep*Data@CV_Dep), betaconv(Data@Dep, Data@Dep*Data@CV_Dep)) ### not using [x]for the non-beta distributions
    Bt_K[is.na(Bt_K)]<-0   # removes nas created in the non-beta distributions
    Bt_K <- sample(Bt_K,1)
   ###########################################################################
    
    Mdb <- round(rlnorm(1,log(Data@Mort[x]),Data@Mort[x]*Data@CV_Mort[x]),2) ## M changed to fit SSS input by AC ##
    Mdb <- Mdb[Mdb < 0.9][1]
    if (is.na(Mdb)) 
      Mdb <- 0.9
    FMSY_M <- trlnorm(1, Data@FMSY_M[x], Data@CV_FMSY_M[x])
    BMSY_K <- rbeta(10000, alphaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
                                       Data@BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
                                                                    Data@BMSY_B0[x]))
    tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1]
    if (is.na(tryBMSY_K)) {
      Min <- min(BMSY_K, na.rm = TRUE)
      Max <- max(BMSY_K, na.rm = TRUE)
      if (Max <= 0.05) 
        BMSY_K <- 0.05
      if (Min >= 0.95) 
        BMSY_K <- 0.95
    }
    
    
    
    if (!is.na(tryBMSY_K)) 
      BMSY_K <- tryBMSY_K
    adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
                            Data@L50[x])), 1)
    scaler <- 1000/mean(C_hist)
    C_hist2 <- scaler * C_hist
    opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist2), 
                                    1000 * mean(C_hist2))), C_hist = C_hist2, nys = length(C_hist2), 
                    Mdb = Mdb, FMSY_M = FMSY_M, BMSY_K = BMSY_K, Bt_K = Bt_K, 
                    adelay = adelay, tol = 0.01)
    Kc <- exp(opt$minimum)/scaler
    BMSYc <- Kc * BMSY_K
    FMSYc <- Mdb * FMSY_M
    UMSYc <- (FMSYc/(FMSYc + Mdb)) * (1 - exp(-(FMSYc + Mdb)))
    MSYc <- Kc * BMSY_K * UMSYc
    TAC[DBSRAcount] <- UMSYc * Kc * Bt_K
    DBSRAcount <- DBSRAcount + 1
  }

  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
  
}

class(DBSRA_new)<-"MP"
sfExport('DBSRA_new')