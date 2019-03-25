#### DCAC adjusted so the Average catch is calculated for the years between Cmax and year 2000
#### when managed was more strict since


DCAC_2000<-function (x, Data, reps = 100) 
  
{
  dependencies = "Data@Cat, Data@Mort, Data@CV_Mort, Data@Dt, Data@CV_Dt,Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0"
  if (is.na(Data@BMSY_B0[x]) | is.na(Data@CV_BMSY_B0[x])) 
    return(NA)
  ##### average catch from year that max catch ir reached until 2000 by AC######
  str<-which(Data@Cat[x,]==max(Data@Cat[x,]))
  end<-which(Data@Year==2000)
  nyears<-length(str:end)
  C_tot<-mean(Data@Cat[x,str:end])*nyears
  #########################################################################
  Mdb <- round(rlnorm(reps,log(Data@Mort[x]),Data@Mort[x]*Data@CV_Mort[x]),2) ## M changed to fit SSS input by AC ##
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
  
  Bt_K <- rbeta(reps, alphaconv(Data@Dep,Data@Dep*Data@CV_Dep), ### not using [x]for the non-beta distributions 
                betaconv(Data@Dep,Data@Dep*Data@CV_Dt))
  Bt_K[is.na(Bt_K)]<-0   # removes nas created in the non-beta distributions
  
  Bt_K[Bt_K > 1] <- 1
  if (any(is.na(c(Data@BMSY_B0, Data@CV_BMSY_B0)))) 
    return(NA)
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x]), 
                  betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x]))
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(C_tot/(nyears + ((1 - Bt_K)/(BMSY_K * 
                                                               FMSY_M * Mdb))))
  Rec
}
class(DCAC_2000)<-"MP"
sfExport("DCAC_2000")





