## code for generic changes of SSS input files from a data object in DLM and running SSS from th sss package
## SSS package available in github
#install.packages("devtools")
#library(devtools)
#devtools::install_github("shcaba/SSS", build_vignettes = TRUE)

SSS_new<-function (x,data, reps) 
  
{
    dependencies = "Data@Mort,Data@Dep, Data@vbLinf,Data@vbK, Data@Year, Data@Cat,Data@FMSY_M, Data@BMSY_B0
  Data@CV_Dep,Data@CV_Mort,Data@MaxAge, Data@L50, Data@Rec, Data@wla, Data@wlb"
  
  Dir.in<-paste("C:/Users/" ,"SS_input_files/",sep="")  #  Add the directory that the files are located
  
  ########### change forecast file for each species #######################
  ########## forecast file #####################
  frcst.file<-SS_readforecast(paste(Dir.in,"/forecast.ss",sep=""), verbose = TRUE)
  frcst.file$FirstYear_for_caps_and_allocations<-data@Year[length(data@Year)]+3
  frcst.file$Yinit<-data@Year[length(data@Year)]+1
  SS_writeforecast(frcst.file, dir=Dir.in,file="forecast.ss", overwrite=TRUE,verbose=TRUE)
  ########################################################
  dat.file<-SS_readdat_3.30(paste(Dir.in,"/simple_data.dat",sep=""), verbose = TRUE, echoall = FALSE, section = NULL)
  dat.file$styr<- data@Year[1] ##start year
  dat.file$endyr<-data@Year[length(data@Year)] ##end year
  dat.file$Nsexes<-2 ### number of sexes 
  Lmax<-data@vbLinf[x] +2  ##Lmax for the species
  dat.file$minimum_size<-2
  dat.file$maximum_size<-Lmax+6
  
  dat.file$Nages<-round(data@MaxAge[x]*0.9) ##Amax 90% of the maximum age 
  all.yrs<-length(data@Year) ##all years 
  ## chanching the catches ###
  dat.file$catch<-as.data.frame(matrix(NA,nrow=all.yrs+1,ncol=5))
  names(dat.file$catch)<-paste(c("V1","V2","V3","V4","V5"))
  dat.file$catch$V1<-c(-999,data@Year)
  dat.file$catch$V2<-rep(1,all.yrs+1)
  dat.file$catch$V3<-rep(1,all.yrs+1)
  dat.file$catch$V4<-c(0,data@Cat) ###catches for all years
  dat.file$catch$V5<-rep(0.05,all.yrs+1)
  nlbins<-length(seq(2,Lmax+4,2))
  dat.file$N_lbins<-nlbins
  dat.file$lbin_vector<-seq(2,Lmax+4,2)
  SS_writedat_3.30(dat.file, outfile=paste(Dir.in,"/simple_data.dat",sep=""), overwrite=TRUE)
  dat.new<-readLines(paste(Dir.in,"/","simple_data.dat",sep=""))
  L<-20+all.yrs+13
  dat.new[L-1]<-paste(data@Year[1],"1 2 1 0.01 #" ,sep=" ")
  dat.new[L]<- paste(data@Year[length(data@Year)],"1 2 1 0.01 # DEPLETION" ,sep=" ")
  dat.new[L+1]<-"-9999 1 1 1 1"
  ### required to make data file operational ####
  line3<-strsplit(dat.new[grep("-9999\t0\t0\t0\t0\t0\t0\t0\t0\t#_terminator",dat.new)], " ")[[1]]
  line3<-"#-9999\t0\t0\t0\t0\t0\t0\t0\t0\t#_terminator"
  dat.new[grep("-9999\t0\t0\t0\t0\t0\t0\t0\t0\t#_terminator",dat.new)]<-paste(line3,collapse=" ")
  write(dat.new,paste(Dir.in,"/","simple_data.dat",sep=""))
  ###################################################################
  ### 
  file changes ###
  ctl.new<-readLines(paste(Dir.in,"/simple_control.ctl",sep=""))
  ctl.new[40]<-paste(1, " #_Growth_Age_for_L1( 999 to use as Linf)", sep="")
  ################# females ####################
  L50.line<-strsplit(ctl.new[grep("Mat50%_Fem",ctl.new)], " ")[[1]]
  L50.line[c(5,7)]<-round(data@L50[x])
  L50.line[c(1,3)]<-c(round(data@L50[x])-1,round(data@L50[x])+1)
  ctl.new[grep("Mat50%_Fem",ctl.new)]<-paste(L50.line,collapse=" ")
  
  #### length at age 1 ####
  L1<-data@vbLinf[x]*(1-exp(-data@vbK[x]*(1-data@vbt0[x])))
  
  L1.line<-strsplit(ctl.new[grep("L_at_Amin_Fem_GP_1",ctl.new)], " ")[[1]]
  L1.line[c(5,6)]<-round(L1)
  L1.line[c(1,3)]<-c(1,round(L1)+5)
  ctl.new[grep("L_at_Amin_Fem_GP_1",ctl.new)]<-paste(L1.line,collapse=" ")
  
  Linf.line<-strsplit(ctl.new[grep("L_at_Amax_Fem_GP_1",ctl.new)], " ")[[1]]
  Linf.line[c(5,6)]<-data@vbLinf[x]
  Linf.line[c(1,3)]<-c(Lmax-7,Lmax+3)
  ctl.new[grep("L_at_Amax_Fem_GP_1",ctl.new)]<-paste(Linf.line,collapse=" ")
  
  k.line<-strsplit(ctl.new[grep("VonBert_K_Fem_GP_1",ctl.new)], " ")[[1]]
  k.line[c(3,5)]<- data@vbK[x]
  k.line[c(1,2)]<-c(0,1)
  ctl.new[grep("VonBert_K_Fem_GP_1",ctl.new)]<-paste(k.line,collapse=" ")
  
  A.line<-strsplit(ctl.new[grep("Wtlen_1_Fem",ctl.new)], " ")[[1]]
  A.line[c(3,5)]<-data@wla[x]
  ctl.new[grep("Wtlen_1_Fem",ctl.new)]<-paste(A.line,collapse=" ")
  
  A.line<-strsplit(ctl.new[grep("Wtlen_2_Fem",ctl.new)], " ")[[1]]
  A.line[c(3,5)]<-data@wlb[x]
  ctl.new[grep("Wtlen_2_Fem",ctl.new)]<-paste(A.line,collapse=" ")
  ################### males ############################ 
  L1.line<-strsplit(ctl.new[grep("L_at_Amin_Mal_GP_1",ctl.new)], " ")[[1]]
  L1.line[c(4,5)]<-round(L1)
  L1.line[c(1,3)]<-c(1,round(L1)+5)
  ctl.new[grep("L_at_Amin_Mal_GP_1",ctl.new)]<-paste(L1.line,collapse=" ")
  
  Linf.line<-strsplit(ctl.new[grep("L_at_Amax_Mal_GP_1",ctl.new)], " ")[[1]]
  Linf.line[c(4,5)]<-Lmax-2
  Linf.line[c(1,3)]<-c(Lmax-7,Lmax+3)
  ctl.new[grep("L_at_Amax_Mal_GP_1",ctl.new)]<-paste(Linf.line,collapse=" ")
  
  k.line<-strsplit(ctl.new[grep("VonBert_K_Mal_GP_1",ctl.new)], " ")[[1]]
  k.line[c(4,6)]<- data@vbK[x]
  k.line[c(1,3)]<-c(0,1)
  ctl.new[grep("VonBert_K_Mal_GP_1",ctl.new)]<-paste(k.line,collapse=" ")
  
  A.line<-strsplit(ctl.new[grep("Wtlen_1_Mal",ctl.new)], " ")[[1]]
  A.line[c(3,5)]<-data@wla[x]
  ctl.new[grep("Wtlen_1_Mal",ctl.new)]<-paste(A.line,collapse=" ")
  
  A.line<-strsplit(ctl.new[grep("Wtlen_2_Mal",ctl.new)], " ")[[1]]
  A.line[c(3,5)]<-data@wlb[x]
  ctl.new[grep("Wtlen_2_Mal",ctl.new)]<-paste(A.line,collapse=" ")
  
  
  
  ctl.new[99]<-paste(data@Year[1]-1, " # first year of main recr_devs; early devs can preceed this era", sep="")
  ctl.new[100]<-paste(2016, " # Last year of main recr_devs; early devs can preceed this era", sep="") # setting it to 2016 improves performance
  
  
  #### selectivity changes #############
  slope<- -1
  intercept<-round(data@L50[x])
  x<-seq(1,50,0.05)
  mat.vec<-1/(1+exp(slope*(x-intercept)))
  mat.vec[mat.vec>0.95]
  width.95per<-(x[mat.vec>0.95]-intercept)[1]
  
  ##############
  
  
  
  Sel_1<-strsplit(ctl.new[grep("SizeSel_P1_Fishery",ctl.new)], " ")[[1]]
  Sel_1[25]<-Lmax
  Sel_1[35]<-data@L50[x]
  Sel_1[47]<-data@L50[x]
  ctl.new[grep("SizeSel_P1_Fishery",ctl.new)]<-paste(Sel_1,collapse=" ")
  
  Sel_2<-strsplit(ctl.new[grep("SizeSel_P2_Fishery",ctl.new)], " ")[[1]]
  Sel_2[22]<-round(Lmax)
  Sel_2[41]<-data@L50[x]
  Sel_2[29]<-width.95per 
  ctl.new[grep("SizeSel_P2_Fishery",ctl.new)]<-paste(Sel_2,collapse=" ")
  write(ctl.new,paste(Dir.in,"/","simple_control.ctl",sep=""))
  
  #########changes in the BMSYBO file #############################
  filepath<-Dir.in
  starter.new<-readLines(paste(filepath,"/starter.ss",sep=""))
  starter.new[7]<-1
  write(starter.new,paste(filepath,"/starter.ss",sep=""))
  RUN.SS(paste(filepath,"/",sep=""), ss.exe="ss") #  ,ss.cmd=" -nohess -nox > out.txt 2>&1")
  
  rep.new<-readLines(paste(filepath,"/Report.sso",sep=""))
  #line<-strsplit(rep.new[grep("AGE_SELEX",rep.new)], " ")
  posF<-grep('Fecund', rep.new)
  Fec<-rep.new[posF[2]]
  Fec<-unlist(strsplit(Fec, split=" "))
  Fec<-Fec[10:length(Fec)]
  posAS<-grep('Asel2', rep.new)
  ASF<-rep.new[posAS[2]]
  ASF<-unlist(strsplit(ASF, split=" "))
  ASF<-ASF[8:length(ASF)]
  
  ASM<-rep.new[posAS[4]]
  ASM<-unlist(strsplit(ASM, split=" "))
  ASM<-ASM[8:length(ASM)]
  posBW<-grep('bodywt', rep.new)
  BWF<-rep.new[posBW[1]]
  BWF<-unlist(strsplit(BWF, split=" "))
  BWF<-BWF[8:length(BWF)]
  BWM<-rep.new[posBW[3]]
  BWM<-unlist(strsplit(BWM, split=" "))
  BWM<-BWM[8:length(BWM)]
  
  
  BMSYdat<-readLines(paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
  BMSYdat[3]<-round(data@MaxAge*0.9)
  BMSYdat[7]<-paste(ASF, collapse=' ' )
  BMSYdat[8]<-paste(ASM, collapse=' ' )
  BMSYdat[10]<-paste(BWF, collapse=' ' )
  BMSYdat[11]<-paste(BWF, collapse=' ' )
  BMSYdat[13]<-paste(Fec, collapse=' ' )
  write(BMSYdat,paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
  starter.new<-readLines(paste(filepath,"/starter.ss",sep=""))
  starter.new[7]<-2
  write(starter.new,paste(filepath,"/starter.ss",sep=""))
  file.remove(paste(filepath,"/Report.sso",sep=""))
  ################################################################################## 
  
  if (length(data@Dep)==1){  ### for beta distributions
  
  SSS_run<-SSS(filepath=Dir.in,
                  file.name=c("simple_data.dat","simple_control.ctl"),
                  reps=reps,
                  seed.in=19,
                  M.in=c(3,data@Mort[x],data@Mort[x]*data@CV_Mort[x],3,data@Mort[x],data@Mort[x]*data@CV_Mort[x]),
                  Dep.in=c(2,data@Dep[x],data@Dep[x]*data@CV_Dep[x]),
                  SR_type=9,
                  h.in=c(-2,0.779,0.152), ### negative value in steepness means the model uses FMSY/M and BMSY/B0
                  FMSY_M.in=c(30,data@FMSY_M[x],round(data@FMSY_M[x]*data@CV_FMSY_M[x],2)),
                  BMSY_B0.in=c(2,data@BMSY_B0[x],data@BMSY_B0[x]*data@CV_BMSY_B0[x]),
                  L1.in=c(0,0,0,0),
                  Linf.in=c(0,0,0,0),
                  k.in=c(0,0,0,0),
                  Zfrac.Beta.in=c(-99,0.2,0.6,-99,0.5,2),
                  R_start=c(0,data@Rec[1,2]),
                  sum_age=0,
                  sb_ofl_yrs=c(tail(data@Year, n=1)+1,tail(data@Year, n=1)+1,tail(data@Year, n=1)+2),
                  f_yr=data@Year[length(data@Year)],
                  year0=data@Year[1],
                  genders=T) 
  } else {   ### for the non-beta priors ##
    SSS_run<-SSS(filepath=Dir.in,
                 file.name=c("simple_data.dat","simple_control.ctl"),
                 reps=reps,
                 seed.in=19,
                 M.in=c(3,data@Mort[x],data@Mort[x]*data@CV_Mort[x],3,data@Mort[x],data@Mort[x]*data@CV_Mort[x]),
                 Dep.in=data@Dep,
                 SR_type=9,
                 h.in=c(-2,0.779,0.152), ### negative value in steepness means the model uses FMSY/M and BMSY/B0
                 FMSY_M.in=c(30,data@FMSY_M[x],round(data@FMSY_M[x]*data@CV_FMSY_M[x],2)),
                 BMSY_B0.in=c(2,data@BMSY_B0[x],data@BMSY_B0[x]*data@CV_BMSY_B0[x]),
                 L1.in=c(0,0,0,0),
                 Linf.in=c(0,0,0,0),
                 k.in=c(0,0,0,0),
                 Zfrac.Beta.in=c(-99,0.2,0.6,-99,0.5,2),
                 R_start=c(0,data@Rec[1,2]),
                 sum_age=0,
                 sb_ofl_yrs=c(tail(data@Year, n=1)+1,tail(data@Year, n=1)+1,tail(data@Year, n=1)+2),
                 f_yr=data@Year[length(data@Year)],
                 year0=data@Year[1],
                 genders=T) 
  }
  
  OFL<-SSS_run$Posteriors[,14]
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(OFL)
  Rec
  
}
class(SSS.AC)<-"MP"

sfExport('SSS_new')

