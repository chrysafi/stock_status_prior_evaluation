### load required packages ###
library(DLMtool)
library(datalimited2) 
library(r4ss) ## required for changing SS input files
library(sss) 

### functions required for the beta distribution ###
alphaconv<-function(m,sd)
  ( m*(((m*(1-m))/(sd^2))-1)) 

betaconv<-function(m,sd)
  ((1-m)*(((m*(1-m))/(sd^2))-1) )
####################################################

## ----Prerequisites_sfinit,echo=TRUE,warning=F----------------------------
sfInit(parallel=T,cpus=2) 
sfExportAll() 

## workign directory
Dir.in<- "C:/Users/"   #  Add the directory that the files are located
Dir.in<- "C:/Users/chrysafi/Documents/data_poor_methods/models_performance/Github functions/"

###### import data ################
## experts' priors
setwd(paste(Dir.in,"DLM_data/",sep=""))
EXP_1<-read.table("EXP_1.txt",header=TRUE)
EXP_2<-read.table("EXP_2.txt",header=TRUE)
EXP_3<-read.table("EXP_3.txt",header=TRUE)
## true stock status (SS derived)
Dep<-read.table("Dep.txt",header=TRUE)

## Species names to create new DLM ojects
S<-c("Aurora","Bocaccio","Cabezon","Canary","Darkblotched","Dover","Petrale","ST","LT","Widow",
     "Arrow","GHL","Aplaice","Mackerel","Pop","Sablefish","Yellow","N_Rockfish")

## PSA derived mean and sd
PSA<-read.table(paste(Dir.in,"DLM_data/PSA_S.txt",sep=""),header=TRUE)

reps<-1000
## Array for saving the model derived OFLs
Results<-array(NA,dim=c(9,6,reps),dimnames=list(c("exp1","exp2","exp3","pool","exp_correcte", "BRT","PSA","0.4","TRUE"),c("DCAC","DBSRA","SPMSY","CMSY","OCOM","SSS"),c(1:reps)))
###################################################################################################################
## for loop to run all models for all tested cases for each species. Can also ( and should) be broken in smaller
## loops to run simultaneously multipple times SSS that requires long time to run
###################################################################################################################
for (i in 1:10){
  ## creating a new data object for each species
    data<-new('Data',paste(Dir.in,"DLM_data/", S[i],"_DLM.csv",sep=""))
  
  #Standard input  
  data@BMSY_B0<-0.4
  data@CV_BMSY_B0<-0.1/data@BMSY_B0
  data@CV_Mort<-0.1/data@Mort
  data@CV_vbLinf<-0.1/data@vbLinf
  data@CV_vbK<-0.1/data@vbK
  data@CV_L50<-0.1/data@L50
  
  ##############################################################################################
  ## mean and sd from experts ## 
  
  M<-c(EXP_1$mean[EXP_1$Species_Scenario==paste(i+1,4,sep="")],EXP_2$mean[EXP_2$Species_Scenario==paste(i+1,4,sep="")],EXP_3$mean[EXP_3$Species_Scenario==paste(i+1,4,sep="")])
  SD<-c(EXP_1$std[EXP_1$Species_Scenario==paste(i+1,4,sep="")],EXP_2$std[EXP_2$Species_Scenario==paste(i+1,4,sep="")],EXP_3$std[EXP_3$Species_Scenario==paste(i+1,4,sep="")])
  
  ### for loop for experts 
  for (j in 1:3){
    
    data@Dep<-M[j]
    data@CV_Dep<-SD[j]/M[j]
    
    ### Did not do paraller execution but can also be done
   
    Results[j,1,]<-DCAC_2000(1,data,reps)@TAC
    Results[j,2,]<-DBSRA_new(1,data,reps)@TAC
    Results[j,3,]<-SPMSY_new(1,data,reps)@TAC
    Results[j,4,]<-CMSY_new(1,data,reps)@TAC
    Results[j,5,]<-OCOM_new(1,data,reps)@TAC
    Results[j,6,]<-SSS_new(1,data,reps)@TAC
  }
  
  #### Opinion pool of experts
  
  A<-rbeta(10000, alphaconv(M[1],SD[1]),betaconv(M[1],SD[1]))
  B<-rbeta(10000, alphaconv(M[2],SD[2]),betaconv(M[2],SD[2]))
  C<-rbeta(10000, alphaconv(M[3],SD[3]),betaconv(M[3],SD[3]))
  
  ## equal weights
  Opinion_pool<-c(sample(A,5000), sample(B,5000), sample(C,5000))

  #### mean and sd for opinion pool
  mean<-sample(Opinion_pool,2000)
  sd<-density(Opinion_pool)$bw
  data@Dep<-mean
  data@CV_Dep<-sd/mean 
  

  Results[4,1,]<-DCAC_2000(1,data,reps)@TAC
  Results[4,2,]<-DBSRA_new(1,data,reps)@TAC
  Results[4,3,]<-SPMSY_new(1,data,reps)@TAC
  Results[4,4,]<-CMSY_new(1,data,reps)@TAC
  Results[4,5,]<-OCOM_new(1,data,reps)@TAC
  Results[4,6,]<-SSS_new(1,data,reps)@TAC
  
  ######### Expert corrected bias from Perälä et al. In revision #####
  
  EBC<-read.table(paste(Dir.in,"EBC_Priors/ebc_",i,".txt",sep=""),header=TRUE)
  EBC<-unlist(EBC)
  mean<-EBC
  sd<-density(EBC)$bw
  ebc<-rbeta(2000,alphaconv(mean,sd),betaconv(mean,sd))
  ebc[is.na(ebc)]<-0
  data@Dep<-mean
  data@CV_Dep<-sd/mean 
  
  Results[5,1,]<-DCAC_2000(1,data,reps)@TAC
  Results[5,2,]<-DBSRA_new(1,data,reps)@TAC
  Results[5,3,]<-SPMSY_new(1,data,reps)@TAC
  Results[5,4,]<-CMSY_new(1,data,reps)@TAC
  Results[5,5,]<-OCOM_new(1,data,reps)@TAC
  Results[5,6,]<-SSS_new(1,data,reps)@TAC
  
  
  ########### BRT derived status prior ################
  

  #### function for bootstrapping####################
  Sdistrib = function(n, s_mean) {
    nv = 0 ; n.redo = 0
    while(nv < n) {
      n.redo = n.redo+1
      if(s_mean<=0.5) {
        si1 = fGarch::rsnorm(n*n.redo, mean=max(s_mean,0)-0.072, sd=0.189, xi=0.763)
        si = si1[si1>0 & si1<1]
      } else if(s_mean>0.5) {
        si1 = fGarch::rsnorm(n*n.redo, mean=max(s_mean,0)+0.179, sd=0.223, xi=0.904)
        si = si1[si1>0 & si1<1] }
      if(length(si)>n) si = sample(si,n);
      nv = length(si) }
    return (si)
  }
  ### Estimating stock status for each species from BRTs with the zbrt function in datalimited2 package
  output<-zbrt(year=data@Year,catch=as.vector(data@Cat))
  
  s_mean<-tail(output$ts[[5]],n=1)
  S_draws<-Sdistrib(2000,s_mean)
  
  S_dist<-S_draws
  mean<-S_dist
  sd<-density(S_dist)$bw
  
  data@Dep<-mean
  data@CV_Dep<-sd/mean
  
  Results[6,1,]<-DCAC_2000(1,data,reps)@TAC
  Results[6,2,]<-DBSRA_new(1,data,reps)@TAC
  Results[6,3,]<-SPMSY_new(1,data,reps)@TAC
  Results[6,4,]<-CMSY_new(1,data,reps)@TAC
  Results[6,5,]<-OCOM_newC(1,data,reps)@TAC
  Results[6,6,]<-SSS_new(1,data,reps)@TAC
  
  
  ### PSA derived priors #####
  
  data@Dep<-PSA[i,2]
  data@CV_Dep<-PSA[i,3]/PSA[i,2]
  
  Results[7,1,]<-DCAC_2000(1,data,reps)@TAC
  Results[7,2,]<-DBSRA_new(1,data,reps)@TAC
  Results[7,3,]<-SPMSY_new(1,data,reps)@TAC
  Results[7,4,]<-CMSY_new(1,data,reps)@TAC
  Results[7,5,]<-OCOM_new(1,data,reps)@TAC
  Results[7,6,]<-SSS_new(1,data,reps)@TAC
  
  #### assuming stock is at B40% ##########
  d<-0.4   
  sd<-0.2
  data@Dep<-d
  data@CV_Dep<-sd/d
  
  Results[8,1,]<-DCAC_2000(1,data,reps)@TAC
  Results[8,2,]<-DBSRA_new(1,data,reps)@TAC
  Results[8,3,]<-SPMSY_new(1,data,reps)@TAC
  Results[8,4,]<-CMSY_new(1,data,reps)@TAC
  Results[8,5,]<-OCOM_new(1,data,reps)@TAC
  Results[8,6,]<-SSS_new(1,data,reps)@TAC
 
  #### Using true (SS derived) stock sttus ####
  data@Dep<-Dep[i,2]
  data@CV_Dep<-0.2/Dep[i,2]
  
  Results[9,1,]<-DCAC_2000(1,data,reps)@TAC
  Results[9,2,]<-DBSRA_new(1,data,reps)@TAC
  Results[9,3,]<-SPMSY_new(1,data,reps)@TAC
  Results[9,4,]<-CMSY_new(1,data,reps)@TAC
  Results[9,5,]<-OCOM_new(1,data,reps)@TAC
  Results[9,6,]<-SSS_new(1,data,reps)@TAC
  
  
  
  
  
  D<-Results
  ##### save and read multidimensional array #####################
  saveRDS(D, file=paste(Dir.in,"Results",i,".rda",sep=""))
}

### END ###