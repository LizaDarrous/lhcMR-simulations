### lHC-MR simulation script 
## Generating n data simulations
## Input: set up the parameters corresponding to the scenario, bolean for more complex scenarios,
# boolean indicating detailed saving or not
## Output: betXY, thS (parameters), sig1 + pi1 (local LD structure), nX, nY, m0 (number of SNPs in local region), 
# piU (constant), if desired: hrX + hrU + hrY

library(data.table)
library(rhdf5)
library(lhcMR)
library(dplyr)
library(PearsonDS)

#### start to change
## set working directory
setwd("/data/sgg3/liza/SEM_Sim/R_trans")

## scenario number
dS = 1;
## data generations for this scenario
nsimi = 50
## boolean set up to indicate detailed save (with hrX/U/Y), kurtosis (non-normal scenario), or presence of two confounders (twoU)
detailedSave = FALSE
kurtosis = FALSE #if true change kurtosis level in first if clause below
twoUs = FALSE #if true change two U properties in firs if clause below
#### end to change

## read in pi1 and sig1 + R matrix
Rmat = fread("LD_vsGM2.txt")
c1 = fread("pi1-sigma1_vsGM2.txt")

pi1 = unlist(c1[,1]) # proportion of non-zero component of rho
sig1 = unlist(c1[,2]) # SD of non-zero component of rho
sim_num = length(sig1);  # number of genotyped SNPs
m0 = 5e3+1;  # n_local true sequence variants in each SNP region - used for pi1/sig1
rm(c1)#to save space

#### start to change
## create bX bY from these parameters
piX = 5e-3;
piU = 5e-2;
piY = 1e-2;
h2X = .25;
h2U = .3;
h2Y = .2;
qX = .3;
qY = .2;
alp = .3;
bet = 0;
M = 1e7;
sigX = sqrt(h2X/(M*piX));
sigU = sqrt(h2U/(M*piU));
sigY = sqrt(h2Y/(M*piY));
tX = qX*sqrt(h2U);
tY = qY*sqrt(h2U);
nX = 5e5;
nY = 5e5;
iX = 1;
iY = 1;
iXY = 0;

if(kurtosis){
  krt = 5; #2 5 10
}
if(twoUs){
  qX1 = .3;
  qY1 = .2;
  qX2 = .4;
  qY2 = .3;
  tX1 = qX1*sqrt(h2U);
  tY1 = qY1*sqrt(h2U);
  tX2 = qX2*sqrt(h2U);
  tY2 = qY2*sqrt(h2U);
}
#### end to change

system("mkdir data-mat")
for(simi in 1:nsimi){
  message(simi)
  
  hX =  matrix(rbinom(sim_num*m0,1,piX),ncol=m0); #binornd(1,piX,sim_num,m0);
  seli = which(hX==1); #find(hX==1);
  tmp = matrix(rnorm(length(seli),0,sigX),ncol=1); #normrnd(0,sigX,length(seli),1);
  hX[seli] = tmp;
  if(kurtosis){
    hX = matrix(rbinom(sim_num*m0,1,piX),ncol=m0);
    sv1 = which(hX==1);
    tmp3 = matrix(PearsonDS::rpearson(length(sv1),moments=c(mean=0,variance=sigX,skewness=0,kurtosis=kurt)),ncol=1);#pearsrnd(0,sigX,0,krt,length(sv1),1);
    hX[sv1] = tmp3;
  }
  
  hU =  matrix(rbinom(sim_num*m0,1,piU),ncol=m0); #binornd(1,piU,sim_num,m0);
  seli = which(hU==1); #find(hX==1);
  tmp = matrix(rnorm(length(seli),0,sigU),ncol=1); #normrnd(0,sigU,length(seli),1);
  hU[seli] = tmp;
  if(kurtosis){
    hU = binornd(1,piU,sim_num,m0);
    sv1 = find(hU==1);
    tmp3 = pearsrnd(0,sigU,0,krt,length(sv1),1);
    hU(sv1) = tmp3;
  }
  if(twoUs){
    hU1 = matrix(rbinom(sim_num*m0,1,piU),ncol=m0);
    seli = which(hU1==1);
    tmp = matrix(rnorm(length(seli),0,sigU),ncol=1)
    hU1[seli] = tmp;
    hU2 = matrix(rbinom(sim_num*m0,1,piU),ncol=m0);
    seli = which(hU2==1);
    tmp = matrix(rnorm(length(seli),0,sigU),ncol=1)
    hU2[seli] = tmp;
  }
  
  hY =  matrix(rbinom(sim_num*m0,1,piY),ncol=m0); #binornd(1,piY,sim_num,m0);
  seli = which(hY==1); #find(hX==1);
  tmp = matrix(rnorm(length(seli),0,sigY),ncol=1); #normrnd(0,sigY,length(seli),1);
  hY[seli] = tmp;
  if(kurtosis){
    hY = binornd(1,piY,sim_num,m0);
    sv1 = find(hY==1);
    tmp3 = pearsrnd(0,sigY,0,krt,length(sv1),1);
    hY(sv1) = tmp3;
  }
  
  hrX = rowSums(Rmat*hX); #sum(R.*hX,2);
  hrU = rowSums(Rmat*hU); #sum(R.*hU,2);
  hrY = rowSums(Rmat*hY); #sum(R.*hY,2);
  if(twoUs){
    hrU1 = rowSums(Rmat*hU1); #sum(R.*hU1,2);
    hrU2 = rowSums(Rmat*hU2); #sum(R.*hU2,2);
  }
  ## Noise is generated from a variance-covariance matrix (X LD-score intercept divided by nX, cross LD score intercept rescaled by dividing by sqrt of nX*nY, and Y LD score intercept divided by nY).
  R0 = matrix(c(iX/nX,iXY/sqrt(nX*nY),iXY/sqrt(nX*nY),iY/nY),ncol=2,nrow=2, byrow = T) # [iX/nX,iXY/sqrt(nX*nY);iXY/sqrt(nX*nY),iY/nY];
  C = chol(R0);
  noi = matrix(rnorm(sim_num*2,0,1),ncol=2) #normrnd(0,1,sim_num,2);
  noi = noi %*% C #noi*C;
  
  ## Generate effects
  betX = (hrX + ((bet*qY)+qX)*hrU + bet*hrY)/(1-(alp*bet)) + noi[,1] #(hrX+(bet*qY+qX)*hrU+bet*hrY)/(1-alp*bet)+noi(:,1);
  betY = (alp*hrX + ((alp*qX)+qY)*hrU + hrY)/(1-(alp*bet)) + noi[,2] #(alp*hrX+(alp*qX+qY)*hrU+hrY)/(1-alp*bet)+noi(:,2);
  if(twoUs){
    betX = hrX + ((bet*qY1)+qX1)*hrU1 + ((bet*qY2)+qX2)*hrU2 + bet*hrY + noi[,1];
    betY = alp*hrX + ((alp*qX1)+qY1)*hrU1 + ((alp*qX2)+qY2)*hrU2 + hrY + noi[,2];
  }
  betXY = cbind(betX,betY) # [betX,betY];
  
  #save parameters set up
  thS = c(piX,piY,sigX,sigY,tX,tY,alp,bet,iX,iY,iXY)
  if(twoUs){
    thS = c(piX,piY,sigX,sigY,tX1,tY1,alp,bet,iX,iY,iXY)
  }
  if(detailedSave){
    if(twoUs){
      save(betXY, hrX, hrU1, hrU2, hrY, thS, sig1, pi1, nX, nY, m0, piU, qX1, qX2, qY1, qY2,
           file = paste0("./data-mat/exampleDataB",dS,"_sim",simi,".RData"))
    }else{
      save(betXY, hrX, hrU, hrY, thS, sig1, pi1, nX, nY, m0, piU, file = paste0("./data-mat/exampleDataB",dS,"_sim",simi,".RData"))
    }
  }else{
    if(twoUs){
      save(betXY, thS, sig1, pi1, nX, nY, m0, piU, qX1, qX2, qY1, qY2, file = paste0("./data-mat/exampleDataB",dS,"_sim",simi,".RData"))
    }else{
      save(betXY, thS, sig1, pi1, nX, nY, m0, piU, file = paste0("./data-mat/exampleDataB",dS,"_sim",simi,".RData"))
    }
  }
}

