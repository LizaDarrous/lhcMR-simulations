### lHC-MR simulation optimisation script 
## Optimising over n different sets of starting points for each data generation
## Input: set up the parameters corresponding to the scenario, number of data generations, number of starting points,
# boolean indicating which parallelisation method to use "rslurm" vs "lapply"
## Output: .csv:detailed optimisation for each data generation, .csv:overall best optimisation for all data generations

library(data.table)
library(rhdf5)
library(lhcMR)
library(dplyr)
library(rslurm)

#### start to change
## set working directory
setwd("/data/sgg3/liza/SEM_Sim/R_trans")

## parameter set up
dS = 1 #scenario number
nsimi = 50 #number of data generations from previous scripts
nSP = 50 #number of starting points for each generation
saveSP = TRUE #to save the starting points used for each data generation
paral_method = "rslurm" #"lapply"
#### end to change

## single trait analysis to estimate pi, h2, i, followed by generating 50 SP for each data generation and optimisation
SP = 3; #number of starting points for single trait analysis to estimate pi, h2 and i
M = 1e7; #number of SNPs
nCores=1; #to set up for optimisation in lapply

write.table(t(c("piX", "piY", "h2X", "h2Y", "tX", "tY", "axy", "ayx", "iXY", "iX", "iY")), 
          paste0("SummaryResult_DataB",dS,"_",simi,".csv"), row.names = FALSE, col.names = FALSE, sep=",")

for(simi in 1:nsimi){
  load(paste0("./data-mat/exampleDataB",dS,"_sim",simi,".RData")) 
  wd = rep(1, length(sig1)) #ones(length(sig1),1);
  piX = thS[1];
  piY = thS[2];
  sigX = thS[3];
  sigY = thS[4];
  tX = thS[5];
  tY = thS[6];
  alp = thS[7];
  bet = thS[8];
  iX = thS[9];
  iY = thS[10];
  iXY = thS[11];
  h2X = piX*M*sigX^2;
  h2Y = piY*M*sigY^2;
  herX = h2X+tX^2;
  herY = h2Y+tY^2;
  
  ## SP calculation
  sp_piX = runif(SP,0,0.01)
  sp_h2X = runif(SP,0,0.5)
  sp_iX = runif(SP,0.5,1.5)
  
  # Get piX, iX (and total_h2) in a separate step - single trait analysis for X and Y
  para=cbind(sp_piX,sp_h2X,sp_iX)
  sp_mat=matrix(unlist(para), ncol=SP, byrow = FALSE)
  colnames(sp_mat)=colnames(para)
  par.df = data.frame(par=I(apply(sp_mat,1,as.list))) #generate a dataframe of lists for each row of parameters - input for rslurm/lapply
  
  test.exp <- parallel::mclapply(par.df[[1]], function(x) {
    theta = unlist(x) #theta = unlist(par.df[[1]][1])
    test1 = optim(theta, lhcMR:::singleTrait_likelihood,
                  betX=betXY[,1], pi1=pi1, sig1=sig1, w8s=wd, M=M,
                  m0=m0, nX=nX, bn=2^7, bins=10,
                  method = "Nelder-Mead",
                  control = list(maxit = 5e3))
    
    list("mLL"=test1$value,"par"=test1$par,"conv"=test1$convergence)
  }, mc.cores = nCores)
  test.exp = as.data.frame(t(matrix(unlist(test.exp), nrow=length(unlist(test.exp[1])))))
  
  test.out <- parallel::mclapply(par.df[[1]], function(x) {
    theta = unlist(x)
    test2 = optim(theta, lhcMR:::singleTrait_likelihood,
                  betX=betXY[,2], pi1=pi1, sig1=sig1, w8s=wd, M=M,
                  m0=m0, nX=nY, bn=2^7, bins=10,
                  method = "Nelder-Mead",
                  control = list(maxit = 5e3))
    
    list("mLL"=test2$value,"par"=test2$par,"conv"=test2$convergence)
  }, mc.cores = nCores)
  test.out = as.data.frame(t(matrix(unlist(test.out), nrow=length(unlist(test.out[1])))))
  
  colnames(test.exp)=c("mLL","piX", "h2X","iX","conv")
  colnames(test.out)=c("mLL","piX", "h2X","iX","conv")
  res_exp_min = test.exp[which(test.exp$mLL == min(test.exp$mLL)), ]
  res_exp = abs(res_exp_min[2:4]) #piX,h2X,iX
  res_out_min = test.out[which(test.out$mLL == min(test.out$mLL)), ]
  res_out = abs(res_out_min[2:4]) #piY,h2Y,iY
  
  pi_X = as.numeric(res_exp[1])
  pi_Y = as.numeric(res_out[1])
  h2_x = as.numeric(res_exp[2]) #total heritability, not used
  h2_y = as.numeric(res_out[2]) #total heritability, not used
  i_X = as.numeric(res_exp[3])
  i_Y = as.numeric(res_out[3])
  
  # Generate the rest of the starting points 
  sp_tX = runif(nSP,0,0.5)
  sp_tY = runif(nSP,-0.5,0.5)
  sp_h2X = max(0,h2_x-(sp_tX^2))
  sp_h2Y = max(0,h2_y-(sp_tY^2))
  sp_axy = replicate(nSP, (alp+runif(1,-0.1,0.1)))
  sp_ayx = replicate(nSP, (bet+runif(1,-0.1,0.1)))
  sp_iXY = replicate(nSP, (iXY+runif(1,-0.05,0.05))) #rep(i_XY,50)
  
  para=cbind(sp_h2X,sp_h2Y,sp_tX,sp_tY,sp_axy,sp_ayx,sp_iXY)
  sp_mat=matrix(unlist(para), ncol=7, byrow = FALSE)
  colnames(sp_mat)=c("sp_h2X","sp_h2Y","sp_tX","sp_tY","sp_axy","sp_ayx","sp_iXY")
  
  SP_list = list("iX"=i_X,"iY"=i_Y,"piX"=pi_X,"piY"=pi_Y,"sp_mat"=sp_mat)
  if(saveSP){
    save(SP_list, file = paste0("./data-mat/exampleDataB",dS,"_",simi,"_sp",nSP,".RData"))
  }
  
  ## optimisation of said set of starting points
  bn = 2^7
  bins = 10
  param="comp" #possibility to develop nesting in later version
  w8s=wd #must be the same name as in likelihood function
  iX=i_X;iY=i_Y;piX=pi_X;piY=pi_Y #must be the same name as in likelihood function
  parscale2 = c(1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1)
  nCores = NA #for lapply: update to number wanted, or leave as NA to use third of available cores
  if(is.na(nCores)){
    nCores = max(1,floor((parallel::detectCores())/3))
  } else {
    if(nCores > parallel::detectCores()){
      cat(print("Core number chosen is greater than cores available\n"))
      stop()
    }
  }
  
  # Generate a dataframe of lists for each row of parameters - input for rslurm/lapply
  par.df = data.frame(par=I(apply(sp_mat,1,as.list)))
  
  if(paral_method=="rslurm"){
    # parameter set up
    partition="sgg"
    account="sgg"
    SP_pair=nrow(sp_mat)
    
    cat(print("Running optimisation"))
    sjob = slurm_apply(f = lhcMR:::slurm_pairTrait_twoStep_likelihood, params = par.df, jobname = paste0(dS,"_",simi,"_optim"), nodes = SP_pair, cpus_per_node = 1,
                       global_objects = c("betXY","pi1","sig1","w8s","m0","M","nX","nY","piU","piX","piY",
                                          "iX","iY","param","bn","bins","parscale2"),
                       slurm_options = list(partition = partition, account = account),
                       submit = TRUE)
    # Keep a loop open till the job is completely done.
    wait_counter = 0
    while (wait_counter < 1) {
      wait_counter = 0
      if(lhcMR:::get_job_status_lhc(sjob)){
        wait_counter = wait_counter + 1
      } else{
        wait_counter = wait_counter
      }
    }
    # Get output of minus log likelihood (mLL) and estimated parameters from rslurm in the form of a table with nrows equal to nrows(par.df)
    res_temp = get_slurm_out(sjob, outtype = 'table')
    res_values = as.data.frame(do.call(rbind, res_temp[[1]]))
    res_values %>%
      dplyr::mutate(h2X = abs(h2X),
                    h2Y = abs(h2Y),
                    tX = abs(tX)) -> res_values
    res_values = cbind("SP"=c(1:nrow(res_values)),"mLL"=res_values[,1],"piX"=piX,"piY"=piY,res_values[,-1],"iX"=iX,"iY"=iY)
    write.csv(res_values, paste0("FullResult_DataB",dS,"_",simi,".csv"), row.names = FALSE)
    
    res_est = res_values[which(res_values$mLL==min(res_values$mLL)),]
    res_est = dplyr::select(res_est, -c(SP,mLL,conv))
    write.table(t(res_est), paste0("SummaryResult_DataB",dS,"_",simi,".csv"), sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    
  }
  if(paral_method=="lapply"){
    cat(print("Running optimisation"))
    test.res <- parallel::mclapply(par.df[[1]], function(x) {
      theta = unlist(x)
      test = optim(theta, lhcMR:::pairTrait_twoStep_likelihood,
                   betX=betXY, pi1=pi1, sig1=sig1, w8s=w8s, M=M,
                   m0=m0, nX=nX, nY=nY, pi_U=piU, pi_X=piX, pi_Y=piY, i_X=iX, i_Y=iY,
                   bn=2^7, bins=10,
                   method = "Nelder-Mead",
                   control = list(maxit = 5e3,
                                  parscale = parscale2))
      
      list("mLL"=test$value,"par"=test$par,"conv"=test$convergence)
    }, mc.cores = nCores)
    
    test.res = as.data.frame(t(matrix(unlist(test.res), nrow=length(unlist(test.res[1])))))
    colnames(test.res)=c("mLL", "h2X","h2Y","tX","tY","axy","ayx","iXY","conv")
    
    test.res %>%
      dplyr::mutate(h2X = abs(h2X),
                    h2Y = abs(h2Y),
                    tX = abs(tX)) -> res_values
    
    res_values = cbind("SP"=c(1:nrow(res_values)),"mLL"=res_values[,1],"piX"=piX,"piY"=piY,res_values[,-1],"iX"=iX,"iY"=iY)
    write.csv(res_values, paste0("FullResult_DataB",dS,"_",simi,".csv"), row.names = FALSE)
    
    res_est = res_values[which(res_values$mLL==min(res_values$mLL)),]
    res_est = dplyr::select(res_est, -c(SP,mLL,conv))
    write.table(res_est, paste0("SummaryResult_DataB",dS,"_",simi,".csv"), sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  
}
