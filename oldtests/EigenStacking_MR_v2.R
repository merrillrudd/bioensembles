rm(list=ls())
main_dir <- "C:\\merrill\\bioensembles\\sim"

library(FishLife)
library(mvQuad)
#library(mvtnorm)

library(LIME)
library(LBSPR)
library(mvtnorm)
library(Hmisc)
library(foreach)
library(doParallel)

#########################
# Step 1:  Generate 4 parameters from 1/2 dimensional quadrature nodes
#########################

# Choose species
sp <- Plot_taxa( Search_species(Genus="Lachnolaimus", Species="maximus")$match_taxonomy, mfrow=c(2,2) )

# Choose number of dimensions for approximation, must be 1 or 2
Dim = 2

# species-level
mp <- sp[[1]]$Mean_pred
cov <- sp[[1]]$Cov_pred

# Extract just the ones we want
Mean = mp[c('M','K','Loo','Lm')]
Cov = cov[c('M','K','Loo','Lm'),c('M','K','Loo','Lm')]

# Bivariate quadrature
myGrid <- createNIGrid(dim=Dim, type="GHe", level=5,  ndConstruction="sparse")
XY_i = getNodes(myGrid)

# Use eigen-decomposition to project quadrature nodes into all variables
Eigen = eigen(Cov)
Diag = function(vec){
  if( length(vec)==1 ) return(vec)
  if( length(vec)>=2 ) return(diag(vec))
}
Param_i = t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(XY_i) )
colnames(Param_i) = names(Mean)

#########################
# Step 2:  Run LBSPR using each row i of Param_i, and extract SPR_i for each
#########################
res_dir <- file.path(main_dir, "results")
dir.create(res_dir, showWarnings=FALSE)

itervec <- 1:20

run_dir <- file.path(res_dir, "2param")
dir.create(run_dir, showWarnings=FALSE)

  ##################################
  ## LBSPR 
  ##################################

  ## Simulate true data with species-level information
  ## Run all models with species-level distributions
  ncores <- 5
  cl <- makeCluster(ncores)
  registerDoParallel(cl)      
  start <- Sys.time()
  foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
    runstack(savedir=run_dir, 
          nodes=Param_i, 
          param=c("K","M","Loo","Lm"), 
          mean=Mean, 
          cov=Cov, 
          modname="2param_Species", 
          simulation=TRUE, 
          iter=loop, 
          seed=loop, 
          Fscenario="equil", 
          rewrite=TRUE, 
          model="LBSPR", 
          sim_model="LIME",
          Nyears=1)
  end <- Sys.time() - start   
  stopCluster(cl)


#########################
# Step 3:  Smooth SPR_i in original 1/2 dimensional space
#########################
re <- rem <- rei <- mean <- var <- tval <- rep(NA, length(itervec))
for(i in 1:length(itervec)){
  true <- readRDS(file.path(run_dir, itervec[i], "True.rds"))
  mres <- readRDS(file.path(run_dir, itervec[i], "2param_Species_res_Means_LBSPR.rds"))
  ires <- readRDS(file.path(run_dir, itervec[i], "res_IterTrue_LBSPR.rds"))
  stack <- readRDS(file.path(run_dir, itervec[i], "2param_Species_res_stacking_LBSPR.rds"))
  SPR_i <- sapply(1:length(stack), function(x) stack[[x]]@SPR)

if( Dim==1 ){
  # Step 3A:  Define axis
  # Grid_j = seq(-4,4,length=100)
  Grid_j = seq(min(Param_i[,1]),max(Param_i[,1]), length=100)
  # Step 3B:  Use approx() to do 1D smoother of SPR_i for each row of Grid_j
  interp_h = approx(x=Param_i[,1], y=SPR_i, n=100)
  SPR_j <- interp_h$y
  # Step 3C:  translate Grid to 4 parameters and calculate FishLife probability for each
  Param_j = t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j) )
  colnames(Param_j) = names(Mean)
  Prob_j = dmvnorm( Param_j, mean=Mean, sigma=Cov )
  # Step 3D:  Calculate mean
  SPRmean = weighted.mean( SPR_j, w=Prob_j, na.rm=TRUE )
}
if( Dim==2 ){
  # Step 3A:  Define axis
  # Grid_j = expand.grid( seq(-4,4,length=100), seq(-4,4,length=100) )
  Grid_j = expand.grid( seq(min(Param_i[,1]),max(Param_i[,1]),length=100), seq(min(Param_i[,2]),max(Param_i[,2]),length=100) )
  # Step 3B:  Use akima:interp to do 2D smoother of SPR_i for each row of Grid_j
  interp_h = akima::interp( x=Param_i[,1], y=Param_i[,2], z=SPR_i, nx=100, ny=100, duplicate="median" )
  SPR_j = interp_h$z
  # Step 3C:  translate Grid to 4 parameters and calculate FishLife probability for each
  Param_j = t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j) )
  colnames(Param_j) = names(Mean)
  Prob_j = dmvnorm( Param_j, mean=Mean, sigma=Cov )
  # Step 3D:  Calculate mean
  SPRmean = weighted.mean( SPR_j, w=Prob_j, na.rm=TRUE )
  SPRvar = wtd.var(SPR_j, w=Prob_j, na.rm=TRUE)
}
  re[i] <- (SPRmean - true$SPR)/true$SPR
  rem[i] <- (mres@SPR - true$SPR)/true$SPR
  rei[i] <- (ires@SPR - true$SPR)/true$SPR
  mean[i] <- SPRmean
  var[i] <- SPRvar
  tval[i] <- true$SPR
}

par(mfrow=c(5,2))
for(i in 1:10){
  hist(rlnorm(100000, log(mean[i]), sqrt(var[i])), xlim=c(0,1))
  abline(v=tval[i])
}

par(mfrow=c(1,1))
boxplot(cbind(re, rem, rei))

par(mfrow=c(2,1))
image(matrix(Prob_j, nrow=100))
image(SPR_j)