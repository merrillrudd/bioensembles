
library(FishLife)
library(mvQuad)
#library(mvtnorm)

#########################
# Step 1:  Generate 4 parameters from 1/2 dimensional quadrature nodes
#########################

# Choose species
sp <- Plot_taxa( Search_species(Genus="Lutjanus", Species="guttatus")$match_taxonomy, mfrow=c(2,2) )

# Choose number of dimensions for approximation, must be 1 or 2
Dim = 4

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
# Step 2:  Run LIME using each row i of Param_i, and extract SPR_i for each
#########################

# Add code

#########################
# Step 3:  Smooth SPR_i in original 1/2 dimensional space
#########################

if( Dim==1 ){
  # Step 3A:  Define axis
  Grid_j = seq(-4,4,length=100)
  # Step 3B:  Use approx() to do 1D smoother of SPR_i for each row of Grid_j
  SPR_j = # ADD CODE
  # Step 3C:  translate Grid to 4 parameters and calculate FishLife probability for each
  Param_j = t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j) )
  colnames(Param_j) = names(Mean)
  Prob_j = dmvnorm( Param_j, mean=Mean, sigma=Cov )
  # Step 3D:  Calculate mean
  SPRmean = weighted.mean( SPR_j, w=Prob_j, na.rm=TRUE )
}
if( Dim==2 ){
  # Step 3A:  Define axis
  Grid_j = expand.grid( seq(-4,4,length=100), seq(-4,4,length=100) )
  # Step 3B:  Use akima:interp to do 2D smoother of SPR_i for each row of Grid_j
  SPR_j = # ADD CODE
  # Step 3C:  translate Grid to 4 parameters and calculate FishLife probability for each
  Param_j = t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j) )
  colnames(Param_j) = names(Mean)
  Prob_j = dmvnorm( Param_j, mean=Mean, sigma=Cov )
  # Step 3D:  Calculate mean
  SPRmean = weighted.mean( SPR_j, w=Prob_j, na.rm=TRUE )
}


