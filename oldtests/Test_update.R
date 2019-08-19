
setwd("C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2017 -- Rudd data-poor ensemble/FishLife update")

library(FishLife)

#Ynew_ij = matrix( c("Loo"=log(91),"K"=log(0.08),"Winfinity"=NA,"tmax"=NA,"tm"=log(24.9),"M"=log(0.13),"Lm"=NA,"Temperature"=NA), nrow=1)
load( "Ynew_ij.Rdata" )
colnames(Ynew_ij) = c("Loo","K","Winfinity","tmax","tm","M","Lm","Temperature")

Cov_gjj = FishLife::database$Cov_gjj
obsCov_jj = FishLife::database$obsCov_jj

Update = Update_prediction( Taxon=Search_species(Genus="Lachnolaimus",Species="maximus",add_ancestors=FALSE)$match_taxonomy, Ynew_ij=Ynew_ij, Cov_gjj=Cov_gjj, obsCov_jj=obsCov_jj )
cbind( "Old"=exp(Update$predMean_j), "New"=exp(Update$updateMean_j) )
colMeans(exp(Ynew_ij), na.rm=TRUE)
cbind( "Old"=sqrt(diag(Update$predCov_jj)), "New"=sqrt(diag(Update$updateCov_jj)) )
