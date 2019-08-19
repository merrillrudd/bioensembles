rm(list=ls())

devtools::install_github("merrillrudd/LIME", dependencies=TRUE)
library(LIME)

library(FishLife)

library(mvQuad)
library(mvtnorm)

## directories
wd <- file.path("C:\\merrill\\bioensembles")
sim <- file.path(wd, "sim")
dir.create(sim, showWarnings=FALSE)


####################################
###  call FishLife               ###
####################################
## Puerto Rico hogfish Lachnolaimu maximus
sp <- Plot_taxa( Search_species(Genus="Lachnolaimus", Species="maximus")$match_taxonomy, mfrow=c(2,2) )

# png(file.path(figs, "FishLife_sim.png"), width=10, height=8, units="in", res=200)
# sp <- Plot_taxa( Search_species(Genus="Lutjanus", Species="guttatus")$match_taxonomy, mfrow=c(2,2) )
# dev.off()


## species-level
mp <- sp[[1]]$Mean_pred
cov <- sp[[1]]$Cov_pred


###############################
## 2 parameters
###############################

## choose parameters
param <- c("K","M")#,"Loo")

Mean <- mp[which(names(mp) %in% param)]
Cov <- cov[which(rownames(cov) %in% param), which(colnames(cov) %in% param)]
myGrid <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
rescale(myGrid, m=Mean, C=(Cov+t(Cov))/2, dec.type=1)

png(file.path(figs, "2param_dist_compareTrue.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(1,1))
plot(myGrid, lwd=3, xlab="log(k)", ylab="log(M)")
# points(x=log(0.21), y=log(0.43), col="blue", pch=19, cex=2)
dev.off()

nodes2 <- getNodes(myGrid)
weights2 <- getWeights(myGrid)

saveRDS(nodes2, file.path(sim, "nodes_2param_hogfish.rds"))
saveRDS(weights2, file.path(sim, "weights_2param_hogfish.rds"))


###############################
## 3 parameters
###############################

## choose parameters
param <- c("K","M","Loo")

Mean <- mp[which(names(mp) %in% param)]
Cov <- cov[which(rownames(cov) %in% param), which(colnames(cov) %in% param)]

############## NEW CODE
# Change #1 :  change dim=2 to dim=3 so that it matches exected syntax
# Change #2 : chanve Cov to (Cov+t(Cov))/2 so that its guarunteed to be symmetric

myGrid <- createNIGrid(dim=3, type="GHe", level=4,  ndConstruction="sparse")
# !(all(Cov == t(Cov)) & all(eigen(Cov)$values > 0))
rescale(myGrid, m=Mean, C=(Cov+t(Cov))/2, dec.type=1)

############## END NEW CODE

plot(myGrid)

nodes3 <- getNodes(myGrid)
weights3 <- getWeights(myGrid)



####################################
###  simulate truth              ###
####################################
## based on Lutjanus guttatus
plist <- create_lh_list(linf=64.6, vbk=0.21, t0=-0.01, 
						lwa=0.0245, lwb=2.79, 
						M=0.43, 
						M50=34, maturity_input="length",
						S50=c(30), S95=c(35), selex_input="length",
						SigmaF=0.2, SigmaR=0.737)

## generate data
itervec <- 1:10
data <- generate_data(modpath=sim, itervec=itervec, 
						Fdynamics="Endogenous", Rdynamics="AR", 
						lh=plist, 
						Nyears=20, Nyears_comp=20, comp_sample=200,
						init_depl=c(0.10,0.90), 
						seed=itervec,
						nburn=50,
						rewrite=TRUE)

##############################
## Run model - 2 parameters
##############################
## single simulated fleet
for(iter in 1:length(itervec)){
	iterpath <- file.path(sim, itervec[iter])

	data <- readRDS(file.path(iterpath, "True.rds"))
	data_list <- list("years"=data$years, "LF"=data$LF)

	res <- lapply(1:nrow(nodes2), function(x){
  		lhinp <- with(plist, 
  				create_lh_list(linf=linf, vbk=exp(nodes2[x,1]), t0=t0,
								lwa=lwa, lwb=lwb,
								M=exp(nodes2[x,2]),
								M50=M50, maturity_input="length",
								S50=SL50, S95=SL95, selex_input="length",
								SigmaF=SigmaF, SigmaR=SigmaR))	

  		out <- run_LIME(modpath=NULL, lh=lhinp, input_data=data_list, est_sigma="log_sigma_R", data_avail="LC", rewrite=FALSE)

		return(out)
	})
	saveRDS(res, file.path(iterpath, "res_2param.rds"))
}

spr2 <- sapply(1:length(itervec), function(y){
	res <- readRDS(file.path(sim, itervec[y], "res_2param.rds"))
	gen <- readRDS(file.path(sim, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$SPR_t[length(res[[x]]$Report$SPR_t)])
	}})
	out <- sum(est * weights2)
	true <- gen$SPR
	re <- (out-true)/true
	return(re)
})

f2 <- sapply(1:length(itervec), function(y){
	res <- readRDS(file.path(sim, itervec[y], "res_2param.rds"))
	gen <- readRDS(file.path(sim, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$F_t[length(res[[x]]$Report$F_t)])
	}})
	out <- sum(est * weights2)
	true <- gen$F_t[length(gen$F_t)]
	re <- (out-true)/true
	return(re)
})

sb2 <- sapply(1:length(itervec), function(y){
	res <- readRDS(file.path(sim, itervec[y], "res_2param.rds"))
	gen <- readRDS(file.path(sim, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$SB_t[length(res[[x]]$Report$SB_t)])
	}})
	out <- sum(est * weights2)
	true <- gen$SB_t[length(gen$SB_t)]
	re <- (out-true)/true
	return(re)
})



