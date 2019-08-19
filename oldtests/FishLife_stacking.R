rm(list=ls())

###############################
### packages                ###
###############################

devtools::install_github("merrillrudd/LIME")
library(LIME)

library(FishLife)

devtools::install_github("adrianhordyk/LBSPR")
library(LBSPR)

library(mvQuad)
library(mvtnorm)

library(dplyr)
library(ggplot2)
library(foreach)
library(doParallel)

library(RColorBrewer)

require(StepwiseLH)

###############################################
###                 directories 			###
###############################################
wd <- file.path("C:\\merrill\\bioensembles")

sim <- file.path(wd, "sim")
dir.create(sim, showWarnings=FALSE)

figs <- file.path(wd, "figures")
dir.create(figs, showWarnings=FALSE)

funs <- list.files(file.path(wd, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(wd,"R", funs[x])))

res_dir <- file.path(sim, "results")
dir.create(res_dir, showWarnings=FALSE)

###############################################
###             call FishLife               ###
###############################################
load(file.path(wd, "Return.Rdata"))

## Puerto Rico hogfish Lachnolaimu maximus
sp <- Search_species(Genus="Lachnolaimus", Species="maximus", ParentChild_gz=Return$ParentChild_gz)$match_taxonomy
png(file.path(figs, "FishLife_sim.png"), width=10, height=8, units="in", res=200)
Plot_taxa( sp, mfrow=c(2,2) )
dev.off()

## find species in database
Which <- grep(sp[1], Return$ParentChild_gz[,"ChildName"])
Mean <- Return$beta_gv[Which,]
Cov <- Return$Cov_gvv[Which,,]

##################################################
###         nodes and weights - M 	   ###
###################################################


# choose parameters
param1 <- "M"

grid1 <- createNIGrid(dim=1, type="GHe", level=4, ndConstruction="sparse")
msub <- Mean[which(names(Mean) %in% param1)]
csub <- Cov[which(rownames(Cov) %in% param1), which(colnames(Cov) %in% param1)]
mvQuad::rescale(grid1, m=msub, C=(csub+t(csub))/2, dec.type=1)

nodes1 <- getNodes(grid1)
colnames(nodes1) <- names(msub)

weights1 <- getWeights(grid1)

##################################################
###         nodes and weights - M and K 	   ###
###################################################

# choose parameters
param2 <- c("K","M")

grid2 <- createNIGrid(dim=2, type="GHe", level=4, ndConstruction="sparse")
msub <- Mean[which(names(Mean) %in% param2)]
csub <- Cov[which(rownames(Cov) %in% param2), which(colnames(Cov) %in% param2)]
mvQuad::rescale(grid2, m=msub, C=(csub+t(csub))/2, dec.type=1)

nodes2 <- getNodes(grid2)
colnames(nodes2) <- names(msub)

weights2 <- getWeights(grid2)

uw <- unique(weights2)[order(unique(weights2))]
cols <- rev(topo.colors(length(uw)))
wcols <- sapply(1:nrow(weights2), function(x) cols[which(uw==weights2[x,1])])
plot(nodes2[,"K"], nodes2[,"M"], col=wcols, pch=19)

##################################################
###         nodes and weights - Lm and Linf	   ###
###################################################
# choose parameters
param2_v2 <- c("Lm","Loo")

grid2_v2 <- createNIGrid(dim=2, type="GHe", level=4, ndConstruction="sparse")
msub <- Mean[which(names(Mean) %in% param2_v2)]
csub <- Cov[which(rownames(Cov) %in% param2_v2), which(colnames(Cov) %in% param2_v2)]
mvQuad::rescale(grid2_v2, m=msub, C=(csub+t(csub))/2, dec.type=1)

nodes2_v2 <- getNodes(grid2_v2)
colnames(nodes2_v2) <- names(msub)

weights2_v2 <- getWeights(grid2_v2)

###############################################################
###         nodes and weights - M and K, Linf, Lmat 	   ###
################################################################

# choose parameters
param4 <- c("K","M","Lm","Loo")

grid4 <- createNIGrid(dim=4, type="GHe", level=4, ndConstruction="sparse")
msub <- Mean[which(names(Mean) %in% param4)]
csub <- Cov[which(rownames(Cov) %in% param4), which(colnames(Cov) %in% param4)]
mvQuad::rescale(grid4, m=msub, C=(csub+t(csub))/2, dec.type=1)

nodes4 <- getNodes(grid4)
colnames(nodes4) <- names(msub)

weights4 <- getWeights(grid4)


#################################################
##       Run models - 1 parameter         	  ###
#################################################

itervec <- 1:100

dir_M <- file.path(res_dir, "M")
dir.create(dir_M, showWarnings=FALSE)

	##################################
	## LBSPR - M
	##################################

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_M, 
					nodes=nodes1, 
					param=param1, 
					mean=Mean, 
					cov=Cov, 
					modname="M_Species", 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	mle_M_LBSPR <- 	stack_mle(savedir=dir_M,
					modname="M_Species",
					model="LBSPR",
					nodes=nodes1,
					weights=weights1,
					vals="SPR",
					itervec=itervec)

	##################################
	## LIME - M
	##################################

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_M, 
					nodes=nodes1, 
					param=param1, 
					mean=Mean, 
					cov=Cov, 
					modname="M_Species", 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	mle_M_LIME <- 	stack_mle(savedir=dir_M,
					modname="M_Species",
					model="LIME",
					nodes=nodes1,
					weights=weights1,
					vals=c("SPR","BB0"),
					itervec=itervec)

M_mle <- rbind.data.frame(mle_M_LBSPR, mle_M_LIME) %>% mutate("RunOption"="Stacking_M")

pspr <- ggplot(mle %>% filter(Value=="SPR")) + 
	geom_violin(aes(x=Model, y=RE, color=Model, fill=Model)) +
	geom_hline(aes(yintercept=0), lwd=1) +
	theme_lsd()

pbb <- ggplot(mle %>% filter(Model=="LIME")) + 
	geom_violin(aes(x=Value, y=RE, color=Value, fill=Value)) +
	geom_hline(aes(yintercept=0), lwd=1) +
	theme_lsd()


#################################################
##       Run models - M and K         	  ###
#################################################

itervec <- 1:100

dir_MK <- file.path(res_dir, "MK")
dir.create(dir_MK, showWarnings=FALSE)

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
		runstack(savedir=dir_MK, 
					nodes=nodes2, 
					param=param2, 
					mean=Mean, 
					cov=Cov, 
					modname="MK_Species", 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	mle_MK_LBSPR <- stack_mle(savedir=dir_MK,
					modname="MK_Species",
					model="LBSPR",
					nodes=nodes2,
					weights=weights2,
					vals="SPR",
					itervec=itervec)

	##################################
	## LIME
	##################################

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_MK, 
					nodes=nodes2, 
					param=param2, 
					mean=Mean, 
					cov=Cov, 
					modname="MK_Species", 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	mle_MK_LIME <- 	stack_mle(savedir=dir_MK,
					modname="MK_Species",
					model="LIME",
					nodes=nodes2,
					weights=weights2,
					vals=c("SPR","BB0"),
					itervec=itervec)

MK_mle <- rbind.data.frame(mle_MK_LBSPR, mle_MK_LIME) %>% mutate("RunOption"="Stacking_M")

pspr <- ggplot(MK_mle %>% filter(Value=="SPR")) + 
	geom_violin(aes(x=Model, y=RE, color=Model, fill=Model)) +
	geom_hline(aes(yintercept=0), lwd=1) +
	theme_lsd()

pbb <- ggplot(MK_mle %>% filter(Model=="LIME")) + 
	geom_violin(aes(x=Value, y=RE, color=Value, fill=Value)) +
	geom_hline(aes(yintercept=0), lwd=1) +
	theme_lsd()


#################################################
##       Run models - Lm and Linf        	  ###
#################################################

itervec <- 1:100

dir_LmLinf <- file.path(res_dir, "LmLinf")
dir.create(dir_LmLinf, showWarnings=FALSE)

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
		runstack(savedir=dir_LmLinf, 
					nodes=nodes2_v2, 
					param=param2_v2, 
					mean=Mean, 
					cov=Cov, 
					modname="LmLinf_Species", 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	mle_LmLinf_LBSPR <- stack_mle(savedir=dir_LmLinf,
					modname="LmLinf_Species",
					model="LBSPR",
					nodes=nodes2_v2,
					weights=weights2_v2,
					vals="SPR",
					itervec=itervec)

	##################################
	## LIME
	##################################

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_LmLinf, 
					nodes=nodes2_v2, 
					param=param2_v2, 
					mean=Mean, 
					cov=Cov, 
					modname="LmLinf_Species", 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	mle_LmLinf_LIME <- 	stack_mle(savedir=dir_LmLinf,
					modname="LmLinf_Species",
					model="LIME",
					nodes=nodes2_v2,
					weights=weights2_v2,
					vals=c("SPR","BB0"),
					itervec=itervec)

LmLinf_mle <- rbind.data.frame(mle_LmLinf_LBSPR, mle_LmLinf_LIME) %>% mutate("RunOption"="Stacking_M")

pspr <- ggplot(LmLinf_mle %>% filter(Value=="SPR")) + 
	geom_violin(aes(x=Model, y=RE, color=Model, fill=Model)) +
	geom_hline(aes(yintercept=0), lwd=1) +
	theme_lsd()

pbb <- ggplot(LmLinf_mle %>% filter(Model=="LIME")) + 
	geom_violin(aes(x=Value, y=RE, color=Value, fill=Value)) +
	geom_hline(aes(yintercept=0), lwd=1) +
	theme_lsd()


#################################################
##       Run models - 4 parameters        	  ###
#################################################

itervec <- 1:100

dir_4param <- file.path(res_dir, "4param")
dir.create(dir_4param, showWarnings=FALSE)

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
		runstack(savedir=dir_4param, 
					nodes=nodes4, 
					param=param4, 
					mean=Mean, 
					cov=Cov, 
					modname="4param_Species", 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	mle_4_LBSPR <- stack_mle(savedir=dir_4param,
					modname="4param_Species",
					model="LBSPR",
					nodes=nodes4,
					weights=weights4,
					vals="SPR",
					itervec=itervec)

	savedir=dir_4param
					modname="4param_Species"
					model="LBSPR"
					nodes=nodes4
					weights=weights4
					vals="SPR"
					itervec=itervec

	##################################
	## LIME
	##################################

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_4param, 
					nodes=nodes4, 
					param=param4, 
					mean=Mean, 
					cov=Cov, 
					modname="4param_Species", 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	mle_4_LIME <- 	stack_mle(savedir=dir_4param,
					modname="4param_Species",
					model="LIME",
					nodes=nodes4,
					weights=weights4,
					vals=c("SPR","BB0"),
					itervec=itervec)

4_mle <- rbind.data.frame(mle_4_LBSPR, mle_4_LIME) %>% mutate("RunOption"="Stacking_M")

pspr <- ggplot(mle_4_LBSPR %>% filter(Value=="SPR")) + 
	geom_violin(aes(x=Model, y=RE, color=Model, fill=Model)) +
	geom_hline(aes(yintercept=0), lwd=1) +
	theme_lsd()

pbb <- ggplot(4_mle %>% filter(Model=="LIME")) + 
	geom_violin(aes(x=Value, y=RE, color=Value, fill=Value)) +
	geom_hline(aes(yintercept=0), lwd=1) +
	theme_lsd()





















#################################################
##       Run models - 2 parameters         	  ###
#################################################


dir_2param <- file.path(res_dir, "param2_LBSPR")
dir.create(dir_2param, showWarnings=FALSE)

	##################################
	## species - LBSPR - 2parameters
	##################################
	type <- "Species"	

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_2param, 
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)



	## calculate densities - LBSPR species
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param, 
					modname=paste0("2param_",type), 
					model="LBSPR", 
					iter=loop, 
					vals="SPR", 
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

	####################
	## genus - LBSPR
	####################
	type <- "Genus"		
	## run LBSPR
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)	

	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_2param, 
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	## calculate densities - LBSPR genus
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param, 
					modname=paste0("2param_",type), 
					model="LBSPR", 
					iter=loop, 
					vals="SPR",
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

	####################
	## Family - LBSPR
	####################
	type <- "Family"		
	## run LBSPR
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_2param, 
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	## calculate densities - LBSPR genus
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param, 
					modname=paste0("2param_",type), 
					model="LBSPR", 
					iter=loop, 
					vals="SPR",
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)


######-------------------- equilibrium, 2 parameters -----------------####

dir_2param2 <- file.path(res_dir, "param2_LIME")
dir.create(dir_2param2, showWarnings=FALSE)

	##################################
	## species - LIME - 4parameters
	##################################
	type <- "Species"	

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
	# for(loop in itervec){
		runstack(savedir=dir_2param2, 
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1,
					LFdist=0)
	# }
	end <- Sys.time() - start		
	stopCluster(cl)

	## calculate densities - LIME species
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param2, 
					modname=paste0("2param_",type), 
					model="LIME", 
					iter=loop, 
					vals=c("Depletion","SPR"), 
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)


	####################
	## genus - LIME
	####################
	type <- "Genus"		
	## run LIME
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR', 'mvtnorm')) %dopar% 
		runstack(savedir=dir_2param2, 
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1,
					LFdist=0)
	end <- Sys.time() - start		
	stopCluster(cl)


	## calculate densities - LIME genus
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param2, 
					modname=paste0("2param_",type), 
					model="LIME", 
					iter=loop, 
					vals=c("Depletion","SPR"), 
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

	####################
	## family - LIME
	####################
	type <- "Family"		
	## run LBSPR
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR', 'mvtnorm')) %dopar% 
		runstack(savedir=dir_2param2, 
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1,
					LFdist=0)
	end <- Sys.time() - start		
	stopCluster(cl)


	## calculate densities - LIME family
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param2, 
					modname=paste0("2param_",type), 
					model="LIME", 
					iter=loop, 
					vals=c("Depletion","SPR"), 
					nodes=nodes2[[type]], 
					param=param2, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

#################################################
##       Run models - 2 parameters - invariant         	  ###
#################################################


dir_2param_BH <- file.path(res_dir, "paramBH_LBSPR")
dir.create(dir_2param_BH, showWarnings=FALSE)

	##################################
	## species - LBSPR - 2parameters
	##################################
	type <- "Species"	

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_2param_BH, 
					nodes=nodes2_v2[[type]], 
					param=param2_v2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2paramv2_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	## calculate densities - LBSPR species
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param_BH, 
					modname=paste0("2paramv2_",type), 
					model="LBSPR", 
					iter=loop, 
					vals="SPR", 
					nodes=nodes2_v2[[type]], 
					param=param2_v2, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

	####################
	## genus - LBSPR
	####################
	type <- "Genus"		
	## run LBSPR
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)	

	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_2param_BH, 
					nodes=nodes2_v2[[type]], 
					param=param2_v2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2paramv2_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	## calculate densities - LBSPR genus
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param_BH, 
					modname=paste0("2paramv2_",type), 
					model="LBSPR", 
					iter=loop, 
					vals="SPR",
					nodes=nodes2_v2[[type]],
					param=param2_v2,  
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

	####################
	## Family - LBSPR
	####################
	type <- "Family"		
	## run LBSPR
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_2param_BH, 
					nodes=nodes2_v2[[type]], 
					param=param2_v2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2paramv2_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	## calculate densities - LBSPR genus
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param_BH, 
					modname=paste0("2paramv2_",type), 
					model="LBSPR", 
					iter=loop, 
					vals="SPR",
					nodes=nodes2_v2[[type]], 
					param=param2_v2, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)


######-------------------- equilibrium, 2 parameters -----------------####

dir_2param_BH2 <- file.path(res_dir, "paramBH_LIME")
dir.create(dir_2param_BH2, showWarnings=FALSE)

	##################################
	## species - LIME - 4parameters
	##################################
	type <- "Species"	

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
	# for(loop in itervec){
		runstack(savedir=dir_2param_BH2, 
					nodes=nodes2_v2[[type]], 
					param=param2_v2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2paramv2_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1,
					LFdist=0)
	# }
	end <- Sys.time() - start		
	stopCluster(cl)

	## calculate densities - LIME species
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param_BH2, 
					modname=paste0("2paramv2_",type), 
					model="LIME", 
					iter=loop, 
					vals=c("Depletion","SPR"), 
					nodes=nodes2_v2[[type]], 
					param=param2_v2,
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)


	####################
	## genus - LIME
	####################
	type <- "Genus"		
	## run LIME
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR', 'mvtnorm')) %dopar% 
		runstack(savedir=dir_2param_BH2, 
					nodes=nodes2_v2[[type]], 
					param=param2_v2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2paramv2_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1,
					LFdist=0)
	end <- Sys.time() - start		
	stopCluster(cl)


	## calculate densities - LIME genus
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param_BH2, 
					modname=paste0("2paramv2_",type), 
					model="LIME", 
					iter=loop, 
					vals=c("Depletion","SPR"), 
					nodes=nodes2_v2[[type]], 
					param=param2_v2,
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

	####################
	## family - LIME
	####################
	type <- "Family"		
	## run LBSPR
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR', 'mvtnorm')) %dopar% 
		runstack(savedir=dir_2param_BH2, 
					nodes=nodes2_v2[[type]], 
					param=param2_v2, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("2paramv2_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1,
					LFdist=0)
	end <- Sys.time() - start		
	stopCluster(cl)


	## calculate densities - LIME family
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_2param_BH2, 
					modname=paste0("2paramv2_",type), 
					model="LIME", 
					iter=loop, 
					vals=c("Depletion","SPR"), 
					nodes=nodes2_v2[[type]], 
					param=param2_v2,
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)


#################################################
##  	Run models  - 4 parameters             ###
#################################################

######-------------------- equilibrium, 4 parameters -----------------####

dir_4param <- file.path(res_dir, "param4_LBSPR")
dir.create(dir_4param, showWarnings=FALSE)

	##################################
	## species - LBSPR - 4parameters
	##################################
	type <- "Species"	

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_4param, 
					nodes=nodes4[[type]], 
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("4param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	## calculate densities - LBSPR species
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_4param, 
					modname=paste0("4param_",type), 
					model="LBSPR", 
					iter=loop, 
					vals="SPR", 
					nodes=nodes4[[type]],
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

	####################
	## genus - LBSPR
	####################
	type <- "Genus"		
	## run LBSPR
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_4param, 
					nodes=nodes4[[type]], 
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("4param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)


	## calculate densities - LBSPR genus
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_4param, 
					modname=paste0("4param_",type), 
					model="LBSPR", 
					iter=loop, 
					vals="SPR",
					nodes=nodes4[[type]],
					param=param4,  
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

	####################
	## Family - LBSPR
	####################
	type <- "Family"		
	## run LBSPR
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
		runstack(savedir=dir_4param, 
					nodes=nodes4[[type]], 
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("4param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LBSPR", 
					sim_model="LIME",
					Nyears=1)
	end <- Sys.time() - start		
	stopCluster(cl)

	## calculate densities - LBSPR genus
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_4param, 
					modname=paste0("4param_",type), 
					model="LBSPR", 
					iter=loop, 
					vals="SPR",
					nodes=nodes4[[type]], 
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)





######-------------------- equilibrium, 4 parameters -----------------####

dir_4param2 <- file.path(res_dir, "param4_LIME")
dir.create(dir_4param2, showWarnings=FALSE)

	##################################
	## species - LIME - 4parameters
	##################################
	type <- "Species"	

	## Simulate true data with species-level information
	## Run all models with species-level distributions
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
	# for(loop in itervec){
		runstack(savedir=dir_4param2, 
					nodes=nodes4[[type]], 
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("4param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1,
					LFdist=0)
	# }
	end <- Sys.time() - start		
	stopCluster(cl)

	## calculate densities - LIME species
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_4param2, 
					modname=paste0("4param_",type), 
					model="LIME", 
					iter=loop, 
					vals=c("Depletion","SPR"), 
					nodes=nodes4[[type]], 
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

	####################
	## genus - LIME
	####################
	type <- "Genus"		
	## run LIME
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR', 'mvtnorm')) %dopar% 
		runstack(savedir=dir_4param2, 
					nodes=nodes4[[type]], 
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("4param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1,
					LFdist=0)
	end <- Sys.time() - start		
	stopCluster(cl)


	## calculate densities - LIME genus
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_4param2, 
					modname=paste0("4param_",type), 
					model="LIME", 
					iter=loop, 
					vals=c("Depletion","SPR"), 
					nodes=nodes4[[type]], 
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)

	####################
	## family - LIME
	####################
	type <- "Family"		
	## run LIME
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)			
	start <- Sys.time()
	foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR', 'mvtnorm')) %dopar% 
		runstack(savedir=dir_4param2, 
					nodes=nodes4[[type]], 
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]], 
					modname=paste0("4param_",type), 
					simulation=TRUE, 
					iter=loop, 
					seed=loop, 
					Fscenario="equil", 
					rewrite=FALSE, 
					model="LIME", 
					sim_model="LIME",
					Nyears=1,
					LFdist=0)
	end <- Sys.time() - start		
	stopCluster(cl)


	## calculate densities - LIME family
	ncores <- 5
	cl <- makeCluster(ncores)
	registerDoParallel(cl)		
	foreach(loop=itervec, .packages=c('dplyr', 'mvtnorm', 'Hmisc')) %dopar% 
		stackdensity(savedir=dir_4param2, 
					modname=paste0("4param_",type), 
					model="LIME", 
					iter=loop, 
					vals=c("Depletion","SPR"), 
					nodes=nodes4[[type]], 
					param=param4, 
					mean=Mean[[type]], 
					cov=Cov[[type]],
					rewrite=TRUE,
					rewrite_summary=TRUE)
	stopCluster(cl)





################################################
### 		Summarize Results 				####
################################################
		combine <- c("1param_Species", "1param_Genus", "1param_Family")

		## bring together summaries across iterations
		summary1 <- lapply(itervec, function(x){
			find <- lapply(1:length(combine), function(y){
				sub <- readRDS(file.path(dir_1param, itervec[x], paste0("results_summary_", combine[y], "_LBSPR.rds")))
				lab <- strsplit(combine[y],"_")
				param <- lab[[1]][1]
				taxa <- lab[[1]][2]
				combine1 <- do.call(rbind, sub) %>% mutate("Iteration"=itervec[x]) %>% mutate("Param"="M") %>% mutate("Taxon"=taxa)
				return(combine1)
			})
			combine2 <- do.call(rbind, find)

			return(combine2)
		})	
		summary1 <- do.call(rbind, summary1) %>% mutate("Model"="LBSPR") %>% mutate("Scenario"="Equilibrium")

		## bring together summaries across iterations
		summary12 <- lapply(itervec, function(x){
			find <- lapply(1:length(combine), function(y){
				if(file.exists(file.path(dir_1param2, itervec[x], paste0("results_summary_", combine[y], "_LIME.rds")))){
					sub <- readRDS(file.path(dir_1param2, itervec[x], paste0("results_summary_", combine[y], "_LIME.rds")))
					lab <- strsplit(combine[y],"_")
					param <- lab[[1]][1]
					taxa <- lab[[1]][2]
					combine1 <- do.call(rbind, sub) %>% mutate("Iteration"=itervec[x]) %>% mutate("Param"="M") %>% mutate("Taxon"=taxa)
					return(combine1)
				}
			})
			combine2 <- do.call(rbind, find)
			return(combine2)
		})	
		summary12 <- do.call(rbind, summary12) %>% mutate("Model"="LIME") %>% mutate("Scenario"="Equilibrium")

		combine <- c("2param_Species", "2param_Genus", "2param_Family")

		## bring together summaries across iterations
		summary2 <- lapply(itervec, function(x){
			find <- lapply(1:length(combine), function(y){
				sub <- readRDS(file.path(dir_2param, itervec[x], paste0("results_summary_", combine[y], "_LBSPR.rds")))
				lab <- strsplit(combine[y],"_")
				param <- lab[[1]][1]
				taxa <- lab[[1]][2]
				combine1 <- do.call(rbind, sub) %>% mutate("Iteration"=itervec[x]) %>% mutate("Param"="MK") %>% mutate("Taxon"=taxa)
				return(combine1)
			})
			combine2 <- do.call(rbind, find)

			return(combine2)
		})	
		summary2 <- do.call(rbind, summary2) %>% mutate("Model"="LBSPR") %>% mutate("Scenario"="Equilibrium")

		## bring together summaries across iterations
		summary22 <- lapply(itervec, function(x){
			find <- lapply(1:length(combine), function(y){
				if(file.exists(file.path(dir_2param2, itervec[x], paste0("results_summary_", combine[y], "_LIME.rds")))){
					sub <- readRDS(file.path(dir_2param2, itervec[x], paste0("results_summary_", combine[y], "_LIME.rds")))
					lab <- strsplit(combine[y],"_")
					param <- lab[[1]][1]
					taxa <- lab[[1]][2]
					combine1 <- do.call(rbind, sub) %>% mutate("Iteration"=itervec[x]) %>% mutate("Param"="MK") %>% mutate("Taxon"=taxa)
					return(combine1)
				}
			})
			combine2 <- do.call(rbind, find)
			return(combine2)
		})	
		summary22 <- do.call(rbind, summary22) %>% mutate("Model"="LIME") %>% mutate("Scenario"="Equilibrium")


		combine <- c("2paramv2_Species", "2paramv2_Genus", "2paramv2_Family")

		## bring together summaries across iterations
		summaryBH <- lapply(itervec, function(x){
			find <- lapply(1:length(combine), function(y){
				sub <- readRDS(file.path(dir_2param_BH, itervec[x], paste0("results_summary_", combine[y], "_LBSPR.rds")))
				lab <- strsplit(combine[y],"_")
				param <- lab[[1]][1]
				taxa <- lab[[1]][2]
				combine1 <- do.call(rbind, sub) %>% mutate("Iteration"=itervec[x]) %>% mutate("Param"="MLinf") %>% mutate("Taxon"=taxa)
				return(combine1)
			})
			combine2 <- do.call(rbind, find)

			return(combine2)
		})	
		summaryBH <- do.call(rbind, summaryBH) %>% mutate("Model"="LBSPR") %>% mutate("Scenario"="Equilibrium")

		## bring together summaries across iterations
		summaryBH2 <- lapply(itervec, function(x){
			find <- lapply(1:length(combine), function(y){
				if(file.exists(file.path(dir_2param_BH2, itervec[x], paste0("results_summary_", combine[y], "_LIME.rds")))){
					sub <- readRDS(file.path(dir_2param_BH2, itervec[x], paste0("results_summary_", combine[y], "_LIME.rds")))
					lab <- strsplit(combine[y],"_")
					param <- lab[[1]][1]
					taxa <- lab[[1]][2]
					combine1 <- do.call(rbind, sub) %>% mutate("Iteration"=itervec[x]) %>% mutate("Param"="MLinf") %>% mutate("Taxon"=taxa)
					return(combine1)
				}
			})
			combine2 <- do.call(rbind, find)
			return(combine2)
		})	
		summaryBH2 <- do.call(rbind, summaryBH2) %>% mutate("Model"="LIME") %>% mutate("Scenario"="Equilibrium")


		combine <- c("4param_Species", "4param_Genus", "4param_Family") 

		## bring together summaries across iterations
		summary4 <- lapply(itervec, function(x){
			find <- lapply(1:length(combine), function(y){
				sub <- readRDS(file.path(dir_4param, itervec[x], paste0("results_summary_", combine[y], "_LBSPR.rds")))
				lab <- strsplit(combine[y],"_")
				param <- lab[[1]][1]
				taxa <- lab[[1]][2]
				combine1 <- do.call(rbind, sub) %>% mutate("Iteration"=itervec[x]) %>% mutate("Param"=param) %>% mutate("Taxon"=taxa)
				return(combine1)
			})
			combine2 <- do.call(rbind, find)

			return(combine2)
		})	
		summary4 <- do.call(rbind, summary4) %>% mutate("Model"="LBSPR") %>% mutate("Scenario"="Equilibrium")


		## bring together summaries across iterations
		summary42 <- lapply(itervec, function(x){
			find <- lapply(1:length(combine), function(y){
				if(file.exists(file.path(dir_4param2, itervec[x], paste0("results_summary_", combine[y], "_LIME.rds")))){
					sub <- readRDS(file.path(dir_4param2, itervec[x], paste0("results_summary_", combine[y], "_LIME.rds")))
					lab <- strsplit(combine[y],"_")
					param <- lab[[1]][1]
					taxa <- lab[[1]][2]
					combine1 <- do.call(rbind, sub) %>% mutate("Iteration"=itervec[x]) %>% mutate("Param"=param) %>% mutate("Taxon"=taxa)
					return(combine1)
				}
			})
			combine2 <- do.call(rbind, find)
			return(combine2)
		})	
		summary42 <- do.call(rbind, summary42) %>% mutate("Model"="LIME") %>% mutate("Scenario"="Equilibrium")


		summary <- rbind(summary1, summary12, summary2, summary22, summaryBH, summaryBH2, summary4, summary42) 

		# ## filter out unnecessary best cases
		best <- unique(summary %>% filter(Type=="BestCase") %>% filter(Param=="4param") %>% select(Year, Variable, Type, Iteration, Model, Scenario, Cover50, RE)) %>%  mutate(Param=NA) %>% mutate(Taxon="True")
		means <- unique(summary %>% filter(Type=="Means") %>% filter(Param=="4param") %>% select(Year, Variable, Type, Iteration, Model, Scenario, Taxon, Cover50, RE)) %>% mutate(Param=NA)
		stack <- unique(summary %>% filter(Type=="Stacking") %>% select(Year, Variable, Type, Param, Iteration, Model, Scenario, Taxon, Cover50, RE))

		summary_adj <- rbind(best, means, stack) %>% mutate(RunOption=ifelse(is.na(Param), as.character(Type), paste(Type, Param, sep="_")))

		re <- summary_adj %>% filter(Year==max(Year)) %>% 
				group_by(Variable, Type, Param, Taxon, Scenario, Model) %>%
				select(Variable, Type, Param, Taxon, Scenario, Model, RE) %>%
				summarise_all(funs(mre=median(.,na.rm=TRUE), mare=median(abs(.),na.rm=TRUE)))

		cover <- summary_adj %>% filter(Year==max(Year)) %>%
				group_by(Variable, Type, Param, Taxon, Scenario, Model, RunOption) %>%
				select(Variable, Type, Param, Taxon, Scenario, Model, RunOption, Cover50) %>%
				summarise_all(funs(Coverage=sum(., na.rm=TRUE), Converge=length(which(is.na(Cover50)==FALSE)))) %>%
				mutate(PropCover = Coverage / Converge)

		results <- full_join(re, cover)
		write.csv(results, file.path(res_dir,"Ensemble_results.csv"))


		cols <- brewer.pal(4, "Set1")
		incols <- c(cols[3:1], cols[4])
		p <- ggplot(summary_adj %>% filter(Variable=="SPR")) +
			geom_violin(aes(x=RunOption, y=RE, color=Taxon, fill=Taxon), scale="width") +
			geom_hline(aes(yintercept=0), lwd=1.5) + 
			facet_grid(Model~.) +
			theme_lsd() +
			ylab("Relative error") + xlab("Modeling approach") +
			scale_fill_manual(values=incols) + scale_color_manual(values=incols)
		ggsave(file.path(figs, "SPR_RE.png"), p, width=10, height=8)

		c <- ggplot(results %>% filter(Variable=="SPR")) +
			geom_hline(aes(yintercept=0.5), lwd=1) +
			geom_point(aes(x=RunOption, y=PropCover, color=Taxon, fill=Taxon, shape=Taxon), cex=6) +
			facet_grid(Model~.) +
			theme_lsd() +
			xlab("Modeling approach") +
			scale_fill_manual(values=incols) + scale_color_manual(values=incols)
		ggsave(file.path(figs, "SPR_Coverage.png"), c, width=10, height=8)

		p2 <- ggplot(summary_adj %>% filter(Variable=="Depletion")) +
			geom_violin(aes(x=RunOption, y=RE, color=Taxon, fill=Taxon), scale="width") +
			geom_hline(aes(yintercept=0), lwd=1.5) + 
			facet_grid(Model~.) +
			theme_lsd() +
			ylab("Relative error") + xlab("Modeling approach") +
			scale_fill_manual(values=incols) + scale_color_manual(values=incols)
		ggsave(file.path(figs, "Depletion_RE.png"), p2)

		c2 <- ggplot(results %>% filter(Variable=="Depletion")) +
			geom_hline(aes(yintercept=0.5), lwd=1.5) +
			geom_point(aes(x=RunOption, y=PropCover, color=Taxon, fill=Taxon, shape=Taxon), cex=5) +
			facet_grid(Model~.) +
			theme_lsd() +
			xlab("Modeling approach") +
			scale_fill_manual(values=incols) + scale_color_manual(values=incols)
		ggsave(file.path(figs, "Depletion_Coverage.png"), c2)



### weights of PDF
Mean2 <- Mean[[1]][which(names(Mean[[1]]) %in% param2)]
Cov2 <- Cov[[1]][which(rownames(Cov[[1]]) %in% param2), which(colnames(Cov[[1]]) %in% param2)]
logdens_i <- dmvnorm(nodes2[[1]], mean=Mean2, sigma=Cov2, log=TRUE)

par(mfrow=c(1,1))
## interpolate function
dens <- readRDS(file.path(dir_2param, 5,"results_density_2param_Species_LBSPR.rds"))
dseq <- as.numeric(rownames(dens[[1]][[1]]))
stack <- readRDS(file.path(dir_2param, 5,"results_stackinterpolation_2param_Species_LBSPR.rds"))
stackdens <- readRDS(file.path(dir_2param, 5,"results_stackinterpolation_dens_2param_Species_LBSPR.rds"))
mres <- readRDS(file.path(dir_2param, 5,"2param_Species_res_Means_LBSPR.rds"))
true <- readRDS(file.path(dir_2param, 5,'True.rds'))


## interpolate PDF
interp_logdens = akima::interp( x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], z=logdens_i, xo=seq(min(nodes2[[1]][,"K"]),max(nodes2[[1]][,"K"]),length.out=100), yo=seq(min(nodes2[[1]][,"M"]),max(nodes2[[1]][,"M"]),length.out=100), nx=100, ny=100, duplicate= "median")
si.zmin <- min(interp_logdens$z, na.rm=TRUE)
si.zmax <- max(interp_logdens$z, na.rm=TRUE)
breaks <- pretty(c(si.zmin, si.zmax),10)
colors <- rev(topo.colors(length(breaks)-1))

png(file.path(figs, "PDF_weights.png"), width=5, height=5.5, units="in", res=200)
image(interp_logdens, breaks=breaks, col=colors, xlab="log(K)", ylab="log(M)", cex.axis=1.3, cex.lab=1.3)
## show points on top of PDF
for(i in c("Species")){
	points(x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], cex=2, xpd=NA, lwd=2, pch=19, col="#AAAAAA90")
	points(log(true$vbk), log(true$M), pch=4, cex=2, col="black", lwd=6)
}
dev.off()

png(file.path(figs, "Function_density_atTruth.png"), width=5, height=5.5, units="in", res=200)
est <- round(stack[[1]]*1000)
tval <- round(true$SPR*1000)
interp_h <- akima::interp(x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], z=dens[[1]][[1]][tval,], xo=seq(min(nodes2[[1]][,"K"]),max(nodes2[[1]][,"K"]),length.out=100), yo=seq(min(nodes2[[1]][,"M"]),max(nodes2[[1]][,"M"]),length.out=100), nx=100, ny=100, duplicate="median")
image(interp_h, col=colors, xlab="log(K)", ylab="log(M)", cex.axis=1.3, cex.lab=1.3)
for(i in c("Species")){
	points(x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], cex=2, xpd=NA, lwd=2, pch=19, col="#AAAAAA90")
	points(log(true$vbk), log(true$M), pch=4, cex=2, col="black", lwd=6)
}
dev.off()

png(file.path(figs, "Weighted_density_atTruth.png"), width=5, height=5.5, units="in", res=200)
cell_area = mean(diff(interp_logdens$x)) * mean(diff(interp_logdens$y))
sumdens = sum( exp(interp_logdens$z), na.rm=TRUE ) * cell_area	
step1 <- exp(interp_logdens$z) * interp_h$z
final <- step1 / sumdens	
image(final, col=colors, axes=F, ann=F)

par(new=TRUE)
for(i in c("Species")){
	plot(x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], cex=2, xpd=NA, lwd=2, pch=19, col="#AAAAAA90", cex.axis=1.3, cex.lab=1.3, xlab="log(K)", ylab="log(M)")
	points(log(true$vbk), log(true$M), pch=4, cex=2, col="black", lwd=6)
}
dev.off()





png(file.path(figs, "Function_density_examples.png"), width=6, height=6, units="in", res=200)
par(mfcol=c(3,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
tval <- round(true$SPR*1000)
low <- round(tval*0.8)
high <- min(round(tval*1.2),1001)
est <- round(stack[[1]]*1000)

pseq <- seq(101,by=100,length.out=9)

for(j in 1:length(pseq)){
	interp_h <- akima::interp(x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], z=dens[[1]][[1]][pseq[j],], xo=seq(min(nodes2[[1]][,"K"]),max(nodes2[[1]][,"K"]),length.out=100), yo=seq(min(nodes2[[1]][,"M"]),max(nodes2[[1]][,"M"]),length.out=100), nx=100, ny=100, duplicate= "median")

	ylim <- c(min(interp_h$y),max(interp_h$y)+0.2*abs(max(interp_h$y)))
	
	image(interp_h, col=colors, xaxt="n", yaxt="n", ylim=ylim)
	if(j %% 3 == 0) axis(1, cex.axis=1.5)
	if(j %in% 1:3) axis(2, cex.axis=1.5)
	## show points on top of PDF
	for(i in c("Species")){
	points(x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], cex=1.5, xpd=NA, lwd=2, pch=19, col="#AAAAAA90")
	points(log(true$vbk), log(true$M), pch=4, cex=2, col="black", lwd=6)
	}	
	mtext(side=3, line=-1.8, paste0("SPR=",dseq[pseq[j]]), cex=1.3)
}
mtext(side=1, "log(K)", line=3, cex=1.3, outer=TRUE)
mtext(side=2, "log(M)", line=3, cex=1.3, outer=TRUE)
dev.off()

png(file.path(figs, "Weighted_density_examples.png"), width=6, height=6, units="in", res=200)
par(mfcol=c(3,3), mar=c(0,0,0,0), omi=c(1,1,1,1))


for(j in 1:length(pseq)){
	interp_h <- akima::interp(x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], z=dens[[1]][[1]][pseq[j],], xo=seq(min(nodes2[[1]][,"K"]),max(nodes2[[1]][,"K"]),length.out=100), yo=seq(min(nodes2[[1]][,"M"]),max(nodes2[[1]][,"M"]),length.out=100), nx=100, ny=100, duplicate= "median")
	# interp_h <- akima::interp(x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], z=dens[[1]][[1]][est,], xo=seq(min(nodes2[[1]][,"K"]),max(nodes2[[1]][,"K"]),length.out=100), yo=seq(min(nodes2[[1]][,"M"]),max(nodes2[[1]][,"M"]),length.out=100), nx=100, ny=100, duplicate= "median")

	cell_area = mean(diff(interp_logdens$x)) * mean(diff(interp_logdens$y))
	sumdens = sum( exp(interp_logdens$z), na.rm=TRUE ) * cell_area	
	step1 <- exp(interp_logdens$z) * interp_h$z
	final <- step1 / sumdens	
	image(final, col=colors, xaxt="n", yaxt="n", ylim=c(0,1.25))
	# show points on top of PDF
	par(new=TRUE)

	plot(x=1, y=1, type="n", xlim=c(min(interp_h$x),max(interp_h$x)), ylim=c(min(interp_h$y),max(interp_h$y)+0.2*abs(max(interp_h$y))), xaxt="n", yaxt="n")
	for(i in c("Species")){
		points(x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], cex=1.5, xpd=NA, lwd=2, pch=19, col="#AAAAAA90")
		points(log(true$vbk), log(true$M), pch=4, cex=2, col="black", lwd=6)
	}
	if(j %% 3 ==0) axis(1, cex.axis=1.5)
	if(j %in% 1:3) axis(2, cex.axis=1.5)
	mtext(side=3, line=-1.8, paste0("SPR=",dseq[pseq[j]]), cex=1.3)
}
	mtext(side=1, "log(K)", line=3, cex=1.3, outer=TRUE)
	mtext(side=2, "log(M)", line=3, cex=1.3, outer=TRUE)
dev.off()



png(file.path(figs, "Node_densities_example.png"), width=10, height=4, units="in", res=200)
par(mfrow=c(1,1), mar=c(4,4,4,4), omi=c(0.1,0.1,0,0))
plot(x=1, y=1, xlim=c(0,1), ylim=c(0,max(dens[[1]][[1]])*1.1), type="n", xlab="SPR", ylab="Density", xaxs="i", yaxs="i", cex.axis=1.3, cex.lab=1.3)

## weights of each node
node_order <- order(logdens_i)
wtcol <- rev(topo.colors(n=length(node_order), alpha=0.5))
for(i in 1:length(node_order)){
	polygon(x=c(dseq, rev(dseq)), y=c(rep(0,length(dseq)), rev(dens[[1]][[1]][,node_order[i]])), col=wtcol[i])
}
dev.off()

png(file.path(figs, "Predictive_stacking_density_example.png"), width=8, height=4, units="in", res=200)
plot(x=1, y=1, type="n", xlim=c(0,1), ylim=c(0, max(stackdens[[1]][[1]])*1.3), xlab="SPR", ylab="Density", xaxs="i", yaxs="i", cex.axis=1.3, cex.lab=1.3)
polygon(x=c(dseq, rev(dseq)), y=c(rep(0, length(dseq)), rev(stackdens[[1]][[1]])), lwd=3, col="#AAAAAADD")
abline(v=true$SPR, lwd=3, col="red")
abline(v=stack[["SPR"]], lwd=3)
dev.off()



par(mfcol=c(5,5), mar=c(0.1,0.1,0.1,0.1))
find <- round(seq(10,nrow(dens[[1]][[1]]),length.out=25))
seq <- as.numeric(rownames(dens[[1]][[1]]))

for(i in find){
interp_h <- akima::interp(x=nodes2[[1]][,"K"], y=nodes2[[1]][,"M"], z=dens[[1]][[1]][i,], xo=seq(min(nodes2[[1]][,"K"]),max(nodes2[[1]][,"K"]),length.out=100), yo=seq(min(nodes2[[1]][,"M"]),max(nodes2[[1]][,"M"]),length.out=100), nx=100, ny=100, duplicate= "median")
si.zmin2 <- min(interp_h$z, na.rm=TRUE)
si.zmax2 <- max(interp_h$z, na.rm=TRUE)
breaks2 <- pretty(c(si.zmin2, si.zmax2),10)
colors2 <- rev(heat.colors(length(breaks2)-1))

cell_area = mean(diff(interp_logdens$x)) * mean(diff(interp_logdens$y))
sumdens = sum( exp(interp_logdens$z), na.rm=TRUE ) * cell_area	
final <- exp(interp_logdens$z) * interp_h$z / sumdens	 
image(interp_h, col=colors2)

mtext(side=3, line=-2, seq[i])
}







###############################################################
###         visualizing nodes 						 	   ###
################################################################
pairs(exp(nodes4[[1]]))

png(file.path(figs, "4param_dist.png"), width=10, height=8, units="in", res=200)
col <- brewer.pal(3,"Set1")
par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(0.5,0.5,0.5,0.5))

for(i in 1:length(param4)){
	for(j in 1:length(param4)){
		if(param4[i]==param4[j]){
			plot(x=1, y=1, type="n", axes=F, ann=F, xlim=c(0,1), ylim=c(0,1))
			text(x=0.5,y=0.5, param4[i], cex=3)
			box()
		}
		if(param4[i]!=param4[j]){
			plot(x=1, y=1, type="n", 
				ylim=c(min(exp(nodes4[["Family"]][,param4[i]])), max(exp(nodes4[["Family"]][,param4[i]]))), 
				xlim=c(min(exp(nodes4[["Family"]][,param4[j]])), max(exp(nodes4[["Family"]][,param4[j]]))),
				cex.lab=2, cex.axis=1.5, axes=F, ann=F)
			box()
			for(z in c("Family","Genus","Species")){
				index <- which(labels==z)
				points(exp(nodes4[[index]][,param4[j]]), exp(nodes4[[index]][,param4[i]]), pch=19, col=paste0(col[index], "90"))
			}
			if(i==1) axis(3, cex.axis=1.3)
			if(j==length(param4)) axis(4, cex.axis=1.3)
			if(j==1) axis(2, cex.axis=1.3)
			if(i==length(param4)) axis(1, cex.axis=1.3)
		}
	}
}

# legend("topright", legend=c("Family", "Genus", "Species"), col=paste0(col,"90"), pch=19, cex=2)
dev.off()


png(file.path(figs, "4param_dist_species.png"), width=10, height=8, units="in", res=200)
col <- brewer.pal(3,"Set1")
par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(0.5,0.5,0.5,0.5))

for(i in 1:length(param4)){
	for(j in 1:length(param4)){
		if(param4[i]==param4[j]){
			plot(x=1, y=1, type="n", axes=F, ann=F, xlim=c(0,1), ylim=c(0,1))
			text(x=0.5,y=0.5, param4[i], cex=3)
			box()
		}
		if(param4[i]!=param4[j]){
			plot(x=1, y=1, type="n", 
				ylim=c(min(exp(nodes4[["Species"]][,param4[i]]))*0.95, max(exp(nodes4[["Species"]][,param4[i]]))*1.05), 
				xlim=c(min(exp(nodes4[["Species"]][,param4[j]]))*0.95, max(exp(nodes4[["Species"]][,param4[j]]))*1.05),
				cex.lab=2, cex.axis=1.5, axes=F, ann=F)
			box()
			for(z in "Species"){
				index <- which(labels==z)
				points(exp(nodes4[[index]][,param4[j]]), exp(nodes4[[index]][,param4[i]]), pch=19, col=paste0(col[index], "90"), cex=2)
			}
			if(i==1) axis(3, cex.axis=1.3)
			if(j==length(param4)) axis(4, cex.axis=1.3)
			if(j==1) axis(2, cex.axis=1.3)
			if(i==length(param4)) axis(1, cex.axis=1.3)
		}
	}
}

# legend("topright", legend=c("Family", "Genus", "Species"), col=paste0(col,"90"), pch=19, cex=2)
dev.off()


png(file.path(figs, "MK_dist.png"), width=10, height=8, units="in", res=200)
col <- brewer.pal(3,"Set1")
# local <- data.frame("M"=M_vec, "K"=vbk_vec)
par(mfrow=c(1,1), mar=c(4,5,2,2))
plot(x=1, y=1, type="n", 
	ylim=c(min(c(exp(nodes2[["Family"]][,"M"])), na.rm=TRUE), max(c(exp(nodes2[["Family"]][,"M"])), na.rm=TRUE)), 
	xlim=c(min(c(exp(nodes2[["Family"]][,"K"])),  na.rm=TRUE), max(c(exp(nodes2[["Family"]][,"K"])), na.rm=TRUE)), 
	cex.lab=2, cex.axis=1.5,
	xlab="von Bertalanffy k",
	ylab="Natural mortality")
for(i in c("Family", "Genus", "Species")){
	index <- which(labels == i)
	points(exp(nodes2[[index]][,"K"]), exp(nodes2[[index]][,"M"]), pch=19, cex=2, col=paste0(col[index],"90"))
}
# rug(side=2, local$M, lwd=2, col=gray(0.3))
# rug(side=1, local$K, lwd=2, col=gray(0.3))
# for(i in 1:nrow(local)){
# 	if(is.numeric(local$M[i]) & is.numeric(local$K[i])) points(local$K[i], local$M[i], pch=19, cex=2, col=paste0(gray(0.3),"90"))
# }
# legend("topleft", legend=c("Family", "Genus", "Species", "Local studies"), col=c(paste0(rev(col),"90"), paste0(gray(0.3),"90")), pch=19, cex=2)
legend("topleft", legend=c("Family", "Genus", "Species"), col=paste0(rev(col),"90"), pch=19, cex=2)
dev.off()

png(file.path(figs, "MK_dist_Species.png"), width=10, height=8, units="in", res=200)
col <- brewer.pal(3,"Set1")
# local <- data.frame("M"=M_vec, "K"=vbk_vec)
par(mfrow=c(1,1), mar=c(4,5,2,2))
plot(x=1, y=1, type="n", 
	ylim=c(min(c(exp(nodes2[["Species"]][,"M"])), na.rm=TRUE), max(c(exp(nodes2[["Species"]][,"M"])), na.rm=TRUE)), 
	xlim=c(min(c(exp(nodes2[["Species"]][,"K"])),  na.rm=TRUE), max(c(exp(nodes2[["Species"]][,"K"])), na.rm=TRUE)), 
	cex.lab=2, cex.axis=1.5,
	xlab="von Bertalanffy k",
	ylab="Natural mortality")
for(i in c("Species")){
	index <- which(labels == i)
	points(exp(nodes2[[index]][,"K"]), exp(nodes2[[index]][,"M"]), pch=19, cex=2, col=paste0(col[index],"90"))
}

dev.off()







