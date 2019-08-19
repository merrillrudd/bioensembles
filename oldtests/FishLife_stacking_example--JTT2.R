
###############################
### packages 				###
###############################
wd <- "C:\\merrill\\bioensembles"
# wd <- "C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2017 -- Rudd data-poor ensemble/2018-03-30/"
setwd(wd)

#devtools::install_github("merrillrudd/LIME", dependencies=TRUE, ref="multifleet")
library(LIME)

library(FishLife)

library(mvQuad)
library(mvtnorm)

####################################
###  call FishLife               ###
####################################
load(file.path(wd, "Return.Rdata"))

sp <- Search_species(Genus="Lachnolaimus", Species="maximus", ParentChild_gz=Return$ParentChild_gz)$match_taxonomy
orders <- c("Class","Order","Family","Genus","Species")
labels <- sapply(1:length(sp), function(x){
	split <- strsplit(sp[x],"_")[[1]]
	find <- which(split=="predictive")[1]
	if(is.na(find)) return(orders[length(orders)])
	if(is.na(find)==FALSE) return(orders[find-1])
})
Which <- sapply(1:length(sp), function(x) grep(sp[x], Return$ParentChild_gz[,"ChildName"]))
Mean <- lapply(1:length(sp), function(x) Return$beta_gv[Which[x],])
names(Mean) <- labels
Cov <- lapply(1:length(sp), function(x) Return$Cov_gvv[Which[x],,])
names(Cov) <- labels

lh_list <- create_lh_list(linf=exp(Mean[["Species"]]["Loo"]), vbk=exp(Mean[["Species"]]["K"]), t0=-1.77,
							lwa=0.0237, lwb=2.96,
							M=exp(Mean[["Species"]]["M"]),
							M50=exp(Mean[["Species"]]["tm"]), maturity_input="age",
							S50=exp(Mean[["Species"]]["tm"]), S95=exp(Mean[["Species"]]["tm"])*1.3, selex_input="age",
							SigmaF=0.2, SigmaR=0.737, rho=0.4,
							AgeMax=exp(Mean[["Species"]]["tmax"]),
							binwidth=1)

###########################
###  grid               ###
###########################

Quadrature = function( myGrid, f, ... ){
  Nodes = getNodes(myGrid)
  Weights = getWeights(myGrid)
  # Compute integral
  Integral = rep(NA, nrow(Nodes))
  for( i in seq_along(Integral) ){
    Integral[i] = f( Nodes[i,], ... ) * Weights[i]
  }
  return( sum(Integral) )
}

Integrate = function( Nodes_i, h_i, Mean, Cov, n_interp=100 ){
  #xy_orig = expand.grid( x=seq(-4,4,by=0.01), y=seq(-4,4,by=0.01) )
  #Chol = chol(Cov)
  #xy_loc = as.matrix(xy_orig) %*% Chol

  # Interpolate PDF and get integral
  logdens_i = dmvnorm( Nodes_i, mean=Mean, sigma=Cov, log=TRUE )
  interp_logdens = akima::interp( x=Nodes_i[,1], y=Nodes_i[,2], z=logdens_i, nx=n_interp, ny=n_interp )
  cell_area = mean(diff(interp_logdens$x)) * mean(diff(interp_logdens$y))
  sumdens = sum( exp(interp_logdens$z), na.rm=TRUE ) * cell_area

  # Interpolate product of function and PDF
  interp_h = akima::interp( x=Nodes_i[,1], y=Nodes_i[,2], z=h_i, nx=n_interp, ny=n_interp )
  interp_z = exp(interp_logdens$z) * interp_h$z / sumdens

  # Calculate interal and return
  Integral = sum( interp_z, na.rm=TRUE ) * cell_area
  return( Integral )
}


paramMK <- c("K","M")

gridMK <- createNIGrid(dim=2, type="GHe", level=4, ndConstruction="sparse")
msub <- Mean[["Species"]][which(names(Mean[["Species"]]) %in% paramMK)]
csub <- Cov[["Species"]][which(rownames(Cov[["Species"]]) %in% paramMK), which(colnames(Cov[["Species"]]) %in% paramMK)]
rescale(gridMK, m=msub, C=(csub+t(csub))/2, dec.type=1)

nodesMK <- getNodes(gridMK)
colnames(nodesMK) <- names(msub)
weightsMK <- getWeights(gridMK)

Quadrature(f=dmvnorm, mean=msub, sigma=csub, myGrid=gridMK)

#** outputs first node in nodesMK
f <- function(x, mean, sigma){
	Like <- dmvnorm(x, mean=mean, sigma=sigma)
	Val <- x[1]
	return(Like * Val)
}
Quadrature(f=f, mean=msub, sigma=csub, myGrid=gridMK)

################################
###  simulate data           ###
################################
seed <- 5
set.seed(seed)
vbk_choose <- rlnorm(1, mean=msub["K"], sd=sqrt(csub["K","K"]))
M_choose <- rlnorm(1, mean=msub["M"], sd=sqrt(csub["M","M"]))
plist <- with(lh_list, create_lh_list(linf=linf, vbk=vbk_choose, t0=t0,
									lwa=lwa, lwb=lwb,
									M=M_choose,
									M50=ML50, maturity_input="length",
									S50=SL50, S95=SL95, selex_input="length",
									SigmaF=0.2, SigmaR=0.737, rho=0.4,
									AgeMax=AgeMax,
									binwidth=binwidth,
									Fequil=1.1,
									theta=10,
									nfleets=nfleets))

data <- generate_data(modpath=wd, itervec=1,
					Fdynamics="Endogenous", Rdynamics="Constant",
					lh=plist,
					Nyears=20, Nyears_comp=20, comp_sample=200,
					init_depl=c(0.10,0.90),
					seed=rep(seed+1000,1))

iterpath <- file.path(wd, 1)
data <- readRDS(file.path(iterpath, "True.rds"))
input_data <- list("years"=data$years, "LF"=data$LF)

plot_LCfits(LFlist=list("LF"=data$LF[,,1]))

################################
###  run LIME at each node   ###
################################
# res <- lapply(1:nrow(nodesMK), function(x){

# 	vbk_inp <- exp(nodesMK[x,"K"])
# 	M_inp <- exp(nodesMK[x,"M"])
# 	lhinp <- with(lh_list,
# 				create_lh_list(linf=linf, vbk=vbk_inp, t0=t0,
# 								lwa=lwa, lwb=lwb,
# 								M=M_inp,
# 								M50=ML50, maturity_input="length",
# 								S50=SL50, S95=SL95, selex_input="length",
# 								SigmaF=SigmaF, SigmaR=SigmaR,
# 								AgeMax=AgeMax,
# 								binwidth=binwidth,
# 								theta=10,
# 								nfleets=nfleets))

# 	# input files and run model
# 	input <- create_inputs(lh=lhinp, input_data=input_data)
# 	out <- run_LIME(modpath=iterpath, input=input, data_avail="LC", newtonsteps=3)

# 		## check_convergence
# 		isNA <- all(is.null(out$df))
# 		if(isNA==TRUE){
# 			gradient <- FALSE
# 			pdHess <- FALSE
# 		}
# 		if(isNA==FALSE){
# 			gradient <- out$opt$max_gradient <= 0.001
# 			pdHess <- out$Sdreport$pdHess
# 		}

# 		## run set of adjustments to achieve convergence
# 		if(isNA == TRUE | gradient == FALSE | pdHess == FALSE){
# 			out <- get_converged(results=out)
# 		}

# 	return(out)
# })
# saveRDS(res, file.path(wd, "results.rds"))

# ## check that runs at all nodes converged
# for(i in 1:length(res)){
# 	if(is.null(out$df)) print(i)
# }

res <- readRDS(file.path(wd, "results.rds"))
##########################################################
###  predictive stacking (relative spawning biomass)   ###
##########################################################
library(dplyr)
library(RColorBrewer)

## summarize relative biomass information for all nodes
find <- lapply(1:length(res), function(x){
	if(is.null(res[[x]]$df)==FALSE){
		sub <- res[[x]]$Report
		mle <- sub$D_t
		sdrep <- summary(res[[x]]$Sdreport)
		sd <- sdrep[which(rownames(sdrep)=="lD_t"),2]
	}
	if(is.null(res[[x]]$df)){
		mle <- rep(NA, 20)
		sd <- rep(NA, 20)
	}
	df <- data.frame("Node"=x, "Year"=1:length(mle), "MLE"=mle, "SE"=sd)
	return(df)
})
find <- do.call(rbind, find)

## calculate annual prediction
yrs <- unique(find$Year)
stack1 <- sapply(yrs, function(x){
	yrsub <- find %>% filter(Year == yrs[x])
	out <- sum( yrsub[,"MLE"] * weightsMK ) 
	return(out)
})

stack2 <- sapply(yrs, function(x){
	yrsub <- find %>% filter(Year == yrs[x])
	out <- exp( sum( log(yrsub[,"MLE"]) * weightsMK ) ) 
	return(out)
})

## means
cols <- rev(brewer.pal(4, "Set1"))
plot(x=1,y=1,type="n",xlim=c(1,20), ylim=c(0,1), ylab="Relative spawning biomass", xlab="Year")
for(i in 1:length(res)){
	rep <- res[[i]]$Report
	lines(x=1:length(rep$D_t), y=rep$D_t, col="#AAAAAA80", lwd=3)
	lines(x=1:length(rep$D_t), y=data$D_t, col=cols[1], lwd=3)
	lines(x=1:length(stack1), y=stack1, col=cols[4], lwd=3)
	lines(x=1:length(stack2), y=stack2, col=cols[4], lwd=3, lty=2)
}

## density
plot(x=1, y=1, type="n", xaxs="i", yaxs="i", xlim=c(0,2), ylim=c(0,8.5), xlab="Relative spawning biomass in terminal year", ylab="Density", cex.lab=2, xaxt="n", las=2, cex.axis=1.5)

xx <- seq(0,2,by=0.001)
stack_dens <- sapply(1:length(res), function(x){
	sub <- find %>% filter(Node==x)
	term <- sub %>% filter(Year==max(Year))

## BEGIN CODE CHANGE
	#dxx <- dlnorm(xx, mean=log(term[,"MLE"]), sd=term[,"SE"])
	dxx <- dlnorm(xx, mean=log(term[,"MLE"])-term[,"SE"]^2/2, sd=term[,"SE"])
## END CODE CHANGE

	polygon(x=c(xx, rev(xx)), y=c(dxx, rep(0, length(dxx))), col="#AAAAAA80", border="#222222")
	return(dxx)
})

######## BEGIN NEW CODE
# Check that each integrates to 1.0
colSums( stack_dens * mean(diff(xx)) )
######## END NEW CODE


wdens <- sapply(1:length(res), function(x){
  return(weightsMK[x] * stack_dens[,x])
})
ysum <- sapply(1:nrow(wdens), function(x){
	sub <- wdens[x,]
	return(sum(sub))
})

############# BEGIN NEW CODE
unique_nodesMK = unique(nodesMK)
Match = match( unique_nodesMK[,1], nodesMK[,1] )
unique_stack_dens = stack_dens[,Match]
ysum = rep(NA,nrow(stack_dens))
for(i in 1:length(ysum)){
  ysum[i] = Integrate(h_i=unique_stack_dens[i,], Nodes_i=unique_nodesMK, Mean=msub, Cov=csub)
}
############## END NEW CODE

###### BEGIN CODE CHANGE
# Rescale to have integral of 1 for predictive distribution
#ysum_rescaled = ysum / sum(ysum) / mean(diff(xx))
ysum_rescaled = ysum
#polygon(x=c(xx, rev(xx)), y=c(ysum, rep(0, length(xx))), col=paste0(cols[4],"99"), border=cols[4])
polygon(x=c(xx, rev(xx)), y=c(ysum_rescaled, rep(0, length(xx))), col=paste0(cols[4],"99"), border=cols[4])
#abline(v=xx[which(ysum==max(ysum))], col=cols[4], lwd=4)
abline(v=weighted.mean(xx,w=ysum), col=cols[4], lwd=4)
###### END CODE CHANGE

abline(v=data$D_t[length(data$D_t)], col=cols[1], lwd=4)

