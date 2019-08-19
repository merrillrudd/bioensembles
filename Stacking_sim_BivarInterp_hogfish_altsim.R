rm(list=ls())

# devtools::install_github("james-thorson/FishLife")
library( FishLife )
library(mvQuad)

# devtools::install_github("merrillrudd/LIME")
library(LIME)
library(LBSPR)

library(mvtnorm)
library(Hmisc)
library(foreach)
library(doParallel)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(RuddR)

main_dir <- "C:\\merrill\\bioensembles"

sim_dir <- file.path(main_dir, "sim_hogfish_alt")
dir.create(sim_dir, showWarnings=FALSE)

res_dir <- file.path(sim_dir, "results_FishLifeDefault")
dir.create(res_dir, showWarnings=FALSE)

fig_dir <- file.path(sim_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

funs <- list.files(file.path(main_dir, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(main_dir,"R", funs[x])))

#########################
# Step 1:  Generate 4 parameters from 1/2 dimensional quadrature nodes
#########################

##################################
## Regional studies - life history
##################################
### Adams and Rios 2015
study_vec <- c(NA, "Puerto Rico", "Cuba", "GOM", "GOM", "GOM", "GOM", NA, "Dry Tortugas", "FL Keys", "S. Florida", "GOM", "Marquesas Key FL", "S. Florida", "GOM", "GOM", "Florida", "W. Florida", "SE US")
loo_vec <- c(566, 913, 850, 966, NA, 381, 896, 566, 651, 336, 428, 917, 397, 448, 939, NA, NA, 849, NA)
vbk_vec <- c(0.19, 0.08, 0.10, 0.08, NA, 0.56, 0.09, 0.19, 0.11, 0.56, 0.26, 0.08, 0.17, 0.23, 0.07, NA, NA, 0.11, NA)
t0_vec <- c(-0.78, -1.78, -1.38, -1.77, NA, -0.16, -1.98, -2.33, -3.34, -0.19, -0.93, -1.84, -3.74, -1.02, -2.01, NA, NA, -1.33, NA)
M_vec <- c(0.25, 0.13, rep(NA, 16), mean(c(0.16,0.29)))
Lm_vec <- c(NA, 249, NA, 152, 152, NA, NA, 198, NA, NA, NA, NA, NA, 163, 166, 177, NA, NA, mean(c(152, 193)))

lhstudies <- data.frame("Location"=study_vec, "Loo"=loo_vec/10, "K"=vbk_vec, "t0"=t0_vec, "M"=M_vec, "Lm"=Lm_vec/10)
ifb_rm <- which(lhstudies$Loo %in% c(91.7, 87.0, 65.1, 42.8, 39.7, 33.6))
fb_rm <- lhstudies[ifb_rm,]
lh_incl <- lhstudies[-ifb_rm,]

a_vec <- c(2.55e-5, 1.52e-2, 2.37e-2, 9.5e-5)
b_vec <- c(2.97, 3.11, 2.95, 2.75)

lhstudies_pr <- lhstudies[which(lhstudies$Location=="Puerto Rico"),]

lh_pr_equil <- create_lh_list(linf=lhstudies_pr[,"Loo"], vbk=lhstudies_pr[,"K"], t0=lhstudies_pr[,"t0"], M=lhstudies_pr[,"M"], M50=lhstudies_pr[,"Lm"], maturity_input="length", S50=15, S95=20, selex_input="length", selex_type="logistic", lwa=mean(a_vec), lwb=mean(b_vec), SigmaF=0.001, SigmaR=0.001, rho = 0, binwidth=1)
lh_pr_var <- create_lh_list(linf=lhstudies_pr[,"Loo"], vbk=lhstudies_pr[,"K"], t0=lhstudies_pr[,"t0"], M=lhstudies_pr[,"M"], M50=lhstudies_pr[,"Lm"], maturity_input="length", S50=15, S95=20, selex_input="length", selex_type="logistic", lwa=mean(a_vec), lwb=mean(b_vec), SigmaF=0.1, SigmaR=0.737, rho=0.4, binwidth=1)


##################################
## FishLife default hogfish
##################################
sp = Plot_taxa( Taxa=Search_species(Genus="Lachnolaimus",Species="maximus")$match_taxonomy, mfrow=c(2,2) )

Mean <- sp[[1]]$Mean_pred[c('M','K','Loo','Lm')]
Cov <- sp[[1]]$Cov_pred[c('M','K','Loo','Lm'),c('M','K','Loo','Lm')]

##################################
## updated FishLife quadrature nodes
##################################
Diag = function(vec){
  if( length(vec)==1 ) return(vec)
  if( length(vec)>=2 ) return(diag(vec))
}

## with observation error
Eigen <- eigen(Cov)

### ----- 1 Dimension ------ ###
Dim = 1

# Bivariate quadrature
myGrid <- createNIGrid(dim=Dim, type="GHe", level=4,  ndConstruction="sparse")
Nodes1_i = getNodes(myGrid)

# Use eigen-decomposition to project quadrature nodes into all variables
Param1_i = t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Nodes1_i) )
colnames(Param1_i) = names(Mean)

### ----- 2 Dimensions ------ ###
Dim = 2

# Bivariate quadrature
myGrid <- createNIGrid(dim=Dim, type="GHe", level=4,  ndConstruction="sparse")
Nodes2_i = getNodes(myGrid)

# Use eigen-decomposition to project quadrature nodes into all variables
Param2_i = t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Nodes2_i) )
colnames(Param2_i) = names(Mean)


param <- names(Mean)
png(file.path(fig_dir, "4param_1D.png"), height=10, width=11, units="in", res=200)
par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(0.5,0.5,0.5,0.5))

for(i in 1:length(param)){
  for(j in 1:length(param)){
    if(param[i]==param[j]){
      plot(x=1, y=1, type="n", axes=F, ann=F, xlim=c(0,1), ylim=c(0,1))
      if(param[i]=="Loo") lplot <- "Linf"
      if(param[i]!="Loo") lplot <- param[i]
      text(x=0.5,y=0.5, lplot, cex=3)
      box()
    }
    if(param[i]!=param[j]){
      plot(x=1, y=1, type="n", 
        ylim=c(min(c(exp(Param2_i[,param[i]])),na.rm=TRUE), max(c(exp(Param2_i[,param[i]])),na.rm=TRUE)), 
        xlim=c(min(c(exp(Param2_i[,param[j]])),na.rm=TRUE), max(c(exp(Param2_i[,param[j]])),na.rm=TRUE)),
        cex.lab=2, cex.axis=2, axes=F, ann=F)
      box()
      # points(exp(test[,param[j]]), exp(test[,param[i]]), pch=19, col="#44444450", cex=1.5)
      # points(exp(Param2_i[,param[j]]), exp(Param2_i[,param[i]]), pch=19, col="#AA000050", cex=1.5)
      points(exp(Param1_i[,param[j]]), exp(Param1_i[,param[i]]), pch=19, col="#0000AA50", cex=4)

      if(i==1) axis(3, cex.axis=1.8)
      if(j==length(param)) axis(4, cex.axis=1.8)
      if(j==1) axis(2, cex.axis=1.8)
      if(i==length(param)) axis(1, cex.axis=1.8)
    }
  }
}
dev.off()

png(file.path(fig_dir, "4param_2D.png"), height=10, width=11, units="in", res=200)
par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(0.5,0.5,0.5,0.5))

for(i in 1:length(param)){
  for(j in 1:length(param)){
    if(param[i]==param[j]){
      plot(x=1, y=1, type="n", axes=F, ann=F, xlim=c(0,1), ylim=c(0,1))
      if(param[i]=="Loo") lplot <- "Linf"
      if(param[i]!="Loo") lplot <- param[i]
      text(x=0.5,y=0.5, lplot, cex=3)
      box()
    }
    if(param[i]!=param[j]){
      plot(x=1, y=1, type="n", 
        ylim=c(min(c(exp(Param2_i[,param[i]])),na.rm=TRUE), max(c(exp(Param2_i[,param[i]])),na.rm=TRUE)), 
        xlim=c(min(c(exp(Param2_i[,param[j]])),na.rm=TRUE), max(c(exp(Param2_i[,param[j]])),na.rm=TRUE)),
        cex.lab=2, cex.axis=2, axes=F, ann=F)
      box()
      # points(exp(test[,param[j]]), exp(test[,param[i]]), pch=19, col="#44444450", cex=1.5)
      points(exp(Param2_i[,param[j]]), exp(Param2_i[,param[i]]), pch=19, col="#AA000050", cex=2)
      points(exp(Param1_i[,param[j]]), exp(Param1_i[,param[i]]), pch=19, col="#0000AA50", cex=4)

      if(i==1) axis(3, cex.axis=1.8)
      if(j==length(param)) axis(4, cex.axis=1.8)
      if(j==1) axis(2, cex.axis=1.8)
      if(i==length(param)) axis(1, cex.axis=1.8)
    }
  }
}
dev.off()

png(file.path(fig_dir, "4param_FishLifeDefault.png"), height=10, width=11, units="in", res=200)
param <- c("M","K","Loo","Lm")
par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(0.5,0.5,0.5,0.5))

for(i in 1:length(param)){
  for(j in 1:length(param)){
    if(param[i]==param[j]){
      plot(x=1, y=1, type="n", axes=F, ann=F, xlim=c(0,1), ylim=c(0,1))
      if(param[i]=="Loo") lplot <- "Linf"
      if(param[i]!="Loo") lplot <- param[i]
      text(x=0.5,y=0.5, lplot, cex=3)
      box()
    }
    if(param[i]!=param[j]){
      if(param[i]=="M") yplot <- M_vec
      if(param[i]=="K") yplot <- vbk_vec
      if(param[i]=="Loo") yplot <- loo_vec/10
      if(param[i]=="Lm") yplot <- Lm_vec/10
      if(param[j]=="M") xplot <- M_vec
      if(param[j]=="K") xplot <- vbk_vec
      if(param[j]=="Loo") xplot <- loo_vec/10
      if(param[j]=="Lm") xplot <- Lm_vec/10

      plot(x=1, y=1, type="n", 
        ylim=c(min(c(exp(Param2_i[,param[i]])),na.rm=TRUE), max(c(exp(Param2_i[,param[i]])),na.rm=TRUE)), 
        xlim=c(min(c(exp(Param2_i[,param[j]])),na.rm=TRUE), max(c(exp(Param2_i[,param[j]])),na.rm=TRUE)),
        cex.lab=2, cex.axis=1.8, axes=F, ann=F)

      box()
      names(xplot) <-study_vec
      names(yplot) <- study_vec
      lwd <- rep(3, length(xplot))
      lwd[which(study_vec=="Puerto Rico")] <- 5
      col <- rep("#AAAAAA", length(xplot))
      col[which(study_vec=="Puerto Rico")] <- "black"
      order <- order(lwd)
      col_new <- col[order]
      lwd_new <- lwd[order]
      xplot_new <- xplot[order]
      yplot_new <- yplot[order]
      # points(xplot_new, yplot_new, pch=4, col=col_new, cex=2, lwd=lwd_new, xpd=NA)
      points(exp(Param2_i[,param[j]]), exp(Param2_i[,param[i]]), pch=19, col="#FF000050", cex=1.5)
      points(exp(Param1_i[,param[j]]), exp(Param1_i[,param[i]]), pch=22, col="#FF0000", bg="#FF000050", cex=3)

      if(i==1) axis(3, cex.axis=1.8)
      if(j==length(param)) axis(4, cex.axis=1.8)
      if(j==1) axis(2, cex.axis=1.8)
      if(i==length(param)) axis(1, cex.axis=1.8)
    }
  }
}
dev.off()





## --------------- begin simulation section ------------------##
itervec <- 1:100

##########################
## Simulate equilibrium
##########################
equil_dir <- file.path(res_dir, "equil")
dir.create(equil_dir, showWarnings=FALSE)

start <- Sys.time()
for(loop in itervec){
  generate_data(modpath=equil_dir,
                itervec=loop,
                Fdynamics="Constant",
                Rdynamics="Constant", 
                lh=lh_pr_equil,
                Nyears=20,
                Nyears_comp=20,
                comp_sample=200,
                init_depl=c(0.1,0.9),
                seed=rep(loop,loop))
}
end <- Sys.time() - start

##########################
## Simulate variability
##########################
var_dir <- file.path(res_dir, "variable")
dir.create(var_dir, showWarnings=FALSE)

start <- Sys.time()
for(loop in itervec){
  generate_data(modpath=var_dir,
                itervec=loop,
                Fdynamics="Endogenous",
                Rdynamics="Constant", 
                lh=lh_pr_var,
                Nyears=20,
                Nyears_comp=20,
                comp_sample=200,
                init_depl=c(0.1,0.9),
                seed=rep(loop,loop))
}
end <- Sys.time() - start

png(file.path(fig_dir, "Simulation_examples.png"), height=8, width=7, units="in", res=200)
col <- brewer.pal(3, "Set1")
par(mfcol=c(3,2), mar=c(0,0,0,0), omi=c(1,0.6,0.4,0.4))
set.seed(123)
choose <- sample(itervec, 3)
plot(x=1,y=1,type="n", cex.axis=1.2, ylim=c(0,1), xlim=c(1,20), xaxt="n")
for(i in 1:length(itervec)){
  true <- readRDS(file.path(equil_dir, itervec[i], "True.rds"))
  lines(true$SPR_t, type="l", col="#AAAAAA80")
}
for(i in 1:length(choose)){
  true <- readRDS(file.path(equil_dir, choose[i], "True.rds"))
  lines(true$SPR_t, type="l", lwd=3, col=col[i])  
}
mtext(side=2, "SPR", cex=1.3, line=2)
mtext(side=3, "Equilibrium", cex=1.3, line=1)
plot(x=1,y=1,type="n", cex.axis=1.2, ylim=c(0,0.5), xlim=c(1,20), xaxt="n")
for(i in 1:length(itervec)){
  true <- readRDS(file.path(equil_dir, itervec[i], "True.rds"))
  lines(true$F_ft[1,], type="l", col="#AAAAAA80")
}
for(i in 1:length(choose)){
  true <- readRDS(file.path(equil_dir, choose[i], "True.rds"))
  lines(true$F_ft[1,], type="l", lwd=3, col=col[i])  
}
mtext(side=2, "Fishing mortality", cex=1.3, line=2)
plot(x=1,y=1,type="n", cex.axis=1.2, ylim=c(0,4), xlim=c(1,20))
for(i in 1:length(itervec)){
  true <- readRDS(file.path(equil_dir, itervec[i], "True.rds"))
  lines(true$R_t, type="l", col="#AAAAAA80")
}
for(i in 1:length(choose)){
  true <- readRDS(file.path(equil_dir, choose[i], "True.rds"))
  lines(true$R_t, type="l", lwd=3, col=col[i])  
}
mtext(side=1, "Year", cex=1.3, line=2.5)
mtext(side=2, "Recruitment", cex=1.3, line=2)
plot(x=1,y=1,type="n", cex.axis=1.2, ylim=c(0,1), xlim=c(1,20), xaxt="n", yaxt="n")
for(i in 1:length(itervec)){
  true <- readRDS(file.path(var_dir, itervec[i], "True.rds"))
  lines(true$SPR_t, type="l", col="#AAAAAA80")
}
for(i in 1:length(choose)){
  true <- readRDS(file.path(var_dir, choose[i], "True.rds"))
  lines(true$SPR_t, type="l", lwd=3, col=col[i])  
}
mtext(side=3, "Variable", cex=1.3, line=1)
plot(x=1,y=1,type="n", cex.axis=1.2, ylim=c(0,0.5), xlim=c(1,20), xaxt="n", yaxt="n")
for(i in 1:length(itervec)){
  true <- readRDS(file.path(var_dir, itervec[i], "True.rds"))
  lines(true$F_ft[1,], type="l", col="#AAAAAA80")
}
for(i in 1:length(choose)){
  true <- readRDS(file.path(var_dir, choose[i], "True.rds"))
  lines(true$F_ft[1,], type="l", lwd=3, col=col[i])  
}
plot(x=1,y=1,type="n", cex.axis=1.2, ylim=c(0,4), xlim=c(1,20), yaxt="n")
for(i in 1:length(itervec)){
  true <- readRDS(file.path(var_dir, itervec[i], "True.rds"))
  lines(true$R_t, type="l", col="#AAAAAA80")
}
for(i in 1:length(choose)){
  true <- readRDS(file.path(var_dir, choose[i], "True.rds"))
  lines(true$R_t, type="l", lwd=3, col=col[i])  
}
mtext(side=1, "Year", cex=1.3, line=2.5)
dev.off()

##########################
## Run at truth
##########################

## ---- equilibrium, LBSPR
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_single(modpath=equil_dir, iter=loop, lh=lh_pr_equil, model="LBSPR", name="BestCase")
end <- Sys.time() - start
stopCluster(cl)

re <- cover <- rep(NA, length(itervec))
for(i in 1:length(itervec)){
  true <- readRDS(file.path(equil_dir, itervec[i], "True.rds"))
  tval <- true$SPR[length(true$SPR)]
  res <- readRDS(file.path(equil_dir, itervec[i], "res_BestCase_LBSPR.rds"))
  eval <- res@SPR[length(res@SPR)]
  esd <- sqrt(res@Vars[,"SPR"][length(res@Vars[,"SPR"])])

  re[i] <- (eval - tval)/tval
  cover[i] <- ifelse(tval >= (eval - 0.674*esd) & tval <= (eval + 0.674*esd), 1, 0)
}

## ---- equilibrium, LIME
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_single(modpath=equil_dir, iter=loop, lh=lh_pr_equil, model="LIME", name="BestCase")
end <- Sys.time() - start
stopCluster(cl)

re <- cover <- rep(NA, length(itervec))
for(i in 1:length(itervec)){
  true <- readRDS(file.path(equil_dir, itervec[i], "True.rds"))
  tval <- true$SPR[length(true$SPR)]
  if(file.exists(file.path(equil_dir, itervec[i], "res_BestCase_LIME.rds"))){
    res <- readRDS(file.path(equil_dir, itervec[i], "res_BestCase_LIME.rds"))
    eval <- res$Report$SPR_t[length(res$Report$SPR_t)]
    sd <- summary(res$Sdreport)
    esd <- sd[which(rownames(sd) == "SPR_t")[length(which(rownames(sd)=="SPR_t"))],2] 

    re[i] <- (eval - tval)/tval
    cover[i] <- ifelse(tval >= (eval - 0.674*esd) & tval <= (eval + 0.674*esd), 1, 0)    
  }
}

## ---- variability, LBSPR
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_single(modpath=var_dir, iter=loop, lh=lh_pr_var, model="LBSPR", name="BestCase")
end <- Sys.time() - start
stopCluster(cl)

re <- cover <- rep(NA, length(itervec))
for(i in 1:length(itervec)){
  true <- readRDS(file.path(var_dir, itervec[i], "True.rds"))
  tval <- true$SPR[length(true$SPR)]
  res <- readRDS(file.path(var_dir, itervec[i], "res_BestCase_LBSPR.rds"))
  eval <- res@SPR[length(res@SPR)]
  esd <- sqrt(res@Vars[,"SPR"][length(res@Vars[,"SPR"])])

  re[i] <- (eval - tval)/tval
  cover[i] <- ifelse(tval >= (eval - 0.674*esd) & tval <= (eval + 0.674*esd), 1, 0)
}


## ---- variability, LIME
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_single(modpath=var_dir, iter=loop, lh=lh_pr_var, model="LIME", name="BestCase")
end <- Sys.time() - start
stopCluster(cl)

re <- cover <- rep(NA, length(itervec))
for(i in 1:length(itervec)){
  true <- readRDS(file.path(var_dir, itervec[i], "True.rds"))
  tval <- true$SPR[length(true$SPR)]
  if(file.exists(file.path(var_dir, itervec[i], "res_BestCase_LIME.rds"))){
    res <- readRDS(file.path(var_dir, itervec[i], "res_BestCase_LIME.rds"))
    eval <- res$Report$SPR_t[length(res$Report$SPR_t)]
    sd <- summary(res$Sdreport)
    esd <- sd[which(rownames(sd) == "SPR_t")[length(which(rownames(sd)=="SPR_t"))],2] 

    re[i] <- (eval - tval)/tval
    cover[i] <- ifelse(tval >= (eval - 0.674*esd) & tval <= (eval + 0.674*esd), 1, 0)    
  }
}

##########################
## Run at FishLife means
##########################

lh_means <-   with(lh_pr_var, create_lh_list(vbk=exp(Mean["K"]), linf=exp(Mean["Loo"]),
                 M=exp(Mean["M"]),
                 M50=exp(Mean["Lm"]), maturity_input="length",
                 lwa=lwa, lwb=lwa,
                 S50=exp(Mean["Lm"]), S95=min(exp(Mean["Loo"])*0.95, exp(Mean["Lm"])*1.2), selex_input="length",
                 SigmaF=0.1, SigmaR=0.737, rho=0, 
                 binwidth=1, theta=10))

## ---- equilibrium, LBSPR
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_single(modpath=equil_dir, iter=loop, lh=lh_means, model="LBSPR", name="Means")
end <- Sys.time() - start
stopCluster(cl)

## ---- equilibrium, LIME
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_single(modpath=equil_dir, iter=loop, lh=lh_means, model="LIME", name="Means")
end <- Sys.time() - start
stopCluster(cl)

## ---- variability, LBSPR
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_single(modpath=var_dir, iter=loop, lh=lh_means, model="LBSPR", name="Means")
end <- Sys.time() - start
stopCluster(cl)

## ---- variability, LIME
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_single(modpath=var_dir, iter=loop, lh=lh_means, model="LIME", name="Means")
end <- Sys.time() - start
stopCluster(cl)

##########################
## Run Stacking - 1D
##########################

start1 <- Sys.time()
run_stacking(modpath=equil_dir, iter=1, lh=lh_means, model="LBSPR", nodes=Param1_i, dim="1D")
end1 <- Sys.time() - start1

start2 <- Sys.time()
run_stacking(modpath=equil_dir, iter=1, lh=lh_means, model="LBSPR", nodes=Param2_i, dim="2D")
end2 <- Sys.time() - start2

start1_2 <- Sys.time()
run_stacking(modpath=equil_dir, iter=1, lh=lh_means, model="LIME", nodes=Param1_i, dim="1D")
end1_2 <- Sys.time() - start1_2

start2_2 <- Sys.time()
run_stacking(modpath=equil_dir, iter=1, lh=lh_means, model="LIME", nodes=Param2_i, dim="2D")
end2_2 <- Sys.time() - start2_2

## ---- equilibrium, LBSPR
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_stacking(modpath=equil_dir, iter=loop, lh=lh_means, model="LBSPR", nodes=Param1_i, dim="1D")
end <- Sys.time() - start
stopCluster(cl)

## ---- equilibrium, LIME
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_stacking(modpath=equil_dir, iter=loop, lh=lh_means, model="LIME", nodes=Param1_i, dim="1D")
end <- Sys.time() - start
stopCluster(cl)

## ---- variability, LBSPR
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_stacking(modpath=var_dir, iter=loop, lh=lh_means, model="LBSPR", nodes=Param1_i, dim="1D")
end <- Sys.time() - start
stopCluster(cl)

## ---- variability, LIME
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_stacking(modpath=var_dir, iter=loop, lh=lh_means, model="LIME", nodes=Param1_i, dim="1D")
end <- Sys.time() - start
stopCluster(cl)

##########################
## Run Stacking - 2D
##########################

## ---- equilibrium, LBSPR
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_stacking(modpath=equil_dir, iter=loop, lh=lh_means, model="LBSPR", nodes=Param2_i, dim="2D")
end <- Sys.time() - start
stopCluster(cl)

## ---- equilibrium, LIME
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_stacking(modpath=equil_dir, iter=loop, lh=lh_means, model="LIME", nodes=Param2_i, dim="2D")
end <- Sys.time() - start
stopCluster(cl)

## ---- variability, LBSPR
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_stacking(modpath=var_dir, iter=loop, lh=lh_means, model="LBSPR", nodes=Param2_i, dim="2D")
end <- Sys.time() - start
stopCluster(cl)

## ---- variability, LIME
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)      
start <- Sys.time()
foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR')) %dopar%
  run_stacking(modpath=var_dir, iter=loop, lh=lh_means, model="LIME", nodes=Param2_i, dim="2D")
end <- Sys.time() - start
stopCluster(cl)

##########################
## Summarize results
##########################

models <- c("LBSPR", "LIME")
scenario <- c("Equilibrium", "Variable")

## non-stacking results
approach <- c("BestCase","Means")
res <- lapply(1:length(scenario), function(ss){
  byApproach <- lapply(1:length(approach), function(aa){
    byIter <- lapply(itervec, function(ii){
      if(scenario[ss]=="Equilibrium") true <- readRDS(file.path(equil_dir, itervec[ii], "True.rds"))
      if(scenario[ss]=="Variable") true <- readRDS(file.path(var_dir, itervec[ii], "True.rds"))
      byMod <- lapply(1:length(models), function(mm){
        if(scenario[ss]=="Equilibrium"){
          if(file.exists(file.path(equil_dir, itervec[ii], paste0("res_", approach[aa], "_", models[mm], ".rds")))==TRUE){
            res <- readRDS(file.path(equil_dir, itervec[ii], paste0("res_", approach[aa], "_", models[mm], ".rds")))
          } else { res <- NA }
        }
        if(scenario[ss]=="Variable"){
          if(file.exists(file.path(var_dir, itervec[ii], paste0("res_", approach[aa], "_", models[mm], ".rds")))==TRUE){
            res <- readRDS(file.path(var_dir, itervec[ii], paste0("res_", approach[aa], "_", models[mm], ".rds")))
          } else { res <- NA }
        }
        if(all(is.na(res))==FALSE){
          if(models[mm]=="LBSPR"){
            mle <- res@SPR
            sd <- sqrt(res@Vars[,"SPR"])
            years <- res@Years
          }
          if(models[mm]=="LIME"){
            mle <- res$Report$SPR_t
            sd1 <- summary(res$Sdreport)
            sd <- sd1[which(rownames(sd1)=="SPR_t"),2]
            years <- res$input$years
          }         
        } else {
          years <- true$years
          mle <- NA
          sd <- NA
        }
          df <- data.frame("Iteration"=itervec[ii], "Population"=scenario[ss], "Approach"=approach[aa], "Model"=models[mm], "Year"=years, "Variable"="SPR", "MLE"=mle, "SD"=sd, "LCL"=mle - 0.674*sd, "UCL"=mle + 0.674*sd, "True"=true$SPR_t)
          return(df)
      })
      byMod <- do.call(rbind, byMod)
      return(byMod)
    })
    byIter <- do.call(rbind, byIter)
    return(byIter)
  })
  byApproach <- do.call(rbind, byApproach)
  return(byApproach)
})
res <- do.call(rbind, res)

## model averaging
bivar_approach <- c("1D", "2D")
models <- c("LBSPR","LIME")
scenario <- c("Equilibrium","Variable")
res_avg <- lapply(1:length(scenario), function(ss){
# for(ss in 1:length(scenario)){
  byApproach <- lapply(1:length(bivar_approach), function(aa){
  # for(aa in 1:length(bivar_approach)){
    byIter <- lapply(1:length(itervec), function(ii){
    # for(ii in 1:length(itervec)){
      if(scenario[ss]=="Equilibrium") iterpath <- file.path(equil_dir, itervec[ii])
      if(scenario[ss]=="Variable") iterpath <- file.path(var_dir, itervec[ii])
      true <- readRDS(file.path(iterpath, "True.rds"))
      tval <- true$SPR_t

      # for(mm in 1:length(models)){
      byMod <- lapply(1:length(models), function(mm){
        ## read files
        if(file.exists(file.path(iterpath, paste0("res_stacking_", models[mm], "_", bivar_approach[aa], ".rds")))){
          stack <- readRDS(file.path(iterpath, paste0("res_stacking_", models[mm], "_", bivar_approach[aa], ".rds")))
          if(models[mm]=="LBSPR"){
            years <- stack[[1]]@Years
              y_i <- sapply(1:length(stack), function(x) stack[[x]]@SPR)
              ysd_i <- sapply(1:length(stack), function(x) sqrt(stack[[x]]@Vars[,"SPR"]))
            NLL_i <- sapply(1:length(stack), function(x) stack[[x]]@NLL)
          }
          if(models[mm]=="LIME"){
            years <- stack[[1]]$input$years
              y_i <- sapply(1:length(stack), function(x) stack[[x]]$Report$SPR_t)
              ysd_i <- sapply(1:length(stack), function(x){
                sdrep <- summary(stack[[x]]$Sdreport)
                sd <- sdrep[which(rownames(sdrep)=="SPR_t"),2]
                return(sd)
              })
          }
        } else {
          stack <- NULL
          y_i <- NULL
          ysd_i <- NULL
          NLL_i <- NULL
          years <- 1:20
        }

        if(all(is.null(y_i))==FALSE){
          avgmean <- sapply(1:nrow(y_i), function(x) mean(y_i[x,], na.rm=TRUE))
        } else {
          avgmean <- rep(NA, length(years))
        }
        if(all(is.null(y_i))==FALSE){
          avgsd <- sapply(1:nrow(ysd_i), function(x) mean(ysd_i[x,],na.rm=TRUE))
          lcl <- avgmean - 0.675*avgsd
          ucl <- avgmean + 0.675*avgsd
        } else { 
          avgsd <- rep(NA, length(years))
          lcl <- rep(NA, length(years))
          ucl <- rep(NA, length(years))
        }
        df <- data.frame("Iteration"=itervec[ii], "Population"=scenario[ss], "Approach"=paste0("Average",bivar_approach[aa]), "Model"=models[mm], "Year"=years, "MLE"=avgmean, "SD"=avgsd, "LCL"=lcl, "UCL"=ucl, "Variable"="SPR", "True"=tval)
        return(df)
      })
      byMod <- do.call(rbind, byMod)
      return(byMod)
    })
    byIter <- do.call(rbind, byIter)
    return(byIter)
  })
  byApproach <- do.call(rbind, byApproach)
  return(byApproach)
})
res_avg <- do.call(rbind, res_avg)

## stacking results
bivar_approach <- c("1D")#, "2D")
models <- c("LBSPR","LIME")
scenario <- c("Equilibrium","Variable")
res_stack <- lapply(1:length(scenario), function(ss){
# for(ss in 1:length(scenario)){
  byApproach <- lapply(1:length(bivar_approach), function(aa){
  # for(aa in 1:length(bivar_approach)){
  byIter <- lapply(1:length(itervec), function(ii){
    # for(ii in 1:length(itervec)){
    if(scenario[ss]=="Equilibrium") iterpath <- file.path(equil_dir, itervec[ii])
    if(scenario[ss]=="Variable") iterpath <- file.path(var_dir, itervec[ii])
    true <- readRDS(file.path(iterpath, "True.rds"))
    byMod <- lapply(1:length(models), function(mm){
    # for(mm in 1:length(models)){
    if(models[mm]=="LBSPR") vals <- "SPR"
    if(models[mm]=="LIME") vals <- c("SPR")
    byVal <- lapply(1:length(vals), function(v){
      # for(v in 1:length(vals)){
      if(vals[v]=="SPR"){
        y_seq <- seq(0.0001,1,by=0.0001)
        tval <- true$SPR_t
      }
      if(vals[v]=="BB0"){
        y_seq <- seq(0,4,by=0.001)
        tval <- true$D_t
      }
        ## read files
        if(file.exists(file.path(iterpath, paste0("res_stacking_", models[mm], "_", bivar_approach[aa], ".rds")))){
          stack <- readRDS(file.path(iterpath, paste0("res_stacking_", models[mm], "_", bivar_approach[aa], ".rds")))
          if(models[mm]=="LBSPR"){
            years <- stack[[1]]@Years
            if(vals[v]=="SPR"){
              y_i <- sapply(1:length(stack), function(x) stack[[x]]@SPR)
              ysd_i <- sapply(1:length(stack), function(x) sqrt(stack[[x]]@Vars[,"SPR"]))
            }
            NLL_i <- sapply(1:length(stack), function(x) stack[[x]]@NLL)
          }
          if(models[mm]=="LIME"){
            years <- stack[[1]]$input$years
            if(vals[v]=="SPR"){
              y_i <- sapply(1:length(stack), function(x) stack[[x]]$Report$SPR_t)
              ysd_i <- sapply(1:length(stack), function(x){
                sdrep <- summary(stack[[x]]$Sdreport)
                sd <- sdrep[which(rownames(sdrep)=="SPR_t"),2]
                return(sd)
              })
            }
            if(vals[v]=="BB0"){
              y_i <- sapply(1:length(stack), function(x) stack[[x]]$Report$D_t)
              ysd_i <- sapply(1:length(stack), function(x){
                sdrep <- summary(stack[[x]]$Sdreport)
                sd <- sdrep[which(rownames(sdrep)=="lD_t"),2]
              })
              NLL_i <- sapply(1:length(stack), function(x) stack[[x]]$Report$jnll)
            }
          }
        } else {
          stack <- NULL
          y_i <- NULL
          ysd_i <- NULL
          NLL_i <- NULL
          years <- 1:20
        }
        ## smooth SPR_i in 1/2D space
        if(all(is.null(y_i))==FALSE){
          if(bivar_approach[aa]=="1D"){
            Dim <- 1
            # Step 3A: Use approx() to do 1D smoother of SPR_i across nodes
            interp_h <- lapply(1:length(years), function(y) approx(x=Nodes1_i[,1], y=y_i[y,], n=100))
            y_j <- t(sapply(1:length(years), function(y) interp_h[[y]]$y))

            interp_k <- lapply(1:length(years), function(y){
            # for(y in 1:length(years)){
              bad <- which(is.na(ysd_i[y,]))
              if(length(bad)==0){
                nodes_inp <- Nodes1_i[,1]
                ysd_inp <- ysd_i[y,]
              }
              if(length(bad)>0){
                nodes_inp <- Nodes1_i[-bad,1]
                ysd_inp <- ysd_i[y,-bad]
              }
              if(length(ysd_inp)>1) return(approx(x=nodes_inp, y=ysd_inp, n=100))
              if(length(ysd_inp)<=1) return(list(x=rep(NA,100), y=rep(NA,100)))
            })
            ysd_j <- t(sapply(1:length(years), function(y) interp_k[[y]]$y))

            # Step 3B: translate Grid to 4 parameters and calculate FishLife probability for each
            Grid_j <- interp_h[[1]]$x
            Param_j <- t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j) )
            colnames(Param_j) <- names(Mean)
            Prob_j <- dmvnorm( Param_j, mean=Mean, sigma=Cov ) 

            # Step 3C: Calculate mean
            wtmean <- sapply(1:length(years), function(y) weighted.mean( y_j[y,], w=Prob_j, na.rm=TRUE)) 

            var_info <- t(sapply(1:length(years), function(y){
              if(all(is.na(ysd_j[y,]))==FALSE){ # wtmean[y] <= max(y_seq) & 
                ## density at each interpolated point - interpolated mean and sd from assessment
                dstack <- sapply(1:length(y_j[y,]), function(i){
                  d <- dnorm(y_seq, mean=y_j[y,i], sd=ysd_j[y,i])
                  return(d)
                })
                ## density at interpolated point weighted by FishLife probability
                wdstack <- dstack * outer( rep(1,nrow(dstack)), Prob_j )
                # summary of weighted distribution
                dnew <- rowSums(wdstack)
                quants <- wtd.quantile(y_seq, w=dnew)       

                ## weighted quantile of resulting distribution
                lcl <-  wtd.quantile(y_seq, w=dnew)[2]
                ucl <- wtd.quantile(y_seq, w=dnew)[4]
                wtsd <- sqrt(wtd.var(y_seq, w=dnew))              
              } else{
                wtsd <- NA
                lcl <- NA
                ucl <- NA
              }        
              return(data.frame("SD"=wtsd, "LCL"=lcl, "UCL"=ucl))      
            }))

          }
          if(bivar_approach[aa]=="2D"){  
            Dim <- 2    

            # Step 3A: Use akima::interp to do 2D smoother of SPR_i across nodes
            interp_h <- lapply(1:length(years), function(y) akima::interp( x=Nodes2_i[,1], y=Nodes2_i[,2], z=y_i[y,], nx=100, ny=100, duplicate="median"))
            y_j <- lapply(1:length(years), function(y) interp_h[[y]]$z) 

            interp_k <- lapply(1:length(years), function(y){
            # for(y in 1:length(years)){
              bad <- which(is.na(ysd_i[y,]))
              if(length(bad)==0){
                nodes_inpx <- Nodes2_i[,1]
                nodes_inpy <- Nodes2_i[,2]
                ysd_inp <- ysd_i[y,]
              }
              if(length(bad)>0){
                nodes_inpx <- Nodes2_i[-bad,1]
                nodes_inpy <- Nodes2_i[-bad,2]
                ysd_inp <- ysd_i[y,-bad]
              }
              if(length(ysd_inp)>1) interp <- akima::interp( x=nodes_inpx, y=nodes_inpy, z=ysd_inp, nx=100, ny=100, duplicate="median")
              if(length(ysd_inp)<=1) interp <- list(x=rep(NA,100),y=rep(NA,100), z=matrix(NA, nrow=100,ncol=100))
              return(interp)
            })
            ysd_j <- lapply(1:length(years), function(y) interp_k[[y]]$z)     

            # Step 3B: translate Grid to 4 parameters and calculate FishLife probability for each
           notNA <- sapply(1:length(interp_h), function(y) all(is.na(interp_h[[y]]$x))==FALSE)
           index <- which(notNA)[1]
           Grid_j <- expand.grid(interp_h[[index]]$x, interp_h[[index]]$y)
           Param_j <- t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j) )
           colnames(Param_j) <- names(Mean)
           Prob_j <- matrix(dmvnorm( Param_j, mean=Mean, sigma=Cov ), nrow=100)           
        
           wtmean <- sapply(1:length(y_j), function(y) weighted.mean(y_j[[y]], w=matrix(Prob_j, nrow=100), na.rm=TRUE))

           var_info <- t(sapply(1:length(years), function(y){
               if(wtmean[y] <= max(y_seq) & all(is.na(ysd_j[[y]]))==FALSE){
                 ## density at each interpolated point - interpolated mean and sd from assessment
                 dstack <- sapply(1:nrow(y_j[[y]]), function(i){
                  sub <- sapply(1:ncol(y_j[[y]]), function(j){
                    d <- dnorm(y_seq, mean=y_j[[y]][i,j], sd=ysd_j[[y]][i,j])
                    return(d)
                  })
                  return(sub)
                 }, simplify="array")
                 dprob <- dstack * outer( rep(1,length(y_seq)), Prob_j )
                 dnew <- apply(dprob, 1, FUN=sum, na.rm=TRUE)
                 quants <- wtd.quantile(y_seq, w=dnew)       

               ## weighted quantile of resulting distribution
                 lcl <-  wtd.quantile(y_seq, w=dnew)[2]
                 ucl <- wtd.quantile(y_seq, w=dnew)[4]
                 wtsd <- sqrt(wtd.var(y_seq, w=dnew))              
               } else{
                 wtsd <- NA
                 lcl <- NA
                 ucl <- NA
               }        
               return(data.frame("SD"=wtsd, "LCL"=lcl, "UCL"=ucl))      
             }))
         }
        } else{
          wtmean <- NA
          wtsd <- NA
          lcl <- NA
          ucl <- NA
          var_info <- data.frame("SD" = wtsd, "LCL" = lcl, "UCL" = ucl)
        } 
      ## return data frame
      df <- data.frame("Iteration"=itervec[ii], "Population"=scenario[ss], "Approach"=paste0("Stacking",bivar_approach[aa]), "Model"=models[mm], "Year"=years, "MLE"=wtmean, "SD"=unlist(var_info[,"SD"]), "LCL"=unlist(var_info[,"LCL"]), "UCL"=unlist(var_info[,"UCL"]), "Variable"=vals[v], "True"=tval)
      
      wtmean <- NULL
      wtsd <- NULL
      lcl <- NULL
      ucl <- NULL
      y_i <- NULL
      ysd_i <- NULL
      y_j <- NULL
      ysd_j <- NULL
      Prob_j <- NULL
      return(df)
    })
    byVal <- do.call(rbind, byVal)
    return(byVal)
    })
    byMod <- do.call(rbind, byMod)
    return(byMod)
  })
  byIter <- do.call(rbind, byIter)
  return(byIter)
  })
  byApproach <- do.call(rbind, byApproach)
  return(byApproach)
})
res_stack <- do.call(rbind, res_stack)


# res_all <- rbind(res, res_stack)

# res_all2 <- rbind(res_all, res_avg)
# saveRDS(res_all2, file.path(res_dir, "summary_RE_Cover.rds"))

res_all <- readRDS(file.path(res_dir, "summary_RE_Cover.rds"))
itervec <- 1:100
equil_dir <- file.path(res_dir, "equil")
var_dir <- file.path(res_dir, "variable")
bivar_approach <- c("1D")#, "2D")
models <- c("LBSPR","LIME")
scenario <- c("Equilibrium","Variable")



## add true initial SPR value
initSPR <- lapply(1:length(scenario), function(ss){
    byIter <- lapply(itervec, function(ii){
      if(scenario[ss]=="Equilibrium") iterpath <- file.path(equil_dir, itervec[ii])
      if(scenario[ss]=="Variable") iterpath <- file.path(var_dir, itervec[ii])
      true <- readRDS(file.path(iterpath, "True.rds"))
          SPR1 <- true$SPR_t[1]
          df <- data.frame("Iteration"=itervec[ii], "Population"=scenario[ss], "InitSPR"=SPR1)
          return(df)
    })
    byIter <- do.call(rbind, byIter)
    return(byIter)
})
initSPR <- do.call(rbind, initSPR)

res_all3 <- inner_join(res_all, initSPR)

summary <- res_all3 %>% 
          mutate(RE = (MLE - True)/True) %>%
          mutate(Converge = ifelse(is.na(SD), 0, 1)) %>%
          mutate(Cover = ifelse(True >= LCL & True <= UCL, 1, 0)) %>%
          mutate(CoverAdj =  ifelse(Converge ==1, Cover, NA)) %>%
          filter(Approach %in% c("Stacking2D", "Average2D") == FALSE)

pinit <- ggplot(summary %>% filter(Year==max(Year))) +
         geom_point(aes(x=InitSPR, y=RE, color=Population), cex=1.2) +
         facet_grid(Model+Population~Approach) +
         guides(color=FALSE) +
         mytheme() +
         xlab("Initial SPR") +
         geom_hline(yintercept=0)
ggsave(file.path(fig_dir, "RE_byInitialSPR.png"), height=8, width=13, pinit)

sub <- summary %>% select(c(Year, CoverAdj, InitSPR, Model, Population, Approach)) %>% na.omit()
pinit2 <- ggplot(sub %>% filter(Year==max(Year))) +
         geom_violin(aes(x=as.factor(CoverAdj), y=InitSPR, color=CoverAdj)) +
         facet_grid(Model+Population~Approach) +
         guides(color=FALSE) +
         mytheme() +
         geom_hline(yintercept=0.5) +
         xlab("True SPR within 50% CI") +
         ylab("Initial SPR")
ggsave(file.path(fig_dir, "Cover_byInitialSPR.png"), height=8, width=13, pinit2)

lastyr <- summary %>% filter(Model=="LIME") %>% filter(Population=="Variable") %>% filter(Approach=="Stacking2D")
plot(x=1,y=1,type="n", xlim=c(min(itervec),max(itervec)), ylim=c(-1,1))
for(i in 1:length(itervec)){
  points(x=itervec[i], y=median(lastyr$RE[which(lastyr$Iteration <= itervec[i])], na.rm=TRUE))
}

re <- summary %>% filter(Year==max(Year, na.rm=TRUE)) %>% filter(Converge == 1) %>%
      group_by(Approach, Population, Model, Variable) %>%
      select(Approach, Population, Model, Variable, RE) %>%
      summarise_all(funs(mre=median(.,na.rm=TRUE), mare=median(abs(.),na.rm=TRUE)))

cover <- summary %>% filter(Year==max(Year, na.rm=TRUE)) %>% 
       group_by(Approach, Population, Model, Variable) %>%
       select(Approach, Population, Model, Variable, CoverAdj) %>%
       summarise_all(funs(Coverage=sum(., na.rm=TRUE), Converge=length(which(is.na(.)==FALSE)))) %>%
       mutate(PropCover = ifelse(Converge == 0, 0, Coverage / Converge))

pre <- ggplot(summary %>% filter(Variable=="SPR") %>% filter(Year==max(Year, na.rm=TRUE)) %>% filter(Converge == 1)) +
      geom_violin(aes(x=Approach, y=RE, fill=Population), scale="width", draw_quantiles=0.5) + 
      geom_hline(aes(yintercept=0), lwd=1.5) + 
      facet_grid(Model~.) +
      mytheme() +
      ylab("Relative error") + xlab("Life history approach")
  ggsave(file.path(fig_dir, "RE_SPR.png"), height=8, width=13, pre)

pcover <- ggplot(cover %>% filter(Variable=="SPR")) +
    geom_point(aes(x=Approach, y=PropCover*100, color=Population, fill=Population, shape=Population), cex=6, alpha=0.7) +
    geom_hline(aes(yintercept=50), lwd=1.5) +
    facet_grid(Model~.) + 
    mytheme() +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
    ylab("Interval Coverage") + xlab("Life history approach") +
    ylim(c(0,100))
ggsave(file.path(fig_dir, "Cover_SPR.png"), height=8, width=10, pcover)

pconverge <- ggplot(cover %>% filter(Variable=="SPR")) +
    geom_point(aes(x=Approach, y=Converge, color=Population, fill=Population, shape=Population), cex=6, alpha=0.7) +
    geom_hline(aes(yintercept=50), lwd=1.5) +
    facet_grid(Model~.) + 
    mytheme() +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
    ylab("Convergence") + xlab("Life history approach") +
    ylim(c(0,100))
ggsave(file.path(fig_dir, "Converge_SPR.png"), height=8, width=10, pconverge)















#######################
## example process
#######################
find <- NULL
for(i in 1:length(itervec)){
  iterpath <- file.path(equil_dir, itervec[i])
  if(file.exists(file.path(iterpath, "res_BestCase_LIME.rds"))==FALSE) find <- c(find, itervec[i])
}
true <- readRDS(file.path(equil_dir, find[1], "True.rds"))



iterpath <- file.path(equil_dir, 50)
stack1 <- readRDS(file.path(iterpath, paste0("res_stacking_LBSPR_1D.rds")))
stack2 <- readRDS(file.path(iterpath, paste0("res_stacking_LBSPR_2D.rds")))
true <- readRDS(file.path(iterpath, "True.rds"))
mres <- readRDS(file.path(iterpath, "res_Means_LBSPR.rds"))

y_seq <- seq(0.0001,0.9999,by=0.001)


#############
## 1D
#############
Dim <- 1
## find SPR estimates at nodes
years <- stack1[[1]]@Years
y_i1 <- sapply(1:length(stack1), function(x) stack1[[x]]@SPR)[length(years),]
ysd_i1 <- sapply(1:length(stack1), function(x) sqrt(stack1[[x]]@Vars[,"SPR"]))[length(years),]

## Use approx() to do 1D smoother of SPR_i across nodes
interp_h1 <- approx(x=Nodes1_i[,1], y=y_i1, n=100)
y_j1 <- interp_h1$y

    bad <- which(is.na(ysd_i1))
    if(length(bad)==0){
      nodes_inp <- Nodes1_i[,1]
      ysd_inp <- ysd_i1
    }
    if(length(bad)>0){
      nodes_inp <- Nodes1_i[-bad,1]
      ysd_inp <- ysd_i1[-bad]
    }
interp_k1 <- approx(x=nodes_inp, y=ysd_inp, n=100)
ysd_j1 <- interp_k1$y

# translate Grid to 4 parameters and calculate FishLife probability for each
Grid_j1 <- interp_h1$x
Param_j1 <- t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j1) )
colnames(Param_j1) <- names(Mean)
Prob_j1 <- dmvnorm( Param_j1, mean=Mean, sigma=Cov ) 

# Calculate mean
wtmean1 <- weighted.mean( y_j1, w=Prob_j1, na.rm=TRUE)

## density at each interpolated point - interpolated mean and sd from assessment

dstack1 <- sapply(1:length(y_j1), function(x) dnorm(y_seq, mean=y_j1[x], sd=ysd_j1[x]))
## density at interpolated point weighted by FishLife probability
wdstack1 <- dstack1 * outer( rep(1,nrow(dstack1)), Prob_j1 )
# summary of weighted distribution
dnew1 <- rowSums(wdstack1)
quants1 <- wtd.quantile(y_seq, w=dnew1)       


up <- unique(Prob_j1)
oup <- up[order(up)]
cup <- rev(topo.colors(length(up)))
wcol <- sapply(1:length(Prob_j1), function(x) cup[which(oup == Prob_j1[x])])

png(file.path(fig_dir, "FishLife_Prob_Scale.png"), height=7, width=2, units="in", res=200)
image(x=1, y=oup, z=t(matrix(oup)), col=cup, axes=FALSE, ann=F, xaxs="i", yaxs="i")
# axis(2, cex.axis=2)
dev.off()

png(file.path(fig_dir, "Example_stacking_1D_Steps1-2.png"), height=7, width=15, units="in", res=200)
par(mfrow=c(1,2), mar=c(5,5,1,1))

plot(x=interp_h1$x, y=Prob_j1, col=wcol, pch=19, xaxt="n", ylab="Likelihood", xlab="Nodes", cex.axis=2, cex.lab=2.5, cex=1.5)
points(x=Nodes1_i, y=dmvnorm(Param1_i, mean=Mean, sigma=Cov), pch=19, cex=3)
axis(1, at=seq(min(interp_h1$x),max(interp_h1$x),length.out=4), labels=1:4, cex.axis=2)

plot(interp_h1$x, y=y_j1, ylim=c(0,1), pch=19, cex=1.5, col=wcol, xaxt="n", xlab="Nodes", ylab="SPR", cex.axis=2, cex.lab=2.5, xaxs="i", yaxs="i")
points(x=Nodes1_i, y=y_i1, cex=3, pch=19, xpd=NA)
axis(1, at=seq(min(interp_h1$x),max(interp_h1$x),length.out=4), labels=1:4, cex.axis=2)

dev.off()

png(file.path(fig_dir, "Example_stacking_1D_Step3.png"), height=7, width=18, units="in", res=200)
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(x=1, y=1, type="n", xlim=c(0,1), ylim=c(0,quantile(wdstack1,1)), cex.axis=2, xlab="SPR", ylab="Density", cex.lab=2.5, xaxs="i", yaxs="i")
for(i in 1:ncol(dstack1)){
  polygon(x=c(y_seq, rev(y_seq)), y=c(rep(0,length(y_seq)),rev(wdstack1[,i])), col=paste0(substring(wcol[i],1,7),"50"))
}
dev.off()

png(file.path(fig_dir, "Example_stacking_1D_Step4.png"), height=7, width=18, units="in", res=200)
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(x=1, y=1, type="n", xlim=c(0,1), ylim=c(0, max(dnew1)*1.1), cex.axis=2, xlab="SPR", ylab="Density", cex.lab=2.5, xaxs="i", yaxs="i")
polygon(x=c(y_seq, rev(y_seq)), y=c(rep(0, length(y_seq)), rev(dnew1)), col=gray(0.7), border=NA)
abline(v=wtmean1, col="blue",lwd=4)
abline(v=true$SPR, lwd=4)
abline(v=quants1[2], col="blue", lwd=2, lty=2)
abline(v=quants1[4], col="blue", lwd=2, lty=2)
abline(v=mres@SPR[length(mres@SPR)], col="red", lwd=2)
legend("topright", legend=c("Truth", "Means", "Stacking", "Stacking 50% interval"), col=c("black", "red", "blue", "blue"), lwd=c(4,2,4,2), lty=c(1,1,1,2), cex=1.7)
dev.off()

#########################
# Step 2:  Run LBSPR using each row i of Param_i, and extract SPR_i for each
#########################

  ##################################
  ## LBSPR - 1D
  ##################################

  ## Simulate true data with species-level information
  ## Run all models with species-level distributions
  ncores <- 5
  cl <- makeCluster(ncores)
  registerDoParallel(cl)      
  start <- Sys.time()
  foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
    runstack(savedir=res_dir, 
          nodes=Param1_i, 
          param=c("K","M","Loo","Lm"), 
          mean=Mean, 
          cov=Cov, 
          taxon="Species",
          dim="1D",
          simulation=TRUE, 
          iter=loop, 
          seed=loop, 
          Fscenario="equil", 
          rewrite=TRUE, 
          model="LBSPR", 
          sim_model="LIME",
          Nyears=20)
  end <- Sys.time() - start   
  stopCluster(cl)

  ##################################
  ## LBSPR - 2D
  ##################################

  ## Simulate true data with species-level information
  ## Run all models with species-level distributions
  ncores <- 5
  cl <- makeCluster(ncores)
  registerDoParallel(cl)      
  start <- Sys.time()
  foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
    runstack(savedir=res_dir, 
          nodes=Param2_i, 
          param=c("K","M","Loo","Lm"), 
          mean=Mean, 
          cov=Cov, 
          taxon="Species",
          dim="2D",
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

  ##################################
  ## LIME - 1D
  ##################################

  ## Simulate true data with species-level information
  ## Run all models with species-level distributions
  ncores <- 5
  cl <- makeCluster(ncores)
  registerDoParallel(cl)      
  start <- Sys.time()
  foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
    runstack(savedir=res_dir, 
          nodes=Param1_i, 
          param=c("K","M","Loo","Lm"), 
          mean=Mean, 
          cov=Cov, 
          taxon="Species",
          dim="1D",
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

  ##################################
  ## LIME - 2D
  ##################################

  ## Simulate true data with species-level information
  ## Run all models with species-level distributions
  ncores <- 5
  cl <- makeCluster(ncores)
  registerDoParallel(cl)      
  start <- Sys.time()
  foreach(loop=itervec, .packages=c('TMB','LIME','LBSPR','mvtnorm')) %dopar% 
    runstack(savedir=res_dir, 
          nodes=Param2_i, 
          param=c("K","M","Loo","Lm"), 
          mean=Mean, 
          cov=Cov, 
          taxon="Species",
          dim="2D",
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

#########################
# Step 3:  Smooth SPR_i in original 1/2 dimensional space
#########################

model <- c("LBSPR","LIME")
dim <- c("1D", "2D")

mods <- expand.grid("Model"=model, "Dim"=dim)

res_stack <- lapply(1:nrow(mods), function(m){
# for(m in 1:nrow(mods)){
  byIter <- lapply(1:length(itervec), function(i){
  # for(i in 1:length(itervec)){
    if(mods[m,"Model"]=="LBSPR") vals <- "SPR"
    if(mods[m,"Model"]=="LIME") vals <- c("SPR","BB0")
    byVal <- lapply(1:length(vals), function(v){
    # for(v in 1:length(vals)){
      if(vals[v]=="SPR") y_seq <- seq(0.0001,0.9999,by=0.001)
      if(vals[v]=="BB0") y_seq <- seq(0,4,by=0.001)
        ## read files
        if(file.exists(file.path(res_dir, itervec[i], paste0("res_stacking_", mods[m,"Model"], "_Species_", mods[m,"Dim"], ".rds")))){
          stack <- readRDS(file.path(res_dir, itervec[i], paste0("res_stacking_", mods[m,"Model"], "_Species_", mods[m,"Dim"], ".rds")))
          if(mods[m,"Model"]=="LBSPR"){
            if(vals[v]=="SPR"){
              y_i <- sapply(1:length(stack), function(x) stack[[x]]@SPR[length(stack[[x]]@SPR)])
              ysd_i <- sapply(1:length(stack), function(x) sqrt(stack[[x]]@Vars[,"SPR"][length(stack[[x]]@Vars[,"SPR"])]))
            }
            NLL_i <- sapply(1:length(stack), function(x) stack[[x]]@NLL)
          }
          if(mods[m,"Model"]=="LIME"){
            if(vals[v]=="SPR"){
              y_i <- sapply(1:length(stack), function(x) stack[[x]]$Report$SPR_t[length(stack[[x]]$Report$SPR_t)])
              ysd_i <- sapply(1:length(stack), function(x){
                sdrep <- summary(stack[[x]]$Sdreport)
                sd <- sdrep[which(rownames(sdrep)=="SPR_t"),2][length(which(rownames(sdrep)=="SPR_t"))]
                return(sd)
              })
            }
            if(vals[v]=="BB0"){
              y_i <- sapply(1:length(stack), function(x) stack[[x]]$Report$D_t[length(stack[[x]]$Report$D_t)])
              ysd_i <- sapply(1:length(stack), function(x){
                sdrep <- summary(stack[[x]]$Sdreport)
                sd <- sdrep[which(rownames(sdrep)=="lD_t"),2][length(which(rownames(sdrep)=="lD_t"))]
              })
              NLL_i <- sapply(1:length(stack), function(x) stack[[x]]$Report$jnll)
            }
          }
        } else {
          stack <- NULL
          y_i <- NULL
          ysd_i <- NULL
          NLL_i <- NULL
        }
        ## smooth SPR_i in 1/2D space
        if(all(is.null(y_i))==FALSE){
          if(mods[m,"Dim"]=="1D"){
            Dim <- 1
            # Step 3A: Use approx() to do 1D smoother of SPR_i across nodes
            interp_h <- approx(x=Nodes1_i[,1], y=y_i, n=100)
            y_j <- interp_h$y

            bad <- which(is.na(ysd_i))
            if(length(bad)==0){
              nodes_inp <- Nodes1_i[,1]
              ysd_inp <- ysd_i
            }
            if(length(bad)>0){
              nodes_inp <- Nodes1_i[-bad,1]
              ysd_inp <- ysd_i[-bad]
            }
            interp_k <- approx(x=nodes_inp, y=ysd_inp, n=100)
            ysd_j <- interp_k$y

            # Step 3B: translate Grid to 4 parameters and calculate FishLife probability for each
            Grid_j <- interp_h$x
            Param_j <- t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j) )
            colnames(Param_j) <- names(Mean)
            Prob_j <- dmvnorm( Param_j, mean=Mean, sigma=Cov ) 

            # Step 3C: Calculate mean
            wtmean <- weighted.mean( y_j, w=Prob_j, na.rm=TRUE) 

            if(wtmean <= max(y_seq) & all(is.na(ysd_j))==FALSE){
            ## density at each interpolated point - interpolated mean and sd from assessment
            dstack <- sapply(1:length(y_j), function(i){
              d <- dnorm(y_seq, mean=y_j[i], sd=ysd_j[i])
              return(d)
            })
            ## density at interpolated point weighted by FishLife probability
            wdstack <- dstack * outer( rep(1,nrow(dstack)), Prob_j )
            # summary of weighted distribution
            dnew <- rowSums(wdstack)
            quants <- wtd.quantile(y_seq, w=dnew)     

              ## weighted quantile of resulting distribution
              lcl <-  wtd.quantile(y_seq, w=dnew)[2]
              ucl <- wtd.quantile(y_seq, w=dnew)[4]
              wtsd <- sqrt(wtd.var(y_seq, w=dnew))              
            } else{
              wtsd <- NA
              lcl <- NA
              ucl <- NA
            }
          }
          if(mods[m,"Dim"]=="2D"){  
            Dim <- 2    

            # Step 3A: Use akima::interp to do 2D smoother of SPR_i across nodes
            interp_h <- akima::interp( x=Nodes2_i[,1], y=Nodes2_i[,2], z=y_i, nx=100, ny=100, duplicate="median")
            y_j <- interp_h$z 

            bad <- which(is.na(ysd_i))
            if(length(bad)==0){
              nodes_inpx <- Nodes2_i[,1]
              nodes_inpy <- Nodes2_i[,2]
              ysd_inp <- ysd_i
            }
            if(length(bad)>0){
              nodes_inpx <- Nodes2_i[-bad,1]
              nodes_inpy <- Nodes2_i[-bad,1]
              ysd_inp <- ysd_i[-bad]
            }
            interp_k <- akima::interp( x=nodes_inpx, y=nodes_inpy, z=ysd_inp, nx=100, ny=100, duplicate="median")
            ysd_j <- interp_k$z     

            # Step 3B: translate Grid to 4 parameters and calculate FishLife probability for each
           Grid_j <- expand.grid(interp_h$x, interp_h$y)
           Param_j <- t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j) )
           colnames(Param_j) <- names(Mean)
           Prob_j <- matrix(dmvnorm( Param_j, mean=Mean, sigma=Cov ), nrow=100)           
        
           wtmean <- weighted.mean(y_j, w=matrix(Prob_j, nrow=100), na.rm=TRUE)

          if(wtmean <= max(y_seq) & all(is.na(ysd_j))==FALSE){
           dstack <- lapply(1:length(y_seq), function(x){
            sub <- sapply(1:nrow(y_j), function(i){
              sub2 <- sapply(1:ncol(y_j), function(j){
                d <- dnorm(y_seq[x], mean=y_j[i,j], sd=ysd_j[i,j])
                return(d)
              })
              return(sub2)
            })
            return(sub)
           })
           darray <- array(NA, dim=c(dim(y_j),length(y_seq)))
           for(i in 1:length(dstack)){
            darray[,,i] <- dstack[[i]]
           }
           darray2 <- darray * outer( Prob_j, rep(1,length(y_seq)) )
           dnew <- apply(darray2, 3, FUN=sum, na.rm=TRUE)       

           wtsd <- sqrt(wtd.var(y_seq, w=dnew))
           lcl <- wtd.quantile(y_seq, w=dnew)[2]
           ucl <- wtd.quantile(y_seq, w=dnew)[4]     
          }  else {
            wtsd <- NA
            lcl <- NA
            ucl <- NA
          }
         }
        } else{
          wtmean <- NA
          wtsd <- NA
          lcl <- NA
          ucl <- NA
        } 
      ## return data frame
      df <- data.frame("Iteration"=itervec[i], "Scenario"=paste0("Ensemble",mods[m,"Dim"]), "Model"=mods[m,"Model"], "MLE"=wtmean, "SD"=wtsd, "LCL"=lcl, "UCL"=ucl, "Variable"=vals[v])
      wtmean <- NULL
      wtsd <- NULL
      lcl <- NULL
      ucl <- NULL
      y_i <- NULL
      ysd_i <- NULL
      y_j <- NULL
      ysd_j <- NULL
      Prob_j <- NULL
      return(df)
    })
    byVal <- do.call(rbind, byVal)
    return(byVal)
  })
  byMod <- do.call(rbind, byIter)
  return(byMod)
})
res_stack <- do.call(rbind, res_stack)
res_stack_sub <- res_stack %>% filter(Scenario!="Ensemble2D")
saveRDS(res_stack_sub, file.path(res_dir, "stack_results_TrueMeans1D.rds"))
saveRDS(res_stack, file.path(res_dir, "stack_results.rds"))

res_stack$Iteration[which(res_stack$Model=="LBSPR" & res_stack$Scenario=="Ensemble2D")] <- 1:100
res_stack$Iteration[which(res_stack$Model=="LIME" & res_stack$Variable=="SPR" & res_stack$Scenario=="Ensemble2D")] <- 1:100
res_stack$Iteration[which(res_stack$Model=="LIME" & res_stack$Variable=="BB0" & res_stack$Scenario=="Ensemble2D")] <- 1:100
saveRDS(res_stack, file.path(res_dir, "stack_results.rds"))


find_true <- lapply(1:length(itervec), function(i){
  true <- readRDS(file.path(res_dir, itervec[i], "True.rds"))
  vals <- c("SPR","BB0")
  byVal <- lapply(1:length(vals), function(v){
    if(vals[v]=="SPR") df <- data.frame("Iteration"=itervec[i], "True"=true$SPR_t[length(true$SPR_t)], "Variable"="SPR")   
    if(vals[v]=="BB0") df <- data.frame("Iteration"=itervec[i], "True"=true$D_t[length(true$D_t)], "Variable"="BB0") 
    return(df)
  })
  byVal <- do.call(rbind, byVal)
  return(byVal)
})
find_true <- do.call(rbind, find_true)

res_means <- lapply(1:length(model), function(m){
# for(m in 1:length(model)){
  byIter <- lapply(1:length(itervec), function(i){
  # for(i in 1:length(itervec)){
    if(model[m]=="LBSPR") vals <- "SPR"
    if(model[m]=="LIME") vals <- c("SPR", "BB0")
      byVal <- lapply(1:length(vals), function(v){
      # for(v in 1:length(vals)){
      if(vals[v]=="SPR") y_seq <- seq(0.0001,0.9999,by=0.001)
      if(vals[v]=="BB0") y_seq <- seq(0,4,by=0.001)

        ## read files
        if(file.exists(file.path(res_dir, itervec[i], paste0("res_Means_", model[m], "_Species.rds")))){
          mres <- readRDS(file.path(res_dir, itervec[i], paste0("res_Means_", model[m], "_Species.rds")))
          if(model[m]=="LBSPR"){
            if(vals[v]=="SPR"){
                mle <- mres@SPR[length(mres@SPR)]
                sd <- sqrt(mres@Vars[,"SPR"][length(mres@Vars[,"SPR"])])
            }
          }
          if(model[m]=="LIME"){
            sdrep <- summary(mres$Sdreport)
            if(vals[v]=="SPR"){
                mle <- mres$Report$SPR_t[length(mres$Report$SPR_t)]
                sd <- sdrep[which(rownames(sdrep)=="SPR_t"),2][length(which(rownames(sdrep)=="SPR_t"))]              
            }
            if(vals[v]=="BB0"){
                mle <- mres$Report$D_t[length(mres$Report$D_t)]
                sd <- sdrep[which(rownames(sdrep)=="lD_t"),2][length(which(rownames(sdrep)=="lD_t"))]              
            }
          }
          if(mle <= max(y_seq) & is.na(sd)==FALSE){
            d <- dnorm(y_seq, mean=mle, sd=sd)
            if(sum(d)>0){
              lcl <- wtd.quantile(y_seq, w=d)[2]
              ucl <- wtd.quantile(y_seq, w=d)[4]
            } else{
              lcl <- NA
              ucl <- NA
            }
          } else{
            lcl <- NA
            ucl <- NA
          }
        } else {
          mle <- NA
          sd <- NA
          lcl <- NA
          ucl <- NA
        } 
        ## return data frame
        df <- data.frame("Iteration"=itervec[i], "Scenario"="Means", "Model"=model[m], "MLE"=mle, "SD"=sd, "LCL"=lcl, "UCL"=ucl, "Variable"=vals[v])
        return(df)   
      }) 
      byVal <- do.call(rbind, byVal)
      return(byVal)     
  })
  byMod <- do.call(rbind, byIter)
  return(byMod)
})
res_means <- do.call(rbind, res_means)

res_best <- lapply(1:length(model), function(m){
  byIter <- lapply(1:length(itervec), function(i){
    if(model[m]=="LBSPR") vals <- "SPR"
    if(model[m]=="LIME") vals <- c("SPR", "BB0")
      byVal <- lapply(1:length(vals), function(v){
        if(vals[v]=="SPR") y_seq <- seq(0.0001,0.9999,by=0.001)
        if(vals[v]=="BB0") y_seq <- seq(0,4,by=0.001)

        ## read files
        if(file.exists(file.path(res_dir, itervec[i], paste0("res_IterTrue_", model[m], ".rds")))){
          ires <- readRDS(file.path(res_dir, itervec[i], paste0("res_IterTrue_", model[m], ".rds")))
          if(model[m]=="LBSPR"){
            if(vals[v]=="SPR"){
                mle <- ires@SPR[length(ires@SPR)]
                sd <- sqrt(ires@Vars[,"SPR"][length(ires@Vars[,"SPR"])])
            }
          }
          if(model[m]=="LIME"){
            sdrep <- summary(ires$Sdreport)
            if(vals[v]=="SPR"){
                mle <- ires$Report$SPR_t[length(ires$Report$SPR_t)]
                sd <- sdrep[which(rownames(sdrep)=="SPR_t"),2][length(which(rownames(sdrep)=="SPR_t"))]              
            }
            if(vals[v]=="BB0"){
                mle <- ires$Report$D_t[length(ires$Report$D_t)]
                sd <- sdrep[which(rownames(sdrep)=="lD_t"),2][length(which(rownames(sdrep)=="lD_t"))]              
            }
          }
          if(mle <= max(y_seq) & is.na(sd)==FALSE){
            d <- dnorm(y_seq, mean=mle, sd=sd)
            if(sum(d)>0){
              lcl <- wtd.quantile(y_seq, w=d)[2]
              ucl <- wtd.quantile(y_seq, w=d)[4]
            } else{
              lcl <- NA
              ucl <- NA
            }
          } else{
            lcl <- NA
            ucl <- NA
          }
        } else {
          mle <- NA
          sd <- NA
          lcl <- NA
          ucl <- NA
        } 
        ## return data frame
        df <- data.frame("Iteration"=itervec[i], "Scenario"="BestCase", "Model"=model[m], "MLE"=mle, "SD"=sd, "LCL"=lcl, "UCL"=ucl, "Variable"=vals[v])
        return(df)   
      }) 
      byVal <- do.call(rbind, byVal)
      return(byVal)     
  })
  byMod <- do.call(rbind, byIter)
  return(byMod)
})
res_best <- do.call(rbind, res_best)

join1 <- rbind(res_means, res_best) %>% full_join(find_true)
join2 <- full_join(res_stack, find_true)
res <- full_join(join1, join2) %>% 
        mutate(RE=(MLE - True)/True) %>% 
        mutate(Cover50=ifelse((True >= LCL) & (True <= UCL), 1, 0))

re <- res %>%
      group_by(Scenario, Model, Variable) %>%
      select(Scenario, Model, Variable, RE) %>%
      summarise_all(funs(mre=median(.,na.rm=TRUE), mare=median(abs(.),na.rm=TRUE)))

cover <- res %>% 
       group_by(Scenario, Model, Variable) %>%
       select(Scenario, Model, Variable, Cover50) %>%
       summarise_all(funs(Coverage=sum(., na.rm=TRUE), Converge=length(which(is.na(Cover50)==FALSE)))) %>%
       mutate(PropCover = Coverage / Converge)

pre <- ggplot(res %>% filter(Variable=="SPR")) +
      geom_violin(aes(x=Scenario, y=RE, color=Scenario, fill=Scenario)) + 
      geom_hline(aes(yintercept=0), lwd=1.5) + 
      facet_grid(Model~.) +
      mytheme() +
      ylab("Relative error") + xlab("Modeling approach")
  ggsave(file.path(fig_dir, "RE_SPR.png"), height=8, width=13, pre)

pre2 <- ggplot(res %>% filter(Variable=="BB0")) +
      geom_violin(aes(x=Scenario, y=RE, color=Scenario, fill=Scenario)) + 
      geom_hline(aes(yintercept=0), lwd=1.5) + 
      facet_grid(Model~.) +
      mytheme() +
      ylab("Relative error") + xlab("Modeling approach")
   ggsave(file.path(fig_dir, "RE_BB0.png"), height=4, width=12, pre2)

pcover <- ggplot(cover %>% filter(Variable=="SPR")) +
    geom_point(aes(x=Scenario, y=PropCover, color=Scenario, fill=Scenario), cex=6) +
    geom_hline(aes(yintercept=0.5), lwd=1.5) +
    facet_grid(Model~.) + 
    mytheme() +
    ylab("Interval Coverage") + xlab("Modeling Scenario") +
    ylim(c(0,1))
  ggsave(file.path(fig_dir, "Cover_SPR.png"), height=8, width=13, pcover)

pcover2 <- ggplot(cover %>% filter(Variable=="BB0")) +
    geom_point(aes(x=Scenario, y=PropCover, color=Scenario, fill=Scenario), cex=6) +
    geom_hline(aes(yintercept=0.5), lwd=1.5) +
    facet_grid(Model~.) + 
    mytheme() +
    ylab("Interval Coverage") + xlab("Modeling Scenario") +
    ylim(c(0,1))
  ggsave(file.path(fig_dir, "Cover_BB0.png"), height=4, width=12, pcover2)


## demonstrate one iter
Dim <- 1
true1 <- readRDS(file.path(res_dir, 1, "True.rds"))
stack1 <- readRDS(file.path(res_dir, 1, "res_stacking_LBSPR_Species_1D.rds"))
SPR_i <- sapply(1:length(stack1), function(x) stack1[[x]]@SPR)
SPRsd_i <- sapply(1:length(stack1), function(x) stack1[[x]]@Vars[,"SPR"])
mres1 <- readRDS(file.path(res_dir, 1, "res_Means_LBSPR_Species.rds"))

 # Step 3A: Use approx() to do 1D smoother of SPR_i across nodes
 interp_h <- approx(x=Nodes1_i[,1], y=SPR_i, n=100)
 SPR_j <- interp_h$y

 interp_k <- approx(x=Nodes1_i[,1], y=SPRsd_i, n=100)
 SPRsd_j <- interp_k$y

 # Step 3B: translate Grid to 4 parameters and calculate FishLife probability for each
 Grid_j <- interp_h$x
 Param_j <- t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j) )
 colnames(Param_j) <- names(Mean)
 Prob_j <- dmvnorm( Param_j, mean=Mean, sigma=Cov )  

 wtmean <- weighted.mean(SPR_j, w=Prob_j)


  up <- unique(Prob_j)
  oup <- up[order(up)]
  coup <- rev(topo.colors(n=length(oup)))
  pl <- Prob_j
  cpl <- sapply(1:length(pl), function(x) coup[which(oup == pl[[x]])]) 
  spr_seq <- seq(0,1,by=0.001)

  uspr <- unique(SPR_j)
  ouspr <- uspr[order(uspr)]
  cuspr <- rev(topo.colors(n=length(ouspr)))
  spr <- SPR_j
  cspr <- sapply(1:length(spr), function(x) cuspr[which(ouspr == spr[[x]])])

  ## smooth across nodes
  png(file.path(fig_dir, "Interpolated_SPR_LBSPR_1D.png"), height=6, width=6, units="in", res=200)
  par(mfrow=c(1,1), mar=c(5,5,2,2))
  plot(SPR_j, col="#AAAAAA", pch=19, cex=1, xlab="Interpolated points between nodes", ylab="SPR", cex.axis=2, cex.lab=2, ylim=c(0,1))
  par(new=TRUE)
  plot(SPR_i, col="black", pch=19, cex=2, axes=F, ann=F, ylim=c(0,1))
  dev.off()

  png(file.path(fig_dir, "FishLife_prob_1D.png"), height=6, width=6, units="in", res=200)
  par(mfrow=c(1,1), mar=c(5,5,2,2))
  plot(pl, col=cpl, pch=19, xlab="Interpolated points", ylab="FishLife probability", cex.axis=2, cex.lab=2)
  dev.off()

  png(file.path(fig_dir, "SPR_weights_1D.png"), height=6, width=6, units="in", res=200)
  par(mfrow=c(1,1), mar=c(5,5,2,2))
  plot(SPR_j, col=cpl, pch=19, cex.axis=2, cex.lab=2, xlab="Interpolated points", ylab="SPR", ylim=c(0,1))
  dev.off()

  png(file.path(fig_dir, "SPR_weights_wEstTrue_1D.png"), height=6, width=6, units="in", res=200)
  par(mfrow=c(1,1), mar=c(5,5,2,2))
  plot(SPR_j, col=cpl, pch=19, cex.axis=2, cex.lab=2, xlab="Interpolated points", ylab="SPR", ylim=c(0,1))
  abline(h=wtmean, col="blue", lwd=2, lty=2)
  abline(h=true1$SPR, col="red", lwd=2, lty=2)
  dev.off()


  ## density at each interpolated point - interpolated mean and sd from assessment
  dstack <- sapply(1:length(SPR_j), function(i){
    d <- dnorm(spr_seq, mean=SPR_j[i], sd=SPRsd_j[i])
    return(d)
  })
  ## density at interpolated point weighted by FishLife probability
  wdstack <- dstack * outer( rep(1,nrow(dstack)), Prob_j )
  # summary of weighted distribution
  dnew <- rowSums(wdstack)
  quants <- wtd.quantile(spr_seq, w=dnew)

png(file.path(fig_dir, "Interpolated_distr_SPR_1D.png"), width=10, height=3.8, units="in", res=200)
  par(mfrow=c(1,1), mar=c(5,5,2,2))
  plot(x=1,y=1,type="n",xlim=c(0,1),ylim=c(0,100), xlab="SPR", ylab="Density", cex.axis=2, cex.lab=2)
  for(i in 1:length(SPR_j)){
    polygon(x=c(spr_seq, rev(spr_seq)), y=c(rep(0,length(spr_seq)), rev(dstack[,i])), col=cpl[i])
  }
 plotrix::color.legend(xl=0.03,xr=0.10,yb=50,yt=90, legend=round(seq(min(up),max(up),length.out=5)), rect.col=rev(topo.colors(5)), gradient="y")
  text(x=0.06, y=97, "Weight", cex=1.2)
dev.off()

png(file.path(fig_dir, "Stacked_distr_SPR_1D.png"), width=10, height=3.8, units="in", res=200)
  par(mfrow=c(1,1), mar=c(5,5,2,2))
  plot(x=1,y=1,type="n",xlim=c(0,1),ylim=c(0,max(dnew)*1.1), xlab="SPR", ylab="Density", cex.lab=2, cex.axis=2)
  polygon(x=c(spr_seq, rev(spr_seq)), y=c(rep(0,length(spr_seq)), rev(dnew)), col="#AAAAAA90")
  abline(v=true1$SPR_t, lty=2, col="red", lwd=3)
  abline(v=wtmean, col="blue", lty=2, lwd=3)
  abline(v = quantile(quants)[2], col="goldenrod", lty=2, lwd=3)
  abline(v=quantile(quants)[4], col="goldenrod", lty=2, lwd=3)
  abline(v=mres1@SPR, col="forestgreen", lty=2, lwd=3)
dev.off()



Dim <- 2
true1 <- readRDS(file.path(res_dir, 1, "True.rds"))
stack2 <- readRDS(file.path(res_dir, 1, "res_stacking_LBSPR_Species_2D.rds"))
SPR_i <- sapply(1:length(stack2), function(x) stack2[[x]]@SPR)
SPRsd_i <- sapply(1:length(stack2), function(x) sqrt(stack2[[x]]@Vars[,"SPR"]))

 interp_h <- akima::interp(x=Nodes2_i[,1], y=Nodes2_i[,2], z=SPR_i, nx=100, ny=100, duplicate="median")
 SPR_j <- interp_h$z

png(file.path(fig_dir, "Interpolated_SPR_LBSPR_2D.png"), height=6, width=6, units="in", res=200)
image(SPR_j, col=rev(heat.colors(100)), cex.axis=2, cex.lab=1.2)
par(new=TRUE)
plot(x=Nodes2_i[,1], y=Nodes2_i[,2], col="#AAAAAA", pch=19, cex=1.3, axes=F, ann=F)
dev.off()

interp_k <- akima::interp(x=Nodes2_i[,1], y=Nodes2_i[,2], z=SPRsd_i, nx=100, ny=100, duplicate="median")
 SPRsd_j <- interp_k$z

  # Step 3B: translate Grid to 4 parameters and calculate FishLife probability for each
 Grid_j <- expand.grid(interp_h$x, interp_h$y)
 Param_j <- t( Mean + Eigen$vectors[,1:Dim,drop=FALSE] %*% Diag(sqrt(Eigen$values[1:Dim])) %*% t(Grid_j) )
 colnames(Param_j) <- names(Mean)
 Prob_j <- matrix(dmvnorm( Param_j, mean=Mean, sigma=Cov ), nrow=100) 

png(file.path(fig_dir, "FishLifeProbs_LBSPR_2D.png"), height=6, width=6, units="in", res=200)
image(Prob_j, col=rev(topo.colors(100)), cex.axis=2, cex.lab=1.2)
dev.off()

 wtmean <- weighted.mean(SPR_j, w=matrix(Prob_j, nrow=100), na.rm=TRUE)

png(file.path(fig_dir, "SPR_Prob_2D.png"), height=6, width=6, units="in", res=200)
image(SPR_j * Prob_j, col=rev(topo.colors(100)), cex.axis=2, cex.lab=2)
par(new=TRUE)
contour(SPR_j, axes=F, ann=F, labcex=1.5)
dev.off()

 ## density at each interpolated point
 ## along SPR
 dstack <- lapply(1:length(spr_seq), function(x){
  sub <- sapply(1:nrow(SPR_j), function(i){
    sub2 <- sapply(1:ncol(SPR_j), function(j){
      d <- dnorm(spr_seq[x], mean=SPR_j[i,j], sd=SPRsd_j[i,j])
      return(d)
    })
    return(sub2)
  })
  return(sub)
 })
 darray <- array(NA, dim=c(dim(SPR_j),length(spr_seq)))
 for(i in 1:length(dstack)){
  darray[,,i] <- dstack[[i]]
 }
 darray2 <- darray * outer( Prob_j, rep(1,length(spr_seq)) )
 dnew <- apply(darray2, 3, FUN=sum, na.rm=TRUE)
 plot(x=spr_seq, y=dnew, ylim=c(0,quantile(dnew, prob=0.99)))
 abline(v=wtmean)
 quants <- wtd.quantile(spr_seq[-length(spr_seq)], w=dnew[-length(spr_seq)])
 abline(v=quants[2], col="goldenrod")
 abline(v=quants[4], col="goldenrod")
 abline(v=mres1@SPR, col="green")

png(file.path(fig_dir, "Interpolated_distr_SPR_2D.png"), width=10, height=3.8, units="in", res=200)
  par(mfrow=c(1,4), mar=c(0,0,0,0), omi=c(1,1,0.3,0.3))
  image(darray2[,,100], col=rev(topo.colors(100)))
  image(darray2[,,300], col=rev(topo.colors(100)), yaxt="n")
  image(darray2[,,600], col=rev(topo.colors(100)), yaxt="n")
  image(darray2[,,900], col=rev(topo.colors(100)), yaxt="n")
dev.off()

png(file.path(fig_dir, "Stacked_distr_SPR_2D.png"), width=10, height=3.8, units="in", res=200)
  par(mfrow=c(1,1), mar=c(5,5,2,2))
  plot(x=1,y=1,type="n",xlim=c(0,1),ylim=c(0,quantile(dnew,0.99)*1.1), xlab="SPR", ylab="Density", cex.lab=2, cex.axis=2)
  polygon(x=c(spr_seq, rev(spr_seq)), y=c(rep(0,length(spr_seq)), rev(dnew)), col="#AAAAAA90")
  abline(v=true1$SPR_t, lty=2, col="red", lwd=3)
  abline(v=wtmean, col="blue", lty=2, lwd=3)
  abline(v = quantile(quants)[2], col="goldenrod", lty=2, lwd=3)
  abline(v=quantile(quants)[4], col="goldenrod", lty=2, lwd=3)
  abline(v=mres1@SPR, col="forestgreen", lty=2, lwd=3)
dev.off()


