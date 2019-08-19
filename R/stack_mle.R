stack_mle <- function(savedir, itervec=NULL, modname, model, nodes, weights, vals){
  iout <- sapply(1:length(itervec), function(i){

	if(all(is.null(itervec)==FALSE)) iterpath <- file.path(savedir, itervec[i])
	if(all(is.null(itervec))) iterpath <- savedir

	pres <- readRDS(file.path(iterpath, paste0(modname, "_res_stacking_", model, ".rds")))
	if(file.exists(file.path(iterpath, paste0(modname, "_res_Means_", model,".rds")))){
		mres <- readRDS(file.path(iterpath, paste0(modname, "_res_Means_", model, ".rds")))
	} else { mres <- NULL}
	if(all(is.null(iter))==FALSE){
		if(file.exists(file.path(iterpath, paste0("res_IterTrue_", model, ".rds")))){
			ires <- readRDS(file.path(iterpath, paste0("res_IterTrue_", model, ".rds")))
		} else { ires <- NULL }
		true <- readRDS(file.path(iterpath, "True.rds"))
	}

	if(model=="LBSPR"){
		smle <- lapply(1:length(vals), function(y){
			if(vals[y]=="SPR"){
				fnodes <- sapply(1:length(pres), function(x) pres[[x]]@SPR)
				tval <- true$SPR
			}
			stack <- sum(fnodes * weights)
			re <- (stack - tval)/tval
			df <- data.frame(Iteration=itervec[i], "Stack"=stack, "True"=tval, "RE"=re, 
			"Value"=vals[y],"Model"="LBSPR")
			return(df)
		})
	}
	if(model=="LIME"){
		smle <- lapply(1:length(vals), function(y){
			if(vals[y]=="SPR"){
				fnodes <- sapply(1:length(pres), function(x) pres[[x]]$Report$SPR_t)
				tval <- true$SPR
			}
			if(vals[y]=="BB0"){
				fnodes <- sapply(1:length(pres), function(x) pres[[x]]$Report$D_t)
				tval <- true$D_t
			}
			stack <- sum(fnodes * weights)
			re <- (stack - tval)/tval
			df <- data.frame(Iteration=itervec[i], "Stack"=stack, "True"=tval, "RE"=re, 
			"Value"=vals[y], "Model"="LIME")
			return(df)
		})
	}

	return(smle)

  })
  iout2 <- do.call(rbind,iout)

  saveRDS(iout2, file.path(savedir, paste0(modname, "_stacking_MLE_", model, ".rds")))
  return(iout2)
}