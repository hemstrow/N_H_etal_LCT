make_AIC_comp <- function(dat, vars){
  mstats <- dat[,c("Basin", "Creek", "Relative.percent.RBT", "year", vars)]
  mstats <- na.omit(mstats)
  infs <- which(is.infinite(mstats[,vars]))
  if(length(infs) > 0){
    mstats <- mstats[-infs,]
  }
  
  AIC_list <- vector("list", 3)
  
  # model 1
  formula <- paste0(vars, "~ (1|Basin/Creek) + Relative.percent.RBT + (1|year)")
  mod_nested_ryear <- try(lmerTest::lmer(formula, mstats), silent = T)
  if(!class(mod_nested_ryear) == "try-error"){
    AIC_list[[1]] <- mod_nested_ryear
    names(AIC_list)[1] <- "nested_ryear"
  }
  
  # model 2
  formula <- paste0(vars, "~ (1|Basin/Creek) + Relative.percent.RBT + year")
  mod_nested_fyear <- try(lmerTest::lmer(formula, mstats), silent = T)
  if(!class(mod_nested_fyear) == "try-error"){
    AIC_list[[2]] <- mod_nested_fyear
    names(AIC_list)[2] <- "nested_fyear"
  }
  
  # model 3
  formula <- paste0(vars, "~ (1|Basin) + Creek  + year + Relative.percent.RBT")
  mod_rbasin_fyear <- try(lmerTest::lmer(formula, mstats), silent = T)
  if(!class(mod_rbasin_fyear) == "try-error"){
    AIC_list[[3]] <- mod_rbasin_fyear
    names(AIC_list)[3] <- "rbasin_fyear"
  }
  
  AIC_list <- rlist::list.clean(AIC_list)
  AICs <- lapply(AIC_list, AIC)
  names(AICs) <- names(AIC_list)
  AICs <- as.data.frame(AICs)
  AICs <- cbind(variable = paste(vars, sep = "_"), AICs)
  
  return(list(AIC = AICs, 
              mod_nested_ryear = mod_nested_ryear,
              mod_nested_fyear = mod_nested_fyear,
              mod_rbasin_fyear = mod_rbasin_fyear,
              dat = mstats))
}






make_pred <- function(AIC_tab){
  mstats <- AIC_tab$dat
  vars <- colnames(AIC_tab$dat)[-c(1:4)]
  
  # figure out which model to use
  ## catch errors
  not.errors <- !unlist(lapply(AIC_tab[-1][-(length(AIC_tab) - 1)], function(x) class(x) == "try-error"))
  if(sum(not.errors) < length(not.errors)){
    err <- names(which(!not.errors))
    AIC_tab[[err]] <- NULL
  }
  
  not.singulars <- !unlist(lapply(AIC_tab[-1][-(length(AIC_tab) - 1)], isSingular))
  use.model <- names(which.min(AIC_tab$AIC[,-1][which(not.singulars)]))
  mod <- AIC_tab[[paste0("mod_", use.model)]]
  
  cat("Best model is: ", use.model)
  
  
  unique_streams <- unique(mstats[,c("Creek", "Basin", "year")])
  unique_streams <- na.omit(cbind.data.frame(unique_streams, Relative.percent.RBT = 0))
  unique_streams$Creek <- as.character(unique_streams$Creek)
  unique_streams$Basin <- as.character(unique_streams$Basin)
  
  pred <- predict(mod, newdata = unique_streams)
  pred <- cbind(unique_streams, pred)
  colnames(pred)[ncol(pred)] <- paste0("pred_", vars)
  rownames(pred) <- 1:nrow(pred)
  pred$Relative.percent.RBT <- NULL
  pred <- merge(pred, AIC_tab$dat[,c("Basin", "Creek", "year", "Relative.percent.RBT", vars)])
  
  return(list(pred = pred, mod = mod, used.model = use.model))
}
