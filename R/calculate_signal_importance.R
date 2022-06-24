#' @title getSiganlImport
#' @description calculate the upstream signal pairs or signals importance in the multilay signal network of cell communication
#'
#' @param SiganlActivity List, The sublist contains the activity of upstream signal pairs and expression of downstream target genes in specific celltype pair.
#' @param Lable Character, Denotes which celltype pair to study.
#' @param ProjectName Character, The project name of running jobs for now. Generate a working directory to save the final result.
#' @param NCores Numercial, set the cores for the parallel process.
#' @param AutoPara Logical, Whether to do automatic optimization parameter.
#' @param NTrees Numercial, number of trees in random forests model, see Seurat package for more details.
#' @param NTrys Numercial, Only relevant if AutoPara=F. Number of variables to possibly split at in each node, see Seurat package for more details.
#' @param TreeMethod Character, Only relevant if AutoPara=F. Splitting rule, see Seurat package for more details.
#' @param NodeSize Numercial, Only relevant if AutoPara=F. Minimal node size, see Seurat package for more details.
#' @param NPert Numercial, Number of shuffle.
#'
#' @return  List, The first sublist (df_im) is the importance of upstream signal pairs (Ligand/Receptor pairs), the second sublist (df_pim) is the importance of upstream signals (ligands/receptors).
#' @export
#' @import dplyr ranger caret doParallel doSNOW foreach
getSiganlImport <- function(SiganlActivity, Lable, ProjectName = NULL, NCores = NULL,
                            AutoPara = TRUE, NTrees = 500, NTrys = 10, TreeMethod = 'variance',
                            NodeSize = 5,  NPert = 10){

  ## library

  # loadNamespace("dplyr")
  # loadNamespace("ranger")
  # loadNamespace("caret")
  # loadNamespace("Metrics")
  # loadNamespace("doParallel")
  # loadNamespace("parallel")
  # loadNamespace("doSNOW")

  ## parallel

  N.cores <- parallel::detectCores()
  if(is.null(NCores)|NCores>=N.cores){

    NCores <- ceiling(N.cores/2)

  }

  ## main

  LRTG_allscore <- SiganlActivity
  n.TG <- length(LRTG_allscore$LRs_score)

  if(NCores == 1){

    rfModelList <- lapply(1:n.TG, function(i){

      trainx = LRTG_allscore$LRs_score[[i]]
      trainy = LRTG_allscore$TGs_expr[[i]]
      runRFModel(trainx = trainx, trainy = trainy, auto_para = AutoPara,
                 n.trees = NTrees, n.trys = NTrys,tree.method = TreeMethod,
                 node.size = NodeSize,  nPert = NPert)

    })
    names(rfModelList) <- names(LRTG_allscore$LRs_score)

  }else{

    # t1 <- Sys.time()
    # message(paste0('Start at ',as.character(t1)))
    cl <- snow::makeSOCKcluster(NCores)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(min=0, max=n.TG, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    rfModelList <- foreach::foreach(i=1:n.TG, .packages=c("dplyr","ranger",'caret'),
                           .options.snow=opts, .errorhandling = "stop"
    ) %dopar% {
      message(paste0("\ngene",i,'\n'))
      trainx = LRTG_allscore$LRs_score[[i]]
      trainy = LRTG_allscore$TGs_expr[[i]]
      runRFModel(trainx = trainx, trainy = trainy, auto_para = AutoPara,
                 n.trees = NTrees, n.trys = NTrys,tree.method = TreeMethod,
                 node.size = NodeSize,  nPert = NPert)
    }
    names(rfModelList) <- names(LRTG_allscore$LRs_score)
    close(pb)
    stopCluster(cl)
    gc()
    # t2 <- Sys.time()
    # message(paste0('End at ',as.character(t2)))

  }

  ## output

  result <- mergeVarsImport(rfModelList = rfModelList, label = Lable, projectName = ProjectName)

  return(result)

}

#' @title runRFModel
#' @description train a random forest model for specific target genes to calculate the importance of upstream signal pairs or signals
#'
#' @param trainx Matrix, the activity of upstream signal pairs.
#' @param trainy Vector, the expression of downstream target genes.
#' @param auto_para Logical, Whether to do automatic optimization parameter.
#' @param n.trees Numercial, number of trees in random forests model, see Seurat package for more details.
#' @param n.trys Numercial, Only relevant if AutoPara=F. Number of variables to possibly split at in each node, see Seurat package for more details.
#' @param tree.method Character, Only relevant if AutoPara=F. Splitting rule, see Seurat package for more details.
#' @param node.size Numercial, Only relevant if AutoPara=F. Minimal node size, see Seurat package for more details.
#' @param nPert Numercial, Number of shuffle.
#'
#' @return List, The first sublist (model) is the trained model of specific target genes. The second sublist (df_IM) is the importance of upstream signal pairs corresponding the specific target genes and the third sublist (df_pIM) is the importance of upstream signals corresponding the specific target genes.
#' @export
#' @import dplyr ranger caret
runRFModel = function(trainx, trainy, auto_para = TRUE, n.trees = 500,
                      n.trys = 10, tree.method = 'variance', node.size = 5,  nPert = 10){

  ## library

  # loadNamespace("dplyr")
  # loadNamespace("ranger")
  # loadNamespace("caret")
  # loadNamespace("Metrics")

  ## train set

  LRpairs = colnames(trainx)
  data = as.data.frame(cbind(trainx, trainy))
  colnames(data) = c(paste0('LR',seq(ncol(trainx))), "Target")

  ## obtain partial variables (ligands/receptors in LR pairs)

  LRTab = data.frame(LRpair = LRpairs)
  LRTab$Ligand = strsplit(LRpairs,"_") %>% do.call('rbind',.) %>% .[,1]
  LRTab$Receptor = strsplit(LRpairs,"_") %>% do.call('rbind',.) %>% .[,2]

  ## filter sample

  cat("\nRemove zero Target")
  zeroTarget = which(trainy==0)
  if(length(zeroTarget)>1) data = data[-zeroTarget,]

  ## filter variables

  cat("\nRemove zero LR")
  zeroLR = which(colSums(data[,-ncol(data)]==0)==nrow(data))
  zeroLR = c(zeroLR,which(is.na(colSums(data[,-ncol(data)]==0))))
  if(length(zeroLR)>1){
    data = data[,-zeroLR]
    LRTab = LRTab[-zeroLR,]
  }

  ## Scale

  cat("\nScale data")
  data = scale(data)
  data = as.data.frame(data)

  ## training parameter

  t1 <- Sys.time()
  cat("\nParameter Tuning")
  cat(paste0("\nStart at: ",as.character(t1)))
  if(auto_para == TRUE){

    fitControl <- trainControl(
      method = "cv",
      number = 5,
      search = 'random',
      allowParallel = FALSE)

    set.seed(2021)
    rfFit <- caret::train(Target ~ ., data = data,
                          method = 'ranger',
                          trControl = fitControl,
                          tuneLength = 20)

    sel_mtry = as.numeric(rfFit$bestTune$mtry)
    sel_splitrule = as.character(rfFit$bestTune$splitrule)
    sel_min.node.size = as.numeric(rfFit$bestTune$min.node.size)
    sel_num.trees = n.trees

  }else{

    sel_mtry = n.trys
    sel_splitrule = tree.method
    sel_min.node.size = node.size
    sel_num.trees = n.trees

  }
  t2 <- Sys.time()
  cat(paste0("\nEnd at: ",as.character(t2)))
  cat(paste0('\nAbout ',signif(t2-t1,digits = 4),' ',units(t2-t1)))

  ## training model

  cat("\nTrain final model")
  cat(paste0("\nThe final parameters used for the model: mtry = ",sel_mtry,
             ', splitrule = ',sel_splitrule, ', min.node.size = ',sel_min.node.size,
             ' ,num.trees = ',sel_num.trees))
  finalFit <- ranger::ranger(formula = Target ~ ., data = data, num.trees = sel_num.trees, seed = 2021,
                     splitrule = sel_splitrule, mtry = sel_mtry, min.node.size = sel_min.node.size,
                     importance = 'permutation', keep.inbag = TRUE, oob.error = TRUE, num.threads = 1)

  ## get variables importance

  cat("\nobtain variables importance")
  df_IM = finalFit$variable.importance
  df_IM = data.frame(IM = df_IM)
  df_IM$LRpair = LRTab$LRpair

  ## get partial variables importance

  t1 <- Sys.time()
  cat("\nshuffle Ligand or Receptor")
  cat(paste0("\nStart at: ",as.character(t1)))
  df_pIM <- lapply(seq(nPert), function(j){

    forest = finalFit
    data_X=data[,-ncol(data)]
    data_y=data$Target
    getPartImport(forest, data_X, data_y, LRTab)

  })
  df_pIM <- Reduce("+", df_pIM)/nPert
  t2 <- Sys.time() # 5.74 mins
  cat(paste0("\nEnd at: ",as.character(t2)))
  cat(paste0('\nAbout ',signif(t2-t1,digits = 4),' ',units(t2-t1)))

  result <- list(model = finalFit,
                 df_IM = df_IM,
                 df_pIM = df_pIM)

  return(result)

}

#' @title mergeVarsImport
#' @description merge the feature importance of different target genes
#'
#' @param rfModelList List, the trained model and variable importance of specific target genes in the multilay signal network of cell communication
#' @param label Character, Denotes which celltype pair to study.
#' @param projectName Character, The project name of running jobs for now. Generate a working directory to save the final result.
#'
#' @return List, The first sublist (df_im) is the importance of upstream signal pairs (Ligand/Receptor pairs), the second sublist (df_pim) is the importance of upstream signals (ligands/receptors).
#' @export
#' @import dplyr
mergeVarsImport <- function(rfModelList, label, projectName){

  ## library

  # loadNamespace("dplyr")

  ## workdir

  WorkDir <- paste0("./getPIM/work_",projectName)
  dir.create(WorkDir, recursive = TRUE,showWarnings = F)
  cat(paste0("WorkDir: ",WorkDir,'\n'))

  ## variable Importance

  df_im = lapply(1:length(rfModelList), function(i){

    rfModel = rfModelList[[i]]
    im = rfModel$df_IM

    if(is.null(im)|sum(im$IM,na.rm = T)==0){

      im = data.frame()

    }else{

      im$Ligand = stringr::str_split(im$LRpair,"_",simplify = T)[,1]
      im$Receptor = stringr::str_split(im$LRpair,"_",simplify = T)[,2]
      im$Target = names(res_ls)[i]
      im$im_norm = im$IM/sum(im$IM)

    }

    im

  })
  df_im = do.call("rbind",df_im)
  df_im = df_im[,c(2:5,1,6)]
  saveRDS(df_im, paste0(WorkDir,'/LRTG_im_clean_',label,'.rds'))

  ## partial variable Importance

  df_pim = lapply(1:length(rfModelList), function(i){

    rfModel = rfModelList[[i]]

    df_LR = data.frame(LRpair = rfModel$df_IM$LRpair)
    df_LR$Ligand = stringr::str_split(df_LR$LRpair,"_",simplify = T)[,1]
    df_LR$Receptor = stringr::str_split(df_LR$LRpair,"_",simplify = T)[,2]

    if(is.null(res$df_pIM)|sum(res$df_pIM$pIM,na.rm = T)==0){

      pim = data.frame()

    }else{

      pim = data.frame(regulator = gsub("shuffle_","",rownames(res$df_pIM)))
      pim$Target = names(res_ls)[i]
      pim$pIM = res$df_pIM$pIM
      pim$type = pim$regulator %in% df_LR$Ligand
      pim$type[pim$type==TRUE] = 'Ligand'
      pim$type[pim$type==FALSE] = 'Receptor'

    }
    pim

  })
  df_pim = do.call('rbind',df_pim)
  saveRDS(df_pim, paste0(WorkDir,'/LRTG_pim_clean_',label,'.rds'))

  ## output

  result <- list(df_im = df_im, df_pim = df_pim)

  return(result)

}

#' @title getPartImport
#' @description get importance of partial variable (signals) by shuffle in OOB data
#'
#' @param rfmodel an ranger object, the trained model of specific target genes.
#' @param data_X Matrix, the activity of upstream signal pairs.
#' @param data_y Vector, the expression of downstream target genes.
#' @param LRTab Matrix, the upstream signal pairs.
#'
#' @return Matrix, the importance of partial variable (signals) corresponding the specific target genes
#' @import dplyr ranger
getPartImport = function(rfmodel, data_X, data_y, LRTab)
{

  ## library

  # loadNamespace("dplyr")
  # loadNamespace("ranger")

  ## get pred value based on OOB data and trained model

  oob_pred = getOOBPred(rfmodel, data_X)
  oob_pred = as.data.frame(oob_pred)

  ## evaluate the trained model based on pred and true value

  # pred value come from the OOB data before shuffle
  # oob_mse_bef: vector (500)
  # MSE of each tree in trained model before shuffle
  oob_mse_bef = lapply(seq(ncol(oob_pred)), function(k){

    oob_pred_k = oob_pred[,k]
    mean((oob_pred_k-data_y)^2,na.rm = T)

  })
  oob_mse_bef = unlist(oob_mse_bef)

  ## get partial variables

  partVars = unique(c(LRTab$Ligand,LRTab$Receptor))

  ## evaluate the trained model based on pred and true value

  # pred value come from the OOB data after shuffle
  # oob_mse_aft: matrix (ix500)
  # oob_mse_aft_i: MSE of each tree in trained model after shuffle i-th partial variables
  oob_mse_aft = lapply(seq(length(partVars)), function(i){

    # shuffle partial variable (i) to produce new trainX
    newX = shufflePartVars(data_X,i,LRTab)

    # pred value based on OOB data after shuffle one Lig/Rec
    oob_pred_i = getOOBPred(rfmodel, newX)
    oob_pred_i = as.data.frame(oob_pred_i)

    # oob_mse_aft_i: vector (500)
    # MSE of each tree in trained model after shuffle
    oob_mse_aft_i = lapply(seq(ncol(oob_pred_i)), function(k){

      oob_pred_i_k = oob_pred_i[,k]
      mean((oob_pred_i_k-data_y)^2,na.rm = T)

    })
    oob_mse_aft_i = unlist(oob_mse_aft_i)

  })
  oob_mse_aft = do.call("cbind",oob_mse_aft)

  ## get mean MSE of each partial variables

  df_pIM = apply(oob_mse_aft, 2, function(oob_mse_aft_i){ mean(oob_mse_aft_i-oob_mse_bef) })
  df_pIM = data.frame(pIM = df_pIM)
  rownames(df_pIM) = paste0('shuffle_',partVars)

  return(df_pIM)

}

#' @title getOOBPred
#' @description get predictions for each tree in trained model based on OOB data
#'
#' @param rfmodel an ranger object, the trained model of specific target genes.
#' @param data_X Matrix, the activity of upstream signal pairs.
#'
#' @return Matrix, the predicted activity of upstream signal pairs.
#' @import dplyr ranger
#' @importFrom stats predict
getOOBPred <- function(rfmodel, data_X)
{

  ## library

  # loadNamespace("dplyr")
  # loadNamespace("ranger")

  ## main

  preds <- stats::predict(rfmodel, data_X, predict.all=TRUE)
  oob <- rfmodel$inbag.counts
  oob <- do.call("cbind",oob)
  oob <- oob==0
  oob[which(!oob)] = NA
  preds.oob = oob*preds$predictions
  return (preds.oob)

}

#' @title shufflePartVars
#' @description shuffle the activity of upstream signal pairs that involve i-th partial variables
#'
#' @param data_X Matrix, the activity of upstream signal pairs.
#' @param i i-th partial variables to be shuffled.
#' @param LRTab Matrix, the upstream signal pairs.
#'
#' @return Matrix, the activity of upstream signal pairs after shuffling
#' @import dplyr
shufflePartVars <- function(data_X,i,LRTab)
{

  ## library

  # loadNamespace("dplyr")

  ## main

  Vars = unique(c(LRTab$Ligand,LRTab$Receptor))
  var <- Vars[i]
  var_ids <- union(which(LRTab$Ligand %in% var),which(LRTab$Receptor %in% var))
  for(id in var_ids) {
    data_X[,id] <- sample(x = unlist(data_X[,id]), size = nrow(data_X), replace = FALSE)
  }
  return(data_X)

}
