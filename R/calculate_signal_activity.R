#' @title getSiganlActivity
#' @description calculate the upstream signals activity in the multilay signal network of cell communication
#'
#' @param ExprMat Matrix, The expression matrix with rows as genes (gene symbols) and columns as cells.
#' @param DistMat Matrix, The distance matrix represents the distance between two cells (e.g., Euclidean distance).
#' @param AnnoMat Matrix, The annotation matrix for cell type annotation. The first column is barcode and the second is cell type.
#' @param MulNetList List, The sublist is the multilayer signal network between of specific celltype pair.
#' @param Receiver Character, cell type that sending signals.
#' @param Sender Vector, cell types that receiving signals.
#' @param ProjectName Character, The project name of running jobs for now. Generate a working directory to save the final result.
#'
#' @return List, The sublist contains the activity of upstream signal pairs and expression of downstream target genes in specific celltype pair.
#' @export
#'
#' @import dplyr
#'
getSiganlActivity <- function(ExprMat, DistMat, AnnoMat, MulNetList,
                              Receiver, Sender = NULL, ProjectName = NULL){

  ## library

  # loadNamespace('dplyr')
  # loadNamespace('Seurat')
  # loadNamespace('SeuratWrappers')

  ## work directory

  if(is.null(ProjectName)){

    ProjectName <- format(Sys.time(),format="%Y-%m-%d_%H-%M-%S")

  }
  WorkDir <- paste0("./runModel/work_",ProjectName)
  dir.create(WorkDir, recursive = TRUE,showWarnings = F)
  cat(paste0("WorkDir: ",WorkDir,'\n'))

  ## main

  MulNetList <- MulNetList[grep(paste0("_",Receiver),names(MulNetList))]
  cellpair <- names(MulNetList)

  signalActivity <- list()
  if(is.null(Sender)){


    MulNetTab <- lapply(MulNetList, function(mlnet){

      ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target)
      rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
      tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)

      res = ligrec %>% merge(., rectf, by = 'Receptor') %>%
        merge(., tftg, by = 'TF') %>%
        dplyr::select(Ligand, Receptor, TF, Target) %>%
        arrange(Ligand, Receptor)

    }) %>% do.call('rbind',.)
    LRpairs <- by(MulNetTab, as.character(MulNetTab$Target), function(x){paste(x$Ligand, x$Receptor, sep = "_")})
    LRpairs <- lapply(LRpairs, function(lrtg){lrtg[!duplicated(lrtg)]})
    LRpairs <- LRpairs[lengths(LRpairs)>1]
    TGs <- names(LRpairs)

    cat(paste0("calculate the regulatory score of LR pairs from microenvironment to ",Receiver,'\n'))
    cp = paste('TME',Receiver,sep = "_")
    signalActivity[[cp]] <- getCPSiganlActivity(ExprMat, DistMat, AnnoMat,
                                                LRpairs, TGs, Receiver, Sender)

  }else{

    cellpair <- intersect(cellpair,paste(Sender,Receiver,sep = "_"))
    if(length(cellpair)==0) return(signalActivity)

    for (cp in cellpair) {

      receiver <- gsub('.*_','',cp)
      sender <- gsub('_.*','',cp)

      MulNet <- MulNetList[[cp]]
      LRpairs <- getSiganlLinks(MulNet)
      LRpairs <- LRpairs[lengths(LRpairs)>1]
      TGs <- names(LRpairs)

      cat(paste0("calculate the regulatory score of LR pairs from ",sender,' to ',receiver,'\n'))
      signalActivity[[cp]] = getCPSiganlActivity(ExprMat, DistMat, AnnoMat,
                                                 LRpairs, TGs, receiver, sender)

    }

  }

  ## save

  for (cp in names(signalActivity)) {

    cpSignalActivity <- signalActivity[[cp]]
    if(length(cpSignalActivity[[1]])>0){
      saveRDS(cpSignalActivity, file = paste0(WorkDir,"/LRTG_allscore_",cp,".rds"))
    }


  }

  return(signalActivity)
}

#' @title getCPSiganlActivity
#' @description calculate the upstream signals activity in the multilay signal network of the specific celltype pair
#'
#' @param exprMat Matrix, The expression matrix with rows as genes (gene symbols) and columns as cells.
#' @param distMat Matrix, The distance matrix represents the distance between two cells (e.g., Euclidean distance).
#' @param annoMat Matrix, The annotation matrix for cell type annotation. The first column is barcode and the second is cell type.
#' @param LRpairs List, upstream signal pairs (ligand/receptor) in the multilay signal network between sender and receiver.
#' @param TGs Vector, downstream target genes in the multilay signal network between sender and receiver.
#' @param receiver Character, cell type that sending signals.
#' @param sender Character, cell type that receiving signals.
#'
#' @return List. The sublist contains the activity of upstream signal pairs (LRs_score) and expression of downstream target genes (TGs_expr) between sender and receiver.
#'
#' @import dplyr stringr
#'
getCPSiganlActivity <- function(exprMat, distMat, annoMat,
                                LRpairs, TGs, receiver, sender = NULL)
{

  ## library

  # loadNamespace('dplyr')

  ## obtain cell ID
  receBars = annoMat %>% dplyr::filter(Cluster == receiver) %>%
    dplyr::select(Barcode) %>% unlist() %>% as.character()
  if(is.character(sender)){
    sendBars = annoMat %>% dplyr::filter(Cluster == sender) %>%
      dplyr::select(Barcode) %>% unlist() %>% as.character()
  }else{
    sendBars = annoMat %>% dplyr::filter(Cluster != receiver) %>%
      dplyr::select(Barcode) %>% unlist() %>% as.character()
  }

  ## obtain ligands and receptors for each downstream target
  Receptors = lapply(LRpairs, function(lr){stringr::str_split(lr,"_",simplify = T)[,2]})
  Ligands = lapply(LRpairs, function(lr){stringr::str_split(lr,"_",simplify = T)[,1]})

  ## obtain ligands expression for each downstream target
  LigMats = lapply(TGs, function(tg){
    # print(tg)
    ligands = Ligands[[tg]]
    if(length(ligands)==1){
      lig_count = exprMat[ligands, sendBars]
      lig_count = matrix(lig_count,nrow = 1)
    }else{
      lig_count = exprMat[ligands, sendBars] %>% as.matrix()
    }
    rownames(lig_count) = LRpairs[[tg]]
    lig_count
  })
  names(LigMats) = TGs

  ## obtain receptors expression for each downstream target
  RecMats = lapply(TGs, function(tg){
    receptors = Receptors[[tg]]
    if(length(receptors)==1){
      rec_count = exprMat[receptors, receBars]
      rec_count = matrix(rec_count,nrow = 1)
    }else{
      rec_count = exprMat[receptors, receBars] %>% as.matrix()
    }
    rownames(rec_count) = LRpairs[[tg]]
    rec_count
  })
  names(RecMats) = TGs

  ## distance matrix between sender and receiver
  distMat = distMat[sendBars, receBars]
  distMatReci = 1/distMat

  ## signal (LRpairs) activity
  # t1 <- Sys.time(); message(paste0('Start at: ',as.character(t1)))
  LRs_score = lapply(TGs, function(tg){

    # print(tg)
    LigMat = LigMats[[tg]]
    RecMat = RecMats[[tg]]
    lr = LRpairs[[tg]]

    LR_score = RecMat*(LigMat%*%distMatReci)
    LR_score = t(LR_score)
    colnames(LR_score) = lr
    rownames(LR_score) = receBars
    LR_score

  })
  names(LRs_score) = TGs
  # t2 <- Sys.time(); message(paste0('End at: ',as.character(t2)))
  # t2-t1

  ## Target expression
  TGs_expr = lapply(TGs, function(tg){ exprMat[tg, receBars] })
  names(TGs_expr) = TGs

  ## output
  result = list(LRs_score = LRs_score, TGs_expr = TGs_expr)

  return(result)

}

## 获取指定TG上游LRpair
#' @title getSiganlLinks
#' @description get the upstream signal pairs of downstream target genes in the multilayer signal network
#'
#' @param mulNet List, the multilayer signal network between of specific celltype pair.
#'
#' @return List. The sublist contains the upstream signal pairs of specific target genes.
#' @import dplyr
getSiganlLinks <- function(mulNet)
{

  ## library

  # loadNamespace('dplyr')

  ## main

  ligrec = data.frame(Ligand = mulNet$LigRec$source, Receptor = mulNet$LigRec$target)
  rectf = data.frame(Receptor = mulNet$RecTF$source, TF = mulNet$RecTF$target)
  tftg = data.frame(TF = mulNet$TFTar$source, Target = mulNet$TFTar$target)

  mulNet_tab = ligrec %>% merge(., rectf, by = 'Receptor') %>%
    merge(., tftg, by = 'TF') %>%
    dplyr::select(Ligand, Receptor, TF, Target) %>%
    arrange(Ligand, Receptor)

  LRTG_link = by(mulNet_tab, as.character(mulNet_tab$Target), function(x){paste(x$Ligand, x$Receptor, sep = "_")})
  LRTG_link = lapply(LRTG_link, function(lrtg){lrtg[!duplicated(lrtg)]})

  return(LRTG_link)
}

#' @title getCPSiganlActivity
#' @description get The imputation of expression matrix using the ALRA method in Seurat
#'
#' @param exprMat Matrix, The expression matrix with rows as genes (gene symbols) and columns as cells.
#' @param use.seed Logical, Whether to set a seed.
#' @param seed Only relevant if use.seed=T, set the seed for the random process.
#'
#' @return Matrix, the imputed expression matrix.
#' @export
#'
#' @import dplyr Seurat
#' @importFrom SeuratWrappers RunALRA
runImputation <- function(exprMat, use.seed = TRUE, seed = 4321)
{

  ## library

  # loadNamespace('dplyr')
  # loadNamespace('Seurat')
  # loadNamespace('SeuratWrappers')

  ## main

  expr.Impute <- CreateSeuratObject(exprMat,verbose=F)
  cat('Using imputation method ALRA wrapped in Seurat\n')
  if(use.seed) set.seed(seed)
  expr.Impute <- suppressMessages(RunALRA(expr.Impute))
  exprMat.Impute <- expr.Impute@assays$alra@data

  return(exprMat.Impute)

}

