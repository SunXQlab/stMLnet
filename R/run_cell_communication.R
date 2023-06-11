#' Title
#'
#' @param ExprMat Matrix, The expression matrix with rows as genes (gene symbols) and columns as cells.
#' @param DistMat Matrix, The distance matrix represents the distance between two cells (e.g., Euclidean distance).
#' @param AnnoMat Matrix, The annotation matrix for cell type annotation. The first column is barcode and the second is cell type.
#' @param LigClus Vector, cell types that sending signals.
#' @param RecClus Vector, cell types that receiving signals.
#' @param Normalize Logical, Whether to do normalization.
#' @param NormMethod Character, Denotes which normalization method to use. Available options are: LogNormalze, SCTransform. See Seurat package for more details.
#' @param logfc.ct Numercial, Screening threshold for FindMarkers in Seurat. The default setting is 0.1.
#' @param pct.ct Numercial, Screening threshold for FindMarkers in Seurat. The default setting is 0.05.
#' @param pval.ct Numercial, Screening threshold for FindMarkers in Seurat. The default setting is 0.05.
#' @param expr.ct Numercial, Screening threshold for high expressed gene in groups of cells. Default is 0.05.
#' @param OutputDir Character, The output path of the currently running job where temporary and final results will be saved.
#' @param Databases List, The prior database used by running jobs for now. Databases includes Ligand-Receptor interactions (LigRec.DB), Receptor-TF interactions (RecTF.DB) and TF-Target interactions (TFTG.DB).
#' @param SigMethod Character, Denotes the strategy for filtering downstream pairing signals (Receptor-TF, TF-Target).  Available options are: Fisher(default, meaning Fisher exact test) and Search (meaning searching in database).
#' @param TGList List, The target genes of interest in groups of cells (RecClus).
#' @param LigList List, The potential ligands in groups of cells (LigClus).
#' @param RecList List, The potential receptors in groups of cells (RecClus).
#' @param NCores Numercial, set the cores for the parallel process.
#' @param AutoPara Logical, Whether to do automatic optimization parameter.
#' @param NTrees Numercial, number of trees in random forests model, see Seurat package for more details.
#' @param NTrys Numercial, Only relevant if AutoPara=F. Number of variables to possibly split at in each node, see Seurat package for more details.
#' @param TreeMethod Character, Only relevant if AutoPara=F. Splitting rule, see Seurat package for more details.
#' @param NodeSize Numercial, Only relevant if AutoPara=F. Minimal node size, see Seurat package for more details.
#' @param NPert Numercial, Number of shuffle.
#'
#' @return List, containning the output folder of the multilayer signal network (MLnetDir), the signal activity (ActivityDir) and the signal importance (ImportDir).
#' @export
#'
runstMLnet <- function(
    ExprMat, AnnoMat, DistMat, LigClus = NULL, RecClus = NULL,
    Normalize = F, NormMethod = 'LogNormalze',
    logfc.ct = 0.1, pct.ct = 0.05, pval.ct = 0.05, expr.ct = 0.1,
    OutputDir = NULL, Databases = NULL, SigMethod = 'Fisher',
    TGList = NULL, LigList = NULL, RecList = NULL,
    NCores = 6, AutoPara = TRUE, NTrees = 500, NTrys = 10,
    TreeMethod = 'variance', NodeSize = 5,  NPert = 10
){

  ## Step0 prepare

  clusters <- ex_inputs$annoMat$Cluster %>% unique() %>% as.character()
  exprMat.Impute <- runImputation(exprMat = ex_inputs$exprMat)

  if(is.null(OutputDir)){

    OutputDir <- format(Sys.time(),format="%Y-%m-%d_%H-%M-%S")
    OutputDir <- paste0(getwd(),'/stMLnet_work_',OutputDir,'/')

  }else{

    OutputDir <- OutputDir

  }

  ## Step1 create mulityayer network

  resMLnet <- runMLnet(ExprMat, AnnoMat, LigClus = LigClus, RecClus = RecClus,
                       Normalize = Normalize, NormMethod = NormMethod,
                       OutputDir = OutputDir, Databases = Databases, SigMethod = SigMethod,
                       logfc.ct = logfc.ct, pct.ct = pct.ct, pval.ct = pval.ct, expr.ct = expr.ct,
                       TGList=TGList, LigList=LigList, RecList=RecList)

  ## Step2 calculate Signal Activity

  ex_mulnetlist <- list()
  for (i in 1:length(resMLnet$mlnets)) {

    mlnets <- resMLnet$mlnets[[i]]
    for (j in 1:length(mlnets)) {
      mlnet <- mlnets[[j]]
      if(nrow(mlnet$LigRec)!=0) ex_mulnetlist[[names(mlnets)[j]]] = mlnet
    }

  }

  resSigActList <- list()
  for (cluster in clusters) {

    Sender <- clusters[clusters!=cluster]
    resSigActList[[cluster]] <- getSiganlActivity(ExprMat = exprMat.Impute,
                                                  DistMat = DistMat,
                                                  AnnoMat = AnnoMat,
                                                  MulNetList = ex_mulnetlist,
                                                  Receiver = cluster,
                                                  Sender = Sender,
                                                  OutputDir = OutputDir)

  }

  resSigActList <- do.call(c,resSigActList)
  names(resSigActList) <- gsub('.*\\.','',names(resSigActList))
  for (cp in names(resSigActList)) {

    cpSignalActivity <- resSigActList[[cp]]
    if(length(cpSignalActivity[[1]])==0){
      resSigActList[[cp]] <- NULL
    }

  }

  ## Step3 calculate signal importance

  time_ls <- c()
  resSigImportList <- list()
  for(cp in names(resSigActList)){

    label <- cp
    LRTG_allscore <- resSigActList[[cp]]
    message(paste0('running jobs: ',label))

    t1 <- Sys.time()
    resSigImportList[[cp]] <- getSiganlImport(
      SiganlActivity = LRTG_allscore, Lable = label, OutputDir = OutputDir,
      NCores = NCores, AutoPara = AutoPara, NTrees = NTrees, NTrys = NTrys,
      TreeMethod = TreeMethod, NodeSize = NodeSize,  NPert = NPert
    )
    t2 <- Sys.time()
    time_ls <- c(time_ls,paste(signif(t2-t1,4),units(signif(t2-t1,4)),sep = ' '))

  }

  ## Step4

  dirpath <- list(
    MLnetDir = paste0(OutputDir,'/runscMLnet/'),
    ActivityDir = paste0(OutputDir,'/runModel/'),
    ImportDir = paste0(OutputDir,'/getPIM/')
  )

  cat(paste0('all output are saved to the current directory\n'))
  cat(paste0('Detail of multilayer signal network are saved to the current directory: ',dirpath$MLnetDir,'\n'))
  cat(paste0('Detail of signals activity are saved to the current directory: ',dirpath$ActivityDir,'/\n'))
  cat(paste0('Detail of signals importance are saved to the current directory: ',dirpath$ImportDir,'/\n'))

  return(dirpath)

}
