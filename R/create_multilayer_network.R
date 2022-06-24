#' @title runMLnet
#' @description Creates an object that stores the multilayer signal network between LigClus and RecClus
#'
#' @param ExprMat Matrix, The expression matrix with rows as genes (gene symbols) and columns as cells.
#' @param AnnoMat Matrix, The annotation matrix for cell type annotation. The first column is barcode and the second is cell type.
#' @param LigClus Vector, cell types that sending signals.
#' @param RecClus Vector, cell types that receiving signals.
#' @param Normalize Logical, Whether to do normalization.
#' @param NormMethod Character, Denotes which normalization method to use. Available options are: LogNormalze, SCTransform. See Seurat package for more details.
#' @param logfc.ct Numercial, Screening threshold for FindMarkers in Seurat. The default setting is 0.1.
#' @param pct.ct Numercial, Screening threshold for FindMarkers in Seurat. The default setting is 0.05.
#' @param pval.ct Numercial, Screening threshold for FindMarkers in Seurat. The default setting is 0.05.
#' @param expr.ct Numercial, Screening threshold for high expressed gene in groups of cells. Default is 0.05.
#' @param ProjectName Character, The project name of running jobs for now. Generate a working directory to save the final result.
#' @param Databases List, The prior database used by running jobs for now. Databases includes Ligand-Receptor interactions (LigRec.DB), Receptor-TF interactions (RecTF.DB) and TF-Target interactions (TFTG.DB).
#' @param TGList List, The target genes of interest in groups of cells (RecClus).
#' @param LigList List, The potential ligands in groups of cells (LigClus).
#' @param RecList List, The potential receptors in groups of cells (RecClus).
#'
#' @return List. The first sublist (mlnets) is the multilayer signal network between LigClus and RecClu, the second sublist (detail) is the detail of cell communication between LigClus and RecClu.
#' @export
#'
#' @import Seurat dplyr
#' @importFrom stats quantile
#'
runMLnet <- function(ExprMat, AnnoMat, LigClus = NULL, RecClus = NULL,
                     Normalize = T, NormMethod = "LogNormalize",
                     logfc.ct = 0.1, pct.ct = 0.05, pval.ct = 0.05, expr.ct = 0.1,
                     ProjectName = NULL, Databases = NULL,
                     TGList=NULL, LigList=NULL, RecList=NULL){

  ## library

  # loadNamespace('Seurat')
  # loadNamespace('dplyr')

  ## databases

  if(is.null(Databases)){

    cat("load default database\n")
    quan.cutoff <- 0.98
    Databases <- ex_databases
    Databases$RecTF.DB <- Databases$RecTF.DB %>%
      .[.$score > quantile(.$score, quan.cutoff),] %>%
      dplyr::distinct(source, target)
    Databases$LigRec.DB <- Databases$LigRec.DB %>%
      dplyr::distinct(source, target) %>%
      dplyr::filter(target %in% Databases$RecTF.DB$source)
    Databases$TFTG.DB <- Databases$TFTG.DB %>%
      dplyr::distinct(source, target) %>%
      dplyr::filter(source %in% Databases$RecTF.DB$target)

  }else{

    cat("load user database\n")
    Databases$RecTF.DB <- Databases$RecTF.DB %>%
      dplyr::distinct(source, target)
    Databases$LigRec.DB <- Databases$LigRec.DB %>%
      dplyr::distinct(source, target) %>%
      dplyr::filter(target %in% Databases$RecTF.DB$source)
    Databases$TFTG.DB <- Databases$TFTG.DB %>%
      dplyr::distinct(source, target) %>%
      dplyr::filter(source %in% Databases$RecTF.DB$target)

  }

  ## work directory

  if(is.null(ProjectName)){

    ProjectName <- format(Sys.time(),format="%Y-%m-%d_%H-%M-%S")

  }
  WorkDir <- paste0("./runscMLnet/work_",ProjectName)
  dir.create(WorkDir, recursive = TRUE,showWarnings = F)
  cat(paste0("WorkDir: ",WorkDir,'\n'))

  ## check LigClu and RecClu

  if(is.null(LigClus)){

    LigClus <- unique(AnnoMat$Cluster) %>% as.character()

  }
  if(is.null(RecClus)){

    RecClus <- unique(AnnoMat$Cluster) %>% as.character()

  }

  ## objects

  inputs = list(
      parameters = list(
        LigClus = LigClus,
        RecClus = RecClus,
        Project = ProjectName,
        logfc.ct = logfc.ct,
        pct.ct = pct.ct,
        pval.ct = pval.ct,
        expr.ct = expr.ct
      ),
      data = list(
        df_data = ExprMat,
        df_anno = AnnoMat,
        ls_clusters = as.character(unique(AnnoMat$Cluster)),
        ls_targets = NA,
        ls_ligands = NA,
        ls_receptors = NA
      )
    )

  ## normalization

  if(Normalize){

    cat("perform normalization\n")
    inputs <- runNormalize(inputs, method = NormMethod)

  }else{

    cat("skip normalization\n")
    inputs$data$df_norm <- inputs$data$df_data

  }

  ## potential signals

  if(any(is.null(TGList),is.null(LigList))){

    inputs = getDiffExpGene(inputs)

  }else{

    inputs$data$df_degs = NA

  }
  if(is.null(TGList)){

    inputs$data$ls_targets <- list()
    for (Clu in inputs$parameters$RecClu) {

      df_degs_clu <- inputs$data$df_degs[[Clu]]
      inputs$data$ls_targets[[Clu]] <- rownames(df_degs_clu)[
        df_degs_clu$p_val_adj <= 0.05 & abs(df_degs_clu$avg_log2FC) >= 1
        # potential target filter: abs(avg_log2FC)>=1 & padj<=0.05
      ]

    }

  }else{

    inputs$data$ls_targets <- TGList

  }
  if(is.null(LigList)){

    inputs$data$ls_ligands <- list()
    for (Clu in inputs$parameters$LigClu) {

      df_degs_clu <- inputs$data$df_degs[[Clu]]
      liglist <- rownames(df_degs_clu)[df_degs_clu$p_val_adj <= 0.05 & df_degs_clu$avg_log2FC >= 0]
      liglist <- intersect(liglist,unique(Databases$LigRec.DB$source))
      inputs$data$ls_ligands[[Clu]] <- liglist

    }

  }else{

    inputs$data$ls_ligands <- LigList

  }
  if(is.null(RecList)){

    inputs = getHighExpGene(inputs)
    inputs$data$ls_receptors <- lapply(inputs$data$ls_receptors,function(ls_recs){

      intersect(ls_recs,unique(Databases$LigRec.DB$target))

    })

  }else{

    inputs$data$ls_receptors <- RecList

  }

  ## get multi-layer

  outputs <- list(mlnets = list(), details = list())
  for(RecClu in inputs$parameters$RecClus){

    ## parameter

    LigClus <- inputs$parameters$LigClus
    LigClus <- LigClus[!LigClus %in% RecClu]
    details <- matrix(ncol = length(LigClus), nrow = 10) %>% as.data.frame()
    mlnets <- list()

    for(i in 1:length(LigClus)){

      LigClu <- LigClus[i]
      cat(paste0(LigClu,"_",RecClu,'\n'))

      resMLnet <- getCellPairMLnet(inputs, LigClu, RecClu, Databases)
      mlnets[[i]] <- resMLnet$mlnet
      details[,i] <- resMLnet$detail


    }
    names(mlnets) <- paste(LigClus,RecClu,sep = "_")
    colnames(details) <- paste(LigClus,RecClu,sep = "_")
    rownames(details) <- c('Lig_bk','Rec_bk','target_bk',
                          "LRpair","RecTFpair","TFTGpair",
                          "Ligand", "Receptor", "TF", "Target")
    write.csv(details, file = paste0(WorkDir,"/TME_",RecClu,".csv"))

    outputs$mlnets[[RecClu]] <- mlnets
    outputs$details[[RecClu]] <- details

  }

  return(outputs)

}

#' @title runNormalize
#' @description Normalize the expression data using function in Seurat package.
#'
#' @param inputs List, the parameters and data inputted by current project.
#' @param method Character, Normalization method, available options are: LogNormalze, SCTransform. See Seurat package for more details.
#'
#' @return List, inputs object contains the normalized data.
#'
#' @import Seurat dplyr
#'
runNormalize <- function(inputs, method = c("LogNormalze",'SCTransform'))
{

  ## library

  # loadNamespace('Seurat')
  # loadNamespace('dplyr')

  ## check

  if (!(method %in% c("LogNormalize", "SCTransform")))
    stop("wrong normalization method!")

  ## mian

  df_data = inputs$data$df_data
  seurobj = CreateSeuratObject(counts = df_data)

  if(method == "SCTransform") {
    seurobj = SCTransform(seurobj, verbose = TRUE)
    df_norm = seurobj@assays$SCT@data
  }else {
    seurobj = NormalizeData(seurobj, normalization.method = method, verbose = TRUE)
    df_norm = seurobj@assays$RNA@data
  }

  inputs$data$df_norm <- df_norm

  return(inputs)
}

#' @title getDiffExpGene
#' @description obtain the differently expressed genes using FindMarkers function in Seurat package.
#'
#' @param inputs List, the parameters and data inputted by current project.
#'
#' @return List, inputs object contains the differently expressed genes.
#'
#' @import Seurat dplyr
#'
getDiffExpGene <- function(inputs)
{

  ## library

  # loadNamespace('Seurat')
  # loadNamespace('dplyr')

  ## parameters

  clusters <- unique(c(inputs$parameters$LigClus,inputs$parameters$RecClus))
  df_norm <- inputs$data$df_norm
  df_anno <- inputs$data$df_anno
  logfc.ct <- inputs$parameters$logfc.ct
  pct.ct <- inputs$parameters$pct.ct
  pval.ct <- inputs$parameters$pval.ct

  ## results

  inputs$data$df_degs <- list()
  assayobj = CreateAssayObject(data = df_norm)
  for (Clu in clusters) {

    ## get barcode
    BarListResult <- getBarList(Clu, df_anno)
    Clus.1 <- BarListResult[[1]]
    Clus.2 <- BarListResult[[2]]

    ## find DEGs(use all other cells for FindMarkers)
    DEGs <- FindMarkers(object = assayobj, cells.1 = Clus.1, cells.2 = Clus.2,
                        logfc.threshold = logfc.ct, min.pct = pct.ct, verbose = T)
    if(!is.null(DEGs)){
      DEGs <- DEGs %>% filter(p_val_adj <= pval.ct)
      inputs$data$df_degs[[Clu]] <- DEGs
    }else{
      inputs$data$df_degs[[Clu]] <- NA
    }

  }

  return(inputs)

}

#' @title getHighExpGene
#' @description obtain the highly expressed genes according the detection and mean expression.
#'
#' @param inputs List, the parameters and data inputted by current project.
#'
#' @return List, inputs object contains the highy expressed genes.
#'
#' @import dplyr
#' @importFrom Matrix rowMeans
#' @importFrom Matrix rowSums
#'
getHighExpGene <- function(inputs)
{

  ## library

  # loadNamespace('dplyr')

  ## parameters

  clusters <- inputs$data$ls_clusters
  df_norm <- inputs$data$df_norm
  df_anno <- inputs$data$df_anno
  pct.ct <- pct.ct
  expr.ct <- expr.ct

  ## pct & expr

  df_mean <- lapply(clusters, function(cluster){

    source_mean <- Matrix::rowMeans(df_norm[,getBarList(cluster,df_anno)[[1]]])
    names(source_mean) <- rownames(df_norm)
    source_mean

  }) %>% do.call('cbind',.) %>% as.data.frame()
  colnames(df_mean) <- clusters

  df_pct <- lapply(clusters, function(cluster){

    dat <- df_norm[,getBarList(cluster,df_anno)[[1]]]
    pct <- Matrix::rowSums(dat>0)/ncol(dat)
    names(pct) <- rownames(df_norm)
    pct

  }) %>% do.call('cbind',.) %>% as.data.frame()
  colnames(df_pct) <- clusters

  ## results

  inputs$data$ls_receptors <- list()
  for (cluster in inputs$parameters$RecClus) {

    inputs$data$ls_receptors[[cluster]] <-
      rownames(df_mean)[df_mean[,cluster] >= expr.ct & df_pct[,cluster] >= pct.ct]

  }

  return(inputs)

}


#' @title getCellPairMLnet
#' @description obtain the multilayer signal network of specific cell pair (LigClu-RecClu).
#'
#' @param inputs List, the parameters and data inputted by current project.
#' @param ligclu Character, cell types that sending signals.
#' @param recclu Character, cell types that receiving signals.
#' @param databases List, The prior database used by running jobs for now. Databases includes Ligand-Receptor interactions (LigRec.DB), Receptor-TF interactions (RecTF.DB) and TF-Target interactions (TFTG.DB).
#'
#' @return List, The first sublist (mlnet) is the multilayer signal network between ligclu and recclu, the second sublist (detail) is the detail of cell communication between ligclu and recclu.
#'
#' @import dplyr
#' @importFrom Matrix rowMeans
#' @importFrom stats fisher.test
#'
getCellPairMLnet <- function(inputs, ligclu, recclu, databases)
{

  ## library

  # loadNamespace('dplyr')

  ## parameter

  df_norm <- inputs$data$df_norm
  ls_ligands <- inputs$data$ls_ligands
  ls_receptors <- inputs$data$ls_receptors
  ls_targets <- inputs$data$ls_targets
  workdir <- paste0("./runscMLnet/work_",inputs$parameters$Project)

  ## database

  LigRec.DB <- databases$LigRec.DB
  TFTG.DB <- databases$TFTG.DB
  RecTF.DB <- databases$RecTF.DB

  ## get LigRec subnetwork

  source_abundant <- ls_ligands[[ligclu]]
  cat("source_background:",length(source_abundant),"\n")
  target_abundant <- ls_receptors[[recclu]]
  cat("target_background:",length(target_abundant),"\n")
  tryCatch({
    LigRecTab <- getLigRec(LigRec.DB, source_abundant, target_abundant)
  }, error = function(e){
    cat(conditionMessage(e),"\n")
  })
  tag1 = exists("LigRecTab")
  if(!tag1) LigRecTab = data.frame()

  ## get TFTarget subnetwork

  target_background <- rownames(df_norm)[Matrix::rowMeans(df_norm)>0]
  cat("target_background:",length(target_background),"\n")
  target_interest <- ls_targets[[recclu]]
  if(ligclu %in% names(target_interest)){
    target_interest <- target_interest[[ligclu]]
  }
  cat("target_interest:",length(target_interest),"\n")
  tryCatch({
    TFTGTab <- getTFTG(TFTG.DB, target_interest, target_background)
  },error = function(e){
    cat(conditionMessage(e),"\n")
  })
  tag2 = exists("TFTGTab")
  if(!tag2) TFTGTab = data.frame()

  ## get RecTF subnetwork

  if(tag1 & tag2){
    Rec.list <- getNodeList(LigRecTab, "target")
    TF.list <- getNodeList(TFTGTab, "source")
    target.tfs <- TF.list
    tryCatch({
      RecTFTab <- getRecTF(RecTF.DB, Rec.list, TF.list)
    },error=function(e){
      cat(conditionMessage(e),"\n")
    })
  }
  tag3 = exists("RecTFTab")
  if(!tag3) RecTFTab = data.frame()

  ## update subnetwork

  if(tag3 & tag1){
    Receptors_in_Tab <- getNodeList(RecTFTab, "source")
    LigRecTab_new <- LigRecTab[LigRecTab[,2] %in% Receptors_in_Tab,]
  }else{
    LigRecTab_new <- data.frame()
  }
  cat("LR pairs:",nrow(LigRecTab_new),"\n")

  if(tag3 & tag2){
    TFs_in_Tab <- getNodeList(RecTFTab, "target")
    TFTGTab_new <- TFTGTab[TFTGTab[,1] %in% TFs_in_Tab,]
  }else{
    TFTGTab_new = data.frame()
  }
  cat("TFTG pairs:",nrow(TFTGTab_new),"\n")

  ## multiayer network

  mlnet <- list("LigRec" = LigRecTab_new,
                 "RecTF" = RecTFTab,
                 "TFTar" = TFTGTab_new)

  ## detail of mlnet

  detail <- c(source_abundant %>% length(),
              target_abundant %>% length(),
              target_interest %>% length(),
              nrow(mlnet$LigRec),nrow(mlnet$RecTF),nrow(mlnet$TFTar),
              ifelse(nrow(mlnet$LigRec)==0,0,mlnet$LigRec$source %>% unique() %>% length()),
              ifelse(nrow(mlnet$LigRec)==0,0,mlnet$LigRec$target %>% unique() %>% length()),
              ifelse(nrow(mlnet$TFTar)==0,0,mlnet$TFTar$source %>% unique() %>% length()),
              ifelse(nrow(mlnet$TFTar)==0,0,mlnet$TFTar$target %>% unique() %>% length()))

  ## output

  result <- list(mlnet=mlnet,detail=detail)

  ## save

  cellpair = paste(ligclu,recclu,sep = "_")
  workdir = paste(workdir,cellpair,sep = "/")
  dir.create(workdir,recursive = TRUE,showWarnings = F)
  saveRDS(mlnet, file = paste0(workdir,"/scMLnet.rds"))

  return(result)

}

#' @title getBarList
#' @description get barcode list of cluster of interset and other cluster.
#'
#' @param Aclu cluster of interset
#' @param BarCluTable Matrix, The annotation matrix for cell type annotation. The first column is barcode and the second is cell type.
#'
#' @return List. The first sublist is the barcodes of cluster of interset, the second sublist is the barcodes of other cluster.
#'
#' @import dplyr
#'
getBarList <- function(Aclu, BarCluTable)
{

  ## library

  # loadNamespace('dplyr')

  ## main

  AcluBar <- BarCluTable %>% dplyr::filter(.,Cluster == Aclu) %>%
    dplyr::select(.,Barcode) %>% unlist() %>% as.character()
  names(AcluBar) <- NULL

  AllBar <- BarCluTable %>% dplyr::select(Barcode) %>% unlist() %>% as.character()
  OtherBar <- setdiff(AllBar,AcluBar)

  result <- list(AcluBar,OtherBar)
  return(result)
}

runFisherTest <- function(subset1,subset2,backgrond)
{
  a=length(intersect(subset1,subset2))
  b=length(subset1)-a
  c=length(subset2)-a
  d=length(backgrond)-a-b-c
  matrix=matrix(c(a,c,b,d),nrow=2)
  fisher.test(matrix,alternative="greater")$p.value
}

getNodeList <- function(Database, Nodetype)
{

  ## library

  # loadNamespace('dplyr')

  ## main

  NodeList <- Database[,Nodetype] %>% unlist() %>% unique()

  return(NodeList)

}

getLigRec <- function(LigRec.DB, source_up, target_up)
{

  ## library

  # loadNamespace('dplyr')

  ## main

  if (!is.data.frame(LigRec.DB))
    stop("LigRec.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(LigRec.DB))
    stop("LigRec.DB must contain a column named 'source'")
  if (!"target" %in% colnames(LigRec.DB))
    stop("LigRec.DB must contain a column named 'target'")

  # get ligand and receptor list
  LigGene <- LigRec.DB %>% dplyr::select(source) %>% unlist() %>% unique()
  RecGene <- LigRec.DB %>% dplyr::select(target) %>% unlist() %>% unique()
  TotLigRec <- paste(LigRec.DB$source, LigRec.DB$target, sep = "_") %>% unique()

  # get high expressed ligand and receptor
  LigHighGene <- intersect(LigGene,source_up)
  RecHighGene <- intersect(RecGene,target_up)

  # get activated LR pairs
  LRList <- paste(rep(LigHighGene,each = length(RecHighGene)),RecHighGene,sep = "_")
  LRList <- intersect(LRList,TotLigRec)

  # check result
  if(length(LRList)==0)
    stop("Error: No significant LigRec pairs")

  # get result
  LRTable <- LRList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(LRTable) <- c("source","target")

  cat(paste0("get ",length(LRList)," activated LR pairs\n"))
  return(LRTable)

}

getTFTG <- function(TFTG.DB, target.degs, target.genes)
{

  ## library

  # loadNamespace('dplyr')

  ## main

  if (!is.data.frame(TFTG.DB))
    stop("TFTG.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(TFTG.DB))
    stop("TFTG.DB must contain a column named 'source'")
  if (!"target" %in% colnames(TFTG.DB))
    stop("TFTG.DB must contain a column named 'target'")

  # get TF list
  TF.list <- TFTG.DB %>% dplyr::select(source) %>% unlist() %>% unique()

  # get Target list
  TG.list <- lapply(TF.list, function(x){
    TFTG.DB %>% filter(source == x)  %>% dplyr::select(target) %>% unlist() %>% unique() %>% intersect(.,target.genes)
  })
  names(TG.list) <- TF.list

  # get target differently expressed genes
  DEGs <- target.degs

  # perform fisher test
  TFs <- lapply(TG.list, function(x){runFisherTest(subset1 = x, subset2 = DEGs, backgrond = target.genes)})
  TFs <- unlist(TFs)
  TFs <- names(TFs)[TFs <= 0.05]
  TFs <- TFs[TFs %in% target.genes]

  # get activated LR pairs
  TFTGList <- TG.list[TFs]
  TFTGList <- lapply(TFTGList, function(x){intersect(x, DEGs)})
  TFTGList <- paste(rep(TFs, times = lengths(TFTGList)), unlist(TFTGList), sep = "_")

  # check result
  if(length(TFTGList)==0)
    stop("Error: No significant TFTG pairs")

  # get result
  TFTGTable <- TFTGList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(TFTGTable) <- c("source","target")

  cat(paste0("get ",length(TFTGList)," activated TFTG pairs\n"))
  return(TFTGTable)

}

getRecTF <- function(RecTF.DB, Rec.list, TF.list)
{

  ## library

  # loadNamespace('dplyr')

  ## main

  if (!is.data.frame(RecTF.DB))
    stop("RecTF.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(RecTF.DB))
    stop("RecTF.DB must contain a column named 'source'")
  if (!"target" %in% colnames(RecTF.DB))
    stop("RecTF.DB must contain a column named 'target'")

  # make sure Rec.list in RecTF.DB
  Rec.list <- Rec.list[Rec.list %in% RecTF.DB$source]
  Rec.list <- as.vector(Rec.list)

  # make sure TF.list in RecTF.DB
  TF.list <- TF.list[TF.list %in% RecTF.DB$target]
  TF.list <- as.vector(TF.list)

  # get TF activated by Receptors
  TFofRec <- lapply(Rec.list, function(x){
    RecTF.DB %>% dplyr::filter(source == x)  %>% dplyr::select(target) %>% unlist() %>% unique()
  })
  names(TFofRec) <- Rec.list

  # get all TF
  TFofALL <- RecTF.DB %>% dplyr::select(target) %>% unlist() %>% unique()

  # perform fisher test
  Recs <- lapply(TFofRec, function(x){
    runFisherTest(subset1 = x, subset2 = TF.list, backgrond = TFofALL)
  })
  Recs <- unlist(Recs)
  Recs <- names(Recs)[Recs <= 0.05]
  # Recs <- Recs[Recs %in% target_gene]

  # get activated RecTF pairs
  RecTFList <- TFofRec[Recs]
  RecTFList <- lapply(RecTFList, function(x){intersect(x, TF.list)})
  RecTFList <- paste(rep(Recs, times = lengths(RecTFList)), unlist(RecTFList), sep = "_")

  # check result
  if(length(RecTFList)==0)
    stop("Error: No significant RecTF pairs")

  # get result
  RecTFTable <- RecTFList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(RecTFTable) <- c("source","target")

  cat(paste0("get ",length(RecTFList)," activated RecTF pairs\n"))
  return(RecTFTable)

}
