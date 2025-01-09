#' @title CheckFeedbackLoop
#' @description select multicellular feedback loop
#'
#' @param celltypes Vector, indicating the cell types in microenvironment.
#' @param workdir Character, the path where the multilayer network is stored.
#' @param savepath Character, the path where the result is stored.
#' @param cts_of_interest Vector or Character, the cell type your interest.
#'
#' @return list, containing the sender-receiver pair, and signals. The signals represent the multicellular feedback cruit between sender and receiver.
#' @export
CheckFeedbackLoop <- function(celltypes, wd_path, savepath, cts_of_interest){
  
  cp_of_inter <- data.frame(matrix(ncol = 4,dimnames = list(c(),c('ct1','ct2','keys_of_ct1','keys_of_ct2'))))
  
  for (i in 1:length(celltypes)) {
    
    ct1 <- celltypes[i]
    for (j in (i+1):length(celltypes)) {
      
      ct2 <- celltypes[j]
      cp1 <- paste0(ct1,'_',ct2)
      cp2 <- paste0(ct2,'_',ct1)
      if(cp1 %in% list.files(paste0(wd_path,'/runscMLnet/')) & cp2 %in% list.files(paste0(wd_path,'/runscMLnet/'))){
        
        cat('check in ',cp1,' and ',cp2,'\n')
        mlnet1 <- readRDS(paste0(wd_path,'/runscMLnet/',cp1,'/scMLnet.rds'))
        mlnet2 <- readRDS(paste0(wd_path,'/runscMLnet/',cp2,'/scMLnet.rds'))
        
        ct2_tgs <- unique(mlnet1$TFTar$target)
        ct1_tgs <- unique(mlnet2$TFTar$target)
        
        ct1_ligs <- unique(mlnet1$LigRec$source)
        ct2_ligs <- unique(mlnet2$LigRec$source)
        
        ct1_keys <- intersect(ct1_ligs,ct1_tgs)
        ct2_keys <- intersect(ct2_ligs,ct2_tgs)
        
        # add by yll 2024-12-04
        ct1_keys_filter <- c()
        ct2_keys_filter <- c()
        for (key in ct1_keys) {
          
          res_ct1_key1 <- process_ml_net(mlnet1, lig_key = key)
          res_ct1_loop <- res_ct1_key1[res_ct1_key1$Target %in% ct2_keys, ]
          
          res_ct2_key2 <- process_ml_net(mlnet2, lig_key = ct2_keys)
          res_ct2_loop <- res_ct2_key2[res_ct2_key2$Target %in% key, ]
          
          if (nrow(res_ct2_loop) > 0) {
            ct1_keys_filter <- c(ct1_keys_filter, key)
          }
          
        }
        
        for (key in ct2_keys) {
          
          res_ct1_key1 <- process_ml_net(mlnet1, lig_key = ct1_keys)
          res_ct1_loop <- res_ct1_key1[res_ct1_key1$Target %in% key, ]
          
          res_ct2_key2 <- process_ml_net(mlnet2, lig_key = key)
          res_ct2_loop <- res_ct2_key2[res_ct2_key2$Target %in% ct1_keys, ]
          
          if (nrow(res_ct1_loop) > 0) {
            ct2_keys_filter <- c(ct2_keys_filter, key)
          }
        }
        
        ct1_keys <- ct1_keys_filter
        ct2_keys <- ct2_keys_filter
        
        
        if(length(ct2_keys)>0|length(ct1_keys)>0){
          cp_of_inter <- rbind(cp_of_inter,c(ct1,ct2,length(ct1_keys),length(ct2_keys)))
        }
      }
    }
  }
  
  cp_of_inter <- na.omit(cp_of_inter)
  cp_of_inter <- cp_of_inter[cp_of_inter$keys_of_ct1 != 0 & cp_of_inter$keys_of_ct2 != 0,]
  rownames(cp_of_inter) <- 1:nrow(cp_of_inter)
  
  ## get genes list ####
  
  key_of_inter <- list()
  
  for (k in 1:nrow(cp_of_inter)) {
    
    ct1 <- cp_of_inter$ct1[k]
    ct2 <- cp_of_inter$ct2[k]
    
    mlnet1 <- readRDS(paste0(wd_path, '/runscMLnet/', paste(ct1, ct2, sep = '_'), '/scMLnet.rds'))
    mlnet2 <- readRDS(paste0(wd_path, '/runscMLnet/', paste(ct2, ct1, sep = '_'), '/scMLnet.rds'))
    
    ct2_tgs <- unique(mlnet1$TFTar$target)
    ct1_tgs <- unique(mlnet2$TFTar$target)
    ct1_ligs <- unique(mlnet1$LigRec$source)
    ct2_ligs <- unique(mlnet2$LigRec$source)

    ct1_keys <- intersect(ct1_ligs, ct1_tgs)
    ct2_keys <- intersect(ct2_ligs, ct2_tgs)
    
    # add by yll 2024-12-04
    ct1_keys_filter <- c()
    ct2_keys_filter <- c()
    for (key in ct1_keys) {
      
      res_ct1_key1 <- process_ml_net(mlnet1, lig_key = key)
      res_ct1_loop <- res_ct1_key1[res_ct1_key1$Target %in% ct2_keys, ]
      
      res_ct2_key2 <- process_ml_net(mlnet2, lig_key = ct2_keys)
      res_ct2_loop <- res_ct2_key2[res_ct2_key2$Target %in% key, ]
      
      if (nrow(res_ct2_loop) > 0) {
        ct1_keys_filter <- c(ct1_keys_filter, key)
      }
      
    }
    
    for (key in ct2_keys) {
      
      res_ct1_key1 <- process_ml_net(mlnet1, lig_key = ct1_keys)
      res_ct1_loop <- res_ct1_key1[res_ct1_key1$Target %in% key, ]
      
      res_ct2_key2 <- process_ml_net(mlnet2, lig_key = key)
      res_ct2_loop <- res_ct2_key2[res_ct2_key2$Target %in% ct1_keys, ]
      
      if (nrow(res_ct1_loop) > 0) {
        ct2_keys_filter <- c(ct2_keys_filter, key)
      }
    }
    
    ct1_keys <- ct1_keys_filter
    ct2_keys <- ct2_keys_filter
    
    key_of_inter[[paste(ct1, ct2, sep = '_')]] <- list(
      ct1_keys = ct1_keys,
      ct2_keys = ct2_keys
    )
  }
  
  ## transform to table ####
  # fbloop <- lapply(1:length(key_of_inter), function(i){
  #   
  #   ls <- key_of_inter[[i]]
  #   cp <- names(key_of_inter)[i]
  #   ct1 <- strsplit(cp,'_')[[1]][1]
  #   ct2 <- strsplit(cp,'_')[[1]][2]
  #   ct1_keys <- paste(ls$ct1_keys,collapse = ' ')
  #   ct2_keys <- paste(ls$ct2_keys,collapse = ' ')
  #   rbind(c(ct1,ct2,ct1_keys),c(ct2,ct1,ct2_keys))
  #   
  # }) %>% do.call('rbind',.) %>% as.data.frame()
  # colnames(fbloop) <- c('Sender','Receiver','Siganl')
  # 
  # if (length(cts_of_interest) == 0){
  #   write.csv(fbloop,paste0(savepath,'feddback_loop_detail.csv'))
  # }else{
  #   cts_of_interest <- cts_of_interest
  #   fbloop <- fbloop[fbloop$Sender %in% cts_of_interest & fbloop$Receiver %in% cts_of_interest,]
  #   write.csv(fbloop,paste0(savepath,'feddback_loop_detail.csv'))
  # }
 
  return(key_of_inter)
}


#' @title SaveFeedbackLoop
#' @description save multicellular feedback loop as table
#'
#' @param key_of_inter data.frame, the result of 'CheckFeedbackLoop' function.
#' @param savepath Character, the path where the result is stored.
#' @param cts_of_interest Vector or Character, the cell type your interest.
#'
#' @return Table, containing the sender, receiver, and signals. The signals represent the multicellular feedback cruit between sender and receiver.
#' @export
SaveFeedbackLoop <- function(key_of_inter,cts_of_interest,savepath){
  
  fbloop <- lapply(1:length(key_of_inter), function(i){
    
    ls <- key_of_inter[[i]]
    cp <- names(key_of_inter)[i]
    ct1 <- strsplit(cp,'_')[[1]][1]
    ct2 <- strsplit(cp,'_')[[1]][2]
    ct1_keys <- paste(ls$ct1_keys,collapse = ' ')
    ct2_keys <- paste(ls$ct2_keys,collapse = ' ')
    rbind(c(ct1,ct2,ct1_keys),c(ct2,ct1,ct2_keys))
    
  }) %>% do.call('rbind',.) %>% as.data.frame()
  colnames(fbloop) <- c('Sender','Receiver','Siganl')
  
  if (length(cts_of_interest) == 0){
    write.csv(fbloop,paste0(savepath,'feddback_loop_detail.csv'))
  }else{
    cts_of_interest <- cts_of_interest
    fbloop <- fbloop[fbloop$Sender %in% cts_of_interest & fbloop$Receiver %in% cts_of_interest,]
    write.csv(fbloop,paste0(savepath,'feddback_loop_detail.csv'))
  }
  return(fbloop)
}


process_ml_net <- function(mlnet, lig_key = NULL, tar_key = NULL) {
  if (!is.null(lig_key)) {
    mlnet$LigRec <- mlnet$LigRec[mlnet$LigRec$source %in% lig_key, ]
  }
  if (!is.null(tar_key)) {
    mlnet$TFTar <- mlnet$TFTar[mlnet$TFTar$target %in% tar_key, ]
  }
  if (!is.null(mlnet$LigRec)) {
    mlnet$RecTF <- mlnet$RecTF[mlnet$RecTF$source %in% mlnet$LigRec$target, ]
  }
  if (!is.null(mlnet$TFTar)) {
    mlnet$RecTF <- mlnet$RecTF[mlnet$RecTF$target %in% mlnet$TFTar$source, ]
  }
  if (!is.null(mlnet$RecTF)) {
    mlnet$LigRec <- mlnet$LigRec[mlnet$LigRec$target %in% mlnet$RecTF$source, ]
  }
  
  ligrec <- data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target)
  rectf <- data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
  tftg <- data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
  res <- ligrec %>% 
    merge(., rectf, by = "Receptor") %>% 
    merge(., tftg, by = "TF") %>% 
    dplyr::select(Ligand, Receptor, TF, Target) %>% 
    arrange(Ligand, Receptor)
  
  return(res)
}


#' @title plot_feedback_loop
#' @description plot multicellular feedback loop
#'
#' @param floop table, containing the sender, receiver, and signals. The signals represent the multicellular feedback cruit between sender and receiver.
#' @param savepath Character, the path where the result is stored.
#' @param vertex.size Numercial, the size of nodes.
#' @param p_width Numercial, the width of plot.
#' @param p_height Numercial, the height of plot.
#' 
#' @export
plot_feedback_loop <- function(floop,save_path,vertex.size,p_width,p_height){
  suppressMessages(library(igraph))
  
  nodes <- c(fbloop$Sender, fbloop$Receiver) %>% unlist(.) %>% unique()
  edges <- list()  
  
  for (i in 1:length(fbloop$Sender)) {
    edge_i <- c(fbloop$Sender[i], fbloop$Receiver[i])
    edges <- append(edges, list(edge_i))  
  }
  edges <- unlist(edges)
  vertex_colors <- mycolor_ct[names(mycolor_ct)  %in% nodes]
  graph <- graph(edges = edges, directed = TRUE)
  
  edge_labels <- list()
  for (i in 1:length(fbloop$Siganl)){
    sig <- strsplit(fbloop$Siganl[i]," ")[[1]]
    if (length(sig)!=1){
      label <- paste0(sig[1],"/",sig[2])
    }else{
      label <- sig
    }
    edge_labels <- append(edge_labels, label)
    
    
  }
  edge_labels <- unlist(edge_labels)
  
  plotdir = paste0(save_path,'/visualize/')
  dir.create(plotdir,recursive = T)
  pdf(paste0(plotdir,"plot_feedback_loop_all_cp.pdf"), width = p_width, height = p_height)  
  p1 <- plot(
    graph, 
    #layout = layout,              
    vertex.color = vertex_colors,   
    vertex.size = vertex.size,             
    vertex.label.color = "black", 
    edge.width = 1.5,
    edge.arrow.size = 1,        
    edge.color = "gray",          
    edge.label = edge_labels,     
    edge.label.cex = 0.8,         
    edge.curved = 0.2,             
    edge.label.cex = 0.8,       
    edge.label.dist = 3,    
    edge.label.angle = 0.5,        
    edge.label.color = 'black'
  )
  dev.off()
}

#' @title plot_feedback_network
#' @description plot multilayer signaling feedback network
#'
#' @param wd_path Character, the path where the stMLnet inference result is stored..
#' @param ct1 Character, sender cell.
#' @param ct2 Character, receiver cell.
#' @param ct1_key Character, the signaling key gene of sender cell.
#' @param ct2_key Character, the signaling key gene of receiver cell.
#' @param vertex.size Numercial, the size of nodes.
#' @param rescale Character, Whether to automatically scale.
#' @param do.check Character, check the result for better visualization.
#' @param p_width Numercial, the width of plot.
#' @param p_height Numercial, the height of plot.
#'
#' @export
plot_feedback_network <- function(wd_path, ct1,ct2,ct1_key,ct2_key,vertex.size,rescale=FALSE,do.check=TRUE,p_width=10,p_height=6) {
  
  mlnet1 <- readRDS(paste0(wd_path, '/runscMLnet/', paste(ct1, ct2, sep = '_'), '/scMLnet.rds'))
  mlnet2 <- readRDS(paste0(wd_path, '/runscMLnet/', paste(ct2, ct1, sep = '_'), '/scMLnet.rds'))
  
  res_ct1_key1 <- process_ml_net(mlnet1, lig_key = ct1_key)
  res_ct1_loop <- res_ct1_key1[res_ct1_key1$Target %in% ct2_key,]
  
  res_ct2_key2 <- process_ml_net(mlnet2, lig_key = ct2_key)
  res_ct2_loop <- res_ct2_key2[res_ct2_key2$Target %in% ct1_key, ]
  
  if (dim(res_ct1_loop)[1]==0){
    cat('Cannot find',ct2_key,'in Target of', ct1,'->',ct2)
    return()  
  }else if(dim(res_ct2_loop)[1]==0){
    cat('Cannot find',ct1_key,'in Target of', ct2,'->',ct1)
    return() 
  }
  
  if (do.check){
    
    if (length(unique(res_ct1_loop$Receptor))>2){
      cut_recs <- unique(res_ct1_loop$Receptor)[1:2]
      res_ct1_loop_cut <- res_ct1_loop[res_ct1_loop$Receptor %in% cut_recs,] 
      res_ct1_loop <- res_ct1_loop_cut
    }
    
    if(length(unique(res_ct1_loop$TF))>2){
      cut_tfs <- unique(res_ct1_loop$TF)[1:2]
      res_ct1_loop_cut <- res_ct1_loop[res_ct1_loop$TF %in% cut_tfs,] 
      res_ct1_loop <- res_ct1_loop_cut
    }
    
    if (length(unique(res_ct2_loop$Receptor))>2){
      cut_recs <- unique(res_ct2_loop$Receptor)[1:2]
      res_ct2_loop_cut <- res_ct2_loop[res_ct2_loop$Receptor %in% cut_recs,] 
      res_ct2_loop <- res_ct2_loop_cut
    }
    
    if(length(unique(res_ct2_loop$TF))>2){
      cut_tfs <- unique(res_ct2_loop$TF)[1:2]
      res_ct2_loop_cut <- res_ct2_loop[res_ct2_loop$TF %in% cut_tfs,] 
      res_ct2_loop <- res_ct2_loop_cut
    }
  }
  
  # if (res_ct1_loop$Ligand==res_ct1_loop$Receptor){
  #   res_ct1_loop$Receptor <- paste0(res_ct1_loop$Receptor,'.rec')
  # }
  
  if (dim(res_ct1_loop)[1]==1){
    df_ct1_loop <- res_ct1_loop
    df_ct1_loop$Receptor <- paste0(df_ct1_loop$Receptor,".rec1")
    df_ct1_loop$TF <- paste0(df_ct1_loop$TF,".tf1")
    df_ct1_loop$Target <- paste0(df_ct1_loop$Target,".tg")
  }else{
    df_ct1_loop <- data.frame(Ligand = NA, Receptor = NA, TF = NA, Target = NA)
    df_ct1_loop$Ligand <- unique(res_ct1_loop$Ligand)
    if (length(unique(res_ct1_loop$Receptor))==1){
      df_ct1_loop$Receptor = unique(res_ct1_loop$Receptor)
      df_ct1_loop$Receptor <- paste0(unique(res_ct1_loop$Receptor), ".rec1")
    }else{
      df_ct1_loop$Receptor = paste0(paste(unique(res_ct1_loop$Receptor), collapse = "/"), ".rec1")
    }
    
    if (length(unique(res_ct1_loop$TF))==1){
      df_ct1_loop$TF = unique(res_ct1_loop$TF)
      df_ct1_loop$TF <- paste0(unique(res_ct1_loop$TF), ".tf1")
    }else{
      df_ct1_loop$TF = paste0(paste(unique(res_ct1_loop$TF), collapse = "/"), ".tf1")
    }
    
    if (length(unique(res_ct1_loop$Target))==1){
      df_ct1_loop$Target = unique(res_ct1_loop$Target)
      df_ct1_loop$Target <- paste0(unique(res_ct1_loop$Target), ".tg")
    }else{
      df_ct1_loop$Target = paste(unique(res_ct1_loop$Target), collapse = "/",".tg")
    }
  }
  
  
  if (dim(res_ct2_loop)[1]==1){
    df_ct2_loop <- res_ct2_loop
    df_ct2_loop$Receptor <- paste0(df_ct2_loop$Receptor,".rec2")
    df_ct2_loop$TF <- paste0(df_ct2_loop$TF,".tf2")
    df_ct2_loop$Ligand <- paste0(df_ct2_loop$Ligand,".tg")
  }else{
    df_ct2_loop <- data.frame(Ligand = NA, Receptor = NA, TF = NA, Target = NA)
    df_ct2_loop$Ligand <- unique(res_ct2_loop$Ligand)
    df_ct2_loop$Ligand <- paste0(unique(res_ct2_loop$Ligand), ".tg")
    if (length(unique(res_ct2_loop$Receptor))==1){
      df_ct2_loop$Receptor = unique(res_ct2_loop$Receptor)
      df_ct2_loop$Receptor <- paste0(unique(res_ct2_loop$Receptor), ".rec2")
    }else{
      df_ct2_loop$Receptor = paste0(paste(unique(res_ct2_loop$Receptor), collapse = "/"), ".rec2")
    }
    
    if (length(unique(res_ct2_loop$TF))==1){
      df_ct2_loop$TF = unique(res_ct2_loop$TF)
      df_ct2_loop$TF <- paste0(unique(res_ct2_loop$TF), ".tf2")
    }else{
      df_ct2_loop$TF = paste0(paste(unique(res_ct2_loop$TF), collapse = "/"), ".tf2")
    }
    
    if (length(unique(res_ct2_loop$Target))==1){
      df_ct2_loop$Target = unique(res_ct2_loop$Target)
    }else{
      df_ct2_loop$Target = paste(unique(res_ct2_loop$Target), collapse = "/")
    }
    
  }
  
  nodes <- c(df_ct1_loop, df_ct2_loop) %>% unlist(.) %>% unique()
  
  edges1 <- c(
    df_ct1_loop$Ligand, df_ct1_loop$Receptor,
    df_ct1_loop$Receptor, df_ct1_loop$TF,
    df_ct1_loop$TF, df_ct1_loop$Target,
    df_ct2_loop$Ligand, df_ct2_loop$Receptor,
    df_ct2_loop$Receptor, df_ct2_loop$TF,
    df_ct2_loop$TF, df_ct1_loop$Ligand
  )
  
  coords <- matrix(c(0, 1,    
                     1, 1,    
                     2, 1,    
                     2, 0,    
                     1, 0,    
                     0, 0),   
                   ncol = 2, byrow = TRUE)
  coords[, 1] <- coords[, 1] * 1
  
  x_offset <- (max(coords[, 1]) - min(coords[, 1])) / 2
  y_offset <- (max(coords[, 2]) - min(coords[, 2])) / 2
  coords[, 1] <- coords[, 1] - min(coords[, 1]) - x_offset
  coords[, 2] <- coords[, 2] - min(coords[, 2]) - y_offset
  
  vertex_colors <- c("#D43F3ACC",   
                     "#EEA236CC",  
                     "#5CB85CCC",  
                     "#D43F3ACC",  
                     "#EEA236CC",  
                     "#5CB85CCC")  
  
  edge_lty <- c(2, 2, 2, 1, 1, 1, 1)  
  
  plotdir = paste0(wd_path,'/visualize/')
  dir.create(plotdir,recursive = T)
  
  pdf(paste0(plotdir,"plot_feedback_",ct1,"_",ct1_key,"_",ct2,"_",ct2_key,"_1.pdf"), width = p_width, height = p_height)  # 设置文件名和图像大小
  plot(
    graph(edges = edges1, directed = TRUE), 
    layout = coords,             
    vertex.size = vertex.size,   
    vertex.label.cex = 1,        
    vertex.color = vertex_colors,
    vertex.frame.color = "black",
    edge.arrow.size = 0.9,       
    edge.color = "#696969",      
    vertex.label.color = "black",
    rescale = rescale,             
    edge.lty = edge_lty          
  )

  legend(
    x = 0.2, y = -0.7,
    legend = c(paste0(ct1,"->",ct2),    
               paste0(ct2,"->",ct1)),
    lty = c(2, 1),                  
    col = "#696969",               
    title = "Edge Types",           
    cex = 0.8                       
  )

  legend(
    x = -0.6, y = -0.7,         
    legend = c("Ligand/Target", "Receptor", "TF"), 
    fill = c("#D43F3ACC", "#EEA236CC", "#5CB85CCC"), 
    title = "Nodes",                  
    cex = 0.8                        
  )
  
  dev.off()
}



