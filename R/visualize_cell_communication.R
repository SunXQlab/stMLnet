#' @title PrepareColorDB
#' @description Prepare color database for visualization
#'
#' @param CellTypes Vector, indicating the cell types in microenvironment.
#' @param ColorDB Vector, indicating the colors of cell types.
#'
#' @return List, containing the color of cell types, key siganls and others.
#' @export
#' @import ggsci
PrepareColorDB <- function(CellTypes, ColorDB = NULL){

  # loadNamespace('ggsci')

  ## celltype

  if(is.null(ColorDB)){

    if(length(CellTypes)>8){

      mycolors_nejm <- pal_d3(palette = "category20", alpha = 0.8)(20)
      celltype <- CellTypes
      mycolor_ct <- mycolors_nejm[1:length(celltype)]
      names(mycolor_ct) <- celltype

    }else{

      mycolors_nejm <- pal_nejm(palette = "default", alpha = 0.8)(8)
      celltype <- CellTypes
      mycolor_ct <- mycolors_nejm[1:length(celltype)]
      names(mycolor_ct) <- celltype

    }

  }else{

    if(length(ColorDB)!=length(CellTypes))
      stop('The length of ColorDB is not equal to the length of CellTypes')
    mycolor_ct <- ColorDB
    names(mycolor_ct) <- CellTypes

  }

  ## nodekey

  mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)
  nodekey <- c("Ligand","Receptor","TF","Target")
  mycolor_key <- mycolors_locus[1:4]
  names(mycolor_key) <- nodekey

  ## nodetype

  mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)
  nodetype <- c("cell","Sender","Receiver")
  mycolor_nt <- mycolors_locus[1:3]
  names(mycolor_nt) <- nodetype

  ## output

  result <- list(Celltypes=mycolor_ct, Keys=mycolor_key, Nodes=mycolor_nt)

  return(result)

}

#' @title DrawNetworkPlot
#' @description Draw Network Plot
#'
#' @param InputDir Character, the path where the signals importance is stored.
#' @param Metric Character, indicating the metrics of edges in Network plot. Available options are: n_LRs, n_TGs, IM,  IM_norm,  mean_IM, mean_IM_norm
#' @param ColorDB Vector, indicating the colors of cell types.
#' @param gtitle Character, the title of plot.
#'
#' @export
#' @import dplyr graphics grDevices
#'
DrawNetworkPlot <- function(InputDir, Metric, ColorDB, gtitle = 'CCI'){

  inputdir <- InputDir
  key <- Metric
  colodb <- ColorDB

  outputdir <- './visualize_CCI/NetworkPlot/'
  dir.create(outputdir,recursive=T,showWarnings = F)

  files <- list.files(inputdir)[grep('LRTG_im_clean_',list.files(inputdir))]
  LRTG_detail <- lapply(files, function(f){

    # f = files[1]
    LRTG_im <- readRDS(paste0(inputdir,f))
    c(length(unique(LRTG_im$LRpair)),length(unique(LRTG_im$Target)),
      sum(LRTG_im$IM),sum(LRTG_im$im_norm),
      sum(LRTG_im$IM)/nrow(LRTG_im),sum(LRTG_im$im_norm)/nrow(LRTG_im))

  }) %>% do.call('rbind',.) %>% as.data.frame()

  df_cellpair <- gsub('LRTG_im_clean_|.rds','',files) %>% strsplit(.,"-") %>% do.call('rbind',.) %>% as.data.frame()
  LRTG_detail <- cbind(df_cellpair,LRTG_detail)
  colnames(LRTG_detail) <- c('cell_from','cell_to','n_LRs','n_TGs','IM','IM_norm','mean_IM','mean_IM_norm')

  tmeTab <- LRTG_detail[,c('cell_from','cell_to',key)]
  colnames(tmeTab) <- c('cell_from','cell_to','n')

  png(paste0(outputdir,"networkPlot_",key,".png"),
      height = 7,width = 9, units = 'in', res = 300)
  DrawCellComm(tmeTab,colodb = colodb, gtitle = gtitle)
  dev.off()

  pdf(paste0(outputdir,"networkPlot_",key,".pdf"),height = 7,width = 9)
  DrawCellComm(tmeTab,colodb = colodb, gtitle = gtitle)
  dev.off()

  DrawCellComm(tmeTab,colodb = colodb, gtitle = gtitle)

}

DrawCellComm <- function(CellTab,colodb,gtitle){

  # loadNamespace("igraph")
  # loadNamespace('plotrix')

  aaa <- CellTab
  g <- igraph::graph.data.frame(aaa,directed=TRUE)
  alltype <- c(aaa$cell_from,aaa$cell_to)
  alltype <- unique(alltype)
  # allcolos <- rainbow(length(alltype))
  ChoColorNum <- alltype# 1:length(alltype) # sample(1:length(colodb),length(alltype), replace = FALSE)

  allcolos <- colodb[ChoColorNum]
  names(allcolos) <- alltype

  edge.start <- igraph::ends(g,es=igraph::E(g),names=FALSE)

  layout <- igraph::in_circle()  #igraph packages
  coords <- igraph::layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }

  loop.angle <- ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))

  vertex.label.color <- 'black'
  igraph::V(g)$size <- 20
  igraph::V(g)$color <- allcolos[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- 0

  label <- FALSE
  if(label){
    igraph::E(g)$label<-igraph::E(g)$n
  }

  edge.max.width = 10
  if(max(igraph::E(g)$n)==min(igraph::E(g)$n)){
    igraph::E(g)$width <- 2
  }else{
    igraph::E(g)$width <- 1+edge.max.width/(max(igraph::E(g)$n)-min(igraph::E(g)$n))*(igraph::E(g)$n-min(igraph::E(g)$n))
  }

  edge.label.color <- "black"
  igraph::E(g)$arrow.width<- 1
  igraph::E(g)$arrow.size <- 0.1
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- 1
  igraph::E(g)$color <- igraph::V(g)$color[edge.start[,1]]

  if(sum(edge.start[,2]==edge.start[,1])!=0)
  {
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }

  #draw
  shape = 'circle'
  margin = 0.2
  edgecurved = 0.2
  vertexlabelcex = 1
  plot(
    g,
    edge.curved=edgecurved,
    vertex.label = "",
    vertex.shape=shape,
    layout=coords_scale,
    margin=margin,
    vertex.shape="fcircle",
    vertex.label.cex = vertexlabelcex,
    axes = FALSE
  )

  #draw legend
  xstart <- 1.5
  ystart <- 0.5
  gap <- 0.02
  cirr <- 0.04
  DrawLeg(alltype,allcolos,xstart,ystart,cirr,gap)

  #title
  title(main = list(gtitle, cex=2))

}

DrawLeg <- function(alltype,allcolos,xstart,ystart,cirr,gap){

  # loadNamespace('igraph')
  # loadNamespace('plotrix')

  ThisX <- xstart
  ThisY <- ystart
  for(i in seq(1,length(alltype)))
  {
    ThisType <- alltype[i]
    plotrix::draw.circle(ThisX,ThisY,cirr,col = allcolos[ThisType])
    text(ThisX+cirr+0.1,ThisY,ThisType,adj = 0)

    ThisY <- ThisY - (2*cirr + gap)
  }
}

#' @title DrawMLnetPlot
#' @description Draw MLnet Plot
#'
#' @param MLnetDir Character, the path where the multilayer network is stored.
#' @param ImportDir Character, the path where the signals importance is stored.
#' @param Signal Character, downstream multilayer network of specific signal are shown.
#' @param ColorDB Vector, indicating the colors of key siganls including Ligands, Recetors, TFs and Target genes.
#' @param Check Logical, Whether to do check visualization.
#' @param top.n Numercial, Only relevant if Check=T. Number of top regulation between LR pairs and target genes to show according to the signals importance.
#' @param gtitle  Character, the title of plot.
#' @param p_height Numercial, the height of plot.
#' @param p_width Numercial, the width of plot.
#'
#' @export
#' @import dplyr ggplot2 grDevices
DrawMLnetPlot <- function(MLnetDir, ImportDir, Signal, ColorDB,
                          Check = TRUE, top.n = 10, gtitle = 'CCI',
                          p_height = 4.5, p_width = 7){

  outputdir <- './visualize_CCI/MLnetPlot/'
  dir.create(outputdir,recursive=T,showWarnings = F)

  MLnet <- readRDS(MLnetDir)
  LRTG_im <- readRDS(ImportDir)
  Key <- Signal

  ligands <- unique(MLnet$LigRec$source)
  if(Key %in% ligands){
    Type <- 'Ligand'
  }else{
    Type <- 'Receptor'
  }

  MLnet_key <- prepareMLnetworkPlotData(mlnet=MLnet,lrtg_im=LRTG_im,Key=Key,Type=Type,do.check=Check,top.n=top.n)
  pt <- drawMLnetworkPlot(mlnet=MLnet_key,colodb=ColorDB)

  ## save

  pdf(paste0(outputdir,'/MLnetPlot-',gtitle,'.pdf'),height = p_height,width = p_width)
  print(pt)
  dev.off()

  png(paste0(outputdir,'/MLnetPlot-',gtitle,'.png'),height = p_height,width = p_width, units = 'in', res = 300)
  print(pt)
  dev.off()

  return(pt)

}


prepareMLnetworkPlotData <- function(mlnet = NULL, lrtg_im = NULL,Key, Type, do.check = TRUE, top.n = 10){

  # loadNamespace('dplyr')
  # loadNamespace('ggraph')

  if(Type == 'Ligand'){

    mlnet$LigRec <- mlnet$LigRec[mlnet$LigRec$source %in% Key,]

  }else{

    mlnet$LigRec <- mlnet$LigRec[mlnet$LigRec$target %in% Key,]

  }
  mlnet$RecTF <- mlnet$RecTF[mlnet$RecTF$source %in% mlnet$LigRec$target,]
  mlnet$TFTar <- mlnet$TFTar[mlnet$TFTar$source %in% mlnet$RecTF$target,]
  lrtg_im = lrtg_im %>% dplyr::select(Ligand,Receptor,Target,Score = im_norm) %>% dplyr::filter(Receptor %in% Key)

  if(do.check & nrow(lrtg_im)>top.n){

    message('using the score of LRTG to narrow down for better visualization')
    tg_im_check <- lrtg_im %>% group_by(Target) %>%
      summarise(sum_Score=sum(Score)) %>% arrange(desc(sum_Score)) %>%
      dplyr::select(Target) %>% unlist() %>% head(top.n)

    mlnet_check <- mlnet
    mlnet_check$TFTar <- mlnet_check$TFTar[mlnet_check$TFTar$target %in% tg_im_check,]
    mlnet_check$RecTF <- mlnet_check$RecTF[mlnet_check$RecTF$target %in% mlnet_check$TFTar$source,]
    mlnet_check$LigRec <- mlnet_check$LigRec[mlnet_check$LigRec$target %in% mlnet_check$RecTF$source,]
    mlnet <- mlnet_check

  }


  return(mlnet)
}


drawMLnetworkPlot <- function(mlnet, colodb){

  # loadNamespace('dplyr')
  # loadNamespace('ggraph')

  Ligs <- unique(mlnet$LigRec$source)
  Recs <- unique(mlnet$LigRec$target)
  TFs <- unique(mlnet$TFTar$source)
  TGs <- unique(mlnet$TFTar$target)
  TGs <- TGs[!TGs %in% c(Ligs,Recs,TFs)]
  df_nodes <- data.frame(node = c(Ligs,Recs,TFs,TGs),
                         key = c(rep('Ligand',length(Ligs)),
                                 rep('Receptor',length(Recs)),
                                 rep('TF',length(TFs)),
                                 rep('Target',length(TGs))))
  df_nodes$color <- colodb[df_nodes$key]
  df_edges <- do.call("rbind",list(mlnet$LigRec,mlnet$RecTF,mlnet$TFTar))

  subnet <- igraph::graph_from_data_frame(d = df_edges, vertices = df_nodes)
  df_nodes <- df_nodes[match(names(igraph::V(subnet)),df_nodes$node),]
  root_index <- grep('Ligand',df_nodes$key)
  set.seed(4)
  coords <- igraph::layout_(subnet,layout = igraph::as_tree(root = root_index)) %>% as.data.frame()
  coords <- cbind(coords,df_nodes$key)
  colnames(coords) <- c('dim_x','dim_y','type')

  dist_TG <- 1;len_TG <- table(coords$type)[['Target']]
  dist_TF <- 1.5;len_TF <- table(coords$type)[['TF']]
  dist_Rec <- 2.5;len_Rec <- table(coords$type)[['Receptor']]
  dist_Lig <- 2.5;len_Lig <- table(coords$type)[['Ligand']]
  if(len_Lig==1){
    coords$dim_x[coords$type == 'Ligand'] = 0
  }else{
    dim_x_1 = seq(to = -dist_Lig/2,by = dist_Lig,length.out = ceiling(len_Lig/2))
    dim_x_2 = seq(from = dist_Lig/2,by = dist_Lig,length.out = len_Lig-ceiling(len_Lig/2))
    coords$dim_x[coords$type == 'Ligand'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  if(len_Rec==1){
    coords$dim_x[coords$type == 'Receptor'] = 0
  }else{
    dim_x_1 = seq(to = -dist_Rec/2,by = dist_Rec,length.out = ceiling(len_Rec/2))
    dim_x_2 = seq(from = dist_Rec/2,by = dist_Rec,length.out = len_Rec-ceiling(len_Rec/2))
    coords$dim_x[coords$type == 'Receptor'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  if(len_TG<len_TF & len_TF>10){

    dist_TG <- 1
    dist_TF <- 0.5

  }else if(len_TG>len_TF & len_TG>10){

    dist_TG <- 0.5
    dist_TF <- 1

  }
  if(len_TF==1){
    coords$dim_x[coords$type == 'TF'] = 0
  }else{
    dim_x_1 = seq(to = -dist_TF/2,by = dist_TF,length.out = ceiling(len_TF/2))
    dim_x_2 = seq(from = dist_TF/2,by = dist_TF,length.out = len_TF-ceiling(len_TF/2))
    coords$dim_x[coords$type == 'TF'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  if(len_TG==1){
    coords$dim_x[coords$type == 'Target'] = 0
  }else{
    dim_x_1 = seq(to = -dist_TG/2,by = dist_TG,length.out = ceiling(len_TG/2))
    dim_x_2 = seq(from = dist_TG/2,by = dist_TG,length.out = len_TG-ceiling(len_TG/2))
    coords$dim_x[coords$type == 'Target'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  coords$dim_y <- lapply(coords$type,switch,'Ligand'=0.9,'Receptor'=0.6,'TF'=0.3,'Target'=0) %>% unlist()
  # plot(subnet, layout = as.matrix(coords[,1:2]), vertex.color=df_nodes$color)

  layout <- ggraph::create_layout(subnet, layout = 'tree')
  layout[,1:2] <- coords[,1:2]
  # head(layout)

  temp <- function(key){

    max(coords$dim_x[coords$type==key])+0.5

  }
  df_anno <- data.frame(x = lapply(unique(layout$key), temp) %>% unlist(),
                        y = c(0.92,0.62,0.32,0.02),
                        lab = unique(layout$key),
                        key = unique(layout$key))

  pt <- ggraph::ggraph(layout) +
    ggraph::geom_edge_link(aes(),color="grey",
                   arrow = arrow(length = unit(1.5, 'mm')),
                   start_cap = ggraph::circle(3, 'mm'),
                   end_cap = ggraph::circle(3, 'mm')) +
    ggraph::geom_node_point(aes(fill = key,color = key),shape=21,size = 8) +
    ggraph::geom_node_text(aes(label=name),size=2) + # ,fontface='bold',family='Arial'
    xlim(c(min(layout$x)-2,max(layout$x)+2)) +
    geom_text(data = df_anno, aes(x,y,label=lab,color=key),
              vjust = 1, hjust = 0, size = 4,fontface='bold') + # ,family='ARL'
    scale_fill_manual(values = colodb) + scale_color_manual(values = colodb) +
    guides(fill='none',color='none') +
    ggraph::theme_graph()
  print(pt)

  return(pt)

}

#' @title DrawEdgeBundlingPlot
#' @description Draw EdgeBundling Plot
#'
#' @param InputDir Character, the path where the signals acitivity is stored.
#' @param Cluster  Character, signals acitivity of upstream signal pairs of specific cluster are shown.
#' @param ColorDB Vector, containing the color of cell types, key siganls and others.
#' @param Check Logical, Whether to do check visualization.
#' @param top.n Numercial, Only relevant if Check=T. Number of top regulation between ligands and receptor to show according to the signals acitivity.
#' @param gtitle  Character, the title of plot.
#' @param p_height Numercial, the height of plot.
#' @param p_width Numercial, the width of plot.
#'
#' @export
#' @import dplyr ggplot2 grDevices
DrawEdgeBundlingPlot <- function(
  InputDir, Cluster, ColorDB,
  Check = TRUE, top.n = 50, gtitle = 'CCI',
  p_height = 7.5, p_width = 7
){

  inputdir <- InputDir
  cluster <- Cluster
  colodb <- ColorDB

  outputdir <- './visualize_CCI/EdgeBundlingPlot/'
  dir.create(outputdir,recursive=T,showWarnings = F)

  files <-  list.files(inputdir)
  files <- files[grep(paste0('-',cluster,'.rds'),files)]
  LRTGscore = lapply(files, function(file){

    print(file)
    LRS_score = readRDS(paste0(inputdir,file))[[1]]
    LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]
    if(is.null(dim(LRS_score_merge))) LRS_score_merge = LRS_score[[1]]

    file <- gsub('-','_',file)
    df_LigRec <- data.frame(
      source = colnames(LRS_score_merge) %>% gsub('_.*','',.),
      target = colnames(LRS_score_merge) %>% gsub('.*_','',.),
      LRpair = colnames(LRS_score_merge),
      count = colMeans(LRS_score_merge),
      source_group = strsplit(file,'[_//.]')[[1]][3],
      target_group = strsplit(file,'[_//.]')[[1]][4]
    )

  }) %>% do.call('rbind',.)

  # plot

  df_input <- prepareEdgeBundlingPlotData(LRTGscore = LRTGscore, do.check = Check, top.n = top.n)
  pt <- drawEdgeBundlingPlot(df_input = df_input, colodb = colodb)

  ## save

  pdf(paste0(outputdir,'EdgeBundlingPlot-',gtitle,'.pdf'),height = p_height, width = p_width)
  print(pt)
  dev.off()

  png(paste0(outputdir,'EdgeBundlingPlot-',gtitle,'.png'),
      height = p_height, width = p_width, units = 'in', res = 300)
  print(pt)
  dev.off()

  return(pt)

}

prepareEdgeBundlingPlotData <- function(LRTGscore, do.check = TRUE, top.n = 50){

  # edge
  LRTGscore$from = paste(LRTGscore$source,LRTGscore$source_group,sep = "_")
  LRTGscore$to = paste(LRTGscore$target,LRTGscore$target_group,sep = "_")
  df_edge <- LRTGscore[,grep("from|to",colnames(LRTGscore))]

  # hierarchy
  d1 <- data.frame(from="cell", to=c('Sender','Receiver'))
  d2 <- data.frame(from=c(rep('Sender',length(unique(LRTGscore$source_group))),'Receiver'),
                   to=c(unique(LRTGscore$source_group),unique(LRTGscore$target_group)))
  d3 <- data.frame(from=c(LRTGscore$source_group,LRTGscore$target_group),
                   to=c(df_edge$from,df_edge$to))
  d3 <- d3[!duplicated(d3),]
  df_hierarchy <- do.call('rbind',list(d1,d2,d3))

  # node
  df_node  <-  data.frame(
    name = unique(c(as.character(df_hierarchy$from), as.character(df_hierarchy$to)))
  )
  df_node$label <- gsub("_.*","",df_node$name)
  df_node$group <- df_hierarchy$from[match(df_node$name, df_hierarchy$to)]
  df_node$group <- factor(df_node$group,levels=unique(df_node$group))
  df_node$score <- LRTGscore$count[match(df_node$name,LRTGscore$from)]
  df_node$score <- replace(df_node$score,is.na(df_node$score),1e-5)

  # check size
  if(do.check & nrow(df_edge)>top.n){

    message('try to narrow down the number of LR for better visualization')
    df_node$rank <- Inf
    df_node$rank[df_node$group %in% d2$to[d2$from == 'Sender']] <- df_node$score[df_node$group %in% d2$to[d2$from == 'Sender']] %>% rank()

    df_node_check <- df_node
    df_node_check <- df_node_check[df_node_check$rank >= max(df_node_check$rank[is.finite(df_node_check$rank)])-top.n+1,]

    df_edge_check <- df_edge[df_edge$from %in% df_node_check$name,]

    df_node_check <- df_node_check[df_node_check$name %in%
                                     c(d1$from,d1$to,
                                       df_node_check$group %>% as.character(),
                                       df_edge_check$from,df_edge_check$to),]

    df_hierarchy_check <- df_hierarchy[df_hierarchy$to %in% df_node_check$name,]

    df_node <- df_node_check
    rownames(df_node) <- seq(nrow(df_node))
    df_edge <- df_edge_check
    rownames(df_edge) <- seq(nrow(df_edge))
    df_hierarchy <- df_hierarchy_check
    rownames(df_hierarchy) <- seq(nrow(df_hierarchy))

  }

  # attr

  index_myleaves <- grep('^cell$|Sender|Receiver',df_node$group,invert = T)[-1]
  index_midleaves <- floor(length(index_myleaves)/2) # floor(nrow(df_node)/2)
  index_myleaves <- c(index_midleaves:max(index_myleaves),min(index_myleaves):(index_midleaves-1))
  nleaves <- length(index_myleaves)

  df_node$id <- 0
  df_node$id[index_myleaves] <- seq(1:length(index_myleaves))

  df_node$angle <- 90 - 360 * df_node$id / nleaves
  df_node$hjust <- lapply(df_node$angle,function(x){
    if(x < -90 & x != 90) hjust <- 0.5
    if(x >= -90 & x != 90) hjust <- 2
    if(x == 90) hjust <- 0
    hjust
  }) %>% unlist()
  df_node$angle <- ifelse(df_node$angle < -90, df_node$angle+180, df_node$angle)

  # merge
  df_input <- list(df_node = df_node, df_edge = df_edge, df_hierarchy = df_hierarchy)
  return(df_input)

}

drawEdgeBundlingPlot <-function(df_input,colodb){

  # graph

  df_node <- df_input$df_node # node attr
  df_hierarchy <- df_input$df_hierarchy # edge
  mygraph <- igraph::graph_from_data_frame(df_hierarchy, vertices = df_node)

  node_color <- colodb[names(colodb) %in% unique(df_node$label)]

  df_edge <- df_input$df_edge #
  from  <-  match(df_edge$from, df_node$name)
  to  <-  match(df_edge$to, df_node$name)

  ## plot

  pt <- ggraph::ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
    #alpha线的透明度，width线的宽度，tension是线的“密集”程度
    ggraph::geom_conn_bundle(data = ggraph::get_con(from = from, to = to),
                     alpha=1, aes(colour=group, width=score),tension=0.5) +
    # 设置边颜色
    ggraph::scale_edge_color_manual(values = node_color,guide='none')+
    # 设置边粗细
    ggraph::scale_edge_width(range = c(0,1)) +
    # 设置节点标签，字体大小，文本注释信息
    ggraph::geom_node_text(aes(x = x*1.3, y=y*1.3, filter = leaf, label=label,
                       angle = angle, hjust=hjust*0.5, colour=group), size=3, alpha=1, show.legend=F) +
    # 设置节点的大小，颜色和透明度
    ggraph::geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07,
                        colour=group, size=score, alpha=1)) +
    # 设置节点颜色
    scale_colour_manual(values = node_color,name="celltype",
                        breaks=names(node_color)[!names(node_color) %in% c('cell','Sender','Receiver')],) +
    # 设置节点大小的范围
    scale_size_continuous(range = c(1,5)) + theme_void() +
    # 去除alpha图例
    scale_alpha(guide = 'none') +
    theme(
      legend.position = 'bottom',
      # legend.direction = 'vertical',
      legend.box = 'vertical'
    )
  print(pt)

  return(pt)

}

#' @title DrawAlluviumPlot
#' @description Draw Alluvium Plot
#'
#' @param InputDir Character, the path where the signals importance is stored.
#' @param Cluster  Character, signals acitivity of upstream signal pairs of specific cluster are shown.
#' @param ColorDB Vector, containing the color of cell types, key siganls and others.
#' @param Check Logical, Whether to do check visualization.
#' @param top.n Numercial, Only relevant if Check=T. Number of top regulation between LR pairs and target genes to show according to the signals importance.
#' @param gtitle  Character, the title of plot.
#' @param p_height Numercial, the height of plot.
#' @param p_width Numercial, the width of plot.
#'
#' @export
#' @import dplyr ggplot2 grDevices
#' @importFrom stats sd
DrawAlluviumPlot <- function(
  InputDir, Cluster, ColorDB,
  Check = TRUE, top.n = 30,
  gtitle = 'CCI', p_height = 9, p_width = 8
){

  inputdir <- InputDir
  cluster <- Cluster
  colodb <- ColorDB

  outputdir <- './visualize_CCI/AlluvialPlot/'
  dir.create(outputdir,recursive=T,showWarnings = F)

  files <- list.files(inputdir)
  files <- files[grep("_im_",files)]
  files <- files[grep(paste0('-',cluster,'.rds'),files)]
  LRTG_im_merge <- lapply(files, function(f){

    LRTG_im <- readRDS(paste0(inputdir,f))
    LRTG_im$Sender <- strsplit(f,'-|_')[[1]][4]
    LRTG_im
    # head(LRTG_im)

  }) %>% do.call('rbind',.)

  ## plot

  df_MLnet_long_check <- prepareAlluviumPlotData(lrtg_im = LRTG_im_merge, color.by = 'Sender',
                                                 do.check = Check, top.n = top.n)

  pt <- drawAlluviumPlot(df_input = df_MLnet_long_check, colodb = colodb)

  ## save

  pdf(paste0(outputdir,'AlluviumPlot-',gtitle,'.pdf'),height = p_height,width = p_width)
  print(pt)
  dev.off()

  png(paste0(outputdir,'AlluviumPlot-',gtitle,'.png'),height = p_height,width = p_width, units = 'in', res = 300)
  print(pt)
  dev.off()

  return(pt)

}

prepareAlluviumPlotData <- function(lrtg_im = NULL, color.by = 'Nodekey', do.check = TRUE,top.n = 30){

  df_lrtg_im = lrtg_im
  colnames(df_lrtg_im)[grep('im_norm',colnames(df_lrtg_im))] <- 'Score'
  df_lrtg_im$ID <- seq(nrow(df_lrtg_im))
  df_plot <- df_lrtg_im

  if(do.check & nrow(df_plot)>top.n){

    if(stats::sd(df_plot$Score)!=0){

      message('using the score of LRTG to narrow down for better visualization')
      df_plot_check <- df_plot
      df_plot_check$rank <- df_plot_check$Score %>% rank()
      df_plot_check <- df_plot_check[df_plot_check$rank >= max(df_plot_check$rank)-top.n+1,]
      df_plot <- df_plot_check

    }else{

      message('using the number of Target to narrow down for better visualization')
      df_plot_check <- df_plot
      df_plot_check <- df_plot_check %>% group_by(Ligand,Receptor,TF) %>%
        summarise(count=n()) %>% arrange(desc(count)) %>% ungroup() %>%
        inner_join(df_plot_check,., by = c('Ligand','Receptor','TF'))
      df_plot_check$rank <- df_plot_check$count %>% rank()
      df_plot_check <- df_plot_check[df_plot_check$rank >= max(df_plot_check$rank)-top.n+1,]
      df_plot <- df_plot_check

    }

  }

  df_plot <- df_plot[,-grep('LRpair|IM|rank',colnames(df_plot))]
  df_plot_long <- reshape2::melt(df_plot, id = c('Score','ID','Sender'))
  colnames(df_plot_long) <- c("Score",'ID','Sender',"Nodekey","Node")

  if(color.by == 'Nodekey'){
    df_plot_long$color <- df_plot_long$Nodekey
  }else if(color.by == 'Node'){
    df_plot_long$color <- df_plot_long$Node
  }else{
    df_plot_long$color <- df_plot[[color.by]][match(df_plot_long$ID,df_plot$ID)]
  }

  return(df_plot_long)

}

drawAlluviumPlot <- function(df_input, colodb){

  pt_color <- colodb
  if(!any(unique(df_input$color) %in% names(pt_color))){
    pt_color <- rainbow(length(unique(df_input$color)))
    names(pt_color) <- unique(df_input$color)
  }else{
    pt_color <- pt_color[names(pt_color) %in% unique(df_input$color)]
  }

  pt <- ggplot(data = df_input,
               aes(x = Nodekey,y = Score,
                   stratum = Node, alluvium = ID,label = Node)) +
    ggalluvial::geom_flow(aes(fill = color)) + ggalluvial::geom_stratum() +
    theme_minimal() + geom_text(stat = "stratum", size = 3,check_overlap = F) +
    scale_fill_manual(values = pt_color,  name="color",breaks=names(pt_color)) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 11),
      plot.title = element_text(size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      panel.grid = element_blank(),
      legend.position = 'none'
    )
  print(pt)

  return(pt)
}

#-----------------------------------------------------------------update 2024-03-13 -------------------------------------------
#' @title DrawHeatmapPlot
#' @description Draw Heatmap Plot
#'
#' @param InputDir Character, the path where the stMLnet's result is stored.
#' @param Sender Character, the sender cell type, default = NULL.
#' @param Receiver Character, the receiver cell type, default = NULL.
#' @param outputdir Character, the path to save plot.
#' @param p_height Numercial, the height of plot.
#' @param p_width Numercial, the width of plot.
#'
#' @export
#' @import dplyr ggplot2 utils grDevices
#' @importFrom stats na.omit
#' @importFrom BiocGenerics toTable
# Heatmap plot pval
DrawHeatmapPlot <- function(InputDir, Sender = NULL, Receiver = NULL, outputdir,p_height = 9,p_width = 8.4){

  inputdir <- InputDir
  sender <- Sender
  receiver <- Receiver

  wd <- paste0(InputDir,"/runscMLnet/")
  files <- list.files(wd)[grep('TME_',list.files(wd),invert = TRUE)]

  LRI_allpval <- lapply(files, function(f) {
    cat("Processing file:", f, "\n")

    LRI_pval <- readRDS(paste0(wd, f, "/cellpair_LRI_pval.rds"))
    print(str(LRI_pval))  # Print the structure of LRI_pval

    if (length(LRI_pval) > 0) {
      LRI_pval$Sender <- strsplit(f, "_")[[1]][1]
      LRI_pval$Receiver <- strsplit(f, "_")[[1]][2]
      return(LRI_pval)
    } else {
      # If LRI_pval is empty, return a placeholder or handle it as needed
      cat("Warning: Empty LRI_pval for file", f, "\n")
      return(NULL)  # or return an empty data frame, depending on your needs
    }
  }) %>% do.call('rbind', .)

  if (length(sender) != 0 ){

    LRI_pval_Inter <- LRI_allpval[LRI_allpval$Sender == sender,]
    df_plot <- data.frame(cellpair = paste0(LRI_pval_Inter$Sender,"_",LRI_pval_Inter$Receiver),
                          LRpair = paste0(LRI_pval_Inter$source,"_",LRI_pval_Inter$target),
                          pval = LRI_pval_Inter$pval)

    df_plot <- df_plot[order(df_plot$pval,decreasing=F),]

    # load LR signaling score
    ct <- sender
    workpath <- paste0(inputdir,'/runModel/')
    files = list.files(workpath)
    files = files[grep(paste0("LRTG_allscore_",ct),files)]

    df_LRTGscore = lapply(files, function(file){

      print(file)
      LRS_score = readRDS(paste0(workpath,file))[[1]]
      LRS_score_merge = do.call('cbind',LRS_score)
      if (length(unique(colnames(LRS_score_merge))) == 1){
        LRpair <- unique(colnames(LRS_score_merge))
        LRS_score_merge = LRS_score_merge[,1]

        # file <- gsub('-','_',file)
        df_LigRec <- data.frame(
          source = LRpair %>% gsub('_.*','',.),
          target = LRpair %>% gsub('.*_','',.),
          LRpair = LRpair,
          count = mean(LRS_score_merge),
          source_group = strsplit(file,'[_\\.]')[[1]][3],
          target_group = strsplit(file,'[_\\.]')[[1]][4])

      }else{
        LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]

        # file <- gsub('-','_',file)
        df_LigRec <- data.frame(
          source = colnames(LRS_score_merge) %>% gsub('_.*','',.),
          target = colnames(LRS_score_merge) %>% gsub('.*_','',.),
          LRpair = colnames(LRS_score_merge),
          count = colMeans(LRS_score_merge),
          source_group = strsplit(file,'[_\\.]')[[1]][3],
          target_group = strsplit(file,'[_\\.]')[[1]][4])
      }

    }) %>% do.call('rbind',.)

    df_LRTGscore$pval <- NA
    for (lr in df_LRTGscore$LRpair){
      pos1 <- which(df_plot$LRpair %in% lr)
      pos2 <- which(df_LRTGscore$LRpair %in% lr)
      df_LRTGscore$pval[pos2] <- df_plot$pval[pos1]
    }

    df <- data.frame(cellpair = paste0(df_LRTGscore$source_group,"_",df_LRTGscore$target_group),
                     LRpair = df_LRTGscore$LRpair,
                     pval = df_LRTGscore$pval,
                     count = df_LRTGscore$count)

  }

  if(length(receiver) != 0 ){

    LRI_pval_Inter <- LRI_allpval[LRI_allpval$Receiver == receiver,]
    df_plot <- data.frame(cellpair = paste0(LRI_pval_Inter$Sender,"_",LRI_pval_Inter$Receiver),
                          LRpair = paste0(LRI_pval_Inter$source,"_",LRI_pval_Inter$target),
                          pval = LRI_pval_Inter$pval)

    df_plot <- df_plot[order(df_plot$pval,decreasing=F),]

    # load LR signaling score
    ct <- receiver
    workpath <- paste0(inputdir,'/runModel/')
    files = list.files(workpath)
    files = files[grep(paste0(ct,".rds"),files)]

    df_LRTGscore = lapply(files, function(file){

      print(file)
      LRS_score = readRDS(paste0(workpath,file))[[1]]
      LRS_score_merge = do.call('cbind',LRS_score)
      if (length(unique(colnames(LRS_score_merge))) == 1){
        LRpair <- unique(colnames(LRS_score_merge))
        LRS_score_merge = LRS_score_merge[,1]

        # file <- gsub('-','_',file)
        df_LigRec <- data.frame(
          source = LRpair %>% gsub('_.*','',.),
          target = LRpair %>% gsub('.*_','',.),
          LRpair = LRpair,
          count = mean(LRS_score_merge),
          source_group = strsplit(file,'[_\\.]')[[1]][3],
          target_group = strsplit(file,'[_\\.]')[[1]][4])

      }else{
        LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]

        # file <- gsub('-','_',file)
        df_LigRec <- data.frame(
          source = colnames(LRS_score_merge) %>% gsub('_.*','',.),
          target = colnames(LRS_score_merge) %>% gsub('.*_','',.),
          LRpair = colnames(LRS_score_merge),
          count = colMeans(LRS_score_merge),
          source_group = strsplit(file,'[_\\.]')[[1]][3],
          target_group = strsplit(file,'[_\\.]')[[1]][4])
      }

    }) %>% do.call('rbind',.)

    df_LRTGscore$pval <- NA
    for (lr in df_LRTGscore$LRpair){
      pos1 <- which(df_plot$LRpair %in% lr)
      pos2 <- which(df_LRTGscore$LRpair %in% lr)
      df_LRTGscore$pval[pos2] <- df_plot$pval[pos1]
    }

    df <- data.frame(cellpair = paste0(df_LRTGscore$source_group,"_",df_LRTGscore$target_group),
                     LRpair = df_LRTGscore$LRpair,
                     pval = df_LRTGscore$pval,
                     count = df_LRTGscore$count)

  }

  p1 <- ggplot(df, aes(x = cellpair, y = LRpair, color = pval, size = count)) +
    geom_point(pch = 16) +
    scale_color_gradient(low = "red", high = "yellow")+
    #theme_linedraw() +
    #theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, size = 10,hjust= NULL, vjust = NULL),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")
  print(p1)

  ## save

  pdf(paste0(plotdir,ct,'_bubble_pval.pdf'),height = p_height,width = p_width)
  print(p1)
  dev.off()

}
#-----------------------------------------------------------------update 2024-03-13 -------------------------------------------

DrawNetworkPlot_v2 <- function(InputDir, Metric, ColorDB, gtitle = "CCI"){
  suppressMessages(library(CellChat))
  inputdir <- InputDir
  key <- Metric

  outputdir <- "./visualize_CCI/NetworkPlot/"
  dir.create(outputdir, recursive = T, showWarnings = F)
  files <- list.files(inputdir)[grep("LRTG_im_clean_", list.files(inputdir))]
  
  LRTG_detail <- lapply(files, function(f) {
    LRTG_im <- readRDS(paste0(inputdir, f))
    c(length(unique(LRTG_im$LRpair)), length(unique(LRTG_im$Target)), 
      sum(LRTG_im$IM), sum(LRTG_im$im_norm), sum(LRTG_im$IM)/nrow(LRTG_im), 
      sum(LRTG_im$im_norm)/nrow(LRTG_im))
  }) %>% do.call("rbind", .) %>% as.data.frame()
  
  df_cellpair <- gsub("LRTG_im_clean_|.rds", "", files) %>% 
    strsplit(., "-") %>% do.call("rbind", .) %>% as.data.frame()
  
  LRTG_detail <- cbind(df_cellpair, LRTG_detail)
  LRTG_detail <- na.omit(LRTG_detail)
  colnames(LRTG_detail) <- c("cell_from", "cell_to", "n_LRs", 
                             "n_TGs", "IM", "IM_norm", "mean_IM", "mean_IM_norm")
  
  celltype = unique(c(LRTG_detail$cell_from,LRTG_detail$cell_to))
  LRTG_detail[[key]] <- as.numeric(LRTG_detail[[key]])
  
  tmeTab <- LRTG_detail[,c('cell_from','cell_to',key)] 
  
  mat2 <- matrix(0,nrow = length(celltype),ncol = length(celltype))
  rownames(mat2) <- celltype
  colnames(mat2) <- celltype
  for (i in 1:nrow(tmeTab)){
    c1 <- tmeTab$cell_from[i]
    c2 <- tmeTab$cell_to[i]
    val <- tmeTab$n_LRs[i]
    mat2[c1,c2] <- val

  }
  
  colordb <- ColorDB[rownames(mat2)]
  
  pdf(paste0(outputdir,gtitle,"_NetworkPlot",".pdf"),height = 6,width =6)
  netVisual_circle(mat2, color.use = colordb,vertex.weight = rowSums(mat2),
                   alpha.edge = 0.8,edge.width.max= 3,
                   weight.scale = T, label.edge= F, title.name = "Number of interactions",
                   arrow.width = 0.8,arrow.size = 0.8,
                   text.x = 10,text.y = 1.5)
  dev.off()
  
}

DrawCircosPlot <- function(InputDir, receiver, ColorDB){
  
  suppressMessages(library(iTALK))
  suppressMessages(library(circlize))
  
  inputdir = InputDir
  
  outputdir <- "./visualize_CCI/CircosPlot/"
  dir.create(outputdir, recursive = T, showWarnings = F)
  
  files = list.files(inputdir)
  files = files[grep(paste0('-',receiver,'.rds'),files)]
  
  df_LRTGscore = lapply(files, function(file){
    
    print(file)
    LRS_score = readRDS(paste0(inputdir,file))[[1]]
    LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]
    
    # file <- gsub('-','_',file)
    df_LigRec <- data.frame(
      source = colnames(LRS_score_merge) %>% gsub('_.*','',.),
      target = colnames(LRS_score_merge) %>% gsub('.*_','',.),
      LRpair = colnames(LRS_score_merge),
      count = colMeans(LRS_score_merge),
      source_group = strsplit(file,"[._-]")[[1]][3],
      target_group = strsplit(file,"[._-]")[[1]][4]
    )
    
  }) %>% do.call('rbind',.)
  
  if(!is.null(df_LRTGscore)){
    
    # input
    res <- data.frame(ligand = df_LRTGscore$source,receptor = df_LRTGscore$target,
                      cell_from = df_LRTGscore$source_group,cell_to = df_LRTGscore$target_group, 
                      count = df_LRTGscore$count, comm_type = "growth factor")
    res <- res[order(res$count,decreasing=T),] 
    # plot
    if (dim(res)[1] > 30){
      # select top 20
      res <- res[1:30,]
    }
    
    min(res$count)
    
    colordb <- ColorDB[which(names(ColorDB) %in% c(unique(res$cell_from),unique(res$cell_to)))]
    
    pdf(paste0(outputdir,"ChordPlot_v2_LRscore_","receiver_",receiver,".pdf"),height = 5,width = 5)
    LRPlot(res,datatype='mean count',
           cell_col=colordb,
           link.arr.lwd=res$count,
           link.arr.col="#696969", # "#808080"
           link.arr.width=0.25,
           track.height_1 = uh(1,"mm"),
           track.height_2 = uh(11,"mm"),
           annotation.height_1 = 0.015,
           annotation.height_2 = 0.01,
           text.vjust = "0.4cm")
    dev.off()
  }
    
}

DrawAlluviumPlot_v2 <- function(MLnetDir = MLnetDir,ImportDir= ImportDir, Sender,Receiver, ColorDB = colordb,
                                Check = TRUE, top.n = 30, p_height = 7.5, p_width = 7){
  suppressMessages(library(ggsci))
  mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)
  
  nodekey <- c("Ligand","Receptor","TF","Target")
  mycolor_key <- mycolors_locus[1:4]
  names(mycolor_key) <- nodekey
  
  mlnetdir <- MLnetDir
  importdir <- ImportDir
  sender <- Sender
  receiver <- Receiver
  colodb <- ColorDB
  outputdir <- "./visualize_CCI/AlluvialPlot/"
  dir.create(outputdir, recursive = T, showWarnings = F)
  
  cp = paste0(sender,"-",receiver)
  
  # load MLnet
  
  files <- list.files(mlnetdir)
  MLnet <- readRDS(paste0(mlnetdir,sender,"_",receiver,"/scMLnet.rds"))
  
  df_ligrec = data.frame(Ligand = MLnet$LigRec$source, Receptor = MLnet$LigRec$target) 
  df_rectf = data.frame(Receptor = MLnet$RecTF$source, TF = MLnet$RecTF$target)
  df_tftar = data.frame(TF = MLnet$TFTar$source, Target = MLnet$TFTar$target)
  
  df_mlnet = df_ligrec %>% merge(., df_rectf, by = 'Receptor') %>% 
    merge(., df_tftar, by = 'TF') %>% 
    dplyr::select(Ligand, Receptor, TF, Target) %>% 
    arrange(Ligand, Receptor)
  df_mlnet$LRpair <- paste0(df_mlnet$Ligand,"_",df_mlnet$Receptor)
  
  # load LRTG_score
  files = list.files(importdir)
  files = files[grep("_im_",files)]
  file_cp = files[grep(paste0(cp,".rds"),files)]
  
  LRTG_im_merge <- lapply(file_cp, function(f){
    
    LRTG_im <- readRDS(paste0(importdir,f))
    LRTG_im$Sender <- strsplit(f,'-|_')[[1]][4]
    LRTG_im
    # head(LRTG_im)
    
  }) %>% do.call('rbind',.)
  
  LRTG_im_merge$LRTG <- paste0(LRTG_im_merge$LRpair,"_",LRTG_im_merge$Target)
  df_mlnet$LRTG <- paste0(df_mlnet$LRpair,"_",df_mlnet$Target)
  LRTG_im_merge$TF <- df_mlnet$TF[which(LRTG_im_merge$LRTG %in% df_mlnet$LRTG)]
  
  LRTG_im_merge <- LRTG_im_merge[,-8]
  
  df_MLnet_long <- prepareAlluviumPlotData(lrtg_im = LRTG_im_merge, 
                                            color.by = 'Sender',
                                            do.check = TRUE,top.n=30)
  df_MLnet_long$Nodekey <- factor(df_MLnet_long$Nodekey,levels = c("Ligand", "Receptor","TF", "Target" ))
  
  pt <-ggplot(df_MLnet_long,
              aes(x = Nodekey, stratum = Node, alluvium = ID,
                  y = Score, fill = Nodekey,label = Node)) +
    scale_x_discrete(expand = c(.1, .1)) +
    scale_fill_manual(values = mycolor_key) + 
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3) +
    theme_minimal()+
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 11),
      plot.title = element_text(size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      panel.grid = element_blank(),
      legend.position = 'none'
    ) +
    ggtitle("Multilayer Network")
  pt
  
  pdf(paste0(outputdir,'AlluviumPlot-LRTFTG_',cp,'.pdf'),height = p_height,width = p_width)
  print(pt)
  dev.off()
  
}




