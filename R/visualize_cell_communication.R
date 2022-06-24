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
  files <- files[grep(paste0('_',cluster,'.rds'),files)]
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
      source_group = strsplit(file,'[_\\.]')[[1]][3],
      target_group = strsplit(file,'[_\\.]')[[1]][4]
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
  df_edge <- LRTGscore[,c('from','to')]

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

#' @title DrawEnrichmentPlot
#' @description Draw Enrichment Plot
#'
#' @param InputDir Character, the path where the multilayer network is stored.
#' @param Cluster Character, enrichment analysis of receptors's downstream target genes in specific cluster are shown.
#' @param top.n Numercial, Only relevant if Check=T. Number of top enrichment term to show according to the p value.
#' @param gtitle  Character, the title of plot.
#' @param p_height Numercial, the height of plot.
#' @param p_width Numercial, the width of plot.
#'
#' @export
#' @import dplyr ggplot2 utils grDevices
#' @importFrom stats na.omit
#' @importFrom BiocGenerics toTable
DrawEnrichmentPlot <- function(InputDir, Cluster, top.n = 3, gtitle = 'Enrichment', p_height = 7.5, p_width = 7){

  inputdir <- InputDir
  cluster <- Cluster
  outputdir <- './visualize_CCI/Heatmap/'
  dir.create(outputdir,recursive=T,showWarnings = F)

  ## input

  files <- list.files(inputdir)
  files <- files[grep(paste0('_',cluster,'$'),files)]
  mulNetDF = lapply(files, function(f){

    mlnet <- readRDS(paste0(inputdir,f,"/scMLnet.rds"))
    if(nrow(mlnet$LigRec)!=0){
      colnames(mlnet$LigRec) <- c('Ligand','Receptor')
      colnames(mlnet$RecTF) <- c('Receptor','TF')
      colnames(mlnet$TFTar) <- c('TF','Target')
      mlnet_df <- merge(mlnet$LigRec,mlnet$RecTF,by='Receptor')
      mlnet_df <- merge(mlnet_df,mlnet$TFTar,by='TF')
    }else{
      mlnet_df <- as.data.frame(matrix(ncol=4,dimnames = list(c(),c('TF', 'Receptor', 'Ligand', 'Target'))))
    }
    mlnet_df[,c('Receptor','Target')]

  }) %>% do.call('rbind',.)
  mulNetDF = stats::na.omit(mulNetDF)
  mulNetList = split(mulNetDF,mulNetDF$Receptor)

  ## enrich

  files <- list.files(outputdir)
  check_res <- grep(paste0('res_Enrich_',gtitle,'.rds'),files)
  if(length(check_res)==1){

    res_Enrich <- readRDS(paste0(outputdir,'res_Enrich_',gtitle,'.rds'))

  }else{

    res_GO <- Perf_Enrich(mulNetList,DB='GO')
    res_KEGG <- Perf_Enrich(mulNetList,DB='KEGG')
    res_Enrich <- list(res_GO = res_GO,
                       res_KEGG = res_KEGG)
    saveRDS(res_Enrich, paste0(outputdir,'res_Enrich_',gtitle,'.rds'))

  }

  ## ORA-KEGG

  res_KEGG <- res_Enrich$res_KEGG
  df_kegg <- res_KEGG[res_KEGG$pvalue <= 0.01,]
  if(nrow(df_kegg)>0){

    df_kegg$geneRatio <- lapply(df_kegg$GeneRatio,function(chr){
      x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
      y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
      x/y
    }) %>% unlist()
    df_kegg$bgRatio <- lapply(df_kegg$BgRatio,function(chr){
      x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
      y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
      x/y
    }) %>% unlist()

    keep_term <- lapply(unique(df_kegg$Regulator), function(key){

      df_kegg[df_kegg$Regulator==key,'ID'] %>% head(.,top.n)

    }) %>% unlist() %>% unique()
    df_kegg$ONTOLOGY <- 'KEGG'
    df_kegg <- df_kegg[df_kegg$ID %in% keep_term,c(13,1:2,5:7,10:11)]
    df_kegg <- df_kegg[order(df_kegg$ONTOLOGY,df_kegg$ID),]

  }

  ## ORA-GO

  res_GO <- res_Enrich$res_GO
  df_go <- res_GO[res_GO$pvalue<=0.01,]
  if(nrow(df_go)>0){

    df_go$geneRatio <- lapply(df_go$GeneRatio,function(chr){

      x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
      y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
      x/y
    }) %>% unlist()
    df_go$bgRatio <- lapply(df_go$BgRatio,function(chr){
      x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
      y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
      x/y
    }) %>% unlist()
    df_go <- df_go[df_go$ONTOLOGY=='BP',]

    keep_term <- lapply(unique(df_go$Regulator), function(key){

      df_go[df_go$Regulator==key,'ID'] %>% head(.,top.n)

    }) %>% unlist() %>% unique()
    df_go <- df_go[df_go$ID %in% keep_term,c(1:3,6:8,11:12),]
    df_go <- df_go[order(df_go$ONTOLOGY,df_go$ID),]

  }

  ## data

  df_plot <- rbind(df_go,df_kegg)
  df_plot <- df_plot[df_plot$ONTOLOGY %in% c('KEGG','BP'),]
  df_plot$Description <- factor(df_plot$Description, levels = unique(df_plot$Description))
  df_plot$ONTOLOGY <- factor(df_plot$ONTOLOGY, levels = c('BP','MF','CC','KEGG'))
  anno_x_loca <- lapply(unique(df_plot$ONTOLOGY), function(key){

    df <- df_plot[!duplicated(df_plot$ID),]
    grep(key,df$ONTOLOGY) %>% max()

  }) %>% unlist()
  anno_x_loca <- anno_x_loca+0.5
  names(anno_x_loca) <- unique(df_plot$ONTOLOGY)

  ## plot

  if(nrow(df_plot)>0){

    pt_merge <- ggplot(data=df_plot,aes(x=Description,y=Regulator,size=geneRatio,col=p.adjust))+
      geom_hline(aes(x=Description,y=Regulator,yintercept = 1:nrow(df_plot)),size= 1.5,colour= "#E4EDF2",alpha= .5)+
      geom_vline(aes(x=Description,y=Regulator,xintercept = anno_x_loca[1]),size=0.5,linetype= "dashed")+
      geom_point()+ coord_flip() +
      scale_color_gradient(low = '#C51162',high = '#FCE4EC') +
      scale_size_continuous(range = c(1,4)) +
      theme_bw() + labs(title='Function Enrichment Analysis of CCI',x='',y='') +
      theme(
        plot.title = element_text(hjust = 0.5,size = 12), # 标题居中
        panel.background = element_blank(), # 去除坐标图的背景色
        legend.key = element_blank(), # 去除图例图案的背景色
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), # x轴刻度方向
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank(), # 去除辅助网格线
        # legend.position = 'bottom',
        # legend.direction = 'horizontal'
        legend.position = 'right',
        legend.direction = 'vertical'
      )
    print(pt_merge)

    png(paste0(outputdir,"Enrichment_BP_KEGG_",gtitle,".png"),
        width = 15, height = 5, units = 'in', res = 1000)
    print(pt_merge)
    dev.off()

    pdf(paste0(outputdir,"Enrichment_BP_KEGG_",gtitle,".pdf"),
        width = 15, height = 5)
    print(pt_merge)
    dev.off()

  }

  return(pt)

}

Perf_ORA <- function(ls_tg,DB=c('GO','KEGG')){

  geneList <- ls_tg
  g2s <- toTable(org.Hs.egSYMBOL)
  geneList <- g2s$gene_id[match(geneList,g2s$symbol)]
  geneList <- stats::na.omit(geneList)

  tryCatch({
    if(DB == 'GO'){

      res_enrich <- clusterProfiler::enrichGO(geneList, 'org.Hs.eg.db', ont = 'ALL', keyType = 'ENTREZID',
                             minGSSize = 1, pvalueCutoff = 0.99)@result
      res_enrich$GeneRatio <- res_enrich$Count/length(enrich_gobp@gene)
      res_enrich <- res_enrich[,c("ID","Description","GeneRatio","p.adjust",'ONTOLOGY')]

    }else{

      res_enrich <- clusterProfiler::enrichKEGG(geneList, minGSSize = 1, pvalueCutoff = 0.99)@result
      res_enrich$GeneRatio <- res_enrich$Count/length(enrich_kegg@gene)
      res_enrich <- res_enrich[,c("ID","Description","GeneRatio","p.adjust")]
      res_enrich$ONTOLOGY <- 'KEGG'

    }
  }, error = function(e){
    res_enrich <- as.data.frame(matrix(
      ncol = 5,dimnames = list(c(),c("ID","Description","GeneRatio","p.adjust",'ONTOLOGY'))))
  })

  if(!exists('res_enrich')){
    res_enrich <- as.data.frame(matrix(
      ncol = 5,dimnames = list(c(),c("ID","Description","GeneRatio","p.adjust",'ONTOLOGY'))))
  }

  return(res_enrich)
}

Perf_Enrich <- function(ls_im, DB=c('GO','KEGG')){

  Recs <- names(ls_im)
  res_enrich <- lapply(1:length(Recs), function(i){

    # rec <- 'ACVR2A'
    print(paste0(i,"/",length(Recs)))
    tgs <- ls_im[[i]]$Target

    res_ORA <- Perf_ORA(tgs,DB=DB)
    if(nrow(res_ORA)==0){
      res_ORA <- as.data.frame(matrix(
        ncol = 6,dimnames = list(c(),
                                 c("ID","Description","GeneRatio","p.adjust","ONTOLOGY","Regulator"))))
      res_ORA
    }else{
      res_ORA$Regulator <- Recs[i]
      res_ORA
    }

  }) %>% do.call('rbind',.)

  return(res_enrich)
}
