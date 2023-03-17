library(Seurat)
library(Rtsne)
library(ggplot2)
library(magrittr) 
library(dplyr)  

# Draw scFseCluster scatter plot
# data: feature subset
# label: cell label
# point_size = 2: Scatter plot point size, default equal to 2
# method = "TSNE": Using the method of dimensionality reduction, the default is to use TSNE, but you can also use UMAP
# Trans = TRUE: Whether to transpose, default transpose
# isSave = TRUE: Whether to save the scatter plot image, default save
plot_scFseCluster <- function(data, label ,point_size=2, method="TSNE" ,Trans=TRUE, isSave=TRUE) {

  zedata <- data
  if (Trans) {
    zedata = t(zedata) 
  }
  print(dim(zedata))
  
  ze <- CreateSeuratObject(counts = zedata, project = "plot")
  ze <- NormalizeData(ze, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(ze)
  ze <- ScaleData(ze, features = all.genes)
  ze <- RunPCA(ze, features = all.genes)
  
  pc1 = NULL
  
  if (method == "TSNE") {
    ze <- RunTSNE(ze, dims = 1:10)
    tsne <- ze@reductions[["tsne"]]@cell.embeddings
    pc1 <- data.frame(var_x=tsne[,1],var_y=tsne[,2])
  } else if (method == "UMAP") {
    ze <- RunUMAP(ze, dims = 1:10)
    umap <- ze@reductions[["umap"]]@cell.embeddings
    pc1 <- data.frame(var_x=umap[,1], var_y=umap[,2])
  }
  
  class.label <- label
  class.label <- as.matrix(class.label)
  class.label <- class.label[,2]
  cluster <- class.label
  cluster <- as.factor(cluster)
  
  Color <- c("#FD2516", "#04A108", "#1452FF", "#FFD700", "#FF61CC", "#55BCC3", "#FC7203", "#C1D245", "#CF7CB0")

  p <- ggplot(data=pc1, aes(x=var_x, y=var_y,color=cluster)) + 
    scale_color_manual(values = Color) +
    geom_point(size=1) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme_classic() + 
    theme(legend.position="none") + 
    labs(x =NULL,y = NULL) + 
    theme(legend.position="none") + 
    geom_point(size = point_size) + 
    theme(axis.text=element_text(size=20),axis.title = element_text(face = "bold")) + 
    theme(legend.title = element_blank())
  
  if (isSave) {
    if (method == "UMAP") {
      save_path = paste("./result/scFseCluster_UMAP.tiff",sep="")
      ggsave(save_path,p,dpi = 300) 
    } else if (method == "TSNE") {
      save_path = paste("./result/scFseCluster_TSNE.tiff",sep="")
      ggsave(save_path,p,dpi = 300)
    }
  }
  
  return(p)
}

args = commandArgs(T)
data_path <- args[1]
label_path <- args[2]

zedata <- read.table(data_path,stringsAsFactors=FALSE, header=TRUE, check.names=FALSE,row.names=1, sep=",")
label <- read.table(label_path, header=T,sep=",",check.names=F)

plot_scFseCluster(zedata,label)