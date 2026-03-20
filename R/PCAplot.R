#' @title PCA分析及绘图
#' @description 对RNA-seq的count数据进行PCA分析，通过DEseq2标准化
#' @param data 输入数据框或矩阵
#' @param coldata 分组信息
#' @return PCA图片
#' @export
#' @examples
#' \dontrun{
#' normalized_data <- normalize_data(my_data, method = "CPM")
#' }


PCAplot <- function(datacount,coldata){
  datacount<-count
  summary(is.numeric(datacount[,1]))
  dds <- DESeq2::DESeqDataSetFromMatrix(datacount, coldata, design = ~group)
  dds <- DESeq2::DESeq(dds)
  vsd<-DESeq2::varianceStabilizingTransformation(dds,blind = FALSE)
  #这里是用DEseq的标准化方法，也可以直接用fpkm
  vstmat<-assay(vsd)
  pca<-FactoMineR::PCA(t(vstmat),scale.unit = FALSE,ncp = 5,graph = FALSE)
  factoextra::fviz_pca_ind(pca, geom.ind  = c("point","text"),
             col.ind = condition, # color by groups
             pointsize = 1.5, pointshape = 16,
             addEllipses = T, repel = F ,legend.title = "Tissue",
             title="PCA analysis")+theme(plot.title = element_text(hjust = 0.5))
}

