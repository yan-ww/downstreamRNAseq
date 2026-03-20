
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(DOSE))

GOplot = function(go.res,color) {
  goBP <- subset(go.res, subset = (ONTOLOGY == "BP"))
  if(nrow(goBP)>10){
    goBP = goBP[1:10,]
  }
  goCC <- subset(go.res, subset = (ONTOLOGY == "CC"))
  if(nrow(goCC)>10){
    goCC = goCC[1:10,]
  }
  goMF <- subset(go.res, subset = (ONTOLOGY == "MF"))
  if(nrow(goMF)>10){
    goMF = goMF[1:10,]
  }
  term <- rbind(goBP, goCC, goMF)
  term$Description=factor(term$Description, levels = term$Description)
  gobar = ggplot(data=term, aes(x=Description, y=Count,color=ONTOLOGY)) + geom_bar(stat="identity",width=0.8,aes(fill=ONTOLOGY)) +
    coord_flip()+ xlab("GO term") + ylab("Num of Genes") + theme_bw() + labs(title = sprintf("GO Enrichment of %s",color)) +
    theme(text=element_text(size=20))+scale_x_discrete(labels = function(x) str_wrap(x,width=50))
  #ggsave(gobar,filename = sprintf("GOafWGCNA/GOplotOf%s.png",color),width=500,height=400,units="mm")
  ggsave(gobar,filename = sprintf("GOafWGCNA/GOplotOf%s.pdf",color),width=500,height=400,units="mm")
}

KEGGplot = function(kegg,color){
  if(nrow(kegg)>15){
    kegg = kegg[1:15,]
  }
  kegg$Description = factor(kegg$Description,levels = rev(kegg$Description))
  kegg <- mutate(kegg, ratio = parse_ratio(GeneRatio))
  p = ggplot(data = kegg,aes(x = ratio, y = reorder(Description,Count)))+
    geom_point(aes(size = Count,color = -log10(p.adjust)))+
    theme_bw()+
    scale_colour_gradient(low = "green",high = "red")+
    scale_y_discrete(labels = function(x) str_wrap(x,width = 40))+
    labs(x = "GeneRatio",y = "",title = sprintf("KEGG pathway of %s",color),
         color = expression(-log10(p.adjust)),size = "Count")+
    theme(text=element_text(size=20))

}

