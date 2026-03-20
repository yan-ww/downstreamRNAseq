
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
options(download.file.method="wget")
DE<-function(a){
  b<-colnames(a)
  condition<-substring(b,nchar(b))
  people<-substr(b,3,nchar(b)-1)
  coldata<-data.frame(row.names=colnames(a),people,condition)
  dds<-DESeqDataSetFromMatrix(a,coldata,design=~people+condition)
  dds<-DESeq(dds)
  res<-results(dds,contrast=c("condition","F","N"))
  res<-res[order(res$pvalue),]
  res$change<-as.factor(ifelse((res$pvalue) < 0.05 & (res$log2FoldChange) > log2(1.5),"up",ifelse(res$pvalue < 0.05 & res$log2FoldChange < (-log2(1.5)),"down","not")))
  return(res)
}


volcano<-function(de.res){
  c<-de.res[c("pvalue","log2FoldChange","change")]
  c$gene<-rownames(de.res)
  c$label<-""
  upgenes<-head(c$gene[which(c$change=="up")],10)
  downgenes<-head(c$gene[which(c$change=="down")],10)
  top10genes<-c(as.character(upgenes),as.character(downgenes))
  c$label[match(top10genes,c$gene)]<-top10genes
  c<-as.data.frame(c)
  volcanoplot<-ggplot(data=c, aes(x=log2FoldChange, y=-log10 (pvalue),color=change)) + geom_point(alpha=0.5, size=1.75) +
    theme_set(theme_bw(base_size=20))+  scale_color_manual(values=c('#2f5688','#BBBBBB','#CC0000'))+
    xlab("log[2] fold change(F/N)") + ylab("-log[10] pvalue") + theme(text=element_text(size=20),plot.title = element_text(size=20,hjust = 0.5),
                                                                      legend.position="right",legend.title = element_blank()) + geom_hline(yintercept = 1.30,linetype="dashed")+
    geom_vline(xintercept = c(-0.6,0.6),linetype="dashed")+ geom_text(aes(label = label),size = 5,color = "black",show.legend = FALSE )
  #改动如下：y的值从pvalue改成padj，geom_vline的xintercept从(-0.6,0.6)改成c(-1,1)
  return(volcanoplot)
}

double_enrich_plot = function(df,cate){
  df$pl = ifelse(df$category == "up",-log10(df$p.adjust),log10(df$p.adjust))
  df = arrange(df,category,pl)
  df$Description = factor(df$Description,levels = unique(df$Description),ordered = TRUE)
  tmp = with(df, labeling::extended(range(pl)[1], range(pl)[2], m = 10))
  lm = tmp[c(1,length(tmp))]
  lm = c(floor(min(df$pl)),ceiling(max(df$pl)))
  ggplot(df, aes(x=Description, y= pl)) +labs(y="p.adjust")+
    geom_bar(stat='identity', aes(fill=category), width=.7)+
    scale_fill_manual(values = c("#2874C5", "#f87669"))+
    coord_flip()+theme_light() +ylim(lm)+
    scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
    scale_y_continuous(breaks = tmp,labels = abs(tmp))+
    theme(text=element_text(size=20),panel.border = element_blank())
  ggsave(filename = sprintf("%sofDEG.pdf",cate))
  ggsave(filename = sprintf("%sofDEG.png",cate),width=500,height=400,units="mm")
}

utils::globalVariables(c("category","pl","Description"))
