#'@title  火山图
#'
#'@description  可视化展示差异分析结果
#'
#'@details  火山图
#'
#'@param DEGALL 差异分析的所有结果，其中symbol为其中一列否则直接使用行名，pvaltype需要指定为相应的列名
#'@param fc 差异倍数
#'@param pval 显著性p值大小
#'@param lab 需要在图上标出的前多少个基因（主要按差异倍数排序），或者使用向量直接传入基因名字（需与symbol或者行名为同种类型）
#'@param pvaltype p值类型
#'@return 火山图
#'@importFrom ggplot2 ggplot
#'@export
#'@examples
#'VolcanoPlot(DEGALL,lab = 10)
#'VolcanoPlot(DEGALL,lab = c("CLEC4M", "CLEC1B", "STAB2", "CLEC4G", "BMPER"))


VolcanoPlot <- function(DEGALL, fc=2, pval=0.05,lab=0,pvaltype = 'PValue') {
  #基因上下调分组
  geneList <- DEGALL[which((!is.na(DEGALL$logFC))&(!(is.na(DEGALL[,match(pvaltype,colnames(DEGALL))])))),]
  geneList$threshold <- 'nosig'
  geneList$threshold[geneList$logFC > log(fc,2) & geneList[,match(pvaltype,colnames(DEGALL))]<pval] <- "up"
  geneList$threshold[geneList$logFC < -log(fc,2) & geneList[,match(pvaltype,colnames(DEGALL))]<pval] <- "down"

  #判断symbol
  if(!'symbol' %in% colnames(DEGALL)) {
    geneList$symbol <- rownames(geneList)
  }

  if(pvaltype %in% colnames(DEGALL)) {
    geneList$logPvalue <- -log10(geneList[,match(pvaltype,colnames(DEGALL))])
  }else {
    stop("pvaltype必须在DEGALL的列名中")
  }
  if(is.numeric(lab)){
    texts<-geneList[which(abs(geneList$logFC)>=1),]
    if(sum(is.na(texts[order(texts[,match(pvaltype,colnames(texts))])[0:lab],'symbol']))>0){
      message('lab: symbol中包含NA或全是NA')
    }
    texts <- texts[order(texts[,match(pvaltype,colnames(texts))])[0:lab],]
  }else if(is.vector(lab) & length(intersect(texts$symbol,lab))>0){
    texts <- texts[intersect(texts$symbol,lab),]
  }else{
    message('lab: 所需标明的symbol不包含在DEGALL中')
  }

  lim <- max(max(geneList$logFC), abs(min(geneList$logFC)))+0.5
  volcano <- ggplot(data=geneList, aes(x=logFC,
                                       y = logPvalue,color=threshold))
  p<-volcano+geom_point(alpha=1,size=2) + #alpha 透明度
    scale_colour_manual(values  = c('#BC3C28','#0072B5','gray'), limits = c("up","down","nosig")) +

    geom_vline(xintercept = c(-log(fc,2), log(fc,2)), color = 'gray', size = 0.3,linetype=4) +
    geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.3,linetype=4) +

    theme(axis.text=element_text(size=30), # raw 13
          axis.title=element_text(size=30))+  # raw 16+
    labs(x = 'log2(FoldChange)', y = '-log10(PValue)')+
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
    theme(legend.key = element_rect(fill = 'transparent'),
          legend.background = element_rect(fill = 'transparent'),
          legend.position = c(0.85, 0.93),legend.title=element_blank(),legend.text = element_text(size = 30))+
    geom_text(data=texts,  aes(x=logFC,y=logPvalue,label=symbol),hjust=0,vjust=0,size=5,show.legend=F)

  p
}


#'@title  热图
#'
#'@description  可视化展示差异分析结果，画的所有差异基因的热图，会花费一些时间
#'
#'@details  热图
#'
#'@param DEGsig 差异分析的显著结果，行名为基因
#'@param exprdata 表达数据，一般为标准化后的数据（log2(tpm+1)）
#'@param group 分组向量
#'@return 热图
#'@importFrom pheatmap pheatmap
#'@importFrom gplots greenred
#'@export
#'@examples
#'exprdata <- log2(LIHCexp.tpm+1)
#'HeatmapPlot(DEGsig,exprdata,group = group)


HeatmapPlot <- function(DEGsig,exprdata,group = group){
  loc <- order(group,colSums(exprdata),decreasing = T)
  exprdata <- exprdata[,loc]
  Group = group[loc]

  degName<-rownames(DEGsig)
  degDa <-as.matrix(exprdata[degName,]) #标准化的数据

  annotation_col <- data.frame('Sample' = Group)
  rownames(annotation_col) <- colnames(degDa)

  degDa[which(degDa > quantile(degDa,0.95,na.rm=T))] <- quantile(degDa,0.95,na.rm=T)
  degDa[which(degDa < quantile(degDa,0.05,na.rm=T))] <- quantile(degDa,0.05,na.rm=T)
  degDa <- t(scale(t(degDa),scale = F))
  p<-pheatmap(degDa,
              color = greenred(75),
              #main = 'heatmap', # 图标题
              scale = 'none', #值集中的方向，“column”，“row” “none”
              annotation_col = annotation_col, #列注释
              #annotation_row = annotation_row, #行注释
              #legend_labels = NA,
              cluster_cols = F,          # 以列聚类
              border_color = NA,
              #cluster_rows = FALSE,         # 以行聚类
              clustering_method = "complete", # 聚类方法 “complete” “average” “median”
              show_rownames = F, #不显示行名
              show_colnames = F, #不显示列名
              #gaps_row = 1169, # 分片
              fontsize = 10,
              angle_col=45) #列名的显示角度
  print(p)
}


#'@title  箱线图
#'
#'@description  箱线图展示部分基因的差异表达
#'
#'@details  箱线图
#'
#'@param data 表达数据，行名为基因，未进行log转化
#'@param genes 要展示的基因
#'@param group 分组向量
#'@param ref 对照组
#'@return 箱线图
#'@importFrom ggpubr ggboxplot
#'@export
#'@examples
#'genes=rownames(LIHCexp.tpm)[1:20]
#'ggboxplotrev(LIHCexp.tpm,genes,group,"Tumor")



ggboxplotrev <- function(data,genes,group,ref=NA){
  if(!is.na(ref)){
    if(!ref %in% group){stop('lab: ref不在group中')}
  }else{
    ref <- unique(group)[1]
  }
  if(length(intersect(genes,rownames(data))) == 0){
    message('全部基因均不在数据集中')
    next
  }else if(length(intersect(genes,rownames(data))) != length(genes)){
    message(paste0('部分基因不在数据集中：',setdiff(genes,rownames(data)),sep='\n'))

  }
  genes <- intersect(rownames(data),genes)
  boxplotdata <- t(data[genes,])
  boxdata <- lapply(colnames(boxplotdata), function(x) data.frame(expression=log2(boxplotdata[,x]+1),gene=rep(x,nrow(boxplotdata)),group))
  boxdata <- do.call(rbind,boxdata)
  boxdata$group <- ifelse(group == ref,"#00ba38","#619cff")
  boxdata$group <- factor(boxdata$group,labels = c('Tumor','Normal'))
  p=ggboxplot(boxdata, x="gene", y="expression", color = "group",
              ylab="Gene expression",
              xlab="",
              add = "jitter",size = 1,axis.line =2)+
    stat_compare_means(label = "p.signif",label.y = max(boxdata[,1]),aes(group=group))+rotate_x_text(60)
  return(p)
}

