
#'@title  GEO数据下载
#'
#'@description  提供GSE和GPL编号下载数据
#'
#'@details  提供GSE和GPL编号下载数据
#'
#'@param GEOnumber 数据集GES编号
#'@param gplnumber 数据集GPL编号
#'@return GEO数据集和样本信息
#'@import GEOquery
#'@import idmap2
#'@import idmap1
#'@import idmap3
#'@export
#'@examples
#'GEOdownloaddata('GSE100927','GPL17077')


GEOdownloaddata <- function(GEOnumber,gplnumber){
  if(!grepl('^GSE',GEOnumber) | !grepl('^GPL',gplnumber)){    stop('输入的GSE编号或者GPL编号有误')  }
  gset <-  getGEO(GEOnumber, GSEMatrix =TRUE, AnnotGPL=F)
  gplloc <- ifelse(any(grepl(gplnumber,names(gset))),grep(gplnumber,names(gset)),1)
  gset <- gset[[gplloc]]

  if(!is.na(match(gplnumber,unique(idmap1::p2s_df$gpl)))){
    ids <- idmap1::getIDs(gplnumber)
  }else if(!is.na(match(gplnumber,idmap2::gpl_list$gpl))){
    ids <- idmap2::get_soft_IDs(gplnumber)
  }else if(!is.na(match(gplnumber,idmap3::gpl_list$gpl))){
    ids <- idmap3::get_pipe_IDs(gplnumber)
  }else if(!is.null(fData(gset))){
    ids <- fData(gset)
    if(any(grepl('symbol',tolower(colnames(ids))))){
      ids <- ids[,c(1,grep('symbol',tolower(colnames(ids))))]
    }else{
      print(colnames(ids))
      stop('在列名中未找到symbol')
    }
  }else{
    stop('请使用另外的方法下载数据')
  }

  expr <- exprs(gset)
  colnames(ids)[1:2] <- c('ID','symbol')
  ids <- ids[match(rownames(expr),ids$ID),]

  #data<- expr[!apply(expr,1,sum)==0,]
  data <- data.frame(genename=ids[,2],expr)
  data$median <- apply(data[,-1],1,median)#计算每行的中位数，添加到 data数据中
  data=data[order(data$genename,data$median,decreasing = T),]#排序
  data=data[!duplicated(data$genename),]#去除重复的基因名
  if(length(which(data$genename=='---'| is.na(data$genename)))>0){
    data <- data[-which(data$genename=='---'| is.na(data$genename)),]
  }
  rownames(data)=data$genename#把基因名变成行名

  rnaexpr <- data[,-c(1,ncol(data))]
  pdata <- pData(gset)
  return(list(rnaexpr,pdata))
}

