#---------------------------------------------------------------------------------------
#
# Miscellanous functions for quantitative analysis
#
#---------------------------------------------------------------------------------------


getAUC <- function(ordered.list, gene.list,plot=TRUE,title="ROC"){
  ordered.list <- tolower(unique(as.character(ordered.list)))
  gene.list <- tolower(unique(as.character(gene.list)))
  if(all(gene.list %in% ordered.list)){
    #print("ALL IN.\n")
    #full.gene.list <- gene.list
  }else{
    #full.gene.list <- gene.list
    #print(paste(length(which(gene.list %in% ordered.list)),"/",length(gene.list),"\n"))
    gene.list <- gene.list[which(gene.list %in% ordered.list)]
  }
  outcome <- rep(0,length(ordered.list))
  outcome[which(ordered.list %in% gene.list)] <- 1
  tpr <- cumsum(outcome) / length(gene.list)
  xoutcome <- (outcome*-1)+1
  fpr <- cumsum(xoutcome) / sum(xoutcome)
  #print(max(tpr))
  AUC <- sum(tpr*c(0,diff(fpr)))
  if(plot){
    plot(c(0,fpr),c(0,tpr), xlab="False Positive Rate",ylab="True Positive Rate",
         main=paste(title,"\nAUC = ",signif(AUC,3),"\n",
           #paste(length(which(full.gene.list %in% ordered.list)),"/",length(full.gene.list)),
           sep=""),type="l",lty=1,lwd=2)
    abline(a=0, b=1, lty=2, col=rgb(0,0,0,0.5))
  }
 
  return(AUC)
}

getPR <- function(ordered.list, gene.list, plot=TRUE){
  nomap <- which(ordered.list=="")
  if(length(nomap!=0)){ordered.list <- tolower(unique(as.character(ordered.list)[-nomap]))}
  ordered.list <- tolower(unique(as.character(ordered.list)))
  gene.list <- tolower(unique(as.character(gene.list)))


  if(all(gene.list %in% ordered.list)){
    #print("ALL IN.\n")
  }else{
    #print(paste(length(which(gene.list %in% ordered.list)),"/",length(gene.list),"\n"))
    gene.list <- gene.list[which(gene.list %in% ordered.list)]
  }

  outcome <- rep(0,length(ordered.list))
  outcome[which(ordered.list %in% gene.list)] <- 1
  precision <- cumsum(outcome) / seq(1,length(outcome))
  recall <- cumsum(outcome)/ length(gene.list)

  if(plot){
    plot(recall, precision, xlab="Recall",ylab="Precision",type='l')
  }
  AveP <- sum(precision*outcome)/length(gene.list)
  return(AveP)
}


geomean <- function(x){exp(mean(log(x),na.rm=TRUE))}

