#----------------------------------------------------------------------
#
# Functions
#
#
#----------------------------------------------------------------------

# load libraries
suppressMessages(library(Biobase))
suppressMessages(library(GEOquery))
suppressMessages(library(limma))
suppressMessages(library(lumi))
suppressMessages(library(gplots))
#suppressMessages(library(qvalue))
suppressMessages(library(RColorBrewer))
suppressMessages(library(stringr))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library(bootstrap))
suppressMessages(library(lattice))

#----------------------------------------------------------------------
# Differential Expression Analysis
#----------------------------------------------------------------------



# compute one sided p values for each gene for both up and down regulated 
one.sided <- function (tT){
  tT[which(tT$t>0),"UP"] <- tT[which(tT$t>0),"P.Value"]/2
  tT[which(tT$t<0),"UP"] <- 1-tT[which(tT$t<0),"P.Value"]/2
  tT$DOWN <- 1-tT$UP
  return(tT)

}

# bonferroni correction for FWER
bonferroni <- function(x){
  corrected.p.val <- length(x)*min(x,na.rm=TRUE)
  if(corrected.p.val==Inf)return(NA)
  else if(corrected.p.val > 1) return(1)
  return(corrected.p.val)
}

# collapse probes to genes by taking bonferroni corrected min one tailed p-value
collapseUp <- function(res){
    aggregate(res$UP,list(as.character(res$NCBIids)),bonferroni)
}

# analogous to collapseUp
collapseDown <- function(res){
   aggregate(res$DOWN,list(as.character(res$NCBIids)),bonferroni)
}

# no collapse
collapseNone <- function(res,up=TRUE){
  if(up){
    aggregate(res$UP,list(row.names(res)),mean)
  }else{
    aggregate(res$DOWN,list(row.names(res)),mean)
  }
    
}


# Eliminate those without unique mappings and without p-values
clean <- function(results, id){
  hasUniqueMapping <- results[which(!(id=="" | grepl("|",id,fixed=TRUE))),]
  nopval <-hasUniqueMapping[is.na(hasUniqueMapping$P.Value),]
  haspval <- hasUniqueMapping[!is.na(hasUniqueMapping$P.Value),]
  return(haspval)
}

# Differential expression results
getDE <- function(res,threshold=0.05,pcol="P.Value",adjpcol="adj.P.Val"){
  res$adj.P.val <- p.adjust(res[,pcol],method='fdr')
  diffex <- res[which(res[,adjpcol]<threshold),]

  if(nrow(diffex)==0){
    upregulated <- 0
    downregulated <- 0
  }
  upTable <- diffex[which(diffex$t>0),]
  downTable <- diffex[which(diffex$t<0),]

  if(nrow(upTable)==0){
    upregulated <- 0
   } else{
        upregulated <- length(unique(upTable$NCBIid))
     }

  if(nrow(downTable)==0){
    downregulated <- 0
  } else {
    downregulated <- length(unique(downTable$NCBIid))
  }
  
  return(c(upregulated+downregulated,upregulated,downregulated))
}

# Upregulated results
getUP <- function(res){
  
#  diffex <- res[which(res$adj.P.Val<0.05),]

 # if(nrow(diffex)==0)return(NULL)

  upTable <- res[which(res$t>0),]
  

  if(nrow(upTable)==0){
     writeLines(paste("Up-regulated: ","0",sep=" "))
     return()
   }

  upGenes <- aggregate(upTable$adj.P.Val,list(upTable$Gene.ID),min)

  writeLines(paste("Up-regulated: ",nrow(upGenes),sep=" "))
row.names(upGenes)<- as.character(upGenes$Group.1)
  return(upGenes)
}

# Downregulated expression results
getDOWN <- function(res){
  
  #diffex <- res[which(res$adj.P.Val<0.05),]

  #if(nrow(diffex)==0)return(NULL)


  downTable <- res[which(res$t<0),]

  if(nrow(downTable)==0){
     writeLines(paste("Down-regulated: ","0",sep=" "))
     return()
  }
  
  
  downGenes <- aggregate(downTable$adj.P.Val,list(downTable$Gene.ID),min)
  
  writeLines(paste("Down-regulated: ",nrow(downGenes),sep=" "))
  return(downGenes)
}



#----------------------------------------------------------------------
# Meta Analysis
#----------------------------------------------------------------------

# combines result sets to a list
rset2metagene <- function(res.dat){
   results <- read.table(res.dat, header=TRUE, sep="\t",row.names="ID",fill=TRUE,quote='"')
   return(results)
}

# merge several datasets change colnames to bypass merge bug
mymerge <- function(x,y){
  #print(colnames(x))
  #print(colnames(y))
  m<-merge(x,y,by="Group.1",all=TRUE)
  row.names(m) <- m$Group.1
  colnames(m)[2:ncol(m)] <- make.names(c(2:ncol(m)))
  return(m)
}


#----------------------------------------------------------------------
# Fisher
#----------------------------------------------------------------------

# arg1: merged: A data frame with first column (Group.1) as gene IDs and the following columns as gene p values for each study.
Fisher <- function(merged,k.threshold=3){

  F.stat <- function(dat){
    pvals <- as.numeric(dat[2:length(dat)])
    k <- length(which(!is.na(pvals)))
    F <- -2 * sum(log(pvals),na.rm=TRUE)
    return(c(dat["Group.1"],"F"=F, "k"=k))
  }

  meta <- data.frame(t(apply(merged,1,F.stat)))
  meta$F <- as.numeric(levels(meta$F))[meta$F]
  meta$k <- as.numeric(levels(meta$k))[meta$k]
  meta$df <- 2 * meta$k
  
  # apply threshold and clean
  meta <- meta[which(meta$k>=k.threshold & meta$Group.1!=""),]

  # generate p values and apply multiple correction
  meta$P <- pchisq(meta$F, df=meta$df,lower.tail=FALSE) 
  meta$Q <- qvalue(meta$P)$qvalues
  meta$BH <- p.adjust(meta$P, method="fdr")

  # order by P values
  meta.ord <- meta[order(meta$BH),]
  return(meta.ord)
}



JK.Fisher <- function(merged, fdr.threshold=0.01, k.threshold=3,correct=TRUE){

   F.stat <- function(pvals){
    k <- length(which(!is.na(pvals)))
    F <- -2 * sum(log(pvals),na.rm=TRUE)
    P <- pchisq(F, df=2*k, lower.tail=FALSE)
    return(P)
  }

    meta <- data.frame(t(apply(merged,1,function(x){return(jackknife(as.numeric(x[2:ncol(merged)]),F.stat)$jack.values)})))
 
   row.names(meta) <- merged$Group.1
   print(dim(meta))
   
   # apply threshold and clean
   k <- apply(merged,1,function(x){length(which(!is.na(as.numeric(x[2:ncol(merged)]))))})
   meta <- meta[which(k>=k.threshold & merged$Group.1!=""),]
   print(dim(meta))
   
   # convert to q values
   #meta.q <- data.frame(apply(meta,2,function(x){qvalue(x)$qvalues}))

   # convert to BH values  
   meta.q <- data.frame(apply(meta,2,function(x){p.adjust(x,method='fdr')}))

   if(!correct) meta.q <- meta

   
   meta.q$LOO <-  apply(meta.q,1,function(x){all(as.numeric(x)<fdr.threshold)})

   colnames(meta.q)[1:ncol(merged)-1] <- names(merged)[2:ncol(merged)]
   return(meta.q)

}


## BS.Fisher <- function(merged, fdr.threshold=0.01, k.threshold=3){

##    F <- function(pvals){
##     k <- length(which(!is.na(pvals)))
##     F <- -2 * sum(log(pvals),na.rm=TRUE)
##     P <- pchisq(F, df=2*k, lower.tail=FALSE)
##     return(P)
##   }


##    pass <- function(x){
##      return(all(x<fdr.threshold))
##    }

   
     

##    meta <- data.frame(apply(merged,1,function(x){return(bootstrap(as.numeric(x[2:ncol(merged)]),1000,F,func=pass)$func.thetastar)}))
## #   row.names(meta) <- merged$Group.1
##    print(dim(meta))
   
##    # apply threshold and clean
##   # k <- apply(merged,1,function(x){length(which(!is.na(as.numeric(x[2:ncol(merged)]))))})
##    #meta <- meta[which(k>=k.threshold & merged$Group.1!=""),]
##    #print(dim(meta))
  

   
##    # convert to q values
##   # meta.q <- data.frame(apply(meta,2,function(x){qvalue(x)$qvalues}))
## #   meta.q$LOO <-  apply(meta.q,1,function(x){all(as.numeric(x)<fdr.threshold)})

## #   colnames(meta.q)[1:ncol(merged)-1] <- names(merged)[2:ncol(merged)]
##    #meta$Gene.ID <- merged[,1]
##    return(meta)

## }




#----------------------------------------------------------------------
# Meta-rank Analysis
#----------------------------------------------------------------------

metarank <- function(merged){

  k <- apply(merged[,c(2:ncol(merged))],1,function(x){length(which(!is.na(x)))})

  # take out those that have k >= 3
  merged <- merged[which(k>=3),]
  
  within.rank <- apply(merged[,c(2:ncol(merged))],2,
                       function(x){rank(x,na.last="keep")}) # rank within study

  merged$Ave.Rank <- rowMeans(within.rank,na.rm=TRUE)

  sqindex <- (within.rank - merged$Ave.Rank)^2
  merged$Heterogeneity <- rowSums(sqindex, na.rm=TRUE)
  merged$Meta.Rank <- rank(merged$Ave.Rank) # meta rank
  
  # order by meta rank
  merged.ord <- merged[order(merged$Meta.Rank),]
  return(merged.ord)
}


aggregrank <- function(merged){
    require(RankAggreg)
    within.rank <- apply(merged[,c(2:ncol(merged))],2,rank) # rank within study
    genelist <- apply(within.rank,2,function(x){
      ord.genes <- row.names(merged)[order(x)]
    })

   agg <- RankAggreg(t(genelist),10,method='CE',distance='Spearman',rho=0.1,N=100,ConvIn=5)
   return(genelist)
}


#-------------------------------------------------
#
# Raw expression values
#
#--------------------------------------------------

# genelist: list of gene ids
# gset: data gset (expression set)
# res: differential expression analysis results data frame
# direction: column name for scores
# name: dataname
genDF <- function(genelist,gset,res,name,direction){

  # get all probes which are mapped to genes in the gene list.
  rein <- res[which(res$NCBIids %in% as.character(genelist)),c("ID","NCBIids",toupper(direction))]
  
  if(nrow(rein)==0){ # if dataset does not have all genes
    print(name)
    return(NULL)
  }

  # collapse probes to genes and add results
  rein <- aggregate(as.formula(paste(toupper(direction),"~ NCBIids")),function(x){min(x,na.rm=TRUE)},data=rein)
  merged <- merge(res,rein,by=c("NCBIids",toupper(direction)))

  g.dat <- as(gset,"data.frame")
  g.ctrl <- g.dat[which(g.dat$Grouping=="Control"),] 
  g.asd <- g.dat[which(g.dat$Grouping=="ASD"),]

  # get expression values for genes
  gc <- data.frame(g.ctrl[,make.names(merged$ID)])
  ga <- data.frame(g.asd[,make.names(merged$ID)])

  if(length(merged$NCBIids)==1)
    merged$NCBIids <-paste("Gene",merged$NCBIids,sep=".")

  all(colnames(gc)==merged$ID)
  all(colnames(ga)==merged$ID)
  colnames(gc)<- merged$NCBIids
  colnames(ga) <- merged$NCBIids 
  
  dat.df <- rbind(cbind(Gene=gc,Disease="Control",Study=name),cbind(Gene=ga,Disease="ASD",Study=name))
  return(dat.df)
}




## sexLinkedReanalysis <-function(genelist,gset,res,name){
##   rein <- res[which(res$NCBIids %in% as.character(genelist)),c("ID","NCBIids")]
##   if(nrow(rein)==0){
##     print(name)
##     return(NULL)
##   }
##   if(name %in% c("GSE18123.2","GSE28521","GSE37772","GSE6575","GSE28475")){
##     gset.ord <-gset[as.character(rein$ID),]
##     print(dim(exprs(gset.ord)))
##     names(pData(gset.ord))[which(tolower(names(pData(gset.ord))) %in% c("sex","gender"))[1]] <- "sex"
##     gset.ord$description <- factor(as.character(gset.ord$Grouping))
##     design <- model.matrix(~description + sex + 0,gset.ord)
##     colnames(design) <- make.names(colnames(design))
##     fit <- lmFit(gset.ord, design)
##     cont.matrix <- makeContrasts(descriptionASD-descriptionControl, levels=design)
##     fit2 <- contrasts.fit(fit, cont.matrix)
##     fit2 <- eBayes(fit2, 0.01)
##     tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
##     tT <- one.sided(tT)
##     tT$UP.adj.P.Val <- p.adjust(tT$UP,method="fdr")
##     tT$DOWN.adj.P.Val <- p.adjust(tT$DOWN,method="fdr")
##     return(tT)
##   }
##   return(res[,c("ID","UP","DOWN")])
 
## }






#########################################################################################
# OLD CODE

#  aggregate(res$UP,list(as.character(res$NCBIids)),function(x){length(x)*min(x,na.rm=TRUE)})
 # aggregate(res$DOWN,list(as.character(res$NCBIids)),function(x){length(x)*min(x,na.rm=TRUE)})

 # merged$Ave.Rank <- apply(within.rank,1,function(x){mean(x,na.rm=TRUE)}) # average rank across studies
 # merged$Heterogeneity1 <- apply(sqindex, 1, function(x){sum(x,na.rm=TRUE)})


##  perm.test <- function(merged,n){
##    set.seed(1)
## #   require(combinat)
##    p.averank <- c()
##    for(i in 1:n){
##      p.merged <- apply(merged[,c(2:ncol(merged))],2,sample)
##      within.rank <- apply(p.merged,2,rank) # rank within study
##      p.averank <-rbind(p.averank, rowMeans(within.rank,na.rm=TRUE))
##    }
##    return(p.averank)
##  }

## ## }

## #Fn <- ecdf(p)
## #pend <- 1-Fn(x)


## http://www.inside-r.org/packages/cran/miRtest/docs/limma.one.sided
## limma.one.sided <- function (fit, lower = FALSE) 
## {
##     se.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled
##     df.total <- fit$df.prior + fit$df.residual
##     pt(fit$t, df = df.total, lower.tail = lower)
## }


## isec <- read.table("/home/cchng/Documents/data/platform/isec_genes.txt")

## # Truncate gset to genes in intersected list
## gset2metagene <- function(gset.dat){
##   load(gset.dat)
##   print(gset.dat)

##   #load platform
##   #gpl <- annotation(gset)
##   #platf <- system(paste("find /home/cchng/Documents/data/platform/",gpl,"*.txt",sep=""),intern=TRUE)
##   #writeLines(paste("Platform found:",platf))
##   #annot <- read.table(platf,header=T,fill=TRUE,sep="\t")

##   geo <- gsub("-gset.RData","",gsub("/home/cchng/Documents/data/expression/ASD-QC-Processed/GSE","",gset.dat))
##   annot <- read.table(paste("/home/cchng/Documents/results/",geo,".txt",sep=""),sep="\t",header=T)
##   annot.in <- annot[which(as.character(annot$Gene.ID) %in% as.character(isec$V1)),]
##   print(nrow(annot.in))
  
##   ex <- exprs(gset)
##   ex <- ex[which(row.names(ex) %in% as.character(annot.in$ID)),]
##   print(dim(ex))


  
##   des <- pData(phenoData(gset))
##   print (all(colnames(ex) == row.names(des)))
##   gset <- new("ExpressionSet", phenoData = as(des, "AnnotatedDataFrame"), exprs = as.matrix(ex))
##   return(gset)

## }

