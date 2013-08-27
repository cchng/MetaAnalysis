#----------------------------------------------------------------------
#
# Description: Functions for plotting graphs, heatmaps, etc.
# 
#----------------------------------------------------------------------


library(ggplot2)



genHeatMapDat <- function(genelist,gset,res,name,direction){
  rein <- res[which(res$NCBIids %in% as.character(genelist)),c("ID","NCBIids",toupper(direction))]
  
  print(all(as.character(genelist) %in% res$NCBIids))
  print(paste(length(genelist),nrow(rein)))
  if(nrow(rein)==0){
    print(name)
    return(NULL)
  }
  rein <- aggregate(as.formula(paste(toupper(direction),"~ NCBIids")),function(x){min(x,na.rm=TRUE)},data=rein)

  merged <- merge(res,rein,by=c("NCBIids",toupper(direction)))
  
  gset.ord <-gset[as.character(merged$ID),sampleNames(gset)[order(gset$Grouping)]]
  gdat <- data.frame(exprs(gset.ord))
  row.names(res) <- res$ID
  labels <- res[row.names(gdat),c("GeneSymbols","NCBIids")]
  all(row.names(gdat)==row.names(labels))
  row.names(gdat) <- labels$NCBIids

  # order by input
  gdat.ord <- gdat[as.character(genelist),]
  missing <- 0
  for(i in 1:nrow(gdat.ord)){
    if(row.names(gdat.ord)[i] %in% labels$NCBIids){
    row.names(gdat.ord)[i]<-labels[which(labels$NCBIids %in% row.names(gdat.ord)[i]),"GeneSymbols"]
  } else{
    missing <- missing +1
    row.names(gdat.ord)[i] <-paste("NA",missing,sep="")
  }

  }
  
  return(gdat.ord)
           
}






# Ref: http://stackoverflow.com/questions/6973394/functions-available-for-tufte-boxplots-in-r
panel.tuftebxp <- 
function (x, y, box.ratio = 1, box.width = box.ratio/(1 + box.ratio), horizontal=FALSE,
    pch = box.dot$pch, col = box.dot$col, 
    alpha = box.dot$alpha, cex = box.dot$cex, font = box.dot$font, 
    fontfamily = box.dot$fontfamily, fontface = box.dot$fontface, 
    fill = box.rectangle$fill, varwidth = FALSE, notch = FALSE, 
    notch.frac = 0.5, ..., levels.fos = if (horizontal) sort(unique(y)) else sort(unique(x)), 
    stats = boxplot.stats, coef = 1.5, do.out = TRUE, identifier = "bwplot") 
{
    if (all(is.na(x) | is.na(y))) 
        return()
    x <- as.numeric(x)
    y <- as.numeric(y)
    box.dot <- trellis.par.get("box.dot")
    box.rectangle <- trellis.par.get("box.rectangle")
    box.umbrella <- trellis.par.get("box.umbrella")
    plot.symbol <- trellis.par.get("plot.symbol")
    fontsize.points <- trellis.par.get("fontsize")$points
    cur.limits <- current.panel.limits()
    xscale <- cur.limits$xlim
    yscale <- cur.limits$ylim
    if (!notch) 
        notch.frac <- 0
    #removed horizontal code
     blist <- tapply(y, factor(x, levels = levels.fos), stats, 
            coef = coef, do.out = do.out)
        blist.stats <- t(sapply(blist, "[[", "stats"))
        blist.out <- lapply(blist, "[[", "out")
        blist.height <- box.width
        if (varwidth) {
            maxn <- max(table(x))
            blist.n <- sapply(blist, "[[", "n")
            blist.height <- sqrt(blist.n/maxn) * blist.height
        }
        blist.conf <- if (notch) 
            sapply(blist, "[[", "conf")
        else t(blist.stats[, c(2, 4), drop = FALSE])
        ybnd <- cbind(blist.stats[, 3], blist.conf[2, ], blist.stats[, 
            4], blist.stats[, 4], blist.conf[2, ], blist.stats[, 
            3], blist.conf[1, ], blist.stats[, 2], blist.stats[, 
            2], blist.conf[1, ], blist.stats[, 3])
        xleft <- levels.fos - blist.height/2
        xright <- levels.fos + blist.height/2
        xbnd <- cbind(xleft + notch.frac * blist.height/2, xleft, 
            xleft, xright, xright, xright - notch.frac * blist.height/2, 
            xright, xright, xleft, xleft, xleft + notch.frac * 
                blist.height/2)
        xs <- cbind(xbnd, NA_real_)
        ys <- cbind(ybnd, NA_real_)
        panel.segments(rep(levels.fos, 2), c(blist.stats[, 2], 
            blist.stats[, 4]), rep(levels.fos, 2), c(blist.stats[, 
            1], blist.stats[, 5]), col = box.umbrella$col, alpha = box.umbrella$alpha, 
            lwd = box.umbrella$lwd, lty = box.umbrella$lty, identifier = paste(identifier, 
                "whisker", sep = "."))

        if (all(pch == "|")) {
            mult <- if (notch) 
                1 - notch.frac
            else 1
            panel.segments(levels.fos - mult * blist.height/2, 
                blist.stats[, 3], levels.fos + mult * blist.height/2, 
                blist.stats[, 3], lwd = box.rectangle$lwd, lty = box.rectangle$lty, 
                col = box.rectangle$col, alpha = alpha, identifier = paste(identifier, 
                  "dot", sep = "."))
        }
        else {
            panel.points(x = levels.fos, y = blist.stats[, 3], 
                pch = pch, col = col, alpha = alpha, cex = cex, 
                 identifier = paste(identifier, 
                  "dot", sep = "."))
        }
        panel.points(x = rep(levels.fos, sapply(blist.out, length)), 
            y = unlist(blist.out), pch = plot.symbol$pch, col = plot.symbol$col, 
            alpha = plot.symbol$alpha, cex = plot.symbol$cex*0.5, 
            identifier = paste(identifier, "outlier", sep = "."))

}



combinedPPlot <- function(data, meta.p=NULL,genes=NULL,line,quantile.plot=TRUE,plot.legend=FALSE){
#    if(genes==NULL) genes <- row.names(data)
    set.seed(1)
    dnum <- ncol(data)
    colors <- rainbow(dnum)
    pvalues <- as.numeric(data[line,1:dnum])
    if(quantile.plot){
      layout(matrix(c(1,2), 1, 2, byrow = TRUE))
      x <- -log10(sort(pvalues))
      qqplot(-log10(runif(1000)),x,main="Q-Q plot", pch=19, #bg=colors[order(pvalues,na.last=NA)],
             ylab="Data p-value quantiles",xlab="Uniform distribution quantiles")
      abline(0,1,col="grey50",lty=2)
      if(plot.legend)
      legend("bottomright",colnames(data),pch=19,col=colors,bty='n')
    }
    plot(-log10(pvalues),seq(1,dnum),main=paste(genes[line],"\n",row.names(data)[line]),axes=FALSE,xlab="-log10(pval)",ylab="",xlim=c(0,15),ylim=c(0,dnum+1),pch=19)
    box()
    axis(1)
    axis(2,at=seq(1,length(pvalues)),labels=names(data)[1:dnum],las=2)
    abline(v=-log10(0.05),col="light blue",lwd=2)
    for(k in 1:length(pvalues)){
      abline(h=k,lty=2,col="grey50")
    }
    if(!is.null(meta.p)) abline(v=-log10(meta.p),col="dark blue",lwd=2.7)
}




# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)

      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)

      numPlots = length(plots)

      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
            # Make the panel
            # ncol: Number of columns of plots
            # nrow: Number of rows needed, calculated from # of cols
            layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                                 ncol = cols, nrow = ceiling(numPlots/cols))
          }

     if (numPlots==1) {
           print(plots[[1]])

         } else {
               # Set up the page
               grid.newpage()
                   pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

                   # Make each plot, in the correct location
                   for (i in 1:numPlots) {
                           # Get the i,j matrix positions of the regions that contain this subplot
                           matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

                                 print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                                                           layout.pos.col = matchidx$col))
                         }
             }
  }

# Quantile-quantile plot for ggplot
myqq <- function(sampleq,gene){
    set.seed(1)
      x <- -log10(runif(1000));y <- -log10(as.numeric(sampleq))
      sx <- sort(x); sy <- sort(y)
      lenx <- length(sx)
      leny <- length(sy)
      if (leny < lenx)sx <- approx(1L:lenx, sx, n = leny)$y
      if (leny > lenx)sy <- approx(1L:leny, sy, n = lenx)$y
      return(data.frame(X=as.numeric(sx),Y=as.numeric(sy),Gene=rep(gene,length(sx))))

  }
