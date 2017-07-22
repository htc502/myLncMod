#' @input: 0
#' @output: exp.rda

PlotGeneExpDist <- function(mtr) {
    ##to find out what is the value that we assume the genes with this exp value is not expressed...
    med <- apply(mtr,1,median)
    pdf('PlotGeneExpDist.pdf')
    hist(med,breaks=100,xlab='Median of gene expression')
    dev.off()
    med
}

zero2na <- function(mtr,zero=log2(0.001),rmLowexp=F) {
    ##fill zero expression value with NA, those NA values will not be taken into consideration in the screening step
    mtr[ mtr <= zero ] <- NA
    if(rmLowexp) {
    apply(mtr, 1, function(e) sum(is.na(e))) -> nacount
    print('# of genes and samples:')
    print(dim(mtr))
    print(paste0('rm ', sum(nacount >= ncol(mtr)*.5), 'genes'))
    mtr1 <- mtr[ nacount < ncol(mtr)*.5,,drop=F ]
        } else {
        mtr1 = mtr
        }
    mtr1
}

rmOutlierSamples <- function(mtr,plot=T) {
    ##Here we define outlier as value > q3+1.5*IQR or < q1-1.5*iqr
    t( apply(mtr,1,function(e) {
        iqr <- IQR(e, na.rm=T)
        q13 <- quantile(e,na.rm=T,probs=c(.25,.75))
        idx1 <- e > q13[2] + 1.5*iqr
        idx2 <- e < q13[1] - 1.5*iqr
        sum(idx1|idx2)
    }) ) -> nbad
if(plot) {
    pdf('rmOutlierSamples.pdf')
    hist(nbad,breaks=100,xlab='Number of samples')
    dev.off()
    }

    t( apply(mtr,1,function(e) {
        iqr <- IQR(e, na.rm=T)
        q13 <- quantile(e,na.rm=T,probs=c(.25,.75))
        idx1 <- e > q13[2] + 1.5*iqr
        idx2 <- e < q13[1] - 1.5*iqr
        e[idx1 | idx2] <- NA
        e}) ) -> res
    res
}

PlotGeneCVDist <- function(mtr) {
    cv <- apply(mtr, 1, function(e) {
        cve <- IQR(e,na.rm=T)/median(e,na.rm=T)
    })
    pdfname <- 'GeneCVDist.pdf'
    pdf(pdfname)
    hist(cv, breaks=100,
         xlab='Coeffecient of variance')
    dev.off()
    return(cv)
}

rmlowCVGenes <- function(mtr,cv_cutoff,cv=NULL) {

    if(is.null(cv)) cv <- apply(mtr, 1, function(e) {
        cve <- IQR(e,na.rm=T)/median(e,na.rm=T)
    })
    mtr1 <- mtr[ cv > cv_cutoff, ]
    print(paste0('rm ',nrow(mtr)-nrow(mtr1),' genes with small cv'))
    mtr1
}

PCAplot <- function(mtr,label=NULL,col=NULL,fname=NULL) {
    ##to check there are substurcture or not...
    tmp <- t(mtr)
    tmp[ is.na(tmp) ] <- 0
    pca.res <- prcomp(tmp)
    if(is.null(col)) col='black'
    if(is.null(fname)) fname='PCAplot.pdf'
    pdf(fname)
    plot(pca.res)
    plot(pca.res$x[,1:2],xlab='PC1',
         ylab='PC2',
         pch=20,
         cex=.6,
         col=col)
    if(!is.null(label)) {
        text(pca.res$x[,1],
             pca.res$x[,2],
             label)
    }

    dev.off()
}


prep <- function(mtr, zero) {
    r1 <- zero2na(mtr,zero=zero)
    r2 <- rmOutlierSamples(r1)
    PlotGeneCVDist(r2) -> r2_cv
    readline(prompt="Input the CV cutoff:") -> cv_cutoff
    r3 <- rmlowCVGenes(r2,cv_cutoff,cv=r2_cv)
    PCAplot(r3)
    r3
}





