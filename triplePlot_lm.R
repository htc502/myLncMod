

triPlot <- function(m,e,tg,fname='tmp.pdf',xlab='Sample index',ylab='Log2(RPKM)',points.cex=.8,Ngrp=.25,groupMin=20,
                   labels=NULL) {
  
  mRNAexp0 <- tg
  TFexp0 <- e
  lncRNAexp0 <- m
  
  mRNAexp=mRNAexp0
  TFexp=TFexp0
  lncRNAexp=lncRNAexp0
  
  if(length(mRNAexp) != length(TFexp) | length(TFexp) != length(lncRNAexp)) stop('input of diff length')
  cutoffs <- quantile(lncRNAexp,na.rm=T,probs=c(Ngrp,1-Ngrp))
  low_grp <- lncRNAexp < cutoffs[1]
  mid_grp <- lncRNAexp >= cutoffs[1] & lncRNAexp <= cutoffs[2]
  high_grp <- lncRNAexp > cutoffs[2]
  grp <- rep('mid',length=length(mRNAexp))
  grp[low_grp] <- 'low'
  grp[high_grp] <- 'high'
  grp <- as.factor(grp) #used for lm
  
  ##remove outliers within each group
  TF_extremValue <- tapply(TFexp, grp, function(e) {
    IQR(e,na.rm=T)->iqr
    q13 <- quantile(e,na.rm=T,probs=c(.25,.75))
    list(lower=q13[1]-1.5*iqr,
         upper=q13[2]+1.5*iqr)
  })
  which((TFexp > TF_extremValue$low$upper) & low_grp) -> tf_lowgrp_outlier1
  which((TFexp < TF_extremValue$low$lower) & low_grp) -> tf_lowgrp_outlier2
  tf_lowgrp_outlierIdx <- c(tf_lowgrp_outlier1,tf_lowgrp_outlier2)
  which((TFexp > TF_extremValue$mid$upper) & mid_grp) -> tf_midgrp_outlier1
  which((TFexp < TF_extremValue$mid$lower) & mid_grp) -> tf_midgrp_outlier2
  tf_midgrp_outlierIdx <- c(tf_midgrp_outlier1,tf_midgrp_outlier2)
  which((TFexp > TF_extremValue$high$upper) & high_grp) -> tf_highgrp_outlier1
  which((TFexp < TF_extremValue$high$lower) & high_grp) -> tf_highgrp_outlier2
  tf_highgrp_outlierIdx <- c(tf_highgrp_outlier1,tf_highgrp_outlier2)
  
  mRNA_extremValue <- tapply(mRNAexp, grp, function(e) {
    IQR(e,na.rm=T)->iqr
    q13 <- quantile(e,na.rm=T,probs=c(.25,.75))
    list(lower=q13[1]-1.5*iqr,
         upper=q13[2]+1.5*iqr)
  })
  which((mRNAexp > mRNA_extremValue$low$upper) & low_grp) -> mrna_lowgrp_outlier1
  which((mRNAexp < mRNA_extremValue$low$lower) & low_grp) -> mrna_lowgrp_outlier2
  mRNA_lowgrp_outlierIdx <- c(mrna_lowgrp_outlier1,mrna_lowgrp_outlier2)
  which((mRNAexp > mRNA_extremValue$mid$upper) & mid_grp) -> mrna_midgrp_outlier1
  which((mRNAexp < mRNA_extremValue$mid$lower) & mid_grp) -> mrna_midgrp_outlier2
  mRNA_midgrp_outlierIdx <- c(mrna_midgrp_outlier1,mrna_midgrp_outlier2)
  which((mRNAexp > mRNA_extremValue$high$upper) & high_grp) -> mrna_highgrp_outlier1
  which((mRNAexp < mRNA_extremValue$high$lower) & high_grp) -> mrna_highgrp_outlier2
  mRNA_highgrp_outlierIdx <- c(mrna_highgrp_outlier1,mrna_highgrp_outlier2)
  
  ##set outliers'expression value to NA
      outlier <- list(TF=NULL,mRNA=NULL,common=NULL)
  if(!length(c(tf_lowgrp_outlierIdx,tf_highgrp_outlierIdx,tf_midgrp_outlierIdx))==0) {
    TFexp[ c(tf_lowgrp_outlierIdx,tf_highgrp_outlierIdx,tf_midgrp_outlierIdx) ] <- NA
    names(TFexp)[ c(tf_lowgrp_outlierIdx,tf_highgrp_outlierIdx,tf_midgrp_outlierIdx)] -> outlierSampleTF
    outlier$TF <- outlierSampleTF
  }
  
  if(!length(c(mRNA_lowgrp_outlierIdx,mRNA_highgrp_outlierIdx,mRNA_midgrp_outlierIdx))==0) {
    mRNAexp[ c(mRNA_lowgrp_outlierIdx,mRNA_highgrp_outlierIdx,mRNA_midgrp_outlierIdx) ] <- NA
    names(mRNAexp)[ c(mRNA_lowgrp_outlierIdx,mRNA_highgrp_outlierIdx,mRNA_midgrp_outlierIdx)] -> outlierSamplemRNA
    outlier$mRNA <- outlierSamplemRNA
  }
      
      if(length(outlier$TF) != 0 & length(outlier$mRNA) != 0) {
        outlier$common <- intersect(outlier$TF,outlier$mRNA)
        outlier$TF <- setdiff(outlier$TF,outlier$mRNA)
        outlier$mRNA <- setdiff(outlier$mRNA,outlier$TF)
      }
      
  ## number of available values in cleaned data within each group:
  nlowgrp <- sum(!is.na(TFexp[low_grp] + mRNAexp[low_grp]))
  nhighgrp <- sum(!is.na(TFexp[high_grp] + mRNAexp[high_grp]))
  nmidgrp <- sum(!is.na(TFexp[mid_grp] + mRNAexp[mid_grp]))
  
  if(nlowgrp < groupMin |
     nhighgrp < groupMin |
     nmidgrp < groupMin) {
    stop(paste0(M,' ',E,' ',Tg,': high/mid/low group size < predefined groupMin(default:20) after remove outlier'))
  }
  
    
    pdf(fname)
 dat0 <- data.frame(lncRNA=lncRNAexp0,TF=TFexp0,gene=mRNAexp0)
 rownames(dat0) <- names(lncRNAexp0)
    dat <- na.omit(dat0)
 idx <- match(rownames(dat),rownames(dat0)) 
    cols <- rep('grey',length(mRNAexp))
    cols[low_grp] <- 'green'
    cols[high_grp] <- 'red'
   cols <- cols[idx];rm(idx) 
    plot(dat$TF,dat$gene,col=cols,pch=20,cex=points.cex,
         xlab='TF expression',
         ylab='Gene expression')
    if(length(labels) != 0) {
      text(dat$TF,dat$gene,labels=labels,cex=points.cex)
      }
   
    if(length(outlier$TF) != 0) {
      outlierdat <- subset(dat,rownames(dat) %in% outlier$TF)
      if(nrow(outlierdat) != 0) { ## nrow == 0 means except TF, gene or lncRNA was set to NA in the early step and removed in line 91
    text(outlierdat$TF,outlierdat$gene,labels = rownames(outlierdat),cex = 0.6,col='brown')
      }
}
if(length(outlier$mRNA) != 0) {     
  outlierdat <- subset(dat,rownames(dat) %in% outlier$mRNA)
  if(nrow(outlierdat) != 0) { ##see line 104
    text(outlierdat$TF,outlierdat$gene,labels = rownames(outlierdat),cex = 0.6,col='blue')
  }
}
 if(length(outlier$common) != 0) {     
  outlierdat <- subset(dat,rownames(dat) %in% outlier$common)
    text(outlierdat$TF,outlierdat$gene,labels = rownames(outlierdat),cex = 0.6,col='skyblue')
  
}   
dat0 <- data.frame(lncRNA=lncRNAexp,TF=TFexp,gene=mRNAexp)
rownames(dat0) <- names(lncRNAexp)
    dat <- na.omit(dat0)
    idx <- match(rownames(dat),rownames(dat0))
   low_grp <- low_grp[ idx ] 
   high_grp <- high_grp[ idx ]
   mid_grp <- mid_grp[idx]
    lm_lq25 <- lm( dat$gene[low_grp] ~ dat$TF[ low_grp ] )
    lm_hq25 <- lm( dat$gene[ high_grp] ~ dat$TF[ high_grp])
    lm_middle <- lm(dat$gene[ mid_grp ] ~ dat$TF[ mid_grp ] )
    lm_all <- lm(dat$gene ~ dat$TF)
    abline(lm_hq25,col='red',lty='dashed')
    abline(lm_lq25,col='green',lty='dashed')
    abline(lm_middle,col='grey',lty='dashed')
    abline(lm_all,col='black',lty='dashed')

    legend('bottomright',
           col=c('red','green'),
           legend=c('High LncExp',
                    'Low LncExp'),
           pch=20
           )
    legend('bottomleft',
           col=c('red','green','grey','black'),
           legend=c('High LncExp',
                    'Low LncExp',
                    'Middle LncExp',
                    'All samples'),
           lty=rep(1,4)
           )

    dev.off()
}

##test it

##triplePlot using ggplot2

triPlot_ggplot2 <- function(m,e,tg,fname='tmp.pdf',xlab='TF expression',ylab='targetGene expression',
                   labels=NULL) {
if(!require(ggplot2)) stop('error loading ggplot2')
tmp=as.data.frame(cbind(`LncRNA(log2RPKM)`=m,`TF(log2RPKM)`=e,`Gene(log2RPKM)`=tg))
if(!is.null(labels)) {
tmp = cbind(tmp,labels=labels)
}
midpoint=median(tmp$`LncRNA(log2RPKM)`,na.rm=T)
pdf(fname)
gp=ggplot(tmp,aes(`TF(log2RPKM)`,`Gene(log2RPKM)`)) + geom_point(aes(colour=`LncRNA(log2RPKM)`)) + scale_colour_gradient2(low = "#d7191c",mid='#ffffbf',high='#1a9641',midpoint=midpoint) 
if(!is.null(labels)) {
gp = gp + geom_text(aes(label=tmp$labels))
}
print(gp)
dev.off()
}

##multiple genes
triPlot_ggplot2_mg <- function(m,e,tgs,fname='tmp.pdf',xlab='TF expression',ylab='targetGene expression',
                   labels=NULL,label.size=5,point.size=1,plot=T) {
if(!require(ggplot2)) stop('error loading ggplot2')
tmp =c()
for(i in 1:nrow(tgs)) {
tmpi=as.data.frame(cbind(`LncRNA(log2RPKM)`=m,`TF(log2RPKM)`=e,GeneName=rep(rownames(tgs)[i],length(m)),`Gene(log2RPKM)`=tgs[i,]))
tmp=rbind(tmp,tmpi)
}
#colnames(tmp) = c('Lnc','TF','GeneName','GeneExp')
tmp$`LncRNA(log2RPKM)`=as.numeric(as.character(tmp$`LncRNA(log2RPKM)`))
tmp$`TF(log2RPKM)`=as.numeric(as.character(tmp$`TF(log2RPKM)`))
tmp$`Gene(log2RPKM)`=as.numeric(as.character(tmp$`Gene(log2RPKM)`))
tmp$GeneName=as.factor(as.character(tmp$GeneName))
if(!is.null(labels)) {
tmp = cbind(tmp,labels=labels)
}
midpoint=median(tmp$`LncRNA(log2RPKM)`,na.rm=T)

gp=ggplot(tmp,aes(`TF(log2RPKM)`,`Gene(log2RPKM)`)) + geom_point(aes(colour=`LncRNA(log2RPKM)`),size=point.size) + scale_colour_gradient2(high= "#d7191c",mid='#ffffbf',low='#1a9641',midpoint=midpoint) 
gp=gp+facet_wrap(~GeneName)
if(!is.null(labels)) {
gp = gp + geom_text(aes(label=tmp$labels),size=label.size)
}
  if(!is.null(xlab)) gp = gp + xlab(xlab)
   if(!is.null(ylab)) gp = gp + ylab(ylab)
  if(plot) {
pdf(fname)
print(gp)
dev.off()
    }
  return(gp)
}
