my.tri.app.lm.swap.grp <-
    function(MET,M.exp,E.exp,T.exp,N = 0.25,groupMin=20,
             correction="BH",p.cutoff=0.01, cores=1){
      data.E = E.exp
      data.M = M.exp
      data.T = T.exp

        index.MET.M<-as.character(MET[,1])%in%rownames(data.M)
        index.MET.E<-as.character(MET[,2])%in%rownames(data.E)
        index.MET.T<-as.character(MET[,3])%in%rownames(data.T)
        index.MET<-index.MET.M&index.MET.E&index.MET.T
        MET<-MET[index.MET,]
        data.M <- as.matrix(data.M)
        data.T <- as.matrix(data.T)
        data.E <- as.matrix(data.E)
        if(dim(MET)[1]==0){
            return('no triples,return empty result')
        } else {

            find_triplets <- function(i) {
                M <- MET[i,1]

                E <- MET[i,2]
                Tg <- MET[i,3]
                 if (!(M %in% rownames(data.M)) | 
                     !(E %in% rownames(data.E)) |
                     !(Tg %in% rownames(data.E))) {
                    return(paste0(M,' ',E,' ',Tg,': not found in filtered expMtrix'))
                }

                mRNAexp <- data.T[ Tg,]
                TFexp <- data.E[ E, ]
                lncRNAexp <- data.M[ M,]

                if(length(mRNAexp) != length(TFexp) | length(TFexp) != length(lncRNAexp)) stop('input of diff length')
                cutoffs <- quantile(lncRNAexp,na.rm=T,probs=c(.25,.75))
                low_grp <- lncRNAexp < cutoffs[1]
                mid_grp <- lncRNAexp >= cutoffs[1] & lncRNAexp <= cutoffs[2]
                high_grp <- lncRNAexp > cutoffs[2]
                grp <- rep('mid',length=length(mRNAexp))
                grp[low_grp] <- 'low'
                grp[high_grp] <- 'high'
                grp <- as.factor(grp) #used for lm

                ##we swap the high and low group to test whether the triple is significant in low group
                grp <- factor(grp,levels = c('low','mid','high'))
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
                if(!length(c(tf_lowgrp_outlierIdx,tf_highgrp_outlierIdx,tf_midgrp_outlierIdx))==0) {
                TFexp[ c(tf_lowgrp_outlierIdx,tf_highgrp_outlierIdx,tf_midgrp_outlierIdx) ] <- NA
                }
                
                if(!length(c(mRNA_lowgrp_outlierIdx,mRNA_highgrp_outlierIdx,mRNA_midgrp_outlierIdx))==0) {
                mRNAexp[ c(mRNA_lowgrp_outlierIdx,mRNA_highgrp_outlierIdx,mRNA_midgrp_outlierIdx) ] <- NA
                }
                ## number of available values in cleaned data within each group:
                nlowgrp <- sum(!is.na(TFexp[low_grp] + mRNAexp[low_grp]))
                nhighgrp <- sum(!is.na(TFexp[high_grp] + mRNAexp[high_grp]))
                nmidgrp <- sum(!is.na(TFexp[mid_grp] + mRNAexp[mid_grp]))

                if(nlowgrp < groupMin |
                   nhighgrp < groupMin |
                   nmidgrp < groupMin) {
                    return(paste0(M,' ',E,' ',Tg,': high/mid/low group size < predefined groupMin(default:20) after remove outlier'))
                }

                model <- lm( mRNAexp ~ TFexp + grp + TFexp:grp)
                model_coef <- coefficients(summary(model))
                if(!(nrow(model_coef)==6 & ncol(model_coef)==4)) {
                    return(paste0(M,' ',E,' ',Tg,': the model dimension is not proper, some group may have empty values'))
                }
                p <- model_coef[,4]
                beta <- model_coef[,1]
                names(p) <- names(beta) <- rownames(model_coef)
                resi <- list(triple=c(M,E,Tg),beta=beta,p=p,
                             outlier=list(TF=list(low=tf_lowgrp_outlierIdx,high=tf_highgrp_outlierIdx,mid=tf_midgrp_outlierIdx),
                                          mRNA=list(low=mRNA_lowgrp_outlierIdx,high=mRNA_highgrp_outlierIdx,mid=mRNA_midgrp_outlierIdx)))
                return(resi)
            }
            library(parallel)
            lmFits <- mclapply(1:nrow(MET),find_triplets,mc.cores=cores)
           ## lmFits <- lapply(1:nrow(MET),find_triplets) ##for debug find_triplets
            badME.idx <- unlist(lapply(lmFits, function(e) is.character(e)))
            badres <- unlist(lmFits[ badME.idx] );badres <- cbind(MET[ badME.idx,],badres)
            goodres <- lmFits[ !badME.idx ]
            res <- list(bad = badres,
                        good = goodres)
            return(res)
        }
    }
