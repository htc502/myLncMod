my.tri.app.lm <-
    function(ms,ET,M.exp,E.exp,T.exp,N = 0.25,groupMin=20,
             correction="BH",cores=1){

        index.ET.E<-as.character(ET[,1])%in%rownames(data.E)
        index.ET.T<-as.character(ET[,2])%in%rownames(data.T)
        index.ET<-index.ET.E&index.ET.T
        ET<-ET[index.ET,]
        data.M <- as.matrix(data.M)
        data.T <- as.matrix(data.T)
        data.E <- as.matrix(data.E)
        if(dim(ET)[1]==0){
            return('no effector-Target pairs,return empty result')
        } else {

            E2T <- tapply(as.character(ET[,2]),as.character(ET[,1]),function(e) unique(e))
            E2T_E <- names(E2T)
            ##M_E <- expand.grid(ms,E2T_E)
            MET <- expand.grid(1:length(ms),1:nrow(ET))

            find_triplets <- function(i) {
                M <- as.character(ms[MET[i,1]])
                if (!(M %in% rownames(data.M)) ) {
                    return(paste0(M,' not found in filtered expMtrix'))
                }

                E <- as.character(ET[MET[i,2],1])
                Tg <- as.character(ET[MET[i,2],2])
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

                ##remove outliers within each group
                TF_extremValue <- tapply(TFexp, grp, function(e) {
                    IQR(e,na.rm=T)->iqr
                    q13 <- quantile(e,na.rm=T,probs=c(.25,.75))
                    list(lower=q13[1]-1.5*iqr,
                         upper=q13[2]+1.5*iqr)
                })
                which(TFexp > TF_extremValue$low$upper)[low_grp] -> tf_lowgrp_outlier1
                which(TFexp < TF_extremValue$low$lower)[low_grp] -> tf_lowgrp_outlier2
                tf_lowgrp_outlierIdx <- c(tf_lowgrp_outlier1,tf_lowgrp_outlier2)
                which(TFexp > TF_extremValue$mid$upper)[mid_grp] -> tf_midgrp_outlier1
                which(TFexp < TF_extremValue$mid$lower)[mid_grp] -> tf_midgrp_outlier2
                tf_midgrp_outlierIdx <- c(tf_midgrp_outlier1,tf_midgrp_outlier2)
                which(TFexp > TF_extremValue$high$upper)[high_grp] -> tf_highgrp_outlier1
                which(TFexp < TF_extremValue$high$lower)[high_grp] -> tf_highgrp_outlier2
                tf_highgrp_outlierIdx <- c(tf_highgrp_outlier1,tf_highgrp_outlier2)

                mRNA_extremValue <- tapply(mRNAexp, grp, function(e) {
                    IQR(e,na.rm=T)->iqr
                    q13 <- quantile(e,na.rm=T,probs=c(.25,.75))
                    list(lower=q13[1]-1.5*iqr,
                         upper=q13[2]+1.5*iqr)
                })
                which(mRNAexp > mRNA_extremValue$low$upper)[low_grp] -> mrna_lowgrp_outlier1
                which(mRNAexp < mRNA_extremValue$low$lower)[low_grp] -> mrna_lowgrp_outlier2
                mRNA_lowgrp_outlierIdx <- c(mRNA_lowgrp_outlier1,mRNA_lowgrp_outlier2)
                which(mRNAexp > mRNA_extremValue$mid$upper)[mid_grp] -> mrna_midgrp_outlier1
                which(mRNAexp < mRNA_extremValue$mid$lower)[mid_grp] -> mrna_midgrp_outlier2
                mRNA_midgrp_outlierIdx <- c(mRNA_midgrp_outlier1,mRNA_midgrp_outlier2)
                which(mRNAexp > mRNA_extremValue$high$upper)[high_grp] -> mrna_highgrp_outlier1
                which(mRNAexp < mRNA_extremValue$high$lower)[high_grp] -> mrna_highgrp_outlier2
                mRNA_highgrp_outlierIdx <- c(mRNA_highgrp_outlier1,mRNA_highgrp_outlier2)

                ##set outliers'expression value to NA
                TFexp[ c(tf_lowgrp_outlierIdx,tf_highgrp_outlierIdx,tf_midgrp_outlierIdx) ] <- NA
                mRNAexp[ c(mRNA_lowgrp_outlierIdx,mRNA_highgrp_outlierIdx,mRNA_midgrp_outlierIdx) ] <- NA
                ## number of available values in cleaned data within each group:
                nlowgrp <- sum(!is.na(TF_exp[low_grp] + mRNAexp[low_grp]))
                nhighgrp <- sum(!is.na(TF_exp[high_grp] + mRNAexp[high_grp]))
                nmidgrp <- sum(!is.na(TF_exp[mid_grp] + mRNAexp[mid_grp]))

                if(nlowgrp < groupMin |
                   nhighgrp < groupMin |
                   nmidgrp <- groupMin) {
                    return(paste0(M,' ',E,' ',Tg,': high/mid/low group size < predefined groupMin(default:20) after remove outlier'))
                }

                model <- lm( mRNAexp ~ TFexp + grp)
                model
            }
            lmFits <- mclapply(1:nrow(MET),find_triplets,mc.cores=cores)
            badME.idx <- unlist(lapply(lmFIts, function(e) is.character(e)))
            badres <- unlist(lmFIts[ badME.idx] );badres <- cbind(M_E[ badME.idx,],badres)
            goodres <- lmFIts[ !badME.idx ]
            res <- list(bad = badres,
                        good = goodres)
            return(res)
        }
    }
