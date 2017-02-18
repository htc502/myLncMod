my.tri.app1 <-
    function(ms,ET,M.exp,E.exp,T.exp,N = 0.25,method="pearson",iqr.filter = c(log2(1.5),log2(1.5),log2(1.5)),
             cor.MvsET = c(0.3,0.3),delta.p.cutoff=0.01,
             ET.fc.filter = log2(1.5),ET.p.filter = 0.01,nrand=100,seed=123,correction="BH",cores=1){

        iqr <- apply(M.exp, 1, IQR, na.rm=TRUE)
        data.M <- M.exp[iqr>iqr.filter[1],]
        iqr <- apply(E.exp, 1, IQR, na.rm=TRUE)
        data.E <- E.exp[iqr>iqr.filter[2],]
        iqr <- apply(T.exp, 1, IQR, na.rm=TRUE)
        data.T <- T.exp[iqr>iqr.filter[3],]
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
            M_E <- expand.grid(ms,E2T_E)

            find_triplets <- function(iM_E) {
              method = 'pearson' ##for the reason, see below the F trans part
                M <- as.character(M_E[iM_E,1])
                if (M %in% rownames(data.M)) {
                    E <- as.character(M_E[iM_E,2])
                    pos <- which(E2T_E == E)
                    Ts <- E2T[[pos]];rm(pos)
                    mRNAsexp <- data.T[ Ts, , drop=F]
                    TFexp <- data.E[ E, ]
                    lncRNAexp <- data.M[ M,]

                    if(ncol(mRNAsexp) != length(TFexp) | length(TFexp) != length(lncRNAexp)) stop('input of diff length')
                    cutoffs <- quantile(lncRNAexp,c(.25,.75))
                    low_grp <- lncRNAexp < cutoffs[1]
                    high_grp <- lncRNAexp > cutoffs[2]
                    if(sum(low_grp) == 0 | sum(high_grp) == 0) {
                        return(paste0(M,' ',E,' high group or low group empty, probabaly modulator expression is extremely biased, check them please, there may be extra bonus'))
                    } else {
                        PCClow <- cor( t(mRNAsexp[ , low_grp,drop=F]),TFexp[ low_grp],use="pairwise.complete.obs",method=method )
                        PCChigh <- cor( t(mRNAsexp[ , high_grp,drop=F]),TFexp[ high_grp] ,use="pairwise.complete.obs",method=method)
                        PCClow <- PCClow[,1];PCChigh <- PCChigh[,1]
                        
                        ##test M E T independence
                        DE <- t.test(TFexp[ low_grp] ,TFexp[ high_grp ])
                        DE.p <- DE$p.value
                        DE.FC <- abs(DE$estimate[1]-DE$estimate[2])
                        DE.res.index <- !(DE.p < ET.p.filter & DE.FC > ET.fc.filter)
                        MvsE.cor <- cor(lncRNAexp,TFexp,use="pairwise.complete.obs",method=method)
                        MvsT.cor <- cor(t(mRNAsexp),lncRNAexp,use="pairwise.complete.obs",method=method)
                        MvsE.res.index <- MvsE.cor < cor.MvsET[1]
                        MvsT.res.index <- MvsT.cor < cor.MvsET[2]

                        ##here PCClow or PCChigh can be NA when there is no variance in TF/mRNA expression of high/low group
                        ##I want to put PCC = 0 for these conditions to make them as failed triples
                        ## when using mclapply, these failed triples may affect others in the same computor core
                        
                        PCClow[is.na(PCClow)] <- 0
                        PCChigh[is.na(PCChigh)] <- 0
                        ##according to Wang Q et al (bioinformatics, 2015), when we use PCC, the following transformed version of PCC will follow normal distribution,
                        ##which means we don't need permutation, that's wanderful
                        PCChigh[ PCChigh == 1] <- .999
                        PCClow[ PCClow == 1] <- .999
                        ftrans <- function(pcc) { 0.2*(log((1+pcc)/(1-pcc),base = exp(1))) }
                        if(length(high_grp) <= 3 | length(low_grp) <= 3) {
                            return(paste0(M,E,'samples in high_grp/low_grp <= 3, unable to do F transformation'))
                        }
                        
                        delta.X <- (ftrans(PCChigh)-ftrans(PCClow))/(sqrt(1/(length(high_grp)-3) + 1/(length(low_grp)-3)))
                        
                        delta.p <- 2*pnorm(-abs(delta.X)) ##we want two-side test for both PCChigh > PCClow as well as PCClow > PCChigh
                        
                        deltR.res.index <- delta.p < delta.p.cutoff

                        if(sum(deltR.res.index) != 0) {
                            tmp <- data.frame(rep(M,sum(deltR.res.index)), rep(E,sum(deltR.res.index)),Ts[deltR.res.index], PCClow[deltR.res.index],
                                              PCChigh[deltR.res.index],rep(DE.res.index,sum(deltR.res.index)),rep(MvsE.res.index,sum(deltR.res.index)),
                                              MvsT.res.index[deltR.res.index])

                            delta.p1 <- delta.p[deltR.res.index]
                            tmp <- cbind(tmp,delta.p1)
                            colnames(tmp) <- c("modulator","effector","target","R_low","R_HIGH","DE","MvsE","MvsT","p-value")
                            return(tmp)
                        } else {
                            return(paste0(M,' ',E,' no triples with deltR threshold passed'))
                        }
                    }
                } else {
                    return(paste0(M,' not found in filtered expMtrix'))
                }
            }
            tmpMET <- mclapply(1:nrow(M_E),find_triplets,mc.cores=cores)
            badME.idx <- unlist(lapply(tmpMET, function(e) is.character(e)))
            badres <- unlist(tmpMET[ badME.idx] );badres <- cbind(M_E[ badME.idx,],badres)
            goodres <- do.call('rbind',tmpMET[ !badME.idx ])
            rownames(goodres) <- NULL
            res <- list(bad = badres,
                        good = goodres)
            return(res)
        }
    }
