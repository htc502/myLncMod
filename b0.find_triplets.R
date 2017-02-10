my.tri.app <-
function(ms,ET,M.exp,E.exp,T.exp,N = 0.25,method="pearson",iqr.filter = c(log2(1.5),log2(1.5),log2(1.5)),
                  cor.MvsET = c(0.3,0.3),cor.EvsT.dif = 0.45,cor.EvsT.abs = 0.4,
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
if(sum(low_grp) == 0 | sum(high_grp) == ) {
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

delta.cor <- PCChigh-PCClow
diff.res.index <- abs(delta.cor) > cor.EvsT.dif
abs.res.index <- abs(PCClow) > cor.EvsT.abs | abs(PCChigh) > cor.EvsT.abs
deltR.res.index <- diff.res.index & abs.res.index

if(sum(deltR.res.index) != 0) {
tmp <- data.frame(rep(M,sum(deltR.res.index)), rep(E,sum(deltR.res.index)),Ts[deltR.res.index], PCClow[deltR.res.index],
		PCChigh[deltR.res.index],rep(DE.res.index,sum(deltR.res.index)),rep(MvsE.res.index,sum(deltR.res.index)),
		MvsT.res.index[deltR.res.index])


n <- length(lncRNAexp)
    deltR <- tmp[,5]-tmp[,4]
    ##assess significant level
    
    ##make sure that the seed only function locally..
    if(exists('.Random.seed')) {
    old <- .Random.seed
    on.exit( { .Random.seed <<- old } )
    }
    set.seed(seed)
    deltR_nulls <- c()
	nlow <- sum(low_grp)
	nhigh <- sum(high_grp)
    for(i in 1:nrand) {
        null_grplow <- sample(1:n,nlow)
        null_grphigh <- sample(setdiff(1:n,null_grplow),nhigh)
      
        nullPCClow <- cor(t(mRNAsexp[as.character(tmp[,3]),null_grplow,drop=F]),TFexp[null_grplow],use="pairwise.complete.obs",method=method)
        nullPCClow <- nullPCClow[,1]
        nullPCChigh <- cor(t(mRNAsexp[as.character(tmp[,3]),null_grphigh,drop=F]),TFexp[null_grphigh],use="pairwise.complete.obs",method=method)
        nullPCChigh <- nullPCChigh[,1]
        deltR_null <- nullPCChigh-nullPCClow
        deltR_nulls <- cbind(deltR_nulls,deltR_null)
    }
    deltR_nulls <- cbind(deltR,deltR_nulls)
    p <- apply(deltR_nulls,1,function(r,nrand) {
	if(r[1] <= 0) {tmpp <- sum(r[-1]<r[1])/nrand}
	if(r[1] > 0) {tmpp <- sum(r[-1]>r[1])/nrand}
	tmpp
      },nrand=nrand)
	tmp <- cbind(tmp,p)
	colnames(tmp) <- c("modulator","effector","target","R_low","R_HIGH","DE","MvsE","MvsT","p-value")
    return(tmp)
	}
} else {
return(paste0(M,' ',E,' no triples with deltR threshold passed'))
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


####test it
##lnc1 <- rnorm(100,5,.4)
####lnc1_1 <- rnorm(100,5,.1)
####
##tf1 <- rnorm(100,3,.5)
##g1 <- tf1+rnorm(100,0.1,0.02)
##g1_1<- tf1+rnorm(100,4,2)
##e <- rbind(g1,g1_1)
##find_triplets(e,tf1,lnc1,nrand=1000)
##find_triplets(g1,tf1,lnc1,nrand=1000)
##find_triplets(g1_1,tf1,lnc1,nrand=1000)
##
