find_triplets <- function(mRNAsexp,
                          TFexp,
                          lncRNAexp,
                          nrand=100,
                          seed=123) {
    if(is.vector(mRNAsexp)){
        tmp <- as.data.frame(matrix(data=mRNAsexp,nrow=1,ncol=length(mRNAsexp)))
        mRNAsexp <- tmp
    }

    if(ncol(mRNAsexp) != length(TFexp) | length(TFexp) != length(lncRNAexp)) stop('input of diff length')
cutoffs <- quantile(lncRNAexp,c(.25,.75))
low_grp <- lncRNAexp < cutoffs[1]
high_grp <- lncRNAexp > cutoffs[2]
    PCClow <- cor( t(mRNAsexp[ , low_grp]),TFexp[ low_grp] )
    PCChigh <- cor( t(mRNAsexp[ , high_grp]),TFexp[ high_grp] )
    PCClow <- PCClow[,1];PCChigh <- PCChigh[,1]

    deltR <- PCChigh-PCClow
    ##assess significant level
    set.seed(seed)
    deltR_nulls <- c()
	nlow <- sum(low_grp)
	nhigh <- sum(high_grp)
    for(i in 1:nrand) {
        null_grplow <- sample(1:n,nlow)
        null_grphigh <- sample(setdiff(1:n,null_grplow),nhigh)
        nullPCClow <- cor(t(mRNAsexp[,null_grplow]),TFexp[null_grplow])
        nullPCClow <- nullPCClow[,1]
        nullPCChigh <- cor(t(mRNAsexp[,null_grphigh]),TFexp[null_grphigh])
        nullPCChigh <- nullPCChigh[,1]
        deltR_null <- nullPCChigh-nullPCClow
        deltR_nulls <- cbind(deltR_nulls,deltR_null)
    }
    deltR_nulls <- cbind(deltR,deltR_nulls)
    p <- apply(deltR_nulls,1,function(r,nrand) {
      tmp <- abs(r)
      p <- sum(tmp[-1] > tmp[1])/nrand
      p},nrand=nrand)
    return(cbind(PCClow,PCChigh, deltR, p))
}




####test it
##lnc1 <- rnorm(100,5,.4)
####lnc1_1 <- rnorm(100,5,.1)
####
##tf1 <- rnorm(100,3,.5)
##g1 <- tf1+rnorm(100,0.1,0.02)
##g1_1<- tf1+rnorm(100,4,2)
##e <- rbind(g1,g1_1)
##find_triplets1(e,tf1,lnc1,nrand=1000)
##find_triplets(g1,tf1,lnc1,nrand=1000)
##find_triplets(g1_1,tf1,lnc1,nrand=1000)
##
