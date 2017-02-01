
find_triplets <- function(mRNAexp,
                          TFexp,
                          lncRNAexp,
                          nrand=100) {
    if(length(mRNAexp) != length(TFexp) | length(TFexp) != length(lncRNAexp)) stop('input of diff length')
    ord <- order(lncRNAexp, decreasing=F)
    n <- length(lncRNAexp)
    n_p25 <- round(n*.25)
    low_grp <- ord[ 1:n_p25 ]
    high_grp <- ord[ (n-n_p25+1):n ]

    ##check lncRNA TF independence
    fc <- abs(median(TFexp[high_grp])-median(TFexp[low_grp]))
    p <- t.test(TFexp[high_grp],TFexp[low_grp])$p.value
    if(2^fc > 1.5 & p < 0.05) return(rep(-1,4))
    PCClow <- cor(TFexp[low_grp],mRNAexp[low_grp])
    PCChigh <- cor(TFexp[high_grp],mRNAexp[high_grp])
    deltR <- PCChigh-PCClow
    ##assess significant level
    deltR_nulls <- c()
    for(i in 1:nrand) {
        null_samples <- sample(ord,2*n_p25)
        null_grplow <- null_samples[1:n_p25]
        null_grphigh <- null_samples[(n_p25+1):(2*n_p25)]
        nullPCClow <- cor(TFexp[null_grplow],mRNAexp[null_grplow])
        nullPCChigh <- cor(TFexp[null_grphigh],mRNAexp[null_grphigh])
        deltR_null <- abs(nullPCChigh-nullPCClow)
        deltR_nulls <- c(deltR_nulls,deltR_null)
    }
    p <- sum(deltR_nulls > abs(deltR))/nrand
    return(c(PCClow,PCChigh, PCChigh-PCClow, p))
}




####test it
##lnc1 <- rnorm(100,5,.4)
##lnc1_1 <- rnorm(100,5,.1)
##
##tf1 <- rnorm(100,3,.5)
##g1 <- tf1+rnorm(100,0.1,0.02)
##g1_1<- tf1+rnorm(100,4,2)
