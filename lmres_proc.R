lmres_postproc <- function(filedir,fdrcutoff=.25,RNAseqdata,swap_grp_p_cutoff=.005) {
    setwd(filedir)
    resFnames = list.files(path=filedir,pattern='lmres.*.rda')
    if(length(resFnames) == 0) stop('data files not found')
    stats = c()
    ntotal = 0
    nbad = 0
    ngood = 0
    nsig = 0
    for(i in 1:length(resFnames)) {
        load(resFnames[i])
        tmplmres = lmres;rm(lmres)
        ibad =  nrow(tmplmres$bad)
        igood = length(tmplmres$good)
        isig = length(tmplmres$sig)
        stats=rbind(stats,c(ibad,igood,isig))
        nbad = nbad + ibad
        ngood = ngood + igood
        nsig = nsig + isig
        rm(tmplmres)
    }
    stats=rbind(stats,c(nbad,ngood,nsig))
    colnames(stats) = c('#triples failed QC','#triples passed QC','#triples P<p.cutoff per run')
    write.csv(stats,file='Triple summary.csv')

    PMat <- c()
    BetaMat <- c()
    TripleMat <- c()
    outlierMat <- c()
    for(i in 1:length(resFnames) ) {
        print(i)
        load(resFnames[i])
        tmplmres = lmres;rm(lmres)
        ilmres <- tmplmres$good
        Pvaluei <- do.call('rbind',lapply(ilmres,function(e) e$p))
        betai <- do.call('rbind',lapply(ilmres, function(e) e$beta))
        triplei <- do.call('rbind',lapply(ilmres,function(e) e$triple))
        PMat <- rbind(PMat,Pvaluei)
        BetaMat <- rbind(BetaMat,betai)
        TripleMat <- rbind(TripleMat,triplei)
        outlieri <- do.call('rbind',lapply(ilmres,function(e) c(paste(e$outlier$TF$low,collapse = ';'),
                                                                paste(e$outlier$TF$mid,collapse = ';'),
                                                                paste(e$outlier$TF$high,collapse = ';'),
                                                                paste(e$outlier$mRNA$low,collapse = ';'),
                                                                paste(e$outlier$mRNA$mid,collapse = ';'),
                                                                paste(e$outlier$mRNA$high,collapse = ';'))))
        outlierMat <- rbind(outlierMat,outlieri)
        rm(tmplmres)
    }

    colnames(TripleMat) <- c('lncRNA','TF','mRNA')
    colnames(outlierMat) <- c('TF_low','TF_mid','TF_high','mRNA_low','mRNA_mid','mRNA_high')

    save(PMat,BetaMat,TripleMat,outlierMat,file='triples_pastQC.rda')
    fdr <- p.adjust(PMat[,5],'BH')
    idx <- fdr < 0.25 ##or .25
    PMat1 = PMat[ idx, ]
    BetaMat1 = BetaMat[ idx, ]
    TripleMat1 = TripleMat[ idx, ]
    outlier1 = outlierMat[ idx, ]
    ## in order to classify the triples into six categories, we need to test whether it is significant for TF and target in low grp,

    load(RNAseqdata)
    source('/data/ghan/AD_GBM_comparison/myLncMod/b0.find_triplets_lm_swap_highlow_grp.R', encoding = 'UTF-8')
    swapres <- my.tri.app.lm.swap.grp(TripleMat1,M.exp = tmpM.exp,E.exp = tmpE.exp,T.exp=tmpE.exp ,cores=4)
    source('/data/ghan/AD_GBM_comparison/myLncMod/parse_lmRes.R', encoding = 'UTF-8')
    parse_lmRes(swapres) -> swapres1
    save(swapres,swapres1,file='swith_high_low_group_lmres.rda')

    ##we have to pick up TF-target significant in lowgrp
    swapres1$PMat[,2] < 0.005 -> idx1
    sigSwapres <- list(PMat=subset(swapres1$PMat,idx1),
                       BetaMat=subset(swapres1$BetaMat, idx1),
                       tripleMat=subset(swapres1$tripleMat,idx1),
                       outlierMat=subset(swapres1$outlierMat,idx1))

    match(paste(sigSwapres$tripleMat[,1],sigSwapres$tripleMat[,2],sigSwapres$tripleMat[,3],sep='_'),
          paste(TripleMat1[,1],TripleMat1[,2],TripleMat1[,3],sep='_')) -> pos

    sigOrigin <- list(PMat=PMat1[ pos, ],
                      BetaMat=BetaMat1[ pos, ],
                      tripleMat=TripleMat1[ pos,] ,
                      outlierMat=outlier1[ pos,])
    ##sigOrigin ##high grp as reference
    sigSwap <- sigSwapres ##low grp as reference
    ##classify six categories
    ##category 1: enhanced inhibition
    ## low  high
    ##  - -> --
    betaLow <- sigSwap$BetaMat[,2]
    betaHigh <- sigOrigin$BetaMat[,2]

    idx.cat1 <- betaLow < 0 & betaHigh < 0 & betaLow > betaHigh
    ##categoy 2: attenuates inhibition
    ## -- -> -
    idx.cat2 <- betaLow < 0 & betaHigh < 0 & betaLow < betaHigh
    ## + -> ++
    idx.cat3 <- betaLow > 0 & betaHigh > 0 & betaLow < betaHigh
    ## ++ -> +
    idx.cat4 <- betaLow > 0 & betaHigh > 0 & betaLow > betaHigh
    ## - -> +
    idx.cat5 <- betaLow < 0 & betaHigh > 0
    ## + -> -
    idx.cat6 <- betaLow > 0 & betaHigh < 0

    ##classified result
    res <-list(c1_ei=cbind(sigOrigin$tripleMat[idx.cat1,,drop=F],betaLow[idx.cat1],betaHigh[idx.cat1]),
               c2_ai=cbind(sigOrigin$tripleMat[idx.cat2,,drop=F],betaLow[idx.cat2],betaHigh[idx.cat2]),
               c3_ea=cbind(sigOrigin$tripleMat[idx.cat3,,drop=F],betaLow[idx.cat3],betaHigh[idx.cat3]),
               c4_aa=cbind(sigOrigin$tripleMat[idx.cat4,,drop=F],betaLow[idx.cat4],betaHigh[idx.cat4]),
               c5_ia=cbind(sigOrigin$tripleMat[idx.cat5,,drop=F],betaLow[idx.cat5],betaHigh[idx.cat5]),
               c6_ii=cbind(sigOrigin$tripleMat[idx.cat6,,drop=F],betaLow[idx.cat6],betaHigh[idx.cat6]))
    write.csv(res$c1_ei,file='c1_ei.csv')
    write.csv(res$c2_ai,file='c2_ai.csv')
    write.csv(res$c3_ea,file='c3_ea.csv')
    write.csv(res$c4_aa,file='c4_aa.csv')
    write.csv(res$c5_ia,file='c5_ia.csv')
    write.csv(res$c6_ii,file='c6_ii.csv')
    save(sigOrigin$PMat,
         sigOrigin$BetaMat,
         sigOrigin$tripleMat,
         sigOrigin$outlierMat,file='significant_triples.rda')
}
