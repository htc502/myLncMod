#' @input: triple.rda
#' @output: triple1.rda
lmres_postproc <- function(filedir,fdrcutoff=.25,RNAseqdata,Ngrp=.25) {
    setwd(filedir)
    resFnames = list.files(path=filedir,pattern='lmres_.*.rda')
    if(length(resFnames) == 0) stop('data files not found')
    stats = c()
    ntotal = 0
    nbad = 0
    ngood = 0
    nsig = 0
	print('get Triple summary table')
    for(i in 1:length(resFnames)) {
	print(i)
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
	print('combine result')
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
    idx <- fdr < fdrcutoff ##or .25
    PMat1 = PMat[ idx, ]
    BetaMat1 = BetaMat[ idx, ]
    TripleMat1 = TripleMat[ idx, ]
    outlier1 = outlierMat[ idx, ]
    ## in order to classify the triples into six categories, we need to test whether it is significant for TF and target in low grp,

	print('calculate beta in low group')
    load(RNAseqdata)
    source('/data/ghan/AD_GBM_comparison/AD_RNAseq_AMP_AD/get_rpkm/myLncMod/b0.find_triplets_lm_swap_highlow_grp.R', encoding = 'UTF-8')
    swapres <- my.tri.app.lm.swap.grp(TripleMat1,M.exp = tmpM.exp,E.exp = tmpE.exp,T.exp=tmpE.exp ,cores=4,Ngrp=Ngrp)
    source('/data/ghan/AD_GBM_comparison/AD_RNAseq_AMP_AD/get_rpkm/myLncMod/parse_lmRes.R', encoding = 'UTF-8')
    parse_lmRes(swapres) -> swapres1
    save(swapres,swapres1,file='swith_high_low_group_lmres.rda')

    lmRes_ref_higgrp <- list(PMat=PMat1,
                      BetaMat=BetaMat1,
                      tripleMat=TripleMat1 ,
                      outlierMat=outlier1)
    lmRes_ref_lowgrp <- swapres1
    ##classify six categories
    ##category 1: enhanced inhibition
    ## low  high
    ##  - -> --

	print('identify six categoires')
    betaLow <- lmRes_ref_lowgrp$BetaMat[,2]
    betaHigh <- lmRes_ref_higgrp$BetaMat[,2]

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
    res0 <-list(c1_ei=cbind(lmRes_ref_higgrp$tripleMat[idx.cat1,,drop=F],betaLow=betaLow[idx.cat1],betaLowP=lmRes_ref_lowgrp$PMat[idx.cat1,2,drop=F],betaHigh=betaHigh[idx.cat1],betaHighP=lmRes_ref_higgrp$PMat[idx.cat1,2,drop=F]),
               c2_ai=cbind(lmRes_ref_higgrp$tripleMat[idx.cat2,,drop=F],betaLow=betaLow[idx.cat2],betaLowP=lmRes_ref_lowgrp$PMat[idx.cat2,2,drop=F],betaHigh=betaHigh[idx.cat2],betaHighP=lmRes_ref_higgrp$PMat[idx.cat2,2,drop=F]),
               c3_ea=cbind(lmRes_ref_higgrp$tripleMat[idx.cat3,,drop=F],betaLow=betaLow[idx.cat3],betaLowP=lmRes_ref_lowgrp$PMat[idx.cat3,2,drop=F],betaHigh=betaHigh[idx.cat3],betaHighP=lmRes_ref_higgrp$PMat[idx.cat3,2,drop=F]),
               c4_aa=cbind(lmRes_ref_higgrp$tripleMat[idx.cat4,,drop=F],betaLow=betaLow[idx.cat4],betaLowP=lmRes_ref_lowgrp$PMat[idx.cat4,2,drop=F],betaHigh=betaHigh[idx.cat4],betaHighP=lmRes_ref_higgrp$PMat[idx.cat4,2,drop=F]),
               c5_ia=cbind(lmRes_ref_higgrp$tripleMat[idx.cat5,,drop=F],betaLow=betaLow[idx.cat5],betaLowP=lmRes_ref_lowgrp$PMat[idx.cat5,2,drop=F],betaHigh=betaHigh[idx.cat5],betaHighP=lmRes_ref_higgrp$PMat[idx.cat5,2,drop=F]),
               c6_ii=cbind(lmRes_ref_higgrp$tripleMat[idx.cat6,,drop=F],betaLow=betaLow[idx.cat6],betaLowP=lmRes_ref_lowgrp$PMat[idx.cat6,2,drop=F],betaHigh=betaHigh[idx.cat6],betaHighP=lmRes_ref_higgrp$PMat[idx.cat6,2,drop=F]))
    res <- lapply(res0, function(e) {colnames(e) <- c("lncRNA","TF","mRNA","betaLow","p","betaHigh","p");e})
	print('output six categoires')
    write.csv(res$c1_ei,file='c1_ei.csv')
    write.csv(res$c2_ai,file='c2_ai.csv')
    write.csv(res$c3_ea,file='c3_ea.csv')
    write.csv(res$c4_aa,file='c4_aa.csv')
    write.csv(res$c5_ia,file='c5_ia.csv')
    write.csv(res$c6_ii,file='c6_ii.csv')
    save(lmRes_ref_higgrp,lmRes_ref_lowgrp,file='significant_triples.rda')
}
