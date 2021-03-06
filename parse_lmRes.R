parse_lmRes <- function(lmres, getswap=T) {
  ilmres <- lmres$good
  PMat <- do.call('rbind',lapply(ilmres,function(e) e$p))
  BetaMat <- do.call('rbind',lapply(ilmres, function(e) e$beta))
  tripleMat <- do.call('rbind',lapply(ilmres,function(e) e$triple))
  outlierMat <- do.call('rbind',lapply(ilmres,function(e) c(paste(e$outlier$TF$low,collapse = ';'),
                                                             paste(e$outlier$TF$mid,collapse = ';'),
                                                             paste(e$outlier$TF$high,collapse = ';'),
                                                             paste(e$outlier$mRNA$low,collapse = ';'),
                                                             paste(e$outlier$mRNA$mid,collapse = ';'),
                                                             paste(e$outlier$mRNA$high,collapse = ';'))))
  res <- list(PMat=PMat,
              BetaMat=BetaMat,
              tripleMat=tripleMat,
              outlierMat=outlierMat)
  res
}
                                       
 parse_lmRes1 <- function(lmres) {
  ilmres <- lmres$good
  PMat <- do.call('rbind',lapply(ilmres,function(e) e$p))
  BetaMat <- do.call('rbind',lapply(ilmres, function(e) e$beta))
  PMat_swap <- do.call('rbind',lapply(ilmres, function(e) e$p_swap))
  BetaMat_swap <- do.call('rbind',lapply(ilmres, function(e) e$beta_swap))
  tripleMat <- do.call('rbind',lapply(ilmres,function(e) e$triple))
  outlierMat <- do.call('rbind',lapply(ilmres,function(e) c(paste(e$outlier$TF$low,collapse = ';'),
                                                             paste(e$outlier$TF$mid,collapse = ';'),
                                                             paste(e$outlier$TF$high,collapse = ';'),
                                                             paste(e$outlier$mRNA$low,collapse = ';'),
                                                             paste(e$outlier$mRNA$mid,collapse = ';'),
                                                             paste(e$outlier$mRNA$high,collapse = ';'))))
  NgrpSample <- do.call('rbind',lapply(ilmres, function(e) {tmp = unlist(e$groupInfo$ngrp);names(tmp) = names(e$groupInfo$ngrp);tmp}))
  res <- list(PMat=PMat,
              PMat_swap=PMat_swap,
              BetaMat=BetaMat,
              BetaMat_swap=BetaMat_swap,
              tripleMat=tripleMat,
              outlierMat=outlierMat,
             NgrpSample = NgrpSample)
  res
}
