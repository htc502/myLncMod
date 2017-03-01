parse_lmRes <- function(lmres) {
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