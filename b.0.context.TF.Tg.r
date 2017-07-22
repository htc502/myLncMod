#' @input: exp.rda
#' @output: TF_Tg.rda

getContext_TFtarget <- function(tmpE.exp, tf_target,cores=6) {
  mytf_target2 = tf_target
  tmpT.exp <- tmpE.exp
  tmpET0 <- as.matrix(mytf_target2)
  mode(tmpET0) <- 'character'
  idx1 <- tmpET0[, 1] %in% rownames(tmpE.exp)
  idx2 <- tmpET0[, 2] %in% rownames(tmpE.exp)
  idx3 <- tmpET0[,1] != tmpET0[,2]
  tmpET <- tmpET0[idx1 & idx2 & idx3,]
  print(paste0('# TFs in expressiodat:', length(unique(as.character(tmpET[,1])))))
  print(paste0('# Targets in expressiodat:', length(unique(as.character(tmpET[,2])))))
##  lm_e_t <- function(tf, target) {
  lm_e_t <- function(i) {
tf=tmpET[i,1]
target=tmpET[i,2]
    exp_tf <- tmpE.exp[tf, ]
    exp_mrna <- tmpE.exp[target, ]
    lm.tmp <- lm(exp_mrna ~ exp_tf)
    res.tmp <- c(
      tf,
      target,
      summary(lm.tmp)$r.squared,
      summary(lm.tmp)$coefficients[2, 1],
      summary(lm.tmp)$coefficients[2, 4]
    )
    res.tmp
  }
  ##apply(tmpET, 1, function(et) {
library(parallel)

  mclapply(1:nrow(tmpET), lm_e_t, mc.cores=cores) -> lmres
  lmres1 <- do.call('rbind',lmres)
  padj <- p.adjust(as.numeric(lmres1[, 5]), 'BH')
  lmres1 <- cbind(lmres1, fdr = padj)
  colnames(lmres1) <- c('tf','target','r2','beta','pvalue','fdr')
  lmres1
}
