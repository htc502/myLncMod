getContext_TFtarget <- function(tmpE.exp, tf_target) {
  mytf_target2 = tf_target
  tmpT.exp <- tmpE.exp
  tmpET0 <- as.data.frame(mytf_target2)
  idx1 <- tmpET0[, 1] %in% rownames(tmpE.exp)
  idx2 <- tmpET0[, 2] %in% rownames(tmpE.exp)
  tmpET <- tmpET0[idx1 & idx2,]
  print(paste0('# TFs in expressiodat:', length(unique(as.character(tmpET[,1])))))
  print(paste0('# Targets in expressiodat:', length(unique(as.character(tmpET[,2])))))
  lm_e_t <- function(tf, target) {
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
  apply(tmpET, 1, function(et) {
    tf <- as.character(et[1])
    mrna <- as.character(et[2])
    tmpres <- lm_e_t(tf, mrna)
    tmpres
  }) -> lmres
  lmres <- t(lmres)
  padj <- p.adjust(as.numeric(lmres[, 4]), 'BH')
  lmres <- cbind(lmres, fdr = padj)
  colnames(lmres) <- c('tf','target','r2','beta','pvalue','fdr')
  lmres
}
