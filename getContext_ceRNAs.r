get_context_ceRNAs = function(ceRNAs,tmpM.exp,tmpE.exp,cores=4) {
lncidx=ceRNAs[,1] %in% rownames(tmpM.exp)
mrnaidx=ceRNAs[,4] %in% tmpE.exp
idx=lncidx & mrnaidx
ceRNAs=ceRNAs[idx,c(1,4)]

lmf <- function(iceRNA) {
lncexp=tmpM.exp[ ceRNAs[iceRNA,1], ]
mrnaexp=tmpE.exp[ ceRNAs[iceRNA,2], ]

model=lm(lncexp ~ mrnaexp)
coefficients(summary(model)) -> coef
beta=coef[2,1]
p=coef[2,4]
res=c(ceRNAs[iceRNA,1],ceRNAs[iceRNA,2],beta,p)
res
}

library(parallel)
mclapply(1:nrow(ceRNAs),FUN=lmf,cores=cores) -> lmres
lmres1=do.call('rbind',lmres)
fdr=p.adjust(as.numeric(lmres[,4]),'BH')
lmres1=cbind(lmres1,fdr=fdr)
colnames(lmres1)=c('lncRNA','mRNA','beta','p','fdr')
lmres1
}