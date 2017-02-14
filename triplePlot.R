triPlot <- function(m,e,t,fname='tmp.pdf',xlab='Sample index',ylab='Log2(RPKM)',points.cex=.6) {
    ord <- order(m,decreasing=F)
    m <- m[ord]
    e <- e[ord]
    t <- t[ord]
    ns <- length(m)
    colfunc <- colorRampPalette(c("green", "red"))
    dat <- data.frame(x=1:length(m),lncRNA=m,TF=e,gene=t,cols=colfunc(ns))
    ylim <- c(min(dat[,2:4]),max(dat[,2:4]))
    pdf(fname)
    plot(dat$x,dat$lncRNA,pch=20,ylim=ylim,
         xlab=xlab,
         ylab=ylab,cex=points.cex)
    for(i in 1:nrow(dat)) {
        tf_x <- dat[i,1];tf_y <- dat[i,3]
        g_x <- dat[i,1];g_y <- dat[i,4]
        if(tf_y > g_y ) {
            col <- 'red'
        } else {
            col <- 'green'
        }

        segments(tf_x,tf_y,
                 g_x,g_y,
                 col=col)
    }

    points(dat[,1],dat[,3],pch=6,cex=points.cex)
    points(dat[,1],dat[,4],pch=5,cex=points.cex)
    legend('bottomright',
           legend=c('TF',
                    'Gene',
                    'LncRNA'),
           pch=c(6,5,20)
           )
    q25 <- quantile(m,c(.25,.75))
    n25_l <-sum(m < q25[1])
    n25_h <-sum(m > q25[2])
    n50 <- length(m)-n25_h-n25_l

    cols <- c(rep('green',n25_l),rep('grey',length(m)-n25_h-n25_l),rep('red',n25_h))
    plot(dat$TF,dat$gene,col=cols,pch=20,cex=points.cex,
         xlab='TF expression',
         ylab='Gene expression')

    lm_lq25 <- lm( dat$gene[ m < q25[1]] ~ dat$TF[ m < q25[1]] )
    lm_hq25 <- lm( dat$gene[ m > q25[2]] ~ dat$TF[ m > q25[2]])
    lm_middle <- lm(dat$gene[ m >= q25[1] & m < q25[2] ] ~ dat$TF[ m>= q25[1] & m <= q25[2] ] )
    lm_all <- lm(dat$gene ~ dat$TF)
    abline(lm_hq25,col='red',lty='dashed')
    abline(lm_lq25,col='green',lty='dashed')
    abline(lm_middle,col='grey',lty='dashed')
    abline(lm_all,col='black',lty='dashed')

    legend('bottomright',
           col=c('red','green'),
           legend=c('High LncExp',
                    'Low LncExp'),
           pch=20
           )
    legend('bottomleft',
           col=c('red','green','grey','black'),
           legend=c('High LncExp',
                    'Low LncExp',
                    'Middle LncExp',
                    'All samples'),
           lty=rep(1,4)
           )

    dev.off()
}

##test it
