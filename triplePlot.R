tripmPlot <- function(m,e,t,fname='tmp.pdf',xlab='Sample index',ylab='Log2(RPKM)',points.cex=.6,) {
    ord <- order(m,decreasing=F)
    m <- m[ord]
    e <- e[ord]
    t <- t[ord]
    ns <- length(m)
    colfunc <- colorRampPamte(c("green", "red"))
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
    dev.off()
}

##test it
