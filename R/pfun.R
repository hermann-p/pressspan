ps.read <- function(fileName) {
    read.csv(file=fileName, sep="\t")
}


ps.plot <- function(data, n, ptitle="") {
    if (ptitle == "") {
        ptitle <- names(data)[n]
    }
    data <- as.data.frame(data)
    tick.num <- 10
    bucket.size <- data[1,n]
    buckets <- data[-1,n]
    len <- length(buckets)
    my.ticks <- seq(0, len, len=tick.num)
    my.labels <- as.integer(my.ticks * bucket.size)
    par(mar=c(6, 4.1, 4.1, 2.1))
    plot(buckets, main=ptitle, axes=FALSE, ylab="events", xlab="", t="l")
    axis(side=2)
    axis(side=1, at=my.ticks, labels=my.labels, las=2)
    mtext("Position", side=1, at=c(0,0), adj=1, line=2)
}


bc.balb.circ <- ps.read("bc01-balb/circulars.csv") 
bc.balb.mult <- ps.read("bc01-balb/multistrand.csv")
bc.mm.circ <- ps.read("Bc_d0_1-mm/circulars.csv")  
bc.mm.mult <- ps.read("Bc_d0_1-mm/multistrand.csv")

c57.balb.circ <- ps.read("c5701-balb/circulars.csv")  
c57.balb.mult <- ps.read("c5701-balb/multistrand.csv")
c57.mm.circ <- ps.read("c57_d0_1-mm/circulars.csv")   
c57.mm.mult <- ps.read("c57_d0_1-mm/multistrand.csv") 
