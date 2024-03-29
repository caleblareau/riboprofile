A <- read.table ("././test.data/read.dist.sample.29.txt", sep="\t")
loc1 <- c(-30:50)
loc2 <- c(-50:30)
B1 <- apply(A[1,2:ncol(A)], 2, sum)
B2 <- apply(A[2,2:ncol(A)], 2, sum)
pdf (file="././test.data/plot.readDist.29.pdf")
plot(loc1, B1, xlab="Around start codon", ylab="RPM", type="h", lwd=2, ylim=c(0,max(c(B1,B2))))
plot(loc2, B2, xlab="Around stop codon", ylab="RPM", type="h", lwd=2, ylim=c(0,max(c(B1,B2))))
dev.off()
