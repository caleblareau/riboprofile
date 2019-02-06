library ("e1071")
load ("ribORF.rda")
A1 <- read.table ("./example.data/input.parameters.txt", sep="\t", header=T)
f1 <- A1[,9]
f2 <- A1[,10]
pme <- A1[,14]
A2 <- cbind(f1,f2,pme)
L1 <- sprintf("%.3f", predict(pred1, A2, decision.values = TRUE))
out <- cbind(A1, pvalue=L1)
write.table (out, "./example.data/pred.pvalue.parameters.txt", sep="\t", quote=F, row.names=F, col.names=T)
