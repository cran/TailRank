cm <- 1/2.32

tr.stats <- getStatistic(tr2)

windows(8.7*cm, 6*cm, pointsize=10)

par(mai=c(0.7, 0.7, 0.1, 0.1))
plot(t.statistics, tr.stats, pch='.',
     xlab='T Statistic', ylab='Tail Rank Statistic')
abline(h=15.5, lwd=2)
abline(v=c(4.25, -4.25), lwd=2)

i.t.weird <- c(22552, 23827, 23981, 24175) # deep magic

x <- ordered(clinical.info$Status, levels=c('N', 'T', 'L'))
opar <- par(mai=c(0.3, 0.4, 0.4, 0.1), mfrow=c(2,2))
for (i in i.t.weird) {
  y <- as.vector(t(expression.data[i,]))
  label <- as.character(gene.info$Gene.Symbol[i])
  if (label=='') {
    label <- as.character(gene.info$Accession[i])
  }
  stripchart(y ~ x, xlab='', main=label, method='jitter')
}
par(opar)

k.weird <- tr.stats > 15 & (abs(t.statistics)<1.25)
sum(k.weird)
i.k.weird <- (1:length(k.weird))[k.weird]

gene.info$Gene.Symbol[i.k.weird]

i.k.weird <- c(17710, 36012, 16038, 5266)
x <- ordered(clinical.info$Status, levels=c('N', 'T', 'L'))
opar <- par(mai=c(0.3, 0.4, 0.4, 0.1), mfrow=c(2,2))
for (i in i.k.weird) {
  y <- as.vector(t(expression.data[i,]))
  label <- as.character(gene.info$Gene.Symbol[i])
  if (label=='') {
    label <- as.character(gene.info$Accession[i])
  }
  stripchart(y ~ x, xlab='', main=label, method='jitter')
}
par(opar)

toogood <- tr.stats > 55
sum(toogood)
i.toogood <- (1:length(toogood))[toogood]
i <- i.toogood[22]
  y <- as.vector(t(expression.data[i,]))
  label <- as.character(gene.info$Gene.Symbol[i])
  if (label=='') {
    label <- as.character(gene.info$Accession[i])
  }
  stripchart(y ~ x, xlab='', main=label, method='jitter')
