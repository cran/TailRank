
# TR-prostate.R
# Copyright: Kevin R. Coombes, 2004.

# This file contains the application of the tail-rank test to
# to a prostate cancer data set. For the data, see the paper by
# Lapointe et al., PNAS 2004; 101:811-816. For a description of
# the test, see our technical report.

############################
# load the processed prostate cancer data from Lapointe's study
# We reprocessed the raw data ourselves as described in our paper,
# and saved the processed data to focus on how the tail-rank test
# works.
f <- 'lapointe.Rda'
if (file.exists(f)) {
  load(f)
} else {
  prostate.dir <- 'f:/imported/StanfordProstate'
  expression.data <- read.table(file.path(prostate.dir, 'normalizedData.txt'),
      header=TRUE, sep='\t', row.names=1)
  clinical.info <- read.table(file.path(prostate.dir, 'clinicalInfo.txt'),
                              header=TRUE, sep='\t', row.names=1)
  gene.info <- read.table(file.path(prostate.dir, 'geneInfo.txt'),
                          header=TRUE, sep='\t', row.names=1,
                        quote='', as.is=TRUE)
  save(expression.data, clinical.info, gene.info, file="lapointe.Rda")
}

require(TailRank)

############################
# prostate data gene-wise mean, variance, and t-statistics
levels(clinical.info$Status)

healthy.n <- sum(clinical.info$Status=='N')
cancer.n  <- sum(clinical.info$Status!='N')

healthy.mean <- matrixMean(expression.data[,clinical.info$Status=='N'])
cancer.mean  <- matrixMean(expression.data[,clinical.info$Status!='N'])

healthy.var <- matrixVar(expression.data[,clinical.info$Status=='N'],
                          healthy.mean)
cancer.var  <- matrixVar(expression.data[,clinical.info$Status!='N'],
                          cancer.mean)

pooled.sd <- sqrt( (healthy.var * (healthy.n - 1) + cancer.var * (cancer.n - 1))/(healthy.n + cancer.n - 2) ) 

t.statistics <- (cancer.mean - healthy.mean)/pooled.sd/sqrt(1/healthy.n + 1/cancer.n)

p.values <- sapply(t.statistics, function(tv, df) {
  2*(1-pt(abs(tv), df))
}, healthy.n + cancer.n - 2)
p.values[p.values==0] <- 1.1e-16

prostate.bum <- Bum(p.values)
image(prostate.bum)

############################
# temporary function to compute the tail-rank statistics
krc.test <- function(tumors, threshold, direction='up') {
  # this only works if 'tumors' is a data.frame, not a matrix
  # that's because the comparison operators behave differently, sigh.
  if (direction == 'up') {
    val <- as.matrix(tumors > threshold) %*% rep(1, ncol(tumors))
  } else {
    val <- as.matrix(tumors < threshold) %*% rep(1, ncol(tumors))
  }
  val
}

############################
# 90% upper tolerance bounds on 95th percentile of healthy

confidence <- 0.95
tolerance <- 0.90
specificity <- 0.95
tolerance.factor <- toleranceBound(specificity, tolerance, healthy.n)
tau <- healthy.mean + sqrt(healthy.var)*tolerance.factor
# also get the lower bound for 5th percentile
rho <- healthy.mean - sqrt(healthy.var)*tolerance.factor

############################
# the actual tail-rank statistics
tumor.samples <- expression.data[, clinical.info$Status!='N']
up.markers <- krc.test(tumor.samples, tau)
down.markers <- krc.test(tumor.samples, rho, 'down')

cutoff <- tailRankCutoff(nrow(tumor.samples),
                           N1=sum(clinical.info$Status=='N'),
                           N2=ncol(tumor.samples),
                           psi=specificity,
                           conf=confidence)

sum(up.markers > cutoff)   # 1766
sum(down.markers > cutoff) # 1930
max.markers <- ((up.markers + down.markers) + abs(up.markers-down.markers))/2
sum(max.markers > cutoff)

############################
# recomputing without tolerance bounds

tolerance.factor <- toleranceBound(specificity, 0.50, healthy.n)
tau1 <- healthy.mean + sqrt(healthy.var)*tolerance.factor
rho1 <- healthy.mean - sqrt(healthy.var)*tolerance.factor
up1.markers <- krc.test(tumor.samples, tau1)
down1.markers <- krc.test(tumor.samples, rho1, 'down')

sum(up1.markers > cutoff)   # 3498 !
sum(down1.markers > cutoff) # 3407 !

############################
# recomputing with nonparametric percentiles

x <- apply(expression.data[,clinical.info$Status=='N'], 1, sort)
tau2 <- x[40,]
rho2 <- x[2,]
rm(x)
up2.markers <- krc.test(tumor.samples, tau2)
down2.markers <- krc.test(tumor.samples, rho2, 'down')

sum(up2.markers > cutoff)   # 3543 !
sum(down2.markers > cutoff) # 3389 !

############################
# The above code was used during development, After we were certain
# that it worked properly, we put together a class and some functions
# to implement it so others could reuse it. In the refined version, we
# simply do the following:

# Get the positive (up) biomarkers. This relies on the default values
# of several parameters, which might be a bad idea.
trup <- TailRankTest(expression.data, clinical.info$Status=='N')
summary(trup)

# The fully specified version of the same command is
tolerance <- 0.90
specificity <- 0.95
confidence <- 0.95
trup <- TailRankTest(expression.data, clinical.info$Status=='N',
                       specificity=specificity, tolerance=tolerance,
                       confidence=confidence, direction='up')

# Now we can get the results
summary(trup)
hist(trup, breaks=70, overlay=TRUE)

# Get the negative (down) biomarkers.
trdn <- TailRankTest(expression.data, clinical.info$Status=='N',
                       direction='down')
summary(trdn)
hist(trdn, breaks=70, prob=TRUE)

# Get the two-sided biomarkers.
tr2 <- TailRankTest(expression.data, clinical.info$Status=='N',
                        direction='two')
summary(tr2)
hist(tr2, breaks=70, prob=TRUE)

#forJim <- TailRankTest(expression.data^2, clinical.info$Status=='N',
#                         direction='two')
#forJim.up <- TailRankTest(expression.data^2, clinical.info$Status=='N',
#                         direction='up')
#forJim.dn <- TailRankTest(expression.data^2, clinical.info$Status=='N',
#                         direction='down')

############################
# simulations

# one to start
nr <- 40000
nc <- 110
fake.data <- matrix(rnorm(nr*nc), ncol=nc)
fake.class <- rep(c(TRUE, FALSE), c(40, 70))
null.tr2 <- TailRankTest(fake.data, fake.class)
summary(null.tr2)
hist(null.tr2, overlay=T)

# cleanup
rm(fake.data, fake.class, nr, nc, null.tr2)

# Now we summarized 100 simulations by counting how many false
# positives we get. Sizes are the same as in the prostate cancer
# data set. This takes a long time, so be prepared to work
# on something else while it runs in the background.
notol.results <- sapply(1:100, function(x) {
  print(x)
  nr <- 42129
  nc <- 112
  fake.data <- matrix(rnorm(nr*nc), ncol=nc)
  fake.class <- rep(c(TRUE, FALSE), c(41, 71))
  null.tr <- TailRankTest(fake.data, fake.class, tolerance=0.5) # no tolerance
  selector <- as.logical(null.tr)
  sum(selector)
})

summary(notol.results)
hist(notol.results, breaks=seq(-0.5, max(notol.results)+0.5, by=1))
hist(notol.results, breaks=seq(-0.5, max(notol.results)+0.5, by=1), plot=F)$counts

sim.results <- sapply(1:100, function(x) {
  print(x)
  nr <- 42129
  nc <- 112
  fake.data <- matrix(rnorm(nr*nc), ncol=nc)
  fake.class <- rep(c(TRUE, FALSE), c(41, 71))
  null.tr <- TailRankTest(fake.data, fake.class, tolerance=0.95)
  selector <- as.logical(null.tr)
  sum(selector)
})

summary(sim.results)
hist(sim.results, breaks=seq(-0.5, max(sim.results)+0.5, by=1))
hist(sim.results, breaks=seq(-0.5, max(sim.results)+0.5, by=1), plot=F)$counts

# In spite of our attempts to be hyper-conservative, we still get a
# small number of false positives in almost every case. The median number
# of false positives is 1 out of 42000 and the maximum number is 4 out of
# 42000, which I think most biologists can live with comfortably.

# When we omit the tolerance bounds and use the theoretical maximum
# likelihood estimate, the number of false positives ranges from 47
# to 79 with a median of 63. This is probably less likely to be acceptable.

############################
# resimulating with nonparametric percentiles

nonpar.results <- sapply(1:100, function(x) {
  cat('...')
  nr <- 42129
  nc <- 112
  fake.data <- matrix(0, ncol=nc, nrow=nr)
  cat('simulating...')
  fake.data <- matrix(rnorm(nr*nc), ncol=nc)
  fake.class <- rep(c(TRUE, FALSE), c(41, 71))
  cat('sorting...')
  y <- apply(fake.data[,1:41], 1, sort)
  cat('selecting...')
  tau2 <- y[40,]
  cat('testing...')
  up2.markers <- krc.test(fake.data[,42:112], tau2)
  cat('counting')
  n <- sum(up2.markers > 15)
  cat(paste(' ', n, '\n'))
  n
})

summary(nonpar.results)
hist(nonpar.results, breaks=seq(-0.5, max(nonpar.results)+0.5, by=1))
hist(nonpar.results, breaks=seq(-0.5, max(nonpar.results)+0.5, by=1), plot=F)$counts

############################
# permutations

# one to start
labeled <- sample(clinical.info$Status=='N', nrow(clinical.info))
sum(labeled & clinical.info$Status=='N')
permed <- TailRankTest(expression.data, labeled)
sum(as.logical(permed)) # only 5 false positives on the first go-round!

# now we do 100 permutations of the labels
perm.results <- sapply(1:100, function(x) {
  print(x)
  labeled <- sample(clinical.info$Status=='N', nrow(clinical.info))
  permed <- TailRankTest(expression.data, labeled)
  list(one=sum(as.logical(permed)),
       two=sum(labeled & clinical.info$Status=='N'),
       three=sum(labeled & clinical.info$Status=='L'))
})
falsepos <- unlist(perm.results[1,])
overlap  <- unlist(perm.results[2,])

# first, we check that the permutations are reasonably balanced
summary(overlap)
hist(overlap, breaks=seq(8.5, 21.5, by=1), prob=TRUE)
lines(8:22, dbinom(8:22, 41, 41/112), col='blue')

# now we look at the distribution of false positives
summary(falsepos)
hist(falsepos, breaks=seq(-0.5, max(falsepos)+0.5, by=1))
hp <- hist(falsepos,
           breaks=seq(-0.5, max(falsepos)+0.5, by=1), plot=F)$counts
hp

# There is a single "large" outlier (with 22 false positives) among the
# 100 permutations. Of course, this number is still really small compared
# to the 1500+ markers found in the unpermuted data.

# The outlier is not associated with a particularly unbalanced set.
plot(overlap, falsepos)


# now we do 100 permutations of the labels
notol.perm.results <- sapply(1:100, function(x) {
  print(x)
  labeled <- sample(clinical.info$Status=='N', nrow(clinical.info))
  permed <- TailRankTest(expression.data, labeled, tolerance=0.5)
  list(one=sum(as.logical(permed)),
       two=sum(labeled & clinical.info$Status=='N'),
       three=sum(labeled & clinical.info$Status=='L'))
})
notol.falsepos <- unlist(notol.perm.results[1,])
notol.overlap  <- unlist(notol.perm.results[2,])
summary(notol.falsepos)
hist(notol.falsepos, breaks=seq(-0.5, max(notol.falsepos)+0.5, by=1))

# If we use the tolerance bounds to set the percentiles, we get between
# 0 and 33 false positives, with a median of only 2. This is very similar
# to the siulation results. However, if we use the MLE estimate without the
# tolerance adjustment, then we get between  and 168 false positives, with
# a median of 30, which is again in line with the simulation results.

nonpar.perm.results <- sapply(1:100, function(x) {
  cat(paste('[', x, ']...', sep=''))
  labeled <- sample(clinical.info$Status=='N', nrow(clinical.info))
  cat('sorting...')
  y <- apply(expression.data[,labeled], 1, sort)
  cat('selecting...')
  tau2 <- y[40,]
  cat('testing...')
  up2.markers <- krc.test(expression.data[,!labeled], tau2)
  cat('counting')
  n <- sum(up2.markers > 15)
  cat(paste(' ', n, '\n'))
  n
})
summary(nonpar.perm.results)
hist(nonpar.perm.results, breaks=seq(-0.5, max(nonpar.perm.results)+0.5, by=1))

# we next compare the observed permutation results with the results
# of the simulation of independent data
sp <- hist(sim.results, breaks=seq(-0.5, max(falsepos)+0.5, by=1),
           plot=F)$counts

windows(width=8, height=4)
par(mai=c(1, 1, 0.3, 0.3), cex=1.2)
barplot(rbind(hp, sp), beside=TRUE, names.arg=(1:length(hp))-1,
        col=c('black', 'gray'),
        xlab='Number of false positives',
        ylab='Number of simulations')
legend(50, 20, c('Permuted', 'Simulated'), adj=0, pch=15,
       col=c('black', 'gray'))

#################################
healthy.median <- apply(expression.data[,clinical.info$Status=='N'], 1, median)

healthy.iqr <- apply(expression.data[,clinical.info$Status=='N'], 1,
                     function(x) {quantile(x, 0.75) - quantile(x, 0.25)})

new.tau <- healthy.median + (qnorm(0.95)*healthy.iqr)/(2*qnorm(0.75))

iqr.sim.results <- sapply(1:2, function(x) {
  print(x)
  nr <- 42129
  nc <- 112
  fake.data <- matrix(rnorm(nr*nc), ncol=nc)
  fake.class <- rep(c(TRUE, FALSE), c(41, 71))
  h.median <- apply(fake.data[,fake.class], 1, median)
  h.iqr <- apply(fake.data[,!fake.class], 1,
                 function(x) {quantile(x, 0.75) - quantile(x, 0.25)})
  new.tau <- h.median + (qnorm(0.95)*h.iqr)/(2*qnorm(0.75))
  sum(new.tau > cutoff)
})

summary(sim.results)
hist(sim.results, breaks=seq(-0.5, max(sim.results)+0.5, by=1))
hist(sim.results, breaks=seq(-0.5, max(sim.results)+0.5, by=1), plot=F)$counts

#################################
t.statistics <- as.vector(t.statistics)
descr <- as.character(clinical.info$Sample)

marker.data <- expression.data[as.logical(tr2),]
dim(marker.data)

tronly.data <- expression.data[as.logical(tr2) & abs(t.statistics) < 4.25,]
dim(tronly.data)

tonly.data <- expression.data[!as.logical(tr2) & abs(t.statistics) > 4.25,]
dim(tonly.data)

dum <- as.dist((1-cor(marker.data))/2)
dumonly <- as.dist((1-cor(tronly.data))/2)
dumt <- as.dist((1-cor(tonly.data))/2)

hd <- hclust(dum, 'complete')
hdonly <- hclust(dumonly, 'complete')
hdt <- hclust(dumt, 'complete')

dev.set(2)
plclust(hd, labels=clinical.info$Status)
plclust(hdonly, labels=clinical.info$Status)
plclust(hdt, labels=clinical.info$Status)

windows(height=6, width=6, pointsize=10)

strapped <- clust.strap(marker.data, hclust.cut, K=10, method='complete',
                        verbose=TRUE, n.boot=100)
strappedonly <- clust.strap(tronly.data, hclust.cut, K=12, method='complete',
                        verbose=TRUE, n.boot=200)

strappedt <- clust.strap(tonly.data, hclust.cut, K=8, method='complete',
                        verbose=TRUE, n.boot=20)

dev.set(3)
colorList <- c('red', 'green', 'cyan',
               'purple', 'orange', 'magenta', 'green', 'blue')
heatmap(strapped, Rowv=as.dendrogram(hd), symm=TRUE, revC=FALSE,
        labRow=descr, labCol=descr,
        ColSideColors=colorList[as.numeric(clinical.info$Status)],
        RowSideColors=colorList[3+as.numeric(clinical.info$Subgroups)],
        col=blueyellow(64)
        )
heatmap(strappedonly, Rowv=as.dendrogram(hdonly), symm=TRUE, revC=FALSE,
        labRow=descr, labCol=descr,
        ColSideColors=colorList[as.numeric(clinical.info$Status)],
        RowSideColors=colorList[3+as.numeric(clinical.info$Subgroups)],
        col=blueyellow(64)
        )
heatmap(strappedt, Rowv=as.dendrogram(hdt), symm=TRUE, revC=FALSE,
        labRow=descr, labCol=descr,
        ColSideColors=colorList[as.numeric(clinical.info$Status)],
        RowSideColors=colorList[3+as.numeric(clinical.info$Subgroups)],
        col=blueyellow(64)
        )



heatmap(strapped, Rowv=as.dendrogram(hd), Colv=as.dendrogram(hdonly),
        symm=TRUE, revC=FALSE,
        labRow=descr, labCol=descr,
        ColSideColors=colorList[as.numeric(clinical.info$Status)],
        RowSideColors=colorList[as.numeric(clinical.info$Status)],
        col=blueyellow(64)
        )

#################################
pos.marker <- marker.data > tau[as.logical(tr2)]
neg.marker <- marker.data < rho[as.logical(tr2)]

marker <- pos.marker | neg.marker
marker <- marker[rev(order(apply(marker, 1, sum))),]

whatever <- apply(marker, 2, cumsum)
whenever <- apply(marker[order(apply(marker, 1, sum)),], 2, cumsum)

plot(whatever[500,]); abline(v=c(9.5, 50.5))
plot(whenever[500,]); abline(v=c(9.5, 50.5))


binary.marker <- matrix(0, nrow=nrow(marker), ncol=ncol(marker))
binary.marker[marker] <- 1

cln <- clinical.info$Status!='N'
sum(cln)

bim <- binary.marker[, cln]
dim(bim)

agreement <- (bim %*% t(bim)) + (1-bim) %*% t(1-bim)
x <- seq(5, 75, by=0.03)
b <- 0.5 + (5:75)
hist(agreement, breaks=b, prob=TRUE)
lines(x, dnorm(x, 41, 6.25), col='red')

hc <- hclust(as.dist(71-agreement), method='average')
plclust(hc)
xc <- cutree(hc, k=2)
ran <- (1:length(xc))[xc==1]

ran <- sample(nrow(agreement), 300)
hc2 <- hclust(as.dist(71-agreement[ran, ran]), method='average')
plclust(hc2, labels=as.character(ran))

plot(apply(bim, 1, sum))
plot(apply(bim, 1, sum), jitter(xc))

weights <- (apply(bim, 1, sum))/71
expected <- (outer(weights, weights)+outer(1-weights, 1-weights))*71
summary((agreement-expected)[expected<71])

x0 <- seq(-25, 40, by=0.3)
length(x0)
hist(agreement-expected, breaks=0.5 + (-25:40), prob=TRUE)
lines(x0, dnorm(x0, 3, 4.75), col='red')

ran <- sample(nrow(agreement), 300)
heatmap((agreement-expected)[ran, ran], col=blueyellow(64))

###########################################
clin.new <- read.table('//mdabam1/bioinfo/Public/prostate-meta/Stanford/clinicalinfo-stanford-new.txt',
           header=TRUE, row.names=NULL, sep='\t')

gleason

samtype <- ordered(clinical.info$Status, levels=c('N', 'T', 'L'))

fool <- function(NN) {
  marker.count <- apply(binary.marker[1:NN,], 2, sum)
  plot(marker.count, samtype, yaxt='n', ylab='',
       xlab=paste('Number of markers (out of ', NN, ')', sep=''))
  mtext(c('N', 'T', 'L'), at=1:3, side=2, las=2, line=1, font=2)
  invisible(NN)
}

opar <- par(mfrow=c(2,2), mai=c(0.8,0.6, 0.1, 0.1))
fool(nrow(binary.marker))
fool(1000)
fool(100)
fool(50)
par(opar)

plot(samtype, marker.count)

plot(gleason, marker.count)

plot(marker.count, gleason)

