harvard.dir <- '//mdabam1/bioinfo/Public/prostate-meta/Harvard'

harvard.data <- read.table(file.path(harvard.dir,
                                     'HarvardUpdatedDataAndAnnotation.txt'),
                           header=TRUE, sep='\t', row.names=2)
harvard.genes <- harvard.data[, 2:3]
harvard.data  <- harvard.data[, 4:ncol(harvard.data)]


harvard.clin <- read.table(file.path(harvard.dir, 'Harvard-clinical-new.txt'),
                           header=TRUE, sep='\t', row.names=NULL,
                           comment='', quote='')
dimnames(harvard.clin)[[1]] <- as.character(harvard.clin$Sample.ID)

tumor <- dimnames(harvard.data)[[2]][51:102]
harvard.info <- data.frame(sample=factor(dimnames(harvard.data)[[2]]),
                           status=factor(rep(c('normal', 'cancer'), times=c(50,52))),
                           gleason=factor(paste('G',
                             c(rep(0, 50),
                             as.character(harvard.clin[tumor, 'GS'])),
                             sep='')))
rm(tumor)

plot(t(harvard.data[1,]))

harvard.tr2 <- TailRankTest(logb(harvard.data, 2),
                              harvard.info$status=='normal',
                              direction='two')

summary(harvard.tr2)

x <- gene.info$Cluster.ID[as.logical(tr2)]
y1 <- as.character(gene.info$Cluster.ID)
z1 <- is.element(y1,x) & y1 != ''

y <- as.character(harvard.genes[,1])
z <- is.element(y, x) & y != ''

w <- harvard.genes[as.logical(harvard.tr2), 1]
z2 <- is.element(y1,w) & y1 != ''
