# Copyright (C) Kevin R. Coombes, 2007-2012

# tail-rank-test.R
# Copyright: Kevin R. Coombes, 2004.

require('methods')

# Implements the tail-rank test described in our paper.

############################
# We start with a class to hold and display the results.

setClass("TailRankTest",
         slots = c(statistic='numeric',	# nonnegative integer vector
                   direction='character',	# up, down, or two-sided
                   N1='numeric',		# number of healthy
                   N2='numeric',		# number of cancer
                   tau='numeric or NULL',  # upper thresholds
                   rho='numeric or NULL',  # lower thresholds
                   specificity='numeric',	# between 0 and 1
                   tolerance='numeric',	# between 0 and 1
                   confidence='numeric',	# between 0 and 1
                   model='character',      # 'bb' or 'binomial'
                   cutoff='numeric'	# integer: max value under null
                   )
         )

setMethod('summary', signature(object='TailRankTest'),
          function(object, ...) {
  cat(paste('A tail-rank test object in the',
            object@direction,'direction.\n'))
  cat(paste('The test was performed using the', object@model, 'model.\n'))
  cat(paste('Specificity:', object@specificity,
            'computed with tolerance', object@tolerance, '\n'))
  cat(paste('Significance cutoff:', object@cutoff,
            'based on a family-wise error rate less than',
            1-object@confidence, '\n'))
  cat(paste('There are', sum(as.logical(object)),
            'tail-rank statistics that exceed the cutoff\n'))
})

# altered by KRC to use beta-binomial on 5 July 2007
setMethod('hist', signature(x='TailRankTest'),
          function(x,
                   overlay=FALSE,
                   xlab='tail-rank statistic',
                   main='',
                   ...) {
            if(overlay) {
              xx <- 0:length(x@statistic)
              qq <- 1 - pnorm(toleranceBound(x@specificity, x@tolerance, x@N1))
              if (x@model=="binomial") {
                dd <- dbinom(xx, x@N2, qq)
              } else {
                W <- x@N1+2
                psi <- x@specificity
                dd <- dbb(xx, x@N2, W*(1-psi), W*psi)
              }
              hist(x@statistic, xlab=xlab, main=main,
                   prob=TRUE, ...)
              lines(xx, dd, col='red')
            } else {
              hist(x@statistic, xlab=xlab, main=main, ...)
            }
            # also add a tick mark at the significance cutoff
            rug(x@cutoff)
            invisible(x)
          })

setMethod('as.logical', signature(x='TailRankTest'),
          function(x, ...) {
  x@statistic > x@cutoff
})

if (!isGeneric('getStatistic'))
  setGeneric('getStatistic', function(object, ...) {
    standardGeneric('getStatistic')})

setMethod('getStatistic', signature(object='TailRankTest'),
          function(object, ...) {
              object@statistic
          })

############################
# Here is the constructor that actually performs the test
#
# We compute the tail rank statistic for each gene (as rows of
# the data matrix). The data is split into two groups; the
# first group is used to estimate a tolerance bound (defaults
# to 90%) on a specific quantile (defaults to 95%) of the
# distribution of each gene. The tail-rank statistic is the
# number of samples in the second group outside the bound.
# The test can be applied in the "up", "down", or "two-sided"
# direction, depending on the kinds of markers being sought.
# Also computes the cutoff for significance based on a confidence
# level that is "1 - FWER" for a desired family-wise error rate.
TailRankTest <-
  function(data,
           classes,
           specificity=0.95,
           tolerance=0.50,
           model=c('bb', 'betabinomial', 'binomial'),
           confidence=0.95,
           direction='up') {

  # check consistency of the arguments
    model <- match.arg(model)
    # start with the direction of the test
    directions <- c('up', 'down', 'two-sided')
    dir <- pmatch(direction, directions)
    if(is.na(dir))
      stop(paste(direction, 'is not a valid direction for comparison'))
    if(dir == -1)
      stop(paste(direction, 'is ambiguous'))

    # handle the special inputr case of an ExpressionSet
    if(inherits(data, 'ExpressionSet')) {
      if(is.character(classes)) {
        classes <- as.factor(pData(data)[,classes])
      }
      data <- exprs(data)
    }
    # now convert classes to a logical vector
    if(is.factor(classes))
      classes <- classes == levels(classes)[[1]]
    classes <- as.logical(classes)
    if(sum(classes)==0 | sum(!classes)==0)
      stop('one of the groups is empty!')

    # the data should be a matrix, with number of columns = length of classes
    data <- as.matrix(data)
    if(ncol(data) != length(classes))
      stop(paste('The number of data columns must equal the',
                 'length of the splitting vector'))

    # make sure specificity, tolerance, and confidence are in the correct range
    if (specificity > 1 | specificity < 0)
      stop('specificity must be between 0 and 1')
    if (tolerance > 1 | tolerance < 0)
      stop('tolerance must be between 0 and 1')
    if (confidence > 1 | confidence < 0)
      stop('confidence must be between 0 and 1')

    # we're finally ready to start
    healthyN <- sum(classes)
    cancerN  <- sum(!classes)
    healthyMean <- matrixMean(data[,classes])
    cancerMean  <- matrixMean(data[,!classes])
    healthySD <- sqrt(matrixVar(data[,classes], healthyMean))
    cancerSD  <- sqrt(matrixVar(data[,!classes], cancerMean))
    toleranceFactor <- toleranceBound(specificity, tolerance, healthyN)

    # Now we actually compute the tail rank statistics
    tumors <- as.data.frame(data[,!classes]) # i really hate this
    # you need to use a data.frame in order to get the comparison
    # operators to work correctly
    tau <- NULL
    rho <-  NULL
    if(dir == 1) { # up
      tau <- healthyMean + healthySD*toleranceFactor
      trs <- as.matrix(tumors > tau) %*% rep(1, ncol(tumors))
      cutoff <- tailRankCutoff(nrow(data), healthyN, cancerN, specificity,
                               confidence, model)
    } else if (dir == 2) { # down
      rho <- healthyMean - healthySD*toleranceFactor
      trs <- as.matrix(tumors < rho) %*% rep(1, ncol(tumors))
      cutoff <- tailRankCutoff(nrow(data), healthyN, cancerN, specificity,
                               confidence, model)
    } else { # two-sided
      tau <- healthyMean + healthySD*toleranceFactor
      ups <- as.matrix(tumors > tau) %*% rep(1, ncol(tumors))
      rho <- healthyMean - healthySD*toleranceFactor
      downs <- as.matrix(tumors < rho) %*% rep(1, ncol(tumors))
      trs <- ((ups + downs) + abs(ups-downs))/2	# maximum of one-sided
      cutoff <- tailRankCutoff(2*nrow(data), healthyN, cancerN, specificity,
                               confidence, model)
    }

    trs <- as.vector(trs)
    # let's build a reasonable object for the return journey
    new("TailRankTest",
        statistic=trs,
        direction=directions[dir],
        N1=healthyN,
        N2=cancerN,
        tau=as.vector(tau), rho=as.vector(rho),
        specificity=specificity,
        tolerance=tolerance,
        confidence=confidence,
        model=model,
        cutoff=cutoff
        )
}

