### R code from vignette source 'analysis.Rnw'

###################################################
### code chunk number 1: analysis.Rnw:136-141
###################################################
require(pbh5)
require(pbls)
require(xtable)
require(ggplot2)
source("utils.R")


###################################################
### code chunk number 2: analysis.Rnw:163-165
###################################################
makeSummaryTable("Synthetic", caption = "Summary of synthetically methylated datasets used in this document.",
                 label = "tbl:synthetic")


###################################################
### code chunk number 3: analysis.Rnw:175-177
###################################################
makeSummaryTable("Lambda", caption = "Summary of lambda datasets used in this document.",
                 label = "tbl:lambda")


###################################################
### code chunk number 4: analysis.Rnw:209-211
###################################################
cmpH5 <- PacBioCmpH5("../Data/Lambda/6mA_dam+_native/data/aligned_reads.cmp.h5") 
cmpH5


###################################################
### code chunk number 5: analysis.Rnw:217-220
###################################################
group <- "ref000001/m110818_122604_42141_c100129202555500000315043109121114_s1_p0"
g <- getH5Group(cmpH5, group)
ls(g)


###################################################
### code chunk number 6: analysis.Rnw:238-239
###################################################
head(alnIndex(cmpH5), 2)


###################################################
### code chunk number 7: analysis.Rnw:244-246
###################################################
alns <- getAlignments(cmpH5, idx = c(1, 200, 3))
lapply(alns, head, n = 2)


###################################################
### code chunk number 8: analysis.Rnw:286-287
###################################################
plotIPDForReads(cmpH5, idx = 1:100, range = c(40, 50), maxReads = 5, matches = FALSE)


###################################################
### code chunk number 9: analysis.Rnw:298-299
###################################################
plotIPDForReads(cmpH5, idx = 1:100, range = c(40, 50), maxReads = 5, matches = TRUE)


###################################################
### code chunk number 10: LoadInCmpH5s
###################################################
cmpH5s <- lapply(Sys.glob("../Data/Lambda/*/data/aligned_reads.cmp.h5"), PacBioCmpH5)
names(cmpH5s) <- sapply(cmpH5s, function(h) basename(dirname(dirname(h@fileName))))


###################################################
### code chunk number 11: analysis.Rnw:318-322
###################################################
plotDensity(lapply(cmpH5s, function(cmpH5) {
  v <- getIPD(cmpH5, sample(1:nrow(cmpH5), size = 100))
  do.call(c, v)
}), legend = T, log = 'x', xlab = "IPD")


###################################################
### code chunk number 12: analysis.Rnw:326-330
###################################################
plotDensity(lapply(cmpH5s, function(cmpH5) {
  v <- getPulseWidth(cmpH5, sample(1:nrow(cmpH5), size = 100))
  do.call(c, v)
}), legend = T,  log = 'x', xlab = "Pulse Width")


###################################################
### code chunk number 13: analysis.Rnw:344-353
###################################################
l <- lapply(cmpH5s, function(cmpH5) {
  w <- sample(1:nrow(cmpH5), size = 1000)
  v <- do.call(c, getIPD(cmpH5, w))
  a <- do.call(c, lapply(getAlignments(cmpH5, w), function(b) b[,1]))
  split(v, factor(a, c("A","C","G","T")))
})
show(ggplot(melt(l), aes(x = L2, y = value, color = L1)) + geom_boxplot(outlier.shape = NA) + 
     scale_y_continuous(limits = c(0, .4)) + xlab("Base") + ylab("IPD") + 
     opts(title="IPD by Base"))


###################################################
### code chunk number 14: analysis.Rnw:357-366
###################################################
l <- lapply(cmpH5s, function(cmpH5) {
  w <- sample(1:nrow(cmpH5), size = 1000)
  v <- do.call(c, getPulseWidth(cmpH5, w))
  a <- do.call(c, lapply(getAlignments(cmpH5, w), function(b) b[,1]))
  split(v, factor(a, c("A","C","G","T")))
})
show(ggplot(melt(l), aes(x = L2, y = value, color = L1)) + geom_boxplot(outlier.shape = NA) + 
     scale_y_continuous(limits = c(0, .5)) + xlab("Base") + ylab("PulseWidth") + 
     opts(title="PulseWidth by Base"))


###################################################
### code chunk number 15: analysis.Rnw:383-384
###################################################
head(getByTemplatePosition(cmpH5, idx = 1:2, f = getIPD))


###################################################
### code chunk number 16: analysis.Rnw:391-392
###################################################
head(makeContextDataTable(cmpH5, idx = 1:2, up = 2, down = 2))


###################################################
### code chunk number 17: analysis.Rnw:395-397
###################################################
s <- summarizeByContext(cmpH5, idx = 1:100, up = 1, down = 1, statF = getPulseWidth)
head(s)


###################################################
### code chunk number 18: analysis.Rnw:406-417
###################################################
getByPositionAndStrand <- function(f = getIPD, s = 20000, e = 20025) {
  pbutils::collapse(lapply(cmpH5s, function(cmpH5) {
    x <- getByTemplatePosition(cmpH5, idx = getReadsInRange(cmpH5, 1, s, e), f = f)
    x <- subset(x, position >= s & position <= e)
    ddply(x, c("strand", "position"), function(a) {
      median(a$elt, na.rm = T)
    })
  }))
}             
byPositionAndStrandIPD <- getByPositionAndStrand()
byPositionAndStrandPW <- getByPositionAndStrand(f=getPulseWidth)


###################################################
### code chunk number 19: analysis.Rnw:424-427
###################################################
show(ggplot(byPositionAndStrandIPD, aes(x = position, y = V1, color = factor(strand), lty = L1, 
                                        group = factor(strand):factor(L1))) +
     geom_line() + ylab("Median IPD"))


###################################################
### code chunk number 20: analysis.Rnw:431-435
###################################################
show(ggplot(byPositionAndStrandPW, 
            aes(x = position, y = V1, color = factor(strand), lty = L1, 
                group = factor(strand):factor(L1))) +
     geom_line() + ylab("Median Pulse Width"))


###################################################
### code chunk number 21: analysis.Rnw:465-468
###################################################
getTemplateStrand(cmpH5)[1:10]
tmp <- getByTemplatePosition(cmpH5, idx = 1:2)
head(tmp[order(tmp$position, tmp$strand),])


###################################################
### code chunk number 22: analysis.Rnw:475-478
###################################################
tmp <- associateWithContext(cmpH5, idx = 1:2, f = getTemplatePosition, collapse = T, 
                            useReference = T)
head(tmp[order(tmp$elt),])


###################################################
### code chunk number 23: analysis.Rnw:489-495
###################################################
contextTable <- associateWithContext(cmpH5, idx = sample(1:nrow(cmpH5), size = 1000), f = getIPD, 
                                     collapse = T, useReference = T, up = 1, down = 1)
par(cex.axis = .65)
boxplot(split(contextTable$elt, contextTable$context), ylim = c(0, .5), las = 2, 
        main = "Context-specific IPD distributions", ylab = "IPD", outline = FALSE, 
        col = rep(1:4, each = 4))


###################################################
### code chunk number 24: analysis.Rnw:507-517
###################################################
par(mfrow=c(2,1), mar = c(3, 5, 1, 1))
lapply(cmpH5s[c("6mA_dam+_native", "6mA_dam-_native")], function(cmp) {
  tmp <- associateWithContext(cmp, idx = sample(1:nrow(cmp), size = 5000), 
                              f = getIPD, collapse = T, useReference = T, up = 2, down = 2)
  contextMedians <- tapply(tmp$elt, tmp$context, median, na.rm = T)
  plot(x <- 1:length(contextMedians), y <- contextMedians, pch = 16, 
       ylim = c(0, 1), xlab = "Context", xaxt = 'n', ylab = "Median")
  w <- grep("^GATC", names(y))
  text(x[w], y[w], names(y)[w], col = "darkblue", cex = 1.3)
})


###################################################
### code chunk number 25: analysis.Rnw:564-588
###################################################
cmpH5s <- lapply(Sys.glob("../Data/Synthetic/*/data/aligned_reads.cmp.h5"), PacBioCmpH5)
names(cmpH5s) <- sapply(cmpH5s, function(h) basename(dirname(dirname(h@fileName))))
modifications <- list("2x_5mC"  = c(55,74),
                      "2x_5hmC" = c(51,74),
                      "2x_4mC"  = c(55,74),
                      "2x_6mA"  = c(57,68,112,123))
syntheticStart <- 20
syntheticEnd <- 150

distributionPlot <- function(positions, nm) {
  getIPDForPosition <- function(p) {
    lapply(cmpH5s[c(nm, "control")], function(cmpH5) {
      subset(getByTemplatePosition(cmpH5, idx = sample(1:nrow(cmpH5), size = 5000)), 
             position == p & strand == 0 & read == ref)$elt
    })
  }
  par(mfrow=c(length(positions), 2), mar = c(5, 5, 4, 1))
  lapply(positions, function(z) {
    title <- paste("Position:", z)
    lst <- lapply(getIPDForPosition(z), function(k) log10(k+1/76))
    plotDensity(lst, legend = TRUE, main = title, xlab = "Started log10 of IPD")
    qqPairs(lst, main = title)
  })
}


###################################################
### code chunk number 26: analysis.Rnw:592-594
###################################################
matplot(sapply(cmpH5s, getCoverageInRange, refSeq = 1), type = 'l', ylab = "coverage", xlab = "position")
legend("topleft", names(cmpH5s), fill = 1:length(cmpH5s), bg = 'white')


###################################################
### code chunk number 27: analysis.Rnw:604-605
###################################################
distributionPlot(c(modifications[["2x_6mA"]][1], 40), "2x_6mA")


###################################################
### code chunk number 28: analysis.Rnw:616-617
###################################################
distributionPlot(c(modifications[["2x_5mC"]][1], 40), "2x_5mC")


###################################################
### code chunk number 29: analysis.Rnw:640-659
###################################################
getIPDRatios <- function() {
  getIPDMeanByPosition <- function(cmpH5) {
    s <- subset(getByTemplatePosition(cmpH5, idx = sample(1:nrow(cmpH5), size = 5000)), 
                read == ref)
    tapply(s$elt, factor(s$position, 1:getRefLength(cmpH5, 1)), mean, na.rm = T)
  }
  ctrl <- getIPDMeanByPosition(cmpH5s[["control"]])
  lapply(cmpH5s[names(cmpH5s) != "control"], function(cmpH5) {
    getIPDMeanByPosition(cmpH5)/ctrl
  })
}
ipdRatios <- getIPDRatios()
par(mfrow=c(4,1), mar = c(4, 4, 3, 1))
lapply(names(ipdRatios), function(a) {
  plot(ipdRatios[[a]], pch = 16, xlab = "Position",  ylim = c(0, 8), ylab = "IPD Ratio", main = a, 
       xlim = c(syntheticStart, syntheticEnd))
  abline(v = modifications[[a]], col = "grey", lwd = 4)
  points(ipdRatios[[a]], pch = 16)
})


###################################################
### code chunk number 30: TestPositions
###################################################
COVGS <- c(c10 = 10, c50 = 50, c100 = 100, c500 = 500, c1000 = 1000)

testPositions <- function(treatmentH5, controlH5, testStatistic = wilcox.test, 
                          targetStrands = c(0, 1), targetPositions = NA, 
                          targetCoverage = COVGS, throwOutFirstSubread = TRUE,
                          getData = getByTemplatePosition, whReference = 1) {
  if (any(is.na(targetPositions))) {
    start <- 1
    end <- refInfo(controlH5)$Length[whReference]
  } else {
    start <- min(targetPositions)
    end <- max(targetPositions)
  }
  ## find which reads to grab and then possibly downsample. 
  tReads <- getReadsInRange(treatmentH5, whReference, start, end)
  cReads <- getReadsInRange(controlH5, whReference, start, end)

  if (throwOutFirstSubread) {
    rlowest <- function(h5, r) {
      do.call(c, lapply(split(r, getMoleculeIndex(h5)[r]), function(i) {
        i[-which.min(h5$rStart[i])]
      }))
    }
    tReads <- rlowest(treatmentH5, tReads)
    cReads <- rlowest(controlH5, cReads)
  }
  
  ## retrieve the data as a data.frame.
  tData <- getData(treatmentH5, idx = tReads)
  cData <- getData(controlH5, idx = cReads)

  ## filter the data.
  tData <- subset(tData, (strand %in% targetStrands) & (position >= start & position <= end) & read == ref)
  cData <- subset(cData, (strand %in% targetStrands) & (position >= start & position <= end) & read == ref)
    
  g <- function(v, n) {
    if (length(v) > n) sample(v, size = n) else v
  }
  mapply(function(tIdxs, cIdxs) {
    lapply(targetCoverage, function(n) {
      testStatistic(tData$elt[g(tIdxs, n)], cData$elt[g(cIdxs, n)])
    })
  }, split(1:nrow(tData), factor(tData$position, start:end)), 
         split(1:nrow(cData), factor(cData$position, start:end)), SIMPLIFY = FALSE)
}

plotResult <- function(nm, g = function(z) -log10(z$p.value), ...) {
  tp <- testResults[[nm]]
  modifiedPositions <- modifications[[nm]]
  par(mfrow=c(5,1), mar = c(2, 2, 1, 1))
  sapply(1:5, function(i) {
    positions <- as.integer(names(tp))
    y <- sapply(tp, function(x) g(x[[i]]))
    plot(positions, y, xlab = "", ylab = "", 
         main = paste("-log10 p-values for coverage:", names(tp[[1]])[i], "(", nm, ")"),
         pch = 16, cex = 1.25)
    abline(v = modifiedPositions, col = 'grey', lty = 3)
  })
}
testResults <- lapply(cmpH5s[1:4], function(tH5) {
  testPositions(tH5, cmpH5s$control, targetStrands = 0)
})


###################################################
### code chunk number 31: analysis.Rnw:749-750
###################################################
plotResult("2x_6mA")


###################################################
### code chunk number 32: analysis.Rnw:758-759
###################################################
plotResult("2x_4mC")


###################################################
### code chunk number 33: analysis.Rnw:767-768
###################################################
plotResult("2x_5hmC")


###################################################
### code chunk number 34: analysis.Rnw:776-777
###################################################
plotResult("2x_5mC")


###################################################
### code chunk number 35: DoDifferentTestProcedures
###################################################
doAcrossProcedures <- function(cmpH5) {
  trimmedSlog <- function(x, trim = .975, alpha = 1/100) {
    log(x[x<quantile(x, trim)] + alpha)
  }
  testFunctions <- list(wilcox.test = wilcox.test,
                        trimmed.slog.t = function(x, y) {
                          t.test(trimmedSlog(x), trimmedSlog(y))
                        }, 
                        lr.test = function(x, y) {
                          z <- c(lx <- trimmedSlog(x), ly <- trimmedSlog(y))
                          m1 <- sum(dnorm(z, mean(z), sd(z), log = T))
                          m2 <- sum(dnorm(lx, mean(lx), sd(lx), log = T)) +
                            sum(dnorm(ly, mean(ly), sd(ly), log = T))
                          stat <- -2*(m1-m2)
                          list(statistic = stat, p.value = 1-pchisq(stat, 2))
                        })
  lapply(testFunctions, function(f) {
    testPositions(cmpH5, cmpH5s$control, targetStrands = 0, testStatistic = f)
  })
}
byTestFunction <- lapply(cmpH5s[1:4], doAcrossProcedures)

plotROC <- function(nm) {
  byT <- byTestFunction[[nm]]
  truePositives <- rep(FALSE, length(byT[[1]]))
  truePositives[modifications[[nm]]] <- TRUE
  truePositives <- factor(truePositives, c(TRUE, FALSE))
  v <- c("c10"=1, "c50"=2, "c100"=3, "c500"=4, "c1000"=5)
  pbutils::collapse(lapply(v, function(i) {
    pbutils::collapse(lapply(byT, function(testRes) {
      x <- sapply(testRes, function(r) r[[i]]$p.value)
      x <- do.call(rbind, lapply(c(sort(x), Inf), function(q) {
        tbl <- table(truth = truePositives, observed = factor(x < q, c(TRUE, FALSE)))
        c(tbl[2,1]/sum(tbl[2,]), tbl[1,1]/sum(tbl[1,]))
      }))
      colnames(x) <- c("FPR", "TPR")
      return(x)
    }))
  }))
}


###################################################
### code chunk number 36: analysis.Rnw:838-845
###################################################
tmp <- lapply(names(byTestFunction), plotROC)
names(tmp) <- names(byTestFunction)
tmp <- pbutils::collapse(tmp)
tmp$L2 <- factor(as.character(tmp$L2), c("c10", "c50", "c100", "c500", "c1000"))
tmp$L3 <- factor(as.character(tmp$L3), names(modifications))
show(ggplot(tmp, aes(x = as.numeric(FPR), y = as.numeric(TPR), color = L1)) + facet_grid(L3 ~ L2) + 
     geom_line(lwd = 1))


###################################################
### code chunk number 37: analysis.Rnw:865-867
###################################################
cmpH5s <- lapply(Sys.glob("../Data/Lambda/*/data/aligned_reads.cmp.h5"), PacBioCmpH5)
names(cmpH5s) <- sapply(cmpH5s, function(h) basename(dirname(dirname(h@fileName))))


###################################################
### code chunk number 38: analysis.Rnw:879-890
###################################################
if (! require(Biostrings)) {
  stop("Unable to execute Lambda testing examples without Biostrings package.")
}
lambda <- read.DNAStringSet("../ReferenceRepository/lambdaNEB/lambdaNEB.fa")[[1]]
matches <- matchPattern("GATC", lambda)
gatcExample <- pbutils::collapse(lapply(cmpH5s[c("6mA_dam+_native", "6mA_dam-_native")], function(cmp) {
  s <- start(matches)[1]
  e <- end(matches)[1]
  subset(getByTemplatePosition(cmp, idx = getReadsInRange(cmp, 1, s, e), f = getIPD),
         position >= s & position <= e & read == ref)
}))


###################################################
### code chunk number 39: analysis.Rnw:895-897
###################################################
show(ggplot(gatcExample, aes(x = factor(position), y = elt, color = L1)) + geom_boxplot() + 
     facet_wrap(~ strand) + scale_y_log10())


###################################################
### code chunk number 40: analysis.Rnw:918-921
###################################################
topTable <- makeTopTable(cmpH5s[["6mA_dam+_native"]], cmpH5s[["6mA_dam-_native"]], 
                         start = 1, end = 5000)
head(topTable)


###################################################
### code chunk number 41: analysis.Rnw:929-932
###################################################
plot(topTable$position, -log10(topTable$fdr), pch = 16, ylab = "-log10(FDR adjusted p.values)", 
     xlab = "position", main = "-log10 FDR adjusted p.values")
abline(v = start(matches) + 1)


###################################################
### code chunk number 42: analysis.Rnw:952-953
###################################################
sessionInfo()


