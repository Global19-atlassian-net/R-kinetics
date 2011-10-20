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
### code chunk number 2: analysis.Rnw:155-157
###################################################
makeSummaryTable("Synthetic", caption = "Summary of synthetically methylated datasets used in this document.",
                 label = "tbl:synthetic")


###################################################
### code chunk number 3: analysis.Rnw:167-169
###################################################
makeSummaryTable("Lambda", caption = "Summary of lambda datasets used in this document.",
                 label = "tbl:lambda")


###################################################
### code chunk number 4: analysis.Rnw:201-203
###################################################
cmpH5 <- PacBioCmpH5("../Data/Lambda/6mA_dam+_native/data/aligned_reads.cmp.h5") 
cmpH5


###################################################
### code chunk number 5: analysis.Rnw:209-212
###################################################
group <- "ref000001/m110818_122604_42141_c100129202555500000315043109121114_s1_p0"
g <- getH5Group(cmpH5, group)
ls(g)


###################################################
### code chunk number 6: analysis.Rnw:230-231
###################################################
head(alnIndex(cmpH5), 2)


###################################################
### code chunk number 7: analysis.Rnw:236-238
###################################################
alns <- getAlignments(cmpH5, idx = c(1, 200, 3))
lapply(alns, head, n = 2)


###################################################
### code chunk number 8: analysis.Rnw:278-279
###################################################
plotIPDForReads(cmpH5, idx = 1:100, range = c(40, 50), maxReads = 5, matches = FALSE)


###################################################
### code chunk number 9: analysis.Rnw:290-291
###################################################
plotIPDForReads(cmpH5, idx = 1:100, range = c(40, 50), maxReads = 5, matches = TRUE)


###################################################
### code chunk number 10: LoadInCmpH5s
###################################################
cmpH5s <- lapply(Sys.glob("../Data/Lambda/*/data/aligned_reads.cmp.h5"), PacBioCmpH5)
names(cmpH5s) <- sapply(cmpH5s, function(h) basename(dirname(dirname(h@fileName))))


###################################################
### code chunk number 11: analysis.Rnw:310-314
###################################################
plotDensity(lapply(cmpH5s, function(cmpH5) {
  v <- getIPD(cmpH5, sample(1:nrow(cmpH5), size = 100))
  do.call(c, v)
}), legend = T, log = 'x', xlab = "IPD")


###################################################
### code chunk number 12: analysis.Rnw:318-322
###################################################
plotDensity(lapply(cmpH5s, function(cmpH5) {
  v <- getPulseWidth(cmpH5, sample(1:nrow(cmpH5), size = 100))
  do.call(c, v)
}), legend = T,  log = 'x', xlab = "Pulse Width")


###################################################
### code chunk number 13: analysis.Rnw:336-345
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
### code chunk number 14: analysis.Rnw:349-358
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
### code chunk number 15: analysis.Rnw:375-376
###################################################
head(getByTemplatePosition(cmpH5, idx = 1:2, f = getIPD))


###################################################
### code chunk number 16: analysis.Rnw:383-384
###################################################
head(makeContextDataTable(cmpH5, idx = 1:2, up = 2, down = 2))


###################################################
### code chunk number 17: analysis.Rnw:387-389
###################################################
s <- summarizeByContext(cmpH5, idx = 1:100, up = 1, down = 1, statF = getPulseWidth)
head(s)


###################################################
### code chunk number 18: analysis.Rnw:398-409
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
### code chunk number 19: analysis.Rnw:416-419
###################################################
show(ggplot(byPositionAndStrandIPD, aes(x = position, y = V1, color = factor(strand), lty = L1, 
                                        group = factor(strand):factor(L1))) +
     geom_line() + ylab("Median IPD"))


###################################################
### code chunk number 20: analysis.Rnw:423-427
###################################################
show(ggplot(byPositionAndStrandPW, 
            aes(x = position, y = V1, color = factor(strand), lty = L1, 
                group = factor(strand):factor(L1))) +
     geom_line() + ylab("Median Pulse Width"))


###################################################
### code chunk number 21: analysis.Rnw:457-460
###################################################
getTemplateStrand(cmpH5)[1:10]
tmp <- getByTemplatePosition(cmpH5, idx = 1:2)
head(tmp[order(tmp$position, tmp$strand),])


###################################################
### code chunk number 22: analysis.Rnw:467-470
###################################################
tmp <- associateWithContext(cmpH5, idx = 1:2, f = getTemplatePosition, collapse = T, 
                            useReference = T)
head(tmp[order(tmp$elt),])


###################################################
### code chunk number 23: analysis.Rnw:480-486
###################################################
contextTable <- associateWithContext(cmpH5, idx = sample(1:nrow(cmpH5), size = 1000), f = getIPD, 
                                     collapse = T, useReference = T, up = 1, down = 1)
par(cex.axis = .65)
boxplot(split(contextTable$elt, contextTable$context), ylim = c(0, .5), las = 2, 
        main = "Context-specific IPD distributions", ylab = "IPD", outline = FALSE, 
        col = rep(1:4, each = 4))


###################################################
### code chunk number 24: analysis.Rnw:498-508
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
### code chunk number 25: analysis.Rnw:542-564
###################################################
cmpH5s <- lapply(Sys.glob("../Data/Synthetic/*/data/aligned_reads.cmp.h5"), PacBioCmpH5)
names(cmpH5s) <- sapply(cmpH5s, function(h) basename(dirname(dirname(h@fileName))))
modifications <- list("2x_5mC"  = c(55,74),
                      "2x_5hmC" = c(51,74),
                      "2x_4mC"  = c(55,74),
                      "2x_6mA"  = c(57,68,112,123))

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
### code chunk number 26: analysis.Rnw:568-569
###################################################
distributionPlot(c(modifications[["2x_6mA"]][1], 40), "2x_6mA")


###################################################
### code chunk number 27: analysis.Rnw:580-581
###################################################
distributionPlot(c(modifications[["2x_5mC"]][1], 40), "2x_5mC")


###################################################
### code chunk number 28: analysis.Rnw:604-622
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
  plot(ipdRatios[[a]], pch = 16, xlab = "Position", ylim = c(0, 8), ylab = "IPD Ratio", main = a)
  abline(v = modifications[[a]], col = "grey", lwd = 4)
  points(ipdRatios[[a]], pch = 16)
})


