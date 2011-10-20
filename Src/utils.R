
makeSummaryTable <- function(dataSet = "Synthetic", caption = "", ...) {
  cmpH5s <- lapply(Sys.glob(paste("../Data/",dataSet,"/*/data/aligned_reads.cmp.h5", sep = "")), PacBioCmpH5)
  names(cmpH5s) <- sapply(cmpH5s, function(h) basename(dirname(dirname(h@fileName))))
  d <- do.call(rbind, lapply(cmpH5s, function(a) {
    data.frame(nAlignments = nrow(a), 
               nMolecules  = length(unique(interaction(a$movieName, a$holeNumber, drop = T))),
               nMovies     = length(unique(a$movieName)),
               nReferences = length(unique(a$fullRefName)))
  }))
  xtable(d, caption = caption, ...)
}

plotIPDForReads <- function(cmpH5, idx, matches = FALSE, range = c(1, 10), whStrand = 0, whPos = floor((range[1]+range[2])/2),
                            maxReads = 10, ...) {
  tpos <- getByTemplatePosition(cmpH5, idx = idx, function(cmpH5, idx) {
    mapply(getIPD(cmpH5, idx), getPulseWidth(cmpH5, idx), FUN = function(ipd, pw) {
      data.frame(ipd = ipd, pw = pw)
    }, SIMPLIFY = FALSE)
  })
  tpos <- subset(tpos, position >= range[1] & position <= range[2] & strand == whStrand)
    
  plotStrip <- function(d) {
    ipd <- d$elt.ipd[!is.na(d$elt.ipd)]
    pw <- d$elt.pw[!is.na(d$elt.pw)]
  
    ipd.end <- cumsum(ipd) + cumsum(pw) - pw
    pw.end <- cumsum(ipd) + cumsum(pw)
    
    colors <- rep(c("black","red","black"), c(2,4,4))
    plot(NA, xlim=c(-0.6, pw.end[length(pw.end)]), ylim = c(-.75,1.5), axes = F, xlab = "", ylab = "")
    segments(c(0, pw.end[-nrow(d)]), 0, ipd.end, 0, col="darkgreen")
    segments(ipd.end, .75, pw.end, 0.75, col="blue")
    segments(ipd.end, 0, ipd.end, 0.75, col="blue")
    segments(pw.end, 0.75, pw.end, 0, col="blue")
       
    ipdG <- d$elt.ipd
    pwG <- d$elt.pw
    ipdG[is.na(ipdG)] <- 0
    pwG[is.na(pwG)] <- 0
    m <- cbind(cumsum(ipdG) + cumsum(pwG) - pwG,
               cumsum(ipdG) + cumsum(pwG))
    text((m[,1]+m[,2])/2, jitter(rep(1, nrow(m))), d$read)
    text((m[,1]+m[,2])/2, jitter(rep(-.5, nrow(m))), d$ref)
    text(-0.5, 1, "read")
    text(-0.5, -.5, "reference")
    msk <- which(whPos == d$position & d$read == d$ref)
    if (pw.end[msk-1] == ipd.end[msk]) {
      points(pw.end[msk-1], 0.25, pch = 'X', col = 'red')
    } else {
      arrows(pw.end[msk-1], 0.25, ipd.end[msk], 0.25, col = "red", length=0.1)
    }
  }
  
  byRead <- split(tpos, tpos$idx)
  if (length(byRead) > maxReads) {
    byRead <- byRead[1:maxReads]
  }
  if (matches) {
    byRead <- lapply(byRead, function(a) subset(a, read == ref))
  }
  par(mfrow = c(length(byRead), 1), mar=c(4, 1, 1, 1))
  invisible(lapply(byRead, plotStrip))
}

makeTopTable <- function(treatmentCmpH5, controlCmpH5, start = 1, end = NA, whichStrand = 0,
                         getData = getByTemplatePosition, whReference = 1,
                         test = wilcox.test) {
  stopifnot(getRefName(treatmentCmpH5, whReference) ==
            getRefName(controlCmpH5, whReference))
  end <- if (is.na(end)) getRefLength(treatmentCmpH5, whReference) else end
  gData <- function(cmpH5) {
    reads <- getReadsInRange(cmpH5, whReference, start, end)
    s <- subset(getData(cmpH5, idx = reads),
           position >= start & position <= end & strand == whichStrand & read == ref)
    split(s, factor(s$position, start:end))
  }
  tDta <- gData(treatmentCmpH5)
  cDta <- gData(controlCmpH5)
  dta <- do.call(rbind, mapply(tDta, cDta, FUN = function(tD, cD) {
    res <- test(tD$elt, cD$elt)
    data.frame(position = tD$position[1],
               reference = as.character(tD$ref[1]),
               read = as.character(tD$read[1]),
               p.value = res$p.value,
               statistic = res$statistic,
               n.control = nrow(cD),
               n.treatment = nrow(tD),
               ipd.ratio = mean(tD$elt, na.rm = T)/mean(cD$elt, na.rm = T))
  }, SIMPLIFY = FALSE))
  dta$fdr <- p.adjust(dta$p.value, "fdr")
  dta[order(dta$p.value), ]
}

trimmedSlog <- function(x, trim = .975, alpha = 1/100) {
  log(x[x<quantile(x, trim)] + alpha)
}

##
## This function demonstrates testing using permutations. The primary
## problem with this current implementation is performance
##
makeTopTablePermutation <- function(treatmentCmpH5, controlCmpH5, start = 1, end = NA, whichStrand = 0,
                                    getData = getByTemplatePosition, whReference = 1, N = 100, 
                                    statistic = function(x, y) {
                                      mean(trimmedSlog(x)) - mean(trimmedSlog(y))
                                    }) {
  stopifnot(getRefName(treatmentCmpH5, whReference) ==
            getRefName(controlCmpH5, whReference))
  end <- if (is.na(end)) getRefLength(treatmentCmpH5, whReference) else end
  gData <- function(cmpH5) {
    reads <- getReadsInRange(cmpH5, whReference, start, end)
    s <- subset(getData(cmpH5, idx = reads),
                position >= start & position <= end & strand == whichStrand & read == ref)
    ## instead of splitting by position we split by read idx. 
    lapply(split(s, s$idx), function(x) cbind(x$elt, x$position))
  }
  tDta <- gData(treatmentCmpH5)
  cDta <- gData(controlCmpH5)
  labels <- c(rep("treatment", length(tDta)), rep("control", length(cDta)))
  jDta <- c(tDta, cDta)
  
  computeStatsForPositions <- function(labels) {
    tD <- do.call(rbind, jDta[labels == "treatment"])
    cD <- do.call(rbind, jDta[labels == "control"])
    mapply(split(tD[,1], factor(tD[,2], start:end)),
           split(cD[,1], factor(cD[,2], start:end)), FUN = function(x, y) {
             statistic(x, y)
           })
  }
  observed <- computeStatsForPositions(labels)
  R <- replicate(N, {
    computeStatsForPositions(sample(labels))
  })
  p.approx <- 1 - rowMeans(observed >= R)
  p.approx <- ifelse(p.approx == 0, 1/(N + 1), p.approx)
  d <- data.frame(p.approx = p.approx, statistic = observed, position = start:end)
  d[order(d$p.approx), ]
}
