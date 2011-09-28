
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
