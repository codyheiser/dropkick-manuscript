suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tictoc))

ppath <- "emptydrops_out"
if(!dir.exists(ppath)){dir.create(ppath)}

plotBarcodes <- function(ranks, totals, fitted=NULL, subset=NULL, ...) {
    xlim <- range(ranks[ranks>0])
    ylim <- range(totals[totals>0])

    # Dropping non-unique points to save space.
    # Subsetting performed after range() to create comparable plots, if desired.
    keep <- !duplicated(totals)
    if (!is.null(subset)) {
        alt.keep <- keep
        keep[] <- FALSE
        keep[subset] <- alt.keep[subset]
    }
    Rx <- ranks[keep]
    Tx <- totals[keep]

    # Actually making the plot, also plotting the fitted values if requested.
    plot(Rx, Tx, log="xy", xlab="Rank", ylab="Total count", xlim=xlim, ylim=ylim, ...)
    if (!is.null(fitted)) {
        Fx <- fitted[keep]
        o <- order(Rx)
        lines(Rx[o], Fx[o], col="red", lwd=2)
    }
    return(invisible(NULL))
}

saveRaster <- function(...) {
    tmp <- tempfile(fileext=".png")
    png(tmp, width=7, height=7, units="in", res=300, pointsize=12)
    par(mar=c(0,0,0,0))
    options(bitmapType="cairo")
    plot(..., xlab="", ylab="", axes=FALSE)
    dev.off()
    tmp
}

loadRaster <- function(fpath, log=FALSE) {
    limits <- par("usr")
    img <- png::readPNG(fpath)
    if (log) { limits <- 10^limits }
    rasterImage(img, limits[1], limits[3], limits[2], limits[4])
    invisible(limits)
}

parser <- ArgumentParser(description="Run EmptyDrops analysis on scRNA-seq counts")
parser$add_argument("counts", type="character", nargs='+', help="Counts file as .csv with no headers")
parser$add_argument("lower", type="integer", nargs='+', help="Lower limit of counts to consider initial empty droplets")
args <- parser$parse_args()

# Specifying all the input files.
ALLFILES <- args$counts
lower <- args$lower

for (i in seq_along(ALLFILES)) {
    fpath <- ALLFILES[i]
    fname <- tail(strsplit(as.character(fpath), "/")[[1]], n=1)
    
    if (!file.exists(fpath)) {
        message("Missing data files for '", fname, "'")
        next
    }
    
    # timing
    tic(paste0("EmptyDrops analysis on '",fname,"'"))
    
    # read in csv file(s)
    if (tail(strsplit(as.character(fname), "\\.")[[1]], 1) == "csv"){
        message(paste0("Reading in data for '",fname,"'"))
        x <- read.csv(fpath, header=F)
        x2 <- as(as.matrix(x), "sparseMatrix")
        sce <- SingleCellExperiment(assays=list(counts=x2))
        outname <- strsplit(fname, "\\.")[[1]][1]
        
    # read in 10X .mtx files
    }else{
        message(paste0("Reading in data for '",fname,"'"))
        sce <- read10xCounts(fpath)
        outname <- strsplit(fname, "/")[[1]][1]
    }
    
    message(paste0("Calculating global counts metrics for '",outname,"'"))
    final <- counts(sce)
    stats <- barcodeRanks(final, lower=lower[i])
    
    # Making a barcode plot for diagnostics.
    message(paste0("Plotting barcode log-rank for '",outname,"'"))
    pdf(file.path(ppath, paste0("rank_", outname, ".pdf")))
    plotBarcodes(stats$rank, stats$total, fitted=stats$fitted, cex.axis=1.2, cex.lab=1.4, main=outname, col="grey70", pch=16)
    abline(h=metadata(stats)$knee, col="dodgerblue", lty=2, lwd=2)
    abline(h=metadata(stats)$inflection, col="forestgreen", lty=2, lwd=2)
    legend("bottomleft", lty=c(1,2,2), col=c("red", "dodgerblue", "forestgreen"),
           legend=c("Fitted spline", "Knee", "Inflection"), lwd=2, cex=1.2)
    dev.off()
    
    # EmptyDrops
    message(paste0("Running EmptyDrops for '",outname,"'"))
    e.out <- emptyDrops(final, lower = metadata(stats)$inflection, test.ambient=T)
    write.csv(e.out, file=file.path(ppath, paste0("emptydrops_",outname,".csv")))
    e.keep <- e.out$FDR <= 0.001
    e.keep[is.na(e.keep)] <- FALSE
    
    ############################
    # Examining the distribution of deviances.
    X <- e.out$Total/1000
    Y <- -e.out$LogProb/1000
    xlim <- range(X[!is.na(e.out$LogProb)]) 
    tmp <- saveRaster(X, Y, pch=16, cex=0.5, col=ifelse(e.keep, "red", "grey80"), log="xy", xlim=xlim)
    
    message(paste0("Plotting deviances for '",outname,"'"))
    pdf(file.path(ppath, paste0("dev_", outname, ".pdf")))
    par(mar=c(5.1, 5.1, 4.1, 1.1))
    
    plot(X, Y, log="xy", 
         xlab=expression("Total count (x "*10^3*")"), 
         ylab=expression("-Log likelihood (x "*10^3*")"), 
         xlim=xlim, type="n", cex.axis=1.2, cex.lab=1.4,
         main=outname, cex.main=1.4)
    loadRaster(tmp, log=TRUE)
    box()
    
    legend("bottomright", col=c("red", "grey80"), pch=16, cex=1.2,
           legend=c("Putative cell", "Empty droplet"))
    dev.off()
    
    # timing
    toc()
}
