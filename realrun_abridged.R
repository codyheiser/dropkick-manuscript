suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(source("../simulations/functions.R"))

ppath <- "pics-real"
if(!dir.exists(ppath)){dir.create(ppath)}

# Specifying all the input files.
ALLFILES <- c(#"2771_S1.csv", # labels added
	      #"2771_S2.csv", # labels added
	      #"3072_S1.csv", # labels added
	      #"3072_S2.csv", # labels added
	      #"3072_S3.csv", # labels added
	      #"3254_S1.csv", # labels added
	      #"3777_S1.csv", # labels added
	      #"3777_S2.csv", # labels added
	      #"3907_S1.csv", # labels added
	      #"3907_S2.csv", # labels added
	      #"4058_colon.csv", # labels added
	      #"MPP_03A.csv", # labels added
	      #"MPP_04A.csv", # labels added
	      #"MPP_07A.csv", # labels added
	      #"MPP_08A.csv", # labels added
	      #"MPP_13A.csv", # labels added
	      #"MPP_18A.csv", # labels added
	      #"MPP_19A.csv", # labels added
	      #"MPP_21A.csv", # labels added
	      #"MPP_48A.csv", # labels added
	      #"MPX_03A.csv", # labels added
	      #"MPX_03B.csv", # labels added
	      #"MPX_04A.csv", # labels added
	      #"MPX_08A.csv", # labels added
	      #"MPX_08B.csv", # labels added
	      #"MPX_09A.csv", # labels added
	      #"MPX_09B.csv" # labels added
	      #"MPX_17A.csv", # labels added
	      #"MPX_19A.csv", # labels added
	      #"MPX_19B.csv", # labels added
	      #"MPX_21A.csv", # labels added
	      #"MPX_21B.csv", # labels added
	      #"3659_colon.csv" 
	      #"293t/matrices_mex/hg19/",
	      #"jurkat/matrices_mex/hg19/",
	      #"neuron_9k/raw_gene_bc_matrices/mm10/",
	      #"neurons_900/raw_gene_bc_matrices/mm10/",
	      #"pbmc4k/raw_gene_bc_matrices/GRCh38",
	      #sprintf("placenta%i/placenta%i_raw_gene_bc_matrices/GRCh38/", seq_len(6), seq_len(6)),
	      #"t_4k/raw_gene_bc_matrices/GRCh38/"
	      #"high-background1_X.csv",
	      #"high-background2_X.csv",
	      #"high-background3_X.csv",
	      #"high-background4_X.csv",
	      #"high-background5_X.csv",
	      #"high-background6_X.csv",
	      #"high-background7_X.csv",
	      #"high-background8_X.csv",
	      #"high-background9_X.csv",
	      #"high-background10_X.csv",
	      "low-background1_X.csv",
	      "low-background2_X.csv",
	      "low-background3_X.csv",
	      "low-background4_X.csv",
	      "low-background5_X.csv",
	      "low-background6_X.csv",
	      "low-background7_X.csv",
	      "low-background8_X.csv",
	      "low-background9_X.csv",
	      "low-background10_X.csv"
          )

expected <- c(#2400, #"2771_S1.csv"
	      #1000, #"2771_S2.csv"
	      #1300, #"2771_S3.csv"
	      #1500, #"3072_S1.csv"
	      #2000, #"3072_S2.csv"
	      #2000, #"3072_S3.csv"
	      #1300, #"3254_S1.csv"
	      #1000, #"3777_S1.csv"
	      #1200, #"3777_S2.csv"
	      #2000, #"3907_S1.csv"
	      #2000, #"3907_S2.csv"
	      #2000, #"4058_colon.csv"
	      #900, #MPP_03A
	      #700, #MPP_04A
	      #1300, #MPP_07A
	      #1500, #MPP_08A
	      #1600, #MPP_13A
	      #700, #MPP_18A
	      #1400, #MPP_19A
	      #2000, #MPP_21A
	      #500, #MPP_48A
	      #1100, #MPX_03A
	      #800, #MPX_03B
	      #900, #MPX_04A
	      #1100, #MPX_08A
	      #900, #MPX_08B
	      #500, #MPX_09A
	      #1200 #MPX_09B
	      #800, #MPX_17A
	      #800, #MPX_19A
	      #2700, #MPX_19B
	      #1600, #MPX_21A
	      #1200, #MPX_21B
	      #5000 #3659_colon
	      #2800, #293t
	      #3200, #jurkat
	      #9000, #neuron_9k
	      #900, #neurons_900
	      #4000, #pbmc_4k
	      #rep(5000, 6), #placenta
	      #4000 #t_4k
	      #3000,
	      #3000,
	      #3000,
	      #3000,
	      #3000,
	      #3000,
	      #3000,
	      #3000,
	      #3000,
	      #3000,
	      3000,
	      3000,
	      3000,
	      3000,
	      3000,
	      3000,
	      3000,
	      3000,
	      3000,
	      3000
          )

lower <- c(#500, #"2771_S1.csv"
	   #500, #"2771_S2.csv"
	   #100, #"2771_S3.csv"
	   #500, #"3072_S1.csv"
	   #1000, #"3072_S2.csv"
	   #2000, #"3072_S3.csv"
	   #500, #"3254_S1.csv"
	   #500, #"3777_S1.csv"
	   #500, #"3777_S2.csv"
	   #1000, #"3907_S1.csv"
	   #1000, #"3907_S2.csv"
	   #1000, #"4058_colon.csv"
	   #1000, #MPP_03A
	   #1000, #MPP_04A
	   #500, #MPP_07A
	   #1000, #MPP_08A
	   #1000, #MPP_13A
	   #1000, #MPP_18A
	   #1000, #MPP_19A
	   #5000, #MPP_21A
	   #1000, #MPP_48A
	   #500, #MPX_03A
	   #500, #MPX_03B
	   #2000, #MPX_04A
	   #2000, #MPX_08A
	   #2000, #MPX_08B
	   #1000, #MPX_09A
	   #1000 #MPX_09B
	   #3000, #MPX_17A
	   #3000, #MPX_19A
	   #1000, #MPX_19B
	   #2000, #MPX_21A
	   #2000, #MPX_21B
	   #5000 #3659_colon
	   #100, #293t
	   #100, #jurkat
	   #100, #neuron_9k
	   #100, #neurons_900
	   #100, #pbmc_4k
	   #100, #placenta
	   #100 #t_4k
	   #5000,
	   #5000,
	   #5000,
	   #5000,
	   #5000,
	   #5000,
	   #5000,
	   #5000,
	   #5000,
	   #5000,
	   1000,
	   1000,
	   1000,
	   1000,
	   1000,
	   1000,
	   1000,
	   1000,
	   1000,
	   1000
	   )

for (i in seq_along(ALLFILES)) {
    fname <- ALLFILES[i]
    fpath <- file.path("..", "data", fname)
    if (!file.exists(fpath)) {
        message("missing data files for '", fname, "'")
        next
    }
    
    # timing
    tic(paste0("EmptyDrops analysis on '",fname,"'"))
    
    # read in csv file(s)
    if (tail(strsplit(fname, "\\.")[[1]], 1) == "csv"){
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
    #stub <- sub("/.*", "", fname)
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
    
    # Testing emptyDrops.
    message(paste0("Running EmptyDrops for '",outname,"'"))
    e.out <- emptyDrops(final, lower = metadata(stats)$inflection, test.ambient=T)
    write.csv(e.out, file=paste0("emptydrops_",outname,".csv"))
    e.keep <- e.out$FDR <= 0.001
    e.keep[is.na(e.keep)] <- FALSE
    #write.csv(e.keep, file=paste0("emptydrops_",outname,".csv"))
    
    # Keeping everything above the knee point.
    message(paste0("Calculating knee-point cutoff for '",outname,"'"))
    k.keep <- stats$total >= stats$knee
    write.csv(k.keep, file=paste0("knee_",outname,".csv"))
    
    # Using the CellRanger approach.
    message(paste0("Running CellRanger v2.1 for '",outname,"'"))
    c.keep <- defaultDrops(final, expected=expected[i], lower.prop=0.2)
    write.csv(c.keep, file=paste0("cellranger_",outname,".csv"))
    
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
