# required packages for emptydrops_run.R
install.packages("argparse")
install.packages("Matrix")
install.packages("tictoc")
# DropletUtils package (containing EmptyDrops function)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DropletUtils")
