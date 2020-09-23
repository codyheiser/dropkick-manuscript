![Alt text](data/dropkick_logo.png)

Scripts and notebooks for reproducing analysis from the 2020 dropkick manuscript.

1. Run [`install.sh`](install.sh) to install `dropkick` and dependencies for completing downstream analysis. It is recommended to initiate a new virtual environment to avoid conflicts.

2. To add EmptyDrops labels to the raw `.h5ad` files, run [`emptydrops.sh`](emptydrops.sh). This will first convert your AnnData files into `.csv` format, read them into R, run EmptyDrops analysis, and then deposit labels into the original `.h5ad` file. Ensure you have [`argparse`](https://www.rdocumentation.org/packages/argparse/versions/2.0.1), [`DropletUtils`](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html), [`Matrix`](https://www.rdocumentation.org/packages/Matrix/versions/1.2-18), and [`tictoc`](https://www.rdocumentation.org/packages/tictoc/versions/1.0) installed in R.

3. [`run_and_test.sh`](run_and_test.sh) will perform CellRanger v2 and dropkick filtering on all `.h5ad` files in [`data/`](data/) and then compare `dropkick` outputs to `CellRanger_2` and `EmptyDrops` labels.

4. [`manualfilter_example.ipynb`](manualfilter_example.ipynb) outlines our manual filtering approach performed on inDrop sequencing samples for comparison to `dropkick`.
