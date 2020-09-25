![Alt text](data/dropkick_logo.png)

Scripts and notebooks for reproducing analysis from the 2020 dropkick manuscript.

---

1. Run [`install.sh`](install.sh) to install `dropkick` and dependencies for completing downstream analysis. It is recommended to initiate a new virtual python environment to avoid conflicts.

2. To add EmptyDrops and CellRanger v2 labels to all raw `.h5ad` files in [`data/`](data/), run [`emptydrops_cellranger.sh`](emptydrops_cellranger.sh). This will deposit labels into the original `.h5ad` file.

3. [`run_and_test.sh`](run_and_test.sh) will perform `dropkick` filtering on all `.h5ad` files in [`data/`](data/) and then compare outputs to `CellRanger_2` and `EmptyDrops` labels.

4. [`manualfilter_example.ipynb`](manualfilter_example.ipynb) outlines our manual filtering approach performed on inDrop sequencing samples for comparison to `dropkick`.

---

**NOTE:** This entire process may take several minutes to complete on a well-equipped machine. The [`emptydrops_cellranger.sh`](emptydrops_cellranger.sh) script typically completes in ~5 min per input file (`.h5ad`). The [`run_and_test.sh`](run_and_test.sh) script may also take a couple minutes for dimension reduction, clustering, and embedding to visualize `dropkick` vs. `EmptyDrops` and `CellRanger` filtering.
