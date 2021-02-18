![Alt text](data/dropkick_logo.png)

Scripts and notebooks for reproducing analysis from the [2020 dropkick manuscript](https://doi.org/10.1101/2020.10.08.332288).

---

1. Run [`install.sh`](install.sh) to install `dropkick` and dependencies for completing downstream analysis. It is recommended to initiate a new virtual python environment to avoid conflicts.

2. (optional) Download six replicates of [human placenta data](https://www.nature.com/articles/s41586-018-0698-6) and perform CellBender remove-background by following instructions in the [`cellbender/` directory](cellbender/).

3. (optional) Build a high-background simulation from human PBMC data by following instructions in the [`simulation/` directory](simulation/).

4. To add EmptyDrops and CellRanger v2 labels to all raw `.h5ad` files in [`data/`](data/), run [`emptydrops_cellranger.sh`](emptydrops_cellranger.sh). This will deposit labels into the original `.h5ad` file.

5. [`run_and_test.sh`](run_and_test.sh) will perform `dropkick` filtering on all `.h5ad` files in [`data/`](data/) and then compare outputs to `CellRanger_2` and `EmptyDrops` labels.

6. [`test_3907.sh`](test_3907.sh) will combine the 3907 human colorectal cancer datasets and reduce dimensions, showing differences between filtering tools.

7. [`test_placenta.sh`](test_placenta.sh) will combine the [human placenta datasets](https://www.nature.com/articles/s41586-018-0698-6) downloaded in **step 2** and reduce dimensions.

8. [`test_simulation.sh`](test_simulation.sh) will reduce dimensions for the high-background simulation made in **step 3** above.

9. [`manualfilter_example.ipynb`](manualfilter_example.ipynb) outlines our manual filtering approach performed on inDrop sequencing samples for comparison to `dropkick`.

---

**NOTE:** This entire process may take several minutes to complete on a well-equipped machine. The [`emptydrops_cellranger.sh`](emptydrops_cellranger.sh) script typically completes in ~5 min per input file (`.h5ad`). The [`run_and_test.sh`](run_and_test.sh) script may also take a couple minutes for dimension reduction, clustering, and embedding to visualize `dropkick` vs. `EmptyDrops` and `CellRanger` filtering.
