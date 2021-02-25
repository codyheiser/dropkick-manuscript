High-Background Simulation of Human PBMCs
---

Using a human 4k PBMC dataset from 10x Genomics, we'd like to add high-counts empty droplets simulated from the ground-truth empty droplets native to the low-background dataset.

To do so, we sample ambient genes from a multinomial distribution up to 5k total counts, creating a counts vector that imitates a high-background barcode. Then these barcodes are added back to the original counts matrix with a label.

The following command will run this simulation and save the result to `../data/pbmc_4k_sim.h5ad`.

```sh
python simulate_empty.py ../data/pbmc/pbmc_4k.h5ad --n-bottom-counts 100 --n-ambient 1000 --n-empties 2000 --max-counts-empty 5000 -o ../data
```
