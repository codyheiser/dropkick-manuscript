# reduce dimensions for simulated PBMC data
kitchen recipe pbmc_4k_sim_dropkick.h5ad -s dropkick_label CellRanger_2 EmptyDrops -p --paga -cc -sa -c simulated CellRanger_2 EmptyDrops dropkick_label dropkick_score arcsinh_total_counts pct_counts_mito -de dotplot
