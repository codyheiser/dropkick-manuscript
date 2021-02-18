# combine two replicates (3907_S1 & 3907_S2) and reduce dimensions to visualize results
printf "\nCombining 3907_S1 and 3907_S2 and reducing dimensions:\n"
kitchen concatenate 3907_S1_dropkick.h5ad 3907_S2_dropkick.h5ad -o 3907_combined_dropkick.h5ad

kitchen recipe 3907_combined_dropkick.h5ad -s dropkick_label CellRanger_2 EmptyDrops -p --paga -cc -sa -c batch CellRanger_2 EmptyDrops dropkick_label dropkick_score arcsinh_total_counts arcsinh_n_genes_by_counts pct_counts_mito pct_counts_ambient -de dotplot
