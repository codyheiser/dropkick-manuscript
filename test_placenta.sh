# combine six replicates of placenta data and reduce dimensions to visualize results
printf "\nCombining placenta replicates and reducing dimensions:\n"
kitchen concatenate placenta1_dropkick.h5ad placenta2_dropkick.h5ad placenta3_dropkick.h5ad placenta4_dropkick.h5ad placenta5_dropkick.h5ad placenta6_dropkick.h5ad -o placenta_combined_dropkick.h5ad

kitchen recipe placenta_combined_dropkick.h5ad -s dropkick_label CellRanger_2 EmptyDrops cellbender -p --paga -cc -sa -c batch CellRanger_2 EmptyDrops cellbender dropkick_label dropkick_score arcsinh_total_counts pct_counts_mito -de dotplot
