FILES=`ls data/*.h5ad 2>/dev/null` # get names of all .h5ad files from data dir

for f in $FILES;
do
	file=`basename $f` # get name of .h5ad file

	printf "\nStarting dropkick qc on $file:\n"
	time dropkick qc $f

	printf "\nRunning dropkick filtering on $file:\n"
	time dropkick run $f -j 5

done

# aggregate stats comparing dropkick to EmptyDrops and CellRanger
printf "\nSummarizing statistics comparing dropkick to EmptyDrops and CellRanger:\n"
python dropkick_agg_stats.py `ls | grep .h5ad` -l CellRanger_2 EmptyDrops

# combine two replicates (3907_S1 & 3907_S2) and reduce dimensions to visualize results
printf "\nCombining 3907_S1 and 3907_S2 and reducing dimensions:\n"
kitchen concatenate 3907_S1_dropkick.h5ad 3907_S2_dropkick.h5ad -o 3907_combined_dropkick.h5ad

kitchen recipe 3907_combined_dropkick.h5ad -p --paga -cc -sa -c batch CellRanger_2 EmptyDrops dropkick_label dropkick_score arcsinh_total_counts arcsinh_n_genes_by_counts pct_counts_mito pct_counts_ambient -de dotplot
