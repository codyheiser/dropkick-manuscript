FILES=`ls data/*.h5ad 2>/dev/null` # get names of all .h5ad files from data dir

for f in $FILES;
do
	file=`basename $f` # get name of .h5ad file

	printf "\nStarting dropkick qc on $file\n"
	time dropkick qc $f

	printf "\nRunning dropkick filtering on $file\n"
	time dropkick run $f -j 5

	printf "\nRunning cNMF on $file\n"
	cnmf_p ${file%.h5ad}_dropkick.h5ad --name ${file%.h5ad}_union -k 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 -n 30 -j 5 --subset dropkick_label CellRanger_2 EmptyDrops --auto-k --cleanup

	# move plots into _union directories with the rest of outputs
	mv ${file%.h5ad}_coef.png ${file%.h5ad}_union/${file%.h5ad}_coef.png
	mv ${file%.h5ad}_score.png ${file%.h5ad}_union/${file%.h5ad}_score.png
	mv ${file%.h5ad}_qc.png ${file%.h5ad}_union/${file%.h5ad}_qc.png

	printf "\nReducing dimensions and plotting embeddings for $file\n"
	kitchen recipe ${file%.h5ad}_union/${file%.h5ad}_union*.h5ad -p --paga -cc -sa --cnmf -o ${file%.h5ad}_union --paga -cc -c CellRanger_2 EmptyDrops dropkick_label dropkick_score arcsinh_total_counts arcsinh_n_genes_by_counts pct_counts_mito pct_counts_ambient -de dotplot

done

# aggregate stats comparing dropkick to EmptyDrops and CellRanger
python dropkick_agg_stats.py `ls | grep .h5ad` -l CellRanger_2 EmptyDrops
