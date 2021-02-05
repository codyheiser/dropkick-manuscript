FILES=`ls data/*.h5ad 2>/dev/null`  # get names of .h5ad files in data dir

# convert .h5ad files to .csv for reading into R
for f in $FILES;
do
	printf "\nConverting $f to .csv for EmptyDrops analysis in R:\n"
	cp $f tmp.h5ad
	kitchen transpose tmp.h5ad
	kitchen to_dense tmp.h5ad
	time kitchen to_csv tmp.h5ad -s
	mv tmp_X.csv ${f%.h5ad}.csv
	rm tmp_obs.csv
	rm tmp_var.csv
	rm tmp.h5ad

	# read in above csv file and perform EmptyDrops filtering
	printf "\nPerforming EmptyDrops analysis on $f in R:\n"
	Rscript emptydrops_run.R ${f%.h5ad}.csv 1000

	# add EmptyDrops labels to original .h5ad file
	printf "\nAdding EmptyDrops labels to $f:\n"
	NAME=`basename $f`
	python emptydrops_add.py $f emptydrops_out/emptydrops_${NAME%.h5ad}.csv
	kitchen label_info $f -l EmptyDrops

	# add CellRanger v2 labels to original .h5ad file
	printf "\nAdding CellRanger_2 labels to $f:\n"
	kitchen cellranger2 $f --expected 2000 -lp 0.22
done

printf "\nDone!"
