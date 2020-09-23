FILES=`ls data/*.h5ad 2>/dev/null`  # get names of .h5ad files in data dir

# convert .h5ad files to .csv for reading into R
for f in $FILES;
do
	printf "\nConverting $f to .csv\n"
	cp $f tmp.h5ad
	kitchen transpose tmp.h5ad
	kitchen to_dense tmp.h5ad
	time kitchen to_csv tmp.h5ad -s
	mv tmp_X.csv ${f%.h5ad}.csv
	rm tmp_obs.csv
	rm tmp_var.csv
	rm tmp.h5ad

	# read in above csv file and perform EmptyDrops filtering
	Rscript emptydrops.R ${f%.h5ad}.csv 2000 1000

	# add EmptyDrops labels to original .h5ad file
	printf "\nAdding EmptyDrops labels to $f"
	python emptydrops_to_h5ad.py $f emptydrops_${f%.h5ad}_X.csv
done

printf "\nDone!"
