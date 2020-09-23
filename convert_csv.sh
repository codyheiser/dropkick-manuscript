FILES=`ls $1/*.h5ad 2>/dev/null`  # get names of .h5ad files in $1 dir

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
done

