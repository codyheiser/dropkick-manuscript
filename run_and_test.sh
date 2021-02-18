FILES=`ls data/*.h5ad 2>/dev/null` # get names of all .h5ad files from data dir

for f in $FILES;
do
	file=`basename $f` # get name of .h5ad file

	printf "\nStarting dropkick qc on $file:\n"
	time dropkick qc $f  # generate dropkick QC report

	printf "\nRunning dropkick filtering on $file:\n"
	time dropkick run $f -j 5  # run dropkick filtering with 5 jobs
	NAME=`basename $f`
	kitchen label_info ${NAME%.h5ad}_dropkick.h5ad -l dropkick_label  # show number of cells identified by dropkick

done

# aggregate stats comparing dropkick to EmptyDrops and CellRanger
printf "\nSummarizing statistics comparing dropkick to EmptyDrops and CellRanger:\n"
python dropkick_agg_stats.py `ls | grep .h5ad` -l CellRanger_2 EmptyDrops
