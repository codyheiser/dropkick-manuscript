printf "\nInstalling dropkick and dependencies in Python:\n"
pip install numpy  # only dependency prior to installing dropkick
pip install dropkick  # install the latest version of the dropkick package
pip install progressbar  # needed for dropkick_agg_stats.py script

printf "\nInstalling the kitchen for easy manipulation of anndata objects from command line (https://github.com/codyheiser/kitchen):\n"
git clone https://github.com/codyheiser/kitchen.git ../kitchen && pip install -e ../kitchen  # install kitchen package from source

printf "\nInstalling EmptyDrops and dependencies in R:\n"
Rscript emptydrops_install.R
