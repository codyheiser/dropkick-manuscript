pip install numpy  # only dependency prior to installing dropkick
pip install dropkick  # install the latest version of the dropkick package
pip install progressbar  # needed for dropkick_agg_stats.py script
git clone https://github.com/codyheiser/cnmf.git ../cnmf && pip install -e ../cnmf  # install cNMF package from source
git clone https://github.com/codyheiser/kitchen.git ../kitchen && pip install -e ../kitchen  # install kitchen package from source
Rscript emptydrops_install.R
