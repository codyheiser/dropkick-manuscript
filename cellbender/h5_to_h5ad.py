from glob import glob
import scanpy as sc

for f in glob("*.h5"):
    print("Reading {}".format(f))
    a = sc.read_10x_h5(f)
    a.var_names_make_unique()
    print("Writing {}ad".format(f))
    a.write("{}ad".format(f), compression="gzip")
