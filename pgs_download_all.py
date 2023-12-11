import subprocess
import argparse
import pandas as pd

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-l", "--list_csv", type = str, required = True, help = "the path of the input list of PRS Catalog scores and their corresponding phenotypes", default = "pgs_all_metadata_scores.csv")
parser.add_argument("-o", "--out_dir", type = str, required = True, help = "the output directory to which the PRS files will be saved", default = ".")
parser.add_argument("-b", "--genome_build", type = str, required = True, help = "the genome build of the output file. Options: GRCh37 or GRCh38")
args = parser.parse_args()

df = pd.read_csv(args.list_csv)
pgs_ids = df.iloc[:,0].tolist()

for id in pgs_ids:
	print(id)
	subprocess.run(["download_scorefiles", "-i", id, "-o", args.out_dir, "-b", args.genome_build])