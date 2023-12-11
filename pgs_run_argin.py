import subprocess
import argparse
import pandas as pd
import glob
import gzip
import re

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-l", "--list_csv", type = str, required = True, help = "the path of the input list of PRS Catalog scores and their corresponding phenotypes", default = "pgs_all_metadata_scores.csv")
parser.add_argument("-i", "--input", type = str, required = True, help = "an input file in plink bed/bim/fam format")
parser.add_argument("-id", "--pgs_id", type = str, required = True, help = "the id of the PGS Catalog score you wish to apply")
parser.add_argument("-po", "--prs_out_dir", type = str, required = True, help = "the output directory to which the PRS files will be saved", default = ".")
parser.add_argument("-so", "--score_out_dir", type = str, required = True, help = "the output score directory", default = "./")
parser.add_argument("-sr", "--score_out_root", type = str, required = False, help = "the output score root name")
parser.add_argument("-b", "--genome_build", type = str, required = True, help = "the genome build of the output file. Options: GRCh37 or GRCh38")
parser.add_argument("-t", "--threads", type = str, required = True, help = "the number of threads to use for PRS calculation", default = 1)
args = parser.parse_args()

#Define variant reformatting function
def format_variant(variant):
    	return re.sub(r'(\d+)\.0_', r'\1_', variant)

df = pd.read_csv(args.list_csv)
#pgs_ids = df.iloc[:,0].tolist()

# for id in pgs_ids:
# 	print(id)
# 	subprocess.run(["download_scorefiles", "-i", id, "-o", args.out_dir, "-b", args.genome_build])

#Download the PRS
subprocess.run(["download_scorefiles", "-i", args.pgs_id, "-o", args.prs_out_dir, "-b", args.genome_build])

#Find the index of the PGS ID in the DataFrame
pgs_id_index = df[df.iloc[:,0] == args.pgs_id].index[0]

#Use the index to find the corresponding full name
full_name = df.loc[pgs_id_index, "Reported Trait"]
full_name = full_name.replace(" ", "_")

#Find the score file with a name starting with the PGS ID
score_file_list = glob.glob(args.prs_out_dir + "/" + args.pgs_id + "*.txt.gz")

#Use the first matching file as the score file
if score_file_list:
    score_file = score_file_list[0]
    
    # Create new file path
    score_df = pd.read_csv(score_file, compression='gzip', sep='\t', comment='#', header=0)
    score_df['varid'] = "chr" + score_df['hm_chr'].astype(str) + "_" + score_df['hm_pos'].astype(str) + "_" + score_df['effect_allele']
    
    #Reformat the variant IDs if they are formatted incorrectly
    score_df['varid'] = score_df['varid'].apply(format_variant)

    out_df = score_df[['varid', 'effect_allele', 'effect_weight']]
    out_df.to_csv((args.prs_out_dir + "/" + args.pgs_id + ".reformat.txt"), sep='\t', index=False)
    
    subprocess.run(["rm", score_file])
    
    if args.score_out_root is None:
    	subprocess.run(["plink2", "--bfile", args.input, "--score", (args.prs_out_dir + "/" + args.pgs_id + ".reformat.txt"), "1", "2", "3", "ignore-dup-ids", "cols=+scoresums", "header", "--threads", args.threads, "--allow-extra-chr", "--out", (args.score_out_dir + "/" + args.pgs_id + "_" + full_name)])
    else:
    	subprocess.run(["plink2", "--bfile", args.input, "--score", (args.prs_out_dir + "/" + args.pgs_id + ".reformat.txt"), "1", "2", "3", "ignore-dup-ids", "cols=+scoresums", "header", "--threads", args.threads, "--allow-extra-chr", "--out", (args.score_out_dir + "/" + args.score_out_root + "_" + args.pgs_id + "_" + full_name)])
else:
    print(f"No score file found starting with the PGS ID: {args.pgs_id}")
    
# open the output file
if args.score_out_root is None:
	out_file = args.score_out_dir + "/" + args.pgs_id + "_" + full_name + ".sscore"
else:
	out_file = args.score_out_dir + "/" + args.score_out_root + "_" + args.pgs_id + "_" + full_name + ".sscore"

df_out = pd.read_csv(out_file, delim_whitespace=True)  # assumes the file is whitespace-delimited

# transform the 'score' column to percentiles
df_out['score_pct'] = df_out['SCORE1_SUM'].rank(pct=True)

# save the transformed data
df_out.to_csv(out_file, sep='\t', index=False)