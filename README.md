# prs-pgs-calculator
A set of scripts and accompanying instruction for how to compute polygenic risk scores (PRS) across all scores in the PGS Catalog and calculate risk percentile relative to gnomAD samples.

**Disclaimer: The code in this repository is intended for purely educational use and is not intended to diagnose or treat any condition. Polygenic risk scores have a number of notable flaws (some of which are detailed here and here) that limit their potential for use in a clinical setting.**

**Additionally, the example workflow depicted here show calculates PRS relative to publicly available, de-identified whole genome sequences from gnomAD, which have unknown phenotypes. Thus, it is impossible to compare the polygenic risk profile of target samples to that of actually affected reference individuals versus healthy controls for a given phenotype.**

## Dependencies
The following dependencies are required to run the PRS calculation workflow: 
* ```Python``` (version 3.10)
* ```R``` (version 4.2)
* ```R optparse``` (version 1.20.4)
* ```R data.table``` (version 1.14.10)
* ```pandas``` (version 1.5.3)
* ```pgscatalog-utils``` (version 0.4.3)
* ```plink``` (version 2.0)
* ```bcftools``` (version 1.8)
* ```tabix``` (version 0.2.5)

These dependencies can be installed together via Anaconda using the ```prscatalog_env.yml``` file with the following command:
```
conda env create -f prscatalog_env.yml
```

## Downloading PRS scoring files
First, download a CSV file of the metadata and accession IDs for all PRSs submitted to the PGS Catalog. A copy of this file is also saved in this repository as ```pgs_all_metadata_scores.csv```:
```
wget https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata_scores.csv
```

The scores can then be downloaded using the included ```pgs_download_all.py``` script: 
```
conda activate prscatalog_env

python ~/pgs_download_all.py \
-l /home/mbetti/pgs_all_metadata_scores.csv \
-o /home/mbetti/personal_genome/pgs_scores_hg38 \
-b GRCh38
```
The following arguments are required:
* ```-l```: the path to the CSV file with accession IDs
* ```-o```: The path of the desired output directory to which score files will be downloaded
* ```-b```: The desired reference genome build (GRCh37 and GRCh38 supported)

## Preparing the target dataset
The target dataset (containing individuals for which you wish to compute PRS) need to undergo a series of quality control steps so that they are compatible with scores from the PGS Catalog.

Remove multiallelic variants from the target VCF:
```
bcftools view --max-alleles 2 --exclude-types indels /home/mbetti/personal_genome/pgs_catalog/chr_mb_vcfs_hg38/60820188479382_GFX0435435_MB.recal.snp.indel.hg38.anno.vcf.gz > /home/mbetti/personal_genome/pgs_catalog/chr_mb_vcfs_hg38/60820188479382_GFX0435435_MB.recal.snp.indel.hg38.anno.no_multiallelic.vcf
```

Use ```plink2``` to convert the target VCF to a set of plink bim/bed/fam files:
```
plink2 \
--vcf /home/mbetti/personal_genome/pgs_catalog/chr_mb_vcfs_hg38/60820188479382_GFX0435435_MB.recal.snp.indel.hg38.anno.no_multiallelic.vcf \
--make-bed \
--out /home/mbetti/personal_genome/pgs_catalog/chr_mb_vcfs_hg38/60820188479382_GFX0435435_MB.recal.snp.indel.hg38.anno.no_multiallelic \
--allow-extra-chr
```

At this point, you can also modify or reformat the FIDs and IIDs in the target fam file, if so desired.

The SNPs in the PGS Catalog score files have the format ```chr1_13813_G``` (```chr_pos_alt```). So the SNPs in the target bim file should be reformatted to use this same format:
```
library("data.table")

#Declare the path of the bim file
bim_path <- "/home/mbetti/personal_genome/pgs_catalog/chr_mb_vcfs_hg38/60820188479382_GFX0435435_MB.recal.snp.indel.hg38.anno.no_multiallelic.bim"

#Open the bim file as a data frame
bim_file <- fread(bim_path, header = FALSE, sep = "\t", quote = "")
bim_df <- as.data.frame(bim_file)

#Generate new variant IDs for each SNP using the coordinate and allele information in the other columns
bim_df[,2] <- paste(paste0("chr", bim_df[,1]), bim_df[,4], bim_df[,5], sep = "_")

#Write the modified bim out
write.table(bim_df, file = bim_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

## Preparing the reference dataset

The largest publicly available dataset of individual-level whole genome sequences is made available by gnomAD, consisting of 2,500 samples from the 1000 Genomes Project and an additional 780 from the Human Genome Diversity Project, representing over 60 global populations. These data are based on genome build GRCh38.

This gnomAD dataset can be downloaded using the following commands:
```
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr1.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr2.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr3.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr4.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr5.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr6.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr7.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr8.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr9.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr10.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr11.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr12.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr13.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr14.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr15.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr16.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr17.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr18.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr19.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr20.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr21.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr22.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chrX.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz

wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr1.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr2.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr3.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr4.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr5.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr6.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr7.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr8.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr9.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr10.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr11.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr12.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr13.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr14.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr15.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr16.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr17.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr18.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr19.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr20.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr21.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr22.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chrX.vcf.bgz.tbi
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz.tbi
```

*Note: These data are very large (i.e. >3 TB in size), so expect them to both take a long time to download and require a large amound of storage.*

After downloading, these files should be converted to plink format:
```
mkdir /mnt/16tb/gnomad/plink

for i in {1..22}; do
sbatch \
--job-name=chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--wrap="plink2 \
--vcf /mnt/16tb/gnomad/gnomad.genomes.v3.1.2.hgdp_tgp.chr$i\.vcf.bgz \
--make-bed \
--out /mnt/16tb/gnomad/plink/gnomad.genomes.v3.1.2.hgdp_tgp.chr$i"
done
```

Just as was done for the target dataset, the variant IDs in these reference data need to be reformatted to following the format ```chr_pos_alt```. Because the gnomAD dataset is so large, the included ```recode_bim_with_snpid_hg38.R``` script can be used to run the reformatting in parallel on a server, like in the following command:
```
for i in {1..22}; do
sbatch \
--job-name=chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="Rscript /mnt/16tb/gnomad/plink/recode_bim_with_snpid_hg38.R \
-b /mnt/16tb/gnomad/plink/gnomad.genomes.v3.1.2.hgdp_tgp.chr$i\.bim \
-o /mnt/16tb/gnomad/plink/gnomad.genomes.v3.1.2.hgdp_tgp.chr$i\.bim"
done
```

...where
* ```-b``` is the path of the input bim file
* ```-o``` is the desired path of the reformated bim file that is output










