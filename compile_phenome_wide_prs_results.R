library("optparse")
library("data.table")

option_list = list(
	make_option(c("-i", "--in_dir"), type = "character", default = ".", help = "directory of the input files"),
	make_option(c("-s", "--sample"), type = "character", default = NULL, help = "the IID of the target individual"),
	make_option(c("-o", "--out_dir"), type = "character", default = ".", help = "target directory for the output file"),
	make_option(c("-p", "--out_prefix"), type = "character", default = NULL, help = "prefix for the output name")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#Open each of the files as a data frame and concatenate
cat_df <- data.frame()

for (score_file in list.files(opt$in_dir, pattern = "*.sscore")) {
	print(score_file)
	pheno <- strsplit(score_file, "[.]")
	pheno <- unlist(lapply(pheno, `[[`, 1))
	file <- read.table(paste(opt$in_dir, score_file, sep = "/"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
	df <- as.data.frame(file)
	df <- df[(df[,2] == opt$sample),]
	df$pheno <- pheno
	cat_df <- rbind(cat_df, df)
}

#Sort by percentile
cat_df <- cat_df[order(cat_df[,7], decreasing = TRUE),]
cat_df <- cat_df[,c(8, 3:7)]
names(cat_df) <- c("PHENO", "ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM", "SCORE1_AVG", "SCORE1_SUM", "SCORE_PERCENTILE")

#Write this new fam file out to a fam file
out_file_name <- paste0(opt$out_prefix, ".txt")
write.table(cat_df, file = paste(opt$out_dir, out_file_name, sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)