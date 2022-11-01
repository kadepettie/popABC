#!/usr/bin/R

library(argparse)
library(tidyverse)

downloadUKBBpheno <- function(mani, pcode, sexfilt='both_sexes', download_dir='.') {

    prow <- mani %>% filter(phenotype_code==pcode, sex==sexfilt)
    pdesc <- prow$phenotype_description
    pfname <- file.path(download_dir, prow$file)
    syscomm <- prow$wget_command

    if (file.exists(pfname) & file.size(pfname) > 1000000) {
        print(paste(pdesc, 'GWAS results already downloaded here:'))
        print(paste0("    ", pfname))

    } else {

        print(paste("Downloading", pdesc, "GWAS results..."))
        pcall <- system(syscomm, ignore.stdout=TRUE, ignore.stderr=TRUE)
        if (pcall==0) print('Success!') else stop("Download failed!")

    }

}


parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--pheno_code',
                    type='character',
                    default=NULL,
                    help='Phenotype code to filter to for downloading relevant UKBB phenotype results.')
parser$add_argument('--ukbbmani_fname',
                    type='character',
                    default=NULL,
                    help='UK biobank manifest file for extracting phenotype description matching code in results file name')
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

pheno_code <- opt$pheno_code
ukbbmani_fname <- opt$ukbbmani_fname
outdir <- opt$outdir

ukbbmani <- read_tsv(ukbbmani_fname,na = c("", "NA","N/A")) %>%
    rename_with(., ~ gsub(' ', '_', tolower(.x)))

downloadUKBBpheno(ukbbmani, pheno_code, sexfilt='both_sexes', download_dir=outdir)
