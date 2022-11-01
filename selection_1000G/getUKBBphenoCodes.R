#!/usr/bin/R

library(argparse)
library(tidyverse)


parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--ukbbmani_fname',
                    type='character',
                    default=NULL,
                    help='UK biobank manifest file for extracting phenotype description matching code in results file name')
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

ukbbmani_fname <- opt$ukbbmani_fname
outdir <- opt$outdir

ukbbmani <- read_tsv(ukbbmani_fname,na = c("", "NA","N/A")) %>%
    rename_with(., ~ gsub(' ', '_', tolower(.x)))

# initial script for filtering to test set phenotypes
## write each line to separate file for parallel processing in nextflow
#### NOTE: adding 'malignant' to phenosubstr adds 106 phenotypes

phenosubstr <- 'Cancer|Non-cancer|metabol|blood|haem|corpuscular|platelet|count|fraction|percentage|pulse|cholest|infect|mycoses|bacterial|viral|colitis|pneumonia|lymph|melanoma|mesothelioma|leukaemia|myeloma|enteritis|cardiac|cardio|allergy|chirrosis|crohn|carcinoma|anaemia|iron|thrombo|angio|diabet|athero|sarco|nutrition|obesity|thyro|endomet|electrolyte|fibro|migraine|headache|parkinson|sclerosis|gout|heart|arter|rheumat|hypotension|hypertension|angina|vascular|stroke|neoplasm|acute|ILD diff|ulcer|cyst'

nonphenosubstr <- 'Illnesses of|refraction|Job SOC|Job cod|predicted percentage|fraction of day|over-the-counter|encountering|cholesterol lowering|cardioplen|influencing'

ukbsub <- ukbbmani %>%
    filter(!is.na(phenotype_code),
           sex=='both_sexes',
           grepl(phenosubstr, phenotype_description, ignore.case=TRUE) | grepl("^[A-X].*", phenotype_code),
           !grepl('^S', phenotype_code),
           !grepl("^(XVIII|XIX|XX|Z).*", phenotype_code),
           !grepl(nonphenosubstr, phenotype_description, ignore.case=TRUE),
           !grepl('_raw', phenotype_code))

ukbsub %>%
   pull(phenotype_code) %>%
   writeLines(file.path(outdir, 'phenotype_codes.txt'))
