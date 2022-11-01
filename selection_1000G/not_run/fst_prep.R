#!/usr/bin/R

library(argparse)
library(plyranges)
library(tidyverse)

# traceback to written code
options(error=function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })


parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--gene_set',
                    type='character',
                    default=NULL,
                    help="Gene set to annotate eQTL with")

parser$add_argument('--snps',
                    type='character',
                    default=NULL,
                    help="SNPs to parse for input vcftools FST calc.")

parser$add_argument('-l', '--liftover',
                    type='character',
                    default=NULL,
                    help="Chain file for liftover from GTEx coordinates to hg19 for overlap")

parser$add_argument('-n', '--name',
                    type='character',
                    default='fst_snps',
                    help="Basename for output files")

parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

parser$add_argument('--script_dir',
                    type='character',
                    default=NULL,
                    help="Directory with 'trans_utils.R' file to source for functions if desired.")

opt <- parser$parse_args()

gs_fname <- opt$gene_set
snps_fname <- opt$snps
chain_fname <- opt$liftover
script_dir <- opt$script_dir
outbase <- file.path(opt$outdir, opt$name)

source(
  paste0(script_dir,
    '/trans_utils.R'
  )
)
source(
  paste0(script_dir,
    '/peak_utils.R'
  )
)

snps <- read_hg19_eqtl_granges(snps_fname, chain_fname, gs_fname)

snps %>%
  as_tibble %>%
  mutate(chrom=sub('chr','',seqnames)) %>%
  dplyr::select(chrom, start) %>%
  arrange(chrom, start) %>%
  distinct() %>%
  write_tsv(., paste0(outbase, '.positions.txt'), col_names=FALSE)
