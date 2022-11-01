#!/usr/bin/R

if(!(require(argparse)))  install.packages("argparse")
if(!(require(tidyverse))) install.packages("tidyverse")

library(argparse)
library(tidyverse)


# traceback to written code
options(error=function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

parser <- ArgumentParser(description="Subset a bedfile to regions linked to genes in a gene set via loops and nearest gene.")

parser$add_argument('--script_dir',
                    type='character',
                    default=NULL,
                    help="Directory with 'trans_utils.R' file to source for functions if desired.")
parser$add_argument('-i', '--inbed',
                    type='character',
                    default=NULL,
                    help="Input bed file of regions of interest.")
parser$add_argument('-n', '--name',
                    type='character',
                    default='gs_linked_regions_of_interest.bed',
                    help="Name for output bed file")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

bed_fname <- opt$inbed
script_dir <- opt$script_dir
outdir <- opt$outdir

# load common functions
if (!is.null(script_dir)) {
  source(
    paste0(script_dir,
      '/trans_utils.R'
    )
  )
}

if (grepl('mumerge', bed_fname, ignore.case=TRUE)) {
  sourcename <- 'mumerge'
  cnames = c('seqnames','start','end')
  bed <- read_tsv(bed_fname, col_names=cnames)
} else {
  # bed file should already have column names
  sourcename <- 'genrichUnion'
  bed <- read_tsv(bed_fname) %>%
    dplyr::rename(seqnames=chr)
}

gff <- bed %>%
    mutate(source=sourcename,
          feature='exon',
          start=start+1,
          score='.',
          strand='.',
          frame='.',
          attribute=paste0('gene_id=',seqnames,'_',start,'_',end)) %>%
    dplyr::select(seqnames,
                 source,
                 feature,
                 start,
                 end,
                 score,
                 strand,
                 frame,
                 attribute) %>%
    arrange(seqnames, start, end)

write_tsv(
  gff,
  paste0(
    outdir,
    '/',
    opt$name
  ),
  col_names=FALSE
)
