#!/usr/bin/R

if(!(require(argparse) )) install.packages("argparse")
if(!(require(tidyverse) )) install.packages("tidyverse")

library(argparse)
library(tidyverse)


# traceback to written code
options(error=function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })


# FUNCTIONS

getBinTotals <- function(m) {

    b1totals <- m %>%
        group_by(bin1) %>%
        summarize(total_count = sum(count)) %>%
        dplyr::rename(bin_idx = bin1)

    b2totals <- m %>%
        dplyr::filter(bin1 != bin2) %>%
        group_by(bin2) %>%
        summarize(total_count = sum(count)) %>%
        dplyr::rename(bin_idx = bin2)


    btotals <- rbind(b1totals, b2totals) %>%
        group_by(bin_idx) %>%
        summarize(bin_total = sum(total_count))

    return(btotals)

}

matrixbed_to_bedpe <- function(matrix_dir,
                               matrix_name,
                               bed_name,
                               bin_totals=FALSE) {

    m <- read_tsv(file.path(matrix_dir,matrix_name), col_names=c('bin1','bin2','count'))
    b <- read_tsv(file.path(matrix_dir,bed_name), col_names=c('chr','start','end','bin_idx'))

    if (!grepl('chr', b[[1,1]])) {
        b <- b %>% mutate(chr=paste0('chr',chr))
    }
    b <- b %>% mutate(chr=if_else(chr=='chrMT','chrM',chr))

    if (bin_totals) {

        btotals <- getBinTotals(m)
        b <- merge(b, btotals, all.x=TRUE) %>%
            mutate(across(bin_total, ~replace_na(.x, 0))) %>%
            arrange(bin_idx) %>%
            dplyr::select(chr,start,end,bin_total)

        bpe1 <- b[m$bin1, 1:4] %>%
            dplyr::rename(chr1=chr,
                          start1=start,
                          end1=end,
                          bin1_total=bin_total)

        bpe2 <- b[m$bin2, 1:4] %>%
            dplyr::rename(chr2=chr,
                          start2=start,
                          end2=end,
                          bin2_total=bin_total)

        bpe <- cbind(bpe1,bpe2) %>%
            mutate(count = m$count) %>%
            dplyr::filter(chr1==chr2) %>%
            dplyr::select(chr1,start1,end1,chr2,start2,end2,count,bin1_total,bin2_total)

    } else {

        bpe1 <- b[m$bin1, 1:3] %>%
            dplyr::rename(chr1=chr,
                          start1=start,
                          end1=end)

        bpe2 <- b[m$bin2, 1:3] %>%
            dplyr::rename(chr2=chr,
                          start2=start,
                          end2=end)

        bpe <- cbind(bpe1,bpe2) %>%
            mutate(count = m$count) %>%
            dplyr::filter(chr1==chr2)

    }

    return(bpe)

}

parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--bin_totals',
                    action='store_true',
                    default=FALSE,
                    help='Include total counts for each bin as 2 metadata columns in output bedpe.')
parser$add_argument('--matrix',
                    type='character',
                    default='sample_5000.matrix',
                    help='Sparse matrix of `<idx1> <idx2> <count>` HiC/HiChIP contact counts.')
parser$add_argument('--bed',
                    type='character',
                    default='sample_5000_abs.bed',
                    help='Coorsponding bed file assinging contact window indices to coordinates.')
parser$add_argument('--matrix_dir',
                    type='character',
                    default='.',
                    help='Directory with hicpro raw, sparse matrix output.')
parser$add_argument('-n', '--name',
                    type='character',
                    default='sample_5000.bedpe.gz',
                    help="Bedpe output file name")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

bin_totals <- opt$bin_totals
matrix_dir <- opt$matrix_dir
matrix_name <- opt$matrix
bed_name <- opt$bed

outname <- opt$name
outdir <- opt$outdir

bpe <- matrixbed_to_bedpe(matrix_dir, matrix_name, bed_name, bin_totals=bin_totals)
bpe %>%
    write_tsv(file.path(outdir, outname), col_names=FALSE)
