#!/usr/bin/R

if(!(require(argparse)))  install.packages("argparse")
if(!(require(tidyverse))) install.packages("tidyverse")

library(argparse)
library(tidyverse)

# traceback to written code
options(error=function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
# fix plot plotting dimensions for scaling
options(repr.plot.width = 18, repr.plot.height = 12, repr.plot.res = 100)


parser <- ArgumentParser(description='Run g:Profiler on dataframe with given parameters.')

parser$add_argument('--custom_bg',
                    type='character',
                    default='yes',
                    help="Significant genes from both ancestries used as background ('yes', 'no').")
parser$add_argument('--ordered',
                    type='character',
                    default='yes',
                    help="Forground genes interpreted as ordered by descending significance ('yes', 'no').")
parser$add_argument('--gene_link',
                    type='character',
                    default='loopednearest',
                    help="Method to link peaks to genes ('eGene','loopednearest','any').")
parser$add_argument('--pthresh',
                    type='double',
                    default=0.005,
                    help="DA p-value threshold for defining foreground (and custom background).")
parser$add_argument('--fg_anc',
                    type='character',
                    default='EUR',
                    help="High CA ancestry to use in g:Profiler foreground.")
parser$add_argument('--script_dir',
                    type='character',
                    default=NULL,
                    help="Directory with 'trans_utils.R' file to source for functions if desired.")
parser$add_argument('--peak_da',
                    type='character',
                    default=NULL,
                    help="File with all tested peaks and their diffATAC p-values.")
parser$add_argument('-n', '--name',
                    type='character',
                    default='peak_snp_aggregate',
                    help="Basename for output files")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

bglog <- TRUE
if (opt$custom_bg=='no') bglog <- FALSE
olog <- TRUE
if (opt$ordered=='no') olog <- FALSE
pv <- opt$pthresh
glink <- opt$gene_link
highca <- opt$fg_anc
peaks <- opt$peak_da
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

print(paste0("Loading ", peaks, " ..."))
peak_df <- read_tsv(peaks)

print("Running g:Profiler...")
gost_df <- gprofiler_diffatac(peak_df,
                              highCA=highca,
                              gene_link=glink,
                              pthresh=pv,
                              ordered=olog,
                              custom_bg=bglog,
                              debug=FALSE,
                              outroot=outbase)

if (is.null(gost_df)) {
	print('writing no results')
	write('no results', 'nores_gprofiler.txt')
	write('no results', 'nores_gprofiler_summary.txt')
} else {
	print('g:Profiler success')
}
