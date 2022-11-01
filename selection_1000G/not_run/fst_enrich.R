#!/usr/bin/R

library(argparse)
library(plyranges)
library(rstatix)
library(tidyverse)
library(ggplot2)

# traceback to written code
options(error=function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

bname <- function(f,pathstrip=TRUE) {
  b <- unlist(strsplit(f,'\\.'))
  b <- b[length(b)-1] # grab element preceding extension
  if (length(unlist(strsplit(b,'\\.')))==1) {
    if (pathstrip) {
      b <- unlist(strsplit(b,'/'))
      return(b[length(b)])
    } else {
      return(b)
    }
  } else {
    return(bname(b,pathstrip=pathstrip))
  }
}

peak_fst_density <- function(df, outdir='.', name='fst') {

  cdf <- df %>%
             group_by(peak_type) %>%
             summarize(mean = mean(FST), count=length(FST))
  print(cdf)
  write_tsv(cdf,
            file.path(opt$outdir,
                      paste0(opt$name,
                             ".summary.txt")))

  ggplot(df, aes(x=FST, fill=peak_type, colour=peak_type)) +
    geom_density(alpha=0.2) +
    geom_vline(data=cdf, aes(xintercept=mean,  colour=peak_type),
               linetype="dashed", size=0.5) +
    theme_classic(15) + theme(legend.position=c(0.8,0.5))

  ggsave(paste0(outdir,"/",name,".density.png"))

}

peak_fst_violin <- function(df, outdir='.', name='fst') {

  ggplot(df, aes(peak_type, FST)) +
    geom_jitter(height = 0, width = 0.1, alpha=0.5) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), aes(fill = peak_type), alpha=0.2) +
    theme_classic(15) +
    guides(fill=FALSE) +
    labs(x = "Peak type",
         y = "FST")

  ggsave(paste0(outdir,"/",name,".violin.png"))

}

parser <- ArgumentParser(description='Test for enrichment of high Fst in SNPs for which Fst has been calculated using vcftools --weir-fst-pop')
parser$add_argument('-b','--background',
                    type='character',
                    default=NULL,
                    help="File with Fst calculations for backgrounds SNPs (e.g., EUR_vs_AFR_welchs_pvalue_gt0.5.weir.fst)")
parser$add_argument('-f', '--foreground',
                    type='character',
                    default=NULL,
                    help="File with Fst calculations for foreground SNPs")
parser$add_argument('-p', '--peak_genes',
                    type='character',
                    default=NULL,
                    help="File with differential peaks corresponding to the provided foreground SNPs linked to candidate target genes.")
parser$add_argument('-d', '--peak_densities',
                    type='character',
                    default=NULL,
                    help="File with all tested peaks and their densities in each population.")
parser$add_argument('-g', '--gene_set',
                    type='character',
                    default=NULL,
                    help="List of genes to use for additional foreground test.")
parser$add_argument('-n', '--name',
                    type='character',
                    default='fst_enrichment',
                    help="Basename for output files")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

cnames <- c('seqnames','start','end','gene_chr','gene_start','gene_end','ensembl_gene_id', 'gene_strand', 'distance')
ctypes <- 'ciiciic_ci'
gp <- read_tsv(opt$peak_genes, col_names=cnames, col_types=ctypes, quote="") %>%
        as_granges()

gs <- readLines(opt$gene_set)
gs_name <- bname(opt$gene_set)

bg <- read_tsv(opt$background, quote='', na = c("", "NA", "NaN", "-nan")) %>%
        drop_na() %>%
        rename(FST = WEIR_AND_COCKERHAM_FST) %>%
        mutate(FST = if_else(FST < 0, 0, FST))
fg <- read_tsv(opt$foreground, quote='', na = c("", "NA", "NaN", "-nan")) %>%
        drop_na() %>%
        rename(FST = WEIR_AND_COCKERHAM_FST) %>%
        mutate(FST = if_else(FST < 0, 0, FST))

pd <- read_tsv(opt$peak_densities, quote="") %>%
        mutate(peak_id=paste(chr,start,end,sep="_")) %>%
        rename(seqnames=chr) %>%
        as_granges()

fst <- as_tibble(merge(fg, bg,
                       by=c('CHROM','POS'),
                       all=TRUE,
                       suffixes=c('_diffATAC','_non-diffATAC'))) %>%
            pivot_longer(cols=starts_with('FST'),
                        names_to='peak_type',
                        names_prefix='FST_',
                        values_to='FST',
                        values_drop_na=TRUE) %>%

            mutate(CHROM=paste0('chr',CHROM)) %>%
            rename(seqnames=CHROM, start=POS) %>%
            as_granges(width=1) %>%

            # overlap with all peaks
            join_overlap_left(., pd) %>%
            as_tibble() %>%
            arrange(peak_id, desc(FST)) %>%
            # take only the SNP with highest Fst per peak
            filter(peak_id != lag(peak_id, default=1)) %>%
            as_granges() %>%

            # overlap with differntial peaks linked to genes
            join_overlap_left(., gp) %>%

            as_tibble() %>%
            # identify diffATAC peaks linked to genes in gene set
            mutate(snp_id=paste(seqnames,start,end,sep='_'),
                   peak_type=if_else(ensembl_gene_id %in% gs,
                                     paste0('diffATAC_',gs_name),
                                     peak_type)) %>%
            mutate(peak_type=fct_relevel(peak_type,
                                         c(paste0('diffATAC_',gs_name),
                                           'diffATAC',
                                           'non-diffATAC'))) %>%
            arrange(snp_id, peak_type) %>%
            # remove duplicate SNPs linked to multiple genes
            # prioritizing those linked to gene set of interest
            filter(snp_id != lag(snp_id, default=1))


write_tsv(fst,
          file.path(opt$outdir,
                    paste0(opt$name,
                           ".txt")))

wfst <- fst %>%
          wilcox_test(FST ~ peak_type, p.adjust.method="none")
print(wfst)
write_tsv(wfst,
          file.path(opt$outdir,
                    paste0(opt$name,
                           ".wilcox.txt")))

peak_fst_density(fst, outdir=opt$outdir, name=opt$name)
peak_fst_violin(fst, outdir=opt$outdir, name=opt$name)
