#!/usr/bin/R

# if(!(require(argparse)  )) install.packages("argparse")
# if(!(require(biomaRt))) {
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install("biomaRt")
# }
if(!(require(plyranges) )) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("plyranges")
}
# if(!(require(rstatix)   )) install.packages("rstatix")
if(!(require(tidyverse) )) install.packages("tidyverse")
# if(!(require(ggplot2)   )) install.packages("ggplot2")

library(argparse)
# library(biomaRt)
library(plyranges)
# library(rstatix)
library(tidyverse)
# library(ggplot2)


# traceback to written code
options(error=function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
# fix plot plotting dimensions for scaling
# options(repr.plot.width = 18, repr.plot.height = 12, repr.plot.res = 100)


# FUNCTIONS

top_peaks <- function(df,topN=150000,group_col='ancestry',count_col='count') {
    df <- df %>%
            group_by(.data[[group_col]]) %>%
            arrange(desc(.data[[group_col]]), .by_group=TRUE) %>%
            dplyr::slice(1:topN) %>%
            ungroup
}

read_group_peaks <- function(peaks_dir,
                             counts_dir,
                             g1name,
                             g2name,
                             psuff='.narrowPeak.sorted',
                             csuff='_merged.bam.Counts.bed',
                             group_type='ancestry',
                             nstrongest=NULL) {

    pcnames <- c('seqnames',
             'start',
             'end',
             'name',
             'score',
             'strand',
             'signalValue',
             'neglog10pval',
             'neglog10qval',
             'summit_pos')
    ccnames <- c('seqnames','start','end','count')

    g1 <- merge(
        read_tsv(file.path(peaks_dir,paste0(g1name,psuff)),
                 col_names=pcnames),
        read_tsv(file.path(counts_dir,paste0(g1name,psuff,'.',g1name,csuff)),
                 col_names=ccnames)
    ) %>%
        mutate(!!group_type := g1name)

    g2 <- merge(
        read_tsv(file.path(peaks_dir,paste0(g2name,psuff)),
                 col_names=pcnames),
        read_tsv(file.path(counts_dir,paste0(g2name,psuff,'.',g2name,csuff)),
                 col_names=ccnames)
    ) %>%
        mutate(!!group_type := g2name)

    gm <- rbind(g1,g2)

    if (!is.null(nstrongest)) {
        gm <- top_peaks(gm,
                        topN=nstrongest,
                        group_col=group_type,
                        count_col='count')
    }

    return(gm)

}

extend_summits <- function(summits, gbuild='hg19', peak_extend=250) {

  ggbuild <- genome_info(gbuild) %>%
     dplyr::filter(seqnames %in% (summits %>% pull(seqnames) %>% unique))
  gseqnames <- ggbuild %>%
     as_tibble() %>%
     pull(seqnames) %>%
     as.character()
  gseqlengths <- ggbuild %>%
      as_tibble() %>%
      pull(end)

  tsumgr <- summits %>%
    mutate(strand = '*') %>%
    as_granges()

  seqlevelorder <- intersect(seqlevels(ggbuild), seqlevels(tsumgr))
  seqlevels(tsumgr) <- seqlevelorder
  tsumgr <- set_genome_info(tsumgr,
                    genome=gbuild,
                    seqnames=gseqnames,
                    seqlengths=gseqlengths,
                    is_circular=rep(FALSE,length(gseqlengths))) %>%
    trim()

  tcan <- tsumgr %>%
    anchor_start() %>%
    mutate(start = start + summit_pos,
           width = 1) %>%
    anchor_center() %>%
    mutate(width=peak_extend*2+1) %>%
    trim()

  return(tcan)

}

merge_interleave_group_candidates <- function(tcan,
                                              g1name,
                                              g2name,
                                              group_col='ancestry',
                                              rank_col='count',
                                              n_enhancers=150000) {
    # granges doesn't handle dynamic column names like tidyverse does
    # need to convert to tibble first, rename columns then back to granges
    # to merge overlapping ranges by group
    # don't need to retain genome info since we are done with range arithmatic
    
    print('DEBUG')
    print(colnames(tcan %>% as_tibble))
    print(dim(tcan %>% as_tibble()))	
    tcan_groupmerge_pre <- tcan %>%
        as_tibble() %>%
        dplyr::rename(group := {{group_col}},
                      rank_var := {{rank_col}}) %>%
        as_granges()

    print(colnames(tcan_groupmerge_pre %>% as_tibble()))
    tcan_groupmerge <- tcan_groupmerge_pre %>% 
        group_by(group) %>%
        reduce_ranges(rank_var = max(rank_var))

    print('Post group-merge')	
    print(dim(tcan_groupmerge))
    tcan_allmerge <- tcan_groupmerge %>% 
        reduce_ranges(rankvar_g1 = max(rank_var[group==g1name]),
                      rankvar_g2 = max(rank_var[group==g2name]),
                      group = paste(sort(group),collapse=',')) %>%
        as_tibble() %>%
        dplyr::rename("{rank_col}_{g1name}" := rankvar_g1,
                      "{rank_col}_{g2name}" := rankvar_g2,
                      !!group_col := group)

    rankcolg1 <- paste(rank_col,g1name,sep='_')
    rankcolg2 <- paste(rank_col,g2name,sep='_')

    tcan_g1 <- tcan_allmerge %>%
        arrange(desc(.data[[rankcolg1]])) %>%
        dplyr::filter(across(all_of(rankcolg1)) != -Inf,
                      nchar(as.character(seqnames)) < 6)

    tcan_g2 <- tcan_allmerge %>%
        arrange(desc(.data[[rankcolg2]])) %>%
        dplyr::filter(across(all_of(rankcolg2)) != -Inf,
                      nchar(as.character(seqnames)) < 6)

    interleaved_ranks <- order(c(1:nrow(tcan_g1), 1:nrow(tcan_g2)))

    tcan_final <- rbind(tcan_g1,tcan_g2)[interleaved_ranks,] %>%
        distinct() %>%
        dplyr::slice(1:n_enhancers)

    return(tcan_final)
}

parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--genome_build',
                    type='character',
                    default='hg19',
                    help='For getting seqlengths via plyranges.')
parser$add_argument('--peaks_suffix',
                    type='character',
                    default='.narrowPeak.sorted')
parser$add_argument('--counts_suffix',
                    type='character',
                    default='_merged.bam.Counts.bed')
parser$add_argument('--group1_name',
                    type='character',
                    default='AFR')
parser$add_argument('--group2_name',
                    type='character',
                    default='EUR')
parser$add_argument('--group_type',
                    type='character',
                    default='ancestry')
parser$add_argument('--peaks_dir',
                    type='character',
                    default='.')
parser$add_argument('--counts_dir',
                    type='character',
                    default='.')
parser$add_argument('--peakExtendFromSummit',
                    type='integer',
                    default=NULL,
                    help="Number of bp to extend in either direction from peak summits.")
parser$add_argument('--nStrongestPeaks',
                    type='integer',
                    default=150000,
                    help="Number of top peaks to output after interleaving.")
parser$add_argument('-n', '--name',
                    type='character',
                    default='combinedTopCcres',
                    help="Basename for output files")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

counts_dir <- opt$counts_dir
peaks_dir <- opt$peaks_dir
gtype <- opt$group_type
g1name <- opt$group1_name
g2name <- opt$group2_name
psuff <- opt$peaks_suffix
csuff <- opt$counts_suffix
peak_extend <- opt$peakExtendFromSummit
gbuild <- opt$genome_build
n_enhancers <- opt$nStrongestPeaks
outbase <- opt$name
outdir <- opt$outdir

gsummits <- read_group_peaks(peaks_dir,
                             counts_dir,
                             g1name,
                             g2name,
                             psuff=psuff,
                             csuff=csuff,
                             group_type=gtype,
                             nstrongest=NULL)

tsummits <- top_peaks(gsummits,
                     topN=n_enhancers,
                     group_col=gtype,
                     count_col='count')

top_cands <- extend_summits(tsummits, gbuild=gbuild, peak_extend=peak_extend)

comb_top_cands <- merge_interleave_group_candidates(top_cands,
                                                    g1name,
                                                    g2name,
                                                    group_col=gtype,
                                                    rank_col='count',
                                                    n_enhancers=n_enhancers)

comb_top_cands %>%
    write_tsv(file.path(outdir,paste0(outbase,".Counts.txt")))

comb_top_cands %>%
    dplyr::select(1:3) %>%
    write_tsv(file.path(outdir,paste0(outbase,".bed")), col_names=FALSE)
