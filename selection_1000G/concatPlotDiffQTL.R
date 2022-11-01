#!/usr/bin/R

library(argparse)
library(tidyverse)
library(ggplot2)

concatForPlotting <- function(fdir,
                             f_patt='fisherEnrichments.txt$') {

    fnames <- list.files(fdir,
                         pattern=f_patt,
                         full.names=FALSE)

    adf <- tibble()
    for (f in fnames) {

        print(paste0('Reading ',f))

        df <- read_tsv(file.path(fdir, f))

        adf <- rbind(adf, df)

    }

    return(adf)

}



parser <- ArgumentParser(description='Aggregate differential test results of ABC Score components for QC/diagnostics and spot checking.')

parser$add_argument('-n', '--name',
                    type='character',
                    default='diffABC_DE_overlap',
                    help="Basename for output files")
parser$add_argument('--plotdir',
                    type='character',
                    default=NULL, # /home/kpettie/code/github/plotting
                    help="Directory with plotting function in 'plotting.R'. Required to make plots")
parser$add_argument('--wilcox_pattern',
                    type='character',
                    default="fstWilcoxData.txt.gz",
                    help="Pattern uniquely identifying dataframes to aggregate from agg_dir for Wilcoxon FST tests (glob with partial matching allowed, not regex, to avoid escape characters on command line).")
parser$add_argument('--fisher_pattern',
                    type='character',
                    default="fisherEnrichments.txt",
                    help="Pattern uniquely identifying fisher enrichments to aggregate from agg_dir (glob with partial matching allowed, not regex, to avoid escape characters on command line).")
parser$add_argument('--agg_dir',
                    type='character',
                    default='.',
                    help="Directory with files from different QTL overlap analyses for aggregating")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

plotdir <- opt$plotdir
w_patt <- opt$wilcox_pattern
f_patt <- opt$fisher_pattern
a_dir <- opt$agg_dir
outdir <- opt$outdir
outbase <- opt$name

# fisher enrichments and binomSignTest
bfe <- concatForPlotting(a_dir, f_patt=f_patt) %>%
    mutate(fst_percentile = as_factor(fst_percentile)) %>%
    distinct()
write_tsv(bfe, file.path(outdir, paste0(outbase, ".diffQTLoverlap.fisherEnrichmentsAgg.txt")))

# Wilcoxon test data
bfw <- concatForPlotting(a_dir, f_patt=w_patt) %>%
    distinct()
write_tsv(bfw, file.path(outdir, paste0(outbase, ".diffQTLoverlap.fstWilcoxDataAgg.txt.gz")))

if (!is.null(plotdir)) {

    source(file.path(plotdir,"plotting.R"))

    ####### Fisher enrichment PLOTS #######

    # bQTL
    # overall
    p1 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(!is.na(chip_type),
                                        test_type=='diffQTL',
                                        cre_type=='all'),
                      yvar_col='chip_type',
                      ylabel='QTL type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar=NULL,
                      groupvarlab=NULL,
                      groupvarord=NULL,
                      colorvals=NULL,
                      sep_groups=FALSE,
                      w=21,
                      h=7,
                      facets='score_type',
                      frows=1,
                      fscales='fixed', # free, free_x, free_y
                      fdir='h',
                      legendpos=NULL, # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Diff-CRE enrichment for QTL") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p1
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.bQTLoverall.png")), width=21, height=7)

    p2 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(!is.na(chip_type),
                                        test_type=='diffQTL',
                                        cre_type!='all'),
                      yvar_col='chip_type',
                      ylabel='QTL type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='cre_type',
                      groupvarord=NULL,
                      colorvals=c('pink3','darkgreen'),
                      sep_groups=FALSE,
                      w=21,
                      h=7,
                      facets='score_type',
                      frows=1,
                      fscales='fixed', # free, free_x, free_y
                      fdir='h',
                      legendpos=c(.3,.2), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Diff-CRE enrichment for QTL") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p2
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.bQTLoverall.byCreType.png")), width=21, height=7)

    # direction
    p3 <- plotFisherEnrichments(bfe %>%
                                    dplyr::filter(!is.na(chip_type),
                                                  test_type=='diffQTLdirection',
                                                  cre_type=='all'),
                                yvar_col='chip_type',
                                ylabel='QTL type',
                                xvar_col='odds_ratio',
                                althypoth='greater',
                                line_intercept=1,
                                groupvar='fst_percentile',
                                groupvarord=NULL,
                                colorvals=c('black','purple2'),
                                sep_groups=FALSE,
                                w=21,
                                h=7,
                                facets='score_type',
                                frows=1,
                                fscales='fixed', # free, free_x, free_y
                                fdir='v',
                                legendpos=c(.3,.2), # c(xpos, ypos)
                                debug=FALSE) +
                  ggtitle("Diff-CRE enrichment for matching QTL direction") +
                  theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p3
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.bQTLdirection.png")), width=21, height=7)

    p4 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(!is.na(chip_type),
                                        test_type=='diffQTLdirection',
                                        cre_type!='all',
                                        fst_percentile=='0'),
                      yvar_col='chip_type',
                      ylabel='QTL type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='cre_type',
                      groupvarord=NULL,
                      colorvals=c('pink3','darkgreen'),
                      sep_groups=FALSE,
                      w=21,
                      h=7,
                      facets='score_type',
                      frows=1,
                      fscales='fixed', # free, free_x, free_y
                      fdir='h',
                      legendpos=c(.3,.2), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Diff-CRE enrichment for matching QTL direction") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p4
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.bQTLdirection.byCreType.png")), width=21, height=7)

    p5 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(!is.na(chip_type),
                                        test_type=='diffQTLdirection',
                                        cre_type!='all'),
                      yvar_col='chip_type',
                      ylabel='QTL type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='fst_percentile',
                      groupvarord=NULL,
                      colorvals=c('black','purple2'),
                      sep_groups=FALSE,
                      w=21,
                      h=14,
                      facets=c('score_type','cre_type'),
                      frows=2,
                      fscales='fixed', # free, free_x, free_y
                      fdir='v',
                      legendpos=c(.3,.2), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Diff-CRE enrichment for matching QTL direction") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p5
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.bQTLdirection.byCreTypeFST.png")), width=21, height=14)

    # DE direction
    # LCL
    p6 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(!is.na(chip_type),
                                        test_type=='DEdiffQTLdirection',
                                        is.na(self_promoter),
                                        de_celltype=='lcl'),
                      yvar_col='chip_type',
                      ylabel='QTL type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='fst_percentile',
                      groupvarord=NULL,
                      colorvals=c('black','purple2'),
                      sep_groups=FALSE,
                      w=21,
                      h=7,
                      facets='score_type',
                      frows=1,
                      fscales='fixed', # free, free_x, free_y
                      fdir='h',
                      legendpos=c(.3,.2), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Matching diff-CRE/QTL enrichment for matching LCL DE direction") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p6
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.bQTLdirectionDE.LCL.png")), width=21, height=7)

    p7 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(!is.na(chip_type),
                                        test_type=='DEdiffQTLdirection',
                                        !is.na(self_promoter),
                                        de_celltype=='lcl') %>%
                          mutate(self_promoter = if_else(self_promoter, 'self-promoter', 'non-self-promoter')),
                      yvar_col='chip_type',
                      ylabel='QTL type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='fst_percentile',
                      groupvarord=NULL,
                      colorvals=c('black','purple2'),
                      sep_groups=FALSE,
                      w=21,
                      h=14,
                      facets=c('score_type','self_promoter'),
                      frows=2,
                      fscales='fixed', # free, free_x, free_y
                      fdir='v',
                      legendpos=c(.2,.8), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Matching diff-CRE/QTL enrichment for matching LCL DE direction") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p7
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.bQTLdirectionDE.bySelfPromoter.LCL.png")), width=21, height=14)

    # PBMC
    p8 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(!is.na(chip_type),
                                        test_type=='DEdiffQTLdirection',
                                        is.na(self_promoter),
                                        de_celltype=='pbmc'),
                      yvar_col='chip_type',
                      ylabel='QTL type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='fst_percentile',
                      groupvarord=NULL,
                      colorvals=c('black','purple2'),
                      sep_groups=FALSE,
                      w=21,
                      h=7,
                      facets='score_type',
                      frows=1,
                      fscales='fixed', # free, free_x, free_y
                      fdir='h',
                      legendpos=c(.3,.2), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Matching diff-CRE/QTL enrichment for matching PBMC DE direction") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p8
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.bQTLdirectionDE.PBMC.png")), width=21, height=7)

    p9 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(!is.na(chip_type),
                                        test_type=='DEdiffQTLdirection',
                                        !is.na(self_promoter),
                                        de_celltype=='pbmc') %>%
                          mutate(self_promoter = if_else(self_promoter, 'self-promoter', 'non-self-promoter')),
                      yvar_col='chip_type',
                      ylabel='QTL type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='fst_percentile',
                      groupvarord=NULL,
                      colorvals=c('black','purple2'),
                      sep_groups=FALSE,
                      w=21,
                      h=14,
                      facets=c('score_type','self_promoter'),
                      frows=2,
                      fscales='fixed', # free, free_x, free_y
                      fdir='v',
                      legendpos=c(.3,.2), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Matching diff-CRE/QTL enrichment for matching PBMC DE direction") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p9
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.bQTLdirectionDE.bySelfPromoter.PBMC.png")), width=21, height=14)

    # eQTL
    # direction
    p10 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(qtl_type=='eqtl_lcl',
                                        test_type=='diffQTLdirection'),
                      yvar_col='cre_type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='fst_percentile',
                      groupvarord=NULL,
                      colorvals=c('black','purple2'),
                      sep_groups=FALSE,
                      w=21,
                      h=7,
                      facets='score_type',
                      frows=1,
                      fscales='fixed', # free, free_x, free_y
                      fdir='h',
                      legendpos=c(.3,.2), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Diff-CRE enrichment for matching LCL eQTL direction") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p10
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.eQTLdirection.LCL.png")), width=21, height=7)

    p11 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(qtl_type=='eqtl_pbmc',
                                        test_type=='diffQTLdirection'),
                      yvar_col='cre_type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='fst_percentile',
                      groupvarord=NULL,
                      colorvals=c('black','purple2'),
                      sep_groups=FALSE,
                      w=21,
                      h=7,
                      facets='score_type',
                      frows=1,
                      fscales='fixed', # free, free_x, free_y
                      fdir='h',
                      legendpos=c(.3,.2), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Diff-CRE enrichment for matching PBMC eQTL direction") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p11
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.eQTLdirection.PBMC.png")), width=21, height=7)


    # DE direction
    p12 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(qtl_type=='eqtl_lcl',
                                        test_type=='DEdiffQTLdirection') %>%
                          mutate(self_promoter=case_when(is.na(self_promoter) ~ 'all',
                                                         self_promoter ~ 'self-promoter',
                                                         !self_promoter ~ 'non-self-promoter')),
                      yvar_col='self_promoter',
                      ylabel='cre_type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='fst_percentile',
                      groupvarord=NULL,
                      colorvals=c('black','purple2'),
                      sep_groups=FALSE,
                      w=21,
                      h=7,
                      facets='score_type',
                      frows=1,
                      fscales='fixed', # free, free_x, free_y
                      fdir='h',
                      legendpos=c(.9,.2), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Matching diff-CRE/LCL eQTL enrichment for matching LCL DE direction") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p12
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.eQTLdirectionDE.bySelfPromoter.LCL.png")), width=21, height=7)

    p13 <- plotFisherEnrichments(bfe %>%
                          dplyr::filter(qtl_type=='eqtl_pbmc',
                                        test_type=='DEdiffQTLdirection') %>%
                          mutate(self_promoter=case_when(is.na(self_promoter) ~ 'all',
                                                         self_promoter ~ 'self-promoter',
                                                         !self_promoter ~ 'non-self-promoter')),
                      yvar_col='self_promoter',
                      ylabel='cre_type',
                      xvar_col='odds_ratio',
                      althypoth='greater',
                      line_intercept=1,
                      groupvar='fst_percentile',
                      groupvarord=NULL,
                      colorvals=c('black','purple2'),
                      sep_groups=FALSE,
                      w=21,
                      h=7,
                      facets='score_type',
                      frows=1,
                      fscales='fixed', # free, free_x, free_y
                      fdir='h',
                      legendpos=c(.9,.2), # c(xpos, ypos)
                      debug=FALSE) +
        ggtitle("Matching diff-CRE/PBMC eQTL enrichment for matching PBMC DE direction") +
        theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p13
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.eQTLdirectionDE.bySelfPromoter.PBMC.png")), width=21, height=7)


    ####### Binomial sign test PLOTS #######
    print('Making binomial sign test plots...')

    # bQTL

    p14 <- plotFisherEnrichments(bfe %>%
                                    dplyr::filter(qtl_type=='bqtl',
                                                  test_type=='binomSignTest',
                                                  n_trials > 0,
                                                  cre_type=='all'),
                                yvar_col='chip_type',
                                ylabel='QTL type',
                                xvar_col='odds_ratio',
                                althypoth='two.sided',
                                line_intercept=1,
                                groupvar='fst_percentile',
                                groupvarord=NULL,
                                colorvals=c('black','purple2'),
                                sep_groups=FALSE,
                                w=21,
                                h=7,
                                facets='score_type',
                                frows=1,
                                fscales='fixed', # free, free_x, free_y
                                fdir='h',
                                legendpos=c(.3,.2), # c(xpos, ypos)
                                debug=FALSE) +
                  ggtitle("Matching diff-CRE/QTL binomial sign test") +
                  theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p14
    ggsave(file.path(outdir, paste0(outbase, ".binomSignTest.bQTL.png")), width=21, height=7)

    p15 <- plotFisherEnrichments(bfe %>%
                                    dplyr::filter(qtl_type=='bqtl',
                                                  test_type=='binomSignTest',
                                                  n_trials > 0,
                                                  cre_type!='all'),
                                yvar_col='chip_type',
                                ylabel='QTL type',
                                xvar_col='odds_ratio',
                                althypoth='two.sided',
                                line_intercept=1,
                                groupvar='fst_percentile',
                                groupvarord=NULL,
                                colorvals=c('black','purple2'),
                                sep_groups=FALSE,
                                w=21,
                                h=14,
                                facets=c('score_type','cre_type'),
                                frows=2,
                                fscales='fixed', # free, free_x, free_y
                                fdir='v',
                                legendpos=c(.3,.2), # c(xpos, ypos)
                                debug=FALSE) +
                  ggtitle("Matching diff-CRE/QTL binomial sign test") +
                  theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p15
    ggsave(file.path(outdir, paste0(outbase, ".binomSignTest.bQTL.byCreType.png")), width=21, height=14)


    # eQTL
    # LCL
    p16 <- plotFisherEnrichments(bfe %>%
                                    dplyr::filter(qtl_type=='eqtl_lcl',
                                                  test_type=='binomSignTest',
                                                  n_trials > 0),
                                yvar_col='cre_type',
                                ylabel='CRE type',
                                xvar_col='odds_ratio',
                                althypoth='two.sided',
                                line_intercept=1,
                                groupvar='fst_percentile',
                                groupvarord=NULL,
                                colorvals=c('black','purple2'),
                                sep_groups=FALSE,
                                w=21,
                                h=7,
                                facets='score_type',
                                frows=1,
                                fscales='fixed', # free, free_x, free_y
                                fdir='h',
                                legendpos=c(.3,.2), # c(xpos, ypos)
                                debug=FALSE) +
                  ggtitle("Matching diff-CRE/LCL eQTL binomial sign test") +
                  theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p16
    ggsave(file.path(outdir, paste0(outbase, ".binomSignTest.eQTL.LCL.png")), width=21, height=7)

    p17 <- plotFisherEnrichments(bfe %>%
                                    dplyr::filter(qtl_type=='eqtl_pbmc',
                                                  test_type=='binomSignTest',
                                                  n_trials > 0),
                                yvar_col='cre_type',
                                ylabel='CRE type',
                                xvar_col='odds_ratio',
                                althypoth='two.sided',
                                line_intercept=1,
                                groupvar='fst_percentile',
                                groupvarord=NULL,
                                colorvals=c('black','purple2'),
                                sep_groups=FALSE,
                                w=21,
                                h=7,
                                facets='score_type',
                                frows=1,
                                fscales='fixed', # free, free_x, free_y
                                fdir='h',
                                legendpos=c(.3,.2), # c(xpos, ypos)
                                debug=FALSE) +
                  ggtitle("Matching diff-CRE/PBMC eQTL binomial sign test") +
                  theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))
    p17
    ggsave(file.path(outdir, paste0(outbase, ".binomSignTest.eQTL.PBMC.png")), width=21, height=7)

    ####### Wilcoxon test PLOTS #######

    p18 <- groupedBoxplot(bfw %>%
                             mutate(diff_status = fct_relevel(diff_status,c('non-diff','diff'))) %>%
                             dplyr::filter(qtl_type=='bqtl'),
             'chip_type',
             'fst_var',
             groupvar='diff_status',
             xlab='ChIP type',
             ylab='QTL FST',
             plotpoints=TRUE,
             show_wilcox=TRUE,
             alt='less',
             logscale=FALSE,
             colorX=FALSE,
             colorvals=c('darkgray','red3'),
             showleg=TRUE,
             legendpos='top',
             horizontal=TRUE,
             facets=c('cre_type','score_type'),
             frows=2,
             fscales='fixed',
             fdir='h',
             w=15,h=10,
             angleHjVj=c(0,.5,0),
             debug=FALSE) +
          theme(legend.title = element_blank())
    p18
    ggsave(file.path(outdir, paste0(outbase, ".fstWilcox.bQTL.png")), width=15, height=10)

    p19 <- groupedBoxplot(bfw %>%
                    dplyr::filter(qtl_type=='bqtl',
                                  diff_status=='diff'),
             'chip_type',
             'fst_var',
             groupvar='cre_type',
             xlab='ChIP type',
             ylab='QTL FST',
             plotpoints=TRUE,
             show_wilcox=TRUE,
             alt='greater',
             logscale=FALSE,
             colorX=FALSE,
             colorvals=c('pink3','darkgreen'),
             showleg=TRUE,
             legendpos='top',
             horizontal=TRUE,
             facets='score_type',
             frows=1,
             fscales='fixed',
             fdir='h',
             w=15,h=7,
             angleHjVj=c(0,.5,0),
             debug=FALSE) +
          labs(fill='Diff-CRE type', alpha='Diff-CRE type')
    p19
    ggsave(file.path(outdir, paste0(outbase, ".fstWilcox.bQTLdiff.png")), width=15, height=7)

    p20 <- groupedBoxplot(bfw %>%
                             mutate(diff_status = fct_relevel(diff_status,c('non-diff','diff'))) %>%
                             dplyr::filter(qtl_type!='bqtl'),
             'score_type',
             'fst_var',
             groupvar='diff_status',
             xlab='Score type',
             ylab='QTL FST',
             plotpoints=TRUE,
             show_wilcox=TRUE,
             alt='less',
             logscale=FALSE,
             colorX=FALSE,
             colorvals=c('darkgray','red3'),
             showleg=TRUE,
             legendpos='top',
             horizontal=TRUE,
             facets=c('cre_type','qtl_type'),
             frows=2,
             fscales='free_x',
             fdir='h',
             w=7,h=7,
             angleHjVj=c(0,.5,0),
             debug=FALSE)  +
          theme(legend.title = element_blank())
    p20
    ggsave(file.path(outdir, paste0(outbase, ".fstWilcox.eQTL.png")), width=7, height=7)


}
