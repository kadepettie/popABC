#!/usr/bin/R

# if(!(require(argparse)  )) install.packages("argparse")
# if(!(require(biomaRt))) {
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install("biomaRt")
# }

if(!(require(argparse) )) install.packages("argparse")
if(!(require(tidyverse) )) install.packages("tidyverse")
if(!(require(rstatix)   )) install.packages("rstatix")
if(!(require(ggplot2)   )) install.packages("ggplot2")
if(!(require(ggrepel)   )) install.packages("ggrepel")

library(argparse)
library(tidyverse)
library(rstatix)
# library(biomaRt)
library(ggplot2)
library(ggrepel)


# traceback to written code
options(error=function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
# fix plot plotting dimensions for scaling
# options(repr.plot.width = 18, repr.plot.height = 12, repr.plot.res = 100)


# FUNCTIONS

positive_ccres_allsamps <- function(pred_dir,
                                    qnorm_suffix,
                                    ethreshold=0.015,
                                    pthreshold=0.1,
                                    score_column='powerlaw.Score',
                                    cnames=TRUE) {

    #' Get cCRE-gene pair predictions that meet the give thresholds in at least one sample
    #'
    #' * `ethreshold` applies to non-promoters
    #' * `pthreshold` applies to promoters and should be more stringent to
    #'    minimize false positives due to transcriptional interference,
    #'    trans effects, and/or promoter competition

    predfnames <- list.files(pred_dir,
                        pattern='EnhancerPredictionsAllPutative.txt.gz',
                        recursive=TRUE,
                        full.names=TRUE) %>%
        grep(pattern=paste0('[A-Z]{3}_rep[1-2]\\.',qnorm_suffix), value=TRUE)

    ccremax <- tibble()
    for (predfname in predfnames) {

        cccremax <- read_tsv(predfname, col_names=cnames) %>%
            dplyr::filter(case_when(class=='promoter' ~ across(all_of(score_column)) >= pthreshold,
                                    TRUE ~ across(all_of(score_column)) >= ethreshold)) %>%
            dplyr::select(chr,start,end,class,TargetGene,TargetGeneTSS,TargetGeneExpression,distance,isSelfPromoter)

        prev_Nccres <- nrow(ccremax)
        ccremax <- rbind(ccremax, cccremax) %>% distinct()
        curr_Nccres <- nrow(ccremax)

        print(paste('Done with', predfname))
        print(paste0('    ',curr_Nccres-prev_Nccres,' new cCREs added'))

    }

    return(ccremax)
}

concat_positive_predictions <- function(pred_dir,
                                        qnorm_suffix,
                                        positive_preds,
                                        cnames=TRUE) {

    #' Get scores in all samples for cCRE-gene pair predictions made in at least one sample

    predfnames <- list.files(pred_dir,
                        pattern='EnhancerPredictionsAllPutative.txt.gz',
                        recursive=TRUE,
                        full.names=TRUE) %>%
        grep(pattern=paste0('[A-Z]{3}_rep[1-2]\\.',qnorm_suffix), value=TRUE)

    allpreds <- tibble()
    for (predfname in predfnames) {

        print(predfname)

        currpred <- merge(positive_preds, read_tsv(predfname, col_names=cnames)) %>%
            dplyr::select(
              -any_of(
                c('hic_contactpromoter.quantile',
                  'hic_contactnonpromoter.quantile')
                )
              )

        allpreds <- rbind(allpreds, currpred)

    }

    return(allpreds)
}

addAncestry <- function(df) {

    afrpops <- c('ESN','GWD','LWK','YRI')
    eurpops <- c('CEU','FIN','IBS','TSI')

    df <- df %>%
        mutate(ancestry = case_when(population %in% afrpops ~ 'AFR',
                                    population %in% eurpops ~ 'EUR'))

    df <- df %>%
            mutate(population=fct_relevel(population,
                                          'CEU',
                                          'FIN',
                                          'IBS',
                                          'TSI',
                                          'ESN',
                                          'GWD',
                                          'LWK',
                                          'YRI'))

    return(df)
}

filtScoresSummPairs <- function(scores,
                                score_column='ABC.Score',
                                cdt_col='ancestry',
                                cdt1='AFR',
                                cdt2='EUR',
                                groupvars=c('chr','start','end','class','TargetGene','TargetGeneTSS',
                                        'TargetGeneExpression','distance','isSelfPromoter','name',
                                        'TargetGeneIsExpressed'),
                                min_nonzero_hic=16) {

    scores <- scores %>%
        dplyr::rename(test.Score := {{score_column}},
                      condition := {{cdt_col}}) %>%
        dplyr::filter(!is.na(test.Score)) %>%
        group_by(across(all_of(groupvars))) %>%
        dplyr::filter(length(test.Score)==min_nonzero_hic,
                      sum(hic_contact > 0) >= min_nonzero_hic) %>%
        dplyr::rename(!!score_column := test.Score) %>%
        # removes downstream LFC outliers driven by a few samples with no ATAC
        # reads or with no HiC contacts in all other enhancers in 5 Mb window
        dplyr::filter(all(ABC.Score < 1),
                      all(activity_base > 0)) %>%
        dplyr::rename(test.Score := {{score_column}})

    pairsumm <- scores %>%
        summarize(test.Score_cdt1 = mean(test.Score[condition==cdt1]),
                  test.Score_cdt2 = mean(test.Score[condition==cdt2]),
                  sd_cdt1 = sd(test.Score[condition==cdt1]),
                  sd_cdt2 = sd(test.Score[condition==cdt2]),
                  n_cdt1 = sum(condition==cdt1),
                  n_cdt2 = sum(condition==cdt2)) %>%
        ungroup() %>%
        mutate(sdadj_cdt1 = if_else(sd_cdt1 < .2*test.Score_cdt1, .2*test.Score_cdt1, sd_cdt1),
               sdadj_cdt2 = if_else(sd_cdt2 < .2*test.Score_cdt2, .2*test.Score_cdt2, sd_cdt2),
               log2FoldChange = log2(test.Score_cdt2/test.Score_cdt1),
               gsea_stat = (test.Score_cdt2 - test.Score_cdt1)/sqrt(sdadj_cdt2/n_cdt2 + sdadj_cdt1/n_cdt1)) %>%
        dplyr::rename("{score_column}_{cdt1}" := test.Score_cdt1,
                      "{score_column}_{cdt2}" := test.Score_cdt2,
                      "sd_{cdt1}" := sd_cdt1,
                      "sd_{cdt2}" := sd_cdt2,
                      "n_{cdt1}" := n_cdt1,
                      "n_{cdt2}" := n_cdt2,
                      "sdadj_{cdt1}" := sdadj_cdt1,
                      "sdadj_{cdt2}" := sdadj_cdt2)

    scores <- scores %>%
        ungroup() %>%
        dplyr::rename(!!score_column := test.Score,
                      !!cdt_col := condition)

    return(list(scores=scores, pairSummary=pairsumm))

}

parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--abc_positive_ccres',
                    action='store_true',
                    default=FALSE,
                    help='Use ABC.Score to define positive cCRE-gene pairs using provided thresholds, even when a different score column is provided for differential analysis.')
parser$add_argument('--no_header',
                    action='store_true',
                    default=FALSE,
                    help='Include if files in prediction_dir lack headers')
parser$add_argument('--min_nonzero_hic',
                    type='double',
                    default=16,
                    help='Number of samples with non-zero HiC/HiChIP contacts required to test E-G pair for diffABC.')
parser$add_argument('--test_activity',
                    action='store_true',
                    default=FALSE,
                    help="Also perform t-test on `activity_base` column for comparing results using normalized cCRE activity with scores normalized by window (e.g., 5 Mb).")
parser$add_argument('--prediction_dir',
                    type='character',
                    default='.',
                    help='Directory with parent directories for each sample containing `EnhancerPredictionsAllPutative.txt.gz` output from predict.py')
parser$add_argument('--qnorm_suffix',
                    type='character',
                    default='CEU_rep1QN',
                    help='Sample used for quantile normalization + `QN`.')
parser$add_argument('--enh_threshold',
                    type='double',
                    default=0.015,
                    help='Score to threshold on for enhancers.')
parser$add_argument('--promoter_threshold',
                    type='double',
                    default=0.1,
                    help='Score to threshold on for promoters.')
parser$add_argument('--score_column',
                    type='character',
                    default='ABC.Score',
                    help='Column name with score to use in thresholding and differential analysis.')
parser$add_argument('-n', '--name',
                    type='character',
                    default='AFR_EUR.diff',
                    help="Basename for output files")
parser$add_argument('--plotdir',
                    type='character',
                    default='/home/kpettie/code/github/plotting',
                    help="Directory with 'plotting.R' plotting functions.")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

abc_positive_ccres <- opt$abc_positive_ccres
no_header <- opt$no_header
min_nonzero_hic <- opt$min_nonzero_hic
test_activity <- opt$test_activity
pred_dir <- opt$prediction_dir
qn_suffix <- opt$qnorm_suffix
ethreshold <- opt$enh_threshold
pthreshold <- opt$promoter_threshold
score_column <- opt$score_column
outbase <- opt$name
outdir <- opt$outdir

source(file.path(opt$plotdir,'plotting.R'))

cdt_col <- 'ancestry'
cdt1 <- 'AFR'
cdt2 <- 'EUR'
groupvars <- c('chr','start','end','class','TargetGene','TargetGeneTSS',
              'TargetGeneExpression','distance','isSelfPromoter','name',
              'TargetGeneIsExpressed')
repcorrgvars <- c(groupvars, 'distance', 'isSelfPromoter', score_column, 'population', 'rep')

if (no_header) {
    # these colnames are out of date but all files should have header now
    cnames <- c(
      'chr',
      'start',
      'end',
      'name',
      'class',
      'TargetGene',
      'TargetGeneTSS',
      'TargetGeneExpression',
      'TargetGeneIsExpressed',
      'distance',
      'isSelfPromoter',
      'powerlaw_contact',
      'powerlaw_contact_reference',
      'hic_contact',
      'hic_enh_total', # check position
      'hic_tss_total', # check position
      'activity_base',
      'normalized_hic_contact',
      'hic_contact_pl_scaled',
      'hic_contact_sum1',
      'hic_contact_pl_scaled_adj',
      'ABC.Score.Numerator',
      'ABC.Score',
      'powerlaw.Score.Numerator',
      'powerlaw.Score',
      'CellType'
    )
} else {
  cnames <- TRUE
}

candidate_def_col <- score_column
if (abc_positive_ccres) candidate_def_col <- 'ABC.Score'

egpairs_allsamps <- positive_ccres_allsamps(pred_dir,
                                            qn_suffix,
                                            ethreshold=ethreshold,
                                            pthreshold=pthreshold,
                                            score_column=candidate_def_col,
                                            cnames=cnames)

print("Extracting positive predictions in at least one sample from all samples...")
scores_allsamps <- concat_positive_predictions(pred_dir, qn_suffix, egpairs_allsamps, cnames=cnames) %>%
    tidyr::separate(CellType, c('sample','qnorm'), sep='\\.', remove=FALSE) %>%
    tidyr::separate(sample, c('population','rep'), sep='_', remove=FALSE)
scores_allsamps <- addAncestry(scores_allsamps)

scores_allsamps %>%
    write_tsv(file.path(outdir,paste0(outbase,'.',score_column,".scoresAllSamps.txt.gz")))

# filter scores and get summary stats to match with t-test preds
fsScores <- filtScoresSummPairs(scores_allsamps,
                                score_column=score_column,
                                cdt_col=cdt_col,
                                cdt1=cdt1,
                                cdt2=cdt2,
                                groupvars=groupvars,
                                min_nonzero_hic=min_nonzero_hic)

fsScores$scores %>%
    write_tsv(file.path(outdir,paste0(outbase,'.',score_column,".scoresAllSampsFilt.txt.gz")))

scores_allsamps <- fsScores$scores %>%
    dplyr::rename(test.Score := {{score_column}},
                  condition := {{cdt_col}})

ttest_preds_all <- scores_allsamps %>%
    group_by(across(all_of(groupvars))) %>%
    t_test(test.Score ~ condition, p.adjust.method='bonferroni') %>%
    mutate(.y. = score_column)

ttest_preds <- merge(fsScores$pairSummary, ttest_preds_all)

ttest_preds %>%
    write_tsv(file.path(outdir,paste0(outbase,'.',score_column,".txt.gz")))

p1 <- groupedDistPlot(ttest_preds,
                      xvar='p',
                      groupvar='class',
                      nbins=150,
                      logscale=FALSE,
                    colorvals=c('red','green','blue'),
                      legendpos=c(.8,.8),
                      w=7,h=7)
p1
ggsave(file.path(outdir,paste0(outbase,'.',score_column,".pvalDist.png")), width=7,height=7)

p5 <- plotVolcano(ttest_preds,
                  lfc_var='log2FoldChange',
                  p_var='p',
                  pthresh='bonf',
                  nominal_p=TRUE,
                  lfcthresh=0.5,
                  label_var='TargetGene',
                  label_above = c(2.5,0),
                  xlabel='Log2 Fold Change',
                  ylabel='-log10(P)',
                  groupvar='class',
                  groupvarlab='cCRE class',
                  colorvals=NULL, # c("#2E3191","#F05A28") for orange, blue
                  legendpos=c(.8,.8),
                  w=7,h=7)
p5
ggsave(file.path(outdir,paste0(outbase,'.',score_column,".volcano.png")), width=7,height=7)

p6 <- plotDirectionByPval(ttest_preds,
                          direction_var='log2FoldChange',
                          p_var='p',
                          pvs=c(1,.5,0.05,0.005,5e-4,5e-5,0),
                          colorvals=c("#2E3191","#F05A28"),
                          updownlabs=c('EUR','AFR') # c(label_for_up_direction, label_for_down_direction)
                         )
p6
ggsave(file.path(outdir,paste0(outbase,'.',score_column,".directionByPval.png")), width=7,height=7)

scores_allsamps <- scores_allsamps %>%
    dplyr::rename(!!score_column := test.Score,
                  !!cdt_col := condition)

sasw <- scores_allsamps %>%
    dplyr::select(all_of(repcorrgvars)) %>%
    dplyr::rename(test.Score := {{score_column}}) %>%
    pivot_wider(names_from = rep,
                values_from = test.Score,
                names_prefix = paste0(score_column,'_'))

p7 <- plotCorrelation(sasw,
                            paste0(score_column,'_rep1'),
                            paste0(score_column,'_rep2'),
                            groupvar='class',
                            colorvals=NULL,
                            addline=c(0,1), # c(intercept, slope)
                            legendpos='bottom',
                            corrpos=c(.1,.9),
                            corrcolor='black',
                            w=22,h=14,
                            logscale=FALSE,
                            facets='population',
                            frows=2)
p7
ggsave(file.path(outdir,paste0(outbase,'.',score_column,".repsCorrelations.png")), width=22, height=14)


if (test_activity) {

    print("Testing `activity_base` column...")
    rm(fsScores)

    score_column <- 'activity_base'

    # filter scores and get summary stats to match with t-test preds
    fsScores <- filtScoresSummPairs(scores_allsamps,
                                    score_column=score_column,
                                    cdt_col=cdt_col,
                                    cdt1=cdt1,
                                    cdt2=cdt2,
                                    groupvars=groupvars,
                                    min_nonzero_hic=min_nonzero_hic)

    scores_allsamps <- fsScores$scores %>%
        dplyr::rename(test.Score := {{score_column}},
                      condition := {{cdt_col}})

    ttest_preds_all_base <- scores_allsamps %>%
        group_by(across(all_of(groupvars))) %>%
        t_test(test.Score ~ condition, p.adjust.method='bonferroni') %>%
        mutate(.y. = score_column)

    ttest_preds_base <- merge(fsScores$pairSummary, ttest_preds_all_base)

    ttest_preds_base %>%
        write_tsv(file.path(outdir,paste0(outbase,'.',score_column,".txt.gz")))

    p2 <- groupedDistPlot(ttest_preds_base,
                          xvar='p',
                          groupvar='class',
                          nbins=150,
                          logscale=FALSE,
                        colorvals=c('red','green','blue'),
                          legendpos=c(.8,.8),
                          w=7,h=7)
    p2
    ggsave(file.path(outdir,paste0(outbase,'.',score_column,".pvalDist.png")), width=7,height=7)

    ###### COMPARE ATAC ALONE TO ABC SCORE #########

    score_column <- opt$score_column

    ttest_comp <- merge(
        ttest_preds_base %>%
            dplyr::select(-c(.y., group1, group2, n1, n2)) %>%
            dplyr::rename(L2FC_ATAC = log2FoldChange,
                          statistic_ATAC = statistic,
                          df_ATAC = df,
                          p_ATAC = p),
        ttest_preds %>%
            dplyr::select(-c(.y., group1, group2, n1, n2)) %>%
            dplyr::rename(L2FC_ABC = log2FoldChange,
                          statistic_ABC = statistic,
                          df_ABC = df,
                          p_ABC = p)

    ) %>% dplyr::filter(!is.infinite(L2FC_ABC),
                        !is.infinite(L2FC_ATAC))

    p3 <- plotCorrelation(ttest_comp,
                xvar='L2FC_ATAC',
                yvar='L2FC_ABC',
                xlabel='Normalized ATAC L2FC',
                ylabel=paste0(score_column,' L2FC'),
                groupvar='class',
                colorvals=NULL,
                addline=c(0,1), # c(intercept, slope)
                legendpos=c(.2,.8),
                corrpos=c(5,3),
                w=7,h=7,
                logscale=FALSE)
    p3
    ggsave(file.path(outdir,paste0(outbase,'.',score_column,"vsATAC.LFC.png")), width=7,height=7)

    p4 <- plotCorrelation(ttest_comp %>%
                              mutate(nlp_ATAC = -log10(p_ATAC),
                                     nlp_ABC = -log10(p_ABC)),
                          xvar='nlp_ATAC',
                          yvar='nlp_ABC',
                          xlabel='Normalized ATAC -log10(P)',
                          ylabel=paste0(score_column, ' -log10(P)'),
                          groupvar='class',
                          colorvals=NULL,
                          addline=c(0,1), # c(intercept, slope)
                          legendpos=c(.9,.2),
                          corrpos=c(12,9),
                          w=7,h=7,
                          logscale=FALSE)

    p4
    ggsave(file.path(outdir,paste0(outbase,'.',score_column,"vsATAC.pvals.png")), width=7,height=7)


}
