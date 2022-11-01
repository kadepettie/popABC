#!/usr/bin/R

library(argparse)
library(tidyverse)

narrow_scores <- function(df,
                          groupvars=c('chr','start','end','class','TargetGene','TargetGeneTSS',
                                      'TargetGeneExpression','distance','isSelfPromoter','name',
                                      'TargetGeneIsExpressed','ABC.Score','atac.Score','chip.Score',
                                      'hic.Score','sample','population','rep','ancestry'),
                          cdt_col='ancestry',
                          subj_col='population',
                          debug=FALSE,
                          Ntestpairs=100,
                          Ntestsamps=16) {

    if (debug) {
        df <- df %>%
            select(any_of(groupvars)) %>%
            arrange(TargetGene,name) %>%
            head(Ntestpairs*Ntestsamps)
    }

    dfn <- df %>%
        rename(cdt := {{cdt_col}},
               subject := {{subj_col}}) %>%
        select(any_of(groupvars),subject,cdt) %>%
        pivot_longer(cols=ends_with('Score'),
                     names_to=c("score_type", ".value"),
                     names_pattern="(.*)\\.(.*)") %>%
        rename(!!cdt_col := cdt,
               !!subj_col := subject)


}

shuffle_cdt_reps_proxy <- function(scores,
                                  cdt_col='ancestry',
                                  subj_col='population',
                                  model_reps=TRUE,
                                  groupvars=c('chr','start','end','class','TargetGene','TargetGeneTSS',
                                                'TargetGeneExpression','distance','isSelfPromoter','name',
                                                'TargetGeneIsExpressed','score_type'),
                                  debug=FALSE,
                                  checkout=FALSE) {

    scores <- scores %>%
        dplyr::rename(cdt := {{cdt_col}},
                      subject := {{subj_col}})

    if (debug) {
        scores <- scores %>%
            dplyr::select(TargetGene,name,Score,sample,subject,rep,cdt) %>%
            arrange(TargetGene,name) %>%
            head(32)
    }

    # check if equal number of samples between conditions
    # to know if reps can be tested directly
    # if not, assign one subject from condition with more
    # samples randomly to one replicate as proxy
    cdtsamps <- scores %>%
        select(subject, sample, cdt) %>%
        distinct() %>%
        group_by(cdt) %>%
        mutate(Nsubjects = n_distinct(subject),
                  Nsamples = n_distinct(sample)) %>%
        ungroup()

    if (debug) return(cdtsamps)

    cdts <- cdtsamps %>% pull(cdt) %>% unique()
    Nsamps_cdt1 <- cdtsamps %>% filter(cdt==cdts[1]) %>% nrow()
    Nsamps_cdt2 <- cdtsamps %>% filter(cdt==cdts[2]) %>% nrow()

    # check if equal number of samples between conditions
    if (Nsamps_cdt1==Nsamps_cdt2) {
        scores <- scores %>%
            group_by(across(any_of(groupvars)), subject) %>%
            mutate(perm_reps = sample(rep)) %>%
            ungroup()
    # if not, assign one subject from condition with more
    # samples randomly to one replicate as proxy
    } else if (abs(Nsamps_cdt1-Nsamps_cdt2)==2) {
        if (Nsamps_cdt1>Nsamps_cdt2) {
            holdoutcdt <- cdts[1]
        } else {
            holdoutcdt <- cdts[2]
        }
        scores <- scores %>%
            group_by(across(any_of(groupvars))) %>%
            mutate(holdoutsub = sample(unique(subject[cdt==holdoutcdt]), 1)) %>%
            group_by(subject, .add=TRUE) %>%
            mutate(perm_reps = if_else(subject==holdoutsub, rep(sample(rep, 1), n()), sample(rep))) %>%
            ungroup()
    } else {
        stop('Check number of samples per condition and handle accordingly!')
    }

    if (checkout) {
        scores <- scores %>%
            dplyr::select(TargetGene,name,Score,sample,subject,rep,cdt,holdoutsub,perm_reps) %>%
            arrange(TargetGene,name)
        return(scores)
    }

    scores <- scores %>%
        dplyr::rename(!!cdt_col := cdt,
                      !!subj_col := subject)

    return(scores)

}

estimateFDR <- function(tdf, ndf, group_col='score_type', pcol='p', pthresh=0.05) {

    tsum <- as_tibble(tdf) %>%
        select(all_of(c(group_col, pcol))) %>%
        group_by(across(all_of(group_col))) %>%
        summarize(Nsig_obs = sum(.data[[pcol]] < pthresh, na.rm=TRUE))

    nsum <- as_tibble(ndf) %>%
        select(all_of(c(group_col, pcol))) %>%
        group_by(across(all_of(group_col))) %>%
        summarize(Nsig_exp = sum(.data[[pcol]] < pthresh, na.rm=TRUE))

    fdrdf <- merge(tsum, nsum) %>%
        mutate(pvthresh = pthresh,
               perm_fdr = Nsig_exp/Nsig_obs)

    return(fdrdf)
}

parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--scores_allsamps_fname',
                    type='character',
                    default='meanQN.16.AFR_EUR.diff.ABC.Score.scoresAllSamps.txt.gz')
parser$add_argument('--d_all_fname',
                    type='character',
                    default='allZerosFilt.meanQN.16.AFR_EUR.diff.allComponents.txt.gz')
parser$add_argument('--plotdir',
                    type='character',
                    default='/home/kpettie/code/github/plotting',
                    help="Directory with 'plotting.R' plotting functions.")
parser$add_argument('--fontdir',
                    type='character',
                    default='/cashew/shared_data/fonts',
                    help="Directory with .ttf font files.")
parser$add_argument('-n', '--name',
                    type='character',
                    default='allZerosFilt.meanQN.16.AFR_EUR.diff.allComponents',
                    help="Basename for output files")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

scores_allsamps_fname <- opt$scores_allsamps_fname
d_all_fname <- opt$d_all_fname
outdir <- opt$outdir

groupvars <- c('chr','start','end','class','TargetGene','TargetGeneTSS',
              'TargetGeneExpression','distance','isSelfPromoter','name',
              'TargetGeneIsExpressed','score_type')

d_all <- read_tsv(d_all_fname)
scores_allsamps <- read_tsv(scores_allsamps_fname)

print("Narrowing scores...")
snarr <- narrow_scores(scores_allsamps)
print("Shuffling reps...")
srep <- shuffle_cdt_reps_proxy(snarr,
                                  cdt_col='ancestry',
                                  subj_col='population',
                                  groupvars=groupvars)

print("Computing null P-values...")
srep_null <- srep %>%
    group_by(across(all_of(groupvars))) %>%
    rstatix::t_test(Score ~ perm_reps, p.adjust.method='bonferroni')

srep_null %>%
    write_tsv(file.path(outdir, paste0(opt$name, ".nullPvals.txt.gz")))

print("Estimating FDR...")
fdrest <- estimateFDR(d_all, srep_null)

fdrest %>%
    write_tsv(file.path(outdir, paste0(opt$name, ".FDRs.txt")))
