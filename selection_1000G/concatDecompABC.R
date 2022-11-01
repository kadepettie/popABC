#!/usr/bin/R

if(!(require(argparse) )) install.packages("argparse")
if(!(require(tidyverse) )) install.packages("tidyverse")
if(!(require(ggplot2)   )) install.packages("ggplot2")
if(!(require(gridExtra) )) install.packages("gridExtra")
if (!(require(msigdbr))) install.packages("msigdbr")

if(!(require(plyranges))) {
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("plyranges")
}
if (!require(fgsea)) {
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("fgsea")

}

library(argparse)
library(tidyverse)
library(plyranges)
library(ggplot2)
library(gridExtra)
library(msigdbr)
library(fgsea)

hgnc2ensembl <- function (df) {

    hgncgenes <- df %>% pull(hgnc_symbol) %>% unique()
    hgncgenes <- hgncgenes[which(!is.na(hgncgenes))]
    merge_col <- "hgnc_symbol"

    if(!(require(EnsDb.Hsapiens.v79))) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("EnsDb.Hsapiens.v79")
    }
    library(EnsDb.Hsapiens.v79)

    gene_info <- ensembldb::select(EnsDb.Hsapiens.v79,
                                   keys = hgncgenes,
                                   keytype = "SYMBOL",
                                   columns = c("SYMBOL", "GENEID", "ENTREZID")) %>%
        as_tibble() %>%
        dplyr::rename(ensembl_gene_id = GENEID,
                      hgnc_symbol = SYMBOL,
                      entrez_id = ENTREZID) %>%
        dplyr::filter(grepl('ENSG', ensembl_gene_id)) # filter out alternate gene symbols

    print("Adding ensembl_gene_ids...")
    df <- merge(df, gene_info, by = merge_col, all.x = TRUE)
    return(df)

}

prepABCforGSEA <- function(df) {
    # converts TargetGenes (hgnc) to unique ensIDs

    # some hgcn's have multiple ensg's but usually lower ID number is in ENSG
    # keep just one
    df <- df %>%
        dplyr::rename(hgnc_symbol = TargetGene) %>%
        hgnc2ensembl(.) %>%
        group_by(hgnc_symbol,name) %>%
        arrange(ensembl_gene_id, entrez_id, .by_group=TRUE) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        arrange(p)

    return(df)

}

direction_matches <- function(x) {
    if (all(x >= 0) | all(x < 0)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

addEnhancerDirectionMatching <- function(df,
                                  dist_thresh=NULL,
                                  p_thresh=0.05,
                                  abc_p_col='p',
                                  abc_lfc_col='log2FoldChange',
                                  gene_col='TargetGene',
                                  nEnh=2,
                                  score_type='Score',
                                  cdt1 = 'AFR',
                                  cdt2 = 'EUR',
                                  allscores = FALSE,
                                  abcrank = FALSE,
                                  debug=FALSE) {

    # gprof_p = threshold for returning gprofiler results

    score_col1 <- paste0(score_type,'_',cdt1)
    score_col2 <- paste0(score_type,'_',cdt2)

    if (allscores & abcrank) {

        meanabc <- df %>%
            dplyr::filter(score_type=='ABC') %>%
            mutate(mean.Score = (.data[[score_col1]] + .data[[score_col2]])/2) %>%
            dplyr::select(TargetGene,
                          seqnames,
                          start,
                          end,
                          class,
                          name,
                          distance,
                          isSelfPromoter,
                          mean.Score)

        dfs <- merge(meanabc, df)

    } else {

        dfs <- df %>%
            mutate(mean.Score := (.data[[score_col1]] + .data[[score_col2]]) / 2)

    }

    dfs <- dfs %>%
        dplyr::rename(abc_p := {{abc_p_col}},
                      abc_lfc := {{abc_lfc_col}},
                      gene := {{gene_col}}) %>%
        dplyr::filter(abc_p < p_thresh) %>%
        group_by(gene)

    if (allscores) dfs <- dfs %>% group_by(score_type, .add=TRUE)

    if (is.numeric(nEnh)) {
        dfs <- dfs %>% dplyr::filter(n()==nEnh)
    } else {
        dfs <- dfs %>% dplyr::filter(n() > 1)
    }

    dfs <- dfs %>%
        arrange(desc(mean.Score), .by_group=TRUE) %>%
        mutate(enh_idx = 1:n(),
               enh_type = if_else(enh_idx==1,'primary','secondary'),
               secondary_matches = if_else(direction_matches(abc_lfc[enh_type=='secondary']), TRUE, FALSE)) %>%
        dplyr::filter(secondary_matches)

    if (debug) return(dfs)

    dfs <- dfs %>%
        mutate(direction = case_when(all(abc_lfc > 0)                              ~ "up_up",
                                     all(abc_lfc[enh_type=='primary'] > 0)
                                         & all(abc_lfc[enh_type=='secondary'] < 0) ~ "up_down",
                                     all(abc_lfc[enh_type=='primary'] < 0)
                                         & all(abc_lfc[enh_type=='secondary'] > 0) ~ "down_up",
                                     all(abc_lfc < 0)                              ~ "down_down"),
               binary_direction = if_else(direction %in% c('up_up','down_down'), 'reinforcing', 'opposing')
              )

    if (!is.null(dist_thresh)) {

        tgenes <- dfs %>%
            group_by(binary_direction, .add=TRUE) %>%
            summarize(enh_distance = max(distance) - min(distance)) %>%
            ungroup() %>%
            dplyr::filter(enh_distance > dist_thresh) %>%
            pull(gene)

        dfs <- dfs %>%
            dplyr::filter(gene %in% tgenes)

    }

    return(dfs)

}

snp_selection_diff_score_wilcox <- function(df,
                                           diff_enh_pairs=FALSE,
                                           snp_stat='iHS',
                                           dabc_thresh=0.05,
                                           dabc_col='p',
                                           lfc_col='log2FoldChange',
                                           perc_thresh=0.99,
                                           perc_col='rank2_percentile',
                                           rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                           by_gene=FALSE, # count by gene number instead of enhancer number
                                           all_or_none=FALSE,
                                           run_wilcox=FALSE,
                                           wilcox_col='max_abs_ihs',
                                           gene_col='TargetGene',
                                           allscores=FALSE,
                                           score_col=NULL,
                                           debug=FALSE) {

    # output dataframe for plotting snp selection stats (iHS or FST) in diff- vs
    # non-diff-E-G pairs with option to output wilcoxon test results
    # includes some legacy options from previous function incorporating
    # fisher's exact test on top selection candidates by percentile
    # vs not top selection candidates

    # all_or_none requires all SNPs overlapping enhancers
    # or none to be top selection candidates for the gene
    # to be included in the test

    # if diff_enh_pairs==FALSE (not restricted to genes with 2 diff-enhancers) and by_gene==TRUE,
    # then test is on top selection SNP enhancer (diff or non-diff) per gene

    dfs <- df %>%
        dplyr::rename(snp_percentile := {{perc_col}},
                      gene := {{gene_col}},
                      log2FoldChange := {{lfc_col}},
                      p := {{dabc_col}},
                      rank_var := {{rank_col}}) %>%
        # results in uneven counting if by_gene==FALSE and diff_enh_pairs==TRUE
        # thus, by_gene should be the preferred method for diff_enh_pairs (reinforcing vs opposing direction)
        dplyr::filter(!is.na(snp_percentile)) %>%
        mutate(selection_candidate = if_else(snp_percentile >= perc_thresh, TRUE, FALSE))

    if (diff_enh_pairs) {
        enhtype <- 'pairs'
        dfs <- dfs %>%
            mutate(diff_status = if_else(direction %in% c('up_up','down_down'), 'diff', 'non-diff'))
    } else {
        enhtype <- 'singletons'
        dfs <- dfs %>%
            mutate(diff_status = case_when(p < dabc_thresh ~ 'diff',
                                           p >= 0.5        ~ 'non-diff')) %>%
            dplyr::filter(!is.na(diff_status))
    }

    tset <- 'enh-level'

    if (by_gene) {

        tset <- 'gene-level'

        dfs <- dfs %>%
            group_by(gene)
        if (allscores) dfs <- dfs %>% group_by(score_type, .add=TRUE)
        if (all_or_none) {
            dfs <- dfs %>%
                mutate(all_or_none = if_else(all(selection_candidate) | !any(selection_candidate), TRUE, FALSE))
        }
        dfs <- dfs %>%
            arrange(desc(rank_var), .by_group=TRUE) %>%
            dplyr::slice(1) %>%
            ungroup()

        if (all_or_none) {
            dfs <- dfs %>%
                dplyr::filter(all_or_none)

        }
    }

    dfs <- dfs %>% dplyr::rename(!!rank_col := rank_var)

    if (allscores) dfs <- dfs %>% group_by(score_type)

    if (run_wilcox) {

        wilcox_df <- dfs %>%
            dplyr::rename(wilcox_test_var := wilcox_col)
        wilcox_summ <- wilcox_df %>%
            group_by(diff_status, .add=TRUE) %>%
            summarize(N = n(),
                      mean_stat = mean(wilcox_test_var)) %>%
            ungroup()

        if (allscores) {

            noTestScoreTypes <- wilcox_summ %>%
                group_by(score_type) %>%
                summarize(N = dplyr::n()) %>%
                dplyr::filter(N < 2) %>%
                pull(score_type)

            wilcox_summf <- wilcox_summ %>%
                dplyr::filter(!(score_type %in% noTestScoreTypes))

            wilcox_res <- wilcox_df %>%
                dplyr::filter(!(score_type %in% noTestScoreTypes)) %>%
                wilcox_test(wilcox_test_var ~ diff_status) %>%
                mutate(wilcox_var = wilcox_col,
                       mean_stat_fg = wilcox_summf[wilcox_summf$diff_status=='diff',][['mean_stat']],
                       mean_stat_bg = wilcox_summf[wilcox_summf$diff_status=='non-diff',][['mean_stat']]) %>%
                dplyr::rename(wilcox_P = p) %>%
                dplyr::select(any_of('score_type'),wilcox_var,mean_stat_fg,mean_stat_bg,wilcox_P)

            wilcox_summnt <- wilcox_summ %>%
                dplyr::filter(score_type %in% noTestScoreTypes)

            msf <- wilcox_summnt[wilcox_summnt$diff_status=='diff',][['mean_stat']]
            msb <- wilcox_summnt[wilcox_summnt$diff_status=='non-diff',][['mean_stat']]
            if (length(msf)==0) msf <- NA
            if (length(msb)==0) msb <- NA

            wilcox_res_nt <- tibble(score_type=wilcox_summnt$score_type,
                                    wilcox_var=wilcox_col,
                                    mean_stat_fg = msf,
                                    mean_stat_bg = msb,
                                    wilcox_P = NA)

            if (debug) return(wilcox_res_nt)

            wilcox_res <- rbind(wilcox_res, wilcox_res_nt)


        } else {

            # check for N_diff_status == 2

            if (nrow(wilcox_summ)==2) {

                wilcox_res <- wilcox_df %>%
                    wilcox_test(wilcox_test_var ~ diff_status) %>%
                    mutate(wilcox_var = wilcox_col,
                           mean_stat_fg = wilcox_summ[wilcox_summ$diff_status=='diff',][['mean_stat']],
                           mean_stat_bg = wilcox_summ[wilcox_summ$diff_status=='non-diff',][['mean_stat']]) %>%
                    dplyr::rename(wilcox_P = p) %>%
                    dplyr::select(wilcox_var,mean_stat_fg,mean_stat_bg,wilcox_P)

            } else {

                msf <- wilcox_summ[wilcox_summ$diff_status=='diff',][['mean_stat']]
                msb <- wilcox_summ[wilcox_summ$diff_status=='non-diff',][['mean_stat']]
                if (length(msf)==0) msf <- NA
                if (length(msb)==0) msb <- NA

                wilcox_res <- tibble(wilcox_var=wilcox_col,
                                        mean_stat_fg = msf,
                                        mean_stat_bg = msb,
                                        wilcox_P = NA)

            }

        }

        wilcox_res <- wilcox_res %>%
            mutate(snp_stat = snp_stat, enh_type = enhtype, test_set = tset)

        if (debug) return(wilcox_res)

    }

    dfs <- ungroup(dfs)

    if ('binary_direction' %in% colnames(dfs)) {
        dfs <- dfs %>% mutate(enh_status = binary_direction)
    } else {
        dfs <- dfs %>% mutate(enh_status = diff_status)
    }
    dfs <- dfs %>%
        mutate(snp_stat = snp_stat, enh_type = enhtype, test_set = tset) %>%
        dplyr::rename(wilcox_var := {{wilcox_col}})

    plotcols <- c('seqnames','start','end','hgnc_symbol','gene','ensembl_gene_id',
                  'entrez_id','class','TargetGeneTSS','distance','isSelfPromoter',
                  'name','log2FoldChange','p','RSNUM','snp_percentile','selection_candidate','enh_status',
                  'wilcox_var','snp_stat','enh_type','test_set','score_type')

    dfs <- dfs %>% dplyr::select(any_of(plotcols))

    if (!is.null(score_col)) {
        dfs <- dfs %>% mutate(score_column = score_col)
        eres <- eres %>% mutate(score_column = score_col)
    }

    print(colnames(dfs))

    if (run_wilcox) {
        return(list(fw=wilcox_res, plotdf=dfs))
    } else {
        return(dfs)
    }


}

snp_selection_diff_score_fisher <- function(df,
                                           diff_enh_pairs=FALSE,
                                           snp_stat='iHS',
                                           dabc_thresh=0.05,
                                           dabc_col='p',
                                           lfc_col='log2FoldChange',
                                           perc_thresh=0.99,
                                           mid_filt=FALSE,
                                           alt='greater',
                                           perc_col='rank2_percentile',
                                           rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                           by_gene=FALSE, # count by gene number instead of enhancer number
                                           all_or_none=FALSE,
                                           gene_col='TargetGene',
                                           allscores=FALSE,
                                           score_col=NULL,
                                           mat_out=FALSE,
                                           debug=FALSE) {

    # fisher's exact test on diff-score/reinforcing enhancer pair enrichment in
    # top selection candidates by percentile vs not top selection candidates

    # all_or_none requires all SNPs overlapping enhancers
    # or none to be top selection candidates for the gene
    # to be included in the test

    # if diff_enh_pairs==FALSE (not restricted to genes with 2 diff-enhancers) and by_gene==TRUE,
    # then test is on top selection SNP enhancer (diff or non-diff) per gene

    dfs <- df %>%
        dplyr::rename(snp_percentile := {{perc_col}},
                      rank_var := {{rank_col}},
                      gene := {{gene_col}},
                      log2FoldChange := {{lfc_col}},
                      p := {{dabc_col}}) %>%
        # results in uneven counting if by_gene==FALSE and diff_enh_pairs==TRUE
        # thus, by_gene should be the preferred method for diff_enh_pairs (reinforcing vs opposing direction)
        dplyr::filter(!is.na(snp_percentile)) %>%
        mutate(selection_candidate = if_else(snp_percentile >= perc_thresh, TRUE, FALSE))

    if (diff_enh_pairs) {
        enhtype <- 'pairs'
        dfs <- dfs %>%
            mutate(diff_status = if_else(direction %in% c('up_up','down_down'), 'diff', 'non-diff'))
    } else {
        enhtype <- 'singletons'
        dfs <- dfs %>%
            mutate(diff_status = case_when(p < dabc_thresh ~ 'diff',
                                           p >= 0.5        ~ 'non-diff')) %>%
            dplyr::filter(!is.na(diff_status))
    }

    tset <- 'enh-level'

    if (by_gene) {

        tset <- 'gene-level'

        dfs <- dfs %>%
            group_by(gene)
        if (allscores) dfs <- dfs %>% group_by(score_type, .add=TRUE)
        if (all_or_none) {
            dfs <- dfs %>%
                mutate(all_or_none = if_else(all(selection_candidate) | !any(selection_candidate), TRUE, FALSE))
        }
        dfs <- dfs %>%
            arrange(desc(rank_var), .by_group=TRUE) %>%
            dplyr::slice(1) %>%
            ungroup()

        if (all_or_none) {
            dfs <- dfs %>%
                dplyr::filter(all_or_none)

        }
    }

    lowstat <- perc_thresh
    if (mid_filt) {
        lowstat <- 0.5
        dfs <- dfs %>% dplyr::filter(selection_candidate | snp_percentile <= lowstat)
    }

    if (allscores) {

        score_types <- dfs %>% pull(score_type) %>% unique()

        eres <- tibble()
        for (st in score_types) {
            dfss <- dfs %>% dplyr::filter(score_type==st)

            up_up <- dfss %>%
                dplyr::filter(diff_status=='diff',
                              selection_candidate) %>%
                nrow()

            up_down <- dfss %>%
                dplyr::filter(diff_status=='diff',
                              !selection_candidate) %>%
                nrow()

            down_up <- dfss %>%
                dplyr::filter(diff_status=='non-diff',
                              selection_candidate) %>%
                nrow()

            down_down <- dfss %>%
                dplyr::filter(diff_status=='non-diff',
                              !selection_candidate) %>%
                nrow()

            mat11 <- up_up
            mat21 <- up_down
            mat12 <- down_up
            mat22 <- down_down
            tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                              dimnames =
                       list(c(paste0('top ',snp_stat), paste0('not top ',snp_stat)),
                            c('diff', 'non-diff')))


            print(tmat)
            ft <- fisher.test(tmat, alternative=alt)

            if (mat_out) {
                return(
                    list(
                        tmat,
                        ft
                    )
                )
            }

            curreres <- list(score_type = st,
                             snp_stat = snp_stat,
                             enh_type = enhtype,
                             test_set = tset,
                             dabc_significance = dabc_thresh,
                             top_stat_percentile = perc_thresh,
                             low_stat_percentile = lowstat,
                             odds_ratio = ft$estimate[[1]],
                             pvalue = ft$p.value,
                             conf_lower = ft$conf.int[1],
                             conf_upper = ft$conf.int[2],
                             up_up = mat11,
                             up_down = mat21,
                             down_up = mat12,
                             down_down = mat22)

            eres <- rbind(eres, curreres)

        }

    } else {

        up_up <- dfs %>%
            dplyr::filter(diff_status=='diff',
                          selection_candidate) %>%
            nrow()

        up_down <- dfs %>%
            dplyr::filter(diff_status=='diff',
                          !selection_candidate) %>%
            nrow()

        down_up <- dfs %>%
            dplyr::filter(diff_status=='non-diff',
                          selection_candidate) %>%
            nrow()

        down_down <- dfs %>%
            dplyr::filter(diff_status=='non-diff',
                          !selection_candidate) %>%
            nrow()

        mat11 <- up_up
        mat21 <- up_down
        mat12 <- down_up
        mat22 <- down_down
        tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                          dimnames =
                   list(c(paste0('top ',snp_stat), paste0('not top ',snp_stat)),
                        c('diff', 'non-diff')))


        print(tmat)
        ft <- fisher.test(tmat, alternative=alt)

        if (mat_out) {
            return(
                list(
                    tmat,
                    ft
                )
            )
        }

        eres <- list(snp_stat = snp_stat,
                     enh_type = enhtype,
                     test_set = tset,
                     dabc_significance = dabc_thresh,
                     top_stat_percentile = perc_thresh,
                     low_stat_percentile = lowstat,
                     odds_ratio = ft$estimate[[1]],
                     pvalue = ft$p.value,
                     conf_lower = ft$conf.int[1],
                     conf_upper = ft$conf.int[2],
                     up_up = mat11,
                     up_down = mat21,
                     down_up = mat12,
                     down_down = mat22)

    }


    dfs <- ungroup(dfs)

    if ('binary_direction' %in% colnames(dfs)) {
        dfs <- dfs %>% mutate(enh_status = binary_direction)
    } else {
        dfs <- dfs %>% mutate(enh_status = diff_status)
    }
    dfs <- dfs %>% mutate(snp_stat = snp_stat, enh_type = enhtype, test_set = tset)

    if (!is.null(score_col)) {
        dfs <- dfs %>% mutate(score_column = score_col)
        eres <- eres %>% mutate(score_column = score_col)
    }

    print(colnames(as_tibble(eres)))

    return(eres)


}

mergeDecompDiffScores <- function(diffdir,
                                  diffpatt='meanQN\\.16\\.AFR_EUR\\.diff\\..*\\.Score\\.txt\\.gz') {

    fnames <- list.files(diffdir,
                         pattern=diffpatt,
                         full.names=TRUE)

    remove_cols <- c('sd_AFR','sd_EUR','sdadj_AFR','sdadj_EUR','.y.','statistic','df')

    adf <- tibble()
    for (f in fnames) {

        print(paste0('Reading ',f))

        df <- read_tsv(f)
        st <- str_split(df[['.y.']][1],pattern='\\.',simplify=TRUE)[1]
        df <- df %>%
            dplyr::select(-all_of(remove_cols)) %>%
            mutate(score_type = st) %>%
            dplyr::rename(Score_AFR := paste0(st,'.Score_AFR'),
                          Score_EUR := paste0(st,'.Score_EUR'))

        adf <- rbind(adf, df)

    }

    return(adf)
}

mergeEnhReinfFisherEnrichments <- function(rdir,
                                           rpatt='.enhReinforcingDirection.fisherEnrichments.txt.gz') {

    fnames <- list.files(rdir, pattern=rpatt, full.names=TRUE)

    adf <- tibble()
    for (f in fnames) {
        print(paste0('Reading ', f))

        adf <- rbind(adf, read_tsv(f))
    }

    return(
        adf %>%
            mutate(condition = if_else(datatype=='ABC_overall',
                                       'ABC_overall',
                                       condition)) %>%
            distinct()
    )

}

mergeSnpSelectionFisherEnrichments <- function(adir,
                                               spatt='.selectionSNPs.fisherEnrichmentsWilcox.txt') {

    fnames <- list.files(adir, pattern=spatt, full.names=TRUE)

    adf <- tibble()
    for (f in fnames) {
        print(paste0('Reading ', f))

        adf <- rbind(adf, read_tsv(f))
    }

    return(
        adf %>% distinct()
    )

}

readTopPairFisherEnrichments <- function(fdir,
                                         fe_patt='fisherEnrichments.txt$') {

    fnames <- list.files(fdir,
                         pattern=fe_patt,
                         full.names=FALSE)

    adf <- tibble()
    for (f in fnames) {

        print(paste0('Reading ',f))

        df <- read_tsv(file.path(fdir, f))

        if (!('score_type' %in% colnames(df))) {
            scoretype <- str_split(f, '\\.')[[1]][2:3] %>% paste(., collapse='.')
            df <- df %>% mutate(score_type = scoretype)
        }

        adf <- rbind(adf, df)

    }

    return(adf)

}

add_ihs <- function(df, ihsf) {

    if (!is.null(ihsf)) {

        include_cols <- c('seqnames',
                          'start',
                          'RSNUM',
                          'max_abs_ihs',
                          'rank2_abs_ihs',
                          'rank2_percentile',
                          'selection_candidate')

        groupcols <- c('score_type','TargetGene','name')

        ihs <- read_tsv(ihsf, col_select=any_of(include_cols)) %>%
            mutate(ihs_pos = start)

        if (!('seqnames' %in% colnames(df))) {
            df <- df %>% dplyr::rename(seqnames=chr)
        }

        df <- join_overlap_left(
                    df %>% as_granges(),
                    ihs %>% mutate(end = start + 1) %>% as_granges()
                ) %>%
                as_tibble() %>%
                group_by(across(all_of(groupcols))) %>%
                arrange(desc(rank2_abs_ihs), .by_group=TRUE) %>%
                dplyr::slice(1) %>%
                ungroup()

    }

    return(df)

}

add_fst <- function(df, fstf) {

    if (!is.null(fstf)) {

        include_cols <- c('seqnames',
                          'start',
                          'FST',
                          'fst_percentile')

        groupcols <- c('score_type','TargetGene','name')

        fst <- read_tsv(fstf, col_select=any_of(include_cols)) %>%
            mutate(FST_pos = start)

        if (!('seqnames' %in% colnames(df))) {
            df <- df %>% dplyr::rename(seqnames=chr)
        }

        df <- join_overlap_left(
                    df %>% as_granges(),
                    fst %>% mutate(end = start + 1) %>% as_granges()
                ) %>%
                as_tibble() %>%
                group_by(across(all_of(groupcols))) %>%
                arrange(desc(FST), .by_group=TRUE) %>%
                dplyr::slice(1) %>%
                ungroup()

    }

    return(df)

}

# also in plotting.R, but required for main gsea function

plotTopPways <- function(res, topN=5, cwidths=c(5, 3, 0.6, .8, .8)) {

    topPathwaysUp <- res$res_df[ES > 0][head(order(pval), n=topN), pathway]
    topPathwaysDown <- res$res_df[ES < 0][head(order(pval), n=topN), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    p1 <- plotGseaTable(res$pwaylist[topPathways],
                        res$ranklist,
                        res$res_df,
                        gseaParam=0.5,
                        colwidths=cwidths,
                        render=FALSE)

    return(p1)

}

gridPlotTopPways <- function(res1, res2, score_column='ABC.Score') {

    p1 <- plotTopPways(res1)
    p2 <- plotTopPways(res2, cwidths=c(7, 3, 0.6, .8, .8))

    options(repr.plot.width = 21, repr.plot.height = 5, repr.plot.res = 200)

    g1 <- grid.arrange(
          p2,
          p1,
          nrow = 1,
          top = grid::textGrob(
               score_column,
               gp = grid::gpar(fontsize = 26)
          ),
          padding = unit(2, "line")
        )

    return(g1)

}

run_gsea <- function(df,
                     genesets=NULL,
                     gene_col='entrez_id',
                     gsea_col='gsea_stat',
                     unique_enh=FALSE,
                     return_plot_data=FALSE,
                     debug=FALSE) {

    if (!(gene_col %in% colnames(df))) {
        df <- prepABCforGSEA(df)
    }

    df <- df %>%
        dplyr::rename(gene := {{gene_col}},
                      gsea_var := {{gsea_col}}) %>%
        dplyr::filter(!is.na(gene)) %>%
        group_by(gene) %>%
        arrange(desc(abs(gsea_var)), .by_group=TRUE) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        arrange(gsea_var)

    if (unique_enh) {
        df <- df %>%
            group_by(name) %>%
            # choose gene semi-arbitrarily for now (lower ens number)
            arrange(gene) %>%
            dplyr::slice(1) %>%
            ungroup()
    }

    if (debug) return(df)

    rankstat <- setNames(df$gsea_var, df$gene)

    if (!require(fgsea)) {
        if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")

        BiocManager::install("fgsea")

    }
    library(fgsea)
    if (is.null(genesets)) {
        print("Using example pathways")
        data(examplePathways)
        msigdbr_list <- examplePathways
    } else {
        if (gene_col=='entrez_id') {
            gscol <- 'entrez_gene'
        } else if (gene_col=='ensembl_gene_id') {
            gscol <- 'ensembl_gene'
        } else if (gene_col %in% c('hgnc_symbol','TargetGene')) {
            gscol <- 'gene_symbol'
        } else {
            stop(paste('Gene col',gene_col,'not supported!'))
        }
        msigdbr_list <- split(x = genesets[[gscol]], f = genesets$gs_name)
    }

    res <- fgsea(pathways = msigdbr_list,
                  stats    = rankstat,
                   minSize = 15, maxSize = 500) %>%
        arrange(pval)

    if (return_plot_data) {
        # pwaylist and ranklist used as input to plotEnrichment
        # for characteristic GSEA single pathway enrichment plots
        return(list(res_df=res, pwaylist=msigdbr_list, ranklist=rankstat))
    }

    return(res)

}

aggregate_gsea <- function(d, allgenesets=NULL, remove_antisense=FALSE) {

    if (is.null(allgenesets)) {

        if (!(require(msigdbr))) install.packages("msigdbr")
        library(msigdbr)

        allgenesets <- msigdbr("Homo sapiens")

    }

    adg <- d %>%
        group_by(score_type, TargetGene) %>%
        arrange(p, desc(abs(gsea_stat)), .by_group=TRUE) %>%
        dplyr::slice(1) %>%
        ungroup()

    if (remove_antisense) {

        # remove antisense genes whos' sense genes have the same value

        print("Removing redundant antisense genes...")
        adg <- adg %>%
            separate(TargetGene,
                     c('senseGene','antisenseNum'),
                     sep='-AS',
                     remove=FALSE,
                     convert=TRUE,
                     fill='right') %>%
            group_by(score_type, gsea_stat) %>%
            mutate(N_AS_pairs = length(senseGene) - n_distinct(senseGene)) %>%
            ungroup() %>%
            dplyr::filter(!(!is.na(antisenseNum) & N_AS_pairs > 0)) %>%
            dplyr::select(-all_of(c('senseGene','antisenseNum')))

    }

    sts <- adg %>% pull(score_type) %>% unique()

    ares_df <- tibble()

    for (i in seq_along(sts)) {

        st <- sts[i]

        print(paste0('Processing diff-', st))

        dg <- adg %>% dplyr::filter(score_type==st)

        print('    Hallmark enrichments...')
        hreslist <- run_gsea(dg,
                            genesets=allgenesets %>%
                                dplyr::filter(gs_cat=='H'),
                            gene_col='TargetGene',
                            gsea_col='gsea_stat',
                            unique_enh=FALSE,
                            return_plot_data=TRUE,
                            debug=FALSE)

        print('    GO:BP enrichments...')
        bpreslist <- run_gsea(dg,
                            genesets=allgenesets %>%
                                dplyr::filter(gs_subcat=='GO:BP'),
                            gene_col='TargetGene',
                            gsea_col='gsea_stat',
                            unique_enh=FALSE,
                            return_plot_data=TRUE,
                            debug=FALSE)

        cres_df <- rbind(
            hreslist$res_df %>%
                mutate(gs_category='Hallmark'),
            bpreslist$res_df %>%
                mutate(gs_category='GO:BP')
        ) %>%
            mutate(score_type=st)

        ares_df <- rbind(ares_df, cres_df)

        g1 <- gridPlotTopPways(hreslist, bpreslist, score_column=st)
        if (i==1) {
            agg_g <- g1
        } else {
            agg_g <- rbind(agg_g, g1)
        }
        i = i + 1

    }

    return(list(res_df=ares_df, plot_obj=agg_g))

}

enhDirTestFromPairs <- function(df,
                                dist_thresh=NULL,
                                mat_out=FALSE,
                                alt='greater',
                                debug=FALSE) {

    dfs <- df %>% group_by(gene)

    if (!is.null(dist_thresh)) {

        tgenes <- dfs %>%
            group_by(binary_direction, .add=TRUE) %>%
            summarize(enh_distance = max(distance) - min(distance)) %>%
            ungroup() %>%
            dplyr::filter(enh_distance > dist_thresh) %>%
            pull(gene)

        dfs <- dfs %>%
            dplyr::filter(gene %in% tgenes)

    }

    up_up <- dfs %>%
        dplyr::filter(direction=='up_up') %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        nrow()

    up_down <- dfs %>%
        dplyr::filter(direction=='up_down') %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        nrow()

    down_up <- dfs %>%
        dplyr::filter(direction=='down_up') %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        nrow()

    down_down <- dfs %>%
        dplyr::filter(direction=='down_down') %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        nrow()

    mat11 <- up_up
    mat21 <- up_down
    mat12 <- down_up
    mat22 <- down_down
    tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                      dimnames =
               list(c('secondary up', 'secondary down'),
                    c('primary up', 'primary down')))

    ft <- fisher.test(tmat, alternative=alt)

    if (mat_out) {
        return(
            list(
                tmat,
                ft
            )
        )
    }

    eres <- list(odds_ratio = ft$estimate[[1]],
                 pvalue = ft$p.value,
                 conf_lower = ft$conf.int[1],
                 conf_upper = ft$conf.int[2],
                 up_up = mat11,
                 up_down = mat21,
                 down_up = mat12,
                 down_down = mat22)

    return(eres)

}

aggregate_enhdirtest <- function(df, dist_thresh=NULL) {

    sts <- df %>% pull(score_type) %>% unique()

    ares <- tibble()

    for (st in sts) {

        cdf <- df %>% dplyr::filter(score_type==st)

        res <- enhDirTestFromPairs(cdf, dist_thresh=dist_thresh) %>%
            as_tibble() %>%
            mutate(score_type=st)

        ares <- rbind(ares, res)

    }

    return(ares)

}

parser <- ArgumentParser(description='Aggregate differential test results of ABC Score components for QC/diagnostics and spot checking.')

parser$add_argument('--dist_thresh',
                    type='double',
                    default=NULL,
                    help='Require reinforcing/opposing enhancers to be at least this far apart (default of 0 means ~500 bp enhancers can be separated by 1 bp)')
parser$add_argument('--FST',
                    type='character',
                    default=NULL,
                    help="File of SNPs with FST values calculate between 4 pops of each ancestry.")
parser$add_argument('--iHS',
                    type='character',
                    default=NULL,
                    help="File of SNPs with iHS scores for at least 2 pops and percentile of |iHS| for second highest |iHS| per SNP (`rank2_percentile`).")
parser$add_argument('--reinf_pattern',
                    type='character',
                    default=".enhReinforcingDirection.fisherEnrichments.txt.gz",
                    help="Pattern uniquely identifying enhancer reinforcing direction/DE matching direction files to aggregate from agg_dir (glob with partial matching allowed, not regex, to avoid escape characters on command line).")
parser$add_argument('-n', '--name',
                    type='character',
                    default='diffABC_DE_overlap',
                    help="Basename for output files")
parser$add_argument('--plotdir',
                    type='character',
                    default=NULL, # /home/kpettie/code/github/plotting
                    help="Directory with plotting function in 'plotting.R'. Required to make plots")
parser$add_argument('--diff_pattern',
                    type='character',
                    default=".Score.txt.gz",
                    help="Pattern uniquely identifying differntial score files to aggregate from agg_dir (glob with partial matching allowed, not regex, to avoid escape characters on command line).")
parser$add_argument('--agg_dir',
                    type='character',
                    default='.',
                    help="Directory with files from different score component and overlap analyses for aggregating")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

dist_thresh <- opt$dist_thresh
ihs_fname <- opt$iHS
fst_fname <- opt$FST
plotdir <- opt$plotdir
r_patt <- opt$reinf_pattern
d_patt <- opt$diff_pattern
a_dir <- opt$agg_dir
outdir <- opt$outdir
outbase <- opt$name

envPlotColors <- tribble(
                    ~condition,  ~hex, ~category,
                    'GARD',  "#CC79A7", 'Immune stimulant',
                    'FSL1',  "#CC79A7", 'Immune stimulant',
                    'IFNG',  "#CC79A7", 'Immune stimulant',
                    'BAFF',  "#CC79A7", 'Immune stimulant',
                    'DEX',   "#56B4E9", 'Hormone',
                    'IGF',   "#56B4E9", 'Hormone',
                    'ACRYL', "#009E73", 'Contaminant or stressor',
                    'PFOA',  "#009E73", 'Contaminant or stressor',
                    'BPA',   "#009E73", 'Contaminant or stressor',
                    'TUNIC', "#009E73", 'Contaminant or stressor',
                    'ETOH',  "#009E73", 'Contaminant or stressor',
                    'H20',   "#000000", 'Vehicle control'
                )

# save this for spot checking
d_all <- mergeDecompDiffScores(a_dir, diffpatt=d_patt)
d_all <- add_ihs(d_all, ihs_fname)
d_all <- add_fst(d_all, fst_fname)
write_tsv(d_all, file.path(outdir, paste0(outbase, ".txt.gz")))

# run hallmark and GO:BP GSEA for each score type ranked by gsea_stat
allgenesets <- msigdbr("Homo sapiens")
gsea_res <- aggregate_gsea(d_all, allgenesets=allgenesets, remove_antisense=TRUE)

grid::grid.draw(gsea_res$plot_obj)
ggsave(file.path(outdir, paste0(outbase, ".fgseaTopPathways.png")), plot=gsea_res$plot_obj, width=21, height=25)

write_tsv(gsea_res$res_df, file.path(outdir , paste0(outbase, ".fgseaHallmarkGOBP.txt.gz")))


# get subset of E-G pairs that have 2 diff enhancers per gene for reinforcing
# direction sign and SNP selection stat enrichment tests
ap <- addEnhancerDirectionMatching(d_all,
                                  dist_thresh=dist_thresh,
                                  p_thresh=0.01,
                                  abc_p_col='p',
                                  abc_lfc_col='log2FoldChange',
                                  gene_col='TargetGene',
                                  nEnh=2,
                                  score_type='Score',
                                  cdt1 = 'AFR',
                                  cdt2 = 'EUR',
                                  allscores = TRUE,
                                  abcrank = TRUE,
                                  debug=FALSE)

print("Wilcox dfs #1")
# format for facetted plotting by enhancer level, gene level,
# enhancer status (diff vs. non-diff, reinforcing vs. opposing), etc.
wp1 <- snp_selection_diff_score_wilcox(d_all,
                                       diff_enh_pairs=FALSE,
                                       snp_stat='iHS',
                                       dabc_thresh=0.05,
                                       dabc_col='p',
                                       lfc_col='log2FoldChange',
                                       perc_thresh=0.95,
                                       perc_col='rank2_percentile',
                                       rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                       by_gene=TRUE, # count by gene number instead of enhancer number
                                       all_or_none=FALSE,
                                       run_wilcox=FALSE,
                                       wilcox_col='max_abs_ihs',
                                       gene_col='TargetGene',
                                       allscores=TRUE,
                                       score_col=NULL)

wp2 <- snp_selection_diff_score_wilcox(d_all,
                                       diff_enh_pairs=FALSE,
                                       snp_stat='iHS',
                                       dabc_thresh=0.05,
                                       dabc_col='p',
                                       lfc_col='log2FoldChange',
                                       perc_thresh=0.95,
                                       perc_col='rank2_percentile',
                                       rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                       by_gene=FALSE, # count by gene number instead of enhancer number
                                       all_or_none=FALSE,
                                       run_wilcox=FALSE,
                                       wilcox_col='max_abs_ihs',
                                       gene_col='TargetGene',
                                       allscores=TRUE,
                                       score_col=NULL)

wp3 <- snp_selection_diff_score_wilcox(ap,
                                       diff_enh_pairs=TRUE,
                                       snp_stat='iHS',
                                       dabc_thresh=0.01,
                                       dabc_col='abc_p',
                                       lfc_col='abc_lfc',
                                       perc_thresh=0.95,
                                       perc_col='rank2_percentile',
                                       rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                       by_gene=TRUE, # count by gene number instead of enhancer number
                                       all_or_none=FALSE,
                                       run_wilcox=FALSE,
                                       wilcox_col='max_abs_ihs',
                                       gene_col='gene',
                                       allscores=TRUE,
                                       score_col=NULL,
                                       debug=FALSE)

wp4 <- snp_selection_diff_score_wilcox(ap,
                                       diff_enh_pairs=TRUE,
                                       snp_stat='iHS',
                                       dabc_thresh=0.01,
                                       dabc_col='abc_p',
                                       lfc_col='abc_lfc',
                                       perc_thresh=0.95,
                                       perc_col='rank2_percentile',
                                       rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                       by_gene=FALSE, # count by gene number instead of enhancer number
                                       all_or_none=FALSE,
                                       run_wilcox=FALSE,
                                       wilcox_col='max_abs_ihs',
                                       gene_col='gene',
                                       allscores=TRUE,
                                       score_col=NULL,
                                       debug=FALSE)

print("Fisher dfs #1")

f1 <- snp_selection_diff_score_fisher(d_all,
                                      diff_enh_pairs=FALSE,
                                      snp_stat='iHS',
                                      dabc_thresh=0.05,
                                      dabc_col='p',
                                      lfc_col='log2FoldChange',
                                      perc_thresh=0.95,
                                      mid_filt=TRUE,
                                      alt='greater',
                                      perc_col='rank2_percentile',
                                      rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                      by_gene=TRUE, # count by gene number instead of enhancer number
                                      all_or_none=FALSE,
                                      gene_col='TargetGene',
                                      allscores=TRUE,
                                      score_col=NULL)

f2 <- snp_selection_diff_score_fisher(d_all,
                                      diff_enh_pairs=FALSE,
                                      snp_stat='iHS',
                                      dabc_thresh=0.05,
                                      dabc_col='p',
                                      lfc_col='log2FoldChange',
                                      perc_thresh=0.95,
                                      mid_filt=TRUE,
                                      alt='greater',
                                      perc_col='rank2_percentile',
                                      rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                      by_gene=FALSE, # count by gene number instead of enhancer number
                                      all_or_none=FALSE,
                                      gene_col='TargetGene',
                                      allscores=TRUE,
                                      score_col=NULL)

f3 <- snp_selection_diff_score_fisher(ap,
                                      diff_enh_pairs=TRUE,
                                      snp_stat='iHS',
                                      dabc_thresh=0.01,
                                      dabc_col='abc_p',
                                      lfc_col='abc_lfc',
                                      perc_thresh=0.95,
                                      mid_filt=TRUE,
                                      alt='two.sided',
                                      perc_col='rank2_percentile',
                                      rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                      by_gene=TRUE, # count by gene number instead of enhancer number
                                      all_or_none=FALSE,
                                      gene_col='gene',
                                      allscores=TRUE,
                                      score_col=NULL)

f4 <- snp_selection_diff_score_fisher(ap,
                                      diff_enh_pairs=TRUE,
                                      snp_stat='iHS',
                                      dabc_thresh=0.01,
                                      dabc_col='abc_p',
                                      lfc_col='abc_lfc',
                                      perc_thresh=0.95,
                                      mid_filt=TRUE,
                                      alt='two.sided',
                                      perc_col='rank2_percentile',
                                      rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                      by_gene=FALSE, # count by gene number instead of enhancer number
                                      all_or_none=FALSE,
                                      gene_col='gene',
                                      allscores=TRUE,
                                      score_col=NULL)

print("Wilcox dfs #2")
# format for facetted plotting by enhancer level, gene level,
# enhancer status (diff vs. non-diff, reinforcing vs. opposing), etc.
wp5 <- snp_selection_diff_score_wilcox(d_all,
                                       diff_enh_pairs=FALSE,
                                       snp_stat='FST',
                                       dabc_thresh=0.05,
                                       dabc_col='p',
                                       lfc_col='log2FoldChange',
                                       perc_thresh=0.95,
                                       perc_col='fst_percentile',
                                       rank_col='FST', # column for taking the maximum per enhancer
                                       by_gene=TRUE, # count by gene number instead of enhancer number
                                       all_or_none=FALSE,
                                       run_wilcox=FALSE,
                                       wilcox_col='FST',
                                       gene_col='TargetGene',
                                       allscores=TRUE,
                                       score_col=NULL)

wp6 <- snp_selection_diff_score_wilcox(d_all,
                                       diff_enh_pairs=FALSE,
                                       snp_stat='FST',
                                       dabc_thresh=0.05,
                                       dabc_col='p',
                                       lfc_col='log2FoldChange',
                                       perc_thresh=0.95,
                                       perc_col='fst_percentile',
                                       rank_col='FST', # column for taking the maximum per enhancer
                                       by_gene=FALSE, # count by gene number instead of enhancer number
                                       all_or_none=FALSE,
                                       run_wilcox=FALSE,
                                       wilcox_col='FST',
                                       gene_col='TargetGene',
                                       allscores=TRUE,
                                       score_col=NULL)

wp7 <- snp_selection_diff_score_wilcox(ap,
                                       diff_enh_pairs=TRUE,
                                       snp_stat='FST',
                                       dabc_thresh=0.01,
                                       dabc_col='abc_p',
                                       lfc_col='abc_lfc',
                                       perc_thresh=0.95,
                                       perc_col='fst_percentile',
                                       rank_col='FST', # column for taking the maximum per enhancer
                                       by_gene=TRUE, # count by gene number instead of enhancer number
                                       all_or_none=FALSE,
                                       run_wilcox=FALSE,
                                       wilcox_col='FST',
                                       gene_col='gene',
                                       allscores=TRUE,
                                       score_col=NULL,
                                       debug=FALSE)

wp8 <- snp_selection_diff_score_wilcox(ap,
                                       diff_enh_pairs=TRUE,
                                       snp_stat='FST',
                                       dabc_thresh=0.01,
                                       dabc_col='abc_p',
                                       lfc_col='abc_lfc',
                                       perc_thresh=0.95,
                                       perc_col='fst_percentile',
                                       rank_col='FST', # column for taking the maximum per enhancer
                                       by_gene=FALSE, # count by gene number instead of enhancer number
                                       all_or_none=FALSE,
                                       run_wilcox=FALSE,
                                       wilcox_col='FST',
                                       gene_col='gene',
                                       allscores=TRUE,
                                       score_col=NULL,
                                       debug=FALSE)

print("Fisher dfs #2")

f5 <- snp_selection_diff_score_fisher(d_all,
                                      diff_enh_pairs=FALSE,
                                      snp_stat='FST',
                                      dabc_thresh=0.05,
                                      dabc_col='p',
                                      lfc_col='log2FoldChange',
                                      perc_thresh=0.95,
                                      mid_filt=TRUE,
                                      alt='greater',
                                      perc_col='fst_percentile',
                                      rank_col='FST', # column for taking the maximum per enhancer
                                      by_gene=TRUE, # count by gene number instead of enhancer number
                                      all_or_none=FALSE,
                                      gene_col='TargetGene',
                                      allscores=TRUE,
                                      score_col=NULL)

f6 <- snp_selection_diff_score_fisher(d_all,
                                      diff_enh_pairs=FALSE,
                                      snp_stat='FST',
                                      dabc_thresh=0.05,
                                      dabc_col='p',
                                      lfc_col='log2FoldChange',
                                      perc_thresh=0.95,
                                      mid_filt=TRUE,
                                      alt='greater',
                                      perc_col='fst_percentile',
                                      rank_col='FST', # column for taking the maximum per enhancer
                                      by_gene=FALSE, # count by gene number instead of enhancer number
                                      all_or_none=FALSE,
                                      gene_col='TargetGene',
                                      allscores=TRUE,
                                      score_col=NULL)

f7 <- snp_selection_diff_score_fisher(ap,
                                      diff_enh_pairs=TRUE,
                                      snp_stat='FST',
                                      dabc_thresh=0.01,
                                      dabc_col='abc_p',
                                      lfc_col='abc_lfc',
                                      perc_thresh=0.95,
                                      mid_filt=TRUE,
                                      alt='two.sided',
                                      perc_col='fst_percentile',
                                      rank_col='FST', # column for taking the maximum per enhancer
                                      by_gene=TRUE, # count by gene number instead of enhancer number
                                      all_or_none=FALSE,
                                      gene_col='gene',
                                      allscores=TRUE,
                                      score_col=NULL)

f8 <- snp_selection_diff_score_fisher(ap,
                                      diff_enh_pairs=TRUE,
                                      snp_stat='FST',
                                      dabc_thresh=0.01,
                                      dabc_col='abc_p',
                                      lfc_col='abc_lfc',
                                      perc_thresh=0.95,
                                      mid_filt=TRUE,
                                      alt='two.sided',
                                      perc_col='fst_percentile',
                                      rank_col='FST', # column for taking the maximum per enhancer
                                      by_gene=FALSE, # count by gene number instead of enhancer number
                                      all_or_none=FALSE,
                                      gene_col='gene',
                                      allscores=TRUE,
                                      score_col=NULL)

print("Binding fisher dfs...")

aggf <- rbind(f1,f2,f3,f4,f5,f6,f7,f8)

print("Binding wilcox dfs...")
aggwp <- rbind(wp1, wp2, wp3, wp4, wp5, wp6, wp7, wp8)
write_tsv(aggwp, file.path(outdir, paste0(outbase, ".selectionCandidatesBackground.txt.gz")))


fe <- readTopPairFisherEnrichments(a_dir, fe_patt='fisherEnrichments.txt$')
fer <- readTopPairFisherEnrichments(a_dir, fe_patt='fisherEnrichmentsDirection.txt$')
topfer <- readTopPairFisherEnrichments(a_dir, fe_patt='fisherEnrichmentsDirectionTopN.txt$')

ed <- mergeEnhReinfFisherEnrichments(a_dir, rpatt=r_patt)
edd <- aggregate_enhdirtest(ap, dist_thresh=10000)

write_tsv(fe, file.path(outdir, paste0(outbase, ".DEoverlap.fisherEnrichments.txt")))
write_tsv(fer, file.path(outdir, paste0(outbase, ".DEoverlap.fisherEnrichmentsDirection.txt")))
write_tsv(topfer, file.path(outdir, paste0(outbase, ".DEoverlap.fisherEnrichmentsDirectionTopN.txt")))
write_tsv(ed, file.path(outdir, paste0(outbase, ".enhReinforcingDirection.fisherEnrichments.txt")))
write_tsv(aggf, file.path(outdir, paste0(outbase, ".selectionSNPs.fisherEnrichments.txt")))

if (!is.null(plotdir)) {

    source(file.path(plotdir,"plotting.R"))

    ####### DIFF-SCORE MA PLOTS #######

    p1 <- facetAbcPlotMA(d_all)
    p1
    ggsave(file.path(outdir, paste0(outbase, ".componentScoreMAplot.png")), width=33.75, height=15)


    ####### DE GENE ENRICHMENT PLOTS #######

    dim1 <- c(15,10)
    dim2 <- c(15,10)
    combdims <- c(dim1[1] + dim2[1], dim1[2] + dim2[2])

    g2 <- envCellTypeGridPlot(fe %>%
                                  dplyr::filter(is.na(self_promoter),
                                                doublecount_enh),
                        envPlotColors,
                        p1filt='bulk',
                        p2filt='pseudobulk',
                        gtitle="Diff-score - DE gene enrichments",
                        dim1=dim1,
                        dim2=dim2)
    grid::grid.draw(g2)
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.DEgene.LeaRandolph.png")), plot=g2, width=combdims[1], height=combdims[2])

    p4 <- plotFisherEnrichments(fe %>%
                                    dplyr::filter(datatype %in% c('nzmedian','nzfraction'),
                                                  is.na(self_promoter),
                                                  doublecount_enh),
                            yvar_col = "celltype",
                            ylabel = "Cell type", althypoth = "greater",
                            groupvar = "condition",
                            colorvals = c('darkred','darkgray'), w = 25, h = 10,
                            facets = c('score_type','datatype'), frows = 2,
                            fscales = "free",
                            fdir='v',
                            legendpos=c(.97,.19))
    p4
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.DEgene.scDecomp.png")), width=25, height=10)


    ####### DE DIRECTION ENRICHMENT PLOTS #######

    g3 <- envCellTypeGridPlot(fer %>%
                                  dplyr::filter(is.na(self_promoter),
                                                doublecount_enh),
                        envPlotColors,
                        p1filt='bulk',
                        p2filt='pseudobulk',
                        gtitle="Diff-score - DE direction enrichments",
                        dim1=dim1,
                        dim2=dim2)
    grid::grid.draw(g3)
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.DEdirection.LeaRandolph.png")), plot=g3, width=combdims[1], height=combdims[2])

    p5 <- plotFisherEnrichments(fer %>%
                                dplyr::filter(datatype %in% c('nzmedian','nzfraction'),
                                              is.na(self_promoter),
                                              doublecount_enh),
                            yvar_col = "celltype",
                            ylabel = "Cell type", althypoth = "greater",
                            groupvar = "condition",
                            colorvals = c('darkred','darkgray'), w = 25, h = 10,
                            facets = c('score_type','datatype'), frows = 2,
                            fscales = "free",
                            fdir='v',
                            legendpos=c(.97,.19))
    p5
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.DEdirection.scDecomp.png")), width=25, height=10)

    ####### TOP ABC - DE DIRECTION ENRICHMENT PLOTS #######
    # to compare with enhancer reinforcing direction DE matching results

    g4 <- envCellTypeGridPlot(topfer %>%
                                  dplyr::filter(is.na(self_promoter),
                                                doublecount_enh),
                        envPlotColors,
                        p1filt='bulk',
                        p2filt='pseudobulk',
                        gtitle="Top diff-score - DE direction enrichments",
                        dim1=dim1,
                        dim2=dim2)
    grid::grid.draw(g4)
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.DEdirection.LeaRandolph.topDiffScore.png")), plot=g4, width=combdims[1], height=combdims[2])

    p6 <- plotFisherEnrichments(topfer %>%
                                    dplyr::filter(datatype %in% c('nzmedian','nzfraction'),
                                                  is.na(self_promoter),
                                                  doublecount_enh),
                            yvar_col = "celltype",
                            ylabel = "Cell type", althypoth = "greater",
                            groupvar = "condition",
                            colorvals = c('darkred','darkgray'), w = 25, h = 10,
                            facets = c('score_type','datatype'), frows = 2,
                            fscales = "free",
                            fdir='v',
                            legendpos=c(.97,.19))
    p6
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.DEdirection.scDecomp.topDiffScore.png")), width=25, height=10)


    ####### ENHANCER MATCHING DIRECTIONALITY ENRICHMENT PLOTS #######

    p2 <- plotFisherEnrichments(ed %>%
                                    dplyr::filter(datatype=='ABC_overall'),
                                yvar_col = "score_type",
                                althypoth = "greater",
                                groupvar = NULL,
                                groupvarord = NULL,
                                colorvals = NULL,
                                w = 6.5, h = 5,
                                facets = NULL,
                                frows = 2,
                                fscales = "fixed")
    p2
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.enhReinforcingDirection.overall.png")), width=6.5, height=5)

    p9 <- plotFisherEnrichments(edd,
                                yvar_col = "score_type",
                                althypoth = "greater",
                                groupvar = NULL,
                                groupvarord = NULL,
                                colorvals = NULL,
                                w = 6.5, h = 5,
                                facets = NULL,
                                frows = 2,
                                fscales = "fixed")
    p9
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.enhReinforcingDirection.overall.10kb.png")), width=6.5, height=5)

    g1 <- envCellTypeGridPlot(ed,
                        envPlotColors,
                        p1filt='bulk',
                        p2filt='pseudobulk',
                        gtitle="Diff-score reinforcing enhancer - DE directionality matching enrichments",
                        dim1=dim1,
                        dim2=dim2)
    grid::grid.draw(g1)
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.enhReinforcingDirection.LeaRandolphDE.png")), plot=g1, width=combdims[1], height=combdims[2])

    p3 <- plotFisherEnrichments(ed %>%
                                dplyr::filter(datatype %in% c('nzmedian','nzfraction')),
                            yvar_col = "celltype",
                            ylabel = "Cell_type", althypoth = "greater",
                            groupvar = "condition",
                            colorvals = c('darkred','darkgray'), w = 25, h = 10,
                            facets = c('score_type','datatype'), frows = 2,
                            fscales = "free",
                            fdir='v',
                            legendpos=c(.93,.13))
    p3
    ggsave(file.path(outdir, paste0(outbase, ".fisherEnrichments.enhReinforcingDirection.scDecompDE.png")), width=25, height=10)


    ####### DIFF-ABC/ENHANCER MATCHING DIRECTIONALITY - SNP SELECTION STAT ENRICHMENT PLOTS #######

    p7 <- plotFisherEnrichments(aggf,
                      yvar_col = "score_type",
                      althypoth = "two.sided",
                      groupvar = "test_set",
                      groupvarlab = NULL,
                      groupvarord = NULL,
                      colorvals = c('blue','gold'), w = 12, h = 12,
                      facets = c('snp_stat','enh_type'), frows = 2,
                      fscales = "free",
                      legendpos=c(.9,.1))
    p7
    ggsave(file.path(outdir, paste0(outbase, ".selectionSNPs.fisherEnrichments.png")), width=12, height=12)

    p8 <- groupedBoxplot(aggwp %>%
                             dplyr::filter(snp_stat=='iHS') %>%
                             mutate(enh_status = fct_relevel(enh_status,
                                                             c('opposing',
                                                               'reinforcing',
                                                               'non-diff',
                                                               'diff'))),
                         'enh_status',
                         'wilcox_var',
                         groupvar=NULL,
                         xlab='Enhancer status',
                         ylab='Max |iHS|',
                         plotpoints=TRUE,
                         show_wilcox=TRUE,
                         logscale=FALSE,
                         colorX=TRUE,
                         colorvals=c('blue','gold','blue','gold'),
                         showleg=FALSE,
                         legendpos='right',
                         horizontal=FALSE,
                         facets=c('score_type','enh_type','test_set'),
                         frows=2,
                         fscales='free',
                         fdir='v',
                         w=22.5,h=12,
                         debug=FALSE)
    p8
    ggsave(file.path(outdir, paste0(outbase, ".iHS.Wilcoxon.png")), width=22.5, height=12)

    p10 <- groupedBoxplot(aggwp %>%
                             dplyr::filter(snp_stat=='FST') %>%
                             mutate(enh_status = fct_relevel(enh_status,
                                                             c('opposing',
                                                               'reinforcing',
                                                               'non-diff',
                                                               'diff'))),
                         'enh_status',
                         'wilcox_var',
                         groupvar=NULL,
                         xlab='Enhancer status',
                         ylab='FST',
                         plotpoints=TRUE,
                         show_wilcox=TRUE,
                         logscale=FALSE,
                         colorX=TRUE,
                         colorvals=c('blue','gold','blue','gold'),
                         showleg=FALSE,
                         legendpos='right',
                         horizontal=FALSE,
                         facets=c('score_type','enh_type','test_set'),
                         frows=2,
                         fscales='free',
                         fdir='v',
                         w=22.5,h=12,
                         debug=FALSE)
    p10
    ggsave(file.path(outdir, paste0(outbase, ".FST.Wilcoxon.png")), width=22.5, height=12)

}
