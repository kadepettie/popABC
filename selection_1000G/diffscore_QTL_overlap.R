#!/usr/bin/R

library(argparse)
library(docstring)
library(readxl)
library(plyranges)
library(ggplot2)
library(tidyverse)

diffQtlFishers <- function(df,
                      qtl_p_thresh=0.05,
                      qtl_p_nonthresh=0.5,
                      qtl_p_col='bqtl_p',
                      abc_p_thresh=0.05,
                      abc_p_nonthresh=0.5,
                      abc_p_col='atac_p',
                      class_filt='none', # 'promoter','non-promoter'
                      class_res='atac',
                      mat_out=FALSE,
                      alt='two.sided') {

    dfs <- df %>%
      dplyr::rename(qtl_p := {{qtl_p_col}},
                    abc_p := {{abc_p_col}})



    if (grepl('chip', abc_p_col) | grepl('hic', abc_p_col)) {
      # only count each HiChIP bin once
      dfs <- dfs %>%
            group_by(hgnc_symbol, name) %>%
            arrange(qtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup() %>%
            group_by(hgnc_symbol, chip_lfc, abc_p) %>%
            mutate(chipbin_has_promoter = if_else(any(class=='promoter'), TRUE, FALSE))

      if (class_res=='chip') {

        if (class_filt=='promoter') {
            dfs <- dfs %>% dplyr::filter(chipbin_has_promoter)
        } else if (class_filt=='non-promoter') {
            dfs <- dfs %>% dplyr::filter(!chipbin_has_promoter)
        }

      } else {
        if (class_filt=='promoter') {
            print('Testing promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class=='promoter')
        } else if (class_filt=='non-promoter') {
            print('Testing non-promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class!='promoter')
        } else {
            print('No CRE class type filtering.')
        }
      }

      # take the CRE with the most significant QTL in the HiChIP bin
      dfs <- dfs %>%
            arrange(qtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup()

    } else {
        if (class_filt=='promoter') {
            print('Testing promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class=='promoter')
        } else if (class_filt=='non-promoter') {
            print('Testing non-promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class!='promoter')
        } else {
            print('No CRE class type filtering.')
        }
    }

    countsdf <- dfs %>%
      group_by(name) %>%
      arrange(qtl_p, abc_p, .by_group=TRUE) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(abc_type = case_when(abc_p < abc_p_thresh ~ 'diff-ABC',
                                  abc_p >= abc_p_nonthresh ~ 'non-diff-ABC',
                                  TRUE ~ 'ambig'),
             qtl_type = case_when(qtl_p < qtl_p_thresh ~ 'QTL',
                                  qtl_p >= qtl_p_nonthresh ~ 'non-QTL',
                                  TRUE ~ 'ambig')) %>%
      dplyr::filter(abc_type != 'ambig',
                    qtl_type != 'ambig') %>%
      group_by(abc_type) %>%
      summarize(
        var_count = sum(!is.na(qtl_p)),
        qtl_count =sum(qtl_p < qtl_p_thresh, na.rm=TRUE)
      )

    fg_numerator <- countsdf %>%
        dplyr::filter(abc_type=='diff-ABC') %>%
        pull(qtl_count)

    fg_denominator <- countsdf %>%
        pull(qtl_count) %>%
        sum()

    bg_numerator <- countsdf %>%
        dplyr::filter(abc_type=='diff-ABC') %>%
        pull(var_count)

    bg_denominator <- countsdf %>%
        dplyr::select(var_count) %>%
        sum()

    fe <- (fg_numerator/fg_denominator) / (bg_numerator/bg_denominator)

    mat11 <- fg_numerator
    mat21 <- fg_denominator-fg_numerator
    mat12 <- bg_numerator-fg_numerator
    mat22 <- bg_denominator-mat11-mat21-mat12
    tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                      dimnames =
               list(c('diff-ABC', 'non-diff-ABC'),
                    c('QTL', 'non-QTL')))

    ft <- fisher.test(tmat, alternative=alt)

    if (mat_out) {
        return(
            list(
                tmat,
                ft
            )
        )
    }

    eres <- list(qtl_significance = qtl_p_thresh,
                 diff_significance = abc_p_thresh,
                 odds_ratio = ft$estimate[[1]],
                 pvalue = ft$p.value,
                 conf_lower = ft$conf.int[1],
                 conf_upper = ft$conf.int[2],
                 up_up = mat11,
                 up_down = mat21,
                 down_up = mat12,
                 down_down = mat22)

    return(eres)

}

diffQTLdirectionFishers <- function(df,
                                      qtl_p_thresh=0.05,
                                      qtl_p_col='bqtl_p',
                                      qtl_dir_col='bqtl_direction',
                                      abc_p_thresh=0.05,
                                      abc_p_col='atac_p',
                                      abc_dir_col='atac_monodirection',
                                      fst_perc_thresh=0.0, # [0.0-1.0]
                                      fst_perc_col='fst_percentile',
                                      class_filt='none', # 'promoter', 'non-promoter', 'none'
                                      class_res='atac',
                                      mat_out=FALSE,
                                      alt='two.sided') {

    dfs <- df %>%
        dplyr::rename(qtl_p := {{qtl_p_col}},
                      abc_p := {{abc_p_col}},
                      qtl_dir := {{qtl_dir_col}},
                      abc_dir := {{abc_dir_col}}) %>%
        dplyr::filter(qtl_p < qtl_p_thresh,
                      abc_p < abc_p_thresh,
                      !is.na(qtl_dir))

    if (grepl('chip', abc_p_col) | grepl('hic', abc_p_col)) {
      # only count each HiChIP bin once
      dfs <- dfs %>%
            group_by(hgnc_symbol, name) %>%
            arrange(qtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup() %>%
            group_by(hgnc_symbol, chip_lfc, abc_p) %>%
            mutate(chipbin_has_promoter = if_else(any(class=='promoter'), TRUE, FALSE))

      if (class_res=='chip') {

        if (class_filt=='promoter') {
            dfs <- dfs %>% dplyr::filter(chipbin_has_promoter)
        } else if (class_filt=='non-promoter') {
            dfs <- dfs %>% dplyr::filter(!chipbin_has_promoter)
        }

      } else {
        if (class_filt=='promoter') {
            print('Testing promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class=='promoter')
        } else if (class_filt=='non-promoter') {
            print('Testing non-promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class!='promoter')
        } else {
            print('No CRE class type filtering.')
        }
      }

      # take the CRE with the most significant QTL in the HiChIP bin
      dfs <- dfs %>%
            arrange(qtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup()

    } else {
        if (class_filt=='promoter') {
            print('Testing promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class=='promoter')
        } else if (class_filt=='non-promoter') {
            print('Testing non-promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class!='promoter')
        } else {
            print('No CRE class type filtering.')
        }
    }

    if (!is.null(fst_perc_thresh)) {
        dfs <- dfs %>%
            dplyr::rename(fst_perc := {{fst_perc_col}}) %>%
            dplyr::filter(fst_perc > fst_perc_thresh)
    }

    dfs <- dfs %>%
      group_by(name) %>%
      arrange(qtl_p, abc_p, .by_group=TRUE) %>%
      slice(1) %>%
      ungroup()

    # # could try majortiy of bQTL direction to match,
    # but for now just take top
    de_up_da_up <- dfs %>%
        dplyr::filter(qtl_dir=='up',
               abc_dir=='up') %>%
        nrow()

    de_up_da_down <- dfs %>%
        dplyr::filter(qtl_dir=='up',
               abc_dir=='down') %>%
        nrow()

    de_down_da_up <- dfs %>%
        dplyr::filter(qtl_dir=='down',
               abc_dir=='up') %>%
        nrow()

    de_down_da_down <- dfs %>%
        dplyr::filter(qtl_dir=='down',
               abc_dir=='down') %>%
        nrow()

    mat11 <- de_up_da_up
    mat21 <- de_up_da_down
    mat12 <- de_down_da_up
    mat22 <- de_down_da_down
    tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                      dimnames =
               list(c('DA up', 'DA down'),
                    c('QTL up', 'QTL down')))

    ft <- fisher.test(tmat, alternative=alt)

    if (mat_out) {
        return(
            list(
                tmat,
                ft
            )
        )
    }

    eres <- list(qtl_significance = qtl_p_thresh,
                 diff_significance = abc_p_thresh,
                 odds_ratio = ft$estimate[[1]],
                 pvalue = ft$p.value,
                 conf_lower = ft$conf.int[1],
                 conf_upper = ft$conf.int[2],
                 up_up = mat11,
                 up_down = mat21,
                 down_up = mat12,
                 down_down = mat22)

    return(eres)

}

DEdiffQTLdirectionFishers <- function(df,
                                      qtl_p_thresh=0.05,
                                      qtl_p_col='bqtl_p',
                                      qtl_dir_col='bqtl_direction',
                                      abc_p_thresh=0.05,
                                      abc_p_col='atac_p',
                                      abc_dir_col='atac_monodirection',
                                      de_dir_col='direction_de_lcl',
                                      fst_perc_thresh=0.0, # [0.0-1.0]
                                      fst_perc_col='bqtl_fst_percentile',
                                      class_res='atac',
                                      class_filt='none', # 'promoter', 'non-promoter', 'none'
                                      self_promoter=NULL,
                                      mat_out=FALSE,
                                      alt='two.sided',
                                      debug=FALSE) {

    dfs <- df %>%
        dplyr::rename(qtl_p := {{qtl_p_col}},
                      abc_p := {{abc_p_col}},
                      qtl_dir := {{qtl_dir_col}},
                      abc_dir := {{abc_dir_col}},
                      de_dir := {{de_dir_col}}) %>%
        dplyr::filter(qtl_p < qtl_p_thresh,
                      abc_p < abc_p_thresh,
                      !is.na(qtl_dir),
                      !is.na(de_dir))

    if (grepl('chip', abc_p_col) | grepl('hic', abc_p_col)) {
      # only count each HiChIP bin once
      dfs <- dfs %>%
            group_by(hgnc_symbol, name) %>%
            arrange(qtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup() %>%
            group_by(hgnc_symbol, chip_lfc, abc_p) %>%
            mutate(chipbin_has_promoter = if_else(any(class=='promoter'), TRUE, FALSE))

      if (class_res=='chip') {

        if (class_filt=='promoter') {
            dfs <- dfs %>% dplyr::filter(chipbin_has_promoter)
        } else if (class_filt=='non-promoter') {
            dfs <- dfs %>% dplyr::filter(!chipbin_has_promoter)
        }

      } else {
        if (class_filt=='promoter') {
            print('Testing promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class=='promoter')
        } else if (class_filt=='non-promoter') {
            print('Testing non-promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class!='promoter')
        } else {
            print('No CRE class type filtering.')
        }
      }

      # take the CRE with the most significant QTL in the HiChIP bin
      dfs <- dfs %>%
            arrange(qtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup()

    } else {
        if (class_filt=='promoter') {
            print('Testing promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class=='promoter')
        } else if (class_filt=='non-promoter') {
            print('Testing non-promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class!='promoter')
        } else {
            print('No CRE class type filtering.')
        }
    }

    if (!is.null(fst_perc_thresh)) {
        dfs <- dfs %>%
            dplyr::rename(fst_perc := {{fst_perc_col}}) %>%
            dplyr::filter(fst_perc > fst_perc_thresh)
    }

    dfs <- dfs %>%
      group_by(hgnc_symbol) %>%
      arrange(abc_p, qtl_p, .by_group=TRUE) %>%
      slice(1) %>%
      ungroup()

    dfs <- dfs %>%
        group_by(name) %>%
        mutate(deqtl_abc_dir = case_when(
                    all(de_dir=='up') & (qtl_dir=='up' & abc_dir=='up') ~ 'up_up',
                    all(de_dir=='up') & (qtl_dir=='down' & abc_dir=='down') ~ 'up_down',
                    all(de_dir=='down') & (qtl_dir=='up' & abc_dir=='up') ~ 'down_up',
                    all(de_dir=='down') & (qtl_dir=='down' & abc_dir=='down') ~ 'down_down'
                )
        ) %>%
        dplyr::filter(!is.na(deqtl_abc_dir)) %>%
        mutate(includes_selfpromoter = if_else(any(isSelfPromoter), TRUE, FALSE)) %>%
        slice(1) %>%
        ungroup()

    if (!is.null(self_promoter)) {
        if (self_promoter) {
            dfs <- dfs %>% dplyr::filter(includes_selfpromoter)
        } else {
            dfs <- dfs %>% dplyr::filter(!includes_selfpromoter)
        }
    }

    if (debug) return(dfs)

    # # could try majortiy of bQTL direction to match,
    # but for now just take top
    de_up_da_up <- dfs %>%
        dplyr::filter(deqtl_abc_dir=='up_up') %>%
        nrow()

    de_up_da_down <- dfs %>%
        dplyr::filter(deqtl_abc_dir=='up_down') %>%
        nrow()

    de_down_da_up <- dfs %>%
        dplyr::filter(deqtl_abc_dir=='down_up') %>%
        nrow()

    de_down_da_down <- dfs %>%
        dplyr::filter(deqtl_abc_dir=='down_down') %>%
        nrow()

    mat11 <- de_up_da_up
    mat21 <- de_up_da_down
    mat12 <- de_down_da_up
    mat22 <- de_down_da_down
    tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                      dimnames =
               list(c('QTL/DA up', 'QTL/DA down'),
                    c('DE up', 'DE down')))

    ft <- fisher.test(tmat, alternative=alt)

    if (mat_out) {
        return(
            list(
                tmat,
                ft
            )
        )
    }

    eres <- list(qtl_significance = qtl_p_thresh,
                 diff_significance = abc_p_thresh,
                 odds_ratio = ft$estimate[[1]],
                 pvalue = ft$p.value,
                 conf_lower = ft$conf.int[1],
                 conf_upper = ft$conf.int[2],
                 up_up = mat11,
                 up_down = mat21,
                 down_up = mat12,
                 down_down = mat22)

    return(eres)

}

diffQTLsignTest <-  function(df,
                              qtl_p_thresh=0.05,
                              qtl_p_col='bqtl_p',
                              qtl_dir_col='bqtl_direction',
                              abc_p_thresh=0.05,
                              abc_p_col='atac_p',
                              abc_dir_col='atac_monodirection',
                              fst_perc_thresh=0.0, # [0.0-1.0]
                              fst_perc_col='bqtl_fst_percentile',
                              class_filt='none', # 'promoter', 'non-promoter', 'none'
                              class_res='atac',
                              alt='two.sided') {

    dfs <- df %>%
        dplyr::rename(qtl_p := {{qtl_p_col}},
                      abc_p := {{abc_p_col}},
                      qtl_dir := {{qtl_dir_col}},
                      abc_dir := {{abc_dir_col}}) %>%
        # dplyr::filter to QTL with directional info
        dplyr::filter(qtl_p < qtl_p_thresh,
                      !is.na(qtl_dir))

    if (grepl('chip', abc_p_col) | grepl('hic', abc_p_col)) {
      # only count each HiChIP bin once
      dfs <- dfs %>%
            group_by(hgnc_symbol, name) %>%
            arrange(qtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup() %>%
            group_by(hgnc_symbol, chip_lfc, abc_p) %>%
            mutate(chipbin_has_promoter = if_else(any(class=='promoter'), TRUE, FALSE))

      if (class_res=='chip') {

        if (class_filt=='promoter') {
            dfs <- dfs %>% dplyr::filter(chipbin_has_promoter)
        } else if (class_filt=='non-promoter') {
            dfs <- dfs %>% dplyr::filter(!chipbin_has_promoter)
        }

      } else {
        if (class_filt=='promoter') {
            print('Testing promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class=='promoter')
        } else if (class_filt=='non-promoter') {
            print('Testing non-promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class!='promoter')
        } else {
            print('No CRE class type filtering.')
        }
      }

      # take the CRE with the most significant QTL in the HiChIP bin
      dfs <- dfs %>%
            arrange(qtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup()

    } else {
        if (class_filt=='promoter') {
            print('Testing promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class=='promoter')
        } else if (class_filt=='non-promoter') {
            print('Testing non-promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class!='promoter')
        } else {
            print('No CRE class type filtering.')
        }
    }


    if (!is.null(fst_perc_thresh)) {
        dfs <- dfs %>%
            dplyr::rename(fst_perc := {{fst_perc_col}}) %>%
            dplyr::filter(fst_perc > fst_perc_thresh)
    }

    # only count each CRE once
    # (remove double counts due to multiple target genes per CRE)
    dfs <- dfs %>%
      group_by(name) %>%
      arrange(qtl_p, abc_p, .by_group=TRUE) %>%
      slice(1) %>%
      ungroup()

    # probability of success as proportion down bQTL & down DE
    # out of all CREs with matching direction
    dfs <- dfs %>%
        mutate(diffqtl_direction = case_when(qtl_dir=='up' & abc_dir=='up' ~ 'up_up',
                                             qtl_dir=='down' & abc_dir=='down' ~ 'down_down'),
               diff_status = if_else(abc_p < abc_p_thresh, 'diff', 'non-diff')) %>%
        dplyr::filter(!is.na(diffqtl_direction))

    gws <- dfs  %>%
        group_by(diffqtl_direction) %>%
        summarize(N = dplyr::n()) %>%
        mutate(prop = N/sum(N))

    print("Background QTL - diff-score direction matching (all CREs):")
    print(gws)

    diffs <- dfs %>%
        dplyr::filter(diff_status=='diff') %>%
        group_by(diffqtl_direction) %>%
        summarize(N = dplyr::n()) %>%
        mutate(prop = N/sum(N))

    psuccess <- gws[gws$diffqtl_direction=='down_down', "prop", drop=TRUE]
    nsuccess <- diffs[diffs$diffqtl_direction=='down_down', "N", drop=TRUE]
    if (length(nsuccess)==0) nsuccess <- 0
    nfail <- diffs[diffs$diffqtl_direction=='up_up', "N", drop=TRUE]
    if (length(nfail)==0) nfail <- 0
    ntrials <- nsuccess + nfail
    if (length(ntrials)==0) ntrials <- 0

    if (ntrials==0) {

        print("No diff-CRE/QTL to test!")

        bres <- list(qtl_significance = qtl_p_thresh,
                     diff_significance = abc_p_thresh,
                     odds_ratio = NA,
                     pvalue = NA,
                     conf_lower = NA,
                     conf_upper = NA,
                     n_success = nsuccess,
                     n_trials = ntrials,
                     null_prob = psuccess)

    } else {

        bt <- binom.test(nsuccess, ntrials, p=psuccess, alternative=alt)
        print(bt)

        bres <- list(qtl_significance = qtl_p_thresh,
                     diff_significance = abc_p_thresh,
                     odds_ratio = bt$estimate[[1]]/psuccess,
                     pvalue = bt$p.value,
                     conf_lower = bt$conf.int[1]/psuccess,
                     conf_upper = bt$conf.int[2]/psuccess,
                     n_success = nsuccess,
                     n_trials = ntrials,
                     null_prob = psuccess)

    }

    return(bres)

}

diffQTLFSTWilcoxon <-  function(df,
                              qtl_p_thresh=0.05,
                              qtl_p_col='bqtl_p',
                              qtl_dir_col='bqtl_direction',
                              abc_p_thresh=0.05,
                              abc_p_nonthresh=abc_p_thresh,
                              abc_p_col='atac_p',
                              abc_dir_col='atac_monodirection',
                              fst_col='bqtl_fst',
                              class_filt='none', # 'promoter', 'non-promoter', 'none'
                              class_res='atac',
                              conf_int=FALSE,
                              return_df=TRUE) {

    outcols <- c('seqnames',
                 'start',
                 'end',
                 'name',
                 'class',
                 'isSelfPromoter',
                 'qtl_p',
                 'abc_p',
                 'qtl_dir',
                 'abc_dir',
                 'atac_lfc',
                 'chip_lfc',
                 'hic_lfc',
                 'fst_var',
                 'diff_status')

    dfs <- df %>%
        dplyr::rename(qtl_p := {{qtl_p_col}},
                      abc_p := {{abc_p_col}},
                      qtl_dir := {{qtl_dir_col}},
                      abc_dir := {{abc_dir_col}}) %>%
        # dplyr::filter to QTL with directional info
        dplyr::filter(qtl_p < qtl_p_thresh)

    if (grepl('chip', abc_p_col) | grepl('hic', abc_p_col)) {
      # only count each HiChIP bin once
      dfs <- dfs %>%
            group_by(hgnc_symbol, name) %>%
            arrange(qtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup() %>%
            group_by(hgnc_symbol, chip_lfc, abc_p) %>%
            mutate(chipbin_has_promoter = if_else(any(class=='promoter'), TRUE, FALSE))

      if (class_res=='chip') {

        if (class_filt=='promoter') {
            dfs <- dfs %>% dplyr::filter(chipbin_has_promoter)
        } else if (class_filt=='non-promoter') {
            dfs <- dfs %>% dplyr::filter(!chipbin_has_promoter)
        }

      } else {
        if (class_filt=='promoter') {
            print('Testing promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class=='promoter')
        } else if (class_filt=='non-promoter') {
            print('Testing non-promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class!='promoter')
        } else {
            print('No CRE class type filtering.')
        }
      }

      # take the CRE with the most significant QTL in the HiChIP bin
      dfs <- dfs %>%
            arrange(qtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup()

    } else {
        if (class_filt=='promoter') {
            print('Testing promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class=='promoter')
        } else if (class_filt=='non-promoter') {
            print('Testing non-promoter CREs only')
            dfs <- dfs %>% dplyr::filter(class!='promoter')
        } else {
            print('No CRE class type filtering.')
        }
    }


    dfs <- dfs %>%
        dplyr::rename(fst_var := {{fst_col}}) %>%
        dplyr::filter(!is.na(fst_var))

    # only count each CRE once
    # (remove double counts due to multiple target genes per CRE)
    dfs <- dfs %>%
      group_by(name) %>%
      arrange(qtl_p, abc_p, .by_group=TRUE) %>%
      slice(1) %>%
      ungroup()

    # probability of success as proportion down bQTL & down DE
    # out of all CREs with matching direction
    dfs <- dfs %>%
        mutate(diff_status = case_when(abc_p < abc_p_thresh ~ 'diff',
                                       abc_p >= abc_p_nonthresh ~ 'non-diff')) %>%
        dplyr::filter(!is.na(diff_status))

    diff_fst <- dfs %>%
        dplyr::filter(diff_status=='diff') %>%
        pull(fst_var)
    nondiff_fst <- dfs %>%
        dplyr::filter(diff_status=='non-diff') %>%
        pull(fst_var)

    wt <- wilcox.test(diff_fst, nondiff_fst, alternative='greater', conf.int=conf_int)

    wres <- list(qtl_significance = qtl_p_thresh,
                 diff_significance = abc_p_thresh,
                 pvalue = wt$p.value)

    if (return_df) {
        outdf <- dfs %>%
            dplyr::select(any_of(outcols)) %>%
            mutate(diff_wilcoxP = wt$p.value)
        return(outdf)
    }

    return(wres)

}

genesetQtlFishers <- function(df,
                                dscore_thresh=0.05,
                                qtl_sig_col='bqtl_p',
                                qtl_sig=0.005,
                                qtl_nonsig=0.5,
                                mat_out=FALSE,
                                alt='greater',
                                debug=FALSE) {

    dfs <- df %>%
        dplyr::rename(qtl_p := {{qtl_sig_col}}) %>%
        dplyr::filter(p < dscore_thresh,
                      qtl_p < qtl_sig | qtl_p >= qtl_nonsig) %>%
        group_by(name) %>%
        mutate(direction = case_when(any(gene_set != 'other') & qtl_p < qtl_sig ~ 'up_up',
                                     any(gene_set != 'other') & qtl_p >= qtl_nonsig ~ 'up_down',
                                     all(gene_set == 'other') & qtl_p < qtl_sig ~ 'down_up',
                                     all(gene_set == 'other') & qtl_p >= qtl_nonsig ~ 'down_down'))

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
               list(c('bQTL enhancer', 'non-bQTL enhancer'),
                    c('target gene in set', 'no target gene in set')))

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

liftover_grange <- function(df, ch_fname, from_b="GRCh38", to_b="hg19") {

    #' Liftover granges object
    #'
    #' @param df Genomic ranges object to liftover
    #' @param ch_fname Filename of liftover chain file (must correspond to `from_b`/`to_b`)
    #' @param from_b Genome of input granges object
    #' @param to_b Genome to liftover to

    genome(df) = from_b
    if(!(require(rtracklayer))) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

      BiocManager::install("rtracklayer")
    }
    library(rtracklayer)
    ch = import.chain(ch_fname)
    seqlevelsStyle(df) = "UCSC"  # necessary
    print(paste0("Lifting over ranges from ", from_b, " to ", to_b, "..."))
    df_new = liftOver(df, ch)
    df_new = unlist(df_new)
    genome(df_new) = to_b
    nc <- length(df_new) - length(df)
    print(paste0("Ranges gained/lost in liftOver = ", if_else(nc>0,'+',''), nc))
    return(df_new)
}

add_sharing_status <- function(df,
                               es_pre='beta',
                               pv_pre='lfsr',
                               gene_col='hgnc_symbol',
                               var_col='eqtl_pos',
                               eqtl=FALSE,
                               sharing_def='es',
                               negate_beta=FALSE,
                               pivot=TRUE,
                               pivotw=FALSE) {

    if (negate_beta) {
        df <- df %>%
            mutate(across(dplyr::starts_with(es_pre), ~ -.x))
    }

    if (pivot) {
        df <- df %>%
            pivot_longer(cols=matches(paste0(es_pre,'_|',pv_pre,'_')),
                     names_to=c(".value","condition","celltype","datatype"),
                     names_pattern="(.*)_(.*)_(.*)_(.*)")
    }

    df <- df %>%
        dplyr::rename(gene := {{gene_col}},
                      beta := {{es_pre}},
                      pval := {{pv_pre}})

    if (sharing_def=='es') {
        df <- df %>%
            group_by(across(any_of(c('gene',var_col)))) %>%
            mutate(
                type_de = case_when(
                    any(pval < 0.1) &
                        all(dplyr::between(beta[pval!=min(pval)],
                                           if_else(median(beta[pval==min(pval)])<0,median(beta[pval==min(pval)])*2,median(beta[pval==min(pval)])/2),
                                           if_else(median(beta[pval==min(pval)])<0,median(beta[pval==min(pval)])/2,median(beta[pval==min(pval)])*2)))
                            ~ 'ubiquitous',
                    any(pval < 0.1) &
                        any(dplyr::between(beta[pval!=min(pval)],
                                           if_else(median(beta[pval==min(pval)])<0,median(beta[pval==min(pval)])*2,median(beta[pval==min(pval)])/2),
                                           if_else(median(beta[pval==min(pval)])<0,median(beta[pval==min(pval)])/2,median(beta[pval==min(pval)])*2)))
                            ~ 'context-dependent',
                    any(pval < 0.1) ~ 'context-specific'

                )
            )
    } else if (sharing_def=='lfsr') {
        df <- df %>%
            group_by(across(any_of(c('gene',var_col)))) %>%
            mutate(
                type_de = case_when(
                    all(pval < 0.1) & (all(beta > 0) | all(beta < 0)) ~ 'ubiquitous',
                    sum(pval < 0.1) == 1 ~ 'context-specific',
                    any(pval < 0.1) ~ 'context-dependent'
                )
            )
    } else {
        stop('Sharing def not recognized')
    }

    df <- df %>%
        mutate(
            direction_de = case_when(
                any(pval < 0.1) & median(beta[pval==min(pval)]) < 0 ~ 'down',
                any(pval < 0.1) & median(beta[pval==min(pval)]) > 0 ~ 'up'
            )
        ) %>%
        ungroup() %>%
        dplyr::rename(!!gene_col := gene,
                      !!es_pre := beta,
                      !!pv_pre := pval)

    if (eqtl) df <- df %>% dplyr::rename(type_eqtl = type_de,
                                         direction_eqtl = direction_de)

    if (pivotw) {
        df <- df %>%
            pivot_wider(names_from=c(condition,celltype,datatype),
                names_glue="{.value}_{condition}_{celltype}_{datatype}",
                values_from=c(es_pre, pv_pre))
    }

    return(df)

}

add_hgnc <- function (df, eGene = FALSE, biomart = FALSE, trans_eGene=FALSE) {
    if (eGene) {
        ensgenes <- df %>% pull(eGene) %>% unique()
        ensgenes <- ensgenes[which(!is.na(ensgenes))]
        merge_col <- "eGene"
    } else if (trans_eGene) {
        ensgenes <- df %>% pull(trans_eGene) %>% unique()
        ensgenes <- ensgenes[which(!is.na(ensgenes))]
        merge_col <- "trans_eGene"
    }
    else {
        ensgenes <- df %>% pull(ensembl_gene_id) %>% unique()
        ensgenes <- ensgenes[which(!is.na(ensgenes))]
        merge_col <- "ensembl_gene_id"
    }
    if (biomart) {
        if(!(require(biomaRt))) {
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
          BiocManager::install("biomaRt")
        }
        library(biomaRt)
        ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
        mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
        gene_attributes = c("ensembl_gene_id", "hgnc_symbol")
        print("Getting hgnc_symbols from biomaRt...")
        gene_info = getBM(attributes = gene_attributes, filters = "ensembl_gene_id",
            values = ensgenes, mart = ensembl)
    }
    else {
        if(!(require(EnsDb.Hsapiens.v79))) {
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")

          BiocManager::install("EnsDb.Hsapiens.v79")
        }
        library(EnsDb.Hsapiens.v79)
        gene_info <- ensembldb::select(EnsDb.Hsapiens.v79,
                                       keys = ensgenes,
                                       keytype = "GENEID",
                                       columns = c("SYMBOL", "GENEID")) %>%
            as_tibble() %>%
            dplyr::rename(ensembl_gene_id = GENEID, hgnc_symbol = SYMBOL)

        # unload namespaces due to dplyr masking (dplyr::filter, etc.)
        # unloadNamespace('EnsDb.Hsapiens.v79')
        # unloadNamespace('ensembldb')
    }
    if (eGene) {
        gene_info <- gene_info %>% as_tibble() %>% dplyr::rename(eGene = ensembl_gene_id,
            eGene_hgnc = hgnc_symbol)
    }
    if (trans_eGene) {
        gene_info <- gene_info %>% as_tibble() %>% dplyr::rename(trans_eGene = ensembl_gene_id,
            trans_eGene_hgnc = hgnc_symbol)
    }

    print("Adding hgnc_symbols...")
    df <- merge(df, gene_info, by = merge_col, all.x = TRUE)

    return(df)
}

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
        dplyr::filter(grepl('ENSG', ensembl_gene_id)) # dplyr::filter out alternate gene symbols

    print("Adding ensembl_gene_ids...")
    df <- merge(df, gene_info, by = merge_col, all.x = TRUE)

    # unload namespaces due to dplyr masking (dplyr::filter, etc.)
    # unloadNamespace('EnsDb.Hsapiens.v79')
    # unloadNamespace('ensembldb')

    return(df)

}

prepABCforGSEA <- function(df, allscores=FALSE) {
    # converts TargetGenes (hgnc) to unique ensIDs

    # some hgcn's have multiple ensg's but usually lower ID number is in ENSG
    # keep just one
    df <- df %>%
        dplyr::rename(hgnc_symbol = TargetGene) %>%
        hgnc2ensembl(.)

    if (allscores) df <- df %>% group_by(score_type)

    df <- df %>%
        group_by(hgnc_symbol,name, .add=TRUE) %>%
        arrange(ensembl_gene_id, entrez_id, .by_group=TRUE) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        arrange(p)

    return(df)

}

widen_diffABC <- function(df) {

    dabc <- prepABCforGSEA(df, allscores=TRUE) %>%
        dplyr::select(-c(group1,group2,n1,n2)) %>%
        mutate(neglog10p_ABC = -log10(p)) %>%
        dplyr::rename(log2FoldChange_ABC = log2FoldChange,
                      p_ABC = p)

    edircols <- c('mean_ABC','enh_idx','enh_type','direction','binary_direction','enh_distance')

    chipatac <- dabc %>%
        dplyr::filter(score_type %in% c('chip','atac','hic')) %>%
        dplyr::rename(lfc = log2FoldChange_ABC,
                      p = p_ABC,
                      nltp = neglog10p_ABC) %>%
        dplyr::select(ensembl_gene_id,
                      entrez_id,
                      hgnc_symbol,
                      TargetGeneTSS,
                      seqnames,
                      start,
                      end,
                      class,
                      name,
                      distance,
                      isSelfPromoter,
                      max_abs_ihs,
                      rank2_percentile,
                      ihs_pos,
                      FST,
                      fst_percentile,
                      FST_pos,
                      lfc,
                      p,
                      nltp,
                      score_type,
                      any_of(edircols))

    valsfrom <- c('lfc','p','nltp')
    if ('binary_direction' %in% colnames(chipatac)) {
        valsfrom <- c(valsfrom, 'enh_idx', 'enh_type', 'direction', 'binary_direction', 'enh_distance')
    }

    chipatac <- chipatac %>%
        pivot_wider(names_from=score_type,
                    names_glue="{score_type}_{.value}",
                    values_from=valsfrom)

    outdf <- chipatac

    # avoids merging to only diff-ABC score enhancers
    # when applying this function to df that has been filtered
    # to reinforcing/opposing diff-enhancers
    if (!('mean_ABC' %in% colnames(chipatac))) {
        meanabc <- dabc %>%
            dplyr::filter(score_type=='ABC') %>%
            mutate(mean_ABC = (Score_AFR + Score_EUR)/2) %>%
            dplyr::select(ensembl_gene_id,
                          entrez_id,
                          hgnc_symbol,
                          TargetGeneTSS,
                          seqnames,
                          start,
                          end,
                          class,
                          name,
                          distance,
                          isSelfPromoter,
                          mean_ABC)

        outdf <- merge(meanabc, outdf)
    }

    return(outdf)

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
                                  return_all = FALSE,
                                  orig_cnames = FALSE,
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

        df <- merge(meanabc, df)

    } else {

        df <- df %>%
            mutate(mean.Score := (.data[[score_col1]] + .data[[score_col2]]) / 2)

    }

    df <- df %>%
        dplyr::rename(abc_p := {{abc_p_col}},
                      abc_lfc := {{abc_lfc_col}},
                      gene := {{gene_col}})

    dfs <- df %>%
        dplyr::filter(abc_p < p_thresh) %>%
        group_by(gene)

    if (allscores) dfs <- dfs %>% group_by(score_type, .add=TRUE)

    if (is.numeric(nEnh)) {
        dfs <- dfs %>% dplyr::filter(dplyr::n()==nEnh)
    } else {
        dfs <- dfs %>% dplyr::filter(dplyr::n() > 1)
    }

    dfs <- dfs %>%
        arrange(desc(mean.Score), .by_group=TRUE) %>%
        mutate(enh_idx = 1:dplyr::n(),
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
              ) %>%
        group_by(binary_direction, .add=TRUE) %>%
        # max distance between diff-enhancers of the same target gene
        mutate(enh_distance = max(abs(distance[enh_type=='primary'] - distance[enh_type=='secondary']))) %>%
        ungroup(binary_direction)

    if (!is.null(dist_thresh)) {

        dfs <- dfs %>%
            dplyr::filter(enh_distance > dist_thresh)

    }

    outdf <- dfs

    if (return_all) outdf <- merge(df, outdf, all.x=TRUE)

    if (orig_cnames) {
        outdf <- outdf %>%
            dplyr::rename(!!abc_p_col := abc_p,
                          !!abc_lfc_col := abc_lfc,
                          !!gene_col := gene) %>%
            ungroup()
    }

    if (abcrank) outdf <- outdf %>% dplyr::rename(mean_ABC=mean.Score)

    return(outdf)

}


parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--class_resolution',
                    type='character',
                    default=NULL)
parser$add_argument('--qtl_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--qtl_name',
                    type='character',
                    default=NULL)
parser$add_argument('--d_all_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--gs_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--frq_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--fst_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--chain_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--ref_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--eqtl_lfsr_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--eqtl_es_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--de_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--rde_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--rdr_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--teqtl_fname',
                    type='character',
                    default=NULL)
parser$add_argument('-n', '--name',
                    type='character',
                    default='diffscore_QTL_overlap',
                    help="Basename for output files")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

qtl_fname <- opt$qtl_fname
qtl_name <- opt$qtl_name
outbase <- file.path(opt$outdir, paste0(opt$name, '.', qtl_name))

if (is.null(qtl_fname)) stop('QTL file required.')

d_all_fname <- opt$d_all_fname
gs_fname <- opt$gs_fname
frq_fname <- opt$frq_fname
fst_fname <- opt$fst_fname


chain_fname <- opt$chain_fname
ref_fname <- opt$ref_fname
eqtl_lfsr_fname <- opt$eqtl_lfsr_fname
eqtl_es_fname <- opt$eqtl_es_fname

de_fname <- opt$de_fname
rde_fname <- opt$rde_fname
rdr_fname <- opt$rdr_fname
teqtl_fname <- opt$teqtl_fname

classres <- opt$class_resolution


########## FIXED ARGS ###############
edircols <- c('mean_ABC','enh_idx','enh_type','direction','binary_direction','enh_distance')
firstcols <- c('hgnc_symbol','name','TargetGeneTSS','distance','mean_ABC','width','rank2_percentile',
               'fst_percentile','ihs_pos','FST_pos','bqtl_pos','bqtl_pos_all','eqtl_pos_lcl','eqtl_pos_pbmc',
               'bqtl_alt_affinity','bqtl_alt_high_freq','bqtl_concordance','type_de_lcl','direction_de_lcl',
               'type_de_pbmc','direction_de_pbmc')

########## READ IN DATA #############

# ancestry relative allele freq data
frq <- read_tsv(frq_fname)
fst <- read_tsv(fst_fname)
# bQTL data (either all tested or sig only (CTCF))
qtl <- read_tsv(qtl_fname)
if (qtl_name=='ctcf') {
    qtl <- qtl %>%
        mutate(bqtl_alt_affinity = if_else(beta < 0, 'low', 'high')) %>%
        dplyr::select(-beta)
} else {
    qtl <- qtl %>%
        mutate(bqtl_alt_affinity = case_when(POSTfreq > prechipfreq ~ 'low',
                                             POSTfreq < prechipfreq ~ 'high',
                                             POSTfreq == prechipfreq ~ 'neutral'))
}
qtl <- qtl %>%
     merge(
     .,
     frq %>%
         dplyr::select(-start) %>%
         dplyr::rename(Chr=seqnames,
                       position=end),
     all.x=TRUE
  ) %>%
  merge(.,
       fst %>%
           dplyr::rename(Chr=seqnames,
                  position=start,
                  bqtl_fst=FST,
                  bqtl_fst_percentile=fst_percentile),
       all.x=TRUE)

# gene set
gs <- readLines(gs_fname)


# Lea et al. eQTL data
es_pre <- 'es'
pv_pre <- 'lfsr'
ref <- read_tsv(ref_fname, col_types=cols_only('c','c'))
eqtl_lfsr <- read_tsv(eqtl_lfsr_fname)
eqtl_es <- read_tsv(eqtl_es_fname)
print('Binding ref to eQTL')
print(dim(ref))
print(dim(eqtl_es))
print(dim(eqtl_lfsr))
eqtl <- rbind(
    cbind(ref, eqtl_es) %>%
        mutate(stat = 'es'),
    cbind(ref, eqtl_lfsr) %>%
        mutate(stat = 'lfsr')
)
print('Pivoting eQTL')
eqtl <- eqtl %>%
    pivot_wider(names_from=stat,
                names_glue="{stat}_{.value}",
                values_from=3:last_col(1)) %>%
    dplyr::rename_with(~paste0(.x, '_LCL_bulk_eqtl'), dplyr::starts_with(c('es_','lfsr_'))) %>%
    separate(SNP, c('seqnames','start'), sep=":", remove=TRUE, convert=TRUE) %>%
    mutate(seqnames = paste0('chr',seqnames)) %>%
    as_granges(width=1)

eqtl_19 <- liftover_grange(eqtl, chain_fname) %>%
    as_tibble %>%
    mutate(eqtl_pos = start) %>%
    dplyr::rename(eGene = gene) %>%
    as_granges

eqtl_19tf <- eqtl_19 %>%
    as_tibble %>%
    pivot_longer(cols=matches(paste0(es_pre,'_|',pv_pre,'_')),
                 names_to=c(".value","condition","celltype","datatype"),
                 names_pattern="(.*)_(.*)_(.*)_(.*)_.*") %>%
    dplyr::rename_with(~paste0(.x, '_eqtl'), dplyr::matches(paste0(es_pre,'|',pv_pre))) %>%
    add_sharing_status(.,
                       es_pre='es_eqtl',
                       pv_pre='lfsr_eqtl',
                       gene_col='eGene',
                       var_col='eqtl_pos',
                       eqtl=TRUE,
                       sharing_def='es',
                       negate_beta=FALSE,
                       pivot=FALSE,
                       pivotw=TRUE) %>%
    rowwise() %>%
    mutate(eqtl_p_min = min(c_across(dplyr::starts_with('lfsr_')))) %>%
    ungroup() %>%
    merge(
        .,
        frq %>%
            mutate(start=start+1) %>%
            dplyr::rename(seqnames_eqtl=seqnames),
        all.x=TRUE
    ) %>%
    merge(.,
         fst %>%
             dplyr::rename(seqnames_eqtl=seqnames,
                    eqtl_fst=FST,
                    eqtl_fst_percentile=fst_percentile),
         all.x=TRUE)

# Lea et al. DE data
es_pre <- 'beta'
pv_pre <- 'pval'
de <- read_tsv(de_fname)
de <- add_hgnc(de)
# remove duplicated gene IDs
# all duplicated gene ids are gene on chrX
ensdups <- de %>%
    group_by(ensembl_gene_id,hgnc_symbol) %>%
    summarize(N = dplyr::n()) %>%
    dplyr::filter(N>1) %>%
    pull(ensembl_gene_id)
de <- de %>%
    as_tibble %>%
    dplyr::filter(!(ensembl_gene_id %in% ensdups)) %>%
    pivot_longer(cols=matches(paste0(es_pre,'_|',pv_pre,'_')),
                 names_to=c(".value","condition","celltype","datatype"),
                 names_pattern="(.*)_(.*)_(.*)_(.*)") %>%
    add_sharing_status(.,
                       es_pre=es_pre,
                       pv_pre=pv_pre,
                       gene_col='ensembl_gene_id',
                       var_col='eqtl_pos',
                       eqtl=FALSE,
                       sharing_def='es',
                       pivot=FALSE)
dew <- de %>%
   pivot_wider(names_from=c(condition,celltype,datatype),
               names_glue="{.value}_{condition}_{celltype}_{datatype}",
               values_from=c(beta, pval))

# Randolph et al. DE data
# Randolph et al. eQTL data
es_pre <- 'beta'
pv_pre <- 'lfsr'
# DE per condition
rde <- read_tsv(rde_fname) %>%
    add_sharing_status(.,
                      es_pre=es_pre,
                       pv_pre=pv_pre,
                       gene_col='hgnc_symbol',
                       var_col='eqtl_pos',
                       eqtl=FALSE,
                       sharing_def='lfsr',
                       negate_beta=TRUE,
                       pivot=TRUE) %>%
    pivot_wider(names_from=c(condition,celltype,datatype),
                names_glue="{.value}_{condition}_{celltype}_{datatype}",
                values_from=c(beta, lfsr))

# differential response across conditions by ancestry
rdr <- read_tsv(rdr_fname) %>%
    pivot_longer(cols=matches(paste0(es_pre,'_|',pv_pre,'_')),
                 names_to=c(".value","celltype"),
                 names_pattern="(.*)_(.*)") %>%
    mutate(condition='fluResponse',
           datatype='pseudobulk') %>%
     add_sharing_status(.,
                       es_pre=es_pre,
                        pv_pre=pv_pre,
                        gene_col='hgnc_symbol',
                        var_col='eqtl_pos',
                        eqtl=FALSE,
                        sharing_def='lfsr',
                        negate_beta=TRUE,
                        pivot=FALSE,
                        pivotw=TRUE)
# top cis-SNP per tested gene across cell types
teqtl_lfsr <- read_excel(teqtl_fname, sheet=1, skip=7) %>%
    separate(gene_SNP,
             into=c('hgnc_symbol', 'seqnames', 'start', 'ref', 'alt'),
             sep='_',
             convert=TRUE) %>%
    mutate(seqnames = paste0('chr', seqnames)) %>%
    rename_with(.,
                ~ paste('lfsr', str_split(.x,'_',simplify=TRUE)[,2], str_split(.x,'_',simplify=TRUE)[,1], sep='_'),
                .cols=6:ncol(.))
teqtl_beta <- read_excel(teqtl_fname, sheet=2) %>%
    separate(gene_SNP,
             into=c('hgnc_symbol', 'seqnames', 'start', 'ref', 'alt'),
             sep='_',
             convert=TRUE) %>%
    mutate(seqnames = paste0('chr', seqnames)) %>%
    rename_with(.,
                ~ paste('beta', str_split(.x,'_',simplify=TRUE)[,2], str_split(.x,'_',simplify=TRUE)[,1], sep='_'),
                .cols=6:ncol(.))
teqtl <- merge(teqtl_beta,teqtl_lfsr) %>%
    rename_with(.,
                ~ paste(.x, 'pseudobulk', sep='_'),
                .cols=6:ncol(.)) %>%
    as_granges(width=1) %>%
    liftover_grange(., chain_fname) %>%
    as_tibble() %>%
    pivot_longer(cols=matches(paste0(es_pre,'_|',pv_pre,'_')),
                 names_to=c(".value","condition","celltype","datatype"),
                 names_pattern="(.*)_(.*)_(.*)_(.*)") %>%
    dplyr::rename_with(~paste0(.x, '_eqtl'), dplyr::matches(paste0(es_pre,'|',pv_pre))) %>%
    dplyr::rename(eGene = hgnc_symbol,
                  seqnames_eqtl = seqnames) %>%
    mutate(eqtl_pos = start) %>%
    add_sharing_status(.,
               es_pre='beta_eqtl',
               pv_pre='lfsr_eqtl',
               gene_col='eGene',
               var_col='eqtl_pos',
               eqtl=TRUE,
               sharing_def='lfsr',
               negate_beta=FALSE,
               pivot=FALSE,
               pivotw=TRUE) %>%
    rowwise() %>%
    mutate(eqtl_p_min = min(c_across(dplyr::starts_with('lfsr_')))) %>%
    ungroup() %>%
    merge(
        .,
        frq %>%
            mutate(start=start+1) %>%
            dplyr::rename(seqnames_eqtl=seqnames,
                          frq_ref=ref,
                          frq_alt=alt),
        all.x=TRUE
    ) %>%
    # some eQTL 'ref' alleles don't match the frequency file confirmed 1000 genomes ref allele
    # ignore these for now, but PBMC format could be minor allele second
    # regression likely still done on number of alternate alleles
    dplyr::select(-c(ref,alt)) %>%
    dplyr::rename(ref=frq_ref,
                  alt=frq_alt) %>%
    merge(.,
         fst %>%
             dplyr::rename(seqnames_eqtl=seqnames,
                    eqtl_fst=FST,
                    eqtl_fst_percentile=fst_percentile),
         all.x=TRUE)

# diff-score CREs
abc <- read_tsv(d_all_fname)

# get reinforcing enhancers to add later
reinfenh <- addEnhancerDirectionMatching(abc,
                                          dist_thresh=NULL,
                                          p_thresh=0.05,
                                          abc_p_col='p',
                                          abc_lfc_col='log2FoldChange',
                                          gene_col='TargetGene',
                                          nEnh='2+',
                                          score_type='Score',
                                          cdt1 = 'AFR',
                                          cdt2 = 'EUR',
                                          allscores = TRUE,
                                          abcrank = TRUE,
                                          return_all = FALSE,
                                          orig_cnames = TRUE,
                                          debug=FALSE) %>%
    widen_diffABC(.) %>%
    dplyr::select(ensembl_gene_id,
                  entrez_id,
                  hgnc_symbol,
                  TargetGeneTSS,
                  seqnames,
                  start,
                  end,
                  class,
                  name,
                  distance,
                  isSelfPromoter,
                  contains(edircols))

# widen and add bQTL
wabc <- join_overlap_left(
    widen_diffABC(abc) %>% as_granges(),
    qtl %>%
        dplyr::rename(seqnames=Chr,
                      start=position,
                      bqtl_p=pvalue) %>%
        mutate(end=start+1,
               bqtl_pos=start) %>%
        as_granges()
) %>%
    as_tibble() %>%
    mutate(coded_alt_affinity = case_when(bqtl_alt_affinity=='high' ~ 1,
                                          bqtl_alt_affinity=='low' ~ -1),
           coded_alt_high_freq = case_when(alt_high_freq=='AFR' ~ 1,
                                           alt_high_freq=='EUR' ~ -1)) %>%
    group_by(hgnc_symbol, name) %>%
    arrange(bqtl_p, .by_group=TRUE) %>%
    mutate(bqtl_N=sum(bqtl_p < 0.05),
           bqtl_pos_all=paste0(bqtl_pos[bqtl_p < 0.05], collapse=';'),
           bqtl_fst_percentile_all=paste0(bqtl_fst_percentile[bqtl_p < 0.05], collapse=';'),
           bqtl_alt_affinity_all=paste0(bqtl_alt_affinity[bqtl_p < 0.05], collapse=';'),
           bqtl_alt_high_freq_all=paste0(alt_high_freq[bqtl_p < 0.05], collapse=';'),
           # concordance = if high frequency allele is the high or low affinity allele consistently throughout a CRE
           bqtl_concordance=case_when(sum(bqtl_p < 0.05)>0 & (paste0(coded_alt_affinity[bqtl_p < 0.05], collapse=';')==paste0(coded_alt_high_freq[bqtl_p < 0.05], collapse=';') |
                                           paste0(-coded_alt_affinity[bqtl_p < 0.05], collapse=';')==paste0(coded_alt_high_freq[bqtl_p < 0.05], collapse=';')) ~ TRUE,
                                       sum(bqtl_p < 0.05)>0 ~ FALSE)) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::select(-c(ref,alt,AFR_ref_frq,EUR_ref_frq),
                  -starts_with('coded_')) %>%
    dplyr::rename(bqtl_AFR_alt_frq=AFR_alt_frq,
                  bqtl_EUR_alt_frq=EUR_alt_frq,
                  bqtl_alt_high_freq=alt_high_freq)

# add DE/eQTL
wabc_olap <- wabc %>%
    # add gene set info
    mutate(gene_set = if_else(hgnc_symbol %in% gs, 'TNFA_NFKB', 'other')) %>%
    # add LCL condition DE
    merge(
        .,
        dew %>%
            dplyr::select(1:4) %>%
            dplyr::rename(type_de_lcl = type_de,
                          direction_de_lcl = direction_de),
        all.x=TRUE
    ) %>%
    as_granges() %>%
    # include eQTL for non-ABC-target eGenes in case the eQTL is for an eGene in a gene set of interest
    join_overlap_left(
        .,
        eqtl_19tf %>%
            dplyr::select(-c(ref,alt,AFR_ref_frq,EUR_ref_frq)) %>%
            dplyr::rename(seqnames=seqnames_eqtl,
                          type_eqtl_lcl=type_eqtl,
                          direction_eqtl_lcl=direction_eqtl) %>%
            dplyr::select(1:9,all_of(c('eqtl_p_min','AFR_alt_frq','EUR_alt_frq','alt_high_freq','eqtl_fst','eqtl_fst_percentile'))) %>%
            as_granges
    ) %>%
    as_tibble %>%
    mutate(coded_eqtl_direction = case_when(direction_eqtl_lcl=='up' ~ 1,
                                          direction_eqtl_lcl=='down' ~ -1),
           coded_alt_high_freq = case_when(alt_high_freq=='AFR' ~ 1,
                                           alt_high_freq=='EUR' ~ -1)) %>%
    group_by(hgnc_symbol,name,eGene,type_eqtl_lcl) %>%
    arrange(eqtl_p_min, .by_group=TRUE) %>%
    mutate(type_eqtl_lcl_N=sum(eqtl_p_min<0.1, na.rm=TRUE),
           type_eqtl_lcl_pos=paste0(eqtl_pos[eqtl_p_min<0.1], collapse=';'),
           eqtl_lcl_fst_percentile_all=paste0(eqtl_fst_percentile[eqtl_p_min<0.1], collapse=';'),
           eqtl_lcl_alt_expr_all=paste0(direction_eqtl_lcl[eqtl_p_min<0.1], collapse=';'),
           eqtl_lcl_alt_high_freq_all=paste0(alt_high_freq[eqtl_p_min<0.1], collapse=';'),
           eqtl_concordance=case_when(sum(eqtl_p_min < 0.1)>0 & (paste0(coded_eqtl_direction[eqtl_p_min < 0.1], collapse=';')==paste0(coded_alt_high_freq[eqtl_p_min < 0.1], collapse=';') |
                                           paste0(-coded_eqtl_direction[eqtl_p_min < 0.1], collapse=';')==paste0(coded_alt_high_freq[eqtl_p_min < 0.1], collapse=';')) ~ TRUE,
                                       sum(eqtl_p_min < 0.1)>0 ~ FALSE)) %>%
    dplyr::slice(1) %>%
    dplyr::select(-starts_with('coded_')) %>%
    ungroup() %>%
    dplyr::rename(eqtl_lcl_alt_expr = direction_eqtl_lcl,
                  eqtl_lcl_AFR_alt_frq=AFR_alt_frq,
                  eqtl_lcl_EUR_alt_frq=EUR_alt_frq,
                  eqtl_lcl_alt_high_freq = alt_high_freq,
                  eqtl_lcl_p=eqtl_p_min,
                  eGene_lcl=eGene,
                  eqtl_pos_lcl=eqtl_pos,
                  eqtl_lcl_fst=eqtl_fst,
                  eqtl_lcl_fst_percentile=eqtl_fst_percentile) %>%
    # add lead cis-SNP per gene from Randolph et al. PBMCs (ideally rerun their eQTL mapping to include other SNPs)
    as_granges() %>%
    join_overlap_left(
        .,
        teqtl %>%
            dplyr::select(-c(ref,alt,AFR_ref_frq,EUR_ref_frq)) %>%
            dplyr::rename(seqnames=seqnames_eqtl,
                          type_eqtl_pbmc=type_eqtl,
                          eqtl_pbmc_alt_expr=direction_eqtl,
                          eqtl_pos_pbmc=eqtl_pos) %>%
            dplyr::select(1:9,all_of(c('eqtl_p_min','AFR_alt_frq','EUR_alt_frq','alt_high_freq','eqtl_fst','eqtl_fst_percentile'))) %>%
            as_granges
    ) %>%
    as_tibble() %>%
    dplyr::rename(eqtl_pbmc_p=eqtl_p_min,
                  eqtl_pbmc_AFR_alt_frq=AFR_alt_frq,
                  eqtl_pbmc_EUR_alt_frq=EUR_alt_frq,
                  eqtl_pbmc_alt_high_freq=alt_high_freq,
                  eGene_pbmc=eGene,
                  eqtl_pbmc_fst=eqtl_fst,
                  eqtl_pbmc_fst_percentile=eqtl_fst_percentile) %>%
    # add PBMC popDE
    merge(
        .,
        rde %>%
            dplyr::select(1:3) %>%
            dplyr::rename(type_de_pbmc = type_de,
                          direction_de_pbmc = direction_de),
        all.x=TRUE
    ) %>%
    # add PBMC popDR
    merge(
        .,
        rdr %>%
            dplyr::select(1:3) %>%
            dplyr::rename(type_dr_pbmc = type_de,
                          direction_dr_pbmc = direction_de),
        all.x=TRUE
    ) %>%
    # add reinforcing enh
    merge(., reinfenh, all.x=TRUE) %>%
    mutate(atac_monodirection = if_else(atac_lfc>0, 'up', 'down'),
           chip_monodirection = if_else(chip_lfc>0, 'up', 'down'),
           hic_monodirection = if_else(hic_lfc>0, 'up', 'down'),
           de_direction_matches_atac = case_when(
               (!is.na(direction_de_lcl) & !is.na(direction_de_pbmc)) & (atac_monodirection==direction_de_lcl & atac_monodirection==direction_de_pbmc) ~ 'both',
               !is.na(direction_de_lcl) & atac_monodirection==direction_de_lcl ~ 'lcl',
               !is.na(direction_de_pbmc) & atac_monodirection==direction_de_pbmc ~ 'pbmc'),
           de_direction_matches_chip = case_when(
               (!is.na(direction_de_lcl) & !is.na(direction_de_pbmc)) & (chip_monodirection==direction_de_lcl & chip_monodirection==direction_de_pbmc) ~ 'both',
               !is.na(direction_de_lcl) & chip_monodirection==direction_de_lcl ~ 'lcl',
               !is.na(direction_de_pbmc) & chip_monodirection==direction_de_pbmc ~ 'pbmc')
           ) %>%
    dplyr::select(all_of(firstcols), everything()) %>%
    mutate(bqtl_direction = case_when(
                ((bqtl_alt_affinity=='high' & bqtl_alt_high_freq=='AFR') |
                    (bqtl_alt_affinity=='low' & bqtl_alt_high_freq=='EUR')) &
                bqtl_concordance ~ 'down',
                ((bqtl_alt_affinity=='high' & bqtl_alt_high_freq=='EUR') |
                    (bqtl_alt_affinity=='low' & bqtl_alt_high_freq=='AFR')) &
                bqtl_concordance ~ 'up'
            )
          ) %>%
    mutate(eqtl_lcl_direction = case_when(
                ((eqtl_lcl_alt_expr=='up' & eqtl_lcl_alt_high_freq=='AFR') |
                    (eqtl_lcl_alt_expr=='down' & eqtl_lcl_alt_high_freq=='EUR')) &
                eqtl_concordance ~ 'down',
                ((eqtl_lcl_alt_expr=='up' & eqtl_lcl_alt_high_freq=='EUR') |
                    (eqtl_lcl_alt_expr=='down' & eqtl_lcl_alt_high_freq=='AFR')) &
                eqtl_concordance ~ 'up'
            ),
           eqtl_pbmc_direction = case_when(
                (eqtl_pbmc_alt_expr=='up' & eqtl_pbmc_alt_high_freq=='AFR') |
                    (eqtl_pbmc_alt_expr=='down' & eqtl_pbmc_alt_high_freq=='EUR') ~ 'down',
                (eqtl_pbmc_alt_expr=='up' & eqtl_pbmc_alt_high_freq=='EUR') |
                    (eqtl_pbmc_alt_expr=='down' & eqtl_pbmc_alt_high_freq=='AFR') ~ 'up'
            )
          )


######### bQTL ENRICHMENT TESTS FOR AGGREGATING ############

sts <- c('chip','atac','hic')
qts <- c('eqtl_lcl', 'eqtl_pbmc', 'bqtl')
cfs <- c('all','promoter','non-promoter')
fts <- c(0.0, 0.95)
dcts <- c('lcl','pbmc')
spts <- c(TRUE, FALSE)

ares <- tibble()
wres <- tibble() # for storing full Wilcoxon dfs for FST plotting

# score types
for (st in sts) {
  # QTL types
  for (qt in qts) {

    # skip overall enrichment for eQTL and ctcf for now (don't have all tested SNPs for all)
    if (qt=='bqtl' & qtl_name != 'ctcf') {

      # CRE class filters (cre_type)
      for (cf in cfs) {

        res <- diffQtlFishers(wabc_olap,
                              qtl_p_thresh=0.05,
                              qtl_p_nonthresh=0.5,
                              qtl_p_col=paste0(qt, '_p'),
                              abc_p_thresh=0.05,
                              abc_p_nonthresh=0.5,
                              abc_p_col=paste0(st, '_p'),
                              class_filt=cf,
                              class_res=classres,
                              mat_out=FALSE,
                              alt='greater') %>%
              as_tibble() %>%
              mutate(test_type='diffQTL',
                     qtl_type=qt,
                     score_type=st,
                     cre_type=cf,
                     fst_percentile=0.0,
                     de_celltype=NA,
                     self_promoter=NA,
                     n_success=NA,
                     n_trials=NA,
                     null_prob=NA)

        ares <- rbind(ares, res)

      }

    }

    # FST thresholds
    for (ft in fts) {

      # CRE class filters (cre_type)
      for (cf in cfs) {

        res <- diffQTLdirectionFishers(wabc_olap,
                                      qtl_p_thresh=0.05,
                                      qtl_p_col=paste0(qt, '_p'),
                                      qtl_dir_col=paste0(qt, '_direction'),
                                      abc_p_thresh=0.05,
                                      abc_p_col=paste0(st, '_p'),
                                      abc_dir_col=paste0(st, '_monodirection'),
                                      class_filt=cf,
                                      class_res=classres,
                                      fst_perc_thresh=ft,
                                      fst_perc_col=paste0(qt, '_fst_percentile'),
                                      mat_out=FALSE,
                                      alt='greater') %>%
                      as_tibble() %>%
                      mutate(test_type='diffQTLdirection',
                             qtl_type=qt,
                             score_type=st,
                             cre_type=cf,
                             fst_percentile=ft,
                             de_celltype=NA,
                             self_promoter=NA,
                             n_success=NA,
                             n_trials=NA,
                             null_prob=NA)

        ares <- rbind(ares, res)

        # sign test for diff-score identification of QTL directional selection
        # with binomial test where 'success' -> diff-score - bQTL both in
        # 'down' direction
        res <- diffQTLsignTest(wabc_olap,
                              qtl_p_thresh=0.05,
                              qtl_p_col=paste0(qt, '_p'),
                              qtl_dir_col=paste0(qt, '_direction'),
                              abc_p_thresh=0.05,
                              abc_p_col=paste0(st, '_p'),
                              abc_dir_col=paste0(st, '_monodirection'),
                              fst_perc_thresh=ft, # [0.0-1.0]
                              fst_perc_col=paste0(qt, '_fst_percentile'),
                              class_filt=cf, # 'promoter', 'non-promoter', 'none'
                              class_res='chip',
                              alt='two.sided') %>%
                      as_tibble() %>%
                      mutate(test_type='binomSignTest',
                            qtl_type=qt,
                            score_type=st,
                            cre_type=cf,
                            fst_percentile=ft,
                            de_celltype=NA,
                            self_promoter=NA,
                            up_up = NA,
                             up_down = NA,
                             down_up = NA,
                             down_down = NA)

        ares <- rbind(ares, res)
      }

      # DE cell types
      for (dct in dcts) {

        if ((grepl('lcl', qt) & dct=='pbmc') |
            (grepl('pbmc', qt) & dct=='lcl')) {
          next
        }

        # no self promoter filtering
        res <- DEdiffQTLdirectionFishers(wabc_olap,
                                        qtl_p_thresh=0.05,
                                        qtl_p_col=paste0(qt, '_p'),
                                        qtl_dir_col=paste0(qt, '_direction'),
                                        abc_p_thresh=0.05,
                                        abc_p_col=paste0(st, '_p'),
                                        abc_dir_col=paste0(st, '_monodirection'),
                                        de_dir_col=paste0('direction_de_', dct),
                                        fst_perc_thresh=ft,
                                        fst_perc_col=paste0(qt, '_fst_percentile'),
                                        class_filt='none',
                                        class_res=classres,
                                        self_promoter=NULL,
                                        mat_out=FALSE,
                                        alt='greater') %>%
                        as_tibble() %>%
                        mutate(test_type='DEdiffQTLdirection',
                               qtl_type=qt,
                               score_type=st,
                               cre_type='all',
                               fst_percentile=ft,
                               de_celltype=dct,
                               self_promoter=NA,
                               n_success=NA,
                               n_trials=NA,
                               null_prob=NA)

          ares <- rbind(ares, res)

        # self promoter only vs non-self promoter only
        for (spt in spts) {

          res <- DEdiffQTLdirectionFishers(wabc_olap,
                                          qtl_p_thresh=0.05,
                                          qtl_p_col=paste0(qt, '_p'),
                                          qtl_dir_col=paste0(qt, '_direction'),
                                          abc_p_thresh=0.05,
                                          abc_p_col=paste0(st, '_p'),
                                          abc_dir_col=paste0(st, '_monodirection'),
                                          de_dir_col=paste0('direction_de_', dct),
                                          fst_perc_thresh=ft,
                                          fst_perc_col=paste0(qt, '_fst_percentile'),
                                          class_filt='none',
                                          class_res=classres,
                                          self_promoter=spt,
                                          mat_out=FALSE,
                                          alt='greater') %>%
                          as_tibble() %>%
                          mutate(test_type='DEdiffQTLdirection',
                                 qtl_type=qt,
                                 score_type=st,
                                 cre_type=NA,
                                 fst_percentile=ft,
                                 de_celltype=dct,
                                 self_promoter=spt,
                                 n_success=NA,
                                 n_trials=NA,
                                 null_prob=NA)

            ares <- rbind(ares, res)

        }

      }

    }

    # Wilcoxon tests on FST (don't loop through FST thresholds)
    # but loop through class filters except 'all'
    for (cf in cfs) {

      if (cf=='all') next

      res <- diffQTLFSTWilcoxon(wabc_olap,
                                  qtl_p_thresh=0.05,
                                  qtl_p_col=paste0(qt, '_p'),
                                  qtl_dir_col=paste0(qt, '_direction'),
                                  abc_p_thresh=0.05,
                                  abc_p_nonthresh=0.05,
                                  abc_p_col=paste0(st, '_p'),
                                  abc_dir_col=paste0(st, '_monodirection'),
                                  fst_col=paste0(qt, '_fst'),
                                  class_filt=cf, # 'promoter', 'non-promoter', 'none'
                                  class_res='chip') %>%
                    mutate(cre_type=cf,
                           qtl_type=qt,
                           score_type=st)

        wres <- rbind(wres, res)

    }

  }
}

ares <- ares %>%
  mutate(chip_type=case_when(qtl_type=='bqtl' ~ qtl_name))
wres <- wres %>%
  mutate(chip_type=case_when(qtl_type=='bqtl' ~ qtl_name))

write_tsv(ares, paste0(outbase, ".fisherEnrichments.txt"))
write_tsv(wres, paste0(outbase, ".fstWilcoxData.txt.gz"))
write_tsv(wabc_olap, paste0(outbase, ".txt.gz"))
