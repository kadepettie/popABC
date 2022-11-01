#!/usr/bin/R

if(!(require(biomaRt))) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
if(!(require(plyranges))) {
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("plyranges")
}

if(!(require(argparse) )) install.packages("argparse")
if(!(require(tidyverse) )) install.packages("tidyverse")
if(!(require(readxl) )) install.packages("readxl")
if(!(require(rstatix)   )) install.packages("rstatix")
if(!(require(ggplot2)   )) install.packages("ggplot2")
if(!(require(ggrepel)   )) install.packages("ggrepel")

library(argparse)
library(tidyverse)
library(readxl)
library(rstatix)
library(plyranges)
library(biomaRt)
library(ggplot2)
library(ggrepel)


###### FUNCTIONS #######
direction_matches <- function(x) {
    if (all(x >= 0) | all(x < 0)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

abc_de_fishers <- function(df,
                      de_p_thresh=0.05,
                      de_p_nonthresh=0.5,
                      de_p_col='padj_DE',
                      abc_p_thresh=0.05,
                      abc_p_nonthresh=0.5,
                      abc_p_col='p_ABC',
                      self_promoter=NULL,
                      doublecount_enh=TRUE,
                      mat_out=FALSE,
                      alt='two.sided') {

    countsdf <- df %>%
      dplyr::rename(de_p := {{de_p_col}},
                    abc_p := {{abc_p_col}}) %>%
      mutate(abc_type = case_when(abc_p < abc_p_thresh ~ 'diff-ABC',
                                  abc_p >= abc_p_nonthresh ~ 'non-diff-ABC',
                                  TRUE ~ 'ambig'),
             de_type = case_when(de_p < de_p_thresh ~ 'DE',
                                  de_p >= de_p_nonthresh ~ 'non-DE',
                                  TRUE ~ 'ambig')) %>%
      dplyr::filter(abc_type != 'ambig',
                    de_type != 'ambig') %>%
      group_by(name)

    if (!doublecount_enh) {
      # count each enh only once and a hit as any enh w/ at least 1 DE target gene

      # for non-self-promoter test, a hit can be an enh that is a self-promoter
      # for a non-DE gene, but not a DE gene
      # a diff-CRE linked only to non-DE genes is counted as non-hit only if it is
      # not a self-promoter for any of those genes
      # a non-diff-CRE linked only to non-DE genes is counted as a 'hit' only
      # if it is not a self-promoter for any of those genes
      countsdf <- countsdf %>%
        arrange(de_type, desc(isSelfPromoter), .by_group=TRUE) %>%
        dplyr::slice(1)
    } else {
      # if doublecount_enh==TRUE, define includes_selfpromoter within DE vs non-DE groups
      countsdf <- countsdf %>%
        group_by(de_type, .add=TRUE)
    }

    # if doublecount_enh==TRUE (aka gene-level analysis), a DE gene is a hit if linked
    # to at least one diff-CRE, but that same diff-CRE may count toward a non-hit if
    # also linked to a non-DE gene, and vice versa
    countsdf <- countsdf %>%
      mutate(includes_selfpromoter = if_else(any(isSelfPromoter), TRUE, FALSE)) %>%
      ungroup() %>%
      group_by(abc_type)

    if (!is.null(self_promoter)) {
      # if doublecount_enh==TRUE for non-self-promoter test, a DE gene is a hit
      # if linked to a diff-CRE that is not a self-promoter for it or any other DE gene,
      # but this diff-CRE may be a self-promoter for a non-DE gene (for which
      # category it would not count toward non-hits)
      if (self_promoter) {
          countsdf <- countsdf %>% dplyr::filter(includes_selfpromoter)
      } else {
          countsdf <- countsdf %>% dplyr::filter(!includes_selfpromoter)
      }

    }

    countsdf <- countsdf %>%
      summarize(
        gene_count = sum(!is.na(de_p)),
        de_gene_count =sum(de_p < de_p_thresh, na.rm=TRUE)
      )

    fg_numerator <- countsdf %>%
        dplyr::filter(abc_type=='diff-ABC') %>%
        pull(de_gene_count)

    fg_denominator <- countsdf %>%
        pull(de_gene_count) %>%
        sum()

    bg_numerator <- countsdf %>%
        dplyr::filter(abc_type=='diff-ABC') %>%
        pull(gene_count)

    bg_denominator <- countsdf %>%
        dplyr::select(gene_count) %>%
        sum()

    fe <- (fg_numerator/fg_denominator) / (bg_numerator/bg_denominator)

    mat11 <- fg_numerator
    mat21 <- fg_denominator-fg_numerator
    mat12 <- bg_numerator-fg_numerator
    mat22 <- bg_denominator-mat11-mat21-mat12
    tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                      dimnames =
               list(c('diff-ABC', 'non-diff-ABC'),
                    c('DE', 'non-DE')))

    ft <- fisher.test(tmat, alternative=alt)

    if (mat_out) {
        return(
            list(
                tmat,
                ft
            )
        )
    }

    eres <- list(de_significance = de_p_thresh,
                 enrichment = fe,
                 odds_ratio = ft$estimate[[1]],
                 pvalue = ft$p.value,
                 conf_lower = ft$conf.int[1],
                 conf_upper = ft$conf.int[2],
                 fg_success=fg_numerator,
                 fg_total=fg_denominator,
                 bg_success=bg_numerator,
                 bg_total=bg_denominator)

    return(eres)

}

abc_de_directional_fishers <- function(df,
                                      de_p_thresh=0.05,
                                      de_p_col='lfsr',
                                      abc_p_thresh=0.05,
                                      abc_p_col='p_ABC',
                                      negate_beta=FALSE, # for matching direction of log2FoldChange_ABC
                                      self_promoter=NULL,
                                      doublecount_enh=TRUE,
                                      mat_out=FALSE,
                                      alt='two.sided') {

    dfs <- df %>%
        dplyr::rename(de_p := {{de_p_col}},
                      abc_p := {{abc_p_col}}) %>%
        dplyr::filter(de_p < de_p_thresh,
                      abc_p < abc_p_thresh)
    if (negate_beta) {
        dfs <- dfs %>% mutate(beta = -beta)
    }

    # rename for compatibility with old function
    dfs <- dfs %>%
        dplyr::rename(de_log2FC = beta,
                      da_log2FC = log2FoldChange_ABC)

    if (!doublecount_enh) {
      # count each diff-CRE only once and a hit as any diff-CRE w/ all DE target
      # genes in the same direction

      # for non-self-promoter test, for each category a CRE must not be a self-promoter
      # for any of its target genes
      # for self-promoter test, for each category a CRE must be a self-promoter
      # for at least one of its target genes

      dfs <- dfs %>%
          group_by(name) %>%
          mutate(de_abc_dir = case_when(
                    all(de_log2FC > 0) & da_log2FC > 0 ~ 'up_up',
                    all(de_log2FC > 0) & da_log2FC < 0 ~ 'up_down',
                    all(de_log2FC < 0) & da_log2FC > 0 ~ 'down_up',
                    all(de_log2FC < 0) & da_log2FC < 0 ~ 'down_down'
                )
          )

    } else {
      # count each gene once, may double count diff-CREs with multiple target genes
      dfs <- dfs %>%
          group_by(ensembl_gene_id) %>%
          mutate(de_abc_dir = case_when(
                    de_log2FC > 0 & sum(da_log2FC > 0, na.rm=TRUE) >= 0.5 ~ 'up_up',
                    de_log2FC > 0 & sum(da_log2FC < 0, na.rm=TRUE) >= 0.5 ~ 'up_down',
                    de_log2FC < 0 & sum(da_log2FC > 0, na.rm=TRUE) >= 0.5 ~ 'down_up',
                    de_log2FC < 0 & sum(da_log2FC < 0, na.rm=TRUE) >= 0.5 ~ 'down_down'
                )
          )

    }

    # still grouped by gene or CRE (from above)
    dfs <- dfs %>%
        dplyr::filter(!is.na(de_abc_dir)) %>%
        mutate(includes_selfpromoter = if_else(any(isSelfPromoter), TRUE, FALSE)) %>%
        dplyr::slice(1) %>%
        ungroup()

    if (!is.null(self_promoter)) {
        if (self_promoter) {
            dfs <- dfs %>% dplyr::filter(includes_selfpromoter)
        } else {
            dfs <- dfs %>% dplyr::filter(!includes_selfpromoter)
        }
    }

    # # upreg genes linked to majority increased ABC score pairs
    de_up_da_up <- dfs %>%
        dplyr::filter(de_abc_dir=="up_up") %>%
        nrow()

    de_up_da_down <- dfs %>%
        dplyr::filter(de_abc_dir=="up_down") %>%
        nrow()

    de_down_da_up <- dfs %>%
        dplyr::filter(de_abc_dir=="down_up") %>%
        nrow()

    de_down_da_down <- dfs %>%
        dplyr::filter(de_abc_dir=="down_down") %>%
        nrow()

    mat11 <- de_up_da_up
    mat21 <- de_up_da_down
    mat12 <- de_down_da_up
    mat22 <- de_down_da_down
    tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                      dimnames =
               list(c('DA up', 'DA down'),
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

    eres <- list(de_significance = de_p_thresh,
                 odds_ratio = ft$estimate[[1]],
                 pvalue = ft$p.value,
                 conf_lower = ft$conf.int[1],
                 conf_upper = ft$conf.int[2],
                 DE_up_DA_up = mat11,
                 DE_up_DA_down = mat21,
                 DE_down_DA_up = mat12,
                 DE_down_DA_down = mat22)

    return(eres)

}

fisher_enrichments <- function(df,
                               de_p_col='p_DE',
                               abc_p_col='p_ABC',
                               conditions=df %>% pull(condition) %>% unique(),
                               celltypes=df %>% pull(celltype) %>% unique(),
                               datatypes=df %>% pull(datatype) %>% unique()) {

    agg_df <- tibble()
    for (cond in conditions) {
        for (ct in celltypes) {
            for (dt in datatypes) {
                for (dce in c(TRUE, FALSE)) {

                    dfs <- df %>%
                        dplyr::filter(condition==cond,
                                      celltype==ct,
                                      datatype==dt)

                    if (nrow(dfs)==0) next

                    curr_df <- as_tibble(
                        abc_de_fishers(dfs,
                              de_p_thresh=0.05,
                              de_p_nonthresh=0.5,
                              de_p_col=de_p_col,
                              abc_p_thresh=0.05,
                              abc_p_nonthresh=0.5,
                              abc_p_col=abc_p_col,
                              self_promoter=NULL,
                              doublecount_enh=dce,
                              mat_out=FALSE,
                              alt='greater')
                    )
                    curr_df['condition'] <- cond
                    curr_df['celltype'] <- ct
                    curr_df['datatype'] <- dt
                    curr_df['self_promoter'] <- NA
                    curr_df['doublecount_enh'] <- dce

                    agg_df <- rbind(agg_df, curr_df)

                    for (spt in c(TRUE, FALSE)) {

                        curr_df <- as_tibble(
                            abc_de_fishers(dfs,
                                  de_p_thresh=0.05,
                                  de_p_nonthresh=0.5,
                                  de_p_col=de_p_col,
                                  abc_p_thresh=0.05,
                                  abc_p_nonthresh=0.5,
                                  abc_p_col=abc_p_col,
                                  self_promoter=spt,
                                  doublecount_enh=dce,
                                  mat_out=FALSE,
                                  alt='greater')
                        )
                        curr_df['condition'] <- cond
                        curr_df['celltype'] <- ct
                        curr_df['datatype'] <- dt
                        curr_df['self_promoter'] <- spt
                        curr_df['doublecount_enh'] <- dce

                        agg_df <- rbind(agg_df, curr_df)

                    }

                }

            }

        }

    }

    return(agg_df)
}

fisher_enrichments_direction <- function(df,
                                de_p_col='p_DE',
                                abc_p_col='p_ABC',
                                negate_beta=FALSE,
                               conditions=df %>% pull(condition) %>% unique(),
                               celltypes=df %>% pull(celltype) %>% unique(),
                               datatypes=df %>% pull(datatype) %>% unique()) {

    agg_df <- tibble()
    for (cond in conditions) {
        for (ct in celltypes) {
            for (dt in datatypes) {
                for (dce in c(TRUE, FALSE)) {

                    dfs <- df %>%
                        dplyr::filter(condition==cond,
                                      celltype==ct,
                                      datatype==dt)

                    if (nrow(dfs)==0) next

                    curr_df <- as_tibble(
                        abc_de_directional_fishers(dfs,
                                                  de_p_thresh=0.05,
                                                  de_p_col=de_p_col,
                                                  abc_p_thresh=0.05,
                                                  abc_p_col=abc_p_col,
                                                  negate_beta=negate_beta, # for matching direction of log2FoldChange_ABC
                                                  doublecount_enh=dce,
                                                  mat_out=FALSE,
                                                  alt='greater')
                    )
                    curr_df['condition'] <- cond
                    curr_df['celltype'] <- ct
                    curr_df['datatype'] <- dt
                    curr_df['self_promoter'] <- NA
                    curr_df['doublecount_enh'] <- dce

                    agg_df <- rbind(agg_df, curr_df)

                    for (spt in c(TRUE, FALSE)) {

                        curr_df <- as_tibble(
                            abc_de_directional_fishers(dfs,
                                                      de_p_thresh=0.05,
                                                      de_p_col=de_p_col,
                                                      abc_p_thresh=0.05,
                                                      abc_p_col=abc_p_col,
                                                      negate_beta=negate_beta, # for matching direction of log2FoldChange_ABC
                                                      self_promoter=spt,
                                                      doublecount_enh=dce,
                                                      mat_out=FALSE,
                                                      alt='greater')
                        )
                        curr_df['condition'] <- cond
                        curr_df['celltype'] <- ct
                        curr_df['datatype'] <- dt
                        curr_df['self_promoter'] <- spt
                        curr_df['doublecount_enh'] <- dce

                        agg_df <- rbind(agg_df, curr_df)

                    }

                }

            }

        }

    }

    return(agg_df)
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
        arrange(ensembl_gene_id, .by_group=TRUE) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        arrange(p)

    return(df)

}

gprofilerTryCatch <- function(fg, bg, sig_thresh=0.05, debug=FALSE) {

    print("Running g:Profiler query...")

    ordsource <- c('GO:BP','KEGG','HP','HPA','REAC','CORUM','TF','GO:CC','GO:MF','MIRNA','WP')

    gostres <- tryCatch(
      {
        gostres <- gost(fg,
                        organism = "hsapiens",
                        ordered_query = FALSE,
                        multi_query = FALSE,
                        significant = TRUE,
                        exclude_iea = TRUE,
                        measure_underrepresentation = FALSE,
                        evcodes = TRUE,
                        user_threshold = sig_thresh,
                        correction_method ="g_SCS",
                        domain_scope = "known",
                        custom_bg = bg,
                        numeric_ns = "",
                        sources = NULL)
      },
      error=function(e) {
        message(e)
        return(NULL)
      }
    )

    if (is.null(gostres)) {
      print("bad request")
      return(NULL)
    }


    if (debug) {
        return(gostres)
    }

    if (is.null(gostres$result)) {
      print("No results to show")
      return(NULL)
    } else {
      gostdf <- gostres$result %>%
          as_tibble() %>%
          dplyr::select(-c(parents)) %>% # remove list column before writing
          arrange(factor(source, levels=ordsource), p_value)

      gostdf_summ <- gostdf %>%
          dplyr::select(source,
                        term_name,
                        p_value,
                        query_size,
                      term_size,
                    intersection_size)


      print(gostdf_summ %>% dplyr::filter(source=='GO:BP'))

      return(gostdf)
    }

}

enhancerDirectionTest <- function(df,
                                  dist_thresh=NULL,
                                  p_thresh=0.05,
                                  abc_p_col='p',
                                  abc_lfc_col='log2FoldChange',
                                  gene_col='TargetGene',
                                  nEnh=2,
                                  score_type='ABC.Score',
                                  cdt1 = 'AFR',
                                  cdt2 = 'EUR',
                                  direction_out = FALSE,
                                  mat_out=FALSE,
                                  alt='two.sided',
                                  run_gsea=FALSE,
                                  gprof_p=0.05,
                                  debug=FALSE) {

    # gprof_p = threshold for returning gprofiler results

    score_col1 <- paste0(score_type,'_',cdt1)
    score_col2 <- paste0(score_type,'_',cdt2)

    dfs <- df %>%
        dplyr::rename(abc_p := {{abc_p_col}},
                      abc_lfc := {{abc_lfc_col}},
                      gene := {{gene_col}}) %>%
        dplyr::filter(abc_p < p_thresh) %>%
        mutate(mean.Score := (.data[[score_col1]] + .data[[score_col2]]) / 2) %>%
        group_by(gene)

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
                binary_direction = if_else(direction %in% c('up_up','down_down'),
                                           'reinforcing',
                                           'opposing')
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

    if (direction_out) return(dfs %>% ungroup())

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

    eres <- list(abc_significance = p_thresh,
                 de_significance = NA,
                 odds_ratio = ft$estimate[[1]],
                 pvalue = ft$p.value,
                 conf_lower = ft$conf.int[1],
                 conf_upper = ft$conf.int[2],
                 up_up = mat11,
                 up_down = mat21,
                 down_up = mat12,
                 down_down = mat22)

    if (run_gsea) {

        if(!(require(gprofiler2))) install.packages("gprofiler2")
        library(gprofiler2)

        if ('ensembl_gene_id' %in% colnames(dfs)) {
            dfs <- dfs %>% dplyr::rename(hgnc_symbol = gene)
        } else {
            dfs <- dfs %>%
                ungroup() %>%
                dplyr::rename(TargetGene = gene,
                              p = abc_p) %>%
                prepABCforGSEA(.)
        }
        dfs <- dfs %>%
            group_by(ensembl_gene_id) %>%
            arrange(mean.Score, .by_group=TRUE)
        print("ensembl_gene_ids added, ready for gprofiler!")

        bg_genes <- dfs %>%
            dplyr::slice(1) %>%
            pull(ensembl_gene_id)

        up_genes <- dfs %>%
            dplyr::filter(direction=='up_up') %>%
            dplyr::slice(1) %>%
            pull(ensembl_gene_id)

        down_genes <- dfs %>%
            dplyr::filter(direction=='down_down') %>%
            dplyr::slice(1) %>%
            pull(ensembl_gene_id)

        print(paste0("Testing ",
                     length(up_genes),
                     " 'up' genes against ",
                     length(bg_genes),
                     " genes background..."))

        gres_up <- gprofilerTryCatch(up_genes, bg_genes, sig_thresh=gprof_p)

        print(paste0("Testing ",
                     length(down_genes),
                     " 'down' genes against ",
                     length(bg_genes),
                     " genes background..."))

        gres_down <- gprofilerTryCatch(down_genes, bg_genes, sig_thresh=gprof_p)

        return(list(eres, gres_up, gres_down))

    } else {

        return(eres)

    }



}

DEreinforcingEnhancerTest <- function(df,
                                      p_thresh=0.05,
                                      alt='greater',
                                      mat_out=FALSE) {

    dfs <- df %>% dplyr::filter(p_DE < p_thresh)

    up_up <- dfs %>%
        dplyr::filter(direction=='up_up',
                      beta > 0) %>%
        nrow()

    up_down <- dfs %>%
        dplyr::filter(direction=='up_up',
                      beta < 0) %>%
        nrow()

    down_up <- dfs %>%
        dplyr::filter(direction=='down_down',
                      beta > 0) %>%
        nrow()

    down_down <- dfs %>%
        dplyr::filter(direction=='down_down',
                      beta < 0) %>%
        nrow()

    mat11 <- up_up
    mat21 <- up_down
    mat12 <- down_up
    mat22 <- down_down
    tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                      dimnames =
               list(c('DE up', 'DE down'),
                    c('reinforcing up', 'reinforcing down')))

    ft <- fisher.test(tmat, alternative=alt)

    if (mat_out) {
        return(
            list(
                tmat,
                ft
            )
        )
    }

    eres <- list(de_significance = p_thresh,
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

DE_reinforcing_enhancers <- function(abc,
                                     de,
                                     dist_thresh=NULL,
                                abc_p_thresh=0.01,
                                de_p_col='p_DE',
                                de_p_thresh=0.05,
                                es_prefix='beta',
                                pv_prefix='lfsr',
                                abc_p_col='p_ABC',
                                score_type='ABC.Score',
                                negate_beta=FALSE,
                                all_enh_out=FALSE,
                                run_gsea=FALSE) {

    agg_df <- tibble()

    fe_gprof <- enhancerDirectionTest(abc,
                  dist_thresh=dist_thresh,
                  p_thresh=abc_p_thresh,
                  abc_p_col='p_ABC',
                  abc_lfc_col='log2FoldChange_ABC',
                  gene_col='hgnc_symbol',
                  nEnh=2,
                  score_type=score_type,
                  cdt1 = 'AFR',
                  cdt2 = 'EUR',
                  direction_out = FALSE,
                  mat_out=FALSE,
                  alt='greater',
                  run_gsea=run_gsea,
                  gprof_p=1,
                  debug=FALSE)

    if (run_gsea) {
        curr_df <- as_tibble(fe_gprof[[1]])
    } else {
        curr_df <- as_tibble(fe_gprof)
    }

    curr_df['celltype'] <- 'ABC_overall'
    curr_df['datatype'] <- 'ABC_overall'

    reinf_abc <- enhancerDirectionTest(abc,
                      dist_thresh=dist_thresh,
                      p_thresh=abc_p_thresh,
                      abc_p_col='p_ABC',
                      abc_lfc_col='log2FoldChange_ABC',
                      gene_col='hgnc_symbol',
                      nEnh=2,
                      score_type=score_type,
                      cdt1 = 'AFR',
                      cdt2 = 'EUR',
                      direction_out = TRUE, # return pairs with reinforcing info
                      mat_out=FALSE,
                      alt='greater',
                      run_gsea=FALSE,
                      debug=FALSE)

    dfa <- merge(reinf_abc %>% dplyr::rename(hgnc_symbol = gene),
                 de,
                 all.x=TRUE) %>%
        pivot_longer(cols=matches(paste0(es_prefix,'_|',pv_prefix,'_')),
                     names_to=c(".value","condition","celltype","datatype"),
                     names_pattern="(.*)_(.*)_(.*)_(.*)") %>%
        dplyr::rename(beta := {{es_prefix}},
                      p_DE := {{pv_prefix}})

    if (negate_beta) dfa <- dfa %>% mutate(beta = -beta)

    df <- dfa %>%
        dplyr::filter(!is.na(p_DE)) %>%
        group_by(hgnc_symbol, condition, celltype, datatype) %>%
        arrange(enh_idx, .by_group=TRUE) %>%
        dplyr::slice(1) %>%
        ungroup()

    conditions=df %>% pull(condition) %>% unique()
    celltypes=df %>% pull(celltype) %>% unique()
    datatypes=df %>% pull(datatype) %>% unique()

    baseline_cdt <- 'ABC_overall'
    if ('NI' %in% conditions) baseline_cdt <- 'NI'
    curr_df['condition'] <- baseline_cdt

    agg_df <- rbind(agg_df, curr_df)

    for (cond in conditions) {
        for (ct in celltypes) {
            for (dt in datatypes) {

              dfs <- df %>%
                  dplyr::filter(condition==cond,
                                celltype==ct,
                                datatype==dt)

              if (nrow(dfs)==0) next

              curr_df <- as_tibble(
                  DEreinforcingEnhancerTest(dfs,
                                            p_thresh=de_p_thresh,
                                            alt='greater')
              )
              curr_df['abc_significance'] <- abc_p_thresh
              curr_df['condition'] <- cond
              curr_df['celltype'] <- ct
              curr_df['datatype'] <- dt

              agg_df <- rbind(agg_df, curr_df)

            }

        }

    }

    agg_df['score_type'] <- score_type

    if (all_enh_out) {
        if (run_gsea) {
            return(list(fe_res=agg_df,
                        direction_df=dfa,
                        gres_up=fe_gprof[[2]],
                        gres_down=fe_gprof[[3]]))
        } else {
            return(list(fe_res=agg_df, direction_df=dfa))
        }

    } else {
        return(agg_df)
    }


}

snp_selection_diff_score_fisherwilcox <- function(df,
                                               diff_enh_pairs=FALSE,
                                               snp_stat='iHS',
                                               mat_out=FALSE,
                                               alt='two.sided',
                                               dabc_thresh=0.05,
                                               perc_thresh=0.99,
                                               mid_filt=FALSE, # filter out SNPs w/ 0.5 < iHS percentile < perc_thresh
                                               perc_col='rank2_percentile',
                                               rank_col='rank2_abs_ihs', # column for taking the maximum per enhancer
                                               by_gene=FALSE, # count by gene number instead of enhancer number
                                               all_or_none=FALSE,
                                               run_wilcox=FALSE,
                                               wilcox_col='max_abs_ihs',
                                               return_df=FALSE,
                                               score_col=NULL,
                                               debug=FALSE) {

    # all_or_none requires all SNPs overlapping enhancers
    # or none to be top selection candidates for the gene
    # to be included in the test

    # if diff_enh_pairs==FALSE (not restricted to genes with 2 diff-enhancers) and by_gene==TRUE,
    # then test is on top selection SNP enhancer (diff or non-diff) per gene

    dfs <- df %>%
        dplyr::rename(snp_percentile := {{perc_col}},
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
            mutate(diff_status = case_when(p_ABC < dabc_thresh ~ 'diff',
                                       p_ABC >= 0.5        ~ 'non-diff')) %>%
            dplyr::filter(!is.na(diff_status))
    }

    if (by_gene) {
        dfs <- dfs %>%
            group_by(hgnc_symbol)
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

    if (run_wilcox) {

        wilcox_df <- dfs %>%
            dplyr::rename(wilcox_test_var := wilcox_col)
        wilcox_summ <- wilcox_df %>%
            group_by(diff_status) %>%
            summarize(N = n(),
                      mean_stat = mean(wilcox_test_var))

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

        if (debug) return(wilcox_res)

    }

    if (mid_filt) {
        dfss <- dfs %>% dplyr::filter(selection_candidate | snp_percentile <= 0.5)
    }

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

    if (by_gene) {
        tset <- 'gene-level'
    } else {
        tset <- 'enh-level'
    }
    if (mid_filt) {
        lowstat <- 0.5
    } else {
        lowstat <- perc_thresh
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

    if (run_wilcox) {
        eres <- cbind(as_tibble(eres) %>%
                          dplyr::rename(fishers_P = pvalue),
                      wilcox_res)
    }

    if ('abc_p' %in% colnames(dfs)) dfs <- dfs %>% dplyr::rename(p_ABC = abc_p)
    if ('abc_lfc' %in% colnames(dfs)) dfs <- dfs %>% dplyr::rename(log2FoldChange_ABC = abc_lfc)
    if ('binary_direction' %in% colnames(dfs)) {
        dfs <- dfs %>% mutate(enh_status = binary_direction)
    } else {
        dfs <- dfs %>% mutate(enh_status = diff_status)
    }
    dfs <- dfs %>% mutate(snp_stat = snp_stat, enh_type = enhtype, test_set = tset)

    plotcols <- c('seqnames','start','end','hgnc_symbol','ensembl_gene_id',
                  'entrez_id','class','TargetGeneTSS','distance','isSelfPromoter',
                  'name','log2FoldChange_ABC','p_ABC','neglog10p_ABC','RSNUM',
                  'rank_var','snp_percentile','selection_candidate','enh_status',
                  wilcox_col,'snp_stat','enh_type','test_set')

    dfs <- dfs %>% dplyr::select(any_of(plotcols))

    if (!is.null(score_col)) {
        dfs <- dfs %>% mutate(score_column = score_col)
        eres <- eres %>% mutate(score_column = score_col)
    }

    if (return_df) {
        return(list(fw=eres, plotdf=dfs))
    } else {
        return(eres)
    }


}

parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--dist_thresh',
                    type='double',
                    default=0,
                    help='Require reinforcing/opposing enhancers to be at least this far apart (default of 0 means ~500 bp enhancers can be separated by 1 bp)')
parser$add_argument('--FST',
                    type='character',
                    default=NULL,
                    help="File of SNPs with FST values calculate between 4 pops of each ancestry.")
parser$add_argument('--iHS',
                    type='character',
                    default=NULL,
                    help="File of SNPs with iHS scores for at least 2 pops and percentile of |iHS| for second highest |iHS| per SNP (`rank2_percentile`).")
parser$add_argument('--score_column',
                    type='character',
                    default='ABC.Score',
                    help="Type of score to process and specify in score_type column for aggregating later.")
parser$add_argument('--facetscales',
                    type='character',
                    default='fixed',
                    help="Fix axes scales across facets or allow to adapt to data with 'free', 'free_x', or 'free_y'")
parser$add_argument('--facetvar',
                    type='character',
                    default=NULL,
                    help="Set to 'datatype' if DE has multiple datatypes to make separate enrichment plots for")
parser$add_argument('--enrichgroupcolors',
                    type='character',
                    default='gray,darkred',
                    help="Colors to distinguish enrichment groups.")
parser$add_argument('--enrichgroup',
                    type='character',
                    default=NULL,
                    help="Variable to group by/color on in enrichment plots.")
parser$add_argument('--enrichvarlab',
                    type='character',
                    default=NULL,
                    help="Y-axis label for enrichment plots.")
parser$add_argument('--enrichvar',
                    type='character',
                    default='condition',
                    help="Variable to plot on y-axis of enrichment plots.")
parser$add_argument('--plotdir',
                    type='character',
                    default=NULL, # /home/kpettie/code/github/plotting
                    help="Directory with plotting function in 'plotting.R'. Required to make plots")
parser$add_argument('--pval_prefix',
                    type='character',
                    default='pval',
                    help='prefix of columns with p-values')
parser$add_argument('--ES_prefix',
                    type='character',
                    default='beta',
                    help='prefix of columns with effect sizes')
parser$add_argument('--DE',
                    type='character',
                    default='./ancestry_DESeq_Geuvadis_LRT.txt.gz',
                    help='File with differential expression results')
parser$add_argument('--diffABC_preds',
                    type='character',
                    default='./meanQN.16.AFR_EUR.diff.ABC.Score.txt.gz',
                    help='File with diff-ABC predictions from t-test or other model')
parser$add_argument('--DE_thresh',
                    type='double',
                    default=0.05,
                    help='DE p-value for defining DE genes')
parser$add_argument('--diffABC_thresh',
                    type='double',
                    default=0.05,
                    help='diffABC P-value threshold for testing DE enrichment')
parser$add_argument('--top_diffABC',
                    action='store_true',
                    default=FALSE,
                    help="Subset analysis to top diff-ABC E-G pair per gene")
parser$add_argument('-n', '--name',
                    type='character',
                    default='diffABC_DE_overlap',
                    help="Basename for output files")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

dist_thresh <- opt$dist_thresh
ihs_fname <- opt$iHS
fst_fname <- opt$FST
de_fname <- opt$DE
dabc_fname <- opt$diffABC_preds
de_thresh <- opt$DE_thresh
dabc_thresh <- opt$diffABC_thresh
outbase <- opt$name
outdir <- opt$outdir
top_dabc <- opt$top_diffABC
pv_pre <- opt$pval_prefix
es_pre <- opt$ES_prefix
plotdir <- opt$plotdir
enrichvar <- opt$enrichvar
enrichvarlab <- opt$enrichvarlab
enrichgroup <- opt$enrichgroup
colorvals <- str_split(opt$enrichgroupcolors, ',')[[1]]
facetvar <- opt$facetvar
if (!is.null(facetvar)) facetvar <- str_split(facetvar, ',')[[1]]
facetscales <- opt$facetscales
score_column <- opt$score_column

negate_beta <- FALSE
if (grepl('RandolphEtAl2021_popDE', de_fname)) negate_beta <- TRUE

cdt1 <- 'AFR'
cdt2 <- 'EUR'

dabc <- prepABCforGSEA(read_tsv(dabc_fname)) %>%
    dplyr::select(-c(.y.,group1,group2,n1,n2)) %>%
    mutate(neglog10p_ABC = -log10(p)) %>%
    dplyr::rename(log2FoldChange_ABC = log2FoldChange,
                  stat_ABC = statistic,
                  df_ABC = df,
                  p_ABC = p)

if (grepl('chip', tolower(score_column)) | grepl('hic', tolower(score_column))) {

    dabc <- dabc %>%
        group_by(hgnc_symbol, log2FoldChange_ABC, p_ABC) %>%
        arrange(desc(isSelfPromoter), class, start, .by_group=TRUE) %>%
        mutate(start = min(start),
               end = max(end),
               name = paste0('HiChIP|',chr,':',start,'-',end)) %>%
        dplyr::slice(1) %>%
        ungroup()

}

de <- read_tsv(de_fname)
if (pv_pre=='p') {
    de <- de %>%
        rename_with(.,
                    ~ gsub("p_","pval_",.x),
                    .cols=tidyselect::starts_with('p_'))
    pv_pre <- 'pval'
}

# for overall enrichment/directionality matching
abcde <- merge(dabc, de) %>%
    pivot_longer(cols=matches(paste0(es_pre,'_|',pv_pre,'_')),
                 names_to=c(".value","condition","celltype","datatype"),
                 names_pattern="(.*)_(.*)_(.*)_(.*)") %>%
    dplyr::rename(beta := {{es_pre}},
                  p_DE := {{pv_pre}})

write_tsv(abcde, file.path(outdir,paste0(outbase,'.txt')))

if (!is.null(ihs_fname)) {

    include_cols <- c('seqnames',
                      'start',
                      'RSNUM',
                      'max_abs_ihs',
                      'rank2_abs_ihs',
                      'rank2_percentile',
                      'selection_candidate')

    ihs <- read_tsv(ihs_fname, col_select=any_of(include_cols))

    dabc <- join_overlap_left(
                dabc %>% dplyr::rename(seqnames=chr) %>% as_granges(),
                ihs %>% mutate(end = start + 1) %>% as_granges()
            ) %>%
            as_tibble() %>%
            group_by(hgnc_symbol,name) %>%
            arrange(desc(rank2_abs_ihs), .by_group=TRUE) %>%
            dplyr::slice(1) %>%
            ungroup()

}

if (!is.null(fst_fname)) {

    fst <- read_tsv(fst_fname)
    ### CHECK FST FILE FORMATTING
    ### ADD PERCENTILE COLUMN

    dabc <- join_overlap_left(
            dabc %>% as_granges(),
            ihs %>%
                dplyr::rename(start = POS) %>%
                mutate(end = start + 1) %>%
                as_granges()
        ) %>%
        as_tibble() %>%
        group_by(hgnc_symbol,name) %>%
        arrange(desc(FST), .by_group=TRUE) %>%
        dplyr::slice(1) %>%
        ungroup()

}

# for directionality matching among reinforcing enhancers
fe_dirdf <- DE_reinforcing_enhancers(dabc,
                                    de,
                                    dist_thresh=dist_thresh,
                                    abc_p_thresh=0.01,
                                    de_p_col='p_DE',
                                    de_p_thresh=0.1,
                                    es_prefix=es_pre,
                                    pv_prefix=pv_pre,
                                    abc_p_col='p_ABC',
                                    score_type=score_column,
                                    negate_beta=negate_beta,
                                    all_enh_out=TRUE,
                                    run_gsea=TRUE)
fe_sign <- fe_dirdf$fe_res
write_tsv(fe_sign,
          file.path(outdir,
                    paste0(outbase,
                           '.enhReinforcingDirection.fisherEnrichments.txt.gz')))
dirdf <- fe_dirdf$direction_df
write_tsv(dirdf,
          file.path(outdir,
                    paste0(outbase,
                           '.enhReinforcingDirection.allPairs.txt.gz')))

if (!is.null(fe_dirdf$gres_up)) {

    write_tsv(fe_dirdf$gres_up,
              file.path(outdir,
                        paste0(outbase,
                               '.enhReinforcingDirection.gprofiler',
                               cdt2,
                               '.txt.gz')))

}

if (!is.null(fe_dirdf$gres_down)) {

  write_tsv(fe_dirdf$gres_down,
           file.path(outdir,
                     paste0(outbase,
                            '.enhReinforcingDirection.gprofiler',
                            cdt1,
                            '.txt.gz')))

}

aggt <- tibble()

if (!is.null(ihs_fname)) {

  t1 <- snp_selection_diff_score_fisherwilcox(dabc,
                                        snp_stat='iHS',
                                        diff_enh_pairs=FALSE,
                                 mat_out=FALSE,
                                 alt='greater',
                                 dabc_thresh=0.05,
                                 perc_thresh=0.95,
                                 mid_filt=TRUE, # filter out SNPs w/ 0.5 < iHS percentile < perc_thresh
                                 perc_col='rank2_percentile',
                                 by_gene=TRUE, # count by gene number instead of enhancer number
                                 all_or_none=FALSE,
                                 run_wilcox=TRUE,
                                 wilcox_col='max_abs_ihs',
                                 return_df=FALSE,
                                 score_col=score_column,
                                 debug=FALSE)

  t2 <- snp_selection_diff_score_fisherwilcox(dabc,
                                        snp_stat='iHS',
                                        diff_enh_pairs=FALSE,
                                 mat_out=FALSE,
                                 alt='greater',
                                 dabc_thresh=0.05,
                                 perc_thresh=0.95,
                                 mid_filt=TRUE, # filter out SNPs w/ 0.5 < iHS percentile < perc_thresh
                                 perc_col='rank2_percentile',
                                 by_gene=FALSE, # count by gene number instead of enhancer number
                                 all_or_none=FALSE,
                                 run_wilcox=TRUE,
                                 wilcox_col='max_abs_ihs',
                                 return_df=FALSE,
                                 score_col=score_column)

  t3 <- snp_selection_diff_score_fisherwilcox(dirdf,
                                      snp_stat='iHS',
                                        diff_enh_pairs=TRUE,
                                 mat_out=FALSE,
                                 alt='greater',
                                 dabc_thresh=0.01,
                                 perc_thresh=0.95,
                                 mid_filt=TRUE, # filter out SNPs w/ 0.5 < iHS percentile < perc_thresh
                                 perc_col='rank2_percentile',
                                 by_gene=TRUE, # count by gene number instead of enhancer number
                                 all_or_none=FALSE,
                                 run_wilcox=TRUE,
                                 wilcox_col='max_abs_ihs',
                                 return_df=FALSE,
                                 score_col=score_column)

  t4 <- snp_selection_diff_score_fisherwilcox(dirdf,
                                      snp_stat='iHS',
                                        diff_enh_pairs=TRUE,
                                 mat_out=FALSE,
                                 alt='greater',
                                 dabc_thresh=0.01,
                                 perc_thresh=0.95,
                                 mid_filt=TRUE, # filter out SNPs w/ 0.5 < iHS percentile < perc_thresh
                                 perc_col='rank2_percentile',
                                 by_gene=FALSE, # count by gene number instead of enhancer number
                                 all_or_none=FALSE,
                                 run_wilcox=TRUE,
                                 wilcox_col='max_abs_ihs',
                                 return_df=FALSE,
                                 score_col=score_column)
  aggt <- rbind(t1,t2,t3,t4)

}

if (!is.null(fst_fname)) {

  t5 <- snp_selection_diff_score_fisherwilcox(dabc,
                                        snp_stat='FST',
                                        diff_enh_pairs=FALSE,
                                 mat_out=FALSE,
                                 alt='greater',
                                 dabc_thresh=0.05,
                                 perc_thresh=0.95,
                                 mid_filt=TRUE, # filter out SNPs w/ 0.5 < iHS percentile < perc_thresh
                                 perc_col='percentile',
                                 by_gene=TRUE, # count by gene number instead of enhancer number
                                 all_or_none=FALSE,
                                 run_wilcox=TRUE,
                                 wilcox_col='FST',
                                 return_df=FALSE,
                                 score_col=score_column,
                                 debug=FALSE)

  t6 <- snp_selection_diff_score_fisherwilcox(dabc,
                                        snp_stat='FST',
                                        diff_enh_pairs=FALSE,
                                 mat_out=FALSE,
                                 alt='greater',
                                 dabc_thresh=0.05,
                                 perc_thresh=0.95,
                                 mid_filt=TRUE, # filter out SNPs w/ 0.5 < iHS percentile < perc_thresh
                                 perc_col='percentile',
                                 by_gene=FALSE, # count by gene number instead of enhancer number
                                 all_or_none=FALSE,
                                 run_wilcox=TRUE,
                                 wilcox_col='FST',
                                 return_df=FALSE,
                                 score_col=score_column)

  t7 <- snp_selection_diff_score_fisherwilcox(dirdf,
                                      snp_stat='FST',
                                        diff_enh_pairs=TRUE,
                                 mat_out=FALSE,
                                 alt='greater',
                                 dabc_thresh=0.01,
                                 perc_thresh=0.95,
                                 mid_filt=TRUE, # filter out SNPs w/ 0.5 < iHS percentile < perc_thresh
                                 perc_col='percentile',
                                 by_gene=TRUE, # count by gene number instead of enhancer number
                                 all_or_none=FALSE,
                                 run_wilcox=TRUE,
                                 wilcox_col='FST',
                                 return_df=FALSE,
                                 score_col=score_column)

  t8 <- snp_selection_diff_score_fisherwilcox(dirdf,
                                      snp_stat='FST',
                                        diff_enh_pairs=TRUE,
                                 mat_out=FALSE,
                                 alt='greater',
                                 dabc_thresh=0.01,
                                 perc_thresh=0.95,
                                 mid_filt=TRUE, # filter out SNPs w/ 0.5 < iHS percentile < perc_thresh
                                 perc_col='percentile',
                                 by_gene=FALSE, # count by gene number instead of enhancer number
                                 all_or_none=FALSE,
                                 run_wilcox=TRUE,
                                 wilcox_col='FST',
                                 return_df=FALSE,
                                 score_col=score_column)

  aggt <- rbind(aggt,t5,t6,t7,t8)

}

if (nrow(aggt)>0) {

    aggt %>%
        write_tsv(., file.path(outdir,paste0(outbase,'.selectionSNPs.fisherEnrichmentsWilcox.txt')))

}

if (top_dabc) {

  # top diff-CRE per gene may be same CRE for multiple nearby genes
  abcde <- abcde %>%
    group_by(ensembl_gene_id, condition, celltype, datatype) %>%
    arrange(p_ABC, p_DE, .by_group=TRUE) %>%
    dplyr::slice(1) %>%
    ungroup()

}

fer <- fisher_enrichments(abcde) %>%
    mutate(score_type = score_column)
write_tsv(fer, file.path(outdir,paste0(outbase,'.fisherEnrichments.txt')))

ferd <- fisher_enrichments_direction(abcde, negate_beta=negate_beta) %>%
    mutate(score_type = score_column)
write_tsv(ferd, file.path(outdir,paste0(outbase,'.fisherEnrichmentsDirection.txt')))

# for testing DE directionality matching among top diff-score genes per condition
topN <- fe_sign %>%
    dplyr::filter(datatype != 'ABC_overall') %>%
    rowwise() %>%
    mutate(nGenes = sum(c_across(contains(c('up_','down_'))), na.rm=TRUE)) %>%
    ungroup() %>%
    pull(nGenes) %>%
    max()
topabcde <- abcde %>%
    dplyr::filter(p_DE < de_thresh) %>%
    group_by(condition, celltype, datatype) %>%
    arrange(p_ABC, .by_group=TRUE) %>%
    dplyr::slice(1:topN) %>%
    ungroup()
topferd <- fisher_enrichments_direction(topabcde, negate_beta=negate_beta) %>%
    mutate(score_type = score_column,
           nTop = topN)
write_tsv(topferd, file.path(outdir,paste0(outbase,'.fisherEnrichmentsDirectionTopN.txt')))

if (!is.null(plotdir)) {

  source(file.path(plotdir, 'plotting.R'))

  if (!is.null(enrichgroup)) {
    if (n_distinct(fer[[enrichgroup]])!=length(colorvals)) {
      colorvals <- NULL
    }
  } else {
    colorvals <- NULL
  }

  if (!is.null(facetvar)) {
    ewidth <- 15
    eheight <- 7
    cwidth <- 32
    cheight <- 26
  } else {
    ewidth <- 7
    eheight <- 7
    cwidth <- 19
    cheight <- 12
  }

  p1 <- plotFisherEnrichments(fer,
                      yvar_col=enrichvar,
                      ylabel=enrichvarlab,
                      althypoth='greater',
                      groupvar=enrichgroup,
                      # groupvarlab=groupvar,
                      groupvarord=NULL,
                      colorvals=colorvals,
                      w=ewidth,
                      h=eheight,
                      facets=facetvar,
                      frows=1,
                      fscales=facetscales)
  ggsave(file.path(outdir,paste0(outbase,'.fisherEnrichments.png')),
         width=ewidth,
         height=eheight)

  p2 <- plotFisherEnrichments(ferd,
                      yvar_col=enrichvar,
                      ylabel=enrichvarlab,
                      althypoth='greater',
                      groupvar=enrichgroup,
                      # groupvarlab=groupvar,
                      groupvarord=NULL,
                      colorvals=colorvals,
                      w=ewidth,
                      h=eheight,
                      facets=facetvar,
                      frows=1,
                      fscales=facetscales)
  ggsave(file.path(outdir,paste0(outbase,'.fisherEnrichmentsDirection.png')),
         width=ewidth,
         height=eheight)

  p4 <- plotFisherEnrichments(fe_sign,
                      yvar_col=enrichvar,
                      ylabel=enrichvarlab,
                      althypoth='greater',
                      groupvar=enrichgroup,
                      # groupvarlab=groupvar,
                      groupvarord=NULL,
                      colorvals=colorvals,
                      w=ewidth,
                      h=eheight,
                      facets=facetvar,
                      frows=1,
                      fscales=facetscales)
  ggsave(file.path(outdir,paste0(outbase,'.enhReinforcingDirection.fisherEnrichments.png')),
         width=ewidth,
         height=eheight)

  corrfacets <- enrichvar
  if (!is.null(enrichgroup)) corrfacets <- c(corrfacets, enrichgroup)
  if (!is.null(facetvar)) corrfacets <- c(corrfacets, facetvar)

  corrfrows <- 3
  if (length(corrfacets)>2) corrfrows <- 4

  if (negate_beta) abcde <- abcde %>% mutate(beta = -beta)

  p3 <- plotCorrelation(abcde %>%
                          dplyr::filter(p_ABC<dabc_thresh,
                                        p_DE<de_thresh),
                        "log2FoldChange_ABC",
                        "beta",
                        xlabel = "ABC L2FC",
                        ylabel = "DE Effect Size",
                        groupvar = "isSelfPromoter",
                        colorvals = c("gray", "red"),
                        addline = c(0,1),
                        legendpos = "right",
                        corrpos = NULL,
                        corrcolor = "black",
                        w = cwidth,
                        h = cheight,
                        logscale = FALSE,
                        facets = corrfacets,
                        frows = corrfrows,
                        orderpval=TRUE,
                        fscales='free',
                        addaxes=TRUE)
  ggsave(file.path(outdir,paste0(outbase,'.effectSizeCorrelations.png')),
         width=cwidth,
         height=cheight)

}
