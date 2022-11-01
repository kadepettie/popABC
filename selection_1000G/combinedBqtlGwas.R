#!/usr/bin/R

library(argparse)
library(ggplot2)
library(tidyverse)



hypergeomDiffQTLgwas <- function(df,
                                 abc_p_col='chip_p',
                                 abc_p_thresh=0.05,
                                 abc_p_nonthresh=0.05,
                                 gwas_thresh=5e-5,
                                 gwas_nonthresh=5e-5,
                                 hits_out=FALSE,
                                 mat_out=FALSE,
                                 alt='greater') {

    dfs <- df %>%
        dplyr::rename(abc_p := {{abc_p_col}}) %>%
        # don't double count CREs
        group_by(seqnames,start,end,name) %>%
        arrange(abc_p) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        # don't double count tested GWAS variant
        # take top diff-CRE per GWAS var
        group_by(seqnames, gwas_pos) %>%
        arrange(abc_p) %>%
        dplyr::slice(1) %>%
        ungroup()

    if (grepl('chip', abc_p_col) | grepl('hic', abc_p_col)) {
      # only count each HiChIP bin once
      dfs <- dfs %>%
            group_by(hgnc_symbol, chip_lfc, abc_p) %>%
            mutate(chipbin_has_promoter = if_else(any(class=='promoter'), TRUE, FALSE))

      # take the CRE linked to a trait, then bQTL in LD most significant GWAS hit in the HiChIP bin
      dfs <- dfs %>%
            arrange(gwas_pval, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup()

    }

    # eliminate middle range pvalues, if applicable
    dfs <- dfs %>%
        filter(abc_p < abc_p_thresh | abc_p >= abc_p_nonthresh,
               gwas_pval < gwas_thresh | gwas_pval >= gwas_nonthresh)

    mat11 <- dfs %>% filter(abc_p < abc_p_thresh, gwas_pval < gwas_thresh) %>% nrow()
    mat21 <- dfs %>% filter(abc_p < abc_p_thresh, gwas_pval >= gwas_nonthresh) %>% nrow()
    mat12 <- dfs %>% filter(abc_p >= abc_p_nonthresh, gwas_pval < gwas_thresh) %>% nrow()
    mat22 <- dfs %>% filter(abc_p >= abc_p_nonthresh, gwas_pval >= gwas_nonthresh) %>% nrow()
    tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                      dimnames =
               list(c('GWAS hit', 'not GWAS hit'),
                    c('diff-CRE bQTL', 'non-diff-CRE bQTL')))

    ft <- fisher.test(tmat, alternative=alt)

    if (mat_out) {
        return(
            list(
                tmat,
                ft
            )
        )
    }

    eres <- list(diff_significance = abc_p_thresh,
                 gwas_significance = gwas_thresh,
                 odds_ratio = ft$estimate[[1]],
                 pvalue = ft$p.value,
                 conf_lower = ft$conf.int[1],
                 conf_upper = ft$conf.int[2],
                 up_up = mat11,
                 up_down = mat21,
                 down_up = mat12,
                 down_down = mat22,
                 method = 'Fisher',
                 alternative = alt)

    if (hits_out) {
        dfsout <- dfs %>%
            filter(gwas_pval < gwas_thresh) %>%
            dplyr::rename(!!abc_p_col := abc_p)
        return(list(as_tibble(eres), dfsout))
    }

    return(eres)

}

loopHypergeomDiffQTLgwas <- function(df) {

    ares <- tibble()
    adfs <- tibble()
    sts <- df %>% pull(score_type) %>% sort() %>% unique()
    for (st in sts) {

        cres <- hypergeomDiffQTLgwas(df %>% filter(score_type==st),
                                 abc_p_col='p',
                                 abc_p_thresh=0.05,
                                 abc_p_nonthresh=0.05,
                                 gwas_thresh=5e-5,
                                 gwas_nonthresh=5e-5,
                                 hits_out=TRUE,
                                 mat_out=FALSE,
                                 alt='greater') %>%
            imap(~ mutate(.x, score_type=st))

        ares <- rbind(ares,
                      cres[[1]])

        adfs <- rbind(adfs,
                      cres[[2]])

    }

    resdfs <- list(ares %>%
                       mutate(test_type='gwasDiffCRE'),
                   adfs)

    return(resdfs)

}

mergeLDbqtl <- function(lddir='.',ldpatt="expandLD.txt.gz") {

    ldbqtl <- tibble()

    for (cbqtl_fname in list.files(lddir, pattern=ldpatt)) {

      qtl_name <- str_split(cbqtl_fname,
                            pattern='\\.',
                            simplify=TRUE)[,8]
      cbqtl <- read_tsv(file.path(lddir, cbqtl_fname)) %>%
          mutate(bqtl_name = qtl_name)

      if ('ld_snps' %in% colnames(cbqtl)) {
          cbqtl <- cbqtl %>%
              select(-ld_snps) %>%
              rename(gwas_pos = ld_pos) %>%
              mutate(gwas_pos = if_else(is.na(gwas_pos), as.double(bqtl_pos), as.double(gwas_pos)))
      }

      ldbqtl <- rbind(ldbqtl, cbqtl)

    }

    extrafreqcols <- c('Depth',
                       'ALTdepth',
                       'REFDepth',
                       'ALTallele',
                       'POSTallele',
                       'POSTfreq',
                       'prechipfreq',
                       'eqtl_lcl_AFR_alt_frq',
                       'eqtl_lcl_EUR_alt_frq',
                       'eqtl_pbmc_AFR_alt_frq',
                       'eqtl_pbmc_EUR_alt_frq',
                       'eqtl_lcl_alt_high_freq',
                       'eqtl_lcl_alt_high_freq_all',
                       'eqtl_pbmc_alt_high_freq')

    # keeps multiple target genes for each CRE for checking potential GWAS hit
    # target genes later
    ldbqtl <- ldbqtl %>%
        group_by(seqnames,
                 start,
                 end,
                 hgnc_symbol,
                 name,
                 chip_lfc,
                 atac_lfc,
                 chip_p,
                 atac_p,
                 bqtl_pos,
                 tag_rsID,
                 gwas_pos,
                 bqtl_direction,
                 bqtl_fst,
                 bqtl_fst_percentile) %>%
        arrange(bqtl_p, .by_group=TRUE) %>%
        mutate(min_p = pmin(chip_p, atac_p, na.rm=T),
               bqtl_names = paste0(bqtl_name, collapse=';')) %>%
        arrange(bqtl_p, min_p, .by_group=TRUE) %>%
        slice(1) %>%
        ungroup() %>%
        dplyr::select(-any_of(extrafreqcols)) %>%
        rename(bqtl_ref = ref,
               bqtl_alt = alt,
               alt_high_freq = bqtl_alt_high_freq,
               cre_max_fst = FST) %>%
        rename_with(.cols = contains(c('frq','freq')), ~ paste0('bqtl_', .x)) %>%
        as_tibble()

    return(ldbqtl)

}

gwasFstAlleleFreqMerge <- function(dfa,
                                   fst_dir='.',
                                   fst_pattern='fst.bed',
                                   fst_prefix='EUR_vs_AFR',
                                   frq_dir='.',
                                   aprefix='AFR',
                                   eprefix='EUR',
                                   debug_chr=NULL) {

    chrnames <- c(1:22, 'X')
    adfm <- tibble()

    if (!is.null(debug_chr)) chrnames <- c(debug_chr)

    for (currchr in chrnames) {

        print(paste0('Processing chr', currchr))
        print('    splitting variant column')
        # tested GWAS variants per chromosome should be fewer than variants with allele freq or FST info
        df <- dfa %>%
            filter(grepl(paste0('^',currchr,':.*'), variant)) %>%
            separateFast(.,
                         'variant',
                         c('seqnames','gwas_pos','ref','alt'),
                         sep=":",
                         remove=FALSE,
                         coltypes=c('cdcc')) %>%
            mutate(seqnames=paste0('chr',seqnames),
                   RiskAllele = if_else(beta > 0, alt, ref))

        currchr <- paste0('chr', currchr)


        # read allele frequencies
        frqa_fname <- file.path(frq_dir, paste(aprefix, currchr, sep='.'))
        frqe_fname <- file.path(frq_dir, paste(eprefix, currchr, sep='.'))

        if (file.exists(frqa_fname) & file.exists(frqe_fname)) {
            # outputs position col as 'gwas_pos' for overlap
            print(paste0('    Reading/merging ', aprefix, ' and ', eprefix, ' allele frequencies'))
            cfrq <- alleleFreqMergeFast(frqa_fname, frqe_fname)

            # merge with GWAS
            print('    Merging AF with GWAS and splitting AF cols')
            dfm <- merge(df, cfrq, all.x=TRUE) %>%
                separateFast(.,
                             'ref_frqa',
                             c('frqref','AFR_ref_frq'),
                             sep=":(?=[0-9])",
                             remove=TRUE,
                             coltypes=c('cd')) %>%
                separateFast(.,
                             'alt_frqa',
                             c('frqalt','AFR_alt_frq'),
                             sep=":(?=[0-9])",
                             remove=TRUE,
                             coltypes=c('cd')) %>%
                separateFast(.,
                             'ref_frqe',
                             c('frqref2','EUR_ref_frq'),
                             sep=":(?=[0-9])",
                             remove=TRUE,
                             coltypes=c('cd')) %>%
                separateFast(.,
                             'alt_frqe',
                             c('frqalt2','EUR_alt_frq'),
                             sep=":(?=[0-9])",
                             remove=TRUE,
                             coltypes=c('cd')) %>%
                mutate(alt_high_freq = case_when(EUR_alt_frq > AFR_alt_frq ~ 'EUR',
                                                 EUR_alt_frq < AFR_alt_frq ~ 'AFR')) %>%
                select(-any_of(c('frqref2','frqalt2')))

        } else {
            print(paste0('    Missing ', currchr, ' allele frequency file'))

            dfm <- df %>%
                mutate(frqref = NA,
                       frqalt = NA,
                       AFR_ref_frq = NA,
                       AFR_alt_frq = NA,
                       EUR_ref_frq = NA,
                       EUR_alt_frq = NA,
                       alt_high_freq = NA)
        }


        # read FST
        cfst_fname <- file.path(fst_dir, paste(fst_prefix, currchr, fst_pattern, sep='.'))

        if (file.exists(cfst_fname)) {
            print(paste0('    Reading ', fst_prefix, ' FST'))
            cfst <- read_tsv(cfst_fname, col_names=c('seqnames','start','end','FST')) %>%
                select(-start) %>%
                rename(gwas_pos = end)

            print('    Merging FST with GWAS')
            dfm <- merge(dfm, cfst, all.x=TRUE)

        } else {
            print(paste0('    Missing ', currchr, ' FST file'))

            dfm <- dfm %>%
                mutate(FST = NA)
        }


        adfm <- rbind(adfm, dfm)

    }

    return(adfm)

}

getPhenoCodeDescription <- function(mani_fname, gwas_fname) {

    if (is.null(mani_fname)) return(list(NA, NA))

    print(paste0('Getting phenotype code and description for ', base::basename(gwas_fname)))

    mani <- read_tsv(mani_fname, na = c("", "NA","N/A")) %>%
        rename_with(., ~ gsub(' ', '_', tolower(.x)))

    fnamesplit <- str_split(base::basename(gwas_fname), "\\.", simplify=TRUE)
    pcode <- fnamesplit[1]
    sexfilt <- fnamesplit[4]

    if (!(sexfilt %in% c('both_sexes', 'female', 'male'))) sexfilt <- 'both_sexes'

    prow <- mani %>% filter(phenotype_code==pcode, sex==sexfilt)
    if (nrow(prow)==0) return(list(pcode, NA))
    pdesc <- prow$phenotype_description

    pdescsplit <- str_split(pdesc, ": ", simplify=TRUE)
    pdesc <- pdescsplit[ncol(pdescsplit)]

    return(list(pcode, pdesc))

}

narrow_diffABC <- function(df, displayformat=FALSE) {

    #' Pivot longer atac/chip-widened diff-ABC df
    #'
    #' @description Removes all atac/chip columns except
    #' `lfc` and `p`, then pivots these longer adding
    #' `score_type` column to distinguish

    df <- df %>%
        select(-contains('atac'), atac_lfc, atac_p) %>%
        select(-contains('chip'), chip_lfc, chip_p) %>%
        pivot_longer(cols=contains(c('chip','atac')),
                     names_to=c("score_type", ".value"),
                     names_pattern="(.*)_(.*)")

    if (displayformat) {
        df <- df %>% mutate(score_type=case_when(score_type=='chip' ~ 'ChIP',
                                                score_type=='atac' ~ 'ATAC',
                                                TRUE ~ score_type))
    }

    return(df)

}



parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--ukbbmani_fname',
                    type='character',
                    default=NULL,
                    help='UK biobank manifest file for extracting phenotype description matching code in results file name')
parser$add_argument('--fst_prefix',
                    type='character',
                    default='EUR_vs_AFR')
parser$add_argument('--eprefix',
                    type='character',
                    default='EUR')
parser$add_argument('--aprefix',
                    type='character',
                    default='AFR')
parser$add_argument('--frq_dir',
                    type='character',
                    default='.')
parser$add_argument('--frq_pattern',
                    type='character',
                    default=NULL,
                    help='Pattern for matching files with allele frequency split by chromosome.')
parser$add_argument('--pval_gwas',
                    type='double',
                    default=0.5e-5)
parser$add_argument('--fst_dir',
                    type='character',
                    default='.')
parser$add_argument('--d_dir',
                    type='character',
                    default='.')
parser$add_argument('--gw_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--d_pattern',
                    type='character',
                    default="expandLD.txt.gz",
                    help="Pattern for matching multiple LD-expanded, CRE-bQTL overlaps in working/output directory")
parser$add_argument('--fst_pattern',
                    type='character',
                    default='fst.bed',
                    help='Pattern for matching files with FST split by chromosome.')
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
                    default='diffscore_QTL_overlap',
                    help="Basename for output files")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

pval_gwas <- opt$pval_gwas
fst_dir <- opt$fst_dir
fst_pattern <- opt$fst_pattern
fst_prefix <- opt$fst_prefix
frq_dir <- opt$frq_dir
aprefix <- opt$aprefix
eprefix <- opt$eprefix
d_dir <- opt$d_dir
d_pattern <- opt$d_pattern
gw_fname <- opt$gw_fname
ukbbmani_fname <- opt$ukbbmani_fname
outdir <- opt$outdir

# set to maximize power/data usage
pval_nongwas <- pval_gwas

pheno_att <- getPhenoCodeDescription(ukbbmani_fname, gw_fname)
pheno_code <- pheno_att[[1]]
pheno_desc <- pheno_att[[2]]

source(file.path(opt$plotdir,'plotting.R'))
source(file.path(this.path::this.dir(default='~/code/github/selection_1000G'), 'variantTools.R'))

# merge LD-expanded bQTL
# multiple target genes per CRE sliced to one before GWAS enrichment test
ldbqtls <- mergeLDbqtl(lddir=d_dir,ldpatt=d_pattern)

# read UKBB GWAS phenotype results
cukba <- read_tsv(gw_fname)

gw_fname <- base::basename(gw_fname)

# merge FST and allele freq with gwas
gwm <- gwasFstAlleleFreqMerge(cukba,
                               fst_dir=fst_dir,
                               fst_pattern=fst_pattern,
                               fst_prefix=fst_prefix,
                               frq_dir=frq_dir,
                               aprefix=aprefix,
                               eprefix=eprefix)

gwm %>%
    write_tsv(
      file.path(
        outdir,
        paste0(
          gw_fname,
          '.',
          fst_prefix,
          'alleleFreqFST.txt.gz'
        )
      )
    )

# merge LD bQTL with gwas and write out for looking up hits later
print("Merging GWAS with LD-expanded bQTL...")
gbq <- merge(gwm %>%
                 rename_with(.cols = -all_of(c('seqnames','gwas_pos')),
                             ~ paste0('gwas_', .x)),
             ldbqtls) %>%
    # remove LD expansions of bQTL in LD with GWAS that did not directly overlap GWAS
    group_by(hgnc_symbol,seqnames,start,end,name,bqtl_pos,tag_rsID) %>%
    arrange(gwas_pval, .by_group=TRUE) %>%
    dplyr::slice(1) %>%
    ungroup()

bqtl_names <- gbq %>% pull(bqtl_name) %>% sort() %>% unique
bqtloutstr <- paste0(bqtl_names, collapse='_')

gbq %>%
    write_tsv(
      file.path(
        outdir,
        paste0(
          gw_fname,
          '.',
          bqtloutstr,
          '_CREs.txt.gz'
        )
      )
    )

gbql <- narrow_diffABC(gbq, displayformat=TRUE)
hglist <- loopHypergeomDiffQTLgwas(gbql) %>%
    imap(~ mutate(.x, phenotype_code=pheno_code,
                      phenotype_description=pheno_desc))


# wilcoxon test for GWAS hit FST enrichment relative to all tested non-hits
# (on full GWAS table irrespective of CREs)
gw_wres <- groupedBoxplot(gwm %>%
                   filter(!is.na(FST),
                          !is.na(pval)) %>%
                   mutate(variant_type=case_when(pval >= pval_nongwas ~ 'non-hit',
                                                 pval < pval_gwas ~ 'hit'),
                          variant_type=fct_relevel(variant_type, rev)),
                           'variant_type',
                           'FST',
                           xlab='Variant type',
                           groupvar=NULL,
#                            xlab=xvar,
#                            ylab=yvar,
                           plotpoints=FALSE,
                           points_thresh=200, # show all points for xvar with N observations below this number
                           show_wilcox=TRUE,
                           paired_t_test=FALSE,
                          testout=TRUE,
                           alt='greater',
                           orderpval=FALSE,
                           logscale=FALSE,
                           colorX=TRUE,
                           colorvals=c('darkgrey','red3'),
                           showleg=FALSE,
                           legendpos='right',
                           horizontal=TRUE,
                           facets=NULL,
                           frows=2,
                           fscales='fixed',
                           fdir='h',
                           w=6,h=5,
#                            angleHjVj=c(20,.5,.5),
                           debug=FALSE)
p1 <- gw_wres[[1]] +
        ggtitle(pheno_desc) +
        theme(plot.title = element_text(size=22, hjust=0, margin=margin(0,0,15,0)))
p1
ggsave(
    file.path(
      outdir,
      paste0(
        gw_fname,
        '.',
        fst_prefix,
        'wilcoxFST.png'
      )
    ),
    width=6,
    height=5
)
# don't save as svg since there will be many outlier points plotted
# ggsave(
#     file.path(
#       outdir,
#       paste0(
#         gw_fname,
#         '.',
#         fst_prefix,
#         'wilcoxFST.svg'
#       )
#     ),
#     width=6,
#     height=5
# )


# hypergeometric test results
hgres <- hglist[[1]]
# GWAS hits in diff- vs non-diff CREs for Wilcoxon test
gbqhits <- hglist[[2]]

# check that is at least one diff and one non-diff hit for ATAC or ChIP
# and not all FST values of either variant type are identical
# to dermine whether to run test
checktest <- gbqhits %>%
    filter(!is.na(gwas_FST),
           !is.na(gwas_pval)) %>%
    mutate(variant_type=case_when(p >= 0.05 ~ 'nondiff',
                                  p < 0.05 ~ 'diff')) %>%
    group_by(score_type, variant_type) %>%
    summarize(N = dplyr::n(),
              Nunique = n_distinct(gwas_FST),
              testableFST = if_else(N > 1 & Nunique==1, FALSE, TRUE)) %>%
    group_by(score_type) %>%
    mutate(vartypeN = dplyr::n())

if (any(checktest$vartypeN==2) & all(checktest$testableFST)) {

    # wilcoxon test for GWAS hit FST enrichment in diff-CREs relative to hits in non-diff-CREs
    wres <- tryCatch(
        {
            wres <- groupedBoxplot(gbqhits %>%
                               filter(!is.na(gwas_FST),
                                      !is.na(gwas_pval)) %>%
                               mutate(variant_type=case_when(p >= 0.05 ~ 'non-diff-CRE\nGWAS hit',
                                                             p < 0.05 ~ 'diff-CRE\nGWAS hit'),
                                      variant_type=fct_relevel(variant_type, rev),
                                      score_type=case_when(score_type=='chip' ~ 'ChIP',
                                                           score_type=='atac' ~ 'ATAC',
                                                           TRUE ~ score_type)),
                           'variant_type',
                           'gwas_FST',
                           xlab='Variant type',
                           ylab='FST',
                           groupvar=NULL,
                           plotpoints=TRUE,
                           points_thresh=200, # show all points for xvar with N observations below this number
                           show_wilcox=TRUE,
                           paired_t_test=FALSE,
                           testout=TRUE,
                           alt='greater',
                           orderpval=FALSE,
                           logscale=FALSE,
                           colorX=TRUE,
                           colorvals=c('darkgrey','red3'),
                           showleg=FALSE,
                           legendpos='right',
                           horizontal=TRUE,
                           facets='score_type',
                           frows=1,
                           fscales='fixed',
                           fdir='h',
                           w=10,h=5,
                           angleHjVj=c(0,.5,.5),
                           debug=FALSE)
        },
        error=function(e) {
            print(e)
            return(NA)
        }
    )

    if (!is.na(wres)) {
        #       wres %>%
    #           select(xvarcol, yvarcol)

        p2 <- wres[[1]] +
               ggtitle(pheno_desc) +
               theme(plot.title = element_text(size=22, hjust=0, margin=margin(0,0,15,0)))
        p2

        ggsave(
            file.path(
              outdir,
              paste0(
                gw_fname,
                '.',
                fst_prefix,
                'wilcoxFST.',
                bqtloutstr,
                '_CREs.png'
              )
            ),
            width=10,
            height=5
        )
        ggsave(
            file.path(
              outdir,
              paste0(
                gw_fname,
                '.',
                fst_prefix,
                'wilcoxFST.',
                bqtloutstr,
                '_CREs.svg'
              )
            ),
            width=10,
            height=5
        )

    } else {
        print("Not enough GWAS hits in LD with CRE bQTL to test or ties!")
    }

} else {
    print("Not enough GWAS hits in LD with CRE bQTL to test!")

    wres <- NA
}


# concat test results tables and write out

if (!is.na(wres)) {
  awres <- bind_rows(wres[[2]], gw_wres[[2]])
} else {
  awres <- gw_wres[[2]]
}

awres <- awres %>%
    mutate(gwas_significance = pval_gwas,
           phenotype_code = pheno_code,
           phenotype_description = pheno_desc)

testres <- bind_rows(hgres, awres)

testres %>%
    write_tsv(
      file.path(
        outdir,
        paste0(
          gw_fname,
          '.',
          bqtloutstr,
          '_CREs.testResults.txt'
        )
      )
    )
