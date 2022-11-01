#!/usr/bin/R

library(argparse)
library(ggplot2)
library(tidyverse)

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

        if (!('gwas_RiskAllele' %in% colnames(dfa))) {
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
        } else {
            df <- dfa %>%
                filter(seqnames==paste0('chr', currchr)) %>%
                # remove gwas freq cols, but not bqtl freq cols
                select(-contains(c('frq','freq')), contains('bqtl_'))
        }

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

        if (!('gwas_FST' %in% colnames(dfm))) {

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

        }

        adfm <- rbind(adfm, dfm)

    }

    return(adfm)

}

concatForPlotting <- function(fdir,
                             f_patt='fisherEnrichments.txt$',
                             checkfirst=FALSE) {

    fnames <- list.files(fdir,
                         pattern=f_patt,
                         full.names=FALSE)

    adf <- tibble()
    for (f in fnames) {

        print(paste0('Reading ',f))

        df <- read_tsv(file.path(fdir, f))

        adf <- rbind(adf, df)

        if (checkfirst) return(adf)

    }

    return(adf)

}

concatForRisk <- function(fdir,
                          phenotype_prevalence,
                          f_patt='_CREs.txt.gz',
                          abc_p_thresh=0.05,
                         abc_p_nonthresh=0.05,
                         gwas_thresh=5e-5,
                         gwas_nonthresh=5e-5,
                          checkfirst=FALSE) {

    fnames <- list.files(fdir,
                         pattern=f_patt,
                         full.names=FALSE)

    abc_cols <- c('seqnames',
                  'start',
                  'end',
                  'name',
                  'hgnc_symbol',
                  'TargetGene',
                  'chip_lfc',
                  'chip_p',
                  'atac_lfc',
                  'atac_p',
                  'min_p')

    pheno_codes <- phenotype_prevalence %>% pull(phenotype_code) %>% unique

    adf <- tibble()
    for (f in fnames) {

        pcode <- str_split(f, '\\.', simplify=TRUE)[,1]

        if (!(pcode %in% pheno_codes)) next

        print(paste0('Reading ',f))

        df <- read_tsv(file.path(fdir, f)) %>%
            filter(gwas_pval < gwas_thresh)

        cpp <- phenotype_prevalence %>%
            filter(phenotype_code==pcode) %>%
            select(contains('pheno'), anc_assoc)

        df <- bind_cols(df, cpp) %>%
            select(any_of(abc_cols), contains(c('gwas','pheno','anc_assoc','bqtl')))

        # use bind_rows instead of rbind, as some dfs will have
        # more columns than others (e.g., gwas_expected_case_minor_AC)
        adf <- bind_rows(adf, df)

        if (checkfirst) return(adf)

    }

    return(adf)

}

addSubstringPhenotypes <- function(df, allphenos, debug=FALSE) {

    allphenos <- allphenos %>%
        filter(test_type=='gwasDiffCRE') %>%
        group_by(phenotype_code) %>%
        arrange(pvalue, .by_group=TRUE) %>%
        slice(1) %>%
        ungroup()

    pcodes <- df %>% pull(phenotype_code)

    adf <- tibble()

    for (pcode in pcodes) {

        cdf <- df %>% filter(phenotype_code==pcode)
        c_anc_assoc <- cdf[['anc_assoc']]
        phenomatch <- cdf[['pheno_substr']]
        phenononmatch <- cdf[['nonpheno_substr']]

        cphenos <- allphenos %>%
            filter(grepl(phenomatch, phenotype_description, ignore.case=TRUE),
                   !grepl(phenononmatch, phenotype_description, ignore.case=TRUE),
                   phenotype_code != pcode) %>%
            mutate(anc_assoc = c_anc_assoc)

        adf <- bind_rows(adf, cdf, cphenos)

    }

    return(adf)

}

riskPrevalenceBinomial <- function(df,
                                   fst_thresh=NULL,
                                   test_type='binom',
                                   dfout=FALSE) {

    if (!is.null(fst_thresh)) df <- df %>% filter(gwas_FST > fst_thresh)

    bgmatch <- df %>%
        filter(anc_assoc==gwas_risk_high_freq) %>%
        nrow()

    psuccess <- bgmatch / nrow(df)

    ntrials <- df %>%
        filter(min_p < 0.05) %>%
        nrow()

    nsuccess <- df %>%
        filter(min_p < 0.05,
               anc_assoc==gwas_risk_high_freq) %>%
        nrow()

    if (test_type=='binom') {
        res <- binom.test(nsuccess, ntrials, p=psuccess, alternative='greater')
    } else if (test_type=='hyper') {
        res <- phyper(nsuccess-1, bgmatch, nrow(df)-bgmatch, ntrials, lower.tail=FALSE)
    } else if (test_type=='fisher') {
        tmat <- matrix(
            c(nsuccess,
              bgmatch-nsuccess,
              ntrials-nsuccess,
              nrow(df)-bgmatch-ntrials+nsuccess),
            2,
            2
        )
        res <- fisher.test(tmat, alternative='greater')
    } else {
        stop('Test type not recognized')
    }

    if (dfout) {

        Nnondiff <- df %>%
            filter(min_p >= 0.05) %>%
            nrow()
        Nnondiffmatch <- df %>%
            filter(min_p >= 0.05,
                   anc_assoc==gwas_risk_high_freq) %>%
            nrow()
        Ndiff <- ntrials
        Ndiffmatch <- nsuccess

        propnondiffsuccess <- tribble(
            ~cre_type, ~prop_risk_matching,
            'non-diff', Nnondiffmatch/Nnondiff,
            'diff',    Ndiffmatch/Ndiff
        )

        return(list(res, propnondiffsuccess))

    }

    return(res)

}

# df with manually curated ancestry-disease phenotype associations
# with full phenotype name and substring to detect any nominally significant
# related phenotype to include variants in risk allele direction test
# some pheno_substr will match phenotypes that only include variants for background
# risk allele matching rate
# nonpheno_substr: for phenotypes including pheno_substr, exclude those that also match this substring
prevdf <- tribble(
    ~ phenotype_description,                       ~ pheno_substr,          ~nonpheno_substr,  ~ anc_assoc,
    'Peripheral angiopathy',                       'angiopathy',            NA,                'AFR',
    'E11 Non-insulin-dependent diabetes mellitus', 'Non-insulin-dependent', NA,                'AFR',
    'Type 2 diabetes with peripheral circulatory complications', 'Type 2 diabetes with peri', NA, 'AFR',
    'Type 2 diabetes',                             'Type 2 diabetes',       'with',            'AFR',
    'I78 Diseases of capillaries',                 'capillar',              NA,                'AFR',
    'heart attack/myocardial infarction',          'heart attack/myo',      NA,                'AFR',
    # unclear if any association for basophil percentage
    'Basophill count',                             'Basophill count',       NA,                'EUR',
    'essential hypertension',                      'hypertension',          'pregnancy|[0-9]', 'AFR',
    'transient ischaemic attack (tia)',            'transient ischaemic',   NA,                'AFR',
    'E04 Other non-toxic goitre',                  'goitre',                NA,                'AFR',
    'K92 Other diseases of digestive system',      'Other diseases of(.*)?digestive system|helicobacter', NA, 'AFR',
    'Lymphocyte count',                            'Lymphocyte count',      NA,                'AFR',
    'sarcoidosis',                                 'sarcoidosis',           NA,                'AFR',
    'I51 Complications and ill-defined descriptions of heart disease', 'descriptions of heart disease', NA, 'AFR',
    'Platelet distribution width',                 'Platelet distribution', NA,                'AFR',
    'low platelets/platelet disorder',             'low platelets',         NA,                'AFR',
    'Diseases of veins, lymphatic vessels and lymph nodes, not elsewhere classified', 'Diseases of veins', NA, 'AFR',
    'abdominal hernia',                            'abdominal hernia',      NA,                'EUR',
    'connective tissue disorder',                  'connective tissue disorder', NA,           'AFR',
    'systemic lupus erythematosis/sle',            'lupus',                 NA,                'AFR',
    'malabsorption/coeliac disease',               'malabsorption|coeliac', NA,                'EUR',
    'Psoriasis',                                   'psoriasis',             NA,                'EUR',
    'Rosacea',                                     'Rosacea',               NA,                'EUR',
    'L71 Rosacea',                                 'L71 Rosacea',           NA,                'EUR',
    'Mean corpuscular haemoglobin concentration',  'Mean corpuscular haemoglobin concentration', NA, 'EUR',
    'ankylosing spondylitis',                      'ankylosing spondylitis', 'M45',            'EUR',
    'herpes simplex',                              'herpes simplex',        NA,                'AFR',
    'thyroxine product',                           'thyroxine product|hyperthyroidism', NA,    'AFR',
    'infective/viral hepatitis',                   'infective/viral hepatitis', NA,            'AFR',
    'Blood clot in the lung',                      'blood clot in the lung|pulmonary embolism', NA, 'AFR',
    'crohns disease',                              'crohn',                 'not crohn',       'EUR',
    'Aortic aneurysm',                             'Aortic aneurysm',       NA,                'EUR'
)

# manually curated lower confidence ancestry-phenotype associations for hits
# with multiple trait associations and conflicting directionality
low_conf <- c('chr1;156571217;crohns disease',
              "chr1;156571217;Crohn's disease of large intestine",
              "chr1;156571217;K50 Crohn's disease [regional enteritis]",
              'chr1;172418268;Platelet distribution width',
              "chr10;115719472;Lymphocyte count",
              'chr11;6624660;abdominal hernia',
              'chr11;9780891;Platelet distribution width',
              'chr11;9780891;hypertension',
              'chr11;118747576;Lymphocyte count',
              'chr12;740009;Platelet distribution width',
              'chr12;9911584;Lymphocyte count',
              'chr12;9911586;Lymphocyte count',
              'chr14;57857162;Lymphocyte count',
              'chr15;41789091;Platelet distribution width',
              'chr15;66076794;Platelet distribution width',
              'chr16;29856082;Rosacea',
              'chr16;29856082;L71 Rosacea',
              'chr17;1675174;I51 Complications and ill-defined descriptions of heart disease',
              'chr19;13076186;Lymphocyte count',
              'chr19;45430280;Platelet distribution width',
              'chr3;49208865;Platelet distribution width',
              'chr6;26090025;Platelet distribution width',
              'chr6;26090025;Diseases of veins, lymphatic vessels and lymph nodes, not elsewhere classified',
              'chr6;26190874;hypertension',
              'chr6;26286839;Diseases of veins, lymphatic vessels and lymph nodes, not elsewhere classified',
              'chr6;26286839;hypertension',
              'chr6;26399615;Diseases of veins, lymphatic vessels and lymph nodes, not elsewhere classified',
              'chr6;27491299;Diseases of veins, lymphatic vessels and lymph nodes, not elsewhere classified',
              'chr6;27859453;Diseases of veins, lymphatic vessels and lymph nodes, not elsewhere classified',
              'chr6;28400295;sarcoidosis')

parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--r_pattern',
                    type='character',
                    default='_CREs.txt.gz')
parser$add_argument('--t_pattern',
                    type='character',
                    default='_CREs.testResults.txt')
parser$add_argument('--t_dir',
                    type='character',
                    default='.',
                    help="Directory with results from combinedBqtlGwas.R")
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

t_dir <- opt$t_dir
t_pattern <- opt$t_pattern
r_pattern <- opt$r_pattern
outdir <- opt$outdir

# fst threshold for additional risk allele matching test
topfst <- 0.189

source(file.path(opt$plotdir,'plotting.R'))
source(file.path(this.path::this.dir(default='~/code/github/selection_1000G'), 'variantTools.R'))

# set font for plotting to Arial
# so illustrator displays .eps/.svg correctly
# and Arial mandated by some journals (e.g., PLOS)
extrafont::font_import(path=opt$fontdir, prompt=FALSE)
extrafont::choose_font('Arial')


# concat hypergeometric test restults for plotting top 10 from ChIP and ATAC
tres <- concatForPlotting(t_dir, t_pattern, checkfirst=FALSE)
tres %>%
    write_tsv(
      file.path(
        outdir,
        'aggTestResults.txt.gz'
      )
    )

p1 <- plotFisherEnrichments(tres %>%
                                filter(test_type=='gwasDiffCRE') %>%
                          group_by(score_type) %>%
                          arrange(pvalue, .by_group=TRUE) %>%
                          dplyr::slice(1:10) %>%
                          ungroup() %>%
                          mutate(phenotype_description = case_when(phenotype_description=='Diseases of veins, lymphatic vessels and lymph nodes, not elsewhere classified'
                                                      ~ 'Diseases of veins, lymphatic vessels and\nlymph nodes, not elsewhere classified',
                                                   TRUE ~ phenotype_description)),
                      yvar_col='phenotype_description',
                      ylabel="GWAS trait",
                      xvar_col="odds_ratio",
                      althypoth='greater',
                      success_total=TRUE,
                      line_intercept=1,
                      groupvar=NULL,
                        groupvarord=NULL,
                        colorvals=c('black','purple2'),
                        sep_groups=TRUE,
                        w=27,
                        h=9,
                        facets='score_type',
                        frows=1,
                        fscales='free', # free, free_x, free_y
                        fdir='v',
                        legendpos=c(.3,.2), # c(xpos, ypos)
                        debug=FALSE) +
                  ggtitle(paste0("Top diff-CRE-enriched, JUND/NFKB/PU.1-linked GWAS trais have ancestry-associated prevalence")) +
                  theme(plot.title = element_text(size=30, hjust=.5, margin=margin(0,0,15,0)))
p1
ggsave(
  file.path(
    outdir,
    'hypergeomTop10.png'
  ),
  width=27,
  height=9
)


# add putative ancestry-prevalence associations for risk allele matching
pres <- merge(prevdf,
              tres %>%
                 filter(test_type=='gwasDiffCRE') %>%
                 group_by(phenotype_code) %>%
                 arrange(pvalue, .by_group=TRUE) %>%
                 slice(1) %>%
                 ungroup())
# add redundant phenotypes that did not come up directly in hypergeom nominal
# enrichments in case they add additional variants for risk allele matching
pres <- addSubstringPhenotypes(pres, tres)
# get variant-level info for only these traits for risk allele test
rdf <- concatForRisk(t_dir,
                      pres,
                      f_patt=r_pattern,
                      abc_p_thresh=0.05,
                     abc_p_nonthresh=0.05,
                     gwas_thresh=5e-5,
                     gwas_nonthresh=5e-5,
                      checkfirst=FALSE)

# should not need to merge fst or frq here as they should have already been
# merged in combinedBqtlGwas.R

rdf <- rdf %>%
    select(seqnames, gwas_pos, gwas_ref, gwas_alt, gwas_RiskAllele, gwas_frqref, gwas_frqalt, phenotype_description, anc_assoc, gwas_FST, everything()) %>%
    group_by(seqnames, gwas_pos, phenotype_code) %>%
    arrange(min_p) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(gwas_risk_refalt = if_else(gwas_RiskAllele==gwas_alt, 'alt', 'ref'),
           gwas_risk_high_freq = case_when(gwas_risk_refalt=='alt' ~ gwas_alt_high_freq,
                                           gwas_risk_refalt=='ref' & gwas_alt_high_freq=='EUR' ~ 'AFR',
                                           gwas_risk_refalt=='ref' & gwas_alt_high_freq=='AFR' ~ 'EUR'))

nonconflicts <- rdf %>%
   group_by(seqnames, gwas_pos) %>%
   filter(n_distinct(gwas_RiskAllele)==1,
          n_distinct(anc_assoc)==1) %>%
   arrange(gwas_pval, .by_group=TRUE) %>%
   slice(1) %>%
   ungroup()

conflicts <- rdf %>%
   mutate(hit_id = paste(seqnames,gwas_pos,phenotype_description,sep=';')) %>%
   group_by(seqnames, gwas_pos) %>%
   filter(n() > 1,
          n_distinct(gwas_RiskAllele) > 1 | n_distinct(anc_assoc) > 1) %>%
   filter(!(hit_id %in% low_conf)) %>%
   arrange(gwas_pval, .by_group=TRUE) %>%
   slice(1) %>%
   ungroup() %>%
   select(-hit_id)

trisk <- rbind(nonconflicts, conflicts)


sink(
  file.path(
    outdir,
    "riskPrevalenceTests.Rout"
  )
)

print('Binomial test all:')
riskPrevalenceBinomial(trisk)
print('\n')
print('Hypergeometric test all:')
riskPrevalenceBinomial(trisk, test_type='hyper')
prin('\n')
print('Binomial test top 5% FST:')
riskPrevalenceBinomial(trisk, fst_thresh=topfst)
prin('\n')
print('Hypergeometric test top 5% FST:')
riskPrevalenceBinomial(trisk, fst_thresh=topfst, test_type='hyper')
print('\n')
print("Fisher's exact test top 5% FST:")
riskPrevalenceBinomial(trisk, fst_thresh=topfst, test_type='fisher')
print('\n')

sink()


rrl <- riskPrevalenceBinomial(trisk, fst_thresh=NULL, test_type='hyper', dfout=TRUE)
rrlfst <- riskPrevalenceBinomial(trisk, fst_thresh=topfst, test_type='hyper', dfout=TRUE)

rpmdf <- rbind(
    rrl[[2]] %>% mutate(fst_perc_thresh='all FST'),
    rrlfst[[2]] %>% mutate(fst_perc_thresh='top 5% FST')
)


options(repr.plot.width = 7, repr.plot.height = 7, repr.plot.res = 200)

titlesize <- 28
textsize <- 26

p2 <- rpmdf %>%
    mutate(cre_type = fct_relevel(factor(cre_type),c('non-diff','diff')),
          fst_perc_thresh = fct_relevel(factor(fst_perc_thresh),c('all FST','top 5% FST'))) %>%
    ggplot(aes(x=cre_type,
                   y=prop_risk_matching,
                   fill=fst_perc_thresh)) +
        geom_bar(stat="identity",
                 width=0.5,
                 position=position_dodge(0.6),
                 show.legend=TRUE) +
    theme_classic() +
    scale_fill_manual(values=c('black','purple2')) +
    labs(y = "proportion of risk alleles matching\nancestry-trait association",
         x = "CRE type") +
    theme(axis.text.x = element_text(size=textsize),
          axis.text.y = element_text(size=textsize),
          axis.title.x = element_text(size=titlesize),
          axis.title.y = element_text(size=titlesize),
          legend.text = element_text(size=textsize),
          legend.title = element_blank(),
          legend.position = 'top',
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent",colour = NA))
p2
ggsave(
  file.path(
    outdir,
    'proportionMatchingBarplot.png'
  ),
  width=7,
  height=7
)


# wilcoxon test for GWAS hit FST enrichment in diff-CREs relative to hits in non-diff-CREs
wres <- groupedBoxplot(narrow_diffABC(trisk, displayformat=TRUE) %>%
                               filter(!is.na(gwas_FST),
                                      !is.na(gwas_pval)) %>%
                               mutate(variant_type=case_when(p >= 0.05 ~ 'non-diff-CRE\nGWAS hit',
                                                             p < 0.05 ~ 'diff-CRE\nGWAS hit'),
                                      variant_type=fct_relevel(variant_type, rev)),
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

wres[[1]]
ggsave(
  file.path(
    outdir,
    'diffNonDiffFstWilcoxAtacChip.png'
  ),
  width=10,
  height=5
)

# wilcoxon test for GWAS hit FST enrichment in diff-CREs relative to hits in non-diff-CREs
wres <- groupedBoxplot(narrow_diffABC(trisk, displayformat=TRUE) %>%
                               filter(!is.na(gwas_FST),
                                      !is.na(gwas_pval),
                                      anc_assoc==gwas_risk_high_freq) %>%
                               mutate(variant_type=case_when(p >= 0.05 ~ 'non-diff-CRE\nGWAS hit',
                                                             p < 0.05 ~ 'diff-CRE\nGWAS hit'),
                                      variant_type=fct_relevel(variant_type, rev)),
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

wres[[1]]
ggsave(
  file.path(
    outdir,
    'diffNonDiffFstWilcoxAtacChipMatching.png'
  ),
  width=10,
  height=5
)

wres <- groupedBoxplot(trisk %>%
                               filter(!is.na(gwas_FST),
                                      !is.na(gwas_pval)) %>%
                               mutate(variant_type=case_when(min_p >= 0.05 ~ 'non-diff-CRE\nGWAS hit',
                                                             min_p < 0.05 ~ 'diff-CRE\nGWAS hit'),
                                      variant_type=fct_relevel(variant_type, rev)),
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
                           facets=NULL,
                           frows=1,
                           fscales='fixed',
                           fdir='h',
                           w=7,h=5,
                           angleHjVj=c(0,.5,.5),
                           debug=FALSE)

wres[[1]]
ggsave(
  file.path(
    outdir,
    'diffNonDiffFstWilcox.png'
  ),
  width=7,
  height=5
)

wres <- groupedBoxplot(trisk %>%
                               filter(!is.na(gwas_FST),
                                      !is.na(gwas_pval),
                                      anc_assoc==gwas_risk_high_freq) %>%
                               mutate(variant_type=case_when(min_p >= 0.05 ~ 'non-diff-CRE\nGWAS hit',
                                                             min_p < 0.05 ~ 'diff-CRE\nGWAS hit'),
                                      variant_type=fct_relevel(variant_type, rev)),
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
                           facets=NULL,
                           frows=1,
                           fscales='fixed',
                           fdir='h',
                           w=7,h=5,
                           angleHjVj=c(0,.5,.5),
                           debug=FALSE)

wres[[1]]
ggsave(
  file.path(
    outdir,
    'diffNonDiffFstWilcoxMatching.png'
  ),
  width=7,
  height=5
)

# write cands for diffinteractR
rdf %>%
    filter(!is.na(gwas_FST),
          !is.na(gwas_pval),
          anc_assoc==gwas_risk_high_freq) %>%
    group_by(seqnames,gwas_pos) %>%
    mutate(Nphenos = n()) %>%
    ungroup() %>%
    arrange(desc(Nphenos), desc(gwas_FST)) %>%
    write_tsv(file.path(outdir, 'riskPrevalenceTopCands.txt'))

rdf %>%
    filter(!is.na(gwas_FST),
           !is.na(gwas_pval),
           gwas_FST > topfst) %>%
    mutate(variant_type=case_when(min_p >= 0.05 ~ 'non-diff',
                                  min_p < 0.05 ~ 'diff')) %>%
    group_by(phenotype_code,phenotype_description,anc_assoc,variant_type) %>%
    summarize(Nvars = n(),
              prop_matching = sum(gwas_risk_high_freq==anc_assoc, na.rm=TRUE)/Nvars) %>%
    group_by(phenotype_code,phenotype_description,anc_assoc) %>%
    filter(any(variant_type=='diff')) %>%
    mutate(prop_matching_diff = prop_matching[variant_type=='diff'],
           Nvars_diff = Nvars[variant_type=='diff']) %>%
    ungroup() %>%
    arrange(desc(prop_matching_diff),desc(Nvars_diff)) %>%
    write_tsv(file.path(outdir, 'diffNonDiffTopFstPropMatching.txt'))

rdf %>%
    filter(!is.na(gwas_FST),
           !is.na(gwas_pval),
           min_p < 0.05) %>%
    group_by(phenotype_code,phenotype_description,anc_assoc) %>%
    summarize(Nvars = n(),
              prop_matching = sum(gwas_risk_high_freq==anc_assoc,na.rm=TRUE)/Nvars,
              prop_fst_matching = sum(gwas_FST[gwas_risk_high_freq==anc_assoc])/sum(gwas_FST),
              mean_matching_fst = mean(gwas_FST[gwas_risk_high_freq==anc_assoc])) %>%
    ungroup() %>%
    arrange(desc(mean_matching_fst),desc(Nvars)) %>%
    write_tsv(file.path(outdir, 'diffMatchingMeanFst.txt'))

plfc <- rdf %>%
    filter(!is.na(gwas_FST),
           !is.na(gwas_pval),
           gwas_risk_high_freq==anc_assoc) %>%
    mutate(variant_type=case_when(min_p >= 0.05 ~ 'non-diff',
                                  min_p < 0.05 ~ 'diff')) %>%
    group_by(phenotype_code,phenotype_description,anc_assoc) %>%
    filter(length(gwas_FST[variant_type=='diff'])>0,
           length(gwas_FST[variant_type=='non-diff'])>0) %>%
    summarize(Nvars=n(),
               mean_diff = mean(gwas_FST[variant_type=='diff']),
               mean_nondiff = mean(gwas_FST[variant_type=='non-diff']),
               fst_lfc = log2(mean_diff/mean_nondiff)) %>%
    ungroup() %>%
    arrange(desc(fst_lfc)) %>%
    filter(Nvars > 3)
plfc %>%
    write_tsv(file.path(outdir, 'diffNonDiffMatchingFstLFC.txt'))

options(repr.plot.width = 19, repr.plot.height = 14, repr.plot.res = 200)

calloutpheno <- 'systemic lupus erythematosis/sle'
calloutcolor <- 'green3'
diffcolor <- '#CD0000'
nondiffcolor <- 'darkgrey'

plfcf <- plfc %>%
    mutate(phenotype_description = case_when(phenotype_description=='Diseases of veins, lymphatic vessels and lymph nodes, not elsewhere classified'
                                                      ~ 'Diseases of veins, lymphatic vessels and\nlymph nodes, not elsewhere classified',
                                             phenotype_description=='Type 2 diabetes with peripheral circulatory complications'
                                                     ~ 'Type 2 diabetes with\nperipheral circulatory complications',
                                                   TRUE ~ phenotype_description),
           phenotype_description = fct_reorder(factor(phenotype_description), fst_lfc),
           barcolor = case_when(phenotype_description==calloutpheno ~ calloutcolor,
                                fst_lfc > 0 ~ diffcolor,
                                fst_lfc < 0 ~ nondiffcolor),
           ymaxcol = fst_lfc,
           yoffset = if_else(fst_lfc>0, .05, -.05)) %>%
    arrange(fst_lfc)

p3 <- plfcf %>%
    ggplot(aes(x=phenotype_description,
               y=fst_lfc,
               color=phenotype_description,
               fill=phenotype_description)) +
            geom_bar(stat="identity",
                     width=0.5,
                     show.legend=FALSE) +
            geom_text(size=9,
                      aes(y=ymaxcol+yoffset,label=Nvars,hjust=if_else(fst_lfc>0,0,1)),
                      show.legend=FALSE) +
        scale_color_manual(values=plfcf$barcolor) +
        scale_fill_manual(values=plfcf$barcolor) +
        theme_classic() +
        labs(y = expression(log[2]*" fold change FST (diff/non-diff)"),
             x = "GWAS trait") +
        coord_flip() +
        theme(axis.text.x = element_text(size=textsize),
              axis.text.y = element_text(size=textsize),
              axis.title.x = element_text(size=titlesize),
              axis.title.y = element_text(size=titlesize),
              panel.background = element_rect(fill = "transparent",colour = NA))
p3
ggsave(
    file.path(
        outdir,
        'Fig5g.png'
    ),
    width=19,
    height=14
)

rdf %>%
    filter(!is.na(gwas_FST),
           !is.na(gwas_pval),
           gwas_risk_high_freq==anc_assoc) %>%
    mutate(variant_type=case_when(min_p >= 0.05 ~ 'non-diff',
                                  min_p < 0.05 ~ 'diff')) %>%
    group_by(phenotype_code,phenotype_description,anc_assoc) %>%
    filter(length(gwas_FST[variant_type=='diff'])>0,
           length(gwas_FST[variant_type=='non-diff'])>0) %>%
    rstatix::wilcox_test(gwas_FST ~ variant_type, alternative='greater') %>%
    arrange(p) %>%
    write_tsv(file.path(outdir, 'diffNonDiffMatchingWilcoxIndividualPhenos.txt'))
