#!/usr/bin/R

library(argparse)
library(ggplot2)
library(docstring)
library(readxl)
library(tidyverse)

# from diffinteractR.R

add_rsIDs <- function(df, genome='hg19', format=FALSE, debug=FALSE) {
    # adds rsID column to granges-formatted df
    if (format) {
      df <- df %>%
        separate(variant_id,
                 c('seqnames','start','ref','alt','build'),
                 sep="_",
                 remove=FALSE,
                 convert=TRUE) %>%
        mutate(end=start)
    }

    h2rs <- df %>%
        mutate(seqnames=gsub('chr','',seqnames)) %>%
        plyranges::as_granges()

    if (genome %in% c('hg19', 'GRCh37')) {
        snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37
    } else if (genome %in% c('hg38', 'GRCh38')) {
        snps <- SNPlocs.Hsapiens.dbSNP151.GRCh38::SNPlocs.Hsapiens.dbSNP151.GRCh38
    } else {
        stop('Genome not recognized!')
    }

    # print('Getting rsIDs with SNPloc...')
    rsids <- as.data.frame(BSgenome::snpsByOverlaps(snps, h2rs)) %>%
        as_tibble() %>%
        dplyr::rename(start=pos) %>%
        mutate(seqnames=paste0('chr',seqnames),
              end=start) %>%
        dplyr::rename(rsID=RefSNP_id) %>%
        dplyr::select(-c(strand,alleles_as_ambig))

    # print('Removing SNPs with no rsID...')
    dfrs <- as_tibble(merge(df, rsids))

    print(paste0(nrow(df)-nrow(dfrs)," variants dropped due to no rsID."))

    return(dfrs)
}

# from trans_utils.R

ld_expand_df <- function(r, plink_bed_dir, plink_pre,
                         tmpdir='./tmp',
                         chrcol='seqnames',
                         poscol='start',
                         rscol='rsID',
                         plink_dir="/home/trevor/Programs",
                         sep="") {

  cursnp = r[rscol]

  curchr = as.character(r[chrcol])
  curpos = as.integer(r[poscol])

  print(paste0(curchr, '_', curpos))

  # print(
  #   paste(curchr,curpos,sep='_')
  # )

  curchrname = paste(plink_pre,curchr,sep=sep)
  # print(curchrname)
  syscomm <- paste0(
      "mkdir -p ",
      tmpdir,
      "; ",
      plink_dir,
      "/plink --bfile ",
      plink_bed_dir,
      "/",
      curchrname,
      " --r2 --ld-snp ",
      cursnp,
      " --ld-window 9999 --ld-window-kb 500 --ld-window-r2 .8 --out ",
      tmpdir,
      "/snpld",
      cursnp
  )
  # print(syscomm)

  cursnpcall = system(syscomm,ignore.stdout=TRUE,ignore.stderr=TRUE)

  curreadfile = paste0(tmpdir,"/snpld",cursnp,".ld")

  if(file.exists(curreadfile)) {

    cursnpdatac <- tryCatch(
      {
        cursnpdata <- read.table(curreadfile,header=TRUE) %>%
            as_tibble() %>%
            mutate(coord=paste0('chr',CHR_B,'_',BP_B))

        cursnpdatac = as.character(cursnpdata$coord)
      },
      error=function(e) {
        message(e)
        return("novalidvar")
      }
    )

  } else {
      print('No LD file found')
  }

  if(!file.exists(curreadfile)) {

    cursnpdatac = "novalidvar"

  }

  #deal with rsIDs separated by semicolon that point to same coordinates
  snpIDs = unlist(strsplit(cursnpdatac,";"))

  if(is.null(snpIDs)) {
    snpIDs = cursnp
  }

  if("novalidvar" %in% snpIDs) {
    snpIDs = cursnp
  }

  return(as.list(snpIDs))

}

get_ld_snps <- function(df,
                        plink_bed_dir="/cashew/users/kade/nextflow/selection_1000G/vcf/output170620/plink",
                        plink_pre='CEU',
                        prefix_sep='.',
                        plink_dir="/home/kpettie/bin/plink_linux_x86_64_20200616",
                        rscol='rsID',
                        chrcol='seqnames',
                        poscol='start',
                        mergecols=c('module','seqnames','start','rsID','variant_id'),
                        replace_ranges=TRUE) {

    rs_df <- df %>%
        dplyr::select(any_of(mergecols)) %>%
        distinct()

    startcol <- 'start'
    if (!('start' %in% colnames(rs_df))) {
        startcol <- grep('pos', colnames(rs_df), value=TRUE)
        rs_df <- rs_df %>%
            dplyr::rename(start := {{startcol}}) %>%
            mutate(end = start)
    }

    if (!('rsID' %in% colnames(df))) {
        print("Adding rsIDs...")
        rs_df <- add_rsIDs(rs_df, format=replace_ranges)
    }


    rs_df$ld_snps <- apply(rs_df, MARGIN=1, FUN=ld_expand_df,
                            plink_bed_dir=plink_bed_dir,
                            plink_pre=plink_pre,
                            plink_dir=plink_dir,
                            rscol=rscol,
                            chrcol=chrcol,
                            poscol=poscol,
                            sep=prefix_sep)

    df <- merge(df,
                rs_df %>%
                    dplyr::select(-end) %>%
                    dplyr::rename(!!startcol := start),
                by=mergecols)

    if (replace_ranges) {
        df <- unnest(df,ld_snps) %>%
            dplyr::select(-c(seqnames,start,end)) %>%
            separate(., ld_snps, c("seqnames","start"), sep = "_", convert=TRUE, remove=FALSE) %>%
            mutate(end=start) %>%
            dplyr::rename(tag_variant_id=variant_id)
    } else {
        df <- unnest(df,ld_snps) %>%
            separate(., ld_snps, c(NA,"ld_pos"), sep = "_", convert=TRUE, remove=FALSE) %>%
            dplyr::rename(tag_rsID=rsID)
    }

    return(df)
}

hypergeomDiffQTLgwas <- function(df,
                                 testtrait,
                                 abc_p_col='chip_p',
                                 abc_p_thresh=0.05,
                                 qtl_p_thresh=0.05,
                                 df_sub=NULL,
                                 mat_out=FALSE,
                                 alt='greater') {

    df <- df %>%
      dplyr::rename(abc_p := {{abc_p_col}})

    if (!is.null(df_sub)) {
        if (df_sub=='cre-level') {
            df <- df %>%
                group_by(seqnames, start, end, name, Trait) %>%
                arrange(bqtl_p, desc(mean_ABC)) %>%
                dplyr::slice(1) %>%
                ungroup()
        } else if (df_sub=='gwas-level') {
            df <- df %>%
                mutate(gwas_pos_proxy = ifelse(is.na(gwas_pos), name, gwas_pos)) %>%
                group_by(seqnames, gwas_pos_proxy, Trait) %>%
                arrange(abc_p, bqtl_p, desc(mean_ABC)) %>%
                dplyr::slice(1) %>%
                ungroup()
        } else {
            stop('DF sub not recognized.')
        }
    }

    dfs <- df %>%
        group_by(seqnames,start,end,name) %>%
        filter(all(is.na(Trait)) | Trait==testtrait) %>%
        ungroup()

    if (grepl('chip', abc_p_col) | grepl('hic', abc_p_col)) {
      # only count each HiChIP bin once
      dfs <- dfs %>%
            group_by(hgnc_symbol, chip_lfc, abc_p) %>%
            mutate(chipbin_has_promoter = if_else(any(class=='promoter'), TRUE, FALSE))

      # take the CRE linked to a trait, then with the most significant QTL in the HiChIP bin
      dfs <- dfs %>%
            arrange(Trait, bqtl_p, .by_group=TRUE) %>%
            slice(1) %>%
            ungroup()

    }

    fg_numerator <- dfs %>%
        filter(Trait==testtrait,
               abc_p < abc_p_thresh) %>%
        nrow()
    fg_denominator <- dfs %>%
        filter(abc_p < abc_p_thresh) %>%
        nrow()
    bg_numerator <- dfs %>%
        filter(Trait==testtrait) %>%
        nrow()
    bg_denominator <- dfs %>%
        nrow()

    mat11 <- fg_numerator
    mat21 <- fg_denominator-fg_numerator
    mat12 <- bg_numerator-fg_numerator
    mat22 <- bg_denominator-mat11-mat21-mat12
    tmat <- matrix(c(mat11, mat21, mat12, mat22), nrow = 2,
                      dimnames =
               list(c('GWAS-linked', 'no GWAS link'),
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

aggHyperGwas <- function(dfa, sts=c('chip','atac'), df_sub=NULL) {

    alltraits <- dfa %>%
        filter(!is.na(Trait)) %>%
        pull(Trait) %>%
        sort %>%
        unique

    hres <- tibble()

    for (st in sts) {

        abcpcol <- paste0(st, '_p')
        df <- dfa %>%
          dplyr::rename(abc_p := {{abcpcol}})

        if (!is.null(df_sub)) {
            if (df_sub=='cre-level') {
                df <- df %>%
                    group_by(seqnames, start, end, name, Trait) %>%
                    arrange(bqtl_p, desc(mean_ABC)) %>%
                    dplyr::slice(1) %>%
                    ungroup()
            } else if (df_sub=='gwas-level') {
                df <- df %>%
                    mutate(gwas_pos_proxy = ifelse(is.na(gwas_pos), name, gwas_pos)) %>%
                    group_by(seqnames, gwas_pos_proxy, Trait) %>%
                    arrange(abc_p, bqtl_p, desc(mean_ABC)) %>%
                    dplyr::slice(1) %>%
                    ungroup()
            } else {
                stop('DF sub not recognized.')
            }
        }

        df <- df %>%
            dplyr::rename(!!abcpcol := abc_p)

        print(paste0('Running ', st, ' tests...'))

        for (currtrait in alltraits) {

            print(paste0('    Testing ', currtrait, ' ...'))

            cres <- hypergeomDiffQTLgwas(df,
                                         currtrait,
                                         abc_p_col=paste0(st, '_p'),
                                         abc_p_thresh=0.05,
                                         mat_out=FALSE,
                                         df_sub=NULL) %>%
                as_tibble() %>%
                mutate(score_type=st,
                       trait=currtrait)

            hres <- rbind(hres, cres)

        }
    }

    return(hres)

}


parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--gw2_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--d_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--qtl_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--qtl_name',
                    type='character',
                    default=NULL)
parser$add_argument('--frq_fname',
                    type='character',
                    default=NULL)
parser$add_argument('--fst_fname',
                    type='character',
                    default=NULL)
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

gw2_fname <- opt$gw2_fname
d_fname <- opt$d_fname
qtl_fname <- opt$qtl_fname
qtl_name <- opt$qtl_name
frq_fname <- opt$frq_fname
fst_fname <- opt$fst_fname

source(file.path(opt$plotdir,'plotting.R'))


gw2 <- read_excel(gw2_fname, sheet=4, skip=2) %>%
    rename_with(., ~ gsub('/', '', gsub('-', '', gsub(' ', '', .x)))) %>%
    rename(gwas_rsID=rsID,
           seqnames=Chr,
           gwas_pos=Position,
           gwas_p=Pvalue) %>%
    mutate(seqnames=paste0('chr',seqnames))

otherbqtlcols <- c('Depth','ALTdepth','REFDepth','ALTallele','POSTallele','POSTfreq','prechipfreq')

# CREs with bQTL already overlapped
d <- read_tsv(d_fname) %>%
    select(-starts_with('bqtl_'),
           -any_of(otherbqtlcols)) %>%
    group_by(hgnc_symbol, name) %>%
    # limit top eQTL per peak from previous overlap
    arrange(eqtl_lcl_p, eqtl_pbmc_p) %>%
    dplyr::slice(1) %>%
    ungroup()


# bQTL

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


# merge CREs with bQTL
# widen and add bQTL
wabc <- plyranges::join_overlap_left(
    d %>% plyranges::as_granges(),
    qtl %>%
        dplyr::rename(seqnames=Chr,
                      start=position,
                      bqtl_p=pvalue,
                      bqtl_alt_high_freq=alt_high_freq) %>%
        mutate(end=start+1,
               bqtl_pos=start) %>%
        plyranges::as_granges()
    ) %>%
    as_tibble() %>%
    mutate(bqtl_direction = case_when(
                ((bqtl_alt_affinity=='high' & bqtl_alt_high_freq=='AFR') |
                    (bqtl_alt_affinity=='low' & bqtl_alt_high_freq=='EUR')) ~ 'down',
                ((bqtl_alt_affinity=='high' & bqtl_alt_high_freq=='EUR') |
                    (bqtl_alt_affinity=='low' & bqtl_alt_high_freq=='AFR')) ~ 'up'
            )
          )

# subset to only significant bQTL for gwas overlap
wabc_sbqtl <- wabc %>%
    filter(bqtl_p < 0.05)

ldbqtl <- get_ld_snps(wabc_sbqtl,
                        plink_bed_dir="/cashew/users/kade/nextflow/selection_1000G/vcf/output170620/plink",
                        plink_pre='CEU',
                        prefix_sep='.',
                        plink_dir="/home/kpettie/bin/plink_linux_x86_64_20200616",
                        rscol='rsID',
                        chrcol='seqnames',
                        poscol='start',
                        mergecols=c('seqnames','bqtl_pos'),
                        replace_ranges=FALSE) %>%
          mutate(bqtl_name=qtl_name)

ldbqtl %>%
    write_tsv(
      file.path(opt$outdir, paste0(opt$name, '.expandLD.txt.gz'))
    )

print('Debug here:')
ldbqtl <- ldbqtl %>%
    select(-ld_snps) %>%
    rename(gwas_pos = ld_pos) %>%
    mutate(gwas_pos = if_else(is.na(gwas_pos), as.double(bqtl_pos), as.double(gwas_pos)))

# above should work so overwrite here
ldbqtl %>%
    write_tsv(
      file.path(opt$outdir, paste0(opt$name, '.expandLD.txt.gz'))
    )

print("Merging GWAS...")
bqg <- merge(ldbqtl, gw2, all.x=TRUE) %>%
    # remove LD expansions of bQTL in LD with GWAS that did not directly overlap GWAS
    group_by(hgnc_symbol,seqnames,start,end,name,bqtl_pos,tag_rsID) %>%
    mutate(gwas_linked = if_else(any(!is.na(Trait)), TRUE, FALSE)) %>%
    ungroup() %>%
    filter(!(gwas_linked & is.na(Trait))) %>%
    group_by(hgnc_symbol,seqnames,start,end,name,bqtl_pos,tag_rsID,Trait) %>%
    mutate(pmids = case_when(any(!is.na(PUBMEDID)) ~ paste0(PUBMEDID[!is.na(PUBMEDID)], collapse=';')),
           sources = case_when(any(!is.na(Source)) ~ paste0(Source[!is.na(Source)], collapse=';'))) %>%
    arrange(gwas_p, desc(RiskAllele), .by_group=TRUE) %>%
    dplyr::slice(1) %>%
    ungroup()



bqg %>%
    write_tsv(
      file.path(opt$outdir, paste0(opt$name, '.expandLDgwas.txt.gz'))
    )

gwasenrich <- aggHyperGwas(bqg, df_sub='gwas-level') %>%
    arrange(pvalue)

gwasenrich %>%
    write_tsv(
      file.path(opt$outdir, paste0(opt$name, '.expandLDgwas.hypergeom.txt'))
    )

# set font for plotting to Arial
# so illustrator displays .eps/.svg correctly
# and Arial mandated by some journals (e.g., PLOS)
extrafont::font_import(path=opt$fontdir, prompt=FALSE)
extrafont::choose_font('Arial')

p1 <- plotFisherEnrichments(gwasenrich %>%
                          group_by(score_type) %>%
                          arrange(pvalue, .by_group=TRUE) %>%
                          dplyr::slice(1:10) %>%
                          ungroup() %>%
                          mutate(trait = if_else(trait=='Pulmonary artery enlargement in chronic obstructive pulmonary disease',
                                                 'Pulmonary artery enlargement in\nchronic obstructive pulmonary disease',
                                                 trait)),
                      yvar_col='trait',
                      ylabel="GWAS trait",
                      xvar_col="odds_ratio",
                      althypoth='greater',
                      success_total=TRUE,
                      line_intercept=1,
                      groupvar=NULL,
                        groupvarord=NULL,
                        colorvals=c('black','purple2'),
                        sep_groups=TRUE,
                        w=19,
                        h=7,
                        facets='score_type',
                        frows=1,
                        fscales='free', # free, free_x, free_y
                        fdir='v',
                        legendpos=c(.3,.2), # c(xpos, ypos)
                        debug=FALSE) +
                  ggtitle(paste0(qtl_name, "-linked GWAS enrichment in diff-CREs (GWAS-level)")) +
                  theme(plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,15,0)))


p1
ggsave(file.path(opt$outdir, paste0(opt$name, '.expandLDgwas.hypergeomTop10.png')), width=19, height=7)
ggsave(file.path(opt$outdir, paste0(opt$name, '.expandLDgwas.hypergeomTop10.eps')), width=19, height=7)
