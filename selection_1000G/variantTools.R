#!/usr/bin/R

if(!(require(ggplot2)   )) install.packages("ggplot2")
if(!(require(tidyverse) )) install.packages("tidyverse")

library(ggplot2)
library(tidyverse)

# traceback to written code
options(error=function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })


separateFast <- function(df,
                         col,
                         into,
                         sep="[^[:alnum:]]+",
                         remove=TRUE,
                         coltypes=NULL) {

    #' Faster dplyr::separate
    #'
    #' @description Wrapper around stringi::stri_split_regex
    #' that allows separating a column on a regular expression
    #' much faster than dplyr::separate, but without some
    #' of the latter's usually extra special case functionality
    #'
    #' @param col character. column name of `df` to split
    #' @param into character vector. names of the columns to split `col` into
    #' @param sep character. regex to split on
    #' @param remove logical. remove original column once split
    #' @param coltypes character vector. column types to set of new columns

    cmat <- stringi::stri_split_regex(df[[col]],
                                       pattern=sep,
                                       simplify=TRUE)
    colnames(cmat) <- into
    ctib <- type_convert(as_tibble(cmat), col_types=coltypes)

    atib <- cbind(df, ctib)

    if (remove) atib <- atib %>% select(-all_of(c(col)))

    return(atib)

}

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
      plink_dir,
      "/plink --bfile ",
      plink_bed_dir,
      "/",
      curchrname,
      " --r2 --ld-snp ",
      cursnp,
      " --ld-window 9999 --ld-window-kb 500 --ld-window-r2 .8 --out /tmp/snpld",
      cursnp
  )
  # print(syscomm)

  cursnpcall = system(syscomm,ignore.stdout=TRUE,ignore.stderr=TRUE)

  curreadfile = paste0("/tmp/snpld",cursnp,".ld")

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


# allele frequency within 2 population groups

alleleFreqMerge <- function(frqa_fname, frqe_fname) {

    frqcnames <- c('seqnames','start','end','ref_frq','alt_frq')

    frqa <- read_tsv(frqa_fname, col_names=frqcnames)
    frqe <- read_tsv(frqe_fname, col_names=frqcnames)
    frq <- merge(
        frqa %>%
            separate(ref_frq,
                 into=c('ref','ref_frq'),
                 sep=":(?=[0-9])",
                 convert=TRUE,
                 extra='merge') %>%
            separate(alt_frq,
                     into=c('alt','alt_frq'),
                     sep=":(?=[0-9])",
                     convert=TRUE,
                     extra='merge') %>%
            dplyr::rename_with(~ paste0('AFR_',.x), dplyr::ends_with('_frq')),

        frqe %>%
            separate(ref_frq,
                 into=c('ref','ref_frq'),
                 sep=":(?=[0-9])",
                 convert=TRUE,
                 extra='merge') %>%
            separate(alt_frq,
                     into=c('alt','alt_frq'),
                     sep=":(?=[0-9])",
                     convert=TRUE,
                     extra='merge') %>%
            dplyr::rename_with(~ paste0('EUR_',.x), dplyr::ends_with('_frq'))
    ) %>%
        mutate(alt_high_freq = case_when(EUR_alt_frq > AFR_alt_frq ~ 'EUR',
                                         EUR_alt_frq < AFR_alt_frq ~ 'AFR'))

    print(

      frq %>%
          group_by(alt_high_freq) %>%
          summarize(N = dplyr::n()) %>%
          mutate(prop = N/sum(N))

    )

    return(frq)

}

alleleFreqMergeFast <- function(frqa_fname, frqe_fname) {

    afrqcnames <- c('seqnames','start1','end1','ref_frqa','alt_frqa')
    efrqcnames <- c('seqnames','start2','end2','ref_frqe','alt_frqe')

    frqa <- read_tsv(frqa_fname, col_names=afrqcnames)
    frqe <- read_tsv(frqe_fname, col_names=efrqcnames)

    cfrq <- cbind(
        frqa,
        frqe %>%
            select(-seqnames)
    )

    if (all(frqa$start1==frqe$start2)) {

        print('Successful allele frequency cbind!')
        cfrq <- cfrq %>%
            select(-c(start1,start2,end2)) %>%
            dplyr::rename(gwas_pos=end1) %>%
            as_tibble()

    } else {
        stop("Failed allele frequency cbind!")
    }

    return(cfrq)

}
