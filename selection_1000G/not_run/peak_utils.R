#!/usr/bin/R

############ READING FUNCTIONS #############

read_ukb_gwas_granges <- function(gwas_fname, fst=NULL, bf=NULL) {
  # for parsing UKB GWAS from Neale Lab
  cnames <- c('variant','minor_allele','minor_AF','expected_case_minor_AC','low_confidence_variant','n_complete_samples','AC','ytx','beta','se','tstat','pval')
  gwas <- read_tsv(gwas_fname, col_names=cnames) %>%
            separate(variant, c('seqnames','start','ref','alt'), sep=":", remove=TRUE, convert=TRUE) %>%
            mutate(seqnames=paste0('chr',seqnames),
                   gwas_variant_id=paste(seqnames,start,ref,alt,'b37',sep='_')) %>%
            dplyr::rename(gwas_pvalue=pval,
                          gwas_beta=beta,
                          gwas_se=se) %>%
            dplyr::select(seqnames,
                          start,
                          gwas_variant_id,
                          gwas_beta,
                          gwas_se,
                          gwas_pvalue) %>%
            as_granges(width=1)

  if (!is.null(fst)) {
    if (class(fst)!="GRanges") {
      fst <- as_granges(fst)
    }
    print("Joining FST with eQTL...")
    gwas <- join_overlap_left(gwas, fst) %>%
      dplyr::mutate(FST_eqtl=FST) %>%
      dplyr::select(-FST)
  }

  if (!is.null(bf)) {
    if (class(bf)!="GRanges") {
      bf <- as_granges(bf)
    }
    print("Joining BF with eQTL...")
    gwas <- join_overlap_left(gwas, bf) %>%
      dplyr::mutate(BF_eqtl=median_BF) %>%
      dplyr::select(-median_BF)
  }

  return(gwas)
}

read_hg19_gwas_granges <- function(gwas_fname, mhc_rm=TRUE, pthresh=5e-6, fst=NULL, bf=NULL) {
  print(paste("Loading GWAS SNPs:", gwas_fname))
  gwas <- read_tsv(gwas_fname, quote="")
  # subset to more 'categorical' phenotypes
  print("Removing eQTL, methylation QTL, and metabolomics GWAS hits")
  gwas <- gwas %>%
            subset(is.na(EqtlMethMetabStudy),
                   select=c(NHLBIkey,
                            `chr(hg19)`,
                            `pos(hg19)`,
                            `SNPid(in paper)`,
                            Pvalue,
                            Phenotype,
                            PaperPhenotypeDescription,
                            PaperPhenotypeCategories,
                            Title,
                            GWASancestryDescription,
                            InGene,
                            NearestGene,
                            dbSNPfxn)) %>%
            filter(Pvalue < pthresh) %>%
            rename(gwas_pval = Pvalue,
                    gwas_rsID = `SNPid(in paper)`) %>%
            mutate(seqnames=`chr(hg19)`,
                   start=`pos(hg19)`,
                   end=`pos(hg19)`+1,
                 gwas_id=paste(seqnames,start,sep='_')) %>%
            as_granges()

  # remove mhc before overlap if specified
  if (mhc_rm) {
    mhc <- data.frame(seqnames = 6,
                     start = 28477797,
                     end = 33448354) %>%
            as_granges
    print("Removing GWAS SNPs in MHC...")
    gwas <- gwas %>%
      filter_by_non_overlaps(mhc)
  }

  gwas <- gwas %>%
    as_tibble() %>%
    mutate(seqnames=paste0('chr',seqnames)) %>%
    as_granges()

  if (!is.null(fst)) {
    if (class(fst)!="GRanges") {
      fst <- as_granges(fst)
    }
    print("Joining FST with GWAS...")
    gwas <- join_overlap_left(gwas, fst) %>%
      dplyr::mutate(FST_gwas=FST) %>%
      dplyr::select(-FST)
  }

  if (!is.null(bf)) {
    if (class(bf)!="GRanges") {
      bf <- as_granges(bf)
    }
    print("Joining BF with GWAS...")
    gwas <- join_overlap_left(gwas, bf) %>%
      dplyr::mutate(BF_gwas=median_BF) %>%
      dplyr::select(-median_BF)
  }

  return(gwas)

}

read_deseq2 <- function(res_fname, pthresh=0.05) {
  res <- read.csv(res_fname, header=TRUE) %>% as_tibble()
  res <- res %>%
              dplyr::select(ensembl_gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,neglog10padj) %>%
              mutate(DE_significant=if_else(padj<pthresh,TRUE,FALSE),
                     DE_direction=if_else(log2FoldChange<0,'AFR','EUR')) %>%
              dplyr::rename(DE_baseMean = baseMean,
                     DE_log2FoldChange = log2FoldChange,
                     DE_lfcSE = lfcSE,
                    DE_stat = stat,
                    DE_pvalue = pvalue,
                    DE_padj = padj,
                    DE_neglog10padj = neglog10padj)
  return(res)
}

read_deseq2_krt <- function(res_fname,
                            pthresh=0.05,
                            gs_dir='/cashew/users/kade/snakemake/atac/selection_1000G/peak_norm/output/gene_sets') {

    krtz_fname           <- 'keratinization.txt'
    krt_epithelial_fname <- 'KRT_epithelial.txt'
    krt_hair_fname       <- 'KRT_hair.txt'
    krt_root_fname       <- 'KRT_root_sheath.txt'

    krtz           <- readLines(file.path(gs_dir,krtz_fname))
    krt_epithelial <- readLines(file.path(gs_dir,krt_epithelial_fname))
    krt_hair       <- readLines(file.path(gs_dir,krt_hair_fname))
    krt_root       <- readLines(file.path(gs_dir,krt_root_fname))

    if(base::endsWith(res_fname,'.csv')) {
      sep=','
    } else {
      sep='\t'
    }

    cnames <- c(
      'dropped_cov',
      'ensembl_gene_id',
      'baseMean',
      'log2FoldChange',
      'lfcSE',
      'stat',
      'pvalue',
      'padj',
      'neglog10padj'
    )

    res <- read_delim(res_fname, sep) %>%
        dplyr::select(any_of(cnames)) %>%
        mutate(
            DE_significant=if_else(padj<pthresh,TRUE,FALSE),
            DE_direction=if_else(log2FoldChange<0,'AFR','EUR'),
            krt_type = case_when(ensembl_gene_id %in% krt_root       ~ 'root sheath',
                                 ensembl_gene_id %in% krt_hair       ~ 'hair',
                                 ensembl_gene_id %in% krt_epithelial ~ 'epithelial'),
            krt_type = fct_relevel(krt_type,'epithelial','root sheath','hair'),
            gene_set = if_else(ensembl_gene_id %in% krtz,
                              'keratinization',
                              'other'),
            gene_set = fct_relevel(gene_set, rev)
        )


    return(res)

}

annotateKRT <- function(df,
                        eGene=FALSE,
                        gs_dir='/cashew/users/kade/snakemake/atac/selection_1000G/peak_norm/output/gene_sets') {

    krtz_fname           <- 'keratinization.txt'
    krt_epithelial_fname <- 'KRT_epithelial.txt'
    krt_hair_fname       <- 'KRT_hair.txt'
    krt_root_fname       <- 'KRT_root_sheath.txt'

    krtz           <- readLines(file.path(gs_dir,krtz_fname))
    krt_epithelial <- readLines(file.path(gs_dir,krt_epithelial_fname))
    krt_hair       <- readLines(file.path(gs_dir,krt_hair_fname))
    krt_root       <- readLines(file.path(gs_dir,krt_root_fname))

    if (eGene) {
        df <- df %>%
            mutate(eKrt_type = case_when(eGene %in% krt_root       ~ 'root sheath',
                                        eGene %in% krt_hair       ~ 'hair',
                                        eGene %in% krt_epithelial ~ 'epithelial'),
                   eKrt_type = fct_relevel(eKrt_type,'epithelial','root sheath','hair'),
                   eGene_set = if_else(eGene %in% krtz,
                                      'keratinization',
                                      'other'),
                   eGene_set = fct_relevel(eGene_set, rev))
    } else {
        df <- df %>%
            mutate(krt_type = case_when(ensembl_gene_id %in% krt_root       ~ 'root sheath',
                                        ensembl_gene_id %in% krt_hair       ~ 'hair',
                                        ensembl_gene_id %in% krt_epithelial ~ 'epithelial'),
                   krt_type = fct_relevel(krt_type,'epithelial','root sheath','hair'),
                   gene_set = if_else(ensembl_gene_id %in% krtz,
                                      'keratinization',
                                      'other'),
                   gene_set = fct_relevel(gene_set, rev))
    }


    return(df)

}

read_hg19_eqtl_granges_new <- function(eqtl_fname, ch_fname, gs_fname, fst=NULL, bf=NULL, de=NULL, lift_over=FALSE) {
  # add eQTLs from GTEx tissue tehaas format
  eqtl <- read_tsv(eqtl_fname, quote="") %>%
            separate(Gene_ID, c('ensembl_gene_id',NA), sep="\\.", remove=TRUE) %>%
            separate(SNP_ID, c('seqnames','start','ref','alt','build'), sep="_", remove=FALSE, convert=TRUE) %>%
            dplyr::rename(variant_id=SNP_ID) %>%
            dplyr::rename(eqtl_pval=`P-Value`) %>%
            as_granges(width=1)
  if (grepl("v8", eqtl_fname) | lift_over) {
    genome(eqtl) = "GRCh38"
    if(!(require(rtracklayer))) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

      BiocManager::install("rtracklayer")
    }
    library(rtracklayer)
    ch = import.chain(ch_fname)
    seqlevelsStyle(eqtl) = "UCSC"  # necessary
    print("Lifting over eQTL coordinates from hg38 to hg19")
    eqtl_19 = liftOver(eqtl, ch)
    eqtl_19 = unlist(eqtl_19)
    genome(eqtl_19) = "hg19"
    print(paste(length(eqtl) - length(eqtl_19),"eQTLs dropped in liftOver"))
  } else {
    eqtl_19 <- eqtl
  }
  gs <- readLines(gs_fname)
  eqtl_19 <- eqtl_19 %>%
                  as_tibble() %>%
                  dplyr::select(c(seqnames,
                                 start,
                                 end,
                                 variant_id,
                                 ensembl_gene_id,
                                 eqtl_pval)) %>%
                  dplyr::rename(eqtl_id_new=variant_id,
                        eGene=ensembl_gene_id) %>%
                  mutate(eGene_in_set=if_else(eGene %in% gs,
                                       TRUE,
                                       FALSE),
                         eqtl_id=paste(seqnames,start,sep='_')) %>%
                  as_granges()

  if (!is.null(fst)) {
    if (class(fst)!="GRanges") {
      fst <- as_granges(fst)
    }
    print("Joining FST with eQTL...")
    eqtl_19 <- join_overlap_left(eqtl_19, fst) %>%
      dplyr::mutate(FST_eqtl=FST) %>%
      dplyr::select(-FST)
  }

  if (!is.null(bf)) {
    if (class(bf)!="GRanges") {
      bf <- as_granges(bf)
    }
    print("Joining BF with eQTL...")
    eqtl_19 <- join_overlap_left(eqtl_19, bf) %>%
      dplyr::mutate(BF_eqtl=median_BF) %>%
      dplyr::select(-median_BF)
  }

  if (!is.null(de)) {
      eqtl_19 <- eqtl_19 %>%
          as_tibble %>%
          merge(.,
                de %>%
                    dplyr::select(ensembl_gene_id, DE_log2FoldChange, DE_padj, DE_significant, DE_direction) %>%
                    dplyr::rename(eGene=ensembl_gene_id),
                by='eGene',
                all.x=TRUE)

  }

  return(eqtl_19)
}

read_celltype_eqtl <- function(fname, res='low', gsfname=NULL, gsname='hair_keratin') {
    df <- read_tsv(fname) %>%
        separate(gene_id, c('ensembl_gene_id',NA), sep="\\.", remove=TRUE) %>%
        dplyr::rename(eGene_hgnc=gene_name,
                      eGene=ensembl_gene_id,
                      eGene_start=start,
                      eGene_end=end) %>%
        dplyr::select(-c(ref,alt,pos,rsid)) %>%
        separate(id, c('seqnames','start','ref','alt'), sep="_", remove=TRUE, convert=TRUE) %>%
        mutate(
            seqnames = paste0('chr',seqnames),
            eqtl_id = paste(seqnames,start,sep='_'),
            eqtl_id_new = paste(seqnames,start,ref,alt,'b37',sep='_'),
            ct_specific = if_else(n_signif==1, TRUE, FALSE),
            ct_associated = if_else(n_signif >= 1, TRUE, FALSE)
        ) %>%
        as_granges(width=1)


    if (grepl('lowres', fname)) {
        df <- df %>%
            mutate(specific_cell_type = case_when(
                    gt.collapsed_epidermal.qval < 0.1 & ct_specific ~ 'collapsed_epidermal',
                    gt.collapsed_inner_bulge.qval < 0.1 & ct_specific ~ 'collapsed_inner_bulge',
                    gt.collapsed_leukocyte.qval < 0.1 & ct_specific ~ 'collapsed_leukocyte')
                  )
    } else {
        df <- df %>%
            mutate(specific_cell_type = case_when(
                    gt.epidermis.qval < 0.1 & ct_specific ~ 'epidermis',
                    gt.outer_bulge.qval < 0.1 & ct_specific ~ 'outer_bulge',
                    gt.inner_bulge.qval < 0.1 & ct_specific ~ 'inner_bulge',
                    gt.leukocyte.qval < 0.1 & ct_specific ~ 'leukocyte',
                    gt.epidermis_stem_cell.qval < 0.1 & ct_specific ~ 'epidermis_stem_cell',
                    gt.epidermis_basal.qval < 0.1 & ct_specific ~ 'epidermis_basal')
                  )
    }

    df <- df %>% dplyr::select(-dplyr::starts_with('gt.'))

    if (!is.null(gsfname)) {
        gs <- read_lines(gsfname)
        df <- df %>% mutate(eGene_in_set=if_else(eGene %in% gs, TRUE, FALSE))
    }

    return(df)

}

read_hg19_eqtl_granges <- function(eqtl_fname, ch_fname, gs_fname, fst=NULL, bf=NULL, lift_over=FALSE) {
  # add eQTLs from GTEx tissue
  if (grepl("tehaas", eqtl_fname)) {
    eqtl <- read_tsv(eqtl_fname, quote="") %>%
            separate(Gene_ID, c('ensembl_gene_id',NA), sep="\\.", remove=TRUE) %>%
            separate(SNP_ID, c('seqnames','start','ref','alt','build'), sep="_", remove=FALSE, convert=TRUE) %>%
            dplyr::rename(variant_id=SNP_ID) %>%
            dplyr::rename(eqtl_pval=`P-Value`) %>%
            as_granges(width=1)
    lift_over <- TRUE
  } else if (grepl("GENOA", eqtl_fname)) {
    cnames <- c('chr','GENE','rs','ps','allele1','allele0','af','beta','se','p_wald')
    eqtl <- read_delim(eqtl_fname, ' ', col_names=cnames) %>%
        mutate(seqnames=paste0('chr',chr),
               variant_id=paste(seqnames,ps,allele1,allele0,'b37',sep="_")) %>%
        dplyr::rename(ensembl_gene_id=GENE,
                      start=ps,
                      slope=beta,
                      slope_se=se,
                      eqtl_pval=p_wald) %>%
        dplyr::select(seqnames,
                           start,
                           variant_id,
                           ensembl_gene_id,
                           eqtl_pval,
                           slope,
                           slope_se) %>%
        as_granges(width=1)
  } else if (grepl("Donovan", eqtl_fname)) {
    eqtl <- read_celltype_eqtl(eqtl_fname, gsfname=gs_fname)
  } else {
    eqtl <- read_tsv(eqtl_fname, quote="") %>%
              separate(gene_id, c('ensembl_gene_id',NA), sep="\\.", remove=TRUE) %>%
              separate(variant_id, c('seqnames','start','ref','alt','build'), sep="_", remove=FALSE, convert=TRUE) %>%
              dplyr::rename(eqtl_pval=pval_nominal) %>%
              as_granges(width=1)
  }

  if (grepl("v8", eqtl_fname) | lift_over) {
    genome(eqtl) = "GRCh38"
    if(!(require(rtracklayer))) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

      BiocManager::install("rtracklayer")
    }
    library(rtracklayer)
    ch = import.chain(ch_fname)
    seqlevelsStyle(eqtl) = "UCSC"  # necessary
    print("Lifting over eQTL coordinates from hg38 to hg19")
    eqtl_19 = liftOver(eqtl, ch)
    eqtl_19 = unlist(eqtl_19)
    genome(eqtl_19) = "hg19"
    print(paste(length(eqtl) - length(eqtl_19),"eQTLs dropped in liftOver"))
  } else {
    eqtl_19 <- eqtl
  }

  gs <- readLines(gs_fname)

  if (grepl("tehaas", eqtl_fname)) {
    eqtl_19 <- eqtl_19 %>%
        as_tibble() %>%
        dplyr::select(c(seqnames,
                       start,
                       end,
                       variant_id,
                       ensembl_gene_id,
                       eqtl_pval)) %>%
        dplyr::rename(eqtl_id_new=variant_id,
              eGene=ensembl_gene_id) %>%
        mutate(eGene_in_set=if_else(eGene %in% gs,
                             TRUE,
                             FALSE),
               eqtl_id=paste(seqnames,start,sep='_')) %>%
        as_granges()
  } else if (grepl("GENOA", eqtl_fname)) {
    eqtl_19 <- eqtl_19 %>%
        as_tibble() %>%
        dplyr::rename(eqtl_id_new=variant_id,
                      eGene=ensembl_gene_id,
                      eqtl_slope=slope,
                      eqtl_slope_se=slope_se) %>%
        mutate(eGene_in_set=if_else(eGene %in% gs,
                             TRUE,
                             FALSE),
               eqtl_id=paste(seqnames,start,sep='_')) %>%
        as_granges()
  } else if (grepl("Donovan", eqtl_fname)) {
    print('Cell-type-specific eQTL loaded')
  } else {
    eqtl_19 <- eqtl_19 %>%
        as_tibble() %>%
        dplyr::select(c(seqnames,
                       start,
                       end,
                       variant_id,
                       ensembl_gene_id,
                       tss_distance,
                       eqtl_pval,
                       slope,
                       slope_se)) %>%
        dplyr::rename(eqtl_id_new=variant_id,
              eGene=ensembl_gene_id,
              eGene_distance=tss_distance,
              eqtl_slope=slope,
              eqtl_slope_se=slope_se) %>%
        mutate(eGene_in_set=if_else(eGene %in% gs,
                             TRUE,
                             FALSE),
               eqtl_id=paste(seqnames,start,sep='_')) %>%
        as_granges()
  }

  if (!is.null(fst)) {
    if (class(fst)!="GRanges") {
      fst <- as_granges(fst)
    }
    print("Joining FST with eQTL...")
    eqtl_19 <- join_overlap_left(eqtl_19, fst) %>%
      dplyr::mutate(FST_eqtl=FST) %>%
      dplyr::select(-FST)
  }

  if (!is.null(bf)) {
    if (class(bf)!="GRanges") {
      bf <- as_granges(bf)
    }
    print("Joining BF with eQTL...")
    eqtl_19 <- join_overlap_left(eqtl_19, bf) %>%
      dplyr::mutate(BF_eqtl=median_BF) %>%
      dplyr::select(-median_BF)
  }

  return(eqtl_19)
}

read_caqtl_granges <- function(caqtl_fname, pop=1, fst=NULL, bf=NULL) {
  if (grepl('shared', caqtl_fname)) {
    qtls <- read_tsv(caqtl_fname,quote="")
    qtls <- qtls %>%
              mutate(seqnames=paste0("chr",chr)) %>%
              dplyr::select(-c(start,chr, peak_len)) %>%
              dplyr::rename(start=pos,
                            open_allele=open_best,
                            closed_allele=closed_best,
                            caqtl_pval=fishers_pval) %>%
              as_granges(width=1)
  } else {
    if(!(require(readxl))) install.packages("readxl")
    library(readxl)
    qtls <- read_excel(caqtl_fname, skip=0,sheet=pop,col_names=T) %>%
      dplyr::rename(seqnames=Chr,
                 start=position,
                 alt=ALTallele,
                 ref=REFallele,
                 open_allele=higherBindingAllele,
                 caqtl_pval=pvalue) %>%
      mutate(seqnames=paste0("chr",seqnames),
             closed_allele=if_else(open_allele==ref, alt, ref)) %>%
      as_granges(width=1)
  }

  qtls <- qtls %>%
    mutate(caqtl_id=paste(seqnames, start, sep="_"))

  if (!is.null(fst)) {
    if (class(fst)!="GRanges") {
      fst <- as_granges(fst)
    }
    print("Joining FST with caQTL...")
    qtls <- join_overlap_left(qtls, fst) %>%
      dplyr::mutate(FST_caqtl=FST) %>%
      dplyr::select(-FST)
  }

  if (!is.null(bf)) {
    if (class(bf)!="GRanges") {
      bf <- as_granges(bf)
    }
    print("Joining BF with caQTL...")
    qtls <- join_overlap_left(qtls, bf) %>%
      dplyr::mutate(BF_caqtl=median_BF) %>%
      dplyr::select(-median_BF)
  }

  return(qtls)

}

read_bed_noheader_granges <- function(bed) {
  cnames = c('seqnames','start','end')
  bedr <- read_tsv(bed, col_names=cnames) %>%
      as_granges()
  return(bedr)
}

read_peak_density_granges <- function(peaks_fname, fst=NULL, bf=NULL) {

  # if FST provided, adds column with maximum FST in peak

  if (str_detect(peaks_fname, '(DESeq|deseq)')) {

    pdf <- read_tsv(peaks_fname) %>%
      separate(peak_id,
             into=c('seqnames','start','end'),
             sep='_',
             remove=FALSE,
             convert=TRUE) %>%
      dplyr::rename(DA_baseMean=baseMean,
                    DA_log2FoldChange=log2FoldChange,
                    DA_lfcSE=lfcSE,
                    DA_stat=stat,
                    DA_pvalue=pvalue,
                    DA_padj=padj) %>%
      mutate(high_CA=if_else(DA_log2FoldChange < 0, 'AFR', 'EUR', missing='none'),
	     start = as.numeric(as.character(start)),
	     end = as.numeric(as.character(end))) %>% 
      dplyr::filter(!is.na(start),
		    !is.na(end)) %>%  
      as_granges()

  } else {

    pdf <- read_tsv(peaks_fname, quote="")  %>%
            mutate(peak_id=paste(chr,start,end,sep="_")) %>%
            dplyr::select(chr,start,end,peak_id,welchs_pvalue,mean_AFR,mean_EUR,high_CA) %>%
            dplyr::rename(seqnames=chr, DA_pvalue=welchs_pvalue) %>%
            as_granges()

  }

  if (!is.null(fst)) {
    if (class(fst)!="GRanges") {
      fst <- as_granges(fst)
    }
    fst <- fst %>% mutate(fst_id=paste(seqnames,start,sep="_"))
    print("Joining FST with peaks and retaining max FST per peak...")
    # SNP coordinates as 'FST_id'
    pdf <- join_overlap_left(pdf, fst) %>%
      as_tibble() %>%
      group_by(peak_id) %>%
      arrange(desc(FST), .by_group=TRUE) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      dplyr::rename(FST_peak_max=FST) %>%
      as_granges()
  }

  # add BF here

  if (!is.null(bf)) {
    if (class(bf)!="GRanges") {
      bf <- as_granges(bf)
    }
    bf <- bf %>% mutate(bf_id=paste(seqnames,start,sep="_"))
    print("Joining BF with peaks and retaining max BF per peak...")
    # SNP coordinates as 'BF_id'
    pdf <- join_overlap_left(pdf, bf) %>%
      as_tibble() %>%
      group_by(peak_id) %>%
      arrange(desc(median_BF), .by_group=TRUE) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      dplyr::rename(BF_peak_max=median_BF) %>%
      as_granges()
  }

  return(pdf)
}

read_bayenv_granges <- function(bayenv_fname, fst=NULL) {
    res <- read_tsv(bayenv_fname) %>%
        dplyr::rename(seqnames=`#CHR`,
                 start=POS,
                 median_BF=`median BF`) %>%
        mutate(seqnames=paste0('chr',seqnames)) %>%
        dplyr::select(seqnames,start,median_BF) %>%
        as_granges(width=1)

    if (!is.null(fst)) {
      if (class(fst)!="GRanges") {
        fst <- as_granges(fst)
      }
      print("Joining FST with Bayes Factors...")
      res <- join_overlap_left(res, fst) %>%
        dplyr::mutate(FST_bf=FST) %>%
        dplyr::select(-FST)
    }

    return(res)
}

read_fst_granges <- function(fst_fname, cnames=FALSE) {
    if (cnames) {
        fst <- read_tsv(fst_fname, quote='', na = c("", "NA", "NaN", "-nan")) %>%
            dplyr::rename(seqnames=CHROM,
                          start=POS,
                          FST = WEIR_AND_COCKERHAM_FST)
    } else {
        cnames <- c('seqnames','start','FST')
        ctypes <- 'cid'
        fst <- read_tsv(fst_fname,
                        quote='',
                        na = c("", "NA", "NaN", "-nan"),
                        col_names=cnames,
                        col_types=ctypes)
    }

    if (!grepl('chr', fst %>% pull(seqnames) %>% head(1))) {
        fst <- fst %>%
            mutate(seqnames=paste0('chr',seqnames))
    }

    fst <- fst %>%
        drop_na() %>%
        mutate(FST = if_else(FST < 0, 0, FST)) %>%
        as_granges(width=1)

    return(fst)
}

peaks2genes <- function(peaks, gene_locs_fname, loops_fname, gs_fname, expansion_dist=24100, collapse=TRUE) {

  gs <- readLines(gs_fname)
  if (class(peaks)=='character') {
    if (grepl('density', peaks)) {
      pdf <- read_peak_density_granges(peaks)
    } else {
      pdf <- read_bed_noheader_granges(peaks)
    }
  } else if (class(peaks)!='GRanges') {
    pdf <- as_granges(peaks)
  } else {
    pdf <- peaks
  }

  pdf_exp <- stretch(anchor_center(pdf), expansion_dist*2)

  cnames <- c('anchor_chr','anchor_start','anchor_end','seqnames','start','end')
  loops <- read_tsv(loops_fname, col_names=cnames) %>%
      as_granges() %>%
      anchor_center() %>%
      stretch(expansion_dist*2)

  cnames <- c('seqnames','start','end','ensembl_gene_id')
  ctypes <- 'ciic__'
  gene_locs <- read_tsv(gene_locs_fname, col_names=cnames, col_types=ctypes) %>%
      as_granges

  loop_genes <- join_overlap_inner(loops, gene_locs) %>%
      as_tibble() %>%
      dplyr::select(-c(seqnames,start,end,width,strand)) %>%
      dplyr::rename(seqnames=anchor_chr,
                   start=anchor_start,
                   end=anchor_end) %>%
      as_granges()

  pdflg <- join_overlap_inner(pdf_exp, loop_genes) %>%
      anchor_center() %>%
      stretch(-expansion_dist*2) %>%
      mutate(gene_linkage='loop')

  pdfng <- join_nearest(pdf, gene_locs) %>%
      mutate(gene_linkage='nearest')

  pg <- bind_ranges(pdflg, pdfng) %>%
      as_tibble() %>%
      dplyr::select(-c(width,strand)) %>%
      mutate(in_gene_set=if_else(ensembl_gene_id %in% gs, TRUE, FALSE))

  if (collapse) {
    pg <- pg %>%
      dplyr::group_by(peak_id, ensembl_gene_id) %>%
      arrange(gene_linkage, .by_group=TRUE) %>%
      dplyr::slice(1) %>%
      ungroup()
  }

  return(pg)

}

merge_new_gwas_peaks <- function(gwas_dir, pgt, rpattern="hair_*", debug=FALSE) {

    i <- 0
    fnames <- list.files(gwas_dir, pattern=rpattern)
    for (fname in fnames) {
        i <- i + 1

        if (debug & i>2) {
            pgt <- pgt %>%
                mutate_if(is.character, list(~na_if(.,'')))
            return(pgt)
        }

        pheno <- str_split(fname, '\\.')[[1]][1]

        gwas_fname <- file.path(gwas_dir,fname)
        print(paste("Reading",pheno,"GWAS..."))
        gwas <- read_ukb_gwas_granges(gwas_fname, fst=NULL, bf=NULL) %>%
            mutate(gwas_pheno=pheno)

        print(paste("Joining",pheno,"GWAS with peaks..."))
        pgt <- join_overlap_left(pgt %>% as_granges, gwas) %>%
            as_tibble %>%
            group_by(peak_id) %>%
            arrange(gwas_pvalue, .by_group=TRUE) %>%
            dplyr::slice(1) %>%
            ungroup() %>%
            mutate(gwas_pvalue=as.character(gwas_pvalue),
                   gwas_beta=as.character(gwas_beta),
                   gwas_se=as.character(gwas_se))

        if (i==1) {
           pgt <- pgt %>%
                dplyr::rename(new_gwas_pvalue=gwas_pvalue,
                          new_gwas_variant_id=gwas_variant_id,
                          new_gwas_beta=gwas_beta,
                          new_gwas_se=gwas_se,
                          new_gwas_pheno=gwas_pheno)
        } else {

            pgt <- pgt  %>%
                unite('new_gwas_pheno', dplyr::ends_with('gwas_pheno'), sep=';', na.rm=TRUE) %>%
                unite('new_gwas_pvalue', dplyr::ends_with('gwas_pvalue'), sep=';', na.rm=TRUE) %>%
                unite('new_gwas_variant_id', dplyr::ends_with('gwas_variant_id'), sep=';', na.rm=TRUE) %>%
                unite('new_gwas_beta', dplyr::ends_with('gwas_beta'), sep=';', na.rm=TRUE) %>%
                unite('new_gwas_se', dplyr::ends_with('gwas_se'), sep=';', na.rm=TRUE) %>%
                mutate_if(is.character, list(~na_if(.,'')))
        }
    }

    return(pgt)

}

merge_new_gwas_snps <- function(gwas_dir, pgt, rpattern="hair_*", debug=FALSE) {

    i <- 0
    fnames <- list.files(gwas_dir, pattern=rpattern)
    for (fname in fnames) {
        i <- i + 1

        if (debug & i>2) {
            pgt <- pgt %>%
                mutate_if(is.character, list(~na_if(.,'')))
            return(pgt)
        }

        pheno <- str_split(fname, '\\.')[[1]][1]

        gwas_fname <- file.path(gwas_dir,fname)
        print(paste("Reading",pheno,"GWAS..."))
        gwas <- read_ukb_gwas_granges(gwas_fname, fst=NULL, bf=NULL) %>%
            mutate(gwas_pheno=pheno)

        print(paste("Joining",pheno,"GWAS with QTL..."))
        pgt <- join_overlap_left(pgt %>% as_granges, gwas) %>%
            as_tibble %>%
            mutate(gwas_pvalue=as.character(gwas_pvalue),
                   gwas_beta=as.character(gwas_beta),
                   gwas_se=as.character(gwas_se))

        if (i==1) {
           pgt <- pgt %>%
                dplyr::rename(new_gwas_pvalue=gwas_pvalue,
                          new_gwas_variant_id=gwas_variant_id,
                          new_gwas_beta=gwas_beta,
                          new_gwas_se=gwas_se,
                          new_gwas_pheno=gwas_pheno)
        } else {

            pgt <- pgt  %>%
                unite('new_gwas_pheno', dplyr::ends_with('gwas_pheno'), sep=';', na.rm=TRUE) %>%
                unite('new_gwas_pvalue', dplyr::ends_with('gwas_pvalue'), sep=';', na.rm=TRUE) %>%
                unite('new_gwas_variant_id', dplyr::ends_with('gwas_variant_id'), sep=';', na.rm=TRUE) %>%
                unite('new_gwas_beta', dplyr::ends_with('gwas_beta'), sep=';', na.rm=TRUE) %>%
                unite('new_gwas_se', dplyr::ends_with('gwas_se'), sep=';', na.rm=TRUE) %>%
                mutate_if(is.character, list(~na_if(.,'')))
        }
    }

    return(pgt)

}

################# REFORMATTING / SUBSETTING FUNCTIONS ################

liftover_grange <- function(df, ch_fname, from_b="GRCh38", to_b="hg19") {
    genome(df) = from_b
    if(!(require(rtracklayer))) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

      BiocManager::install("rtracklayer")
    }
    library(rtracklayer)
    ch = import.chain(ch_fname)
    seqlevelsStyle(df) = "UCSC"  # necessary
    print(paste0("Lifting over eQTL coordinates from ", from_b, " to ", to_b, "..."))
    df_new = liftOver(df, ch)
    df_new = unlist(df_new)
    genome(df_new) = to_b
    nc <- length(df_new) - length(df)
    print(paste0("Ranges gained/lost in liftOver = ", if_else(nc>0,'+',''), nc))
    return(df_new)
}

top_genes_snps_per_peak <- function(df, top_snps=TRUE) {
  if (top_snps) {
    dfs <- df %>%
        group_by(peak_id) %>%
        arrange(desc(eGene_in_set),
                desc(in_gene_set),
                desc(gwas_Hair),
                desc(gwas_Hair_morphology),
                caqtl_pval,
                eqtl_pval,
                gwas_pval,
                DE_padj,
                .by_group=TRUE) %>%
        dplyr::slice(1) %>%
        ungroup()
  } else {
    dfs <- df %>%
        group_by(peak_id) %>%
        arrange(desc(in_gene_set)) %>%
        # mutate to indicator for any candidate gene linked to peak in gene set
        mutate(in_gene_set=dplyr::first(in_gene_set)) %>%
        ungroup() %>%
        dplyr::select(-c(ensembl_gene_id,gene_linkage,starts_with('DE_'))) %>%
        distinct()
  }

  return(dfs)
}

gwas_subset <- function(df, pheno="Hair", pthresh=0.05) {
    dfp <- df %>%
        filter(str_detect(Phenotype, pheno)
           | str_detect(PaperPhenotypeDescription, pheno)
           | str_detect(PaperPhenotypeCategories, pheno),
           DA_pvalue<pthresh) %>%
        arrange(DA_pvalue)

    return(dfp)
}

add_gwas_col <- function(df, pheno='Hair') {
    if (class(df)=='GRanges') df <- df %>% as_tibble()
    phenostr <- sub(' ', '_', pheno)
    cname <- paste('gwas', phenostr, sep='_')
    df[cname] <- if_else(str_detect(df$Phenotype, pheno)
                        | str_detect(df$PaperPhenotypeDescription, pheno)
                        | str_detect(df$PaperPhenotypeCategories, pheno),
                        TRUE,
                        FALSE)
    df <- df %>% as_granges()
    return(df)
}

max_fst_per_peak <- function(df) {

  dfm <- df %>%
    group_by(peak_id) %>%
    arrange(desc(FST), .by_group=TRUE) %>%
    mutate(FST=dplyr::first(FST)) %>%
    group_by(ensembl_gene_id, add=TRUE) %>%
    dplyr::slice(1) %>%
    ungroup()

  return(dfm)
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
                                   columns = c("SYMBOL", "GENEID")) %>%
        as_tibble() %>%
        dplyr::rename(ensembl_gene_id = GENEID,
                      hgnc_symbol = SYMBOL) %>%
        dplyr::filter(grepl('ENSG', ensembl_gene_id)) # filter out alternate gene symbols

    print("Adding hgnc_symbols...")
    df <- merge(df, gene_info, by = merge_col, all.x = TRUE)
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
        gene_info <- ensembldb::select(EnsDb.Hsapiens.v79, keys = ensgenes,
            keytype = "GENEID", columns = c("SYMBOL", "GENEID")) %>%
            as_tibble() %>% dplyr::rename(ensembl_gene_id = GENEID,
            hgnc_symbol = SYMBOL)
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

caqtl_peaks_by_pval <- function(df, pvs=c(1,.5,0.05,0.005,5e-4,5e-5,0), longform=FALSE) {
    dfs <- df %>%
        # top_genes_snps_per_peak(., top_snps=FALSE) %>%
        mutate(peak_width=end-start) %>%
        dplyr::filter(!is.na(DA_pvalue)) %>%
        group_by(pvalue=cut(DA_pvalue, breaks=pvs, include.lowest=TRUE)) %>%
        summarize(peaks=sum(peak_width[!duplicated(peak_id)])/1000,
                  caqtl_peaks=sum(peak_width[!is.na(caqtl_pval) & !duplicated(peak_id)])/1000,
                  eqtl_peaks=sum(peak_width[!is.na(eqtl_pval) & !duplicated(peak_id)])/1000,
                  caqtl_eqtl_peaks=sum(peak_width[!is.na(caqtl_pval) & !is.na(eqtl_pval) & !duplicated(peak_id)])/1000,
                  ker_peaks=sum(peak_width[(in_gene_set | eGene_in_set) & !duplicated(peak_id)], na.rm=TRUE)/1000,
                  looped_near_ker_peaks=sum(peak_width[in_gene_set & !duplicated(peak_id)], na.rm=TRUE)/1000,
                  eGene_ker_peaks=sum(peak_width[eGene_in_set & !duplicated(peak_id)], na.rm=TRUE)/1000,
                 peaks_caQTL=sum(!is.na(caqtl_pval)),
                 peaks_caQTL_norm=peaks_caQTL/peaks,
                 peaks_eQTL=sum(!is.na(eqtl_pval)),
                 peaks_eQTL_norm=peaks_eQTL/peaks,
                 peaks_caQTL_eQTL=sum(!is.na(caqtl_pval) & !is.na(eqtl_pval)),
                 peaks_caQTL_eQTL_norm=peaks_caQTL_eQTL/peaks,
                 peaks_ker=sum(in_gene_set | eGene_in_set, na.rm=TRUE),
                 peaks_looped_near_ker=sum(in_gene_set, na.rm=TRUE),
                 peaks_caQTL_ker=sum(!is.na(caqtl_pval) & (in_gene_set | eGene_in_set), na.rm=TRUE),
                 peaks_caQTL_ker_norm=peaks_caQTL_ker/caqtl_peaks,
                 peaks_caQTL_ker_norm_ker=peaks_caQTL_ker/ker_peaks,
                 peaks_caQTL_looped_near_ker=sum(!is.na(caqtl_pval) & in_gene_set, na.rm=TRUE),
                 peaks_caQTL_looped_near_ker_norm=peaks_caQTL_looped_near_ker/caqtl_peaks,
                 peaks_caQTL_looped_near_ker_norm_ker=peaks_caQTL_looped_near_ker/looped_near_ker_peaks,
                 peaks_eGene_ker=sum(eGene_in_set, na.rm=TRUE),
                 peaks_eGene_ker_norm=peaks_eGene_ker/eqtl_peaks,
                 peaks_caQTL_eGene_ker=sum(!is.na(caqtl_pval) & eGene_in_set, na.rm=TRUE),
                 peaks_caQTL_eGene_ker_norm=peaks_caQTL_eGene_ker/caqtl_eqtl_peaks,
                 peaks_caQTL_eGene_ker_norm_ker=peaks_caQTL_eGene_ker/eGene_ker_peaks) %>%
        mutate(pvalue=fct_relevel(pvalue, rev))

    if (longform) {
        dfs <- dfs %>%
            pivot_longer(cols=dplyr::starts_with('peaks_'),
                         names_to='peak_type',
                         names_prefix='peaks_',
                         values_to='count')
    }

    return(dfs)
}

gwas_peaks_by_pval <- function(df, pvs=c(1,.5,0.05,0.005,5e-4,5e-5,0), longform=FALSE) {
    dfs <- df %>%
        group_by(peak_id) %>%
        arrange(desc(in_gene_set)) %>%
        # mutate to indicator for any candidate gene linked to peak in gene set
        mutate(in_gene_set=dplyr::first(in_gene_set)) %>%
        ungroup() %>%
        dplyr::select(-c(ensembl_gene_id,gene_linkage,starts_with('DE_'))) %>%
        distinct() %>%
        mutate(peak_width=end-start) %>%
        group_by(pvalue=cut(welchs_pvalue, breaks=pvs, include.lowest=TRUE)) %>%
        summarize(peaks=sum(peak_width)/1000,
                  gwas_peaks=sum(peak_width[!is.na(gwas_pval)])/1000,
                  gwas_Hair_peaks=sum(peak_width[!is.na(gwas_pval) & gwas_Hair])/1000,
                  peaks_gwas_all=sum(!is.na(gwas_pval)),
                  peaks_gwas_all_norm=peaks_gwas_all/peaks,
                  peaks_gwas_Hair=sum(gwas_Hair, na.rm=TRUE),
                  peaks_gwas_Hair_norm=peaks_gwas_Hair/gwas_peaks,
                  peaks_gwas_Hair_morphology=sum(gwas_Hair_morphology, na.rm=TRUE),
                  peaks_gwas_Hair_morphology_norm=peaks_gwas_Hair_morphology/gwas_peaks,
                  peaks_gwas_Hair_morphology_norm_hair=peaks_gwas_Hair_morphology/gwas_Hair_peaks) %>%
        mutate(pvalue=fct_relevel(pvalue, rev))

    if (longform) {
        dfs <- dfs %>%
            pivot_longer(cols=dplyr::starts_with('peaks_'),
                         names_to='peak_type',
                         names_prefix='peaks_',
                         values_to='count')
    }

    return(dfs)
}

add_pval_bins <- function(df, pvs=c(1,.5,0.05,0.005,5e-4,5e-5,5e-6,5e-7,0)) {
    dfb <- df %>%
        dplyr::filter(!is.na(DA_pvalue)) %>%
        group_by(pvalue_bin=cut(DA_pvalue, breaks=pvs, include.lowest=TRUE)) %>%
        ungroup() %>%
        mutate(pvalue_bin=fct_relevel(pvalue_bin, rev),
              any_gene_in_set=if_else(eGene_in_set | in_gene_set, TRUE, FALSE, missing=FALSE))
    return(dfb)

}

filter_bin_peaks <- function(df, pvs=c(1,.5,0.05,0.005,0), sub='', top_snps=TRUE, stat='BF') {
  # sub should be vector containing any combination of 'caQTL', 'eQTL'
  # plot all bayes_factors in each peak

  #should be subset before
  # df <- top_genes_snps_per_peak(df, top_snps=top_snps)

  if (stat=='BF') {
    df <- df %>% dplyr::filter(!is.na(BF_peak_max))
  }
  if (stat=='FST') {
    df <- df %>% dplyr::filter(!is.na(FST_peak_max))
  }
  dfb <- add_pval_bins(df, pvs=pvs) # adds 'any_gene_in_set' column (peaks linked to eGenes/looped/nearest in GS)

  if (('caQTL' %in% sub)  & ('eQTL' %in% sub)) {
    dfb <- dfb %>%
      dplyr::filter(!is.na(caqtl_pval) & !is.na(eqtl_pval)) %>%
      # for peaks with both caQTL and eQTL, take whichever FST is greater
      mutate(FST=base::pmax(FST_caqtl, FST_eqtl, na.rm=TRUE),
             bayes_factor=base::pmax(BF_caqtl, BF_eqtl, na.rm=TRUE)) # for FST plots
  } else if ('caQTL' %in% sub) {
    dfb <- dfb %>%
      dplyr::filter(!is.na(caqtl_pval)) %>%
      mutate(FST=FST_caqtl,
             bayes_factor=BF_caqtl)
  } else if ('eQTL' %in% sub) {
    dfb <- dfb %>%
      dplyr::filter(!is.na(eqtl_pval)) %>%
      mutate(FST=FST_eqtl,
             bayes_factor=BF_eqtl)
  } else {
    dfb <- dfb %>%
      mutate(FST=FST_peak_max,
             bayes_factor=BF_peak_max)
  }

  return(dfb)
}

wilcox_gslink_in_peak_bins <- function(df, stat='BF', gene_set_link='any', top_snps=TRUE, pvs=c(1,.5,0.05,0.005,5e-4,5e-5,0), sub='') {
  # stat = 'BF' or 'FST'
  # gene_set_link = 'any', 'eGene', or 'loopednear'
  # ... = arguments passed to filter_bin_peaks
  if(!(require(rstatix)   )) install.packages("rstatix")
  library(rstatix)

  df <- filter_bin_peaks(df, pvs=pvs, sub=sub, top_snps=top_snps, stat=stat)

  dfs <- df %>% dplyr::filter(DA_pvalue < pvs[length(pvs)-1])
  if (gene_set_link=='any') {
      ngsp <- dfs %>% pull(any_gene_in_set) %>% sum(.,na.rm=TRUE)
  } else if (gene_set_link=='eGene') {
      ngsp <- dfs %>% pull(eGene_in_set) %>% sum(.,na.rm=TRUE)
  } else {
      ngsp <- dfs %>% pull(in_gene_set) %>% sum(.,na.rm=TRUE)
  }

  if (ngsp < 1) {
    return(NULL)
  }

  df <- df %>%
    group_by(pvalue_bin)

  w <- tryCatch(
    {
      if (stat=='BF') {

        if (gene_set_link=='any') {
          w <- df %>% wilcox_test(bayes_factor ~ any_gene_in_set, p.adjust.method='none')
        } else if (gene_set_link=='eGene') {
          w <- df %>% wilcox_test(bayes_factor ~ eGene_in_set, p.adjust.method='none')
        } else {
          w <- df %>% wilcox_test(bayes_factor ~ in_gene_set, p.adjust.method='none')
        }

      } else {

        if (gene_set_link=='any') {
          w <- df %>% wilcox_test(FST ~ any_gene_in_set, p.adjust.method='none')
        } else if (gene_set_link=='eGene') {
          w <- df %>% wilcox_test(FST ~ eGene_in_set, p.adjust.method='none')
        } else {
          w <- df %>% wilcox_test(FST ~ in_gene_set, p.adjust.method='none')
        }

      }
    },
    error=function(e) {
      message(e)
      return(NULL)
    }
  )

  return(w)

}

gprofiler_diffatac <- function(df, highCA='EUR', gene_link='loopednearest', pthresh=0.05, ordered=TRUE, custom_bg=FALSE, debug=FALSE, outroot='diffATAC') {
    # custom_bg = TRUE combines EUR and AFR high CA-linked sig genes into background
    # gene_set_link = 'any', 'eGene', 'loopednearest'

    outbase = paste(outroot,
                    highCA,'highCA',
                    'p',pthresh,
                    gene_link,'gene_link',
                    if_else(ordered, 'ordered',''),
                    if_else(custom_bg, 'bgAFREUR', 'bgALL'),
                    sep="_")
    if(!(require(gprofiler2))) install.packages("gprofiler2")
    library(gprofiler2)

    if ('DA_padj' %in% colnames(df)) {
        dfs <- df %>%
            dplyr::filter(DA_padj < pthresh) %>%
            arrange(DA_padj)
    } else if ('DA_pvalue' %in% colnames(df)) {
        dfs <- df %>%
            dplyr::filter(DA_pvalue < pthresh) %>%
            arrange(DA_pvalue)
    } else {
        dfs <- df %>%
            dplyr::filter(welchs_pvalue < pthresh) %>%
            arrange(welchs_pvalue)
    }

    if (gene_link=='loopednearest') {
        fg_genes <- dfs %>%
            dplyr::filter(high_CA==highCA) %>%
            pull(ensembl_gene_id) %>%
            unique()
        bg_genes <- dfs %>%
            dplyr::filter(high_CA=='EUR' | high_CA=='AFR') %>%
            pull(ensembl_gene_id) %>%
            unique()
    } else if (gene_link=='eGene') {
        fg_genes <- dfs %>%
            dplyr::filter(high_CA==highCA,
                          !is.na(eqtl_pval)) %>%
            pull(eGene) %>%
            unique()
        bg_genes <- dfs %>%
            dplyr::filter(high_CA=='EUR' | high_CA=='AFR',
                          !is.na(eqtl_pval)) %>%
            pull(eGene) %>%
            unique()
    } else {
        # combine eGene and ensembl_gene_id (loopednearest) into single column
        dfs <- dfs %>%
            dplyr::rename(gene_eGene=eGene,
                          gene_loopednearest=ensembl_gene_id) %>%
            pivot_longer(cols=c(gene_eGene, gene_loopednearest),
                         names_to='gene_link',
                         names_prefix='gene_',
                         values_to='ensembl_gene_id',
                         values_drop_na=TRUE)

        fg_genes <- dfs %>%
            dplyr::filter(high_CA==highCA) %>%
            pull(ensembl_gene_id) %>%
            unique()
        bg_genes <- dfs %>%
            dplyr::filter(high_CA=='EUR' | high_CA=='AFR') %>%
            pull(ensembl_gene_id) %>%
            unique()
    }

    print(paste0(highCA, ' high CA foreground genes at P<', pthresh,':'))
    print(length(fg_genes))

    if (custom_bg) {
        custom_bg <- bg_genes
        print(paste0('AFR + EUR high CA background genes at P<', pthresh,':'))
        print(length(bg_genes))
    	if (length(bg_genes) > 15000 | length(fg_genes) > 11000) {
    		print("Too many genes for g:Profiler to test quickly")
    		write('To many genes for efficient testing',      paste0(outbase,'_too_many_genes_gprofiler.txt'))
    		write('To many genes for efficient testing',      paste0(outbase,'_too_many_genes_gprofiler_summary.txt'))
    		return(NULL)
    	}
    } else {
        custom_bg <- NULL
        print('Background = all genes')
    }

    print("Running g:Profiler query...")

    gostres <- tryCatch(
      {
        gostres <- gost(fg_genes, organism = "hsapiens", ordered_query = ordered,
                      multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                      measure_underrepresentation = FALSE, evcodes = TRUE,
                      user_threshold = 0.05, correction_method ="g_SCS",
                      domain_scope = "known", custom_bg = custom_bg,
                      numeric_ns = "", sources = NULL)
      },
      error=function(e) {
        message(e)
        return(NULL)
      }
    )

    if (is.null(gostres)) {
      print(paste(outbase, "bad request"))
      return(NULL)
    }


    if (debug) {
        return(gostres)
    }

    if (is.null(gostres$result)) {
      write('No results to show',      paste0(outbase,'_no_sig_res_gprofiler.txt'))
      write('No results to show',      paste0(outbase,'_no_sig_res_gprofiler_summary.txt'))
    } else {
      gostdf <- gostres$result %>%
          as_tibble() %>%
          dplyr::select(-c(parents)) %>% # remove list column before writing
          arrange(source,p_value)

      gostdf_summ <- gostdf %>%
          dplyr::select(source,
                        term_name,
                        p_value,
                        query_size,
                      term_size,
                    intersection_size)


      print(gostdf_summ %>% dplyr::filter(source=='GO:BP'))

      write_tsv(gostdf,      paste0(outbase,'_gprofiler.txt'))
      write_tsv(gostdf_summ, paste0(outbase,'_gprofiler_summary.txt'))
    }

}

##################### STAT FUNCTIONS ##################

keratinization_keratinocyte_DA_enrich <- function(df, peak_thresh=5e-8, highca=FALSE) {

    df <- df %>% dplyr::filter(DA_pvalue < peak_thresh)

    if(highca=="EUR" | highca=="AFR") {

        fg_numerator <- df %>%
                        dplyr::filter(keratinocyte_peak,
                              in_gene_set,
                              high_CA==highca) %>%
                        nrow()

        fg_denominator <- df %>%
                            dplyr::filter(keratinocyte_peak,
                              high_CA==highca) %>%
                            nrow()

        bg_numerator <- df %>%
                            dplyr::filter(in_gene_set,
                                  high_CA==highca) %>%
                             nrow()

        bg_denominator <- df %>%
                            dplyr::filter(high_CA==highca) %>%
                            nrow()

    } else {

        fg_numerator <- df %>%
                        dplyr::filter(keratinocyte_peak,
                              in_gene_set) %>%
                        nrow()

        fg_denominator <- df %>%
                            dplyr::filter(in_gene_set) %>%
                            nrow()

        bg_numerator <- df %>%
                            dplyr::filter(in_gene_set) %>%
                             nrow()

        bg_denominator <- df %>%
                            nrow()

    }


    fe <- (fg_numerator/fg_denominator) / (bg_numerator/bg_denominator)
    pv <- phyper(fg_numerator, bg_numerator, bg_denominator-bg_numerator, fg_denominator, lower.tail=FALSE)

    return(
        list(enrichment = fe, pvalue = pv,
          fg_success=fg_numerator, fg_total=fg_denominator, bg_success=bg_numerator, bg_total=bg_denominator)
    )

}

################## PLOTTING FUNCTIONS #################

cor_label <- function(df, col1=2, col2=3, method='pearson') {
    require(grid)
    require(gridExtra)
    cres <- cor.test(df[[col1]], df[[col2]], method=method)
    groblab <- grobTree(
        textGrob(
            paste0(
                paste0(method," corr. = "),
                round(
                    cres$estimate,
                    4
                ),
                ",\nP = ",
                formatC(cres$p.value, format = "e", digits = 4)
            ),
            x = 0.3,
            y = 0.1,
            hjust = 0,
            gp = gpar(col = "gray",
                      fontsize = 20,
                      fontface = "italic")
        )
    )

    return(groblab)
}

plot_peaks_by_pval <- function(df=NULL, summ_df=NULL, pvs=c(1,.5,0.05,0.005,5e-4,5e-5,0), ptypes=c('caQTL','eQTL')) {

    if (is.null(summ_df)) {
        if (any(str_detect(ptypes, 'gwas'))) {
            peak_summ <- gwas_peaks_by_pval(df, pvs=pvs, longform=TRUE)
        } else {
            peak_summ <- caqtl_peaks_by_pval(df, pvs=pvs, longform=TRUE)
        }
    } else {
        peak_summ <- summ_df
    }

    if (all(str_detect(ptypes, 'norm'))) {

        p <- peak_summ %>% dplyr::filter(peak_type %in% ptypes) %>%
            ggplot( aes(x=pvalue, y=count, color=peak_type, group=peak_type))

    } else {

        p <- peak_summ %>% dplyr::filter(peak_type %in% ptypes) %>%
            ggplot( aes(x=pvalue, y=(count/peaks), color=peak_type, group=peak_type))
    }

    if (any(str_detect(ptypes, 'gwas'))) {
        colorlabel <- 'Phenotype'
        ylabel <- 'GWAS hits per 1 Kb peak'
    } else {
        colorlabel <- 'QTL type'
        ylabel <- 'QTL per 1 Kb peak'
    }

    p +
        geom_line() +
        geom_point( aes(size=peaks) ) +
        theme_classic(15) +
        labs(x='p-value bin',
            y=ylabel,
            size='1 Kb peaks in bin',
            color=colorlabel) +
        scale_color_manual(values=c('black','dark blue','dark orange', 'green', 'yellow', 'pink')) +
        theme(title = element_text(size=28)) +
        theme(axis.text.x = element_text(size=18,angle=45,hjust=1)) +
        theme(axis.text.y = element_text(size=20)) +
        theme(legend.text = element_text(size=18)) +
        theme(legend.title = element_text(size=22)) +
        theme(axis.title.x = element_text(size=26)) +
        theme(axis.title.y = element_text(size=26))
}

plot_selection_stat_peaks_by_pval <- function(df,
                                              pvs=c(1,.5,0.05,0.005,0),
                                              sub='allsnps',
                                              gene_set_link='loopednearest',
                                              stat='BF',
                                              top_snps=TRUE) {
    # sub should be vector containing any combination of 'caQTL', 'eQTL'
    # to plot only peaks having these features

    # gene_set_link can be 'any', 'eGene', or 'loopednear'

    dfb <- filter_bin_peaks(df, pvs=pvs, sub=sub, top_snps=top_snps, stat=stat)

    if (stat=='BF') {
      ylabel = "Association with latitude (Log Median Bayes Factor)"
    } else {
      ylabel = "FST"
    }

    if (gene_set_link=='any') {

      if (stat=='BF') {
        p <- dfb %>%
            ggplot(aes(pvalue_bin, log10(bayes_factor), fill=any_gene_in_set))
      } else {
        p <- dfb %>%
            ggplot(aes(pvalue_bin, FST, fill=any_gene_in_set))
      }

      p <- p +
#         geom_boxplot(outlier.size = 0, outlier.stroke=0, outlier.alpha=0) +
        geom_boxplot(outlier.alpha=0.2) +
        geom_point(pch=19,
                   position = position_jitterdodge(jitter.width=0.3,dodge.width=0.75),
                   aes(color=any_gene_in_set, alpha=any_gene_in_set))

    } else if (gene_set_link=='eGene') {

      if (stat=='BF') {
        p <- dfb %>%
            ggplot(aes(pvalue_bin, log10(bayes_factor), fill=eGene_in_set))
      } else {
        p <- dfb %>%
            ggplot(aes(pvalue_bin, FST, fill=eGene_in_set))
      }

      p <- p +
#         geom_boxplot(outlier.size = 0, outlier.stroke=0, outlier.alpha=0) +
        geom_boxplot(outlier.alpha=0.2) +
        geom_point(pch=19,
                   position = position_jitterdodge(jitter.width=0.3,dodge.width=0.75),
                   aes(color=eGene_in_set, alpha=eGene_in_set))

    } else {

      if (stat=='BF') {
        p <- dfb %>%
            ggplot(aes(pvalue_bin, log10(bayes_factor), fill=in_gene_set))
      } else {
        p <- dfb %>%
            ggplot(aes(pvalue_bin, FST, fill=in_gene_set))
      }

      p <- p +
#         geom_boxplot(outlier.size = 0, outlier.stroke=0, outlier.alpha=0) +
        geom_boxplot(outlier.alpha=0.2) +
        geom_point(pch=19,
                   position = position_jitterdodge(jitter.width=0.3,dodge.width=0.75),
                   aes(color=in_gene_set, alpha=in_gene_set))
    }

    p +
        theme_classic(14) +
        guides(fill=FALSE, color=FALSE) +
        scale_fill_manual(values=c("light gray", "pink")) +
        scale_color_manual(values=c("black", "red")) +
        scale_alpha_manual(values=c(0, 1)) +
        labs(x = "diffATAC p-value bin",
             y = ylabel) +
            theme(axis.text.x = element_text(size=18,angle=45,hjust=1)) +
            theme(axis.text.y = element_text(size=20)) +
            theme(legend.text = element_text(size=18)) +
            theme(legend.title = element_text(size=22)) +
            theme(axis.title.x = element_text(size=26)) +
            theme(axis.title.y = element_text(size=26))

}

plot_peaks_by_pval_loop <- function(df,
                                    stat='QTL',
                                    sub='',
                                    gene_set_link='loopednearest',
                                    top_snps=TRUE,
                                    ptypes=c('caQTL_looped_near_ker','caQTL_eGene_ker','caQTL_ker'),
                                    outroot='peaks_by_pval') {

    # stat = 'BF', 'FST', 'QTL'

    if (stat=='BF' | stat=='FST') {
      outbase <- paste(outroot,
                       if_else(top_snps, 'max', 'all'),
                       stat,
                       paste(sub, collapse="_"),
                       gene_set_link,
                       'gslink',
                       sep='_')
    } else {
      outbase <- paste(outroot, paste(ptypes, collapse='_'), sep='_')
    }


    binlist <- list(
        c(1,.5,.05,0),
        c(1,.5,.05,.005,0),
        c(1,.5,0.05,0.005,5e-4,0),
        c(1,.5,0.05,0.005,5e-4,5e-5,0),
        c(1,.5,0.05,0.005,5e-4,5e-5,5e-6,0),
        c(1,.5,0.05,0.005,5e-4,5e-5,5e-6,5e-7,0)
    )

    i <- 0
    for (bins in binlist) {
        i <- i + 1
        if (stat=='BF' | stat=='FST') {
          plot_selection_stat_peaks_by_pval(df,
                                            pvs=bins,
                                            sub=sub,
                                            gene_set_link=gene_set_link,
                                            stat=stat,
                                            top_snps=top_snps)
          wilcox_res <- wilcox_gslink_in_peak_bins(df,
                                                   stat=stat,
                                                   gene_set_link=gene_set_link,
                                                   pvs=bins,
                                                   sub=sub,
                                                   top_snps=top_snps)

          if (!is.null(wilcox_res)) {
              out_fname <- tryCatch(
                {
                  wilcox_res %>%
                    write_tsv(paste0(outbase, '_wilcox_', i, '.txt'))
                  out_fname <- paste0(outbase, '_wilcox_', i, '.txt')
                },
                error=function(e) {
                  message(e)
                  return(paste0(outbase, '_wilcox_', i, '.txt'))
                }
              )

              print(out_fname)

          }

        } else {
          plot_peaks_by_pval(df=df, pvs=bins, ptypes=ptypes)
        }
        ggsave(paste0(outbase,'_',i,'.png'), scale=1.5)
    }

}
