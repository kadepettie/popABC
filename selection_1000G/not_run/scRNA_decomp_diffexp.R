#!/usr/bin/R

if (!(require(limma))) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("limma")
}
if (!(require(limma))) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("edgeR")
}


if(!(require(argparse) )) install.packages("argparse")
if(!(require(tidyverse) )) install.packages("tidyverse")
if(!(require(reshape2) )) install.packages("reshape2")
if(!(require(grid) )) install.packages("grid")
if(!(require(gridExtra) )) install.packages("gridExtra")
if(!(require(Seurat) )) install.packages("Seurat")
if(!(require(textTinyR) )) install.packages("textTinyR")
if(!(require(pbapply) )) install.packages("pbapply")

library(argparse)
library(limma)
library(edgeR)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(Seurat)
library(textTinyR)
library(pbapply)


# functions

listToTibble <- function(x, cnames=NULL) {
    x <- as.data.frame(t(as_tibble(x)))
    names(x) <- c('V2')
    x$V1 <- rownames(x)
    x <- as_tibble(x[,c(2,1)])
    if (!is.null(cnames)) {
        colnames(x) <- cnames
    }
    return(x)
}

nonzeroMedian <- function(x) {
    N_nz <- sum(x > 0)
    if (N_nz==0) {
        m <- NA
    } else {
        m <- median(x[x>0])
    }
    return(m)
}

matchMeta <- function(r,m) {
    correct_names <- colnames(r)
    m <- m[rownames(m) %in% correct_names,]
    return(m)
}
matchExpr <- function(r,m) {
    ordn <- rownames(m)
    if(length(ordn) == dim(r)[2]){
        r <- r[ordn]
    }else{
        correct_names <- colnames(r)
        m <- m[rownames(m) %in% correct_names,]
        ordn <- rownames(m)
        r <- r[ordn]
    }
    return(r)
}

removeLowCounts <- function(ct, r, vr, decomp=NULL) {
    if (is.null(decomp)) {
        decomp <- r
    }
    if(ct == "CD4_T" || ct == "monocytes" || ct == "monocytes_combined"){
        tab = data.frame(genes = rownames(r), medians=apply(vr$E, 1, median), order = 1:nrow(r))
        tab = tab[order(-tab$medians), ]
        tab$order_by_median = 1:nrow(tab)
        tab = tab[order(tab$order), ]
        decomp <- decomp[which(tab$medians > 1.5), ]
    }else if(ct == "B" || ct == "CD8_T"){
        tab = data.frame(genes = rownames(r), medians=apply(vr$E, 1, median), order = 1:nrow(r))
        tab = tab[order(-tab$medians), ]
        tab$order_by_median = 1:nrow(tab)
        tab = tab[order(tab$order), ]
        decomp <- decomp[which(tab$medians > 2.5), ]
    }else{
        tab = data.frame(genes = rownames(r), medians=apply(vr$E, 1, median), order = 1:nrow(r))
        tab = tab[order(-tab$medians), ]
        tab$order_by_median = 1:nrow(tab)
        tab = tab[order(tab$order), ]
        decomp <- decomp[which(tab$medians > 4.0), ]
        }
    return(decomp)
}

voomCaptureCorrect <- function(r, m) {

    dge <- DGEList(counts = r)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~ 0 + capture, data = m)
    v <- voom(dge, design, plot = FALSE)

    fit <- lmFit(v, design)
    fit <- eBayes(fit)

    ## get residuals (intercept == 0) to regress out capture effect
    residuals <- residuals.MArrayLM(object = fit, v)

    ## get average capture effect
    avg_capture_effect <- rowMeans(fit$coefficients)

    ## add average capture effect back into residuals
    corrected_expression <- apply(residuals,2,function(x){x + avg_capture_effect})
    weights <- v$weights
    colnames(weights) <- colnames(corrected_expression)
    rownames(weights) <- rownames(corrected_expression)

    print(length(which(colnames(corrected_expression)!=rownames(meta_data_i))))

    return(list(E=corrected_expression, weights=weights))

}



nzMedianRowDE <- function(gene,o,ce,nzf,m) {

    # gene = gene name to extract row from o, ce, nzf
    # o = observation to model ancestry's effect on
    # ce = corrected pseudobulk expression
    # nzf = nonzero cell fraction
    # m = meta data

    curr_o <- as.vector(as.matrix(o)[gene,])
    curr_ce <- as.vector(as.matrix(ce)[gene,])
    curr_nzf <- as.vector(as.matrix(nzf)[gene,])

    m$pseudobulk <- curr_ce
    m$nzfrac <- curr_nzf
    m$nzmed <- curr_o

    fit <- lm(nzmed ~ capture + age_Scale + lib.size + pseudobulk + nzfrac + YRI_Scale, m)
    coeffs <- setNames(coef(summary(fit))['YRI_Scale',], c('beta','se','t_stat','p'))

    return(coeffs)

}

nzMedianDEmodel <- function(nzm, ce, nzf, m, keep_zeros=TRUE, use_weights=FALSE, base='NI', stim='flu', debug=FALSE) {

    # nzm = nonzero medians df
    # ce = corrected expression
    # nzf = nonzero fraction
    # m = meta data

    print('    Replacing NA nzmedians with 0')
    nzm <- nzm %>%
        dplyr::mutate(dplyr::across(tidyselect::everything(), ~tidyr::replace_na(.x, 0)))
        # dplyr::mutate(dplyr::across(tidyselect::everything(), ~ if_else(is.na(.x),0,.x)))
    print('    Getting counts')
    dge <- DGEList(counts=as.matrix(nzm))
    print('    Getting norm factors')
    dge <- calcNormFactors(dge)

    if (keep_zeros) {
        pseudocount <- 0.5
    } else {
        dge$counts <- as.data.frame(dge$counts) %>%
            dplyr::mutate(dplyr::across(tidyselect::everything(), ~replace(.x, .x==0, NA))) %>%
            drop_na() %>%
            as.matrix()
        pseudocount <- 0
    }

    # voom for correction of nzmedian 'lib sizes' (log CPM 'median counts')
    # and estimation of inverse variance weights of nzmedians
    if (use_weights) {
        v <- voom(dge)
        e <- v$E
    } else {
        e <- log2((t(t(dge$counts+pseudocount)/dge$samples$lib.size))*1e6)
    }

    ebase <- as.data.frame(e) %>% dplyr::select(contains(base))

    if (debug) ebase <- head(ebase)
    basegenes <- rownames(ebase)

    curr_m <- matchMeta(ebase, m)
    curr_o <- matchExpr(ebase, curr_m)
    curr_ce <- matchExpr(as.data.frame(ce), curr_m)
    curr_nzf <- matchExpr(as.data.frame(nzf), curr_m)

    print('Baseline regression...')

    baseres <- as.data.frame(
        setNames(
            pblapply(basegenes,
                     nzMedianRowDE,
                     o=curr_o,
                     ce=curr_ce,
                     nzf=curr_nzf,
                     m=curr_m),
            basegenes
        )
    ) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(hgnc_symbol = rownames(.)) %>%
        as_tibble()


    estim <- as.data.frame(e) %>% dplyr::select(contains(stim))

    if (debug) estim <- head(estim)
    stimgenes <- rownames(estim)

    curr_m <- matchMeta(estim, m)
    curr_o <- matchExpr(estim, curr_m)
    curr_ce <- matchExpr(as.data.frame(ce), curr_m)
    curr_nzf <- matchExpr(as.data.frame(nzf), curr_m)

    print('Stimulation regression...')

    stimres <- as.data.frame(
        setNames(
            pblapply(stimgenes,
                     nzMedianRowDE,
                     o=curr_o,
                     ce=curr_ce,
                     nzf=curr_nzf,
                     m=curr_m),
            stimgenes
        )
    ) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(hgnc_symbol = rownames(.)) %>%
        as_tibble()

    combres <- rbind(
        baseres %>% dplyr::mutate(condition = base),
        stimres %>% dplyr::mutate(condition = stim)
    )


    return(combres)

}

nzFractionRowDE <- function(gene,o,ce,nzm,m,betafit=FALSE) {

    # gene = gene name to extract row from o, ce, nzf
    # o = observation to model ancestry's effect on
    # ce = corrected pseudobulk expression
    # nzm = nonzero log CPM median per cell
    # m = meta data

    curr_o <- as.vector(as.matrix(o)[gene,])
    curr_ce <- as.vector(as.matrix(ce)[gene,])
    curr_nzm <- as.vector(as.matrix(nzm)[gene,])

    m$pseudobulk <- curr_ce
    m$nzmed <- curr_nzm
    m$nzfrac <- curr_o

    if (betafit) {
        fit <- NULL
        try(
            fit <- betareg(nzfrac ~ capture + age_Scale + lib.size + pseudobulk + nzmed + YRI_Scale, m)
        )

        if (is.null(fit)) {
            coeffs <- setNames(rep(NA,4), c('beta','se','z_value','p'))
        } else {
            coeffs <- setNames(coef(summary(fit))$mean['YRI_Scale',], c('beta','se','z_value','p'))
        }

    } else {
        fit <- lm(nzfrac ~ capture + age_Scale + lib.size + pseudobulk + nzmed + YRI_Scale, m)
        coeffs <- setNames(coef(summary(fit))['YRI_Scale',], c('beta','se','t_stat','p'))
    }

    print(paste(gene, 'processed'))

    return(coeffs)

}

nzFractionDEmodel <- function(nzm, ce, nzf, m, keep_zeros=TRUE, use_weights=FALSE, betafit=FALSE, base='NI', stim='flu', debug=FALSE) {

    # nzm = nonzero medians df
    # ce = corrected expression
    # nzf = nonzero fraction
    # m = meta data

    if (betafit) {
        # approximate extreme values for fitting beta distribution
        nzf <- nzf %>%
            dplyr::mutate(dplyr::across(tidyselect::everything(), ~replace(.x, .x==1, 0.9999))) %>%
            dplyr::mutate(dplyr::across(tidyselect::everything(), ~replace(.x, .x==0, 0.0001)))
    }

    nzm <- nzm %>%
        dplyr::mutate(dplyr::across(tidyselect::everything(), ~tidyr::replace_na(.x, 0)))
    dge <- DGEList(counts=nzm)
    dge <- calcNormFactors(dge)

    if (keep_zeros) {
        pseudocount <- 0.5
    } else {
        dge$counts <- as.data.frame(dge$counts) %>%
            dplyr::mutate(dplyr::across(tidyselect::everything(), ~replace(.x, .x==0, NA))) %>%
            drop_na() %>%
            as.matrix()
        pseudocount <- 0
    }

    # voom for correction of nzmedian 'lib sizes' (log CPM 'median counts')
    # and estimation of inverse variance weights of nzmedians
    if (use_weights) {
        v <- voom(dge)
        e <- v$E
    } else {
        e <- log2((t(t(dge$counts+pseudocount)/dge$samples$lib.size))*1e6)
    }

    nzfbase <- as.data.frame(nzf) %>% dplyr::select(contains(base))
    ebase <- as.data.frame(e) %>% dplyr::select(contains(base))

    if (debug) ebase <- head(ebase)
    basegenes <- rownames(ebase) # use nzmedian here in case zeros removed

    curr_m <- matchMeta(ebase, m)
    curr_o <- matchExpr(nzfbase, curr_m)
    curr_ce <- matchExpr(as.data.frame(ce), curr_m)
    curr_nzm <- matchExpr(ebase, curr_m)

    print('Baseline regression...')

    baseres <- as.data.frame(
        setNames(
            pblapply(basegenes,
                     nzFractionRowDE,
                     o=curr_o,
                     ce=curr_ce,
                     nzm=curr_nzm,
                     m=curr_m,
                     betafit=betafit),
            basegenes
        )
    ) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(hgnc_symbol = rownames(.)) %>%
        as_tibble()


    nzfstim <- as.data.frame(nzf) %>% dplyr::select(contains(stim))
    estim <- as.data.frame(e) %>% dplyr::select(contains(stim))

    if (debug) estim <- head(estim)
    stimgenes <- rownames(estim)

    curr_m <- matchMeta(estim, m)
    curr_o <- matchExpr(nzfstim, curr_m)
    curr_ce <- matchExpr(as.data.frame(ce), curr_m)
    curr_nzm <- matchExpr(estim, curr_m)

    print('Stimulation regression...')

    stimres <- as.data.frame(
        setNames(
            pblapply(stimgenes,
                     nzFractionRowDE,
                     o=curr_o,
                     ce=curr_ce,
                     nzm=curr_nzm,
                     m=curr_m,
                     betafit=betafit),
            stimgenes
        )
    ) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(hgnc_symbol = rownames(.)) %>%
        as_tibble()

    combres <- rbind(
        baseres %>% dplyr::mutate(condition = base),
        stimres %>% dplyr::mutate(condition = stim)
    )


    return(combres)

}

# args

parser <- ArgumentParser(description='Decompose single cell RNA-seq data into non-zero fraction and non-zero median expression for each cell type, then do differential expression analysis on ancetry proportion per cell-type per condition.')

parser$add_argument('--betafit',
                    action='store_true',
                    default=FALSE,
                    help="Fit non-zero fraction diff-exp with beta regression model")
parser$add_argument('--keep_zeros',
                    action='store_true',
                    default=FALSE,
                    help="Keep genes with one or more zeros (add pseudocount) for 'non-zero median' expression.")
parser$add_argument('--scCounts_dir',
                    type='character',
                    default='.',
                    help="Path to file with sample metadata for covariates, etc.")
parser$add_argument('--metadata',
                    type='character',
                    default='.',
                    help="Path to file with sample metadata for covariates, etc.")
parser$add_argument('--pseudobulk',
                    type='character',
                    default='.',
                    help="Directory with pre-processsed psuedobulk data")
parser$add_argument('-n', '--name',
                    type='character',
                    default='scDecomp_DE',
                    help="Basename for output files")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

betafit <- opt$betafit
keep_zeros <- opt$keep_zeros
sccount_dir <- opt$scCounts_dir
meta_data_file <- opt$metadata
pseudobulk_dir <- opt$pseudobulk
outbase <- opt$name
outdir <- opt$outdir

if (keep_zeros) {
  outbase <- paste0(outbase, '_nzmedPseudocount')
}
if (betafit) {
  outbase <- paste0(outbae, '_nzfracBetafit')
}


## read in meta data from individuals
meta_data <- read.table(meta_data_file, header = TRUE, sep = ",")
meta_data$capture <- as.factor(meta_data$capture)
meta_data$indiv_ID <- as.factor(meta_data$indiv_ID)
meta_data$infection_status = factor(meta_data$infection_status,
                                    levels=c("NI","flu"))

# initalize df for rbinding results
aggres <- data.frame()

for (pb in list.files(pseudobulk_dir)) {

    meta_data_i <- meta_data

    pbsplit <- str_split(pb, '_', simplify=TRUE)
    cell_type_i <- paste(pbsplit[1:ncol(pbsplit)-1], collapse='_')

    if (str_detect(cell_type_i,'_combined')) cell_type_i <- gsub('_combined','',cell_type_i)

    if (str_detect(cell_type_i,'all')) {
        cell_type_i <- 'PBMCs'
        rawdata_fname <- list.files(sccount_dir, 'merged*', full.names=TRUE)
        if (length(rawdata_fname)>1) stop('Multiple merged scCount files found')
        print(rawdata_fname[1])
    } else {
        rawdata_fname <- file.path(sccount_dir, paste0(cell_type_i, '_cluster_singlets.rds'))
        print(rawdata_fname)
    }

    # skip PBMC merged for now
    if (cell_type_i=='PBMCs') {
      print('Skipping PBMCs')
      next
    }

    print(paste0('Processing ', cell_type_i, '...'))

    reads <- readRDS(file.path(pseudobulk_dir, pb))

    meta_data_i <- matchMeta(reads, meta_data_i)
    reads <- matchExpr(reads,meta_data_i)

    ## remove lowly-expressed genes
    dge <- DGEList(counts = reads)
    dge <- calcNormFactors(dge)
    design = model.matrix(~ 0 + capture, data = meta_data_i)

    ## remove columns that are all 0s
    design <- design[, colSums(design != 0) > 0]
    # voom for removal of lowly expressed genes in nzfrac/nzmedian
    v <- voom(dge, design, plot = FALSE)

    # get cell counts, non-zero fraction, and non-zero median
    if (cell_type_i=='PBMCs') {
      # check that processing is the same for merged
    } else {
      ## read in single cell data (raw UMI counts stored in RNA slot) for each
      print('    Reading in SC data...')
      dat <- readRDS(rawdata_fname)
      ## get raw UMI counts
      raw_data_sparse <- GetAssayData(dat, assay = "RNA", slot = "counts")
      meta_data_sc <- dat@meta.data
      sample_colname <- "sample_condition"

      IDs <- as.data.frame(meta_data_sc)[, sample_colname]
      unique_ID_list <- as.list(unique(IDs))

      cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(raw_data_sparse[,IDs == x, drop = FALSE])})
      names(cellcount) <- unique_ID_list
      cellcount_df <- listToTibble(cellcount, cnames=c('infection_ID','cellcount'))

      print('    Calculating non-zero fractions...')
      # get non-zero fractions
      ## calculate fraction of cells expressing each gene across all cells identified in each individual,cdt pair
      cellfrac <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){unname(rowSums(raw_data_sparse[,IDs == x, drop = FALSE] > 0))/ncol(raw_data_sparse[,IDs == x, drop = FALSE])}))
      colnames(cellfrac) <- unique_ID_list
      rownames(cellfrac) <- rownames(x = dat)

      print('    Calculating non-zero medians...')
      # get non-zero medians
      ## calculate mean expression among cells with nonzero expression
      ## for each gene across cells identified in each individual,cdt pair
      nzmedian <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){floor(unname(apply(raw_data_sparse[,IDs == x, drop = FALSE], 1, nonzeroMedian)))}))
      colnames(nzmedian) <- unique_ID_list
      rownames(nzmedian) <- rownames(x = dat)

      # remove same low counts as pseudobulk
      cellfrac <- removeLowCounts(cell_type_i, reads, v, decomp=cellfrac)
      nzmedian <- removeLowCounts(cell_type_i, reads, v, decomp=nzmedian)
    }

    # remove low counts pre-library size calculation
    reads <- removeLowCounts(cell_type_i, reads, v)
    dge <- DGEList(counts = reads)

    # add libsizes to metadata
    libsizes <- dge$samples
    libsizes$infection_ID <- rownames(libsizes)
    libsizes <- libsizes %>% select(infection_ID, lib.size)
    meta_data_i <- merge(meta_data_i, libsizes, by='infection_ID')
    rownames(meta_data_i) <- meta_data_i$infection_ID

    ## subset correct meta_data bimodal proportion
    counts <- subset(meta_data_i, select = c(paste0(cell_type_i,"_geneProp")))
    colnames(counts)[1] <- "bimodal_prop"
    meta_data_i <- meta_data_i %>% select(-ends_with('_geneProp'))
    meta_data_i <- cbind(meta_data_i, counts)

    # add cell counts to metadata
    meta_data_i <- merge(meta_data_i, cellcount_df, by='infection_ID')
    rownames(meta_data_i) <- meta_data_i$infection_ID
    meta_data_i <- matchMeta(reads,meta_data_i)

    reads <- matchExpr(reads,meta_data_i)

    print('    Modeling capture effect...')
    # capture correct expression
    voomreads <- voomCaptureCorrect(reads, meta_data_i)
    corrected_expression <- voomreads$E
    weights <- voomreads$weights

    ## write capture-corrected expression and weights for later modeling
    write.table(corrected_expression,
                file.path(outdir,
                          paste0(cell_type_i,"_corrected_expression.txt")),
                          quote = FALSE,
                          sep = ",")
    write.table(weights,
                file.path(outdir,
                          paste0(cell_type_i,"_weights.txt")),
                quote = FALSE,
                sep = ",")

    cellfrac <- cellfrac %>%
        dplyr::mutate(hgnc_symbol = rownames(.)) %>%
        as_tibble() %>%
        dplyr::filter(!is.na(hgnc_symbol)) %>%
        as.data.frame()
    rownames(cellfrac) <- cellfrac %>% pull(hgnc_symbol)
    cellfrac <- cellfrac %>% dplyr::select(-hgnc_symbol)

    nzmedian <- nzmedian %>%
        dplyr::mutate(hgnc_symbol = rownames(.)) %>%
        as_tibble() %>%
        dplyr::filter(!is.na(hgnc_symbol)) %>%
        as.data.frame()
    rownames(nzmedian) <- nzmedian %>% pull(hgnc_symbol)
    nzmedian <- nzmedian %>% dplyr::select(-hgnc_symbol)

    # save nzfracs and nzmedians
    write_tsv(
      cellfrac %>%
          dplyr::mutate(hgnc_symbol = rownames(.)) %>%
          as_tibble() %>%
          dplyr::select(hgnc_symbol, everything()),
      file.path(outdir,
                paste0(cell_type_i,"_nonzeroFraction.txt"))
    )
    if (keep_zeros) {
      nzmedsuffix <- '_nonzeroMedianPseudocounts.txt'
    } else {
      nzmedsuffix <- '_nonzeroMedian.txt'
    }
    write_tsv(
        nzmedian %>%
            dplyr::mutate(hgnc_symbol = rownames(.)) %>%
            as_tibble() %>%
            dplyr::select(hgnc_symbol, everything()),
        file.path(outdir,
                  paste0(cell_type_i,nzmedsuffix))
    )

    # rename celltype for downstream processing
    cell_type_i <- gsub('_','',cell_type_i)

    print('    Modeling non-zero median DE')
    # process non-zero medians
    nzmedres <- nzMedianDEmodel(nzmedian,
                            corrected_expression,
                            cellfrac,
                            meta_data_i,
                            keep_zeros=keep_zeros,
                            use_weights=FALSE,
                            base='NI',
                            stim='flu',
                            debug=FALSE) %>%
          dplyr::mutate(celltype = cell_type_i,
                 datatype = 'nzmedian')

     print('    Modeling non-zero fraction DE')
     nzfracres <- nzFractionDEmodel(nzmedian,
                             corrected_expression,
                             cellfrac,
                             meta_data_i,
                             keep_zeros=keep_zeros,
                             use_weights=FALSE,
                             betafit=betafit,
                             base='NI',
                             stim='flu',
                             debug=FALSE) %>%
          dplyr::mutate(celltype = cell_type_i,
                 datatype = 'nzfraction')

     aggres <- rbind(aggres, nzmedres)
     aggres <- rbind(aggres, nzfracres)
}

## adding genes call and converting to tibble pre-base/stim rbind
## should fix auto renaming rows (dataframes can't have duplicate rownames)
# aggres <- aggres %>%
#     dplyr::mutate(hgnc_symbol=rownames(.),
#            hgnc_symbol=if_else(grepl('*1$',hgnc_symbol) & condition=='flu',
#                                gsub('1$','',hgnc_symbol),
#                                hgnc_symbol)) %>%
#     as_tibble()

aggres %>%
    write_tsv(file.path(outdir,paste0(outbase,'_allStats.txt.gz')))

aggres %>%
    dplyr::select(-c(se,t_stat)) %>%
    pivot_wider(names_from=c(condition,celltype,datatype),
                values_from=c(beta,p),
                names_glue="{.value}_{condition}_{celltype}_{datatype}") %>%
    drop_na() %>%
    write_tsv(file.path(outdir,paste0(outbase,'_forOverlap.txt.gz')))
