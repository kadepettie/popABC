#!/usr/bin/R

if(!(require(argparse) )) install.packages("argparse")
if(!(require(tidyverse) )) install.packages("tidyverse")
if(!(require(ggplot2)   )) install.packages("ggplot2")
if(!(require(caret)   )) install.packages("caret")
if(!(require(lme4)   )) install.packages("lme4")
if(!(require(tidymodels)   )) install.packages("tidymodels")
if(!(require(relaimpo)   )) install.packages("relaimpo")
if(!(require(docstring)   )) install.packages("docstring")

if(!(require(plyranges))) {
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("plyranges")
}

library(argparse)
library(tidyverse)
library(ggplot2)
library(plyranges)
library(caret)
library(lme4)
library(tidymodels)
library(relaimpo)
library(docstring)


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

calc_relimp_mm <- function(model,type='lmg') {

    # adapted from https://gist.github.com/BERENZ/e9b581a4b7160357934e

    if (!isLMM(model) & !isGLMM(model)) {
        stop('Currently supports only lmer/glmer objects', call. = FALSE)
    }
    X <- getME(model,'X')
    X <- X[,-1]
    Y <- getME(model,'y')
    s_resid <- sigma(model)
    s_effect <- getME(model,'theta')*s_resid
    s2 <- sum(s_resid^2,s_effect^2)
    V <- Diagonal(x = s2,n=nrow(X))
    YX <- cbind(Y,X)
    cov_XY <- solve( t(YX) %*% solve(V) %*% as.matrix(YX))
    colnames(cov_XY) <- rownames(cov_XY) <- colnames(YX)
    importances <- calc.relimp(as.matrix(cov_XY),rela=FALSE,type=type)
    return(importances)
}

modelPrep <- function(df,
                      scale_vars=NULL,
                      top_enh=NULL,
                      resp_var='beta',
                      chip_var=NULL,
                      atac_var=NULL,
                      int_var='mean_ABC',
                      ...) {

    #' Prep dataframe for modeling
    #'
    #' Subsets a gene expression-enhancer-matched data set
    #' for modeling with templated interaction formula
    #'
    #' @param df Dataframe with model variables to be subset
    #' @param scale_vars Variables to center and scale before subsetting data for modeling; e.g., c("atac_lfc","atac_nltp","chip_lfc","chip_nltp",'mean_ABC')
    #' @param top_enh Only use top enhancer defined by this metric; e.g., 'abc', 'chip', 'atac', 'self', NULL (default)
    #' @param resp_var Response variable in formula; e.g., 'de_sig', 'snltp_DE', 'beta' (beta (default) = beta among de_sig=1 genes)
    #' @param chip_var ChIP component-related (or arbitrary var 1) variable to include in template
    #' @param atac_var ATAC component-related (or arbitrary var 2) variable to include in template
    #' @param int_var Include interaction effect, but no main effect for this variable; e.g., 'mean_ABC' (default)
    #' @param ... For accepting additional args from wrapper function
    #' @return Named list of `df` (subset df , `form`, and `meth`

    if (is.null(chip_var) | is.null(atac_var)) {

        if (resp_var=='de_sig') {
            atac_var <- 'atac_nltp'
            chip_var <- 'chip_nltp'
        } else if (resp_var=='snltp_DE') {
            atac_var <- 'snltp_atac'
            chip_var <- 'snltp_chip'
        } else if (resp_var=='beta') {
            df <- df %>% dplyr::filter(de_sig==1)
            atac_var <- 'atac_lfc'
            chip_var <- 'chip_lfc'
        } else {
            stop(paste0("Response variable not recognized: ", resp_var))
        }

    }

    if (!is.null(scale_vars)) {
        # possible this should be done separately on promoters and enhancers (i.e., after below steps)
        # but for now scaling/centering only sometimes suggested during lmer() runs
        df <- df %>%
            mutate(across(any_of(scale_vars), ~ as.numeric(scale(.x))))

    }

    mod_base <- paste0(resp_var,' ~ ',chip_var,'*',atac_var)
    mixmod_base <- paste0(resp_var,' ~ ',chip_var,'*',atac_var,' + ',int_var,':(',chip_var,'*',atac_var,')')
    mod_meth <- 'lmer'
    if (resp_var=='de_sig') mod_meth <- 'glmer'

    if (is.null(top_enh)) {

        mod_form <- paste0(mixmod_base," + (1 | name) + (1 | hgnc_symbol)")

    } else {

        mod_form <- paste0(mixmod_base," + (1 | name)")

        if (top_enh=='abc') {

            df <- df %>% dplyr::filter(isSelfPromoter==0) %>%
                dplyr::group_by(hgnc_symbol) %>%
                arrange(desc(mean_ABC)) %>%
                dplyr::slice(1) %>%
                ungroup()

        } else if (top_enh=='chip') {

            df <- df %>% dplyr::filter(isSelfPromoter==0) %>%
                dplyr::group_by(hgnc_symbol) %>%
                arrange(desc(chip_nltp)) %>%
                dplyr::slice(1) %>%
                ungroup()

        } else if (top_enh=='atac') {

            df <- df %>% dplyr::filter(isSelfPromoter==0) %>%
                dplyr::group_by(hgnc_symbol) %>%
                arrange(desc(atac_nltp)) %>%
                dplyr::slice(1) %>%
                ungroup()

        } else if (top_enh=='self') {

            df <- df %>% dplyr::filter(isSelfPromoter==1)

            mod_form <- mod_base
            mod_meth <- 'lm'

        } else {

            stop(paste0("Top enh type not recognized: ", top_enh))

        }

    }

    return(list(df=df, form=mod_form, meth=mod_meth))

}

makeModel <- function(df,
                        stratum=NULL,
                        modform = 'beta ~ .',
                        method='lm',
                        center_scale = NULL,
                        random_train = FALSE,
                        random_control = 'mean_ABC',
                        ...
                       ) {

    #' Train model for testing
    #'
    #' Splits data into training (3/4) and testing (1/4)
    #' and trains model using the provided formula and method.
    #'
    #' @param df Dataframe with model variables included in `modform`
    #' @param stratum Variable to ensure equal representation of in train and test data
    #' @param modform Formula specifying model as character string.
    #' @param method Model method to use; e.g., 'lm' (default), 'lmer', 'glmer'
    #' @param center_scale Character vector of variables to center and scale before modeling
    #' @param random_train Randomize training data response vector before training model
    #' @param random_control Variable to control for when random_train=TRUE
    #' @param ... For accepting additional args from wrapper function
    #' @return Named list of `modelfit` (`lm` or `merMod` object with fit model) and `testdata` (hold out data to test model on)

    if (!is.null(center_scale)) {
        print('Centering and scaling...')
        df <- df %>%
            mutate(across(all_of(center_scale), ~ as.numeric(scale(.x))))
    }

    dsplit <- initial_split(df, prop = 3/4, strata=stratum)
    dtrain <- training(dsplit)
    dtest <- testing(dsplit)

    if (random_train) {

        rdtrain <- permute_df(dtrain,
                             c('snltp_DE','hgnc_symbol','de_sig','beta'),
                             strata_var=random_control,
                             breaks=10)

    }

    print(paste0('Running model: ', modform))

    if (method=='lm') {

        dfit <- lm(eval(parse(text=modform)),
                   data=dtrain)

        if (random_train) {
            rdfit <- lm(eval(parse(text=modform)),
                        data=rdtrain)
        }

    } else if (method=='lmer') {

        dfit <- lmer(eval(parse(text=modform)),
                     data=dtrain,
                     REML=TRUE,
                     control=lmerControl(boundary.tol = 1e-7))

        if (random_train) {
            rdfit <- lmer(eval(parse(text=modform)),
                          data=rdtrain,
                          REML=TRUE,
                          control=lmerControl(boundary.tol = 1e-7))
        }

    } else if (method=='glmer') {

        dfit <- glmer(eval(parse(text=modform)),
                      data=dtrain,
                      family=binomial,
                      control=glmerControl(boundary.tol = 1e-7))

        if (random_train) {
            rdfit <- glmer(eval(parse(text=modform)),
                           data=rdtrain,
                           family=binomial,
                           control=glmerControl(boundary.tol = 1e-7))
        }

    } else {

        stop(paste0("Method ", method, " not supported"))

    }

    print(summary(dfit))

    if (random_train) {
        return(list(modelfit=dfit, rand_modelfit=rdfit, testdata=dtest))
    }
    return(list(modelfit=dfit, testdata=dtest))

}

permute_df <- function(df, perm_vars, strata_var=NULL, breaks=10) {

    # perm_vars: vars to permute (e.g., response var)
    # strata_var: restrict permutation to within `breaks` percentile
    #     increments of this var

    #' Permute selected dataframe columns
    #'
    #' Permutes selected dataframe columns optionally restricted
    #' to within strata of a variable to preserve some association with.
    #'
    #' @param df Dataframe with variables to permuate relative to others
    #' @param perm_vars Variables to permute
    #' @param strata_var Variable to stratify permutation by
    #' @param breaks Number of percentile invervals to break `strata_var` into

    library(rsample)

    df <- df %>%
        mutate(sbins := make_strata(.data[[strata_var]], breaks=breaks))

    hdf <- df %>%
        dplyr::select(-any_of(perm_vars)) %>%
        arrange(sbins) %>%
        dplyr::select(-sbins)

    pvdf <- df %>%
        dplyr::select(any_of(perm_vars), sbins) %>%
        arrange(sbins) %>%
        group_by(sbins) %>%
        dplyr::slice_sample(prop=1) %>%
        ungroup()

    return(cbind(pvdf, hdf))

}

RSQ <- function(y, x, zero_int=TRUE) {
        if (zero_int) {
            lms <- summary(lm(y ~ 0 + x))
            pvidx <- 1
        } else {
            lms <- summary(lm(y ~ x))
            pvidx <- 2
        }
        rsq <- lms$r.squared
        rsqadj <- lms$adj.r.squared
        pv <- lms$coefficients[pvidx,'Pr(>|t|)']
        return(list(Rsq=rsq, Rsq.adj=rsqadj, P=pv))
    }

trainTestModel <- function(df,
                           top_enh=NULL,
                           resp_var='beta',
                           celltype='LCL',
                           condition='H20',
                           debug=FALSE, ...) {

    #' Prep, train, and test model
    #'
    #' Wrapper for `modelPrep` and `makeModel` that also
    #' tests the model on a holdout set and compares performance
    #' to model trained on permuted response variable
    #'
    #' @param df Dataframe to pass to `modelPrep`
    #' @param ... Other args to pass to `modelPrep` and `makeModel`
    #' \subsection{Args to pass}{
    #'     @param df Dataframe with model variables to be subset
    #'     @param scale_vars Variables to center and scale before subsetting data for modeling; e.g., c("atac_lfc","atac_nltp","chip_lfc","chip_nltp",'mean_ABC')
    #'     @param top_enh Only use top enhancer defined by this metric; e.g., 'abc', 'chip', 'atac', 'self', NULL (default)
    #'     @param resp_var Response variable in formula; e.g., 'de_sig', 'snltp_DE', 'beta' (beta (default) = beta among de_sig=1 genes)
    #'     @param chip_var ChIP component-related (or arbitrary var 1) variable to include in template
    #'     @param atac_var ATAC component-related (or arbitrary var 2) variable to include in template
    #'     @param int_var Include interaction effect, but no main effect for this variable; e.g., 'mean_ABC' (default)
    #'
    #'     @param df Dataframe with model variables included in `modform`
    #'     @param modform Formula specifying model as character string.
    #'     @param center_scale Character vector of variables to center and scale before modeling
    #'     @return Named list of `modelfit` (`lm` or `merMod` object with fit model) and `testdata` (hold out data to test model on)
    #'     @param ... For accepting additional args from wrapper function

    #' }

    df_form <- modelPrep(df,
                         int_var='mean_ABC',
                         top_enh=top_enh,
                         resp_var=resp_var,
                         ...)

    fit_test <- makeModel(df_form$df,
                          modform = df_form$form,
                          method = df_form$meth,
                          random_train = TRUE,
                          ...)

    dfit <- fit_test$modelfit
    rdfit <- fit_test$rand_modelfit
    dtest <- fit_test$testdata

    dpred <- predict(dfit, dtest, type='response', allow.new.levels=TRUE)
    rdpred <- predict(rdfit, dtest, type='response', allow.new.levels=TRUE)

    predobs <- cbind(dtest, predicted = dpred)
    rpredobs <- cbind(dtest, predicted = rdpred)

    predobs_pdf <- rbind(
        predobs %>%
            separate(name, c('class',NA), sep='\\|', remove=FALSE, convert=TRUE) %>%
            dplyr::rename(observed := {{resp_var}}) %>%
            dplyr::select(class, observed, predicted) %>%
            mutate(mod_type='Model'),
        rpredobs %>%
            separate(name, c('class',NA), sep='\\|', remove=FALSE, convert=TRUE) %>%
            dplyr::rename(observed := {{resp_var}}) %>%
            dplyr::select(class, observed, predicted) %>%
            mutate(mod_type='Random model')
    ) %>%
        mutate(class = if_else(class=='promoter', class, 'non-promoter'))

    if (df_form$meth=='glmer') {

        p1 <- plotRoc(predobs_pdf %>%
                         group_by(mod_type) %>%
                         roc_curve(truth=observed, predicted),
                      facets='mod_type')

        perf_df <- predobs_pdf %>%
            group_by(mod_type) %>%
            roc_auc(truth=observed, predicted)

        # use add1() internal to makeModel() (not implemented yet)
        # to assess variable importance
        # for now use placeholder for processing aggregate results
        relimp_df <- as_tibble(coef(summary(dfit)), rownames='predictor') %>%
            mutate(importance=NA) %>%
            dplyr::select(predictor, importance)

    } else {

         p1 <- plotCorrelation(predobs_pdf %>%
                                 dplyr::arrange(desc(class)),
                            'observed',
                            'predicted',
                            groupvar = "class",
                            colorvals = NULL,
                            addline = c(0,1),
                            legendpos = c(.82,.9),
                            corrpos = c(0,0.03),
                            corrcolor = "black",
                            w = 14, h = 7,
                            logscale = FALSE,
                            addaxes=TRUE,
                            facets = 'mod_type',
                            frows = 1,
                            debug=FALSE)

        trsq <- RSQ(predobs[[resp_var]], predobs$predicted)
        rtrsq <- RSQ(rpredobs[[resp_var]], rpredobs$predicted)
        trsq$mod_type <- 'Model'
        rtrsq$mod_type <- 'Random model'
        perf_df <- rbind(as_tibble(trsq), as_tibble(rtrsq))

        if (df_form$meth=='lm') {
            relimp_dfit <- calc.relimp(dfit,rela=FALSE,type='lmg')
        } else {
            relimp_dfit <- calc_relimp_mm(dfit)
        }
        print(relimp_dfit)
        relimp_df <- as_tibble(relimp_dfit$lmg, rownames='predictor') %>%
            dplyr::rename(importance=value) %>%
            arrange(desc(importance))

    }

    perf_df <- perf_df %>%
        mutate(enh_type=top_enh,
               condition=condition,
               celltype=celltype,
               response=resp_var,
               formula=df_form$form,
               method=df_form$method)

    relimp_df <- relimp_df %>%
        mutate(enh_type = top_enh,
               condition = condition,
               celltype = celltype,
               response = resp_var,
               formula = df_form$form,
               method = df_form$meth)

    fixef <- as_tibble(coef(summary(dfit)), rownames='predictor') %>%
        dplyr::rename(estimate = Estimate,
                      stderr = `Std. Error`) %>%
    mutate(conf_lower = estimate - stderr,
           conf_upper = estimate + stderr,
           enh_type = top_enh,
           condition = condition,
           celltype = celltype,
           response = resp_var,
           formula = df_form$form,
           method = df_form$meth)

    if (df_form$meth=='lm') {
        ranef <- NULL
    } else {
        ranef <- as.data.frame(VarCorr(dfit)) %>%
            mutate(enh_type = top_enh,
                   condition = condition,
                   celltype = celltype,
                   response = resp_var,
                   formula = df_form$form,
                   method = df_form$meth)
    }

    return(list(relimp=relimp_df, fixef=fixef, ranef=ranef, performance=perf_df, pplot=p1))

}

writeTrainTestResults <- function(res,
                                  outdir='.',
                                  outprefix='modelRes') {

    write_tsv(res$relimp, file.path(outdir, paste0(outprefix, '.relativeImportance.txt')))
    write_tsv(res$fixef, file.path(outdir, paste0(outprefix, '.fixedEffects.txt')))
    if (!is.null(res$ranef)) write_tsv(res$ranef, file.path(outdir, paste0(outprefix, '.randomEffects.txt')))
    write_tsv(res$performance, file.path(outdir, paste0(outprefix, '.performance.txt')))

    res$pplot
    ggsave(file.path(outdir,paste0(outprefix,'.performance.png')),
           width=14,
           height=7)

}


parser <- ArgumentParser(description='Aggregate differential test results of ABC Score components for QC/diagnostics and spot checking.')

parser$add_argument('--response_var',
                    type='character',
                    default='beta',
                    help="Type of DE to model as response variable")
parser$add_argument('--top_enh',
                    type='character',
                    default=NULL,
                    help="Select only one enhancer per gene by this metric to model DE with.")
parser$add_argument('--celltype',
                    type='character',
                    default='LCL',
                    help="Subset to this cell type ('celltype' col) for modeling.")
parser$add_argument('--condition',
                    type='character',
                    default=NULL,
                    help="Subset to this condition ('condition' col) for modeling.")
parser$add_argument('--facetvar',
                    type='character',
                    default=NULL,
                    help="Set to 'datatype' if DE has multiple datatypes to make separate plots for")
parser$add_argument('--pval_prefix',
                    type='character',
                    default='pval',
                    help='prefix of columns with p-values')
parser$add_argument('--ES_prefix',
                    type='character',
                    default='beta',
                    help='prefix of columns with effect sizes')
parser$add_argument('-n', '--name',
                    type='character',
                    default='diffABC_DE_overlap',
                    help="Basename for output files")
parser$add_argument('--plotdir',
                    type='character',
                    default=NULL, # /home/kpettie/code/github/plotting
                    help="Directory with plotting function in 'plotting.R'. Required to make plots")
parser$add_argument('--DE',
                    type='character',
                    default='./LeaEtAl2021_popDE.txt.gz',
                    help="File with differential expression data to predict.")
parser$add_argument('--diff_all',
                    type='character',
                    default='./allZerosFilt.meanQN.16.AFR_EUR.diff.activity.allComponents.txt.gz',
                    help="File with all diff-ABC score components to use as predictor.")
parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

mod_resp <- opt$response_var
etype <- opt$top_enh
cdt <- opt$condition
ct <- opt$celltype
pv_pre <- opt$pval_prefix
es_pre <- opt$ES_prefix
facetvar <- opt$facetvar
if (!is.null(facetvar)) facetvar <- str_split(facetvar, ',')[[1]]

plotdir <- opt$plotdir

de_fname <- opt$DE
d_all_fname <- opt$diff_all
outdir <- opt$outdir
outbase <- opt$name

negate_beta <- FALSE
if (grepl('RandolphEtAl2021_popDE', de_fname)) negate_beta <- TRUE

source(file.path(plotdir,"plotting.R"))

########### PREP DATA #############

dabc <- prepABCforGSEA(read_tsv(d_all_fname), allscores=TRUE) %>%
    dplyr::select(-c(group1,group2,n1,n2)) %>%
    mutate(neglog10p_ABC = -log10(p)) %>%
    dplyr::rename(log2FoldChange_ABC = log2FoldChange,
                  p_ABC = p)

meanabc <- dabc %>%
    dplyr::filter(score_type=='ABC') %>%
    mutate(mean_ABC = (Score_AFR + Score_EUR)/2) %>%
    dplyr::select(ensembl_gene_id,
                  entrez_id,
                  hgnc_symbol,
                  seqnames,
                  start,
                  end,
                  class,
                  name,
                  distance,
                  isSelfPromoter,
                  mean_ABC)

chipatac <- dabc %>%
    dplyr::filter(score_type %in% c('chip','atac')) %>%
    dplyr::rename(lfc = log2FoldChange_ABC,
                  p = p_ABC,
                  nltp = neglog10p_ABC) %>%
    dplyr::select(ensembl_gene_id,
                  entrez_id,
                  hgnc_symbol,
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
                  score_type) %>%
    pivot_wider(names_from=score_type,
                names_glue="{score_type}_{.value}",
                values_from=c(lfc, p, nltp))

wabc <- merge(meanabc, chipatac)

de <- read_tsv(de_fname)
if (pv_pre=='p') {
    de <- de %>%
        rename_with(.,
                    ~ gsub("p_","pval_",.x),
                    .cols=tidyselect::starts_with('p_'))
    pv_pre <- 'pval'
}

abcde_all <- merge(wabc, de) %>%
    pivot_longer(cols=matches(paste0(es_pre,'_|',pv_pre,'_')),
                 names_to=c(".value","condition","celltype","datatype"),
                 names_pattern="(.*)_(.*)_(.*)_(.*)") %>%
    dplyr::rename(beta := {{es_pre}},
                  p_DE := {{pv_pre}})

if (negate_beta) abcde_all <- abcde_all %>% mutate(beta = -beta)

abcde <- abcde_all %>%
    dplyr::filter(condition==cdt,
                  celltype==ct) %>%
    mutate(snltp_DE = -log10(p_DE)*sign(beta),
            snltp_chip = chip_nltp*sign(chip_lfc),
            snltp_atac = atac_nltp*sign(atac_lfc))

mdf <- abcde %>%
    mutate(de_sig = if_else(p_DE < 0.05, 1, 0),
           isSelfPromoter = if_else(isSelfPromoter, 1, 0),
           de_sig = as_factor(de_sig),
           isSelfPromoter = as_factor(isSelfPromoter)) %>%
    dplyr::select(snltp_DE,
                  mean_ABC,
                  snltp_atac,
                  snltp_chip,
                  atac_nltp,
                  chip_nltp,
                  hgnc_symbol,
                  name,
                  isSelfPromoter,
                  de_sig,
                  beta,
                  atac_lfc,
                  chip_lfc) %>%
    mutate(hgnc_symbol = as_factor(hgnc_symbol),
           name = as_factor(name))

################ RUN MODEL #################

set.seed(333)

ttstratum <- 'de_sig'
if (mod_resp=='beta') ttstratum <- NULL
ttout <- trainTestModel(mdf,
                        top_enh=etype,
                        resp_var=mod_resp,
                        stratum=ttstratum,
                        scale_vars=NULL,
                        center_scale=NULL,
                        random_control='mean_ABC',
                        celltype=ct,
                        condition=cdt)

writeTrainTestResults(ttout,outdir=outdir,outprefix=outbase)
