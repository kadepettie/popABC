#!/usr/bin/R

if(!(require(ggplot2)   )) install.packages("ggplot2")
if(!(require(tidyverse) )) install.packages("tidyverse")

library(ggplot2)
library(tidyverse)

# traceback to written code
options(error=function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
# fix plot plotting dimensions for scaling
# options(repr.plot.width = 18, repr.plot.height = 12, repr.plot.res = 100)

popPlotColors <- function(pops=c('CEU','FIN','IBS','TSI','ASW','ESN','GWD','LWK','YRI','CHB'),
                          reps=TRUE,
                          rep_colors=NULL) {

    popco_orig <- tribble(
                    ~population,  ~hex,
                    'CEU', '#EAD78D',
                    'FIN', '#E1AB31',
                    'IBS', '#D67F30',
                    'TSI', '#9A4C28',
                    'ASW', '#BEE1F8',
                    'ESN', '#7FCDF3',
                    'GWD', '#6CA7DA',
                    'LWK', '#4263AE',
                    'YRI', '#2A3D88',
                    'CHB', '#7AC09F'
                ) %>%
                arrange(desc(row_number()))

    if (reps) {
        pops2 <- c(pops,pops) %>% sort()
        reps <- rep(c('rep1','rep2'),length(pops2)/2)
        samps <- tibble(population=pops2,rep=reps) %>%
            mutate(sample = paste(population,rep,sep='_')) %>%
            dplyr::select(-rep)

        popco <- merge(popco_orig,samps)

        if (!is.null(rep_colors)) {
            colnames(rep_colors) <- c('sample','alt_hex')
            popco <- merge(popco,rep_colors,all.x=TRUE) %>%
                mutate(hex = case_when(!is.na(alt_hex) ~ alt_hex,
                                       TRUE ~ hex)) %>%
                dplyr::select(-alt_hex)
        }

        popco <- popco %>%
            arrange(match(population,popco_orig$population),sample)
    }

    popco <- popco %>%
        dplyr::filter(population %in% pops) %>%
        separate(sample,c(NA,'replicate'),remove=FALSE)

    return(popco)

}

addAncestry <- function(df) {

    afrpops <- c('ESN','GWD','LWK','YRI')
    eurpops <- c('CEU','FIN','IBS','TSI')

    df <- df %>%
        mutate(ancestry = case_when(population %in% afrpops ~ 'AFR',
                                    population %in% eurpops ~ 'EUR'))

    df <- df %>%
            mutate(population=fct_relevel(population,
                                          'CEU',
                                          'FIN',
                                          'IBS',
                                          'TSI',
                                          'ESN',
                                          'GWD',
                                          'LWK',
                                          'YRI'))

    return(df)
}

plotPrepFromMerge <- function(df, colprefix='activity_base', filt_zeros=TRUE, add_ancestry=FALSE) {

    # filt_zeros should be true in most cases to avoid memory overload in
    # pivot_longer when there are many zero elements
    # filters elements with zero in any sample
    # purely for visualization

    if (!('class' %in% colnames(df))) {

        if (!('name' %in% colnames(df))) {
            stop('No promoter column info found...')
        }

        df <- df %>%
            separate(name,c('class',NA),sep='\\|',remove=FALSE)

    }

    if (filt_zeros) df <- df %>% dplyr::filter(if_all(starts_with(colprefix), ~ . > 0))

    df <- df %>%
        mutate(class=if_else(class=='promoter','promoter','non-promoter')) %>%
        pivot_longer(cols=dplyr::starts_with(colprefix),
                     names_to=c("data_type","sample"),
                     names_pattern="(.*)\\.(.*)") %>%
        separate(sample, c('pop','replicate'), sep='_', remove=FALSE)

    if (length(colprefix)==1) df <- df %>% dplyr::rename({{colprefix}} := value)

    pclrs <- popPlotColors(pops=unique(df$pop))
    r1colors <- pclrs %>% dplyr::filter(replicate=='rep1')

    df <- df %>%
        mutate(sample = as_factor(sample),
               sample = fct_relevel(sample,pclrs$sample),
               pop = as_factor(pop),
               pop = fct_relevel(pop, r1colors$population),
               replicate = as_factor(replicate),
               replicate = fct_relevel(replicate, c('rep2','rep1')))

    if (add_ancestry) {
        df <- df %>%
            dplyr::rename(population=pop) %>%
            addAncestry(.) %>%
            dplyr::rename(pop=population)
    }

    return(df)

}

groupedDistPlot <- function(df,
                      xvar='count',
                      groupvar='ancestry',
                      density=FALSE,
                      position='stack',
                      bwidth=NULL,
                      nbins=150,
                      logscale=FALSE,
                      colorvals=NULL, # c("#2E3191","#F05A28")
                      legendpos=c(.8,.8),
                      w=7,h=7) {

    options(repr.plot.width = w, repr.plot.height = h, repr.plot.res = 200)

    distp <- df %>%
        ggplot(aes(x=.data[[xvar]],color=.data[[groupvar]],fill=.data[[groupvar]]))

    if (density) {
        distp <- distp + geom_density(alpha=0.25)
    } else {
        if (!is.null(bwidth)) {
            binseq <- seq(min(df[[xvar]])-1,
                          max(df[[xvar]])+1,
                          by=bwidth)
            distp <- distp + geom_histogram(alpha=0.25,breaks=binseq,position=position)
        } else {
            distp <- distp + geom_histogram(alpha=0.25,bins=nbins,position=position)
        }
    }

    if (logscale) {
        distp <- distp + scale_x_log10()
    }

    if (is.null(colorvals)) {
        distp <- distp +
            scale_color_viridis_d(option = "D") +
            scale_fill_viridis_d(option = "D")
    } else {
        distp <- distp +
            scale_fill_manual(values = colorvals) +
            scale_color_manual(values = colorvals)
    }

    distp <- distp + theme_classic(20) + theme(legend.position = legendpos)

    return(distp)
}

groupedBoxplot <- function(df,
                           xvar,
                           yvar,
                           groupvar='replicate',
                           xlab=xvar,
                           ylab=yvar,
                           plotpoints=FALSE,
                           points_thresh=200, # show all points for xvar with N observations below this number
                           show_wilcox=FALSE,
                           paired_t_test=FALSE,
                           testout=FALSE,
                           alt='two.sided',
                           pvalsize=4,
                           nominal_p_only=FALSE,
                           ylimadj=0,
                           orderpval=FALSE,
                           logscale=FALSE,
                           colorX=FALSE,
                           colorvals=NULL,
                           showleg=FALSE,
                           legendpos='right',
                           legtextsize=12,
                           horizontal=TRUE,
                           facets=NULL,
                           frows=2,
                           fscales='fixed',
                           fdir='h',
                           w=7,h=7,
                           angleHjVj=c(20,.5,.5),
                           debug=FALSE) {

    options(repr.plot.width = w, repr.plot.height = h, repr.plot.res = 200)

    if (colorX) {
        colorvar <- xvar
    } else {
        colorvar <- groupvar
    }
    outliers <- 0.2
    if (plotpoints) outliers <- 0.0

    if (horizontal) {
      ytextjust <- 0.9
    } else {
      ytextjust <- 1.1
    }

    if (show_wilcox & is.null(groupvar)) {

        if (!is.null(facets)) {

            testTypes <- df %>%
                group_by(across(any_of(facets))) %>%
                dplyr::rename(xvarvals := {{xvar}}) %>%
                summarize(Nxvar = n_distinct(xvarvals)) %>%
                dplyr::filter(Nxvar == 2) %>%
                ungroup() %>%
                dplyr::select(any_of(facets)) %>%
                distinct()

            wilcox_df <- merge(df, testTypes)

        } else {
            wilcox_df <- df
        }

#         wilcox_df <- df %>%
#             group_by(across(any_of(c(facets,xvar)))) %>%
#             summarize(N = dplyr::n())

        wilcox_df <- wilcox_df %>%
            group_by(across(any_of(facets))) %>%
            dplyr::rename(wilcox_test_var := {{xvar}},
                          ycol := {{yvar}}) %>%
            mutate(maxy = max(ycol)*ytextjust) %>%
            group_by(maxy, .add=TRUE)

        if (is.null(facets)) {

            wilcox_df <- wilcox_df %>% mutate(wilcox_test_var = as.factor(wilcox_test_var))
            wvs <- levels(wilcox_df$wilcox_test_var)

            wv1 <- wilcox_df %>% filter(wilcox_test_var==wvs[1]) %>% pull(ycol)
            wv2 <- wilcox_df %>% filter(wilcox_test_var==wvs[2]) %>% pull(ycol)
            wres <- wilcox.test(wv2, wv1, alternative=alt, conf.int=TRUE)

            wilcox_df <- as_tibble(list(estimate = wres$estimate[[1]],
                                        group1 = wvs[2],
                                        group2 = wvs[1],
                                        n1 = length(wv2),
                                        n2 = length(wv1),
                                        statistic = wres$statistic[[1]],
                                        p = wres$p.value,
                                        conf.low = wres$conf.int[1],
                                        conf.high = wres$conf.int[2],
                                        method = 'Wilcoxon',
                                        alternative = alt,
                                        maxy = wilcox_df %>% pull(maxy) %>% unique))

        } else {

            if (paired_t_test) {
                wilcox_df <- wilcox_df %>%
                    rstatix::t_test(ycol ~ wilcox_test_var,
                                    paired=TRUE,
                                    p.adjust.method='none')
            } else {
                wilcox_df <- wilcox_df %>%
                    mutate(wilcox_test_var = fct_relevel(wilcox_test_var, rev))

                wvs <- levels(wilcox_df$wilcox_test_var)

                print('running wilcox tryCatch')

                wilcox_df <- tryCatch(
                    {
                        rstatix::wilcox_test(wilcox_df,
                                                ycol ~ wilcox_test_var,
                                                 alternative=alt,
                                                     detailed=TRUE,
                                                 p.adjust.method='none')
                    },
                    error=function(e) {
                        print(e)
                        return(
                            wilcox_df %>%
                                summarize(p = Inf) %>%
                                mutate(group1=wvs[2],
                                       conf.low=NA,
                                       conf.high=NA)
                        )
                    }
                )
            }

        }

        wilcox_df <- wilcox_df %>%
            mutate(!!xvar := group1,
                   xposn = 1.5,
                   yposn=maxy,
                   wilcox_lab = paste0('P = ', signif(p, digits=2)))

        if (orderpval) {
            wilcox_df <- wilcox_df %>%
                arrange(p) %>%
                mutate(across(all_of(facets), ~ fct_reorder(as_factor(.x), p)))

            df <- df %>%
                mutate(across(all_of(facets), ~ factor(.x, levels=unique(wilcox_df$.x))))

        }

    } else if (show_wilcox) {
        wilcox_df <- df %>%
            group_by(across(any_of(facets))) %>%
            dplyr::rename(wilcox_test_var := {{groupvar}},
                          ycol := {{yvar}}) %>%
            mutate(maxy = max(ycol)*ytextjust) %>%
            group_by(maxy, across(any_of(xvar)), .add=TRUE)

        if (paired_t_test) {
            wilcox_df <- wilcox_df %>%
                rstatix::t_test(ycol ~ wilcox_test_var,
                                paired=TRUE,
                                p.adjust.method='none')
        } else {
            wilcox_df <- tryCatch(
                {
                    wilcox_df <- rstatix::wilcox_test(wilcox_df,
                                                ycol ~ wilcox_test_var,
                                             alternative=alt,
                                                 detailed=TRUE,
                                             p.adjust.method='none')
                },
                error=function(e) {
                    print(e)
                    return(
                        wilcox_df %>%
                            summarize(p = Inf) %>%
                            mutate(group1=wvs[2],
                                   conf.low=NA,
                                   conf.high=NA)
                    )
                }
            )

        }

        wilcox_df <- wilcox_df %>%
            mutate(yposn=maxy,
                   wilcox_lab = paste0('P=', signif(p, digits=2)))
    }

    if (show_wilcox) {

        wilcox_df <- wilcox_df %>%
            mutate(padj = p*dplyr::n(),
                   wilcox_lab = case_when(padj < 0.00005 ~ paste0(wilcox_lab, '****'),
                                          padj < 0.0005 ~ paste0(wilcox_lab, '***'),
                                          padj < 0.005 ~ paste0(wilcox_lab, '**'),
                                          padj < 0.05 ~ paste0(wilcox_lab, '*'),
                                          TRUE ~ wilcox_lab))
        if (nominal_p_only) {

            wilcox_df <- wilcox_df %>%
                mutate(wilcox_lab = case_when(p>=0.05 ~ "",
                                              p<0.05 ~ wilcox_lab))

        }

        if (debug) return(wilcox_df)

    }

    if (debug) return(wilcox_df)


    if (!is.null(groupvar)) {
        p1 <- df %>%
            ggplot(aes(x=.data[[xvar]],
                       y=.data[[yvar]],
                       fill=.data[[colorvar]],
                       alpha=.data[[groupvar]]))
    } else {
        p1 <- df %>%
            ggplot(aes(x=.data[[xvar]],
                       y=.data[[yvar]],
                       fill=.data[[colorvar]]))
    }



    p1 <- p1 +
        geom_boxplot(outlier.alpha=outliers,
                     show.legend=showleg)

    if (plotpoints) {
        posd <- position_jitter()
        if (!is.null(groupvar)) posd <- position_jitterdodge(jitter.width=0.3,dodge.width=0.75)

        dfo <- df %>%
            group_by(across(all_of(c(facets,xvar)))) %>%
            dplyr::filter(dplyr::n() >= points_thresh) %>%
            ungroup()
        p1 <- p1 +
            geom_boxplot(data=dfo,
                         outlier.alpha=0.2,
                         show.legend=showleg)

        dfp <- df %>%
            group_by(across(all_of(c(facets,xvar)))) %>%
            dplyr::filter(dplyr::n() < points_thresh) %>%
            ungroup()
        p1 <- p1 +
            geom_point(data=dfp,
                       pch=21,
                       position = posd,
                       show.legend=FALSE)
    }


    if (show_wilcox) {

        # dummy df for plottin blank space for text labels
        ylim_blank <- df %>%
            group_by(across(all_of(facets))) %>%
            mutate(or_ylim = max(.data[[yvar]])*(1+ylimadj))

        if (is.null(groupvar)) {
            ylim_blank <- ylim_blank %>%
                dplyr::select(any_of(c(facets, xvar, yvar,'or_ylim')))

            p1 <- p1 +
                geom_text(data=wilcox_df,
                          aes(x=xposn,y=yposn,label=wilcox_lab), size=pvalsize, fontface='italic')
        } else {
            ylim_blank <- ylim_blank %>%
                dplyr::select(any_of(c(facets, groupvar, xvar, yvar,'or_ylim')))

            p1 <- p1 +
                geom_text(data=wilcox_df,
                          inherit.aes=FALSE,
                          aes(x=.data[[xvar]],y=yposn,label=wilcox_lab), size=pvalsize, fontface='italic')
        }

        if (ylimadj>0) p1 <- p1 + geom_blank(data=ylim_blank, aes(x=.data[[xvar]], y=or_ylim, color=NULL))

    }



    if (!is.null(facets)) {

       p1 <- p1 +
            facet_wrap(facets, nrow=frows, dir=fdir, scales=fscales)

    }


    if (!is.null(colorvals)) {

        p1 <- p1 +
            scale_fill_manual(values=colorvals) +
            scale_color_manual(values=colorvals)

    }

    if (logscale) {
        p1 <- p1 + scale_y_log10()
    }
    if (horizontal) p1 <- p1 + coord_flip()

    p1 <- p1 +
        scale_alpha_manual(values=rep(1, df[[colorvar]] %>% n_distinct())) +
        labs(y = ylab,
             x = xlab) +
        theme_classic() +
        theme(axis.text.x = element_text(size=16,angle=angleHjVj[1],hjust=angleHjVj[2],vjust=angleHjVj[3]),
              axis.text.y = element_text(size=16),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              plot.margin = unit(c(0,1,0,0),"line"),
              legend.position = legendpos,
              legend.title = element_text(size=legtextsize),
              legend.text = element_text(size=legtextsize),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_text(size=18),
              panel.spacing = unit(2, "line"))

    if (testout) {

        wilcoxout <- wilcox_df %>%
            select(-any_of(c('maxy','.y.','variant_type','xposn','yposn','wilcox_lab','padj'))) %>%
            rename(conf_lower = conf.low,
                   conf_upper = conf.high,
                   pvalue = p)

        return(list(p1, wilcoxout))

    }

    return(p1)
}

groupedBarPlot <- function(df,
                          xvar,
                          yvar,
                          xlab=xvar,
                          ylab=yvar,
                          xlaborder=NULL,
                          xlaborderval=NULL,
                          xlabfacetord=NULL,
                          xcondcolor=NULL,
                          groupvar=NULL,
                          groupvarlab=groupvar,
                          groupvarord=NULL,
                          colorvals=NULL,
                          horizontal=TRUE,
                          logscale=FALSE,
                          w=7,
                          h=7,
                          facets=NULL,
                          frows=2,
                          fscales='fixed', # free, free_x, free_y
                          fdir='h',
                          legendpos=NULL, # c(xpos, ypos)
                          debug=FALSE
                         ) {

    if (!is.null(groupvar)) {
        df <- df %>%
            dplyr::rename(gvc := {{groupvar}})
        if (!is.null(groupvarord)) {
            if (length(groupvarord)>1) {
                df <- df %>%
                    mutate(gvc = fct_relevel(as_factor(gvc), groupvarord))
            } else {
                df <- df %>%
                    dplyr::rename(gvord := {{groupvarord}}) %>%
                    mutate(gvc = as.factor(gvc),
                           gvc = fct_reorder(gvc, gvord)) %>%
                    dplyr::rename(!!groupvarord := gvord)
            }
        } else {
            df <- df %>%
                mutate(gvc = as.factor(gvc),
                       gvc = fct_relevel(gvc, rev))
        }
    } else {
        if ('hex' %in% colnames(df)) {
            df <- df %>%
                dplyr::mutate(gvc = hex)
        } else {
            gvc <- NULL
        }
    }

    df <- df %>%
        dplyr::rename(xcol := {{xvar}},
                      ycol := {{yvar}})

    if (is.null(xlaborder)) {
        df <- df %>%
            mutate(xcol = as_factor(xcol),
                   xcol = fct_reorder(xcol, desc(ycol)))
    } else {
        # subset df and order by specific grouping/facet vars
        # then merge df with new plot order column
        dfs <- df %>%
            dplyr::select(xcol,ycol,gvc,any_of(c(facets,xlaborderval))) %>%
            distinct() %>%
            dplyr::filter(gvc==xlaborder)

        if (!is.null(xlabfacetord)) {
            dfs <- dfs %>%
                dplyr::rename(fc1 := facets[1]) %>%
                dplyr::filter(fc1 == xlabfacetord)
        }

        if (is.null(xlaborderval)) {
            dfs <- dfs %>% arrange(desc(ycol))
        } else {
            dfs <- dfs %>%
                dplyr::rename(xlov := {{xlaborderval}})
            if (xlaborderval=='plot_order') {
                dfs <- dfs %>% arrange(xlov)
            } else {
                dfs <- dfs %>% arrange(desc(xlov))
            }
        }

        if (debug) return(dfs)

        dfs <- dfs %>%
            mutate(plot_order = seq(1,length(ycol))) %>%
            dplyr::select(-c(ycol,gvc))

        df <- merge(df, dfs) %>%
            mutate(xcol = fct_reorder(as_factor(xcol), plot_order))
    }

    options(repr.plot.width = w, repr.plot.height = h, repr.plot.res = 200)

    p1 <- df %>%
            ggplot(aes(x=xcol,
                       y=ycol,
                       fill=gvc)) +
        geom_bar(stat="identity",
                 position="dodge",
                 width=.75)

    if (horizontal) p1 <- p1 + coord_flip()

    p1 <- p1 +
        theme_bw() +
        labs(x=xlab,
             y=ylab)

    if (logscale) {
        p1 <- p1 + scale_y_log10()
    }

    if ('hex' %in% colnames(df)) {
        p1 <- p1 +
            scale_fill_identity(breaks=df$hex, labels=df$category, guide='legend')
    } else {
        if (!is.null(colorvals)) p1 <- p1 + scale_fill_manual(values=colorvals)
        if (!is.null(groupvar)) p1 <- p1 + labs(fill=groupvarlab) + guides(color=guide_legend(reverse=TRUE))
    }

    if (!is.null(facets)) {
        p1 <- p1 +
            facet_wrap(facets, nrow=frows, dir=fdir, scales=fscales)
    }

    if (is.null(legendpos)) {
        if (!is.null(facets)) {
            legendpos <- c(.4,.2)
        } else {
            legendpos <- c(.8,.2)
        }
    }

    if (is.null(groupvarlab)) {
        p1 <- p1 +
            theme(legend.title = element_blank(),
                  legend.text = element_text(size=16))
    } else {
        p1 <- p1 +
            labs(fill=groupvarlab) +
            theme(legend.title = element_text(size=20),
                  legend.text = element_text(size=18))
    }

    if (!is.null(xcondcolor)) {
        xcolor <- ifelse(grepl(xcondcolor, levels(df$xcol)), 'black', 'darkgray')

    } else {
        xcolor <- 'gray'
    }

    p1 <- p1 +
        theme(axis.text.x = element_text(size=22,angle=20,hjust=0.9,vjust=0.9,color=xcolor),
              axis.text.y = element_text(size=22),
              axis.title.x = element_text(size=28),
              axis.title.y = element_text(size=28),
              plot.margin = unit(c(0,1,0,0),"line"),
              legend.position=legendpos,
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_text(size=22),
              panel.spacing = unit(2, "line"))

    return(p1)

}

plotCorrelation <- function(df,
                            xvar,
                            yvar,
                            xlabel=xvar,
                            xangle=0,
                            ylabel=yvar,
                            groupvar='ancestry',
                            groupvarlab=groupvar,
                            colorvals=c("#2E3191","#F05A28"),
                            addline=NULL, # c(intercept, slope)
                            legendpos=c(.8,.2),
                            corrpos=NULL,
                            corrcolor='black',
                            w=7,h=7,
                            logscale=FALSE,
                            xlogonly=FALSE,
                            facets=NULL,
                            fscales='fixed', # free, free_x, free_y
                            frows=2,
                            fdir='h',
                            facetorder=TRUE,
                            orderpval=FALSE,
                            addaxes=FALSE,
                            debug=FALSE) {

    RSQ <- function(y, x, zero_int=TRUE) {
        if (zero_int) {
            lms <- summary(lm(y ~ 0 + x))
            pvidx <- 1
        } else {
            lms <- summary(lm(y ~ x))
            pvidx <- 2
        }
        rsq <- signif(lms$r.squared, 2)
        pv <- signif(lms$coefficients[pvidx,'Pr(>|t|)'], 2)
        return(list(Rsq=rsq, P=pv))
    }

    Npre <- base::nrow(df)
    df <- df %>% dplyr::filter(if_all(c(xvar,yvar), ~ is.finite(.)))
    Npost <- base::nrow(df)
    print(paste0(Npre-Npost," potentially repeated, non-finite elements dropped"))

    if (is.null(facets)) {
        if(logscale) {
            if (xlogonly) {
                lmform <- paste(yvar,'~',paste0("log10",xvar))
                df[paste0("log10",xvar)] <- log10(df[[xvar]])
                dfsub <- df[(df[[paste0("log10",xvar)]]<=0)
                             & is.finite(df[[paste0("log10",xvar)]]),]

            } else {
                lmform <- paste(paste0("log10",yvar),'~',paste0("log10",xvar))
                df[paste0("log10",yvar)] <- log10(df[[yvar]])
                df[paste0("log10",xvar)] <- log10(df[[xvar]])
                dfsub <- df[(df[[paste0("log10",yvar)]]<=0)
                             &(df[[paste0("log10",xvar)]]<=0)
                             & is.finite(df[[paste0("log10",xvar)]])
                             & is.finite(df[[paste0("log10",yvar)]]),]
            }

            m <- lm(eval(parse(text=lmform)),
                    dfsub
                   )
        } else {
            lmform <- paste(yvar,'~',xvar)
            m <- lm(eval(parse(text=lmform)), df)
        }
        aR <- signif(summary(m)$r.squared, 2)
        aP <- signif(summary(m)$coefficients[2,'Pr(>|t|)'], 2)
        annotlab <- paste0("Rsq=",aR,", P=",aP)
    } else {
        if (logscale) {
            if (xlogonly) {
                df[paste0("log10",xvar)] <- log10(df[[xvar]])
                dfsub <- df[(df[[paste0("log10",xvar)]]<=0)
                             & is.finite(df[[paste0("log10",xvar)]]),]
                sdf <- dfsub %>%
                    group_by(across(all_of(facets))) %>%
                    summarize(Rsq=RSQ(.data[[yvar]],.data[[paste0("log10",xvar)]])$Rsq,
                              P=RSQ(.data[[yvar]],.data[[paste0("log10",xvar)]])$P,
                              nObs=dplyr::n(),
                              xpos=mean(.data[[xvar]]),
                              ypos=min(.data[[yvar]]) - 0.15*abs(min(.data[[yvar]]))) %>%
                    ungroup()
                if (debug) return(sdf)

            } else {
                df[paste0("log10",yvar)] <- log10(df[[yvar]])
                df[paste0("log10",xvar)] <- log10(df[[xvar]])
                dfsub <- df[(df[[paste0("log10",yvar)]]<=0)
                             &(df[[paste0("log10",xvar)]]<=0)
                             & is.finite(df[[paste0("log10",xvar)]])
                             & is.finite(df[[paste0("log10",yvar)]]),]
                sdf <- dfsub %>%
                    group_by(across(all_of(facets))) %>%
                    summarize(Rsq=RSQ(.data[[paste0("log10",yvar)]],.data[[paste0("log10",xvar)]])$Rsq,
                              P=RSQ(.data[[paste0("log10",yvar)]],.data[[paste0("log10",xvar)]])$P,
                              nObs=dplyr::n(),
                              xpos=mean(.data[[xvar]]),
                              ypos=min(.data[[yvar]]) - 0.15*abs(min(.data[[yvar]]))) %>%
                    ungroup()
            }
            sdf <- sdf %>%
                mutate(RsqL=paste0("Rsq=",Rsq,", P=",P),
                       xposn=0,
                       yposn=Inf)


        } else {
            sdf <- df %>%
                group_by(across(all_of(facets))) %>%
                summarize(Rsq=RSQ(.data[[yvar]],.data[[xvar]])$Rsq,
                          P=RSQ(.data[[yvar]],.data[[xvar]])$P,
                          nObs=dplyr::n(),
                          xpos=mean(.data[[xvar]]),
                          ypos=min(.data[[yvar]]) - 0.15*abs(min(.data[[yvar]]))) %>%
                mutate(RsqL=paste0("Rsq=",Rsq,", P=",P),
                       xposn=-Inf,
                       yposn=Inf) %>%
                ungroup()
        }

        if (facetorder) {
            if (orderpval) {
                sdf <- sdf %>%
                    dplyr::arrange(P) %>%
                    mutate(across(all_of(facets), ~ fct_reorder(as_factor(.x), P)))

            } else {
                sdf <- sdf %>%
                    dplyr::arrange(desc(Rsq)) %>%
                    mutate(across(all_of(facets), ~ fct_reorder(as_factor(.x), Rsq, .desc=TRUE)))
            }
        }

        df <- df %>%
            mutate(across(all_of(facets), ~ factor(.x, levels=unique(sdf$.x))))

    }

    if (is.null(corrpos)) {
        meany <- df %>% pull(!!yvar) %>% min()
        meanx <- df %>% pull(!!xvar) %>% mean()
        corrpos <- c(meanx, meany - 0.1*abs(meany))
    }

    options(repr.plot.width = w, repr.plot.height = h, repr.plot.res = 200)

    p1 <- ggplot(df, aes(.data[[xvar]], .data[[yvar]])) +
        geom_point(aes(fill=.data[[groupvar]]), shape=21, size=2) +
        geom_smooth(method=lm,color='gray',alpha=0.1)

    if (!is.null(addline)) {
        p1 <- p1 +
            geom_abline(intercept=addline[1],slope=addline[2],linetype='dashed')
    }
    if (addaxes) {
        p1 <- p1 +
            geom_hline(yintercept=0,size=0.5) +
            geom_vline(xintercept=0,size=0.5)
    }

    p1 <- p1 +
        theme_classic() +
        labs(x = xlabel,
             y = ylabel,
             fill=groupvarlab) +
        theme(axis.text.x = element_text(size=18, angle=xangle,vjust=0.5),
              axis.text.y = element_text(size=20),
              axis.ticks.length.x = unit(0.75, "line"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              legend.position = legendpos,
              legend.background = element_rect(fill = "transparent",colour = NA),
              legend.box.background = element_rect(fill = "transparent",colour = "black"),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA),
              plot.margin = unit(c(0,0,0,0),"line"))

    if (is.null(facets)) {
        p1 <- p1 +
            annotate("text",
                     label=annotlab,
                     x=corrpos[1],
                     y=corrpos[2],
                     color=corrcolor) +
            theme(legend.title = element_text(size=18),
                  legend.text = element_text(size=16),
                  legend.key.size = unit(1.0, "line"))
    } else {
        p1 <- p1 +
            facet_wrap(as.formula(paste0('~',paste(facets,collapse='+'))), nrow=frows, dir=fdir, scales=fscales) +
            geom_text(data=sdf, aes(x=xpos,y=ypos,label=RsqL), size=6, color=corrcolor) +
            geom_text(data=sdf, aes(x=xposn,y=yposn,hjust=-0.1,vjust=1/.8,label=paste0('N=',nObs)), size=6) +
            theme(legend.title = element_text(size=18),
                  legend.text = element_text(size=16),
                  legend.key.size = unit(2.0, "line"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size=22)) +
            guides(color = guide_legend(title.position='left', nrow=frows, byrow=TRUE))
    }

    if(!is.null(colorvals)) {
       p1 <- p1 + scale_fill_manual(values=colorvals)
    } else if (!is.null(groupvar) & n_distinct(df[[groupvar]])>1) {
       p1 <- p1 + scale_fill_viridis_d(option = "D", begin=0.5, direction=-1)
    }

    if (logscale) {
        if (xlogonly) {
            p1 <- p1 + scale_x_log10()
        } else {
            p1 <- p1 + scale_x_log10() +  scale_y_log10()
        }

    }

    return(p1)

}

plotVolcano <- function(df,
                        lfc_var='log2FoldChange',
                        p_var='p',
                        pthresh='bonf',
                        nominal_p=FALSE,
                        lfcthresh=NULL,
                        label_var='TargetGene',
                        label_above=NULL, # c(xval, yval) treats as absolute value
                        xlabel='Log2 Fold Change',
                        ylabel='-log10(P)',
                        groupvar=NULL,
                        groupvarlab=groupvar,
                        colorvals=NULL, # c("#2E3191","#F05A28") for blue, orange
                        legendpos=c(.8,.2),
                        w=7,h=7) {

    options(repr.plot.width = w, repr.plot.height = h, repr.plot.res = 200)

    if (pthresh=='bonf') {
        pthresh = 0.05/nrow(df)
    }

    df <- df %>%
        dplyr::rename(lfc := {{lfc_var}},
                      pval := {{p_var}},
                      label_col := {{label_var}}) %>%
        mutate(nlp = -log10(pval))

    p1 <- ggplot(df, aes(lfc, nlp, label=label_col)) +
        geom_point(aes(color=.data[[groupvar]]), shape=21, size=2) +
        geom_vline(xintercept=0,linetype='dashed',color='gray') +
        geom_hline(yintercept=-log10(pthresh),linetype='dashed',color='gray') +
        labs(x = xlabel,
             y = ylabel,
             color=groupvarlab) +
        theme_classic() +
        theme(axis.text.x = element_text(size=18, angle=90,vjust=0.5),
              axis.text.y = element_text(size=20),
              axis.ticks.length.x = unit(0.75, "line"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              legend.position = legendpos,
              legend.title = element_text(size=18),
              legend.text = element_text(size=16),
              legend.key.size = unit(1.0, "line"),
              legend.background = element_rect(fill = "transparent",colour = NA),
              legend.box.background = element_rect(fill = "transparent",colour = NA),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA),
              plot.margin = unit(c(0,0,0,0),"line"))

    if(!is.null(colorvals)) {
       p1 <- p1 + scale_fill_manual(values=colorvals)
    }
    if(!is.null(lfcthresh)) {
        p1 <- p1 +
            geom_vline(xintercept=lfcthresh,linetype='dashed',color='gray') +
            geom_vline(xintercept=-lfcthresh,linetype='dashed',color='gray')
    }
    if (nominal_p) {
        p1 <- p1 +
            geom_hline(yintercept=-log10(0.05),linetype='dashed',color='gray')
    }
    if (!is.null(label_above)) {
        p1 <- p1 +
            ggrepel::geom_text_repel(
                data = df %>% dplyr::filter((lfc > label_above[1]) | (lfc < -label_above[1]),
                                            nlp > label_above[2])
            )
    }


    return(p1)

}

plotDirectionByPval <- function(df,
                                direction_var='log2FoldChange',
                                p_var='p',
                                pvs=c(1,.5,0.05,0.005,5e-4,5e-5,0),
                                colorvals=c("#2E3191","#F05A28"),
                                updownlabs=NULL # c(label_for_up_direction, label_for_down_direction)
                               ) {

    dfs <- df %>%
        dplyr::rename(dvar := {{direction_var}},
                      pval := {{p_var}}) %>%
        group_by(pvalue=cut(pval, breaks=pvs, include.lowest=TRUE)) %>%
        summarize(N_up=sum(dvar > 0),
                  N_down=sum(dvar < 0)) %>%
        pivot_longer(cols=dplyr::starts_with("N_"),
                     names_to="direction",
                     names_prefix="N_",
                     values_to="count")

    if (!is.null(updownlabs)) {
        dfs <- dfs %>%
            mutate(direction = case_when(direction=='up' ~ updownlabs[1],
                                         direction=='down' ~ updownlabs[2]))
    }

    p1 <- ggplot(dfs, aes(x=pvalue, y=count, fill=direction)) +
        geom_bar(stat="identity", position="dodge", width=0.7) +
        geom_text(position=position_dodge(width=.7),aes(y=count+10,label=count,hjust=0)) +
        scale_y_continuous(expand = c(.1, .1)) +
        scale_fill_manual(values=c("#2E3191","#F05A28")) +
        coord_flip() +
        labs(
            y = "count",
            x = "P-value window",
            fill = "Direction") +
        theme_classic() +
        theme(axis.text.x = element_text(size=20),
              axis.text.y = element_text(size=16),
        #           axis.ticks.length.x = unit(0.75, "line"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              legend.position = c(.75,.22),
              legend.title = element_text(size=20),
              legend.text = element_text(size=16),
              legend.key.size = unit(1.75, "line"),
              legend.background = element_rect(fill = "transparent",colour = NA),
        #           legend.box.background = element_rect(fill = "transparent",colour = NA),
              legend.box.background = element_rect(colour = "black"),
        #           panel.background = element_rect(fill = "transparent",colour = NA)
                  plot.margin = unit(c(0,0,0,0),"line")
        #           plot.background = element_rect(fill = "transparent",colour = NA)
             )

    return(p1)
}

plotFisherEnrichments <- function(df,
                                  yvar_col='condition',
                                  ylabel=yvar_col,
                                  xvar_col='odds_ratio',
                                  xlabel=xvar_col,
                                  althypoth='greater',
                                  success_total=FALSE, # display total success number instead of total elements in test
                                  line_intercept=1,
                                  groupvar=NULL,
                                  groupvarlab=groupvar,
                                  groupvarord=NULL,
                                  colorvals=NULL,
                                  sep_groups=TRUE,
                                  dodgewidth=0.4,
                                  w=7,
                                  h=7,
                                  facets=NULL,
                                  frows=2,
                                  fscales='fixed', # free, free_x, free_y
                                  fdir='h',
                                  flabpos='top',
                                  ngenessize=6,
                                  ngenes_lim_adj=6.4,
                                  ngenes_or_adj=.3,
                                  legendpos=NULL, # c(xpos, ypos)
                                  legtextsize=22,
                                  debug=FALSE
                                 ) {

    if ('fg_total' %in% colnames(df)) {
        df <- df %>% mutate(nGenes = fg_total + bg_total)
    }
    if ('DE_up_DA_up' %in% colnames(df)) {
        df <- df %>% dplyr::rename(up_up = DE_up_DA_up,
                                   up_down = DE_up_DA_down,
                                   down_up = DE_down_DA_up,
                                   down_down = DE_down_DA_down)
    }
    if ('up_up' %in% colnames(df)) {
        if (success_total) {
            df <- df %>% mutate(nGenes = up_up + down_up)
        } else {
            df <- df %>% mutate(nGenes = up_up + up_down + down_up + down_down)
        }

    }
    if ('n_trials' %in% colnames(df)) {
        if (any(!is.na(df['n_trials']))) {
            df <- df %>% mutate(nGenes = n_trials)
        }
    }

    # get bonferroni P to display asterisk significance
    df <- df %>%
        mutate(padj = pvalue*dplyr::n(),
               nGenes = as.character(nGenes),
               nGenes = case_when(padj < 0.00005 ~ paste0(nGenes, '****'),
                                  padj < 0.0005 ~ paste0(nGenes, '***'),
                                  padj < 0.005 ~ paste0(nGenes, '**'),
                                  padj < 0.05 ~ paste0(nGenes, '*'),
                                  TRUE ~ nGenes))

    df <- df %>% dplyr::rename(yvar := {{yvar_col}},
                               odds_ratio := {{xvar_col}})
    if (!is.null(groupvar)) {
        df <- df %>%
            dplyr::rename(gvc := {{groupvar}})
        if (is.null(facets)) {
            if (!is.null(groupvarord)) {
                df <- df %>%
                    mutate(gvc = as.factor(gvc),
                           gvc = fct_reorder(gvc, groupvarord))
            } else {
                df <- df %>%
                    mutate(gvc = as.factor(gvc),
                           gvc = fct_relevel(gvc, rev))
            }
        }

        posd <- position_dodge(dodgewidth)
        ewidth <- 0.2
    } else {
        if ('hex' %in% colnames(df)) {
            df <- df %>%
                dplyr::mutate(gvc = hex)
        } else {
            gvc <- NULL
        }
        posd <- 'identity'
        ewidth <- 0.1
    }
    if (althypoth=='greater') {
        df <- df %>% mutate(ymaxcol = odds_ratio)
    } else {
        df <- df %>% mutate(ymaxcol = conf_upper)
    }

    options(repr.plot.width = w, repr.plot.height = h, repr.plot.res = 200)

    df <- df %>%
        group_by(across(all_of(facets))) %>%
        mutate(yvar = as.factor(yvar),
               yvar = fct_reorder(yvar, conf_lower),
               or_inf = max(ymaxcol[is.finite(ymaxcol)]))

    # if x-axis fixed across facets, perform same adjustment to accommodate ngenes labels in each
    if (fscales != 'free_x') df <- df %>% ungroup()
    if (max(df[is.finite(df$ymaxcol),]$ymaxcol) - min(df[is.finite(df$ymaxcol),]$ymaxcol) > 100) {
        logscale <- TRUE
        df <- df %>%
            mutate(or_offset = if_else(is.finite(ymaxcol), ngenes_or_adj*ymaxcol, ngenes_or_adj*or_inf))
    } else {
        logscale <- FALSE
        df <- df %>%
            mutate(or_offset = max((.15*ngenes_or_adj)*(max(ymaxcol[is.finite(ymaxcol)]) - min(ymaxcol)), (.15*ngenes_or_adj)))
    }

    # amount of space to add within plotting range for points/significance labels to be visible
    ylimadj <- 2
    if ('nGenes' %in% colnames(df)) {
        if (max(nchar(df$nGenes)) > 2) ylimadj <- ngenes_lim_adj
    }
    df <- df %>%
        mutate(or_ylim = max(ymaxcol[is.finite(ymaxcol)]) + max(or_offset)*ylimadj)

    if (debug) return(df)

    if (!is.null(facets) & sep_groups) {
        labexp <- 0.5
        df <- df %>%
            ungroup() %>%
            dplyr::rename(yvar_names = yvar) %>%
            group_by(across(all_of(facets))) %>%
            arrange(conf_lower, .by_group=TRUE) %>%
            ungroup() %>%
            mutate(yvar = row_number()) %>%
            group_by(across(all_of(facets))) %>%
            mutate(labs_xlim_lower = min(yvar) - labexp,
                   labs_xlim_upper = max(yvar) + labexp) %>%
            ungroup()
    }

    # dummy df for plottin blank space for text labels
    ylim_blank <- df %>%
        group_by(across(all_of(facets))) %>%
        dplyr::select(yvar,any_of(c(facets,'odds_ratio','or_ylim')))

    p1 <- ggplot(df, aes(x=yvar, y=odds_ratio, color=gvc)) +
        geom_errorbar(aes(ymin=conf_lower, ymax=ymaxcol), width=ewidth, position=posd) +
        geom_point(position=posd, size=4)

    if ('nGenes' %in% colnames(df)) {

        p1 <- p1 + geom_text(position=posd,
                             size=ngenessize,
                             aes(y=ymaxcol+or_offset,label=nGenes,hjust=0),
                             show.legend=FALSE)

    }

    p1 <- p1 +
        geom_blank(data=ylim_blank, aes(x=yvar, y=or_ylim, color=NULL)) +
        geom_hline(yintercept=line_intercept, linetype='dashed', color='gray') +
        coord_flip() +
        theme_classic() +
        labs(x=ylabel,
             y=xlabel)

    if ('nGenes' %in% colnames(df)) {

        if (!all(is.finite(df$odds_ratio))) {
            p1 <- p1 +
                geom_text(data=df[!is.finite(df$odds_ratio),],
                          position=posd,size=ngenessize,
                          aes(y=or_inf+or_offset,label=nGenes,hjust=0,vjust=-0.5))
        }

    }



    if (logscale) {
        p1 <- p1 + scale_y_log10() + theme(plot.margin = unit(c(0,2,0,0),"line"))
    } else {
        p1 <- p1 + theme(plot.margin = unit(c(0,1.5,0,0),"line"))
    }

    if ('hex' %in% colnames(df)) {
        p1 <- p1 +
            scale_color_identity(breaks=df$hex, labels=df$category, guide='legend')
    } else {
        if (!is.null(colorvals)) p1 <- p1 + scale_color_manual(values=colorvals)
        if (!is.null(groupvar)) p1 <- p1 + labs(color=groupvarlab) + guides(color=guide_legend(reverse=TRUE))
    }

    if (!is.null(facets)) {
        p1 <- p1 +
            facet_wrap(facets, nrow=frows, dir=fdir, scales=fscales, strip.position=flabpos)
        if (sep_groups) {
            p1 <- p1 +
                geom_blank(data=df, aes(x=labs_xlim_lower, y=or_ylim, color=NULL)) +
                geom_blank(data=df, aes(x=labs_xlim_upper, y=or_ylim, color=NULL)) +
                scale_x_continuous(breaks=df$yvar, labels=df$yvar_names, expand=c(0,0))
        }

    }

    if (is.null(legendpos)) {
        if (!is.null(facets)) {
            legendpos <- c(.4,.2)
        } else {
            legendpos <- c(.8,.2)
        }
    }

    if (is.null(groupvarlab)) {
        p1 <- p1 +
            theme(legend.title = element_blank(),
                  legend.text = element_text(size=legtextsize))
    } else {
        p1 <- p1 +
            labs(color=groupvarlab) +
            theme(legend.title = element_text(size=legtextsize+2),
                  legend.text = element_text(size=legtextsize))
    }

    p1 <- p1 +
        theme(axis.text.x = element_text(size=22),
              axis.text.y = element_text(size=22),
              axis.title.x = element_text(size=28),
              axis.title.y = element_text(size=28),
              legend.position=legendpos,
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_text(size=22),
              panel.spacing = unit(2, "line"))

    return(p1)

}

envCellTypeGridPlot <- function(df,
                                p1color_table,
                                p1filt='bulk',
                                p2filt='pseudobulk',
                                gtitle="Diff-score - DE gene enrichments",
                                dim1=c(10,10),
                                dim2=c(10,10),
                                ngenes_lim_adj=6.4,
                                debug=FALSE) {

    p1 <- plotFisherEnrichments(df %>%
                                    dplyr::filter(datatype==p1filt) %>%
                                    merge(., p1color_table),
                            yvar_col = "condition",
                            ylabel = "Condition",
                            althypoth = "greater",
                            groupvar = NULL,
                            legendpos = 'top',
                            groupvarord = NULL,
                            colorvals = NULL,
                            w = dim1[1], h = dim1[2],
                            facets = 'score_type',
                            frows = 2,
                            fscales = "free_y",
                               ngenes_lim_adj=ngenes_lim_adj,
                            debug=debug)

    p2 <- plotFisherEnrichments(df %>%
                                dplyr::filter(datatype==p2filt),
                            yvar_col = "celltype",
                            ylabel = "Cell type", althypoth = "greater",
                            groupvar = "condition",
                            groupvarlab = NULL,
                            colorvals = c('darkred','darkgray'),
                            legendpos = 'top', w = dim2[1], h = dim2[2],
                            facets = "score_type", frows = 2,
                            fscales = "free_y",
                            ngenes_lim_adj=ngenes_lim_adj,
                            debug=debug)

    if (debug) return(p1)

    # if(!(require(gridExtra) )) install.packages("gridExtra")

    options(repr.plot.width = dim1[1] + dim2[1], repr.plot.height = dim1[2] + dim2[2], repr.plot.res = 200)

    g1 <- gridExtra::grid.arrange(
      p1 + theme(plot.margin = unit(c(0,3,0,0),"line")),
      p2 + theme(plot.margin = unit(c(0,1.5,0,0),"line")),
      nrow = 1,
      top = grid::textGrob(
           gtitle,
           gp = grid::gpar(fontsize = 26)
      ),
      padding = unit(3, "line")
    )

    return(g1)

}

plotABCscores <- function(df, tgname=NULL) {

    # tgname is c(TargetGene,name)

    if (!is.null(tgname)) {
        df <- df %>% dplyr::filter(TargetGene==tgname[1], name==tgname[2])
    }

    if (!('ancestry' %in% colnames(df))) {
        afrpops <- c('ESN','GWD','LWK','YRI')
        eurpops <- c('CEU','FIN','IBS','TSI')

        df <- df %>%
            mutate(ancestry = case_when(population %in% afrpops ~ 'AFR',
                                        population %in% eurpops ~ 'EUR'))
    }

    df <- df %>%
        mutate(population=fct_relevel(population,'CEU','FIN','IBS','TSI','ESN','GWD','LWK','YRI'),
               sample = fct_relevel(sample,
                                    'CEU_rep1',
                                    'CEU_rep2',
                                    'FIN_rep1',
                                    'FIN_rep2',
                                    'IBS_rep1',
                                    'IBS_rep2',
                                    'TSI_rep1',
                                    'TSI_rep2',
                                    'ESN_rep1',
                                    'ESN_rep2',
                                    'GWD_rep1',
                                    'GWD_rep2',
                                    'LWK_rep1',
                                    'LWK_rep2',
                                    'YRI_rep1',
                                    'YRI_rep2'))

    global_labeller <- labeller(
      name = label_value,
      TargetGene = label_value,
      .default = label_both
    )

    p1 <- ggplot(df, aes(sample,ABC.Score)) +
        geom_point(aes(fill=ancestry,size=hic_contact_pl_scaled_adj), shape=21) +
        theme_classic() +
        scale_fill_manual(values=c("#2E3191","#F05A28")) +
#         guides(color = guide_legend(reverse = TRUE)) +
        labs(x = "Sample",
             y = "ABC Score",
             fill='Ancestry',
             size='Normalized HiChIP') +
        ggrepel::geom_text_repel(
                aes(label=hic_contact),
                size=4
#                 hjust='left',
#                 direction="y",
#                 segment.alpha=0
            ) +
        theme(axis.text.x = element_text(size=18, angle=90,vjust=0.5),
              axis.text.y = element_text(size=20),
              axis.ticks.length.x = unit(0.75, "line"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_text(size=18),
#               legend.position = c(.95,.65),
              legend.position = 'right',
              legend.title = element_text(size=22),
              legend.text = element_text(size=20),
              legend.key.size = unit(2.0, "line"),
              legend.background = element_rect(fill = "transparent",colour = NA),
              legend.box.background = element_rect(fill = "transparent",colour = NA),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA),
              plot.margin = unit(c(0,0,0,0),"line")) +
        facet_wrap(vars(name,
                        TargetGene,
                        distance),
                   labeller=global_labeller)

    return(p1)

}

abcPlotMA <- function(df,
                      prefix=NULL,
                      depthcol1='Score_AFR',
                      depthcol2='Score_EUR',
                      depthtype='Score',
                      lfc_col='log2FoldChange',
                      sig_col='p',
                      sig_thresh=0.05,
                      facets="isSelfPromoter",
                      fscales="free_x",
                      frows=1,
                      fdir='h',
                      w=14,
                      h=7,
                      debug=FALSE) {

    if (!is.null(prefix)) {
        depthcol1 <- paste0(prefix,'.Score_AFR')
        depthcol2 <- paste0(prefix,'.Score_EUR')
        depthtype <- paste0(prefix,'.Score')
        lfc_col <- paste0(prefix,'.log2FoldChange')
        sig_col <- paste0(prefix,'.p')

        df <- df %>%
            dplyr::filter(across(starts_with(prefix), ~ (!is.na(.x) & is.finite(.x))))
    } else {
        df <- df %>%
            dplyr::filter(across(all_of(c(lfc_col,sig_col)), ~ (!is.na(.x) & is.finite(.x))))
    }

    dfp <- df %>%
        mutate("mean_{depthtype}" := (.data[[depthcol1]] + .data[[depthcol2]])/2,
               isSelfPromoter = if_else(isSelfPromoter,'self-promoter','non-self-promoter'),
               significant = if_else(.data[[sig_col]] < sig_thresh, TRUE, FALSE)) %>%
        arrange(desc(.data[[sig_col]]))

    p1 <- plotCorrelation(dfp,
                    paste0('mean','_',depthtype),
                    lfc_col,
                    groupvar = "significant",
                    groupvarlab=paste0(sig_col,' < ',sig_thresh),
                    colorvals = c("gray","red"),
                    addline = c(0,0),
                    legendpos = c(.97,.35),
                    corrpos = c(.25,2.5),
                    corrcolor = "black",
                    w = w, h = h,
                    logscale = TRUE,
                    xlogonly=TRUE,
                    facets = facets,
                    fscales=fscales,
                    frows = frows,
                    fdir = fdir,
                    facetorder=FALSE,
                    orderpval=FALSE,
                    debug=FALSE)

    return(p1)

}

facetAbcPlotMA <- function(df, w=27, h=15) {

    p1 <- abcPlotMA(df,
                  prefix=NULL,
                  depthcol1='Score_AFR',
                  depthcol2='Score_EUR',
                  depthtype='Score',
                  lfc_col='log2FoldChange',
                  sig_col='p',
                  sig_thresh=0.05,
                  facets=c('score_type','isSelfPromoter'),
                  fscales='free',
                  frows=2,
                  fdir='v',
                  w=w,
                  h=h,
                  debug=FALSE) +
        theme(plot.margin = unit(c(0,2,0,0),"line"))

    return(p1)

}

gridAbcPlotMA <- function(df, grid_col='score_type', debug=FALSE) {

    # need to specify name of plot when save grid.arrange with ggsave,
    # otherwise it will only save the last plot in the grid, e.g.,
    # ggsave(file.path(outdir, paste0(outbase, ".componentScoreMAplot.png")), plot=p1, width=14, height=28)

    if(!(require(gridExtra) )) install.packages("gridExtra")
    library(gridExtra)

    df <- df %>%
        dplyr::rename(grid_var := {{grid_col}})
    pfxs <- df %>% pull(grid_var) %>% unique()

    dfw <- df %>%
        pivot_wider(names_from=grid_var,
                    names_glue="{grid_var}.{.value}",
                    values_from=c(Score_AFR, Score_EUR, log2FoldChange, p))

    plotlist <- lapply(pfxs, abcPlotMA, df=dfw, sig_thresh=0.05)
    pmargin <- theme(plot.margin = unit(c(1,0,1,0),"line"))
    plotlist <- lapply(plotlist, "+", pmargin)

    if (debug) return(plotlist)

    options(repr.plot.width = 14, repr.plot.height = 28, repr.plot.res = 200)
    p2 <- grid.arrange(grobs=plotlist, ncol=1)

    return(p2)

}

plotTopPways <- function(res, topN=5, cwidths=c(5, 3, 0.6, .8, .8)) {

    library(fgsea)

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

plotRoc <- function(df, facets=NULL) {

    #' Plot ROC curve
    #'
    #' Plot ROC curve for each optionally included grouping (facet) var
    #'
    #' @param df Dataframe produced by yardstick::roc_curve()
    #' @param facets Plot ROC curve for each of these groupings (e.g., for df produced by roc_curve on grouped data)
    #' @return ggplot2 plot object

    p1 <- df %>%
        ggplot(aes(x = 1 - specificity, y = sensitivity)) +
            geom_path() +
            geom_abline(lty = 3) +
            coord_equal() +
            theme_classic() +
            theme(axis.text.x = element_text(size=20),
                  axis.text.y = element_text(size=20),
                  axis.ticks.length.x = unit(0.75, "line"),
                  axis.title.x = element_text(size=22),
                  axis.title.y = element_text(size=22),
                  panel.background = element_rect(fill = "transparent",colour = NA),
                  plot.background = element_rect(fill = "transparent",colour = NA),
                  plot.margin = unit(c(0,0,0,0),"line"))

    if (!is.null(facets)) {
        p1 <- p1 +
            facet_wrap(as.formula(paste0('~',paste(facets,collapse='+'))), nrow=1, dir='h', scales='fixed') +
            theme(strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size=22))
    }

}
