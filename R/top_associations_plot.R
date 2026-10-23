#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Plot the abundances of the top associated features against the label
#'
#' @description This function plots the abundances of the top \code{n}
#' features (ranked by p-value) from \link{check.associations} against the
#' label. Binary labels are shown as boxplots, continuous labels as
#' scatterplots with a linear fit.
#'
#' @usage top.associations.plot(siamcat, n=5, feature.type=NULL, fn.plot=NULL,
#' ncol=NULL, font.size=14, feat.renamer=NULL, title=NULL, width=12, height=7)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param n integer, number of top features to plot, defaults to \code{5}
#'
#' @param feature.type string, which feature values to plot. One of
#' \code{"original"}, \code{"filtered"}, or \code{"normalized"}. If
#' \code{NULL}, the feature type used in \link{check.associations} is used,
#' so that the plotted values are the ones that were tested.
#'
#' @param fn.plot string, filename for the plot (any extension supported by
#' ggsave is allowed). If \code{fn.plot} is \code{NULL}, the plot will only
#' be returned as a ggplot object.
#'
#' @param ncol integer, number of facet columns, passed to facet_wrap
#'
#' @param font.size integer, base font size for the plot
#'
#' @param feat.renamer function to use to rename the features for plotting.
#' Defaults to \code{NULL} where no renaming is applied.
#'
#' @param title string or boolean or NULL, title for the plot.
#' If \code{NULL}, the label name will be used.
#' If \code{FALSE}, no title will be added.
#'
#' @param width numeric, width of the plot in inches, defaults to 12
#'
#' @param height numeric, height of the plot in inches, defaults to 7
#'
#' @return Returns the ggplot plot object
#'
#' @keywords SIAMCAT top.associations.plot
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
#' # Example data
#' data(siamcat_example)
#'
#' top.associations.plot(siamcat_example, n = 5, fn.plot = "./top5.pdf")
top.associations.plot <- function(
    siamcat, n = 5, feature.type = NULL, fn.plot = NULL, ncol = NULL,
    font.size = 14, feat.renamer = NULL, title = NULL,
    width = 12, height = 7
) {
    associations <- associations(siamcat, verbose = 0)
    if (is.null(associations)) {
        stop(
            "SIAMCAT object does not contain association testing results! ",
            "Exiting..."
        )
    }
    assoc.param <- assoc_param(siamcat)
    if (is.null(feature.type)) feature.type <- assoc.param$feature.type

    feat <- switch(feature.type,
        "original" = get.orig_feat.matrix(siamcat),
        "filtered" = get.filt_feat.matrix(siamcat),
        "normalized" = get.norm_feat.matrix(siamcat),
        stop("feature.type must be one of 'original', 'filtered', or 'normalized'")
    )

    associations <- associations[order(associations$p.val), ]
    associations <- head(associations[!is.na(associations$p.val), ], n)
    if (nrow(associations) == 0) stop("No features with a valid p-value to plot")
    missing.feats <- setdiff(rownames(associations), rownames(feat))
    if (length(missing.feats) > 0) {
        stop(sprintf(
            "Features not found in the %s feature matrix: %s",
            feature.type, paste(missing.feats, collapse = ", ")
        ))
    }

    feat.names <- rownames(associations)
    if (!is.null(feat.renamer)) {
        if (!is.function(feat.renamer)) stop("feat.renamer must be a function")
        feat.names <- feat.renamer(feat.names)
        if (!is.character(feat.names) || length(feat.names) != nrow(associations)) {
            stop("feat.renamer must return a character vector of the same length as its input")
        }
    }

    label <- label(siamcat)
    samples <- intersect(names(label$label), colnames(feat))

    if (label$type == "BINARY") {
        eff <- if (!is.null(associations$rank.biserial)) associations$rank.biserial else associations$fc
        eff.name <- if (!is.null(associations$rank.biserial)) "r[rb]" else "log[10]~FC"
    } else {
        eff <- associations$spearman
        eff.name <- "rho"
    }
    # plotmath facet labels: feature name on the first line, stats on the second
    facets <- sprintf(
        "atop(%s, %s == %s ~ ',' ~ italic(q) == %s)",
        deparse_strings(feat.names), eff.name,
        deparse_strings(sprintf("%.2f", eff)),
        deparse_strings(sprintf("%.2g", associations$p.adj))
    )

    df <- data.frame(
        facet = factor(
            rep(facets, each = length(samples)), levels = facets
        ),
        value = as.vector(t(feat[rownames(associations), samples, drop = FALSE])),
        label = rep(label$label[samples], times = nrow(associations))
    )

    ylab <- if (feature.type == "normalized") "Normalized abundance" else "Abundance"

    if (label$type == "BINARY") {
        df$label <- factor(
            names(label$info)[match(df$label, label$info)],
            levels = names(sort(label$info))
        )
        col <- pos_neg_neut_palette[c(2, 1)]
        plot <- ggplot(df, aes(x = label, y = value, fill = label)) +
            geom_boxplot(outlier.shape = NA, alpha = 0.5) +
            geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 1) +
            scale_fill_manual(values = col, guide = "none") +
            labs(x = NULL, y = ylab)
    } else {
        plot <- ggplot(df, aes(x = label, y = value)) +
            geom_point(alpha = 0.6) +
            geom_smooth(method = "lm", formula = y ~ x, color = pos_neg_neut_palette[1]) +
            labs(x = label$name, y = ylab)
    }

    plot <- plot +
        facet_wrap(~facet, scales = "free_y", ncol = ncol, labeller = label_parsed) +
        theme_siamcat(font.size)

    if (!isFALSE(title)) {
        if (!is.character(title)) {
            title <- label$name
        }
        plot <- plot + ggtitle(title)
    }

    if (!is.null(fn.plot)) {
        ggsave(
            fn.plot, plot, bg = "white", height = height, width = width
        )
    }

    return(plot)
}

# quote strings so that they are safe to use in plotmath expressions
deparse_strings <- function(x) {
    vapply(x, function(s) deparse(s), character(1), USE.NAMES = FALSE)
}
