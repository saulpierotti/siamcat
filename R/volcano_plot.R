#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Visualize associations between features and classes as volcano plot
#'
#' @description This function creates a volcano plot to visualize the
#' association between features and the label
#'
#' @usage volcano.plot(siamcat)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param fn.plot string, filename for the plot (any extension supported by
#' ggsave is allowed). If \code{fn.plot} is \code{NULL}, the plot will only
#' be returned as a ggplot object.
#'
#' @param color.scheme valid R color scheme or vector of valid R colors (must
#' be of length 3 for positive, negative, and non-significant associations),
#' defaults to \code{c('red', 'blue', 'gray')}
#'
#' @param annotate integer, number of features to annotate with the name
#'
#' @param annot.size integer, size of annotation text
#'
#' @param annot.y.exp float, percent vertical expansion of the plot area
#' to accommodate annotations
#'
#' @param annot.y.shift float, upward shift of the annotations from the
#' respective datapoints
#'
#' @param font.size integer, base font size for the plot
#'
#' @return Returns the ggplot plot object
#'
#' @keywords SIAMCAT volcano.plot
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
#' # Example data
#' data(siamcat_example)
#'
#' # Simple example
#' volcano.plot(siamcat_example, fn.plot = "./volcano.pdf")
volcano.plot <- function(
    siamcat, fn.plot = NULL, color.scheme = c("red", "blue", "gray"),
    annotate = 3,
    annot.size = 4, annot.y.exp = 0.2, annot.y.shift = NULL,
    font.size = 14
) {
    associations <- associations(siamcat, verbose = 0)
    if (is.null(associations)) {
        stop(
            "SIAMCAT object does not contain association testing results! ",
            "Exiting..."
        )
    }
    associations$label <- rownames(associations)
    assoc.param <- assoc_param(siamcat)

    if (length(color.scheme) == 1) {
        tryCatch(
            {
                col <- brewer.pal(3, color.scheme)[c(1, 3, 2)]
            },
            error = function(e) {
                stop(
                    "color.scheme contains 1 element, so it is interpreted ",
                    "as a ColorBrewer palette, but the palette name is invalid."
                )
            }
        )
    } else if (length(color.scheme) == 3) {
        tryCatch(col2rgb(color.scheme), error = function(e) {
            stop(
                "color.scheme contains 3 elements, so it is interpreted ",
                "as containing individual colors, but the color names are",
                " invalid."
            )
        })
        col <- color.scheme
    } else {
        stop(
            "color.scheme must contain 3 R color values or",
            " the name of 1 ColorBrewer palette."
        )
    }

    ns.label <- "n. s."
    if (label(siamcat)$type == "BINARY") {
        xlab <- bquote(log[10] ~ "fold change")
        fill_lab <- "Enrichment"
        control.label <- names(label(siamcat)$info)[label(siamcat)$info == -1]
        case.label <- names(label(siamcat)$info)[label(siamcat)$info == 1]
        names(col) <- c(case.label, control.label, ns.label)
        associations$class <- ifelse(
            associations$p.adj < assoc.param$alpha,
            ifelse(associations$fc > 0, case.label, control.label),
            ns.label
        )
        associations$eff <- associations$fc
    }

    if (label(siamcat)$type == "CONTINUOUS") {
        xlab <- "Effect size"
        fill_lab <- "Effect direction"
        names(col) <- c("Positive", "Negative", "n. s.")
        associations$class <- ifelse(
            associations$p.adj < assoc.param$alpha,
            ifelse(associations$fc > 0, case.label, control.label),
            ns.label
        )
        associations$eff <- associations$beta
    }

    associations.to.label <- associations[
        associations$p.adj < assoc.param$alpha,
    ]
    associations.to.label <- associations.to.label[
        order(associations.to.label$p.adj),
    ]
    associations.to.label <- rbind(
        head(associations.to.label[associations.to.label$eff > 0, ], annotate),
        head(associations.to.label[associations.to.label$eff < 0, ], annotate)
    )

    mult.corr.label <- c(
        "holm" = "q",
        "hochberg" = "q",
        "hommel" = "q",
        "bonferroni" = "q",
        "BH" = "q",
        "BY" = "q",
        "fdr" = "q",
        "none" = "p"
    )[[assoc.param$mult.corr]]

    line.size <- 0.5
    half.line <- font.size / 2
    rel.small <- 12 / 14
    small.size <- rel.small * font.size

    plot <- ggplot(
        associations,
        aes(
            x = eff, y = -log10(p.adj),
            size = pr.all, fill = class, label = label
        )
    ) +
        geom_point(shape = 21, alpha = 0.3) +
        ggrepel::geom_text_repel(
            data = associations.to.label,
            ylim = c(
                max(-log10(associations$p.adj)),
                max(-log10(associations$p.adj)) * (1 + annot.y.exp)
            ),
            show.legend = FALSE, size = annot.size, force = 20, force_pull = 0
        ) +
        geom_hline(
            yintercept = -log10(assoc.param$alpha),
            color = "gray", lty = "dashed", lwd = 0.5
        ) +
        annotate(
            "text",
            x = Inf, y = -log10(assoc.param$alpha),
            label = deparse(bquote(alpha ~ "=" ~ .(assoc.param$alpha))),
            color = "gray40", parse = TRUE, hjust = 1, vjust = -1
        ) +
        scale_fill_manual(
            values = col, breaks = names(col),
            guide = guide_legend(override.aes = list(size = 6))
        ) +
        scale_size(guide = guide_legend(reverse = TRUE)) +
        labs(
            x = xlab, y = bquote(-log[10](italic(.(mult.corr.label)))),
            size = "Prevalence", fill = fill_lab
        ) +
        theme_grey(font.size) +
        theme(
            line = element_line(
                color = "black", linetype = 1,
                lineend = "butt", linewidth = line.size
            ),
            text = element_text(
                size = font.size, hjust = 0.5, vjust = 0.5, angle = 0,
                lineheight = .9, margin = margin(), debug = FALSE
            ),
            axis.line = element_blank(),
            axis.text = element_text(color = "black", size = small.size),
            axis.text.x = element_text(
                margin = margin(t = small.size / 4), vjust = 1
            ),
            axis.text.y = element_text(
                margin = margin(r = small.size / 4), hjust = 1
            ),
            axis.ticks = element_line(color = "black", linewidth = line.size),
            axis.ticks.length = unit(half.line / 2, "pt"),
            axis.title.x = element_text(
                margin = margin(t = half.line / 2), vjust = 1
            ),
            axis.title.y = element_text(
                angle = 90, margin = margin(r = half.line / 2), vjust = 1
            ),
            panel.background = element_blank(),
            panel.border = element_rect(color = "black", linewidth = line.size),
            panel.grid = element_blank(),
            legend.background = element_blank(),
            legend.spacing = unit(font.size, "pt"),
            legend.margin = margin(0, 0, 0, 0),
            legend.key = element_blank(),
            legend.key.size = unit(1.1 * font.size, "pt"),
            legend.text = element_text(size = rel(rel.small)),
            legend.title = element_text(hjust = 0),
            legend.box.background = element_blank(),
            legend.box.spacing = unit(font.size, "pt"),
            plot.margin = margin(half.line, half.line, half.line, half.line),
            complete = TRUE
        )

    if (nrow(associations.to.label) != 0) {
        plot <- plot + scale_y_continuous(
            expand = expansion(mult = c(0.05, annot.y.exp))
        )
    }

    # save the plot
    if (!is.null(fn.plot)) {
        ggsave(
            fn.plot, plot, bg = "white", height = 6, width = 7, device = "png"
        )
    }

    return(plot)
}
