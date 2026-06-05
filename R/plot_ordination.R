#' @title Plot ordination of a SIAMCAT object
#'
#' @description A thin wrapper around \code{\link[phyloseq]{plot_ordination}}
#' and \code{\link[phyloseq]{ordinate}} with sensible defaults and boilerplate
#' code for SIAMCAT objects. Computes and plots an ordination of the
#' phyloseq object contained within a \link{siamcat-class} object.
#'
#' @usage plot.ordination(siamcat, method = "PCoA", distance = "bray",
#' color.by = NULL, name.color.by = NULL, palette = NULL,
#' font.size = 14, fn.plot = NULL)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param method string, ordination method passed to
#' \code{\link[phyloseq]{ordinate}}. Supported methods include
#' \code{c("PCoA", "MDS", "NMDS", "DPCoA", "CAP", "RDA", "CCA", "DCA")},
#' defaults to \code{"PCoA"}
#'
#' @param distance string, distance metric passed to
#' \code{\link[phyloseq]{ordinate}}, defaults to \code{"bray"}
#' (Bray-Curtis dissimilarity)
#'
#' @param color.by string, name of a column in the sample metadata to use
#' for coloring points. If \code{NULL} (default), points are not colored.
#' Numeric columns are colored with a continuous color scale; categorical
#' columns are colored with a discrete palette.
#'
#' @param name.color.by string, label for the color legend. If \code{NULL}
#' (default), the value of \code{color.by} is used.
#'
#' @param palette for continuous \code{color.by}: a ColorBrewer palette name,
#' defaults to \code{"RdBu"}. For categorical \code{color.by}: a vector of
#' valid R colors, defaults to \code{okabe_palette}. If \code{NULL},
#' the appropriate default is used.
#'
#' @param font.size integer, base font size for the plot, defaults to
#' \code{14}
#'
#' @param fn.plot string, filename for the plot (any extension supported by
#' \code{\link[ggplot2]{ggsave}} is allowed). If \code{NULL} (default),
#' the plot is only returned as a ggplot object and not saved to disk.
#'
#' @return Returns the ggplot plot object invisibly
#'
#' @keywords SIAMCAT plot.ordination
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
#' # Example data
#' data(siamcat_example)
#'
#' # Simple PCoA with Bray-Curtis dissimilarity
#' plot.ordination(siamcat_example)
#'
#' # Color by a metadata column
#' plot.ordination(siamcat_example, color.by = "disease")
#'
#' # NMDS with Jaccard distance, colored by a continuous variable
#' plot.ordination(siamcat_example, method = "NMDS", distance = "jaccard",
#'     color.by = "age", name.color.by = "Age (years)")

plot.ordination.siamcat <- function(
    siamcat, method="PCoA", distance="bray", color.by=NULL, name.color.by=NULL, palette=NULL, font.size=14, fn.plot = NULL
) {
    ord <- phyloseq::ordinate(siamcat@phyloseq, method = method, distance = distance)
    xlab <- sprintf("%s1 (%s)", method, distance)
    ylab <- sprintf("%s2 (%s)", method, distance)
    
    meta <- sample_data(siamcat@phyloseq)

    if (!is.null(color.by) && color.by %in% colnames(meta)) {
        if (is.null(name.color.by)) name.color.by <- color.by
        p <- phyloseq::plot_ordination(siamcat@phyloseq, ord, color=color.by) +
            labs(x=xlab, y=ylab, color=name.color.by)
        meta_col <- meta[[color.by]]
        if (is.numeric(meta_col)) {
            if (is.null(palette)) palette <- "RdBu"
            p <- p + scale_color_distiller(palette = palette, labels = label_number(scale = 1e-6, suffix = "M"))
        } else {
            if (is.null(palette)) palette <- okabe_palette
            if (length(unique(meta_col)) > length(okabe_palette)) {
                stop("Number of groups in color.by exceeds the number of colors in the palette.")
            }
            p <- p + scale_color_manual(values = palette)
        }
    } else if (!is.null(color.by)) {
        stop("color.by column not found in sample data.")
    } else {
        p <- phyloseq::plot_ordination(siamcat@phyloseq, ord) +
                labs(x=xlab, y=ylab) 
    }

    p <- p + theme_siamcat(font.size)

    if (!is.null(fn.plot)) {
         ggsave(fn.plot, p, bg="white", width=7, height=6)
    }

    return(p)
}