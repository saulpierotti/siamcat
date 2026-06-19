
#' @title Generate ordination of a SIAMCAT object
#'
#' @description A thin wrapper around 
#' \code{\link[phyloseq]{ordinate}}.
#' Computes and stores an ordination of the
#' phyloseq object contained within a \link{siamcat-class} object
#' in the ordination(siamcat) slot.
#' 
#' @usage make.ordination(siamcat, distance="bray", method="Pcoa")
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
#' @param feature.type string, on which type of features should the function
#' work? Can be either \code{"original"}, \code{"filtered"}.
#' Please only change this paramter if you know what
#' you are doing!

make.ordination <- function(siamcat, distance="bray", method="PCoA", feature.type="filtered"){
    # get the right features
    if (feature.type == 'original'){
        feat <- get.orig_feat.matrix(siamcat)
        if (verbose > 1) message('+ using original features')
    } else if (feature.type == 'filtered'){
        if (is.null(filt_feat(siamcat, verbose=0))){
            stop('Features have not yet been filtered, exiting...\n')
        }
        feat <- get.filt_feat.matrix(siamcat)
    } else if (feature.type == 'normalized'){
        stop("Normalised features are not allowed for ordination.")
    }
    temp_phyloseq <- phyloseq(otu_table=otu_table(feat, taxa_are_rows=TRUE))
    ordination(siamcat) <- list(
        ord = phyloseq::ordinate(temp_phyloseq, method = method, distance = distance),
        distance = distance,
        method = method
    )
    return(siamcat)
}