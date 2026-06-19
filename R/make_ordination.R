
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

make.ordination <- function(siamcat, distance="bray", method="PCoA"){
    ordination(siamcat) <- list(
        ord = phyloseq::ordinate(siamcat@phyloseq, method = method, distance = distance),
        distance = distance,
        method = method
    )
    return(siamcat)
}
    
