#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Split a dataset into training and a test sets.
#'
#' @name create.data.split
#'
#' @description This function prepares the cross-validation by splitting the 
#' data into \code{num.folds} training and test folds for 
#' \code{num.resample} times.
#'
#' @usage create.data.split(siamcat, num.folds = 2, num.resample = 1, 
#' stratify = TRUE, inseparable = NULL, verbose = 1)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param num.folds integer number of cross-validation folds (needs to be 
#' \code{>=2}), defaults to \code{2}, set to Inf to perform
#' leave-one-(group)-out CV
#'
#' @param num.resample integer, resampling rounds (values \code{<= 1} 
#' deactivate resampling), defaults to \code{1}
#'
#' @param stratify boolean, should the splits be stratified so that an equal 
#' proportion of classes are present in each fold?, will be ignored for 
#' regression tasks, defaults to \code{TRUE}
#'
#' @param inseparable string, name of metadata variable to be inseparable,
#' defaults to \code{NULL}, see Details below
#'
#' @param verbose integer, control output: \code{0} for no output at all, 
#' \code{1} for only information about progress and success, \code{2} for
#' normal level of information and \code{3} for full debug information,
#' defaults to \code{1}
#'
#' @keywords SIAMCAT create.data.split
#'
#' @return object of class \link{siamcat-class} with the \code{data_split}-slot
#' filled
#'
#' @details This function splits the labels within a \link{siamcat-class} 
#' object and prepares the internal cross-validation for the model training 
#' (see \link{train.model}). 
#' 
#' The function saves the training and test instances for the different 
#' cross-validation folds within a list in the \code{data_split}-slot of the 
#' \link{siamcat-class} object, which is a list with four entries: \itemize{ 
#' \item \code{num.folds} - the number of cross-validation folds
#' \item \code{num.resample} - the number of repetitions for the 
#' cross-validation
#' \item \code{training.folds} - a list containing the indices for the 
#' training instances
#' \item \code{test.folds} - a list containing the indices for the 
#' test instances }
#'
#' If provided, the data split will take into account a metadata variable
#' for the data split (by providing the \code{inseparable} argument). For
#' example, if the data contains several samples for the same individual,
#' it makes sense to keep data from the same individual within the
#' same fold.
#' 
#' If \code{inseparable} is given, the \code{stratify} argument will be
#' ignored.
#' 
#' @export
#' 
#' @encoding UTF-8
#'
#' @examples
#' data(siamcat_example)
#'
#' # simple working example
#' siamcat_split <- create.data.split(siamcat_example, num.folds=10, 
#' num.resample=5, stratify=TRUE)
create.data.split <- function(
    siamcat, num.folds=2, num.resample=1,
    stratify=TRUE, inseparable=NULL, verbose=1
) {

    if (verbose > 1)
        message("+ starting create.data.split")
    s.time <- proc.time()[3]

    label <- label(siamcat)
    if (label$type == 'CONTINUOUS'){
        stratify <- FALSE
    } else if (label$type=='BINARY') {
        # check that there are enough samples in each class for the CV to work
        group.numbers <- vapply(label$info,
                                FUN = function(x){
                                    sum(label$label == x)},
                                FUN.VALUE = integer(1))
        if (any(group.numbers <= 5)){
            msg <- paste0("Data set has only:\n",
                paste0(names(group.numbers)[1], "\t", group.numbers[1]),
                "\n",
                paste0(names(group.numbers)[2], "\t", group.numbers[2]),
                "\nThis is not enough for SIAMCAT to proceed!")
            stop(msg)
        }
    } else if (label$type == 'TEST'){
        stop("Cannot create data split for TEST object!")
    }

    labelNum <- as.numeric(label$label)
    names(labelNum) <- names(label$label)

    ### check arguments
    if (num.resample < 1) {
        if (verbose > 1){
            msg <- paste0("+++ Resetting num.resample = 1 (", 
                num.resample, 
                " is an invalid number of resampling rounds)")
            message(msg)
        }
        num.resample <- 1
    }
    if (num.folds < 2) {
        if (verbose > 1){
            msg <- paste0("+++ Resetting num.folds = 2 (", 
                num.folds, " is an invalid number of folds)")
            message(msg)
        }
        num.folds <- 2
    }
    if (!is.null(inseparable) && is.null(meta(siamcat))) {
        stop("Meta-data must be provided if the inseparable parameter is not
            NULL")
    }
    if (!is.null(inseparable)) {
        if (is.numeric(inseparable) && length(inseparable) == 1) {
            stopifnot(inseparable <= ncol(meta(siamcat)))
            stopifnot(inseparable >= 1)
            inseparable <- colnames(meta(siamcat))[inseparable]
        } else if (!(is.character(inseparable) &&
                length(inseparable) == 1)) {
            stop(
                "Inseparable parameter must be either a single column index",
                " or a single column name of metadata matrix"
            )
        }
        stopifnot(inseparable %in% colnames(meta(siamcat)))
    }

    stopifnot(all(names(labelNum) == rownames(meta(siamcat))))
    if ("labelNum" %in% colnames(meta(siamcat))) {
        stop(
            "The name 'labelNum' is reserved for the internal data split and cannot be",
            " used as a column name in the metadata. Please rename this column and try again."
        )
    }
    data <- cbind(
        data.frame(labelNum = labelNum),
        meta(siamcat)[names(labelNum),]
    )
    if (label$type == 'BINARY') {
        tsk <- as_task_classif(data, target="labelNum", id="siamcat", positive=label$info[2])
    } else if (label$type == 'CONTINUOUS') {
        tsk <- as_task_regr(data, target="labelNum", id="siamcat")
    } else {
        stop("Unknown label type: ", label$type)
    }
    stopifnot(!is.null(tsk))

    if (!is.null(inseparable)) {
        stopifnot(is.character(inseparable))
        stopifnot(inseparable %in% colnames(meta(siamcat)))
        iseparableVal <- meta(siamcat)[names(labelNum), inseparable]
        if (is.numeric(iseparableVal)) {
            iseparableNum <- iseparableVal
        } else if (is.character(iseparableVal) || is.factor(iseparableVal)) {
            iseparableNum <- as.numeric(as.factor(meta(siamcat)[names(labelNum), inseparable]))
        } else {
            stop(
                "The 'inseparable' column must be either numeric, character or factor.",
                "Inseparable column: ", inseparable, " of class: ", class(iseparableVal),
                ". Please choose a different column or remove the inseparable parameter.",
                " Aborting."
            )
        }
        if (length(unique(iseparableNum)) > length(labelNum)*0.9) {
            stop(
                "Too many unique values in the 'inseparable' column. Did you select a continuous variable?"
            )
        }
        cor_val <- cor(iseparableNum, labelNum)
        if (abs(cor_val) > 0.9) {
            stop(
                "The 'inseparable' column cannot be highly correlated with the label column.",
                " Inseparable column: ", inseparable, " with correlation: ", cor_val,
                ". Please choose a different column or remove the inseparable parameter.",
                " Aborting."
            )
        }
        tsk$set_col_roles(inseparable, "group")
    }

    # if num.folds is bigger than number of samples,
    # reset to leave-one-out CV and ignore stratification
    if (num.folds == Inf) {
        do_loo <- TRUE
    } else if (num.folds >= length(labelNum) && is.null(inseparable)) {
        if (verbose > 1)
            message(
                "+++ Performing unstratified leave-one-out (LOO) cross-validation."
            )
        do_loo <- TRUE
    } else if (!is.null(inseparable) && num.folds >= length(unique(iseparableNum))) {
        if (verbose > 1)
            message(
                "+++ Performing unstratified leave-one-group-out (LOGO) cross-validation."
            )
        num.folds <- length(unique(iseparableNum))
        do_loo <- TRUE
    } else {
        do_loo <- FALSE
    }
    stopifnot(!is.null(do_loo))
    if (do_loo) {
        if (num.resample != 1 || stratify) {
            warning(
                "Performing leave-one-(group)-out (LOO) cross-validation. Ignoring stratification and num.resample parameters."
            )
        }
        num.folds <- NA
        num.resample <- NA
        stratify <- FALSE
    }

    if (stratify) {
        # If stratify is TRUE for classification, make sure that num.folds does not exceed the
        # maximum number of examples for the class with
        # the fewest training examples.
        if (label$type == 'BINARY' && any(as.data.frame(table(label))[, 2] < num.folds)) {
            stop(
                "+++ Number of CV folds is too large for this data set to
                maintain stratification. Reduce num.folds or turn
                stratification off. Exiting."
            )
        }
        tsk$set_col_roles("labelNum", c("target", "stratum"))
    }

    if (do_loo) {
        the_rsmp <- rsmp("loo")
    } else {
        stopifnot(num.folds >= 2 && num.folds < length(labelNum))
        if (!is.null(inseparable)) {
            stopifnot(num.folds < length(unique(iseparableNum)))
        }
        if (num.resample > 1) {
            the_rsmp <- rsmp("repeated_cv", folds = num.folds, repeats = num.resample)
        } else if (num.resample == 1) {
            the_rsmp <- rsmp("cv", folds = num.folds)
        } else {
            stop("num.resample must be >= 1")
        }
    }
    stopifnot(!is.null(the_rsmp))

    data_split(siamcat) <- list(
        task = tsk,
        resampling = the_rsmp,
        num.resample = num.resample,
        num.folds = num.folds,
        loo = do_loo,
        stratify = stratify,
        inseparable = inseparable
    )

    e.time <- proc.time()[3]
    if (verbose > 1){
        msg <- paste("+ finished create.data.split in",
            formatC(e.time - s.time, digits = 3),"s")
        message(msg)
    }
    if (verbose == 1)
        message("Features split for cross-validation successfully.")
    return(siamcat)
}