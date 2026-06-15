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
#' stratify = NULL, inseparable = NULL, verbose = 1)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param num.folds integer number of cross-validation folds (needs to be 
#' \code{>=2}), defaults to \code{2}. Set to \code{Inf} to perform
#' leave-one-(group)-out CV; in that case \code{stratify} and 
#' \code{num.resample} are ignored and stored as \code{NA}.
#'
#' @param num.resample integer, number of resampling rounds.
#' Set to \code{1} to perform a single CV run 
#' (no resampling). Ignored when leave-one-out CV is performed.
#' Defaults to \code{1}.
#'
#' @param stratify boolean, should the splits be stratified so that an equal 
#' proportion of classes are present in each fold? Ignored for regression 
#' tasks and when leave-one-(group)-out CV is performed. Defaults to 
#' \code{NULL}, which enables stratification only for classification
#' tasks.
#'
#' @param inseparable string or integer, name or column index of a metadata 
#' variable whose values should be kept within the same fold (e.g. to avoid 
#' splitting repeated measurements from the same individual across folds).
#' Defaults to \code{NULL}. See Details below.
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
#' The function stores the cross-validation configuration in the 
#' \code{data_split}-slot of the \link{siamcat-class} object as a list with 
#' the following entries: \itemize{
#' \item \code{task} - the \code{mlr3} task object used for resampling
#' \item \code{resampling} - the \code{mlr3} resampling object encoding the 
#' CV strategy
#' \item \code{num.folds} - the number of cross-validation folds (\code{NA} 
#' for leave-one-out CV)
#' \item \code{num.resample} - the number of CV repetitions (\code{NA} for 
#' leave-one-out CV)
#' \item \code{loo} - logical, whether leave-one-(group)-out CV is performed
#' \item \code{stratify} - logical, whether stratification was applied
#' \item \code{inseparable} - the name of the inseparable metadata column, 
#' or \code{NULL} if not used }
#'
#' If \code{num.folds} is set to \code{Inf}, or if \code{num.folds} is 
#' greater than or equal to the number of samples (or groups, when 
#' \code{inseparable} is provided), leave-one-(group)-out CV is performed. 
#' In this case \code{stratify} is set to \code{FALSE} and \code{num.resample}
#' is ignored; both are stored as \code{NA} in the \code{data_split} slot. A 
#' warning is issued if either parameter was set to a non-default value.
#' 
#' If \code{inseparable} is provided, all samples sharing the same value of 
#' that metadata variable are guaranteed to end up in the same fold. This is 
#' useful when the data contains repeated measurements per individual. 
#' The \code{inseparable} column must be categorical and must not be 
#' highly correlated with the label (Pearson |r| <= 0.9).
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
    stratify=NULL, inseparable=NULL, verbose=1
) {

    if (verbose > 1)
        message("+ starting create.data.split")
    s.time <- proc.time()[3]

    label <- label(siamcat)
    labelNum <- as.numeric(label$label)
    names(labelNum) <- names(label$label)

    ##########################################################################
    # check arguments
    ##########################################################################

    # stratify
    if (!is.null(stratify) && !is.logical(stratify)) {
        stop("stratify must be a boolean value (TRUE or FALSE) or NULL")
    }
    if (label$type == 'CONTINUOUS'){
        if (isTRUE(stratify)) {
            stop("'stratify' is not allowed for regression tasks.")
        }
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
        if (is.null(stratify)) stratify <- TRUE
    } else if (label$type == 'TEST'){
        stop("Cannot create data split for TEST object!")
    }
    stopifnot(is.logical(stratify))

    # num.resample
    if(!(is.numeric(num.resample) && num.resample >= 1)) {
        stop("num.resample must be a numeric value >= 1")
    }

    #num.folds
    if(!(is.numeric(num.folds) && num.folds >= 2)) {
        stop("num.folds must be a numeric value >= 2")
    }

    # inseparable
    if (!is.null(inseparable)) {
        if (is.null(meta(siamcat))) {
            stop("Meta-data must be provided if the inseparable parameter is not NULL")
        }
        if (length(inseparable) != 1) {
            stop("Inseparable parameter must be a single column index or name of metadata matrix")
        }
        if (is.numeric(inseparable)) {
            if (inseparable < 1 || inseparable > ncol(meta(siamcat))) {
                stop(
                    "Inseparable parameter is out of bounds.",
                    " It must be a single column index of the metadata matrix."
                )
            }
            inseparable <- colnames(meta(siamcat))[inseparable]
        }
        if (!(is.character(inseparable))) {
            stop(
                "Inseparable parameter must be either a single column index",
                " or a single column name of metadata matrix"
            )
        }
        if(!(inseparable %in% colnames(meta(siamcat)))) {
            stop(
                "Inseparable parameter does not match any column name of the metadata matrix."
            )
        }
    }

    ##########################################################################
    
    # generate data (without features)
    if (!is.null(meta(siamcat))){
        if(!all(names(labelNum) == rownames(meta(siamcat)))) {
            stop(
                "The names of the label vector and the rownames of the metadata matrix do not match.",
                " Please raise an issue with the developers."
            )
        }
        if ("labelNum" %in% colnames(meta(siamcat))) {
            stop(
                "The name 'labelNum' is reserved for the internal data split and cannot be",
                " used as a column name in the metadata. Please rename this column and try again."
            )
        }
        data <- cbind(
            data.frame(labelNum=labelNum),
            meta(siamcat)[names(labelNum),]
        )
    } else {
        data <- data.frame(labelNum=labelNum)
    }
    rownames(data) <- names(labelNum)
    
    # generate task
    if (label$type == 'BINARY') {
        tsk <- as_task_classif(data, target="labelNum", id="siamcat", positive=label$info[2])
    } else if (label$type == 'CONTINUOUS') {
        tsk <- as_task_regr(data, target="labelNum", id="siamcat")
    } else {
        stop("Unknown label type: ", label$type)
    }

    # set inseparable correctly
    if (!is.null(inseparable)) {
        inseparable_val <- meta(siamcat)[names(labelNum), inseparable]
        if (is.numeric(inseparable_val)) {
            inseparable_num <- inseparable_val
        } else if (is.character(inseparable_val) || is.factor(inseparable_val)) {
            inseparable_num <- as.numeric(as.factor(inseparable_val))
        } else {
            stop(
                "The 'inseparable' column must be either numeric, character or factor.",
                " Inseparable column: ", inseparable, " of class: ", class(inseparable_val),
                ". Please choose a different column or remove the inseparable parameter.",
                " Aborting."
            )
        }
        if (length(unique(inseparable_num)) > length(labelNum)*0.9) {
            stop(
                "Too many unique values in the 'inseparable' column. Did you select a continuous variable?"
            )
        }
        cor_val <- cor(inseparable_num, labelNum)
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

    # if num.folds is bigger than number of samples or groups,
    # reset to leave-one-out CV
    if (num.folds == Inf) {
        loo <- TRUE
    } else if (num.folds >= length(labelNum) && is.null(inseparable)) {
        if (verbose > 1)
            message(
                "+++ Performing unstratified leave-one-out (LOO) cross-validation."
            )
        loo <- TRUE
    } else if (!is.null(inseparable) && num.folds >= length(unique(inseparable_num))) {
        if (verbose > 1)
            message(
                "+++ Performing unstratified leave-one-group-out (LOGO) cross-validation."
            )
        loo <- TRUE
    } else {
        loo <- FALSE
    }

    # reset num.resample and num.folds if doing LOO
    if (loo) {
        if (num.resample != 1 || stratify) {
            warning(
                "Performing leave-one-(group)-out (LOO) cross-validation. Ignoring stratification and num.resample parameters."
            )
        }
        num.folds <- NA
        num.resample <- NA
        stratify <- FALSE
    }

    # If stratify is TRUE for classification, make sure that num.folds does not exceed the
    # maximum number of examples for the class with
    # the fewest training examples.
    if (stratify) {
        if (label$type == 'BINARY' && any(as.data.frame(table(label))[, 2] < num.folds)) {
            stop(
                "+++ Number of CV folds is too large for this data set to
                maintain stratification. Reduce num.folds or turn
                stratification off. Exiting."
            )
        }
        if (label$type == 'CONTINUOUS') {
            stop(
                "Stratified CV is not supported for regression tasks."
            )
        }
        tsk$set_col_roles("labelNum", c("target", "stratum"))
    }

    if (loo) {
        the_rsmp <- rsmp("loo")
    } else {
        if(!(num.folds >= 2 && num.folds < length(labelNum))){
            stop("num.folds has illegal value:", num.folds, ", please raise an issue with the developers.")
        }
        if (!is.null(inseparable) && num.folds >= length(unique(inseparable_num))) {
            stop(
                "num.folds has illegal value for grouped CV:",
                num.folds, ", please raise an issue with the developers."
            )
        }
        if (num.resample > 1) {
            the_rsmp <- rsmp("repeated_cv", folds=num.folds, repeats=num.resample)
        } else if (num.resample == 1) {
            the_rsmp <- rsmp("cv", folds=num.folds)
        } else {
            stop("num.resample must be >= 1. Please raise an issue with the developers.")
        }
    }

    data_split(siamcat) <- list(
        task = tsk,
        resampling = the_rsmp,
        num.resample = num.resample,
        num.folds = num.folds,
        loo = loo,
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