
check.confounders.ggplot <- function(siamcat, meta.vars=NULL,
    feature.type='normalized', verbose=1, font.size=14) {

    if (verbose > 1) message("+ starting check.confounders")
    s.time <- proc.time()[3]
    label <- label(siamcat)
    if (label$type!='BINARY'){
        stop("Confounder check is currently only possible for",
            " classification tasks")
    }
    meta <- meta(siamcat)

    if (!is.null(meta.vars)) {
        if (!(all(meta.vars %in% colnames(meta)))) {
            stop("meta.vars contains invalid names which are not columns of the metadata table.")
        }
        meta <- meta[,colnames(meta) %in% meta.vars]
    }

    # get features
    if (feature.type == 'original'){
        feat <- get.orig_feat.matrix(siamcat)
    } else if (feature.type == 'filtered'){
        if (is.null(filt_feat(siamcat, verbose=0))){
            stop('Features have not yet been filtered, exiting...\n')
        }
        feat <- get.filt_feat.matrix(siamcat)
    } else if (feature.type == 'normalized'){
        if (is.null(norm_feat(siamcat, verbose=0))){
            stop('Features have not yet been normalized, exiting...\n')
        }
        feat <- get.norm_feat.matrix(siamcat)
    }
    if (is.null(meta)) {
        stop('SIAMCAT object does not contain any metadata.\nExiting...')
    }
    meta <- factorize.metadata(meta, verbose) # creates data.frame

    # remove nested variables
    indep <- vapply(colnames(meta), FUN=function(x) {
        return(condentropy(label$label, discretize(meta[,x])))},
        FUN.VALUE = numeric(1))
    if ((verbose >= 1) & (length(names(which(indep == 0))) > 0)){
        message("++ metadata variables:\n\t",
                paste(c(names(which(indep == 0))), collapse = " & "),
                "\n++ are nested inside the label and ",
                "have been removed from this analysis")
    }
    meta <- meta[,names(which(indep != 0)), drop=FALSE]

    # remove metavariables with less than 2 levels
    n.levels <- vapply(meta,
        FUN = function(x){length(unique(x))},
        FUN.VALUE = integer(1))
    if (any(n.levels < 2)){
        s.name <- names(which(n.levels < 2))
        if (verbose >= 1){
            message("++ remove metadata variables, since all ",
                "subjects have the same value\n\t", s.name)
        }
        meta <- meta[,which(n.levels > 1), drop=FALSE]
    }

    if (ncol(meta) > 10){
        msg <- paste0("The maximum recommended number of metadata variables is 10.\n",
            "Please be aware that some visualizations may not work.")
        message(msg)
    }

    meta <- cleanup.colnames(meta)

    # THIRD PLOT(S) - original confounder check descriptive stat plots
    confounders.descriptive.plots.ggplot(meta, label, font.size, verbose)

    e.time <- proc.time()[3]
    if (verbose > 1) {
        msg <- paste("+ finished check.confounders in", 
            formatC(e.time - s.time, digits = 3), "s")
        message(msg)
    }
    if (verbose == 1) {
        msg <- paste("Finished checking metadata for confounders,",
            "results plotted to:", fn.plot)
        message(msg)
    }
}

#'@keywords internal
confounders.descriptive.plots.ggplot <- function(meta, label, font.size, verbose) {
    cases <- which(label$label == max(label$info))
    controls <- which(label$label == min(label$info))
    p.lab <- names(which(label$info == max(label$info)))
    n.lab <- names(which(label$info == min(label$info)))
    
    colnames(meta) <- gsub("[_.-]", " ", colnames(meta))
    colnames(meta) <- paste(
        toupper(substring(colnames(meta), 1, 1)),
        substring(colnames(meta), 2),
        sep = ""
    )
    
    colors <- brewer.pal(6, "Spectral")
    histcolors <- brewer.pal(9, "YlGnBu")
    label.vec <- names(label$info)[match(label$label, label$info)]

    for (mname in colnames(meta)) {
        tmp <- data.frame(
            mvar=meta[[mname]],
            label=label.vec
        )
        if (verbose > 1){
            msg <- paste("+++ checking",mname,"as a potential confounder")
            message(msg)
        }
        u.val <- unique(tmp$mvar)

        if (length(u.val) == 1) {
            if (verbose > 1) {
                message("+++ skipped because all subjects have the",
                    "same value")}
        } else if (length(u.val) <= 6) {
            if (verbose > 1) message("++++ discrete variable, using a bar plot")
            if (verbose > 2) message("++++ plotting barplot")
            
            plot <- ggplot(tmp, aes(x=label, fill=mvar)) +
                geom_bar(position="fill") +
                scale_fill_manual(values=colors) +
                scale_y_continuous(expand = c(0,0)) +
                theme_siamcat_hist(font.size) +
                labs(x="", y="", fill="") +
                ggtitle(mname)
            ggsave("tmp.png", plot)

            layout(matrix(c(1, 1, 2)))

            # barplot
            par(mar = c(4.1, 9.1, 4.1, 9.1))
            vps <- baseViewports()
            pushViewport(vps$figure)
            vp1 <- plotViewport()
            bar.plot <- barplot(freq, ylim = c(0, 1), main = mname,
                                names.arg = names(label$info), col = colors)
            legend(2.5, 1, legend = var.level.names[[m]],
                    xpd = NA, lwd = 2, col = colors,
                    inset = 0.5, bg = "grey96", cex = 0.8)
            ifelse(length(u.val) > 4,
                    p.val <- fisher.test(ct, simulate.p.value = TRUE,
                                        B = 2000)$p.value,
                    p.val <- fisher.test(ct)$p.value)
            mtext(
                paste("Fisher Test P Value:", format(p.val, digits = 4)),
                cex = 0.8, side = 1, line = 2)
            popViewport()

            if (verbose > 2)
                message("++++ drawing contingency table")

            # contingency table
            plot.new()
            vps <- baseViewports()
            pushViewport(vps$figure)
            niceLabel <- factor(label$label, levels = label$info,
                                labels = names(label$info))
            vp1 <- plotViewport()
            t <- addmargins(table(mvar, niceLabel, dnn = c(mname, "Label")))
            grid.table(t, theme = ttheme_minimal())
            popViewport()
            par(mfrow = c(1, 1), bty = "o")
        } else {
            if (verbose > 1)
                message("++++ continuous variable, using a Q-Q plot")

            layout(rbind(c(1, 2), c(3, 4)))

            if (verbose > 2)
                message("++++ panel 1/4: Q-Q plot")
            par(mar = c(4.5, 4.5, 2.5, 1.5), mgp = c(2.5, 1, 0))
            ax.int <- c(min(mvar, na.rm = TRUE), max(mvar, na.rm = TRUE))
            qqplot(mvar[controls], mvar[cases], xlim = ax.int,
                    ylim = ax.int, pch = 16, cex = 0.6, xlab = n.lab,
                    ylab = p.lab, main = paste("Q-Q plot for", mname))
            abline(0, 1, lty = 3)
            p.val <- wilcox.test(mvar[controls], mvar[cases],
                                exact = FALSE)$p.value
            text(ax.int[1] + 0.9 * (ax.int[2] - ax.int[1]),
                ax.int[1] + 0.1 * (ax.int[2] - ax.int[1]),
                cex = 0.8,paste("MWW Test P Value:",
                                format(p.val, digits = 4)),
                pos = 2)

            if (verbose > 2)
                message("++++ panel 2/4: X histogram")
            par(mar = c(4, 2.5, 3.5, 1.5))
            hist(mvar[controls], main = n.lab, xlab = mname,
                col = histcolors, breaks = seq(min(mvar, na.rm = TRUE),
                                                max(mvar, na.rm = TRUE),
                                                length.out = 10))
            mtext(paste("N =", length(mvar[controls])), cex = 0.6,
                    side = 3, adj = 1, line = 1)

            if (verbose > 2)
                message("++++ panel 3/4: X boxplot")
            par(mar = c(2.5, 4.5, 2.5, 1.5))
            boxplot(mvar ~ label$label, range = 1.5,
                    use.cols = TRUE, names = names(label$info),
                    ylab = mname, main = paste("Boxplot for", mname),
                    col = histcolors, outpch = NA)
            stripchart(mvar ~ label$label, vertical = TRUE, add = TRUE,
                    method = "jitter", pch = 20)

            if (verbose > 2)
                message("++++ panel 4/4: Y histogram")
            par(mar = c(4.5, 2.5, 3.5, 1.5))
            hist(mvar[cases], main = p.lab, xlab = mname,
                    col = histcolors, breaks = seq(min(mvar, na.rm = TRUE),
                                                max(mvar, na.rm = TRUE),
                                                length.out = 10))
            mtext(paste("N =", length(mvar[cases])), cex = 0.6,
                    side = 3, adj = 1, line = 1)
            par(mfrow = c(1, 1))
        }
            }
}

#'@keywords internal
cleanup.colnames <- function(meta) {
    return(meta)
}

#'@keywords internal
confounders.continuous.var.plot <- function(mvar, label, verbose) {
    u.val <- sort(unique(mvar)[!is.na(unique(mvar))])
    colors <- brewer.pal(6, "Spectral")
    histcolors <- brewer.pal(9, "YlGnBu")

    if (length(u.val) == 1) {
        if (verbose > 1) {
            message("+++ skipped because all subjects have the",
                "same value")}
    } else if (length(u.val) <= 6) {
        if (verbose > 1) message("++++ discrete variable, using a bar plot")

        # create contingency table
        ct <- vapply(u.val, FUN = function(x) {
            # deal with NAs...?
            return(c(length(intersect(which(mvar == x), controls)),
                    length(intersect(which(mvar == x), cases))))},
            USE.NAMES = FALSE,
            FUN.VALUE = integer(2))

        freq <- t(ct / rowSums(ct))
        mvar <- factor(mvar, levels = sort(unique(na.omit(mvar))),
                        labels = var.level.names[[m]])

        if (verbose > 2)
            message("++++ plotting barplot")

        layout(matrix(c(1, 1, 2)))

        # barplot
        par(mar = c(4.1, 9.1, 4.1, 9.1))
        vps <- baseViewports()
        pushViewport(vps$figure)
        vp1 <- plotViewport()
        bar.plot <- barplot(freq, ylim = c(0, 1), main = mname,
                            names.arg = names(label$info), col = colors)
        legend(2.5, 1, legend = var.level.names[[m]],
                xpd = NA, lwd = 2, col = colors,
                inset = 0.5, bg = "grey96", cex = 0.8)
        ifelse(length(u.val) > 4,
                p.val <- fisher.test(ct, simulate.p.value = TRUE,
                                    B = 2000)$p.value,
                p.val <- fisher.test(ct)$p.value)
        mtext(
            paste("Fisher Test P Value:", format(p.val, digits = 4)),
            cex = 0.8, side = 1, line = 2)
        popViewport()

        if (verbose > 2)
            message("++++ drawing contingency table")

        # contingency table
        plot.new()
        vps <- baseViewports()
        pushViewport(vps$figure)
        niceLabel <- factor(label$label, levels = label$info,
                            labels = names(label$info))
        vp1 <- plotViewport()
        t <- addmargins(table(mvar, niceLabel, dnn = c(mname, "Label")))
        grid.table(t, theme = ttheme_minimal())
        popViewport()
        par(mfrow = c(1, 1), bty = "o")
    } else {
        if (verbose > 1)
            message("++++ continuous variable, using a Q-Q plot")

        layout(rbind(c(1, 2), c(3, 4)))

        if (verbose > 2)
            message("++++ panel 1/4: Q-Q plot")
        par(mar = c(4.5, 4.5, 2.5, 1.5), mgp = c(2.5, 1, 0))
        ax.int <- c(min(mvar, na.rm = TRUE), max(mvar, na.rm = TRUE))
        qqplot(mvar[controls], mvar[cases], xlim = ax.int,
                ylim = ax.int, pch = 16, cex = 0.6, xlab = n.lab,
                ylab = p.lab, main = paste("Q-Q plot for", mname))
        abline(0, 1, lty = 3)
        p.val <- wilcox.test(mvar[controls], mvar[cases],
                            exact = FALSE)$p.value
        text(ax.int[1] + 0.9 * (ax.int[2] - ax.int[1]),
            ax.int[1] + 0.1 * (ax.int[2] - ax.int[1]),
            cex = 0.8,paste("MWW Test P Value:",
                            format(p.val, digits = 4)),
            pos = 2)

        if (verbose > 2)
            message("++++ panel 2/4: X histogram")
        par(mar = c(4, 2.5, 3.5, 1.5))
        hist(mvar[controls], main = n.lab, xlab = mname,
            col = histcolors, breaks = seq(min(mvar, na.rm = TRUE),
                                            max(mvar, na.rm = TRUE),
                                            length.out = 10))
        mtext(paste("N =", length(mvar[controls])), cex = 0.6,
                side = 3, adj = 1, line = 1)

        if (verbose > 2)
            message("++++ panel 3/4: X boxplot")
        par(mar = c(2.5, 4.5, 2.5, 1.5))
        boxplot(mvar ~ label$label, range = 1.5,
                use.cols = TRUE, names = names(label$info),
                ylab = mname, main = paste("Boxplot for", mname),
                col = histcolors, outpch = NA)
        stripchart(mvar ~ label$label, vertical = TRUE, add = TRUE,
                method = "jitter", pch = 20)

        if (verbose > 2)
            message("++++ panel 4/4: Y histogram")
        par(mar = c(4.5, 2.5, 3.5, 1.5))
        hist(mvar[cases], main = p.lab, xlab = mname,
                col = histcolors, breaks = seq(min(mvar, na.rm = TRUE),
                                            max(mvar, na.rm = TRUE),
                                            length.out = 10))
        mtext(paste("N =", length(mvar[cases])), cex = 0.6,
                side = 3, adj = 1, line = 1)
        par(mfrow = c(1, 1))
    }
}

#'@keywords internal
get.names <- function(meta) {
    temp <- lapply(meta, FUN=function(x) {
        if (is.numeric(x) | is.character(x)) {
            if (length(unique(x)) <= 6) {return(levels(as.factor(x)))}
            else {return(NULL)}}
        else if (is.factor(x)) {return(levels(x))}
        else {return(levels(x))}})
}

#'@keywords internal
factorize.metadata <- function(meta, verbose) {

    if ('BMI' %in% toupper(colnames(meta))) {
        idx <- match('BMI', toupper(colnames(meta)))
        meta[,idx] <- factorize.bmi(meta[,idx])}

    factorized <- as.data.frame(lapply(meta, FUN=function(x) {
        if (is.numeric(x) & (length(unique(x)) > 5)) {
            quart <- quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE)
            temp <- cut(x, unique(quart), include.lowest = TRUE)
            return(factor(temp, labels = seq_along(levels(temp))))}
        else {(return(as.factor(x)))}}))
    rownames(factorized) <- rownames(meta)

    # check for IDs and other metavariable with too many levels
    n.levels <- vapply(colnames(factorized), FUN=function(x){
        length(levels(factorized[[x]]))}, FUN.VALUE = integer(1))
    if (any(n.levels > 0.9*nrow(meta))){
        remove.meta <- names(which(n.levels > 0.9*nrow(meta)))
        if (verbose > 1){
            message("++ metadata variables:\n\t",
                paste(remove.meta, collapse = " & "),
                "\n++ have too many levels and ",
                "have been removed from this analysis")
        }
        factorized <- factorized[,-which(colnames(factorized) %in% remove.meta)]
    }
    return(factorized)
}

#'@keywords internal
factorize.bmi <- function(bmi) {
    # ranges taken from CDC
    # https://www.cdc.gov/healthyweight/assessing/bmi/adult_bmi/index.html

    if (!is.matrix(bmi)) bmi <- as.matrix(bmi)
    temp <- vapply(bmi, FUN=function(x) {
        if (is.na(x)) {return(as.character(NA))}
        else if (x < 18.5) {return("Underweight")}
        else if ((x >= 18.5) & (x <= 24.9)) {return("Healthy")}
        else if ((x > 24.9) & (x <= 29.9)) {return("Overweight")}
        else if (x > 29.9) {return("Obese")}},
        FUN.VALUE = character(1), USE.NAMES = TRUE)
    #names(temp) <- rownames(bmi)
    return(as.factor(temp))
}
