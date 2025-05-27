# ============================================================
#   This script contains the necessary functions for analysis.
# ============================================================
clean_probe_ids <- function(ids) unique(sapply(strsplit(ids, "_"), `[`, 1))

cleanMatrixForClusterW <- function(mtx, f_row = 0.5, f_col = 0.5) {
    cat("Cleaning...")
    cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
                f_row, f_col))
    cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
    namtx = is.na(mtx)
    good_row = rowSums(namtx) <= ncol(mtx) * f_row
    good_col = colSums(namtx) <= nrow(mtx) * f_col
    cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
    cleaned_mtx <- mtx[good_row, good_col]
    cat("Cleaning completed.\n")
    return(cleaned_mtx)
}

bSubMostVariable <- function(betas, n=3000) {
    print(paste("Original dimensions:", nrow(betas), "x", ncol(betas)))
    std <- apply(betas, 1, sd, na.rm=TRUE)
    result <- betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
    print(paste("New dimensions:", nrow(result), "x", ncol(result)))
    return(result)
}

# wrapper for genomic neighbors
imputeBetasByGenomicNeighborsMatrix <- function(beta_mat, platform = NULL,
                                                BPPARAM = SerialParam(),
                                                max_neighbors = 3, max_dist = 10000) {
    if (!is.matrix(beta_mat)) {
        stop("Input must be a matrix with probes as rows and samples as columns.")
    }

    imputed_mat <- do.call(cbind, bplapply(seq_len(ncol(beta_mat)), function(i) {
        betas_col <- beta_mat[, i]
        imputed_col <- imputeBetasByGenomicNeighbors(
            betas = betas_col,
            platform = platform,
            BPPARAM = SerialParam(),  # inner function is not parallelized
            max_neighbors = max_neighbors,
            max_dist = max_dist
        )
        imputed_col
    }, BPPARAM = BPPARAM))

    colnames(imputed_mat) <- colnames(beta_mat)
    rownames(imputed_mat) <- rownames(beta_mat)
    return(imputed_mat)
}

# imputation
impute <- function(betas, platform = "EPICv2") {
    cat("Imputation started...\n")
    original_nas <- sum(is.na(betas))
    cat("Initial NAs:", original_nas, "\n")

    # 1. Genomic neighbor imputation
    cat("Step 1: Genomic neighbor imputation...\n")
    betas_imp2 <- imputeBetasByGenomicNeighborsMatrix(betas, platform)
    nas_after_2 <- sum(is.na(betas_imp2))
    cat("NAs after genomic neighbor imputation:", nas_after_2, "\n")

    # 2. KNN imputation for remaining NAs
    if (nas_after_2 > 0) {
        cat("Step 2: KNN imputation...\n")
        if (!is.matrix(betas_imp2)) {
            betas_imp2 <- as.matrix(betas_imp2)
        }
        betas_final <- impute.knn(betas_imp2, k = 10, rng.seed = 123)$data
        nas_final <- sum(is.na(betas_final))
        cat("NAs after KNN imputation:", nas_final, "\n")
    } else {
        betas_final <- betas_imp2
        nas_final <- 0
    }

    cat("Imputation complete.\n")
    return(betas_final)
}

#' Create SummarizedExperiment with Quality Control
#'
#' @param betas Beta value matrix
#' @param sample_sheet Sample metadata
#' @param factors Character vector of factor column names to check
#' @return Filtered SummarizedExperiment object
create_se_with_qc <- function(betas, sample_sheet, factors) {
    cat("Creating SummarizedExperiment...\n")
    cat("Input dimensions - Betas:", dim(betas), "Sample sheet:", dim(sample_sheet), "\n")

    # Create initial SE
    se <- SummarizedExperiment(betas, colData = sample_sheet)

    # Quality control checks
    cat("Performing quality control checks...\n")
    se_ok <- rep(TRUE, nrow(se))

    for (factor_name in factors) {
        if (factor_name %in% colnames(colData(se))) {
            factor_check <- checkLevels(assay(se), colData(se)[[factor_name]])
            se_ok <- se_ok & factor_check
            cat("Factor", factor_name, "- Valid probes:", sum(factor_check), "\n")
        } else {
            warning("Factor '", factor_name, "' not found in sample data")
        }
    }

    cat("Total probes passing QC:", sum(se_ok), "out of", length(se_ok), "\n")

    # Filter to good probes
    se_filtered <- se[se_ok, ]

    cat("Final SE dimensions:", dim(se_filtered), "\n")
    return(se_filtered)
}

#' Run Differential Methylation Analysis
#'
#' @param se SummarizedExperiment object
#' @param formula Model formula as string
#' @param analysis_name Name for progress messages
#' @param n_workers Number of parallel workers
#' @return DML results
run_dml_analysis <- function(se, formula, analysis_name, n_workers = N_WORKERS) {
    cat("=== RUNNING", toupper(analysis_name), "ANALYSIS ===\n")
    cat("Formula:", formula, "\n")
    cat("Samples:", ncol(se), "Probes:", nrow(se), "\n")

    start_time <- Sys.time()

    tryCatch({
        cat("Starting DML analysis...\n")
        smry <- DML(se, as.formula(formula),
                    BPPARAM = BiocParallel::MulticoreParam(workers = n_workers))

        cat("Extracting test results...\n")
        res <- summaryExtractTest(smry)

        end_time <- Sys.time()
        cat("Analysis completed in:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
        cat("Results dimensions:", dim(res), "\n")
        cat("First few results:\n")
        print(head(res, 3))
        cat("\n")

        return(res)

    }, error = function(e) {
        cat("ERROR in", analysis_name, "analysis:", e$message, "\n")
        cat("Trying with fewer workers...\n")

        # Retry with fewer workers
        smry <- DML(se, as.formula(formula),
                    BPPARAM = BiocParallel::MulticoreParam(workers = max(1, n_workers/2)))
        res <- summaryExtractTest(smry)

        end_time <- Sys.time()
        cat("Analysis completed (with reduced workers) in:",
            round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

        return(res)
    })
}


preparePlotFDR <- function(df, n_min, n_max, max_fdr) {
    df <- df[df$nD >0,]
    df$FDR[df$FDR==0] <- .Machine$double.xmin
    df <- df[order(df$FDR),]
    if (sum(df$FDR < max_fdr) < n_min) {
        df1 <- head(df, n=n_min)
    } else {
        df_filtered <- df[df$estimate > 0,] # enrichment only, exclude depletion
        df_filtered_sig <- df_filtered[df_filtered$FDR < max_fdr,]
        if (nrow(df_filtered_sig) < n_min) {
            df1 <- head(df_filtered, n=n_min)  # Take top n_min from filtered set
        } else {
            df1 <- head(df_filtered_sig, n=n_max)
        }
    }
    gp <- vapply(str_split(df1$group, "\\."), function(x) x[3], character(1))
    if ("Target" %in% colnames(df1)) {
        df1$db1 <- sprintf("%s: %s (%s)", gp, df1$Target, df1$dbname)
    } else {
        df1$db1 <- sprintf("%s: %s", gp, df1$dbname)
    }
    df1$db1 <- factor(df1$db1, levels=rev(df1$db1))
    df1$dbname <- factor(df1$dbname, levels=rev(df1$dbname))
    df1
}

plotDotFDR <- function(df, n_min = 10, n_max = 10, max_fdr = 0.05) {
    db1 <- FDR <- overlap <- estimate <- NULL
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))
    df1 <- preparePlotFDR(df, n_min, n_max, max_fdr)
    ggplot(df1) +
        geom_point(aes(-log10(FDR), dbname, size=overlap, color=estimate)) +
        scale_color_gradientn(colors = brewer.pal(9, "Purples")[4:9]) +
        xlab("-log10(FDR)") + ylab("") +
        theme_minimal()
}







# sv_waterfall2 <- function(object, row_id = 1L, max_display = 10L,
#                           order_fun = function(s) order(abs(s), decreasing = TRUE), # Changed to decreasing=TRUE
#                           fill_colors = c("#f7d13d", "#a52c60"),
#                           format_shap = getOption("shapviz.format_shap"),
#                           format_feat = getOption("shapviz.format_feat"),
#                           contrast = TRUE, show_connection = TRUE,
#                           show_annotation = TRUE, annotation_size = 3.2, ...) {
#     stopifnot(
#         "Exactly two fill colors must be passed" = length(fill_colors) == 2L,
#         "format_shap must be a function" = is.function(format_shap),
#         "format_feat must be a function" = is.function(format_feat),
#         "order_fun must be a function" = is.function(order_fun)
#     )
#     object <- object[row_id, ]
#     b <- get_baseline(object)
#     dat <- .make_dat(object, format_feat = format_feat, sep = " = ")
#     if (ncol(object) > max_display) {
#         dat <- .collapse(dat, max_display = max_display)
#     }
#     m <- nrow(dat)
#
#     # Add order dependent columns
#     dat <- dat[order_fun(dat$S), ]
#     dat$i <- seq_len(m)
#     dat$to <- cumsum(dat$S) + b
#     dat$from <- .lag(dat$to, default = b)
#
#     # Make a waterfall plot
#     height <- grid::unit(1 / (1 + 2 * m), "npc")
#
#     p <- ggplot2::ggplot(
#         dat,
#         ggplot2::aes(
#             xmin = from,
#             xmax = to,
#             y = stats::reorder(label, i),
#             fill = factor(to < from, levels = c(FALSE, TRUE))
#         )
#     ) +
#         gggenes::geom_gene_arrow(
#             show.legend = FALSE,
#             arrowhead_width = grid::unit(2, "mm"),
#             arrowhead_height = height,
#             arrow_body_height = height
#         ) +
#         ggfittext::geom_fit_text(
#             ggplot2::aes(label = paste0(ifelse(S > 0, "+", ""), format_shap(S))),
#             show.legend = FALSE,
#             contrast = contrast,
#             ...
#         ) +
#         ggplot2::scale_fill_manual(values = fill_colors, drop = FALSE) +
#         ggplot2::theme_bw() +
#         ggplot2::theme(
#             panel.border = ggplot2::element_blank(),
#             panel.grid.minor = ggplot2::element_blank(),
#             panel.grid.major.x = ggplot2::element_blank(),
#             axis.line.x = ggplot2::element_line(),
#             axis.ticks.y = ggplot2::element_blank()
#         ) +
#         ggplot2::labs(y = ggplot2::element_blank(), x = "Prediction")
#
#     if (show_connection) {
#         p <- p +
#             ggplot2::geom_segment(
#                 ggplot2::aes(x = to, xend = to, y = i, yend = .lag(i, lead = TRUE, default = m)),
#                 linewidth = 0.3,
#                 linetype = 2
#             )
#     }
#
#     if (show_annotation) {
#         full_range <- c(dat[m, "to"], dat[1L, "from"])
#         p <- p +
#             ggplot2::annotate(
#                 "segment",
#                 x = full_range,
#                 xend = full_range,
#                 y = c(m, 1),
#                 yend = c(m, 1) + m * c(0.075, -0.075) + 0.13 * c(1, -1),
#                 linewidth = 0.3,
#                 linetype = 2
#             ) +
#             ggplot2::annotate(
#                 "text",
#                 x = full_range,
#                 y = c(m, 1) + m * c(0.1, -0.1) + 0.15 * c(1, -1),
#                 label = paste0(c("f(x)=", "E[f(x)]="), format_shap(full_range)),
#                 size = annotation_size
#             ) +
#             ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
#             ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.3, mult = 0.2)) +
#             ggplot2::coord_cartesian(clip = "off")
#     }
#     p
# }


# Helper functions for sv_waterfall() and sv_force()
.lag <- function(z, default = NA, lead = FALSE) {
    n <- length(z)
    if (n < 2L) {
        return(rep(default, times = n))
    }
    if (isTRUE(lead)) {
        return(c(z[2L:n], default))
    }
    c(default, z[1L:(n - 1L)])
}

# Turns "shapviz" object into a two-column data.frame
.make_dat <- function(object, format_feat, sep = " = ") {
    X <- get_feature_values(object)
    S <- get_shap_values(object)
    if (nrow(object) == 1L) {
        S <- drop(S)
        label <- paste(colnames(X), format_feat(X), sep = sep)
    } else {
        message("Aggregating SHAP values over ", nrow(object), " observations")
        S <- colMeans(S)
        J <- vapply(X, function(z) length(unique(z)) <= 1L, FUN.VALUE = TRUE)
        label <- colnames(X)
        if (any(J)) {
            label[J] <- paste(label[J], format_feat(X[1L, J]), sep = sep)
        }
    }
    data.frame(S = S, label = label)
}

# Used to combine unimportant rows in dat. dat has two columns: S and label
# Note: rownames(dat) = colnames(object)
.collapse <- function(dat, max_display) {
    m_drop <- nrow(dat) - max_display + 1L
    drop_cols <- rownames(dat)[order(abs(dat$S))[seq_len(m_drop)]]
    keep_cols <- setdiff(rownames(dat), drop_cols)
    rbind(
        dat[keep_cols, ],
        data.frame(
            S = sum(dat[drop_cols, "S"]),
            label = paste(length(drop_cols), "other features"),
            row.names = "other"
        )
    )
}

# Adds non-null titles "nms" to list of ggplots
add_titles <- function(a_list, nms = NULL) {
    if (is.null(nms)) {
        return(a_list)
    }
    mapply(function(p, nm) p + ggplot2::ggtitle(nm), a_list, nms, SIMPLIFY = FALSE)
}


sv_importance.shapviz <- function(object, kind = c("bar", "beeswarm", "both", "no"),
                                  max_display = 15L, fill = "#fca50a", bar_width = 2/3,
                                  bee_width = 0.4, bee_adjust = 0.5,
                                  viridis_args = getOption("shapviz.viridis_args"),
                                  color_bar_title = "Feature value",
                                  show_numbers = FALSE, format_fun = format_max,
                                  number_size = 3.2, sort_features = TRUE, ...) {
    stopifnot("format_fun must be a function" = is.function(format_fun))
    kind <- match.arg(kind)
    imp <- .get_imp(get_shap_values(object), sort_features = sort_features)
    if (kind == "no") {
        return(imp)
    }
    # Deal with too many features
    if (ncol(object) > max_display) {
        imp <- imp[seq_len(max_display)]
    }
    ord <- names(imp)
    object <- object[, ord]  # not required for kind = "bar"
    # ggplot will need to work with data.frame
    imp_df <- data.frame(feature = factor(ord, rev(ord)), value = imp)
    is_bar <- kind == "bar"

    # Create minimal theme
    minimal_theme <- ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(color = "gray90"),
            panel.background = ggplot2::element_blank(),
            plot.background = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_line(color = "black", size = 0.2),
            legend.background = ggplot2::element_blank(),
            legend.box.spacing = grid::unit(0, "pt")
        )

    if (is_bar) {
        p <- ggplot2::ggplot(imp_df, ggplot2::aes(x = value, y = feature)) +
            ggplot2::geom_bar(fill = fill, width = bar_width, stat = "identity", ...) +
            ggplot2::labs(x = "mean(|SHAP value|)", y = ggplot2::element_blank()) +
            minimal_theme
    } else {
        # Prepare data.frame for beeswarm plot
        S <- get_shap_values(object)
        X <- .scale_X(get_feature_values(object))
        df <- transform(
            as.data.frame.table(S, responseName = "value"),
            feature = factor(Var2, levels = rev(ord)),
            color = as.data.frame.table(X)$Freq
        )
        p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = feature))
        if (kind == "both") {
            p <- p +
                ggplot2::geom_bar(
                    data = imp_df, fill = fill, width = bar_width, stat = "identity"
                )
        }
        p <- p +
            ggplot2::geom_vline(xintercept = 0, color = "darkgray") +
            ggplot2::geom_point(
                ggplot2::aes(color = color),
                position = position_bee(width = bee_width, adjust = bee_adjust),
                size = 1,
                ...
            ) +
            .get_color_scale(
                viridis_args = viridis_args,
                bar = !is.null(color_bar_title),
                ncol = length(unique(df$color))   # Special case of constant feature values
            ) +
            ggplot2::labs(
                x = "SHAP value", y = ggplot2::element_blank(), color = color_bar_title
            ) +
            minimal_theme
    }
    if (show_numbers) {
        p <- p +
            ggplot2::geom_text(
                data = imp_df,
                ggplot2::aes(
                    x = if (is_bar) value + max(value) / 60 else
                        min(df$value) - diff(range(df$value)) / 20,
                    label = format_fun(value)
                ),
                hjust = !is_bar,
                size = number_size
            ) +
            ggplot2::scale_x_continuous(
                expand = ggplot2::expansion(mult = 0.05 + c(0.12 *!is_bar, 0.09 * is_bar))
            )
    }
    p
}

# Helper functions
.min_max_scale <- function(z, na.rm = TRUE) {
    r <- range(z, na.rm = na.rm)
    d <- diff(r)
    if (is.na(d) || d == 0) {
        z[!is.na(z)] <- 0.5
        return(z)
    }
    (z - r[1L]) /(r[2L] - r[1L])
}

.get_imp <- function(z, sort_features = TRUE) {
    if (is.matrix(z)) {
        imp <- colMeans(abs(z))
        if (sort_features) {
            imp <- sort(imp, decreasing = TRUE)
        }
        return(imp)
    }
    # list/mshapviz
    imp <- sapply(z, function(x) colMeans(abs(x)))
    if (sort_features) {
        imp <- imp[order(-rowSums(imp)), ]
    }
    return(imp)
}

.scale_X <- function(X) {
    X_scaled <- apply(data.matrix(X), 2L, FUN = .min_max_scale)
    if (nrow(X) == 1L) t(X_scaled) else X_scaled
}

# ncol < 2 treats the special case of constant feature values (e.g., if n = 1)
.get_color_scale <- function(viridis_args, bar = TRUE, ncol = 2L) {
    if (bar) {
        viridis_args_plus <-
            list(
                breaks = if (ncol >= 2L) 0:1 else 0.5,
                labels = if (ncol >= 2L) c("Low", "High") else "Avg",
                guide = ggplot2::guide_colorbar(
                    barwidth = 0.4,
                    barheight = 8,
                    title.theme = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0),
                    title.position = "left"
                )
            )
    } else {
        viridis_args_plus <- list(guide = "none")
    }
    return(do.call(ggplot2::scale_color_viridis_c, c(viridis_args, viridis_args_plus)))
}


#' Number Formatter
#'
#' Formats a numeric vector in a way that its largest absolute value determines
#' the number of digits after the decimal separator. This function is helpful in
#' perfectly aligning numbers on plots. Does not use scientific formatting.
#'
#' @param x A numeric vector to be formatted.
#' @param digits Number of significant digits of the largest absolute value.
#' @param ... Further arguments passed to [format()], e.g., `big.mark = "'"`.
#' @returns A character vector of formatted numbers.
#' @export
#' @examples
#' x <- c(100, 1, 0.1)
#' format_max(x)
#'
#' y <- c(100, 1.01)
#' format_max(y)
#' format_max(y, digits = 5)
format_max <- function(x, digits = 4L, ...) {
    mx <- trunc(log10(max(abs(x), na.rm = TRUE))) + 1L
    x_rounded <- round(x, pmax(0L, digits - mx))
    format(x_rounded, scientific = FALSE, trim = TRUE, ...)
}


# Beeswarm position
position_bee <- function(width = NULL, adjust = NULL) {
    ggplot2::ggproto(NULL, PositionBee, width = width, adjust = adjust)
}

PositionBee <- ggplot2::ggproto(
    "PositionBee",
    ggplot2::Position,
    required_aes = c("x", "y"),

    setup_params = function(self, data) {
        list(
            width = if (!is.null(self$width)) self$width else
                ggplot2::resolution(data$y, zero = FALSE) * 0.4,
            adjust = if (!is.null(self$adjust)) self$adjust else 0.5
        )
    },

    compute_panel = function(self, data, params, scales) {
        data <- ggplot2::flip_data(data, params$flipped_aes)
        y_jit <- ave2(data$x, g = data$y, FUN = shifter, adjust = params$adjust)
        data <- ggplot2::transform_position(
            data, trans_y = function(y) y + y_jit * params$width
        )
        ggplot2::flip_data(data, params$flipped_aes)
    }
)

#===========================================================================
# Helper functions to produce the beeswarm plots
#===========================================================================

# Example
# ggplot(iris, aes(Species, Sepal.Width)) +
#   geom_point(position = position_bee(), aes(color = Species))

# Beeswarm position
position_bee <- function(width = NULL, adjust = NULL) {
    ggplot2::ggproto(NULL, PositionBee, width = width, adjust = adjust)
}

PositionBee <- ggplot2::ggproto(
    "PositionBee",
    ggplot2::Position,
    required_aes = c("x", "y"),

    setup_params = function(self, data) {
        list(
            width = if (!is.null(self$width)) self$width else
                ggplot2::resolution(data$y, zero = FALSE) * 0.4,
            adjust = if (!is.null(self$adjust)) self$adjust else 0.5
        )
    },

    compute_panel = function(self, data, params, scales) {
        data <- ggplot2::flip_data(data, params$flipped_aes)
        y_jit <- ave2(data$x, g = data$y, FUN = shifter, adjust = params$adjust)
        data <- ggplot2::transform_position(
            data, trans_y = function(y) y + y_jit * params$width
        )
        ggplot2::flip_data(data, params$flipped_aes)
    }
)

# Shift values according to their density in the unit interval by quasi-random numbers
shifter <- function(y, ...) {
    if (length(y) == 1L) {
        return(0)
    }
    dens <- stats::density(y, ...)
    dens_y <- dens[["y"]] / max(dens[["y"]])
    shift <- halton_sequence(length(y))[rank(y, ties.method = "first")] - 0.5
    2 * shift * stats::approx(dens[["x"]], dens_y, xout = y)[["y"]]
}

# "stats::ave" for grouping variable "g" and additional arguments ...
ave2 <- function(x, g = NULL, FUN = mean, ...) {
    if (is.null(g)) {
        x[] <- FUN(x, ...)
    } else {
        split(x, g) <- lapply(split(x, g), FUN, ...)
    }
    x
}

# First n values of the 1-dimensional Halton sequence (= van der Corput sequence)
# https://en.wikipedia.org/wiki/Halton_sequence
halton_sequence <- function(n, b = 2) {
    vapply(seq_len(n), halton, FUN.VALUE = 0.0)
}

# i-th element of above sequence
halton <- function(i, b = 2) {
    f <- 1
    r <- 0
    while (i > 0) {
        f <- f / b
        r <- r + f * (i %% b)
        i <- trunc(i / b)
    }
    r
}




