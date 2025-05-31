# ============================================================
#   FIG 4.
#   This script contains the analysis for t-SNE
#   analysis of classifier results, including:
#   4D Invasiveness results evaluation
#   4E Driver results evaluation
#   4F LN + primary t-SNE evaluation
# ============================================================
# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- "/Users/jennyzli/Documents/HPC_share/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
FIG_DIR <- file.path(BASE_DIR, "figures")
DATA_DIR <- file.path(BASE_DIR, "data")

CL_DIR <- file.path(DATA_DIR, "loocv_driver_model")
SUM_DIR <- file.path(CL_DIR, "summary")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# LOADING DATA
# ========================
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_ped <- ss %>% filter(Lymph_Node == "F")
betas <- readRDS(file.path(DATA_DIR, "ped_betas_imputed.rds"))

# ========================
# LOADING + MERGING
# ========================
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx")) %>% filter(Lymph_Node == "F")
tsne <- readRDS(file.path(DATA_DIR, "ped88_tsne_coords.rds"))
inv_res <- read.csv(file.path(DATA_DIR, "loocv_invasiveness_model/summary/loocv_all_predictions.csv"))
driver_res <- read.csv(file.path(SUM_DIR, "loocv_all_predictions.csv"))

## UPDATE COLUMN NAMES
inv_conf <- inv_res %>%
    transmute(Sample_ID, Predicted_Invasiveness = Predicted,
              Predicted_Invasiveness_Confidence = pmax(Prob_High, Prob_Low),
              Invasiveness_Accuracy = Correct)

driver_conf <- driver_res %>%
    transmute(Sample_ID, Driver_Group = Predicted,
              LOOCV_Driver_Confidence = pmax(Prob_BRAF.V600E, Prob_DICER1, Prob_Kinase.Fusion, Prob_Ras.like),
              Driver_Accuracy = Correct)

merged <- ss %>%
    left_join(tsne, by = "IDAT") %>%
    left_join(inv_conf, by = "Sample_ID") %>%
    left_join(driver_conf, by = "Sample_ID")

# ========================
# PLOTTING FUNCTION
# ========================
create_tsne_plot <- function(data, color_var, color_values, shape_var, shape_values, alpha_var) {
    ggplot(data, aes(x = tSNE1, y = tSNE2)) +
        geom_point(aes_string(color = color_var, shape = shape_var, alpha = alpha_var), size = 2, na.rm = TRUE) +
        scale_color_manual(values = color_values) +
        scale_shape_manual(values = shape_values) +
        scale_alpha_continuous(range = c(0.4, 1)) +
        coord_fixed() + common_theme + theme(legend.position = "none")
}

# ========================
# PLOT CONFIGS
# ========================
plot_configs <- list(
    list(name = "invasiveness_confidence", color_var = "Invasiveness", color_values = invasiveness_colors,
         shape_var = "Invasiveness_Accuracy", alpha_var = "LOOCV_Invasiveness_Confidence"),
    list(name = "drivergroup_confidence", color_var = "Driver_Group", color_values = driver_colors,
         shape_var = "Driver_Accuracy", alpha_var = "LOOCV_Driver_Confidence")
)

for (cfg in plot_configs) {
    p <- create_tsne_plot(merged, cfg$color_var, cfg$color_values, cfg$shape_var, custom_shapes, cfg$alpha_var)
    ggsave(file.path(PLOT_PATH, sprintf("%s_tsne88_%s.pdf", PLOT_DATE, cfg$name)),
           plot = p, width = PLOT_SIZE$width, height = PLOT_SIZE$height)
}

# # ========================
# # DIFF METH - INVASIVENESS
# # ========================
# ss <- read_excel(file.path(SS_DIR, "adult_master.xlsx")) %>%
#     filter(Clinical_Invasiveness != "NA")
# betas_adult <- readRDS(file.path(DATA_DIR, "adult_betas_imputed.rds"))
# betas <- readRDS(file.path(DATA_DIR, "adult_betas_imputed.rds"))
#
# common <- intersect(ss$Sample_ID, colnames(betas))
# betas <- betas[, common, drop = FALSE]
# cat("Beta dimensions:", dim(betas), "\n")
#
# ss$Age <- as.numeric(ss$Age)
# ss$Sex <- as.factor(ss$Sex)
# ss$Invasiveness <- as.factor(ss$Invasiveness)

# ========================
# 4D. LOOCV INVASIVENESS RESULTS
# ========================
# load t-SNE coordinates
tsne <- readRDS(file.path(DATA_DIR, "ped88_tsne_coords.rds"))
inv_res <- read.csv(file.path(DATA_DIR, "loocv_invasiveness_model", "summary", "loocv_all_predictions.csv"))
driver_res <- read.csv(file.path(DATA_DIR, "loocv_driver_model", "summary", "loocv_all_predictions.csv"))

# ========================
# INVASIVENESS CONFIDENCE
# ========================
inv_probs <- inv_res[, c("Prob_High", "Prob_Low")]
conf_scores <- apply(inv_probs, 1, max)

inv_conf <- data.frame(
    Sample_ID = inv_res$Sample_ID,
    Predicted_Invasiveness = inv_res$Predicted,
    Predicted_Invasiveness_Confidence = conf_scores,
    Predicted_Invasiveness_Accuracy = inv_res$Correct
)

# reorder inv_conf based on sample ID in ss, then merge it

# ========================
# DRIVER CONFIDENCE
# ========================
driver_probs <- driver_res[, c("Prob_BRAF.V600E", "Prob_DICER1", "Prob_Kinase.Fusion", "Prob_Ras.like")]
conf_scores <- apply(driver_probs, 1, max)

driver_conf <- data.frame(
    Sample_ID = driver_res$Sample_ID,
    Predicted_Driver = driver_res$Predicted,
    Predicted_Driver_Confidence = conf_scores,
    Predicted_Driver_Accuracy = driver_res$Correct
)
# reorder driver_conf based on sample ID in ss, then merge it

# Create tSNE plotting function with NA handling
create_tsne_plot <- function(data, color_var, color_values = NULL, shape_var = NULL,
                             shape_values = NULL, continuous = FALSE, alpha_var = NULL) {
    # Basic plot setup
    p <- ggplot(data = data, aes(x = tSNE1, y = tSNE2)) +
        coord_fixed() +
        common_theme

    # Create a complete data frame excluding NA values in the alpha variable if needed
    if (!is.null(alpha_var)) {
        # For alpha points (with confidence values)
        plot_data <- data[!is.na(data[[alpha_var]]), ]

        # For points without confidence (NAs)
        na_data <- data[is.na(data[[alpha_var]]), ]

        # Add points with confidence (using alpha)
        if (!is.null(shape_var)) {
            p <- p + geom_point(data = plot_data,
                                aes(color = !!sym(color_var),
                                    shape = !!sym(shape_var),
                                    alpha = !!sym(alpha_var)),
                                size = 2) +
                scale_shape_manual(values = shape_values) +
                scale_alpha_continuous(range = c(0.4, 1))
        } else {
            p <- p + geom_point(data = plot_data,
                                aes(color = !!sym(color_var),
                                    alpha = !!sym(alpha_var)),
                                size = 2, shape = 16) +
                scale_alpha_continuous(range = c(0.4, 1))
        }

        # Add NA points separately with fixed alpha
        if (nrow(na_data) > 0) {
            if (!is.null(shape_var)) {
                p <- p + geom_point(data = na_data,
                                    aes(color = !!sym(color_var),
                                        shape = !!sym(shape_var)),
                                    size = 2, alpha = 0.2) # Low alpha for NA points
            } else {
                p <- p + geom_point(data = na_data,
                                    aes(color = !!sym(color_var)),
                                    size = 2, shape = 16, alpha = 0.2) # Low alpha for NA points
            }
        }

    } else {
        # No alpha variable, regular plot
        if (!is.null(shape_var)) {
            p <- p + geom_point(aes(color = !!sym(color_var),
                                    shape = !!sym(shape_var)),
                                size = 2) +
                scale_shape_manual(values = shape_values)
        } else {
            p <- p + geom_point(aes(color = !!sym(color_var)),
                                size = 2, shape = 16)
        }
    }

    # Apply color scale
    if (continuous) {
        p <- p + scale_color_gradientn(colors = parula(20))
    } else {
        p <- p + scale_color_manual(values = color_values)
    }

    return(p)
}

# Configure plots
plot_configs <- list(
    list(
        name = "invasiveness_confidence",
        color_var = "Invasiveness",
        color_values = invasiveness_colors,
        shape_var = "Invasiveness_Accuracy",
        shape_values = custom_shapes,
        continuous = FALSE,
        alpha_var = "LOOCV_Invasiveness_Confidence"
    ),
    list(
        name = "drivergroup_confidence",
        color_var = "Driver_Group",
        color_values = driver_colors,
        shape_var = "Driver_Accuracy",
        shape_values = custom_shapes,
        continuous = FALSE,
        alpha_var = "LOOCV_Driver_Confidence"
    )
)

for (config in plot_configs) {
    tryCatch({
        plot <- create_tsne_plot(
            data = ss,
            color_var = config$color_var,
            color_values = config$color_values,
            shape_var = config$shape_var,
            shape_values = config$shape_values,
            continuous = config$continuous,
            alpha_var = config$alpha_var
        )

        # Remove all legends
        plot <- plot + theme(legend.position = "none")

        # Create directory if it doesn't exist
        if (!dir.exists(PLOT_PATH)) {
            dir.create(PLOT_PATH, recursive = TRUE)
        }

        filename <- file.path(PLOT_PATH, sprintf("%s_tsne88_%s_no_legend.pdf", PLOT_DATE, config$name))

        # Save the plot
        ggsave(
            filename = filename,
            plot = plot,
            width = PLOT_SIZE$width,
            height = PLOT_SIZE$height,
            units = "in",
            device = "pdf"
        )

        cat("Successfully created plot:", filename, "\n")
    }, error = function(e) {
        cat("Error creating plot for", config$name, ":", conditionMessage(e), "\n")
    })
}

create_confidence_legend <- function(confidence_var_name) {
    # Create a dummy plot with only the confidence scale
    dummy_data <- data.frame(
        x = 1:100,
        y = 1:100,
        confidence = seq(min(ss[[confidence_var_name]], na.rm = TRUE),
                         max(ss[[confidence_var_name]], na.rm = TRUE),
                         length.out = 100)
    )

    # Create the dummy plot with only confidence represented
    p <- ggplot(dummy_data, aes(x = x, y = y, size = confidence)) +
        geom_point(alpha = 0) +  # Invisible points
        scale_size_continuous(range = c(1, 3), name = confidence_var_name) +
        theme_void() +  # Remove all plot elements
        theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.width = unit(2, "cm"),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8)
        )

    # Extract just the legend
    legend <- cowplot::get_legend(p)

    # Create a blank plot with just the legend
    legend_plot <- cowplot::ggdraw() +
        cowplot::draw_grob(legend)

    return(legend_plot)
}

confidence_legend <- create_confidence_legend(config$alpha_var)

legend_filename <- file.path(PLOT_PATH,
                             sprintf("%s_confidence_legend_%s.pdf",
                                     PLOT_DATE, config$name))

ggsave(
    filename = legend_filename,
    plot = confidence_legend,
    width = 4,
    height = 1,
    units = "in",
    device = "pdf"
)

# ========================
# 4D. LOOCV DRIVER RESULTS
# ========================



# ========================
# 4F. TSNE COORDS 136
# ========================
# subset most variable probes
mtx = t(bSubMostVariable(betas, 30000))

pca_result <- prcomp(mtx, center = TRUE, scale. = TRUE)

# variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
cumvar_explained <- cumsum(var_explained)

# scree plot - how many PCs to use
scree_data <- data.frame(
    PC = 1:length(var_explained),
    Variance = var_explained,
    Cumulative = cumvar_explained
)

p1 <- ggplot(scree_data[1:30, ], aes(x = PC, y = Variance)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 2) +
    labs(title = "Scree Plot",
         x = "Principal Component",
         y = "Variance Explained (%)") +
    theme_minimal()
plot(p1)

# 20 as cutoff from visually examining scree plot
pca_data <- pca_result$x[, 1:20]

set.seed(1234)
tsne <- Rtsne(
    pca_data,
    dims = 2,
    perplexity = 9,
    max_iter = 2000,
    pca = FALSE
)

df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
ss$tSNE1 <- as.numeric(df$tSNE1)
ss$tSNE2 <- as.numeric(df$tSNE2)

rownames(df) <- ss$IDAT
saveRDS(df, file.path(DATA_DIR, "ped98_tsne_coords.rds"))

ggplot(df, aes(x = tSNE1, y = tSNE2)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_minimal() +
    labs(title = "t-SNE plot of samples",
         x = "t-SNE 1",
         y = "t-SNE 2")

# ========================
# 4F. PAIRED t-SNE PLOT
# ========================
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
df <- readRDS(file.path(DATA_DIR, "ped98_tsne_coords.rds"))

ss$tSNE1 <- df$tSNE1
ss$tSNE2 <- df$tSNE2

lymph_nodes <- ss %>%
    filter(Lymph_Node == "T")

pairs_data <- data.frame()
n_pairs <- 10

# Find paired samples and assign a pair ID to each
pair_id <- 1
for (i in 1:nrow(lymph_nodes)) {
    ln <- lymph_nodes[i, ]
    if (is.null(ln$Paired_Primary) || is.na(ln$Paired_Primary) || ln$Paired_Primary == "") {
        next
    }

    primary <- ss %>%
        filter(Sample_ID == ln$Paired_Primary)

    # If found, add both to pairs_data with the same pair_id
    if (nrow(primary) > 0) {
        ln_row <- ln %>%
            mutate(pair_id = paste0("Pair ", pair_id),
                   sample_type = "Lymph Node")

        primary_row <- primary %>%
            mutate(pair_id = paste0("Pair ", pair_id),
                   sample_type = "Primary Tumor")

        pairs_data <- rbind(pairs_data, ln_row, primary_row)
        pair_id <- pair_id + 1
    }
}

pairs_data$Pair_ID <- substr(pairs_data$Sample_ID, 1, 7)

connections <- data.frame()
for (current_pair_id in unique(pairs_data$pair_id)) {
    pair_samples <- pairs_data %>% filter(pair_id == current_pair_id)
    if (nrow(pair_samples) == 2) {
        ln <- pair_samples %>% filter(sample_type == "Lymph Node")
        primary <- pair_samples %>% filter(sample_type == "Primary Tumor")

                connections <- rbind(connections, data.frame(
            pair_id = current_pair_id,
            Pair_ID = ln$Pair_ID,
            x1 = ln$tSNE1,
            y1 = ln$tSNE2,
            x2 = primary$tSNE1,
            y2 = primary$tSNE2
        ))
    }
}

n_colors_needed <- n_pairs

pair_color_mapping <- setNames(
    pair_colors[1:n_pairs],
    unique(pairs_data$Pair_ID)
)

ss <- ss %>%
    mutate(
        is_paired = Sample_ID %in% pairs_data$Sample_ID
    )

p <- ggplot() +
    # background points for non paired samples
    geom_point(data = ss %>% filter(!is_paired),
               aes(x = tSNE1, y = tSNE2, text = paste("Sample ID:", Sample_ID)),
               color = "grey80", alpha = 0.4, size = 1.5) +
    # connection lines
    geom_segment(data = connections,
                 aes(x = x1, y = y1, xend = x2, yend = y2, color = Pair_ID,  # Use Pair_ID here
                     text = paste("Pair:", Pair_ID)),
                 linetype = "dashed", size = 0.7) +
    # points for paired samples
    geom_point(data = pairs_data,
               aes(x = tSNE1, y = tSNE2, color = Pair_ID,
                   shape = sample_type, size = sample_type,
                   text = paste("Sample ID:", Sample_ID,
                                "\nPair:", Pair_ID,
                                "\nType:", sample_type)),
               alpha = 0.8) +
    scale_color_manual(values = pair_color_mapping, name = "Sample ID") +
    scale_shape_manual(values = c("Lymph Node" = 17, "Primary Tumor" = 1),
                       name = "Sample Type") +
    scale_size_manual(values = c("Lymph Node" = 3, "Primary Tumor" = 3.5),
                      name = "Sample Type") +
    guides(
        color = guide_legend(override.aes = list(size = 3)),
        shape = guide_legend(override.aes = list(size = 3))
    ) +
    common_theme +
    theme(
        legend.position = "right",
        legend.box = "vertical"
    )

output_file <- file.path(FIG_DIR, "ped_tsne_paired_ln.pdf")
ggsave(output_file, p, width = 7, height = 5)











# # ========================
# # PLOTTING
# # ========================
# create_tsne_plot_with_labels <- function(data, color_var, color_values = NULL, shape_var = NULL,
#                                          shape_values = NULL, continuous = FALSE, title, legend_title = NULL) {
#     p <- ggplot(data = data, aes(x = tSNE1, y = tSNE2)) +
#         coord_fixed() +
#         common_theme_with_legend
#     if (!is.null(shape_var)) {
#         p <- p + geom_point(aes(color = !!sym(color_var),
#                                 shape = !!sym(shape_var)),
#                             size = 2) +
#             scale_shape_manual(values = shape_values)
#     } else {
#         p <- p + geom_point(aes(color = !!sym(color_var)),
#                             size = 2, shape = 16)
#     }
#     if (continuous) {
#         p <- p + scale_color_gradientn(colors = parula(20))
#     } else {
#         p <- p + scale_color_manual(values = color_values)
#     }
#     if (is.null(legend_title)) {
#         legend_title <- prettify_name(color_var)
#     }
#     p <- p + labs(title = title, color = legend_title)
#     return(p)
# }
#
# prettify_name <- function(name) {
#     name <- gsub("_", " ", name)
#     name <- tools::toTitleCase(name)
#     return(name)
# }
#
# for (config in plot_list) {
#     plot_title <- paste("t-SNE plot colored by", prettify_name(config$name))
#     legend_title <- prettify_name(config$color_var)
#     plot <- create_tsne_plot_with_labels(
#         data = config$data,
#         color_var = config$color_var,
#         color_values = config$color_values,
#         shape_var = config$shape_var,
#         shape_values = config$shape_values,
#         continuous = if (!is.null(config$continuous)) config$continuous else FALSE,
#         title = plot_title,
#         legend_title = legend_title
#     )
#     filename <- file.path(FIG_DIR,
#                           sprintf("tsne88_%s.pdf", config$name))
#     ggsave(filename = filename,
#            plot = plot,
#            width = PLOT_SIZE$width,
#            height = PLOT_SIZE$height,
#            units = "in",
#            device = "pdf")
# }
#
# # ==========================
# # TSNE WITHOUT LEGEND (figure creation)
# # ==========================
# create_tsne_plot <- function(data, color_var, color_values = NULL, shape_var = NULL,
#                              shape_values = NULL, continuous = FALSE) {
#     p <- ggplot(data = data, aes(x = tSNE1, y = tSNE2)) +
#         coord_fixed() +
#         common_theme
#
#     if (!is.null(shape_var)) {
#         p <- p + geom_point(aes(color = !!sym(color_var),
#                                 shape = !!sym(shape_var)),
#                             size = 2) +
#             scale_shape_manual(values = shape_values)
#     } else {
#         p <- p + geom_point(aes(color = !!sym(color_var)),
#                             size = 2, shape = 16)
#     }
#
#     if (continuous) {
#         p <- p + scale_color_gradientn(colors = parula(20))
#     } else {
#         p <- p + scale_color_manual(values = color_values)
#     }
#
#     return(p)
# }
#
# for (config in plot_list) {
#     plot <- create_tsne_plot(
#         data = config$data,
#         color_var = config$color_var,
#         color_values = config$color_values,
#         shape_var = config$shape_var,
#         shape_values = config$shape_values,
#         continuous = if (!is.null(config$continuous)) config$continuous else FALSE
#     )
#
#     filename <- file.path(FIG_DIR,
#                           sprintf("tsne88_%s.pdf",
#                                   config$name))
#
#     ggsave(filename = filename,
#            plot = plot,
#            width = PLOT_SIZE$width,
#            height = PLOT_SIZE$height,
#            units = "in",
#            device = "pdf")
# }
#
#
#
#
#
#
#
# ss$
#     # configure for loop
#     plot_list <- list(
#         list(name = "invasiveness", data = ss_ped, color_var = "Clinical_Invasiveness", shape = accuracy_shapes, color_values = invasiveness_colors),
#         list(name = "driver_group", data = ss_ped, color_var = "Driver_Group", shape = accuracy_shapes, color_values = driver_colors)
#     )
#
# PLOT_SIZE <- list(width = 5.5, height = 4.5)
#
# common_theme <- theme(
#     legend.title = element_blank(),
#     panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
#     panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray"),
#     aspect.ratio = 1,
#     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
#     legend.position = "none"
# )
#
# common_theme_with_legend <- theme(
#     legend.title = element_blank(),
#     panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
#     panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray"),
#     aspect.ratio = 1,
#     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
#     legend.position = "right",   # or "bottom" if you prefer
#     legend.text = element_text(size = 10),
#     legend.key.size = unit(0.5, "cm")
# )
#
# create_tsne_plot_with_labels <- function(data, color_var, color_values = NULL, shape_var = NULL,
#                                          shape_values = NULL, continuous = FALSE, title, legend_title = NULL) {
#     p <- ggplot(data = data, aes(x = tSNE1, y = tSNE2, shape = shape_var)) +
#         coord_fixed() +
#         common_theme_with_legend
#     if (!is.null(shape_var)) {
#         p <- p + geom_point(aes(color = !!sym(color_var),
#                                 shape = !!sym(shape_var)),
#                             size = 2) +
#             scale_shape_manual(values = shape_values)
#     } else {
#         p <- p + geom_point(aes(color = !!sym(color_var)),
#                             size = 2, shape = 16)
#     }
#     # If legend_title is not given, use color_var prettified as default
#     if (is.null(legend_title)) {
#         legend_title <- prettify_name(color_var)
#     }
#     p <- p + labs(title = title, color = legend_title)
#     return(p)
# }
#
# prettify_name <- function(name) {
#     name <- gsub("_", " ", name)
#     name <- tools::toTitleCase(name)
#     return(name)
# }
#
# for (config in plot_list) {
#     plot_title <- paste("t-SNE plot colored by", prettify_name(config$name))
#     legend_title <- prettify_name(config$color_var)
#     plot <- create_tsne_plot_with_labels(
#         data = config$data,
#         color_var = config$color_var,
#         color_values = config$color_values,
#         shape_var = config$shape_var,
#         shape_values = config$shape_values,
#         continuous = if (!is.null(config$continuous)) config$continuous else FALSE,
#         title = plot_title,
#         legend_title = legend_title
#     )
#     filename <- file.path(FIG_DIR,
#                           sprintf("tsne88_%s.pdf", config$name))
#     ggsave(filename = filename,
#            plot = plot,
#            width = PLOT_SIZE$width,
#            height = PLOT_SIZE$height,
#            units = "in",
#            device = "pdf")
# }
#
#
