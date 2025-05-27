# ============================================================
#   This script creates the t-SNEs for Fig. 1A-K and
#   SFig. 1A-C.
# ============================================================
# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- here()
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")
IDAT_DIR <- file.path(BASE_DIR, "IDAT")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# LOADING DATA
# ========================
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss <- ss %>% filter(Lymph_Node == "F")
betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))

# ========================
# t-SNE COORDINATES
# ========================
# subset to 3000 most variable probes
mtx = t(bSubMostVariable(betas, 3000))

set.seed(12345678)
pca_result <- prcomp(mtx, center = TRUE, scale. = TRUE)
var_explained <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
n_components <- which(var_explained >= 0.8)[1]
pca_data <- pca_result$x[, 1:20]

# Compare t-SNE results with 20, 30, 42 PCs
for (n in c(20, 30, 42)) {
    tsne <- Rtsne(pca_result$x[, 1:n], dims = 2, perplexity = 9.5, max_iter = 2000, pca = FALSE)
    plot(tsne$Y, main = paste("t-SNE with", n, "PCs"))
}

set.seed(123467)
tsne <- Rtsne(
    mtx,
    dims = 2,
    perplexity = 14,
    max_iter = 2000,
    pca = FALSE
)

df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
ss$tSNE1 <- as.numeric(df$tSNE1)
ss$tSNE2 <- as.numeric(df$tSNE2)

# saveRDS(df, file.path(DATA_DIR, "ped88_tsne_coords.rds"))
ggplot(ss, aes(x = tSNE1, y = tSNE2, color = Methylation_Clusters)) +
    geom_point(alpha = 0.8) +
    theme_minimal() +
    labs(title = "t-SNE with 42 PCs", x = "t-SNE 1", y = "t-SNE 2")

ggplot(ss, aes(x = tSNE1, y = tSNE2, color = Driver_Group)) +
    geom_point(alpha = 0.8) +
    theme_minimal() +
    labs(title = "t-SNE with 42 PCs", x = "t-SNE 1", y = "t-SNE 2")


# ========================
# t-SNE EXPERIMENTATION
# ========================
# subset to 3000 most variable probes
mtx = t(bSubMostVariable(betas, 3000))

betas_centered <- scale(t(mtx), center = TRUE, scale = FALSE)
cov_matrix <- cov(t(betas_centered))  # 88 × 88
eig_result <- eigs(cov_matrix, k = 100)
betas_random <- apply(betas_subset, 1, sample)
eig_random <- eigs(cov_random, k = 100)
threshold <- max(eig_random$values)
nontrivial_indices <- which(eigenvalues > threshold)
pc_scores <- eig_result$vectors[, nontrivial_indices]

tsne_result <- Rtsne(pc_scores, pca = FALSE, theta = 0, max_iter = 2500)


set.seed(12345678)
pca_result <- prcomp(mtx, center = TRUE, scale. = TRUE)
var_explained <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
n_components <- which(var_explained >= 0.8)[1]
pca_data <- pca_result$x[, 1:20]

# Compare t-SNE results with 20, 30, 42 PCs
for (n in c(20, 30, 42)) {
    tsne <- Rtsne(pca_result$x[, 1:n], dims = 2, perplexity = 9.5, max_iter = 2000, pca = FALSE)
    plot(tsne$Y, main = paste("t-SNE with", n, "PCs"))
}

set.seed(123467)
tsne <- Rtsne(
    mtx,
    dims = 2,
    perplexity = 14,
    max_iter = 2000,
    pca = FALSE
)

df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
ss$tSNE1 <- as.numeric(df$tSNE1)
ss$tSNE2 <- as.numeric(df$tSNE2)

# saveRDS(df, file.path(DATA_DIR, "ped88_tsne_coords.rds"))
ggplot(ss, aes(x = tSNE1, y = tSNE2, color = Methylation_Clusters)) +
    geom_point(alpha = 0.8) +
    theme_minimal() +
    labs(title = "t-SNE with 42 PCs", x = "t-SNE 1", y = "t-SNE 2")

ggplot(ss, aes(x = tSNE1, y = tSNE2, color = Driver_Group)) +
    geom_point(alpha = 0.8) +
    theme_minimal() +
    labs(title = "t-SNE with 42 PCs", x = "t-SNE 1", y = "t-SNE 2")



# ==========================
# 1A-H, S1A-C, 4D TSNE PLOTS
# ==========================
# Load t-SNE coordinates and MIR200 predictions
tsne <- readRDS(file.path(DATA_DIR, "ped88_tsne_coords.rds"))
mir <- read.csv(file.path(DATA_DIR, "ped98_mir200.csv"))
mir88 <- mir[match(ss$IDAT, rownames(mir)), ]

ss$tSNE1 <- as.numeric(tsne$tSNE1)
ss$tSNE2 <- as.numeric(tsne$tSNE2)
ss$Age <- as.numeric(ss$Age)
ss$Confidence <- as.factor(ss$Confidence)
ss$MIR200C <- as.numeric(mir$score)

# Load RNA expression data for ZEB1/2
mtx <- as.data.frame(read.csv(file.path(DATA_DIR, "20250113_rna_log2cpmfiltered.csv"), row.names = 1) %>%
                         .[, colnames(.) %in% ss$Sample_ID])

# Filter samples to match RNA data availability
ss_filtered <- ss %>%
    filter(Sample_ID %in% colnames(mtx))

# Ensure matrix columns match sample order
mtx <- mtx[, ss_filtered$Sample_ID]

# Add ZEB1/2 expression data
ss_filtered$ZEB1 <- as.numeric(mtx["ZEB1", ])
ss_filtered$ZEB2 <- as.numeric(mtx["ZEB2", ])

# Update t-SNE coordinates for filtered dataset
df_filtered <- readRDS(file.path(DATA_DIR, "ped88_tsne_coords.rds")) %>%
    .[rownames(.) %in% ss_filtered$IDAT, ]
ss_filtered$tSNE1 <- df_filtered$tSNE1
ss_filtered$tSNE2 <- df_filtered$tSNE2

# Plot configuration
PLOT_PATH <- here("figures")
PLOT_SIZE <- list(width = 5.5, height = 4.5)

common_theme <- theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray"),
    aspect.ratio = 1,
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "none"
)

create_tsne_plot <- function(data, color_var, color_values = NULL, shape_var = NULL,
                             shape_values = NULL, continuous = FALSE) {
    p <- ggplot(data = data, aes(x = tSNE1, y = tSNE2)) +
        coord_fixed() +
        common_theme

    if (!is.null(shape_var)) {
        p <- p + geom_point(aes(color = !!sym(color_var),
                                shape = !!sym(shape_var)),
                            size = 2) +
            scale_shape_manual(values = shape_values)
    } else {
        p <- p + geom_point(aes(color = !!sym(color_var)),
                            size = 2, shape = 16)
    }

    if (continuous) {
        p <- p + scale_color_gradientn(colors = parula(20))
    } else {
        p <- p + scale_color_manual(values = color_values)
    }

    return(p)
}

plot_configs <- list(
    # Original plots (using full dataset)
    list(name = "invasiveness", data = ss, color_var = "Invasiveness", color_values = invasiveness_colors),
    list(name = "drivergroup", data = ss, color_var = "Driver_Group", color_values = driver_colors),
    list(name = "sex", data = ss, color_var = "Sex", color_values = sex_colors),
    list(name = "leukocyte_fraction", data = ss, color_var = "Leukocyte_Fraction", continuous = TRUE),
    list(name = "T", data = ss, color_var = "T", color_values = t_colors),
    list(name = "N", data = ss, color_var = "N", color_values = n_colors),
    list(name = "M", data = ss, color_var = "M", color_values = m_colors),
    list(name = "age", data = ss, color_var = "Actual_Age", continuous = TRUE),
    list(name = "mir200", data = ss, color_var = "MIR200C", continuous = TRUE),

    # ZEB plots (using filtered dataset with RNA data)
    list(name = "zeb1", data = ss_filtered, color_var = "ZEB1", continuous = TRUE),
    list(name = "zeb2", data = ss_filtered, color_var = "ZEB2", continuous = TRUE)
)

for (config in plot_configs) {
    plot <- create_tsne_plot(
        data = config$data,
        color_var = config$color_var,
        color_values = config$color_values,
        shape_var = config$shape_var,
        shape_values = config$shape_values,
        continuous = if (!is.null(config$continuous)) config$continuous else FALSE
    )

    filename <- file.path(PLOT_PATH,
                          sprintf("%s_tsne88_%s.pdf",
                                  DATE,
                                  config$name))

    ggsave(filename = filename,
           plot = plot,
           width = PLOT_SIZE$width,
           height = PLOT_SIZE$height,
           units = "in",
           device = "pdf")
}

