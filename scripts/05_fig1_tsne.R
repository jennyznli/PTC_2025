# ============================================================
#   FIG. 1
#   This script creates the t-SNEs plots for
#   1A-K: Invasiveness, Driver Groups, Leukocyte Fraction, Age,
#         T, N, M, Sex.
#   S1A-C: MIR200C, ZEB1, ZEB2
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
# AGE
# ========================
# p <- ggplot(ss, aes(x = as.numeric(Age))) +
#     geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
#     labs(title = "Histogram of Age", x = "Age", y = "Count") +
#     theme_minimal()
# pdf(file.path(FIG_DIR, "age_histogram.pdf"), width = 6, height = 4)
# print(p)
# dev.off()

mean_age <- mean(as.numeric(ss$Age), na.rm = TRUE)
# 15.27819
sd_age <- sd(as.numeric(ss$Age), na.rm = TRUE)
# 2.783226
range_age <- range(as.numeric(ss$Age), na.rm = TRUE)
#  6.240    21.425

# ========================
# t-SNE COORDINATES
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

# 20 as cutoff
pca_data <- pca_result$x[, 1:20]

set.seed(123)
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
saveRDS(df, file.path(DATA_DIR, "ped88_tsne_coords.rds"))

ggplot(df, aes(x = tSNE1, y = tSNE2)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_minimal() +
    labs(title = "t-SNE plot of samples",
         x = "t-SNE 1",
         y = "t-SNE 2")

# ==========================
# 1A-H, S1A-C TSNE PLOTS
# ==========================
# load t-SNE coordinates, MIR200, leukocyte fractions
tsne <- readRDS(file.path(DATA_DIR, "ped88_tsne_coords.rds"))

mir <- read.csv(file.path(DATA_DIR, "ped_mir200.csv"))
indices <- match(ss$IDAT, mir[,1])
mir88 <- mir[indices, , drop = FALSE]

leuko <- read.csv(file.path(DATA_DIR, "ped_leuko.csv"))
indices <- match(ss$IDAT, leuko[,1])
leuko88 <- leuko[indices, , drop = FALSE]

# extract T, N, M
ss$'T' = str_extract(ss$TNM, "T[^NMR]+")
ss$'N' = str_extract(ss$TNM, "N[^TMR]+")
ss$'M' = str_extract(ss$TNM, "M[^TNR]+")

ss$tSNE1 <- tsne$tSNE1
ss$tSNE2 <- tsne$tSNE2

ss$Age <- as.numeric(ss$Age)
ss$Confidence <- as.factor(ss$Confidence)
ss$MIR200C <- as.numeric(mir88[,2])
ss$Leukocyte_Fraction <- as.numeric(leuko88[,2])

# RNA expression for ZEB1/2
rna <- as.data.frame(read.csv(file.path(DATA_DIR, "rna_log2cpmfiltered.csv"), row.names = 1) %>%
                         .[, colnames(.) %in% ss$Sample_ID])

ss_filtered <- ss %>%
    filter(Sample_ID %in% colnames(rna))

rna <- rna[, ss_filtered$Sample_ID]

ss_filtered$ZEB1 <- as.numeric(rna["ZEB1", ])
ss_filtered$ZEB2 <- as.numeric(rna["ZEB2", ])

# t-SNE coordinates for filtered ss
df_filtered <- readRDS(file.path(DATA_DIR, "ped88_tsne_coords.rds")) %>%
    .[rownames(.) %in% ss_filtered$IDAT, ]
ss_filtered$tSNE1 <- df_filtered$tSNE1
ss_filtered$tSNE2 <- df_filtered$tSNE2

# configure loop
plot_list <- list(
    # Original plots (using full dataset)
    list(name = "invasiveness", data = ss, color_var = "Clinical_Invasiveness", color_values = invasiveness_colors),
    list(name = "driver_group", data = ss, color_var = "Driver_Group", color_values = driver_colors),
    list(name = "sex", data = ss, color_var = "Sex", color_values = sex_colors),
    list(name = "leukocyte_fraction", data = ss, color_var = "Leukocyte_Fraction", continuous = TRUE),
    list(name = "T", data = ss, color_var = "T", color_values = t_colors),
    list(name = "N", data = ss, color_var = "N", color_values = n_colors),
    list(name = "M", data = ss, color_var = "M", color_values = m_colors),
    list(name = "age", data = ss, color_var = "Age", continuous = TRUE),
    list(name = "mir200", data = ss, color_var = "MIR200C", continuous = TRUE),
    list(name = "methylation_cluster", data = ss, color_var = "Methylation_Clusters", color_values = cluster_colors),

    # ZEB plots (filtered dataset)
    list(name = "ZEB1", data = ss_filtered, color_var = "ZEB1", continuous = TRUE),
    list(name = "ZEB2", data = ss_filtered, color_var = "ZEB2", continuous = TRUE)
)

PLOT_SIZE <- list(width = 5.5, height = 4.5)

common_theme <- theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray"),
    aspect.ratio = 1,
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "none"
)

common_theme_with_legend <- theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray"),
    aspect.ratio = 1,
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "right",   # or "bottom" if you prefer
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm")
)

# ==========================
# TSNE WITH LABELS
# ==========================
create_tsne_plot_with_labels <- function(data, color_var, color_values = NULL, shape_var = NULL,
                                         shape_values = NULL, continuous = FALSE, title, legend_title = NULL) {
    p <- ggplot(data = data, aes(x = tSNE1, y = tSNE2)) +
        coord_fixed() +
        common_theme_with_legend
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
    if (is.null(legend_title)) {
        legend_title <- prettify_name(color_var)
    }
    p <- p + labs(title = title, color = legend_title)
    return(p)
}

for (config in plot_list) {
    plot_title <- paste("t-SNE plot colored by", prettify_name(config$name))
    legend_title <- prettify_name(config$color_var)
    plot <- create_tsne_plot_with_labels(
        data = config$data,
        color_var = config$color_var,
        color_values = config$color_values,
        shape_var = config$shape_var,
        shape_values = config$shape_values,
        continuous = if (!is.null(config$continuous)) config$continuous else FALSE,
        title = plot_title,
        legend_title = legend_title
    )
    filename <- file.path(FIG_DIR,
                          sprintf("tsne88_%s.pdf", config$name))
    ggsave(filename = filename,
           plot = plot,
           width = PLOT_SIZE$width,
           height = PLOT_SIZE$height,
           units = "in",
           device = "pdf")
}

# ==========================
# TSNE WITHOUT LEGEND (figure creation)
# ==========================
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

for (config in plot_list) {
    plot <- create_tsne_plot(
        data = config$data,
        color_var = config$color_var,
        color_values = config$color_values,
        shape_var = config$shape_var,
        shape_values = config$shape_values,
        continuous = if (!is.null(config$continuous)) config$continuous else FALSE
    )

    filename <- file.path(FIG_DIR,
                          sprintf("tsne88_%s.pdf",
                                  config$name))

    ggsave(filename = filename,
           plot = plot,
           width = PLOT_SIZE$width,
           height = PLOT_SIZE$height,
           units = "in",
           device = "pdf")
}

