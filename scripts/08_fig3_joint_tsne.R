# ============================================================
#   FIG 3.
#   This script contains the analysis for the
#   joint t-SNE embedding of adult and pediatric cohorts:
#   3-D Invasiveness, Driver Group, Leukocyte Fraction, Sex
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
betas <- readRDS(file.path(DATA_DIR, "joint_betas_imputed.rds"))
ss <- read_xlsx(file.path(SS_DIR, "joint_master.xlsx"))

# ========================
# JOINT SAMPLESHEET CREATION
# ========================
ss_ped <- read_xlsx(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_adult <- read_xlsx(file.path(SS_DIR, "adult_master.xlsx"))

ss_ped$'T' = str_extract(ss_ped$TNM, "T[^NMR]+")
ss_ped$'N' = str_extract(ss_ped$TNM, "N[^TMR]+")
ss_ped$'M' = str_extract(ss_ped$TNM, "M[^TNR]+")

sel_cols_ped <- c("Source", "IDAT", "Methylation_Clusters", "Clinical_Invasiveness", "Driver_Group", "Sex","Age", "N", "M")
sel_cols_adult <- c("Source", "Sample_ID", "Methylation_Clusters", "Clinical_Invasiveness", "Driver_Group", "Sex","Age", "N", "M")

ss_ped <- ss_ped[, sel_cols_ped]
ss_adult <- ss_adult[, sel_cols_adult]
colnames(ss_ped) <- c("Source", "Sample_ID", "Methylation_Clusters", "Clinical_Invasiveness", "Driver_Group", "Sex","Age", "N", "M")
ss_joint <- rbind(ss_ped, ss_adult)
write.xlsx(ss_joint, file.path(SS_DIR, "joint_master.xlsx"))

# ========================
# TSNE COORDINATES
# ========================
ss <- read_xlsx(file.path(SS_DIR, "joint_master.xlsx"))
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

set.seed(123)
tsne <- Rtsne(
    pca_data,
    dims = 2,
    perplexity = 24,
    max_iter = 2000,
    pca = FALSE
)

df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
ss$tSNE1 <- as.numeric(df$tSNE1)
ss$tSNE2 <- as.numeric(df$tSNE2)

rownames(df) <- ss$Sample_ID
saveRDS(df, file.path(DATA_DIR, "joint_tsne_coords.rds"))

ggplot(df, aes(x = tSNE1, y = tSNE2)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_minimal() +
    labs(title = "t-SNE plot of samples",
         x = "t-SNE 1",
         y = "t-SNE 2")

# ========================
# 3A-D. TSNE WITH LABELS
# ========================
leuko <- read.csv(file.path(DATA_DIR, "joint_leuko.csv"), stringsAsFactors = FALSE)
indices <- match(ss$Sample_ID, leuko[,1])
leuko <- leuko[indices, , drop = FALSE]
ss$Leukocyte_Fraction <- as.numeric(leuko[,2])

mir <- read.csv(file.path(DATA_DIR, "joint_mir200.csv"))
indices <- match(ss$Sample_ID, mir[,1])
mir <- mir[indices, , drop = FALSE]
ss$MIR200C <- as.numeric(mir[,2])

plot_list <- list(
    list(name = "joint_tsne_invasiveness", color = "Clinical_Invasiveness", colors = invasiveness_colors),
    list(name = "joint_tsne_drivers", color = "Driver_Group", colors = driver_colors),
    list(name = "joint_tsne_sex", color = "Sex", colors = sex_colors),
    list(name = "joint_tsne_leuko", color = "Leukocyte_Fraction", gradient = TRUE),
    list(name = "joint_tsne_mir200c", color = "MIR200C", gradient = TRUE)
)

for (plt in plot_list) {
    file_name <- paste0(plt$name, ".pdf")
    pdf(file.path(FIG_DIR, file_name), width = 6, height = 5, onefile = FALSE)

    aes_args <- aes(x = tSNE1, y = tSNE2, color = !!sym(plt$color), shape = Source)
    if (!is.null(plt$tooltip) && plt$tooltip) {
        aes_args <- modifyList(aes_args, aes(text = paste("Sample ID:", Sample_ID)))
    }

    p <- ggplot() +
        geom_point(data = ss, mapping = aes_args) +
        scale_shape_manual(values = custom_shapes) +
        (if (!is.null(plt$gradient) && plt$gradient) {
            scale_color_gradientn(colors = parula(20))
        } else {
            scale_color_manual(values = plt$colors)
        }) +
        theme(
            legend.title = element_blank(),
            panel.background = element_rect(fill = "#ffffff", color = NA),
            panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray")
        ) +
        labs(title = paste("t-SNE plot colored by", tools::toTitleCase(gsub("_", " ", plt$color))))


    plot(p)
    dev.off()
}

# ========================
# NO LEGENDS (figure creation)
# ========================

for (plt in plot_list) {
    file_name <- paste0(plt$name, ".pdf")
    pdf(file.path(FIG_DIR, file_name), width = 5, height = 5, onefile = FALSE)

    aes_args <- aes(x = tSNE1, y = tSNE2, color = !!sym(plt$color), shape = Source)
    if (!is.null(plt$tooltip) && plt$tooltip) {
        aes_args <- modifyList(aes_args, aes(text = paste("Sample ID:", Sample_ID)))
    }

    p <- ggplot() +
        geom_point(data = ss, mapping = aes_args) +
        scale_shape_manual(values = custom_shapes) +
        (if (!is.null(plt$gradient) && plt$gradient) {
            scale_color_gradientn(colors = parula(20))
        } else {
            scale_color_manual(values = plt$colors)
        }) +
        theme(
            legend.position = "none",
            panel.background = element_rect(fill = "#ffffff", color = NA),
            panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray")
        )
    plot(p)
    dev.off()
}
