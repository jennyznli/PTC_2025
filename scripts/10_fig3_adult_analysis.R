# ============================================================
#   FIG S3.
#   This script contains the adult only analysis:
#   TCGA-THCA t-SNE and differential methylation analysis
#   S3A-E. ERK, Differentiation, BRS, Invasiveness, Driver group
#   S3F. DMP enrichment dot plot
# ============================================================

# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- "/Users/jennyzli/Documents/HPC_share/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")
FIG_DIR <- file.path(BASE_DIR, "figures")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# LOADING DATA
# ========================
ss <- read_excel(file.path(SS_DIR, "adult_master.xlsx")) %>%
    filter(Clinical_Invasiveness != "NA")
betas <- readRDS(file.path(DATA_DIR, "adult_betas_imputed.rds"))
common <- intersect(ss$Sample_ID, colnames(betas))
betas <- betas[, common, drop = FALSE]

# ========================
# TSNE COORDINATES
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

set.seed(123)
tsne <- Rtsne(
    pca_data,
    dims = 2,
    perplexity = 22,
    max_iter = 2000,
    pca = FALSE
)

df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
ss$tSNE1 <- as.numeric(df$tSNE1)
ss$tSNE2 <- as.numeric(df$tSNE2)

rownames(df) <- ss$Sample_ID
saveRDS(df, file.path(DATA_DIR, "adult_tsne_coords.rds"))

ggplot(df, aes(x = tSNE1, y = tSNE2)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_minimal() +
    labs(title = "t-SNE plot of samples",
         x = "t-SNE 1",
         y = "t-SNE 2")

# ========================
# 3A-E. TSNE PLOTTING
# ========================
PLOT_WIDTH <- 6
PLOT_HEIGHT <- 5

ss$Differentiation_Score <- as.numeric(ss$Differentiation_Score)
ss$ERK_Score <- as.numeric(ss$ERK_Score)
ss$BRS_Score <- as.numeric(ss$BRS_Score)

plot_list <- list(
    list(var = "Clinical_Invasiveness",  name = "invasiveness",   title = "Clinical Invasiveness", is_numeric = FALSE, color_map = invasiveness_colors),
    list(var = "Sex",                    name = "sex",            title = "Sex",                   is_numeric = FALSE, color_map = sex_colors),
    list(var = "Driver_Group",           name = "driver",         title = "Driver Group",          is_numeric = FALSE, color_map = driver_colors),
    list(var = "Differentiation_Score",  name = "differentiation",title = "Differentiation Score", is_numeric = TRUE,  palette = brewer.pal(9, "BuPu")),
    list(var = "ERK_Score",              name = "erk",            title = "ERK Score",             is_numeric = TRUE,  palette = brewer.pal(9, "BuGn")),
    list(var = "BRS_Score",              name = "brs",            title = "BRS Score",             is_numeric = TRUE,  palette = brewer.pal(9, "Blues"))
)

for (plt in plot_list) {
    file_path <- file.path(FIG_DIR, paste0("adult_tsne_", plt$name, ".pdf"))
    pdf(file_path, width = PLOT_WIDTH, height = PLOT_HEIGHT, onefile = FALSE)

    p <- ggplot(ss, aes(x = tSNE1, y = tSNE2, color = .data[[plt$var]])) +
        geom_point(na.rm = TRUE) +
        theme_minimal() +
        ggtitle(plt$title) +
        labs(title = paste("t-SNE plot colored by", tools::toTitleCase(gsub("_", " ", plt$var))))

    if (plt$is_numeric) {
        p <- p + scale_color_gradientn(colors = plt$palette, na.value = "lightgray")
    } else if (!is.null(plt$color_map)) {
        p <- p + scale_color_manual(values = plt$color_map, na.value = "lightgray")
    }

    print(p)
    dev.off()
}

# ========================
# W/0 LEGENDS (figure creation)
# ========================
PLOT_WIDTH <- 5
PLOT_HEIGHT <- 5

for (plt in plot_list) {
    file_path <- file.path(FIG_DIR, paste0("adult_tsne_", plt$name, ".pdf"))
    pdf(file_path, width = PLOT_WIDTH, height = PLOT_HEIGHT, onefile = FALSE)

    p <- ggplot(ss, aes(x = tSNE1, y = tSNE2, color = .data[[plt$var]])) +
        geom_point(na.rm = TRUE) +
        theme_minimal() +
        theme(
            legend.position = "none",
            plot.title = element_blank()
        )

    if (plt$is_numeric) {
        p <- p + scale_color_gradientn(colors = plt$palette, na.value = "lightgray")
    } else if (!is.null(plt$color_map)) {
        p <- p + scale_color_manual(values = plt$color_map, na.value = "lightgray")
    }

    print(p)
    dev.off()
}

# ========================
# DIFF METH - INVASIVENESS
# ========================
ss <- read_excel(file.path(SS_DIR, "adult_master.xlsx")) %>%
    filter(Clinical_Invasiveness != "NA")

ss$Age <- as.numeric(ss$Age)
ss$Sex <- as.factor(ss$Sex)
ss$Clinical_Invasiveness <- as.factor(ss$Clinical_Invasiveness)

se <- SummarizedExperiment(betas, colData = ss)

se_ok <- (checkLevels(assay(se), colData(se)$Sex) &
              checkLevels(assay(se), colData(se)$Clinical_Invasiveness))

cat("Probes passing QC:", sum(se_ok), "out of", length(se_ok), "\n")
se <- se[se_ok, ]
cat("Final SE dimensions:", dim(se), "\n")

colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
colData(se)$Clinical_Invasiveness <- relevel(factor(colData(se)$Clinical_Invasiveness), "Low")

cat("Reference levels - Sex:", levels(colData(se)$Sex)[1],
    "Clinical_Invasiveness:", levels(colData(se)$Clinical_Invasiveness)[1], "\n")

smry <- DML(se, ~Clinical_Invasiveness + Sex + Age,
            BPPARAM = BiocParallel::MulticoreParam(workers = N_WORKERS))
res_inv <- summaryExtractTest(smry)

cat("Results dimensions:", dim(res_inv), "\n")
cat("First few results:\n")
print(head(res_inv, 3))

saveRDS(res_inv, file = file.path(DATA_DIR, "adult_dm_invasiveness.rds"))
cat("Invasiveness differential methylation completed:", dim(res_inv), "\n")

# ========================
# S3F. ENRICHMENT PLOT
# ========================
res <- readRDS(file.path(DATA_DIR, "adult_dm_invasiveness.rds"))
res$BH_Clinical_InvasivenessHigh <- p.adjust(res$Pval_Clinical_InvasivenessHigh, method = "BH")
saveRDS(res, file.path(DATA_DIR, "adult_dm_invasiveness_bh.rds"))

# HIGH V LOW INVASIVE TFBS HYPO
x <- testEnrichment(res$Probe_ID[res$Est_Clinical_InvasivenessHigh < -0.2], "TFBS", platform = "HM450", universe = res$Probe_ID)
pdf(file.path(FIG_DIR, "adult_low_v_high_hypo_tfbs_enrichment.pdf"),
    width = 2.8, height = 3.9, onefile = FALSE)
plotDotFDR(x, n_max = 30, n_min = 30)
dev.off()

# HIGH V LOW INVASIVE TFBS HYPER
# a <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh > 0.2], "TFBS", platform = "HM450", universe = res$Probe_ID)
# pdf(file.path("figures", "thca_invasiveness_TFBS_enrichment_hyper.pdf"),
#     width = 4, height = 5, onefile = FALSE)
# plotDotFDR(a, n_max = 20, n_min = 20)
# dev.off()


