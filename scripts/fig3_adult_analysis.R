# ============================================================
#   This script contains the analysis for Fig. 3, including:
#   TCGA-THCA only t-SNE and differential methylation analysis
#
#   S3A
#   S3B
#   S3C
#   S3D
#   S3E
#   S3F
#
# ============================================================

source(here::here("functions", "load_packages.R"))
source(here::here("functions", "color_keys.R"))
source(here::here("functions", "functions.R"))

# ========================
# TSNE COORDINATES
# ========================
ss = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx")
betas <- readRDS(here("data", "20250329_thca496_betas_processed_exc_sex.rds")) #395737    496
betas <- betas[, colnames(betas) %in% ss$Sample_ID] #395737    496
betas <- betas[, match(ss$Sample_ID, colnames(betas))]

mtx = bSubMostVariable(betas, 3000)
mtx = t(mtx) #496 3000

set.seed(12345)
tsne = Rtsne(mtx, dims=2, perplexity=22)
df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
rownames(df) = rownames(mtx)
saveRDS(df, here("data", "20250329_thca496_tsne_coords_exc_sex.rds"))

# ========================
# TSNE PLOTS
# ========================
ss <- read_excel(here("ss", "20241023_thca_master.xlsx"))
df <- readRDS(here("data", "20250329_thca496_tsne_coords_exc_sex.rds"))
ss$tSNE1 <- df$tSNE1
ss$tSNE2 <- df$tSNE2
DATE <- "20250329"

# Define Colors
invasiveness_colors <- c(
    "High" = "#f8766dff",
    "Low" = "#00bfc4ff"
)
driver_colors <- c(
    "BRAF V600E" = "#ff6cc3ff",
    "Kinase Fusion" = "#20bb20ff",
    "Ras-like" = "#00b4f0ff",
    "DICER1" = "#b47cffff"
)

# Define consistent plot dimensions
plot_width <- 6  # Increased width to accommodate legend
plot_height <- 5

# Function to create consistent plots
create_tsne_plot <- function(data, color_var, color_mapping = NULL, gradient_palette = NULL,
                             title = NULL, is_numeric = FALSE) {
    p <- ggplot() +
        geom_point(data = data, aes(x = tSNE1, y = tSNE2, color = {{color_var}}), na.rm = TRUE) +
        theme_minimal() +
        theme(legend.position = "right",
              plot.margin = margin(10, 10, 10, 10))

    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }

    if (is_numeric) {
        # For numeric variables, use gradient
        p <- p + scale_color_gradientn(colors = gradient_palette, na.value = "lightgray")
    } else if (!is.null(color_mapping)) {
        # For categorical variables, use manual colors
        p <- p + scale_color_manual(values = color_mapping, na.value = "lightgray")
    }

    return(p)
}

# Invasiveness Plot (With Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_invasiveness.pdf")),
    width = plot_width, height = plot_height, onefile = FALSE)
p <- create_tsne_plot(ss, color_var = Invasiveness, color_mapping = invasiveness_colors,
                      title = "Invasiveness")
plot(p)
dev.off()

# Sex Plot (With Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_sex.pdf")),
    width = plot_width, height = plot_height, onefile = FALSE)
p <- create_tsne_plot(ss, color_var = as.factor(Sex),
                      title = "Sex")
plot(p)
dev.off()

# Driver Group Plot (With Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_driver.pdf")),
    width = plot_width, height = plot_height, onefile = FALSE)
p <- create_tsne_plot(ss, color_var = Driver_Group, color_mapping = driver_colors,
                      title = "Driver Group")
plot(p)
dev.off()

# Differentiation Score Plot (With Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_differentiation.pdf")),
    width = plot_width, height = plot_height, onefile = FALSE)
ss$Differentiation_Score <- as.numeric(ss$Differentiation_Score)
p <- create_tsne_plot(ss, color_var = Differentiation_Score,
                      gradient_palette = brewer.pal(9, "BuPu"),
                      title = "Differentiation Score", is_numeric = TRUE)
plot(p)
dev.off()

# ERK Score Plot (With Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_erk.pdf")),
    width = plot_width, height = plot_height, onefile = FALSE)
ss$ERK_Score <- as.numeric(ss$ERK_Score)
p <- create_tsne_plot(ss, color_var = ERK_Score,
                      gradient_palette = brewer.pal(9, "BuGn"),
                      title = "ERK Score", is_numeric = TRUE)
plot(p)
dev.off()

# BRS Score Plot (With Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_brs.pdf")),
    width = plot_width, height = plot_height, onefile = FALSE)
ss$BRS_Score <- as.numeric(ss$BRS_Score)
p <- create_tsne_plot(ss, color_var = BRS_Score,
                      gradient_palette = brewer.pal(9, "Blues"),
                      title = "BRS Score", is_numeric = TRUE)
plot(p)
dev.off()

# Clusters Plot (With Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_clusters.pdf")),
    width = plot_width, height = plot_height, onefile = FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.margin = margin(10, 10, 10, 10),
          legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray"))
plot(p)
dev.off()

# ========================
# TSNE PLOTS W/O LEGEND
# ========================

# Invasiveness Plot (No Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_invasiveness_no_legend.pdf")), width=5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Invasiveness), na.rm = TRUE) +
    scale_color_manual(values = invasiveness_colors, na.value = "lightgray") +  # Change NA color
    theme_minimal() +
    theme(legend.position = "none")
plot(p)
dev.off()

# Sex Plot (No Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_sex_no_legend.pdf")), width=5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = as.factor(Sex)), na.rm = TRUE) +
    # scale_color_manual(values = invasiveness_colors, na.value = "lightgray") +  # Change NA color
    theme_minimal() +
    theme(legend.position = "none")
plot(p)
dev.off()

# Driver Group Plot (No Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_driver_no_legend.pdf")), width=5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Driver_Group), na.rm = TRUE) +
    scale_color_manual(values = driver_colors, na.value = "lightgray") +  # Change NA color
    theme_minimal() +
    theme(legend.position = "none")
plot(p)
dev.off()

# Differentiation Score Plot (No Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_differentiation_no_legend.pdf")), width=5, height=5, onefile=FALSE)
ss$Differentiation_Score <- as.numeric(ss$Differentiation_Score)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Differentiation_Score), na.rm = TRUE) +
    scale_color_gradientn(colors = brewer.pal(9, "BuPu"), na.value = "lightgray") +  # Adjust NA color
    theme_minimal() +
    theme(legend.position = "none")
plot(p)
dev.off()

# ERK Score Plot (No Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_erk_no_legend.pdf")), width=5, height=5, onefile=FALSE)
ss$ERK_Score <- as.numeric(ss$ERK_Score)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = ERK_Score), na.rm = TRUE) +
    scale_color_gradientn(colors = brewer.pal(9, "BuGn"), na.value = "lightgray") +  # Adjust NA color
    theme_minimal() +
    theme(legend.position = "none")
plot(p)
dev.off()

# BRS Score Plot (No Legend)
pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_brs_no_legend.pdf")), width=5, height=5, onefile=FALSE)
ss$BRS_Score <- as.numeric(ss$BRS_Score)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = BRS_Score), na.rm = TRUE) +
    scale_color_gradientn(colors = brewer.pal(9, "Blues"), na.value = "lightgray") +  # Adjust NA color
    theme_minimal() +
    theme(legend.position = "none")
plot(p)
dev.off()

# test clusters - correct
# pdf(file.path(here(), "figures", paste0(DATE, "_thca496_tsne_clusters.pdf")), width=5, height=5, onefile=FALSE)
# p <- ggplot() +
#     geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
#     theme(legend.title = element_blank(),
#           panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
#           panel.grid.major = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgray")
#     )
# plot(p)
# dev.off()

# ========================
# DM INVASIVENESS
# ========================
# ss = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx")
# ss = filter(ss, ss$Cluster != "NA")
#     # %>% filter(!(Invasiveness == "NA"))
# betas <- readRDS(here("data", "20250204_thca_betas_processed.rds"))
# betas <- betas[, colnames(betas) %in% ss$Sample_ID] #395737    449
#
# se = SummarizedExperiment(betas, colData = ss)
#
# se_ok = checkLevels(assay(se), colData(se)$Invasiveness)
# colData(se)$Invasiveness <- relevel(factor(colData(se)$Invasiveness), "Low")
# se = se[se_ok,] #395737    496
# smry = DML(se, ~Invasiveness)
# res = summaryExtractTest(smry)
# head(res)
#
# saveRDS(smry, here("diff_meth", "20250204_thca_smry_invasiveness.rds"))
# saveRDS(res, here("diff_meth", "20250204_thca_res_invasiveness.rds"))

# ========================
# DM INVASIVENESS BH ADJ
# ========================
res <- readRDS(here("diff_meth", "20250204_thca_res_invasiveness.rds"))
res$Pval_InvasivenessHigh <- p.adjust(res$Pval_InvasivenessHigh, method = "BH")
saveRDS(res, file.path(here(), "diff_meth", "20250324_thca_res_invasiveness_adj.rds"))

# ========================
# DM INVASIVENESS ENRICHMENT
# ========================
DATE = "20250325"
res <- readRDS(file.path(here(), "diff_meth", "20250324_thca_res_invasiveness_adj.rds"))

x <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh > 0.2], platform = "HM450", universe=res$Probe_ID)
y <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh < -0.2], platform = "HM450",  universe=res$Probe_ID)

# HIGH V LOW ALL ENRICHMENT
pdf(here("figures", paste0(DATE, "_thca_invasiveness_covar_all_hyper.pdf")),
    width = 8, height = 5,  onefile = FALSE)
KYCG_plotEnrichAll(x)
dev.off()

pdf(here("figures", paste0(DATE, "_thca_invasiveness_covar_all_hypo.pdf")),
    width = 8, height = 5,  onefile = FALSE)
KYCG_plotEnrichAll(y)
dev.off()

# HIGH V LOW INVASIVE TFBS HYPER
# a <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh > 0.2], "TFBS", platform="HM450", universe=res$Probe_ID)
# pdf(here("figures", paste0(DATE, "_thca_invasiveness_TFBS_enrichment_hyper.pdf")),  width = 4, height = 5, onefile = FALSE)
# plotDotFDR(a, n_max = 20, n_min = 20)
# dev.off()

# HIGH V LOW INVASIVE TFBS HYPO
b <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh < -0.2], "TFBS", platform="HM450", universe=res$Probe_ID)
pdf(here("figures", paste0(DATE, "_thca_invasiveness_TFBS_enrichment_hypo.pdf")), width = 4, height = 5,  onefile = FALSE)
plotDotFDR(b, n_max = 20, n_min = 20)
dev.off()

# HIGH V LOW INVASIVE TISSUE HYER
# c <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh > 0.2], "tissueSignature", platform="HM450", universe=res$Probe_ID)
# pdf(here("figures", paste0(DATE, "_thca_invasiveness_tissue_enrichment_hyper.pdf")), width = 4, height = 5,  onefile = FALSE)
# plotDotFDR(c, n_max = 20, n_min = 20)
# dev.off()

# HIGH V LOW INVASIVE TISSUE HYPO
d <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh < -0.2], "tissueSignature", platform="HM450", universe=res$Probe_ID)
pdf(here("figures", paste0(DATE, "_thca_invasiveness_tissue_enrichment_hypo.pdf")), width = 4, height = 5,  onefile = FALSE)
plotDotFDR(d, n_max = 20, n_min = 20)
dev.off()

# ========================
# VOLCANO PLOT INVASIVENESS COVAR
# ========================
res <- readRDS(file.path(here(), "diff_meth", "20250324_thca_res_invasiveness_adj.rds"))

res$threshold <- "NS"
res$threshold[res$Pval_InvasivenessHigh < 0.05 & res$Est_InvasivenessHigh > 0.2] <- "Up"
res$threshold[res$Pval_InvasivenessHigh < 0.05 & res$Est_InvasivenessHigh < -0.2] <- "Down"

# length(res$threshold[res$Pval_InvasivenessHigh < 0.05 & res$Est_InvasivenessHigh > 0.2]) #0
# length(res$threshold[res$Pval_InvasivenessHigh < 0.05 & res$Est_InvasivenessHigh < -0.2]) #142
#
# length(res$threshold[res$Est_InvasivenessHigh > 0.2]) #0
# length(res$threshold[res$Est_InvasivenessHigh < -0.2]) #142

sig_points <- subset(res, threshold != "NS")
nonsig_points <- subset(res, threshold == "NS")

p <- ggplot() +
    geom_bin2d(data = nonsig_points,
               aes(x = Est_InvasivenessHigh, y = -log10(Pval_InvasivenessHigh)),
               bins = 250) +  # Reduced number of bins to make each bin more visible
    scale_fill_gradient(low = "lightgray",
                        high = "black",  # Darker grey for better contrast
                        trans = "log10") + # Add log transform to better show density differences
    geom_point(data = sig_points,
               aes(x = Est_InvasivenessHigh, y = -log10(Pval_InvasivenessHigh),
                   color = threshold),
               size = 0.1,
               alpha = 0.3) +
    scale_color_manual(values = c(
        "Down" = "#0c66bcff",
        "Up" = "#d61525ff"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "#575757ff") +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "#575757ff") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#575757ff") +
    theme_minimal() +
    theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey90"),
        legend.position = "right",
        legend.title = element_blank()
    ) +
    labs(
        x = "Estimate",
        y = "-log10(FDR)"
    ) +
    guides(fill = "none")

pdf(here("figures", paste0(DATE, "_thca_invasiveness_covar_adj_volcano.pdf")), width=5, height=5, onefile=FALSE)
plot(p)
dev.off()









