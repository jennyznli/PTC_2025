# ============================================================
#   This script contains the analysis for Fig. 2, including:
#   Differential methylation between high v. low Clinical_Invasiveness
#
#   2A Volcano plot high v. low Clinical_Invasiveness
#   2B Heatmap 88 primary tumor samples
#   2C Dot plot
#   2D Dot plot
#
#   S2A Heatmap 111 LN & primary tumor samples
#   S2B DMP composition
#   S2C Barplot global methylation of cluster
#   S2D Barplot global methylation of driver groups
#   S2E Dot plot
#   S2F Dot plot
#   S2G Dot plot

# ============================================================
# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- here()
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")
FIG_DIR <- file.path(BASE_DIR, "figures", "final")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# 2A. INVASIVENESS VOLCANO
# ========================
res <- readRDS(file.path(DATA_DIR, "ped88_dm_invasiveness_bh.rds"))
res$threshold <- "NS"
res$threshold[res$Pval_Clinical_InvasivenessHigh < 0.05 & res$Est_Clinical_InvasivenessHigh > 0.2] <- "Up"
res$threshold[res$Pval_Clinical_InvasivenessHigh < 0.05 & res$Est_Clinical_InvasivenessHigh < -0.2] <- "Down"

sig_points <- subset(res, threshold != "NS")
nonsig_points <- subset(res, threshold == "NS")

p <- ggplot() +
    geom_bin2d(data = nonsig_points,
               aes(x = Est_Clinical_InvasivenessHigh, y = -log10(Pval_Clinical_InvasivenessHigh)),
               bins = 250) +
    scale_fill_gradient(low = "lightgray",
                        high = "black",
                        trans = "log10") + # log transform to better show density differences
    geom_point(data = sig_points,
               aes(x = Est_Clinical_InvasivenessHigh, y = -log10(Pval_Clinical_InvasivenessHigh),
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

pdf(file.path(FIG_DIR, "ped88_invasiveness_volcano.pdf"), width=5, height=5, onefile=FALSE)
plot(p)
dev.off()

# ========================
# 2B. HEATMAP 88
# ========================
betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_primary <- ss %>% filter(Lymph_Node == "F")
res <- readRDS(file.path(DATA_DIR, "ped88_dm_invasiveness_bh.rds"))

# Filter for significant probes
sig_probes <- res %>%
    dplyr::filter(Pval_Clinical_InvasivenessHigh < 0.05,
                  abs(Est_Clinical_InvasivenessHigh) > 0.2) %>%
    mutate(Methylation = ifelse(Est_Clinical_InvasivenessHigh > 0.2, "Hyper", "Hypo"))

# Create annotations
annotation_row <- data.frame(
    Methylation = factor(sig_probes$Methylation),
    row.names = sig_probes$Probe_ID
)

annotation_col <- data.frame(
    Clinical_Invasiveness = factor(ss_primary$Clinical_Invasiveness),
    Methylation_Clusters = factor(ss_primary$Methylation_Clusters),
    row.names = ss_primary$Sample_ID
) %>% arrange(Clinical_Invasiveness)

colnames(betas) <- ss_primary$Sample_ID[match(colnames(betas), ss_primary$IDAT)]
sel_betas <- betas[sig_probes$Probe_ID, rownames(annotation_col) ]

annotation_colors <- list(
    Methylation = methylation_colors,
    Clinical_Invasiveness = invasiveness_colors,
    Methylation_Clusters = cluster_colors
)

pdf(file.path(FIG_DIR, "ped88_dm_invasiveness_heatmap.pdf"), width = 13, height = 12, onefile = FALSE)
pheatmap(sel_betas,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#0070ffff", "white", "#f70000ff"))(50),
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_legend = FALSE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         legend = FALSE)
dev.off()

# ========================
# S2A. HEATMAP 98
# ========================
betas <- readRDS(file.path(DATA_DIR, "ped_betas_imputed.rds"))
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
res <- readRDS(file.path(DATA_DIR, "ped88_dm_invasiveness_bh.rds"))

sig_probes <- res %>%
    dplyr::filter(Pval_Clinical_InvasivenessHigh < 0.05,
                  abs(Est_Clinical_InvasivenessHigh) > 0.2) %>%
    mutate(Methylation = ifelse(Est_Clinical_InvasivenessHigh > 0.2, "Hyper", "Hypo"))

colnames(betas) <- ss$Sample_ID[match(colnames(betas), ss$IDAT)]
sel_betas <- betas[sig_probes$Probe_ID, rownames(annotation_col)]

# Create annotations
annotation_row <- data.frame(
    Methylation = factor(sig_probes$Methylation),
    row.names = sig_probes$Probe_ID
)

annotation_col <- data.frame(
    Clinical_Invasiveness = factor(ss$Clinical_Invasiveness),
    Cluster_Group = factor(ss$Methylation_Clusters),
    Lymph_Node = factor(ss$Lymph_Node),
    row.names = ss$Sample_ID
) %>% arrange(Clinical_Invasiveness)

annotation_colors <- list(
    Cluster_Group = cluster_colors,
    Clinical_Invasiveness = Clinical_Invasiveness_colors,
    Lymph_Node = lymph_node_colors,
    Methylation = methylation_colors
)

pdf(file.path(FIG_DIR, "ped98_dm_invasiveness_heatmap.pdf"), width = 15, height = 12, onefile = FALSE)
pheatmap(sel_betas,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#0070ffff", "white", "#f70000ff"))(50),
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_legend = FALSE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         legend = FALSE)
dev.off()

# ========================
# S2B. TOTAL DMP BARPLOT
# ========================
res1 <- readRDS(file.path(DATA_DIR, "ped88_dm_cluster_li_bh.rds"))
res2 <- readRDS(file.path(DATA_DIR, "ped88_dm_cluster_hil_bh.rds"))

# LI vs HI (GroupHI in res1)
cat(dim(filter(res1, BH_GroupHI < 0.05, Est_Methylation_ClustersHI > 0.2))[1], "hyper (LI vs HI)\n")
cat(dim(filter(res1, BH_GroupHI < 0.05, Est_Methylation_ClustersHI < -0.2))[1], "hypo (LI vs HI)\n\n")

# LI vs HIL (GroupHIL in res1)
cat(dim(filter(res1, BH_GroupHIL < 0.05, Est_Methylation_ClustersHIL > 0.2))[1], "hyper (LI vs HIL)\n")
cat(dim(filter(res1, BH_GroupHIL < 0.05, Est_Methylation_ClustersHIL < -0.2))[1], "hypo (LI vs HIL)\n\n")

# HIL vs HI (GroupHI in res2)
cat(dim(filter(res2, BH_GroupHI < 0.05, Est_Methylation_ClustersHI > 0.2))[1], "hyper (HIL vs HI)\n")
cat(dim(filter(res2, BH_GroupHI < 0.05, Est_Methylation_ClustersHI < -0.2))[1], "hypo (HIL vs HI)\n")

get_sig_probes <- function(res, bh_col, est_col) {
    res %>%
        dplyr::filter(!!sym(bh_col) < 0.05,
                      abs(!!sym(est_col)) > 0.2) %>%
        mutate(Methylation = ifelse(!!sym(est_col) > 0.2, "Hyper", "Hypo"))
}

sig_results <- list( # will have to fix the name of columns
    "LI / HIL" = get_sig_probes(res1, "BH_GroupHIL", "Est_Methylation_ClustersHIL"),
    "LI / HI" = get_sig_probes(res1, "BH_GroupHI", "Est_Methylation_ClustersHI"),
    "HIL / HI" = get_sig_probes(res2, "BH_GroupHI", "Est_Methylation_ClustersHI")
)

data <- bind_rows(
    lapply(names(sig_results), function(group) {
        data.frame(
            Group = group,
            table(sig_results[[group]]$Methylation)
        ) %>%
            rename(DMPs = Freq) %>%
            arrange(desc(Var1))
    })
)
total_dmps <- data.frame(
    Group = names(sig_results),
    Total = sapply(sig_results, nrow),
    stringsAsFactors = FALSE
) %>%
    mutate(Label = format(Total, big.mark = ","))

data$Methylation <- factor(data$Var1, levels = c("Hyper", "Hypo"))

pdf(file.path(FIG_DIR, "ped88_dm_cluster_composition.pdf"), width=4, height=4, onefile=FALSE)
p <- ggplot(data, aes(x = Group, y = DMPs/10000)) +
    geom_bar(aes(fill = Var1), stat = "identity", position = "stack") +
    scale_fill_manual(values = methylation_colors) +
    geom_text(data = total_dmps,
              aes(y = Total/10000, label = Label),
              vjust = -0.3,
              size = 3.5) +
    scale_y_continuous(name = expression("DMPs [" * 10^4 * "]"),
                       limits = c(0, max(total_dmps$Total/10000) * 1.1),
                       breaks = seq(0, ceiling(max(total_dmps$Total/10000)), by = 2)) +
    theme_minimal() +
    theme(legend.title = element_text(size = 10),
          panel.grid.major.x = element_blank()) +
    labs(fill = "Methylation",
         x = "Cluster Comparison",
         y = expression("DMPs [" * 10^4 * "]"))
plot(p)
dev.off()

# ========================
# S2C. GLOBAL MEANS CLUSTERS
# ========================
betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_primary <- ss %>% filter(Lymph_Node == "F")

ss_primary$Global_Means = colMeans(betas)

wcc <- pairwise.wilcox.test(ss_primary$Global_Means, ss_primary$Methylation_Clusters,
                            p.adjust.method = "BH")
print(wcc)
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction
#
# data:  ss_primary$Global_Means and ss_primary$Methylation_Clusters
#
# HI      HIL
# HIL 0.00045 -
#     LI  0.00012 0.57581

pdf(file.path(FIG_DIR, "ped88_cluster_global_means_wc.pdf"), width=4.5, height=4.5, onefile=FALSE)
p <- ggplot(ss_primary, aes(x = Methylation_Clusters, y = Global_Means, fill = Methylation_Clusters)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
    stat_compare_means(method = "wilcox.test",
                       comparisons = list(c("HI", "HIL"),
                                          c("HI", "LI")),
                       label = "p") +
    scale_fill_manual(values = cluster_colors) +
    labs(x = "Cluster",
         y = "Global Mean Methylation",
         fill = "Cluster") +
    theme_minimal()
plot(p)
dev.off()

# ========================
# S2D. GLOBAL MEANS DRIVERS
# ========================
betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_primary <- ss %>% filter(Lymph_Node == "F")

ss_primary$Global_Means = colMeans(betas)

wcc <- pairwise.wilcox.test(ss_primary$Global_Means, ss_primary$Driver_Group,
                            p.adjust.method = "BH")
print(wcc)
# Pairwise comparisons using Wilcoxon rank sum exact test
#
# data:  ss_primary$Global_Means and ss_primary$Driver_Group
#
# BRAF V600E DICER1  Indeterminate Kinase Fusion
# DICER1        0.98780    -       -             -
#     Indeterminate 0.03971    0.14141 -             -
#     Kinase Fusion 0.98780    0.98780 0.03971       -
#     Ras-like      0.00091    0.06550 0.98780       0.00091
#
# P value adjustment method: BH

# Filter out Indeterminate category
ss_filtered <- ss_primary %>%
    filter(Driver_Group != "Indeterminate")

pdf(file.path(FIG_DIR, "ped88_driver_global_means_wc.pdf"), width=5.5, height=5, onefile=FALSE)
p <- ggplot(ss_filtered, aes(x = Driver_Group, y = Global_Means, fill = Driver_Group)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
    stat_compare_means(method = "wilcox.test",
                       comparisons = list(c("BRAF V600E", "Ras-like"),
                                          c("Kinase Fusion", "Ras-like"),
                                          c("DICER1", "Ras-like")),
                       label = "p",
                       hide.ns = TRUE,
                       p.cutoff = 0.01) +
    scale_fill_manual(values = driver_colors) +
    labs(x = "Cluster",
         y = "Global Mean Methylation",
         fill = "Cluster") +
    theme_minimal()
plot(p)
dev.off()

# ========================
# 2D. ENRICHMENT - CLUSTERS, TISSUES
# ========================
res1 <- readRDS(file.path(DATA_DIR, "ped88_dm_cluster_li_bh.rds"))
res2 <- readRDS(file.path(DATA_DIR, "ped88_dm_cluster_hil_bh.rds"))

# HIL vs HI — TISSUE HYPER
c <- testEnrichment(
    res2$Probe_ID[res2$Est_Methylation_ClustersHI > 0.2],
    "tissueSignature", platform = "EPIC",
    universe = res2$Probe_ID
)
pdf(file.path(FIG_DIR, "hil_vs_hi_hyper_tissue_enrichment.pdf"), width = 4.75, height = 2, onefile = FALSE)
plotDotFDR(c, n_max = 10, n_min = 10)
dev.off()

# HIL vs HI — TISSUE HYPO
d <- testEnrichment(
    res2$Probe_ID[res2$Est_Methylation_ClustersHI < -0.2],
    "tissueSignature", platform = "EPIC",
    universe = res2$Probe_ID
)
pdf(file.path(FIG_DIR, "hil_vs_hi_hypo_tissue_enrichment.pdf"), width = 3.5, height = 2, onefile = FALSE)
plotDotFDR(d, n_max = 10, n_min = 10)
dev.off()

# HIL vs LI — TISSUE HYPER
a <- testEnrichment(
    res2$Probe_ID[res2$Est_Methylation_ClustersLI > 0.2],
    "tissueSignature", platform = "EPIC",
    universe = res2$Probe_ID
)
pdf(file.path(FIG_DIR, "hil_vs_li_hyper_tissue_enrichment.pdf"), width = 4.75, height = 2, onefile = FALSE)
plotDotFDR(a, n_max = 10, n_min = 10)
dev.off()

# HIL vs LI — TISSUE HYPO
b <- testEnrichment(
    res2$Probe_ID[res2$Est_Methylation_ClustersLI < -0.2],
    "tissueSignature", platform = "EPIC",
    universe = res2$Probe_ID
)
pdf(file.path(FIG_DIR, "hil_vs_li_hypo_tissue_enrichment.pdf"), width = 4.75, height = 2, onefile = FALSE)
plotDotFDR(b, n_max = 10, n_min = 10)
dev.off()

# ========================
# S2F,G. ENRICHMENT - CLUSTERS, TFBS
# ========================
# HIL vs LI — TFBS HYPER
a <- testEnrichment(
    res2$Probe_ID[res2$Est_Methylation_ClustersLI > 0.2],
    "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path(FIG_DIR, "hil_vs_li_hyper_tfbs_enrichment.pdf"), width = 2.8, height = 3.9, onefile = FALSE)
plotDotFDR(a, n_min = 30, n_max = 30)
dev.off()

# HIL vs HI — TFBS HYPER
b <- testEnrichment(
    res2$Probe_ID[res2$Est_Methylation_ClustersHI > 0.2],
    "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path(FIG_DIR, "hil_vs_hi_hyper_tfbs_enrichment.pdf"), width = 2.8, height = 3.9, onefile = FALSE)
plotDotFDR(b, n_min = 30, n_max = 30)
dev.off()

# HIL vs LI — TFBS HYPO
c <- testEnrichment(
    res2$Probe_ID[res2$Est_Methylation_ClustersLI < -0.2],
    "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path(FIG_DIR, "hil_vs_li_hypo_tfbs_enrichment.pdf"), width = 2.8, height = 3.9, onefile = FALSE)
plotDotFDR(c, n_min = 30, n_max = 30)
dev.off()

# HIL vs HI — TFBS HYPO
d <- testEnrichment(
    res2$Probe_ID[res2$Est_Methylation_ClustersHI < -0.2],
    "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path(FIG_DIR, "hil_vs_hi_hypo_tfbs_enrichment.pdf"), width = 2.8, height = 3.9, onefile = FALSE)
plotDotFDR(d, n_min = 30, n_max = 30)
dev.off()

# ========================
# 2D-E, S2E-G. ENRICHMENT PLOTS - INVASIVENESS
# ========================
res <- readRDS(file.path(DATA_DIR, "ped88_dm_invasiveness_bh.rds"))

# LOW vs HIGH — TFBS HYPER
x <- testEnrichment(
    res$Probe_ID[res$Est_Clinical_InvasivenessHigh > 0.2],
    "TFBSconsensus", platform = "EPIC", universe = res$Probe_ID
)
pdf(file.path(FIG_DIR, "low_vs_high_hyper_tfbs_enrichment.pdf"),
    width = 2.8, height = 3.9, onefile = FALSE)
plotDotFDR(x, n_min = 30, n_max = 30)
dev.off()

# LOW vs HIGH — TFBS HYPO
y <- testEnrichment(
    res$Probe_ID[res$Est_Clinical_InvasivenessHigh < -0.2],
    "TFBSconsensus", platform = "EPIC", universe = res$Probe_ID
)
pdf(file.path(FIG_DIR, "low_vs_high_hypo_tfbs_enrichment.pdf"),
    width = 2.8, height = 3.9, onefile = FALSE)
plotDotFDR(y, n_min = 30, n_max = 30)
dev.off()


