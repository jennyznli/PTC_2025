# ============================================================
#   This script performs the Fig. 2 RNA analysis and plotting.
# ============================================================
# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- here()
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# LOADING DATA
# ========================
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx")) %>%
    dplyr::filter(Lymph_Node == "F")
mtx <- read.csv(file.path(DATA_DIR, "rna_log2cpmfiltered.csv"), row.names = 1) %>%
    .[, colnames(.) %in% ss$Sample_ID] %>%
    unique()  # 22516 genes × 84 samples
ss <- ss %>%
    filter(Sample_ID %in% colnames(mtx)) #84 x 20
mtx <- mtx[,match(ss$Sample_ID, colnames(mtx))]

invasiveness <- factor(ss$Clinical_Invasiveness, levels = c("Low", "High"))

# ========================
# RNA-SEQ QC
# ========================
gene_means <- rowMeans(mtx, na.rm = TRUE)
gene_vars <- apply(mtx, 1, var, na.rm = TRUE)
gene_zeros <- rowSums(mtx == 0, na.rm = TRUE) / ncol(mtx)

# Filter criteria
low_expr_threshold <- quantile(gene_means, 0.1)  # Bottom 10% by mean
high_zero_threshold <- 0.8  # Remove genes with >80% zeros
low_var_threshold <- quantile(gene_vars, 0.05, na.rm = TRUE)  # Bottom 5% by variance

genes_to_keep <- gene_means > low_expr_threshold &
    gene_zeros < high_zero_threshold &
    gene_vars > low_var_threshold &
    !is.na(gene_vars)

cat("Filtering genes:\n")
cat("- Low expression (<", round(low_expr_threshold, 3), "):", sum(gene_means <= low_expr_threshold), "\n")
cat("- High zeros (>", high_zero_threshold*100, "%):", sum(gene_zeros >= high_zero_threshold), "\n")
cat("- Low variance:", sum(gene_vars <= low_var_threshold, na.rm = TRUE), "\n")
cat("- Genes kept:", sum(genes_to_keep), "of", nrow(mtx), "\n")

mtx <- mtx[genes_to_keep, ]
saveRDS(mtx, file.path(DATA_DIR, "rna_seq_filtered.rds"))

# ========================
# DIFFERENTIAL EXPRESSION
# ========================
design <- model.matrix(~ invasiveness)
colnames(design) <- c("Intercept", "InvasivenessHigh")

fit <- mtx %>%
    lmFit(design) %>%
    eBayes()

deg_all <- topTable(fit, coef="InvasivenessHigh", n=Inf, sort.by="P")
saveRDS(deg_all, file.path(DATA_DIR, "diff_exp_genes.rds"))

# Add basic statistics
deg_all$Significant <- deg_all$adj.P.Val < 0.05
deg_all$HighLogFC <- abs(deg_all$logFC) > 1
deg_all$Stringent <- deg_all$Significant & deg_all$HighLogFC

cat("Results summary:\n")
cat("- Total genes tested:", nrow(deg_all), "\n")
cat("- Significant (FDR < 0.05):", sum(deg_all$Significant), "\n")
cat("- High effect size (|logFC| > 1):", sum(deg_all$HighLogFC), "\n")
cat("- Stringent (both criteria):", sum(deg_all$Stringent), "\n")

# Direction of significant changes
if (sum(deg_all$Significant) > 0) {
    sig_up <- sum(deg_all$Significant & deg_all$logFC > 0)
    sig_down <- sum(deg_all$Significant & deg_all$logFC < 0)
    cat("- Upregulated:", sig_up, "\n")
    cat("- Downregulated:", sig_down, "\n")
}

deg_stringent <- deg_all %>%
    filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
    arrange(deg_all$adj.P.val)
saveRDS(deg_stringent, file.path(DATA_DIR, "diff_exp_genes_stringent.rds"))

# ========================
# DIFFERENTIAL METHYLATION TFBS
# ========================
res <- readRDS(file.path(DATA_DIR, "ped88_dm_invasiveness_bh.rds"))

effect_threshold <- 0.2
fdr_threshold <- 0.05

hyper_probes <- res$Probe_ID[res$Est_Clinical_InvasivenessHigh > effect_threshold &
                                 res$BH_Clinical_InvasivenessHigh < fdr_threshold]
hypo_probes <- res$Probe_ID[res$Est_Clinical_InvasivenessHigh < -effect_threshold &
                                res$BH_Clinical_InvasivenessHigh < fdr_threshold]

hyper_enrichment <-
    testEnrichment(
        query = hyper_probes,
        databases = "TFBSconsensus",
        platform = "EPIC",
        universe = res$Probe_ID
    )
tf_hyper <- hyper_enrichment %>%
    filter(FDR < fdr_threshold) %>%
    pull(dbname) %>%
    unique()

hypo_enrichment <-
    testEnrichment(
        query = hypo_probes,
        databases = "TFBSconsensus",
        platform = "EPIC",
        universe = res$Probe_ID
    )
tf_hypo <- hypo_enrichment %>%
    filter(FDR < fdr_threshold) %>%
    pull(dbname) %>%
    unique()

saveRDS(tf_hyper, file.path(DATA_DIR, "dm_invasiveness_hyper_tfbs.rds"))
saveRDS(tf_hypo, file.path(DATA_DIR, "dm_invasiveness_hypo_tfbs.rds"))

# ========================
# 2E. TFBS VOLCANO PLOT
# ========================
# hyper_enrichment <- readRDS(file.path(DATA_DIR, "invasiveness_enrichment_hyper_tfbs.rds")
# hypo_enrichment <- readRDS(file.path(DATA_DIR, "invasiveness_enrichment_hypo_tfbs.rds"))

deg_all <- readRDS(file.path(DATA_DIR, "diff_exp_genes.rds"))
deg_all <- deg_all %>%
    mutate(
        significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                             "Significant", "Not Significant"),
        gene = rownames(.)
    )
deg_all <- deg_all %>%
    mutate(
        tf_status = case_when(
            gene %in% tf_hyper ~ "Hypermeth",
            gene %in% tf_hypo ~ "Hypometh",
            TRUE ~ "Other"
        )
    )

pdf(file.path(FIG_DIR, "diff_exp_meth_tf_volcano.pdf"), width=6.5, height=6, onefile=FALSE)
volcano_plot <- ggplot(deg_all, aes(x = logFC, y = -log10(adj.P.Val))) +
    # plot background dots
    geom_point(data = subset(deg_all, tf_status == "Other"),
               aes(color = tf_status), alpha = 0.25, size = 0.3) +
    # plot TF points on top
    geom_point(data = subset(deg_all, tf_status != "Other"),
               aes(color = tf_status), alpha = 1, size = 0.5) +
    # add lines as thresholds
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#757575ff") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#757575ff") +
    scale_color_manual(values = c("Hypermeth" = "#d61525ff", "Hypometh" = "#0c66bcff", "Other" = "grey")) +
    theme_minimal() +
    labs(x = "log2(FC)",
         y = "-log10(Adj. P-Val)",
         color = "TF Binding")
significant_tfs <- deg_all[deg_all$tf_status != "Other" &
                               deg_all$adj.P.Val < 0.05 &
                               abs(deg_all$logFC) > 1, ]
volcano_plot = volcano_plot +
    geom_text_repel(data = significant_tfs,
                    aes(label = gene),
                    max.overlaps = 6,
                    color = "#4d4d4dff")
plot(volcano_plot)
dev.off()

# ========================
# FIG 2F. TF HEATMAP
# ========================
mtx <- readRDS(file.path(DATA_DIR, "rna_seq_filtered.rds"))
tf_hyper <- readRDS(file.path(DATA_DIR, "dm_invasiveness_hyper_tfbs.rds"))
tf_hypo <- readRDS(file.path(DATA_DIR, "dm_invasiveness_hypo_tfbs.rds"))
deg_stringent <- readRDS(file.path(DATA_DIR, "diff_exp_genes_stringent.rds"))
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx")) %>%
    filter(Lymph_Node == "F") %>%
    filter(Sample_ID %in% colnames(mtx)) #84 x 20

invasiveness <- factor(ss$Clinical_Invasiveness, levels = c("Low", "High"))

deg_hyper_expressed <- deg_stringent %>%
    filter(rownames(.) %in% tf_hyper)
deg_hypo_expressed <- deg_stringent %>%
    filter(rownames(.) %in% tf_hypo)

cat("DE genes in hypermethylated TFs:", nrow(deg_hyper_expressed), "\n")
cat("DE genes in hypomethylated TFs:", nrow(deg_hypo_expressed), "\n")

deg_tf_expr <- mtx[unique(c(rownames(deg_hypo_expressed), rownames(deg_hyper_expressed))), ]

annotation_col <- data.frame(
    Cluster = as.factor(ss$Methylation_Clusters),
    Invasiveness = as.factor(ss$Clinical_Invasiveness),
    row.names = colnames(deg_tf_expr)
)

annotation_colors <- list(
    Cluster = cluster_colors,
    Invasiveness = invasiveness_colors
)

pdf(file.path(FIG_DIR, "diff_exp_meth_tf_heatmap.pdf"), width=8, height=5, onefile=FALSE)
pheatmap(deg_tf_expr,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         color = colorRampPalette(rev(brewer.pal(11, "PiYG")))(50),
         border_color = NA
)
dev.off()

# ========================
# FIG S2G. GSEA
# ========================
deg_all <- readRDS(file.path(DATA_DIR, "diff_exp_genes.rds"))

ranked_genes <- deg_all %>%
    mutate(
        ranking_metric = -log10(P.Value) * sign(logFC),
        gene_id = rownames(.)
    ) %>%
    arrange(desc(ranking_metric))

gene_ranks <- setNames(ranked_genes$ranking_metric, ranked_genes$gene_id)
h_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    mutate(collection = "H")

msig_list <- split(h_sets$gene_symbol, h_sets$gs_name)

gsea_results <- fgsea(
    pathways = msig_list,
    stats = gene_ranks,
    scoreType = 'std',
    minSize = 15,
    maxSize = 500,
    nPermSimple = 10000,
    eps = 0
) %>%
    as_tibble() %>%  #
    arrange(padj)
print(dim(gsea_results))
# 50  8

gsea_results_filtered <- gsea_results %>%
    filter(padj < 0.1) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = 20) %>%
    mutate(
        pathway = pathway %>%
            str_remove("^HALLMARK_") %>%
            str_replace_all("_", " ")
    )
# saveRDS(gsea_results_filtered, file.path(DATA_DIR, "gsea_results_filtered.rds"))
# gsea_results_filtered <- readRDS(file.path(DATA_DIR, "gsea_results_filtered.rds"))

pdf(file.path(FIG_DIR, "gsea_enrichment_dot_hallmark.pdf"), width=6, height=4, onefile=FALSE)
p <- ggplot(gsea_results_filtered,
            aes(x = NES,
                y = reorder(pathway, NES),
                size = -log10(padj),
                color = NES)) +
    geom_point() +
    scale_color_gradientn(colors = brewer.pal(9, "PuRd")[3:9]) +
    scale_size_continuous(name = "-log10(padj)") +
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        panel.grid.minor = element_blank()
    ) +
    labs(x = "Normalized Enrichment Score (NES)",
         y = "Pathway")
plot(p)
dev.off()
