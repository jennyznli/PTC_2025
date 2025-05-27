# ============================================================
#   FIG 4.
#   This script contains the analysis for final LOOCV driver
#   classifier development analysis, including:
#
#   S4B. One-v-all ROC curves
#   S4E-H. One-v-all SHAP analysis
# ============================================================

# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20
DATE <- format(Sys.DATE(), "%Y%m%d")

BASE_DIR <- "/home/lijz/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

CL_DIR <- file.path(BASE_DIR, "loocv_driver_model")
MODEL_DIR <- file.path(CL_DIR, "models")
FEAT_DIR <- file.path(CL_DIR, "features")
SUM_DIR <- file.path(CL_DIR, "summary")
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FEAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUM_DIR, recursive = TRUE, showWarnings = FALSE)

models_list <- list()
pred_results <- data.frame()
imp_results <- data.frame()
feature_list <- c()

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

print("Configuration complete.")

# ========================
# Model Analysis
# ========================

date <- "20250316"
dir <- file.path(here(), "data", "20250316_loocv_driver")
sum_dir <- file.path(dir)
fig_dir <- file.path(here(), "figures")

model_list <- readRDS(file.path(sum_dir, paste0(date, "_loocv_model_list.rds"))) #264000
pred_results <- read.csv(file.path(sum_dir, paste0(date, "_loocv_all_predictions.csv")))
imp_results <- readRDS(file.path(sum_dir, paste0(date, "_loocv_imp_results.rds"))) #264000
feature_list <- readRDS(file.path(sum_dir, paste0(date, "_loocv_feature_list.rds"))) #264000
class_metrics <- readRDS(file.path(sum_dir, paste0(date, "_loocv_per_class_metrics.rds"))) #264000
conf_matrix <- readRDS(file.path(sum_dir, paste0(date, "_loocv_conf_matrix.rds"))) #264000
overall_metrics <- readRDS(file.path(sum_dir, paste0(date, "_loocv_overall_metrics.rds"))) #264000

# Accuracy Overall_Kappa AUC_BRAF V600E AUC_DICER1 AUC_Kinase Fusion AUC_Ras-like
# Kappa 0.9411765     0.9159082      0.9931457          1         0.9702857    0.9761538
# Mean_AUC Macro_Sensitivity Macro_Specificity Macro_Precision Macro_Recall  Macro_F1
#  0.9848963          0.950487         0.9783394            0.95     0.950487 0.9499797
# Macro_BalancedAccuracy
# Kappa              0.9644132


print(conf_matrix)
# Accuracy : 0.9412
#
# Reference
# Prediction      BRAF V600E DICER1 Kinase Fusion Ras-like
# BRAF V600E            20      0             1        0
# DICER1                 0      8             0        0
# Kinase Fusion          1      0            33        1
# Ras-like               1      0             1       19

# Class:                        BRAF V600E      DICER1          Kinase Fusion        Ras-like
# Sensitivity                     0.9091       1.00000               0.9429          0.9500
# Specificity                     0.9841       1.00000               0.9600          0.9692
# Pos Pred Value                  0.9524       1.00000               0.9429          0.9048
# Neg Pred Value                  0.9688       1.00000               0.9600          0.9844
# Prevalence                      0.2588       0.09412               0.4118          0.2353
# Detection Rate                  0.2353       0.09412               0.3882          0.2235
# Detection Prevalence            0.2471       0.09412               0.4118          0.2471
# Balanced Accuracy               0.9466       1.00000               0.9514          0.9596

#                       Class Sensitivity Specificity Precision    Recall        F1
# Sensitivity     BRAF V600E   0.9090909   0.9841270 0.9523810 0.9090909 0.9302326
# Sensitivity1        DICER1   1.0000000   1.0000000 1.0000000 1.0000000 1.0000000
# Sensitivity2 Kinase Fusion   0.9428571   0.9600000 0.9428571 0.9428571 0.9428571
# Sensitivity3      Ras-like   0.9500000   0.9692308 0.9047619 0.9500000 0.9268293
# 1            Macro-Average   0.9504870   0.9783394 0.9500000 0.9504870 0.9499797
#               BalancedAccuracy
# Sensitivity         0.9466089
# Sensitivity1        1.0000000
# Sensitivity2        0.9514286
# Sensitivity3        0.9596154
# 1                   0.9644132

# ========================
# ROC Curve Analysis
# ========================
roc_obj <- readRDS(file.path(sum_dir, paste0(date, "_loocv_roc_curves.rds")))
roc_data <- data.frame()
class_names <- names(roc_obj)
n_curves <- length(roc_obj)
driver_colors <- c(
    "BRAF V600E" = "#ff6cc3ff",
    "Kinase Fusion" = "#20bb20ff",
    "Ras-like" = "#00b4f0ff",
    "DICER1" = "#b47cffff"
)

for (i in 1:length(roc_obj)) {
    current_roc <- roc_obj[[i]]
    current_auc <- as.numeric(current_roc$auc)
    print(current_auc)
    temp_data <- data.frame(
        Specificity = 1 - current_roc$specificities,
        Sensitivity = current_roc$sensitivities,
        Class = class_names[i],
        AUC = current_auc
    )
    print(dim(temp_data))
    roc_data <- rbind(roc_data, temp_data)
}
dim(roc_data) #239   4, since it's one v all

p <- ggplot(roc_data, aes(x = Specificity, y = Sensitivity, color = Class)) +
    geom_line(size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = driver_colors,
                       labels = paste0(class_names, " (AUC = ",
                                       round(tapply(roc_data$AUC, roc_data$Class, unique), 3), ")")) +
    labs(x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_blank()
    ) +
    coord_equal()

pdf(file.path(fig_dir, paste0(date, "_thyroid88_loocv_drivers_roc.pdf")),
    height = 6,
    width = 6)
print(p)
dev.off()


# ========================
# SHAP Analysis - CHECK THIS
# ========================
set.seed(123)
train_data <- as.data.frame(t(train_betas3k))
train_labels <- as.factor(ss$Driver_Group)

# model <- readRDS(file.path(dir, "20250316_final_subset_model_3k.rds"))
X_train <- as.data.frame(t(train_betas3k))

unique_classes <- levels(train_labels)

set.seed(123)
for (class_name in unique_classes) {
    binary_target <- as.numeric(labels == class_name)

    set.seed(123)
    binary_model <- randomForest(
        x = X_train,
        y = as.factor(binary_target),
        mtry = 2,
        nodesize = 5,
        importance = TRUE
    )
    saveRDS(binary_model, file.path(dir, paste0(date, "_", class_name, "_model.rds")))

    set.seed(123)
    unified_model <- randomForest.unify(binary_model, X_train)
    shap <- treeshap(unified_model, train_data, interactions = TRUE)
    shp <- shapviz::shapviz(shap, X = train_data)

    saveRDS(shap, file.path(dir, paste0(date, "_", class_name, "_shap.rds")))
    saveRDS(shp, file.path(dir, paste0(date, "_", class_name, "_shpviz.rds")))

    # most important features
    most_imp <- sort(colMeans(abs(shap$shaps)), decreasing = TRUE)
    print(head(most_imp))

    saveRDS(most_imp, file.path(dir, paste0(date, "_", class_name, "_shap_importance.rds")))
    write.csv(data.frame(Feature = names(most_imp),
                         SHAP_Value = most_imp),
              file.path(dir, paste0(date, "_", class_name, "_shap_importance.csv")))
}


## PLOTTING=========================================================================
set.seed(123)
fig_dir <- file.path(here(), "figures")
date <- "20250404"
dir <- file.path(here(), "data", "20250316_fmodel_drivers")

# model <- readRDS(file.path(dir, "20250316_Kinase Fusion_model.rds"))
# shp <- readRDS(file.path(dir, "20250316_BRAF V600E_shpviz.rds"))
# shap <- readRDS(file.path(dir, "20250316_BRAF V600E_shap.rds"))


# beeswarm
p <- sv_importance.shapviz(shp, kind = "beeswarm", max_display = 30)
pdf(file.path(fig_dir, paste0(date, "_", "BRAF V600E", "_thyroid85_beeswarm.pdf")), height = 6, width = 5)
plot(p)
dev.off()

# contribution
p <- plot_contribution(shap, obs = 1, max_vars = 30)
pdf(file.path(fig_dir, paste0(date, "_", "BRAF V600E", "_thyroid85_contribution.pdf")), height = 6, width = 5)
plot(p)
dev.off()

for (class_name in unique_classes) {
    shp <- readRDS(file.path(dir, paste0("20250316_", class_name, "_shpviz.rds")))
    shap <- readRDS(file.path(dir, paste0("20250316_", class_name, "_shap.rds")))

    # beeswarm
    p <- sv_importance.shapviz(shp, kind = "beeswarm", max_display = 30)
    pdf(file.path(fig_dir, paste0(date, "_", class_name, "_thyroid85_beeswarm.pdf")), height = 6, width = 5)
    plot(p)
    dev.off()

    # contribution
    p <- plot_contribution(shap, obs = 1, max_vars = 30)
    pdf(file.path(fig_dir, paste0(date, "_", class_name, "_thyroid85_contribution.pdf")), height = 7, width = 10)
    plot(p)
    dev.off()
}

# ========================
# Enrichment Analysis
# ========================
res <- testEnrichment(sel_probes, platform = "EPIC")

# Low v High TFBS HYPO
pdf(here("figures", "20250226_thyroid88_rfmodel_feat_enrichment_all.pdf"), width=5, height=7, onefile=FALSE)
KYCG_plotEnrichAll(res)
dev.off()

plotDot2 <- function(df, n_min = 20, n_max = 20, max_fdr = 0.05) {
    db1 <- FDR <- overlap <- estimate <- NULL
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))

    df1 <- preparePlotDF(df, n_min, n_max, max_fdr)
    ggplot(df1) +
        geom_point(aes(dbname, estimate, size=overlap, color=-log10(FDR))) +
        coord_flip() +
        scale_color_gradientn(colors = parula(20)) +
        ylab("Estimate (OR)") + xlab("") +
        theme_minimal()  # Adjust the base font size if needed
}
x <- testEnrichment(sel_probes, "TFBSconsensus", platform = "EPIC")
pdf(here("figures", "20250226_thyroid88_rfmodel_feat_enrichment_tfbs.pdf"), width=3.5, height=4, onefile=FALSE)
plotDot2(x, n_min = 40)
dev.off()



