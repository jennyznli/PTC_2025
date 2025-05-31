# ============================================================
#   FIG 4.
#   This script contains the analysis for final LOOCV driver
#   classifier development analysis, including:
#   4C. Confusion matrix
#   S4B. One-v-all ROC curves
# ============================================================

# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- here()
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

CL_DIR <- file.path(DATA_DIR, "loocv_driver_model")
MODEL_DIR <- file.path(CL_DIR, "models")
FEAT_DIR <- file.path(CL_DIR, "features")
SUM_DIR <- file.path(CL_DIR, "summary")
FIG_DIR <- file.path(BASE_DIR, "figures")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# Model Analysis
# ========================
overall_metrics <- readRDS(file.path(SUM_DIR, "loocv_overall_metrics.rds"))
class_metrics <- readRDS(file.path(SUM_DIR, "loocv_per_class_metrics.rds"))
pred_results <- read.csv(file.path(SUM_DIR, "loocv_all_predictions.csv"))
conf_matrix <- readRDS(file.path(SUM_DIR, "loocv_conf_matrix.rds"))
model_list <- readRDS(file.path(SUM_DIR, "loocv_model_list.rds"))
roc_obj <- readRDS(file.path(SUM_DIR, "loocv_roc_curves.rds"))
imp_results <- readRDS(file.path(SUM_DIR, "loocv_imp_results.rds"))

print(overall_metrics)
# Accuracy Overall_Kappa AUC_BRAF V600E AUC_DICER1 AUC_Kinase Fusion AUC_Ras-like
# Kappa 0.9411765     0.9159082      0.9931457          1         0.9702857    0.9761538
# Mean_AUC Macro_Sensitivity Macro_Specificity Macro_Precision Macro_Recall  Macro_F1
#  0.9848963          0.950487         0.9783394            0.95     0.950487 0.9499797
# Macro_BalancedAccuracy
# Kappa              0.9644132

print(class_metrics)
# Class Sensitivity Specificity Precision    Recall              F1       Balanced         Accuracy
# Sensitivity     BRAF V600E   0.9545455   0.9841270 0.9545455 0.9545455 0.9545455        0.9693362
# Sensitivity1        DICER1   1.0000000   0.9870130 0.8888889 1.0000000 0.9411765        0.9935065
# Sensitivity2 Kinase Fusion   0.9428571   0.9600000 0.9428571 0.9428571 0.9428571        0.9514286
# Sensitivity3      Ras-like   0.8500000   0.9692308 0.8947368 0.8500000 0.8717949        0.9096154
# 1            Macro-Average   0.9368506   0.9750927 0.9202571 0.9368506 0.9275935        0.9559717

print(conf_matrix)
# Reference
# Prediction      BRAF V600E DICER1 Kinase Fusion Ras-like
# BRAF V600E            21      0             1        0
# DICER1                 0      8             0        1
# Kinase Fusion          0      0            33        2
# Ras-like               1      0             1       17
#
# Overall Statistics
#
# Accuracy : 0.9294
# 95% CI : (0.8527, 0.9737)
# No Information Rate : 0.4118
# P-Value [Acc > NIR] : < 2.2e-16
# Statistics by Class:
#
#     Class: BRAF V600E Class: DICER1 Class: Kinase Fusion Class: Ras-like
# Sensitivity                     0.9545       1.00000               0.9429          0.8500
# Specificity                     0.9841       0.98701               0.9600          0.9692
# Pos Pred Value                  0.9545       0.88889               0.9429          0.8947
# Neg Pred Value                  0.9841       1.00000               0.9600          0.9545
# Prevalence                      0.2588       0.09412               0.4118          0.2353
# Detection Rate                  0.2471       0.09412               0.3882          0.2000
# Detection Prevalence            0.2588       0.10588               0.4118          0.2235
# Balanced Accuracy               0.9693       0.99351               0.9514          0.9096

# ========================
# S4B. ROC CURVES
# ========================
roc_data <- data.frame()
class_names <- names(roc_obj)
n_curves <- length(roc_obj)

for (i in 1:length(roc_obj)) {
    current_roc <- roc_obj[[i]]
    current_auc <- as.numeric(current_roc$auc)
    print(current_auc)
    temp_data <- data.frame(
        FPR = 1 - current_roc$specificities,
        TPR = current_roc$sensitivities,
        Class = class_names[i],
        AUC = current_auc
    )
    print(dim(temp_data))
    roc_data <- rbind(roc_data, temp_data)
}
dim(roc_data)

p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Class)) +
    geom_line(size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = driver_colors,
                       labels = paste0(class_names, " (AUC = ",
                                       round(tapply(roc_data$AUC, roc_data$Class, unique), 3), ")")) +
    labs(x = "False Positive Rate",
         y = "True Positive Rate") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_blank()
    ) +
    coord_equal()

pdf(file.path(FIG_DIR, "loocv_driver_roc.pdf"),
    height = 5,
    width = 5)
print(p)
dev.off()
