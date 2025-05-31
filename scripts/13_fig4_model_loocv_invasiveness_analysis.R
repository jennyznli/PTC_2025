# ============================================================
#   FIG 4.
#   This script contains the analysis for Fig. 4, including:
#   Final LOOCV driver classifier development analysis
#   4D. Confusion matrix
#   S4A. ROC curve
# ============================================================

# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- here()
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

CL_DIR <- file.path(DATA_DIR, "loocv_invasiveness_model")
MODEL_DIR <- file.path(CL_DIR, "models")
FEAT_DIR <- file.path(CL_DIR, "features")
SUM_DIR <- file.path(CL_DIR, "summary")
FIG_DIR <- file.path(here(), "figures")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

print("Configuration complete.")

# ========================
# 4B. CONFUSION MATRIX
# ========================
pred_results <- read.csv(file.path(SUM_DIR, "loocv_all_predictions.csv"))
conf_matrix <- readRDS(file.path(SUM_DIR, "loocv_conf_matrix.rds"))
overall_metrics <- readRDS(file.path(SUM_DIR, "loocv_overall_metrics.rds"))
model_list <- readRDS(file.path(SUM_DIR, "loocv_model_list.rds"))
roc_obj <- readRDS(file.path(SUM_DIR, "loocv_roc_curve.rds"))
imp_results <- readRDS(file.path(SUM_DIR, "loocv_imp_results.rds"))

print(conf_matrix)
# Reference
# Prediction High Low
# High   55   3
# Low     2  28
# Accuracy : 0.9432
# 95% CI : (0.8724, 0.9813)
# No Information Rate : 0.6477
# P-Value [Acc > NIR] : 5.272e-11

print(overall_metrics)
# Accuracy Accuracy_Threshold Balanced_Accuracy     Kappa
# Balanced Accuracy 0.9431818          0.9431818          0.934069 0.8745724
# AUC Precision    Recall Specificity        F1
# Balanced Accuracy 0.9932088 0.9482759 0.9649123   0.9032258 0.9565217

# ========================
# S4A. ROC CURVE
# ========================
roc_data <- data.frame(
    Specificity = 1 - roc_obj$specificities,
    Sensitivity = roc_obj$sensitivities
)

p <- ggplot(roc_data, aes(x = Specificity, y = Sensitivity)) +
    geom_line(color = "#f8766dff", size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    annotate("text", x = 0.75, y = 0.25,
             label = paste("AUC =", round(as.numeric(overall_metrics$AUC), 3)),
             size = 5, color = "#f8766dff") +
    labs(x = "False Positive Rate",
         y = "True Positive Rate") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
    ) +
    coord_equal()

pdf(file.path(FIG_DIR, "loocv_invasiveness_roc.pdf"), height = 4, width = 4)
print(p)
dev.off()

