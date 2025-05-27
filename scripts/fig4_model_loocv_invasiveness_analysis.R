# ============================================================
#   This script contains the analysis for Fig. 4, including:
#   Final LOOCV driver classifier development analysis
#
#   S4A ROC curve
#   S4D SHAP analysis
# ============================================================

source(here::here("functions", "load_packages.R"))
source(here::here("functions", "color_keys.R"))
source(here::here("functions", "functions.R"))

# ========================
# Model Analysis
# ========================

date <- "20250316"
dir <- file.path(here(), "data", "20250316_loocv_invasiveness")
sum_dir <- file.path(dir)
fig_dir <- file.path(here(), "figures")

model_list <- readRDS(file.path(sum_dir, paste0(date, "_loocv_model_list.rds"))) #264000
pred_results <- read.csv(file.path(sum_dir, paste0(date, "_loocv_all_predictions.csv")))
imp_results <- readRDS(file.path(sum_dir, paste0(date, "_loocv_imp_results.rds"))) #264000
feature_list <- readRDS(file.path(sum_dir, paste0(date, "_loocv_feature_list.rds"))) #264000
class_metrics <- readRDS(file.path(sum_dir, paste0(date, "_loocv_per_class_metrics.rds"))) #264000
conf_matrix <- readRDS(file.path(sum_dir, paste0(date, "_loocv_conf_matrix.rds"))) #264000
overall_metrics <- readRDS(file.path(sum_dir, paste0(date, "_loocv_overall_metrics.rds"))) #264000

# Driver
# Accuracy     Kappa AUC_BRAF V600E AUC_Indeterminate AUC_Kinase Fusion AUC_Ras-like
# Kappa 0.9090909 0.8632213      0.9938017         0.6509804         0.9638814    0.9696429

# ========================
# Confusion Matrix
# ========================

ss <- readxl::read_excel(file.path(here(), "ss", "20231102_thyroid_master.xlsx")) |>
    dplyr::filter(Include_In_Analysis == "1")
labels <- as.factor(ss$Driver_Group_2)

accuracy <- mean(pred_results$Correct)
# 0.9090909

# X BRAF.V600E Indeterminate Kinase.Fusion Ras.like
# 1    BRAF V600E         20             0             1        0
# 2 Indeterminate          0             0             0        0
# 3 Kinase Fusion          1             2            33        1
# 4      Ras-like          1             1             1       27

# ========================
# ROC Curve Analysis
# ========================
roc_obj <- readRDS(file.path(sum_dir, paste0(date, "_loocv_roc_curve.rds")))
imp_results <- readRDS(file.path(sum_dir, paste0(date, "_loocv_imp_results.rds"))) # 264000      6

roc_data <- data.frame(
    Specificity = 1 - roc_obj$specificities,
    Sensitivity = roc_obj$sensitivities
)

p <- ggplot(roc_data, aes(x = Specificity, y = Sensitivity)) +
    geom_line(color = "blue", size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    annotate("text", x = 0.75, y = 0.25,
             label = paste("AUC =", round(as.numeric(overall_metrics$AUC), 3)),
             size = 5, color = "darkblue") +
    labs(x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
    ) +
    coord_equal()

pdf(file.path(fig_dir, "20250226_thyroid88_loocv_driver_roc.pdf"), height = 4, width = 4)
print(p)
dev.off()

# Load required packages
library(ggplot2)
library(pROC)
library(dplyr)
library(tidyr)

# Load your prediction results
pred_results <- readRDS(file.path(sum_dir, paste0(date, "_loocv_pred_results.rds")))

# Get unique classes
classes <- unique(pred_results$True_Label)
num_classes <- length(classes)

# Initialize lists to store ROC objects and data
roc_objects <- list()
roc_curve_data <- data.frame()
auc_values <- numeric(num_classes)
names(auc_values) <- classes

# Create a color palette for multiple classes
colors <- rainbow(num_classes)

# Calculate ROC curves for each class (one-vs-rest)
for (i in 1:num_classes) {
    current_class <- classes[i]

    # Create binary labels (1 for current class, 0 for others)
    binary_true <- ifelse(pred_results$True_Label == current_class, 1, 0)

    # Get probabilities for the current class
    # Assuming your prediction results have probability columns for each class
    # If you only have one probability column, you may need to adjust this
    class_probs <- pred_results[[paste0("Probability_", current_class)]]

    # If you only have binary data with one probability column:
    if (!exists("class_probs") || is.null(class_probs)) {
        if (current_class == classes[2]) {  # Assuming this is the positive class in binary case
            class_probs <- pred_results$Probability
        } else {
            class_probs <- 1 - pred_results$Probability
        }
    }

    # Calculate ROC
    roc_objects[[i]] <- roc(binary_true, class_probs, quiet = TRUE)
    auc_values[i] <- as.numeric(auc(roc_objects[[i]]))

    # Extract coordinates for plotting
    roc_data <- data.frame(
        Specificity = 1 - roc_objects[[i]]$specificities,
        Sensitivity = roc_objects[[i]]$sensitivities,
        Class = current_class
    )

    roc_curve_data <- rbind(roc_curve_data, roc_data)
}

# Calculate macro-average AUC
macro_auc <- mean(auc_values)

# Create the plot
p <- ggplot(roc_curve_data, aes(x = Specificity, y = Sensitivity, color = Class)) +
    geom_line(size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = colors) +
    labs(
        title = "ROC Curves (One-vs-Rest)",
        x = "False Positive Rate (1 - Specificity)",
        y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right"
    ) +
    coord_equal() +
    # Add AUC annotations for each class
    annotate(
        "text",
        x = rep(0.75, num_classes),
        y = seq(0.25, 0.15, length.out = num_classes),
        label = paste(classes, "AUC =", round(auc_values, 3)),
        color = colors,
        size = 3.5,
        hjust = 0
    ) +
    # Add macro-average AUC
    annotate(
        "text",
        x = 0.75,
        y = 0.05,
        label = paste("Macro-avg AUC =", round(macro_auc, 3)),
        size = 4,
        fontface = "bold",
        color = "black"
    )

# Save the plot
pdf(file.path(fig_dir, "20250226_thyroid88_loocv_multiclass_roc.pdf"), height = 6, width = 6)
print(p)
dev.off()

# For comparison, you might also want to save the original binary ROC curve if applicable
if (num_classes == 2) {
    binary_roc <- roc_objects[[2]]  # Assuming the second class is the positive class
    binary_data <- data.frame(
        Specificity = 1 - binary_roc$specificities,
        Sensitivity = binary_roc$sensitivities
    )

    p_binary <- ggplot(binary_data, aes(x = Specificity, y = Sensitivity)) +
        geom_line(color = "blue", size = 1) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        annotate("text", x = 0.75, y = 0.25,
                 label = paste("AUC =", round(auc_values[2], 3)),
                 size = 5, color = "darkblue") +
        labs(title = "Binary ROC Curve",
             x = "False Positive Rate (1 - Specificity)",
             y = "True Positive Rate (Sensitivity)") +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10)
        ) +
        coord_equal()

    pdf(file.path(fig_dir, "20250226_thyroid88_loocv_binary_roc.pdf"), height = 4, width = 4)
    print(p_binary)
    dev.off()
}

