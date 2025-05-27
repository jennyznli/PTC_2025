# ============================================================
#   FIG 4.
#   This script performs leave-one-out cross validation (LOOCV)
#   to develop and validate a random forest classifier for
#   predicting clinical invasiveness in thyroid cancer samples
#   using DNA methylation data.
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

CL_DIR <- file.path(BASE_DIR, "loocv_invasiveness_model")
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
# Load Packages
# ========================
required_packages <- c(
    "randomForest", "pROC", "caret", "mltools",
    "readxl", "dplyr", "BiocParallel", "batchtools", "here")
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
        stop("Package ", pkg, " not available")
    }
}

# ========================
# Core Functions
# ========================
subset_feature_selection <- function(betas, batch_size, n_feat, train_idx, labels, n_trees = 500, fold = NULL) {
    n_probes <- nrow(betas)
    n_batches <- ceiling(n_probes/batch_size)
    fold_models <- list()
    fold_importance <- data.frame()

    for(s in 0:(n_batches-1)) {
        # Calculate indices for current batch
        start_idx <- s * batch_size + 1
        end_idx <- min((s + 1) * batch_size, n_probes)
        sbetas <- betas[start_idx:end_idx, , drop=FALSE]
        print(paste("Batch dimensions:", dim(sbetas)[1], "x", dim(sbetas)[2]))

        model <- randomForest(
            x = t(sbetas[, train_idx]),
            y = as.factor(labels[train_idx]),
            ntree = n_trees,
            importance = TRUE
        )
        fold_models[[s+1]] <- model

        importance <- data.frame(
            Fold = if(is.null(fold)) NA else fold,
            Model = s,
            Feature = rownames(model$importance),
            MeanDecreaseAccuracy = model$importance[,"MeanDecreaseAccuracy"]
        )
        fold_importance <- rbind(fold_importance, importance)
    }

    features <- fold_importance %>%
        arrange(desc(MeanDecreaseAccuracy)) %>%
        head(n_feat) %>%
        pull(Feature)

    return(features)
}

process_sample <- function(i, all_samples, betas, labels) {
    print(paste("Processing sample", i, "of", length(all_samples)))

    # Create LOOCV fold: one test sample, all others for training
    valid_idx <- all_samples[i]
    train_idx <- all_samples[-i]

    # Feature selection
    sel_probes <- subset_feature_selection(
        betas = betas,
        batch_size = 10000,
        n_feat = 3000,
        train_idx = train_idx,
        labels = labels,
        fold = i
    )

    # Model training
    model <- randomForest(
        x = t(betas[sel_probes, train_idx]),
        y = labels[train_idx],
        ntree = 500,
        importance = TRUE
    )

    # IMPORTANCE
    importance <- data.frame(
        Feature = rownames(model$importance),
        MeanDecreaseAccuracy = model$importance[,"MeanDecreaseAccuracy"],
        MeanDecreaseGini = model$importance[,"MeanDecreaseGini"]
    )

    # Add class-specific importance
    for (level in levels(labels)) {
        if(level %in% colnames(model$importance)) {
            importance[,level] <- model$importance[,level]
        }
    }

    importance$Fold <- i

    # PREDICTIONS
    pred_probs <- predict(model, t(betas[sel_probes, valid_idx, drop=FALSE]), type = "prob")
    pred_class <- predict(model, t(betas[sel_probes, valid_idx, drop=FALSE]))

    preds <- data.frame(
        Sample_ID = valid_idx,
        True_Label = labels[valid_idx],
        Predicted = pred_class,
        Correct = (pred_class == labels[valid_idx]),
        Fold = i
    )

    # Add probability columns for each class
    for (j in 1:length(levels(labels))) {
        level_name <- levels(labels)[j]
        # Create a safe column name that matches R's automatic name conversion
        safe_level_name <- gsub(" ", ".", gsub("-", ".", level_name))
        col_name <- paste0("Prob_", safe_level_name)
        preds[,col_name] <- pred_probs[,j]

        # Store the mapping between original level and column name
        if (j == 1) {
            preds$LevelToColMapping <- paste(level_name, "->", col_name)
        }
    }
    preds$LevelToColMapping <- NULL

    # ========================
    # SAVING RESULTS
    # ========================
    saveRDS(model, file.path(MODEL_DIR, paste0(DATE, "_model_sample", i, ".rds")))
    saveRDS(importance, file.path(FEAT_DIR, paste0(DATE, "_importance_sample", i, ".rds")))

    return(list(
        sample = i,
        prediction = pred_class,
        prediction_probs = pred_probs,
        true_label = labels[valid_idx],
        preds = preds,
        importance = importance,
        selected_features = sel_probes,
        model = model,
        status = "complete"
    ))
}

# ========================
# Main Execution
# ========================
# Load data
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx")) %>%
    dplyr::filter(Lymph_Node == "F")

betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))
valid_ids <- ss$IDAT[ss$IDAT %in% colnames(betas)]
betas <- betas[, valid_ids]
ss <- ss[match(valid_ids, ss$IDAT), ]

# Create sample indices for LOOCV
all_samples <- 1:length(labels)
print(paste("Total samples for LOOCV:", length(all_samples)))

# Run LOOCV
print("Starting LOOCV")
results <- bplapply(1:length(all_samples), function(i) process_sample(i, all_samples, betas, labels), BPMulticoreParam(N_WORKERS) = MulticoreParam(N_WORKERS))
print("Completed all samples")

# Process and combine results
for (i in 1:length(results)) {
    # Store the model
    models_list[[i]] <- results[[i]]$model

    # Accumulate predictions
    pred_results <- rbind(pred_results, results[[i]]$preds)

    # Accumulate importance scores
    imp_results <- rbind(imp_results, results[[i]]$importance)

    # Accumulate selected features
    feature_list <- c(feature_list, results[[i]]$selected_features)
}

# Calculate confusion matrix metrics
conf_matrix <- caret::confusionMatrix(
    data = as.factor(pred_results$Predicted),
    reference = as.factor(pred_results$True_Label)
)

# Calculate accuracy from predictions
accuracy <- mean(pred_results$Correct)

# Initialize overall metrics dataframe
overall_metrics <- data.frame(
    Accuracy = accuracy,
    Overall_Kappa = conf_matrix$overall["Kappa"]
)

# Get the name of the probability column for the second class (usually the positive class)
positive_class <- levels(labels)[2]
safe_positive_class <- gsub(" ", ".", gsub("-", ".", positive_class))
prob_col_name <- paste0("Prob_", safe_positive_class)

# Create ROC object
roc_obj <- pROC::roc(response = pred_results$True_Label,
                     predictor = pred_results[[prob_col_name]],
                     levels = levels(pred_results$True_Label))

# Get AUC value
auc_value <- pROC::auc(roc_obj)

# Calculate accuracy with threshold
all_class_preds <- ifelse(pred_results[[prob_col_name]] > 0.5, levels(labels)[2], levels(labels)[1])
all_class_preds <- factor(all_class_preds, levels = levels(labels))
accuracy_threshold <- mean(all_class_preds == pred_results$True_Label)

overall_metrics <- data.frame(
    Accuracy = accuracy,
    Accuracy_Threshold = accuracy_threshold,
    Balanced_Accuracy = conf_matrix$byClass["Balanced Accuracy"],
    Kappa = conf_matrix$overall["Kappa"],
    AUC = as.numeric(auc_value),
    Precision = conf_matrix$byClass["Pos Pred Value"],
    Recall = conf_matrix$byClass["Sensitivity"],
    Specificity = conf_matrix$byClass["Specificity"],
    F1 = conf_matrix$byClass["F1"]
)

# ========================
# Save Final Results
# ========================
write.csv(pred_results, file.path(SUM_DIR, paste0(DATE, "_loocv_all_predictions.csv")))
saveRDS(conf_matrix, file.path(SUM_DIR, paste0(DATE, "_loocv_conf_matrix.rds")))
saveRDS(feature_list, file.path(SUM_DIR, paste0(DATE, "_loocv_feature_list.rds")))
saveRDS(overall_metrics, file.path(SUM_DIR, paste0(DATE, "_loocv_overall_metrics.rds")))
saveRDS(models_list, file.path(SUM_DIR, paste0(DATE, "_loocv_model_list.rds")))
saveRDS(roc_obj, file.path(SUM_DIR, paste0(DATE, "_loocv_roc_curve.rds")))
saveRDS(imp_results, file.path(SUM_DIR, paste0(DATE, "_loocv_imp_results.rds")))

# Print final metrics
print("LOOCV Complete")
print(overall_metrics)
print(conf_matrix)

