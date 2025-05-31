# ============================================================
#   FIG 4.
#   This script trains the final random forest classifier for
#   predicting driver group in thyroid cancer samples
#   using all 85 primary samples (excluding 3 indeterminate).
#   S4E-H. SHAP beeswarm plot for one vs. all final models.
#   4H. Enrichment plot of 3k most important probes.
# ============================================================

# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- here()
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

MODEL_DIR <- file.path(DATA_DIR, "final_driver_model")
SHAP_DIR <- file.path(MODEL_DIR, "shap")
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SHAP_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# PARAMETERS
# ========================
BATCH_SIZE = 10000
N_TREES <- 500

# ========================
# TRAIN/TEST SPLIT
# (test is indeterminate)
# ========================
ss_all <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx")) %>%
    dplyr::filter(Lymph_Node == "F")

betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))

valid_ids <- intersect(ss_all$IDAT, colnames(betas))
ss_all <- ss_all[match(valid_ids, ss_all$IDAT), ]
betas <- betas[, valid_ids]

ss_train <- ss_all %>% filter(Driver_Group != "Indeterminate")
ss_test  <- ss_all %>% filter(Driver_Group == "Indeterminate")

betas_train <- betas[, ss_train$IDAT]
betas_test  <- betas[, ss_test$IDAT]

train_labels <- as.factor(ss_train$Driver_Group)
test_labels <- as.factor(ss_test$Driver_Group)

# ========================
# FEATURE SELECTION
# ========================
n_probes <- nrow(betas)
n_batches <- ceiling(n_probes/BATCH_SIZE)
fold_models <- list()
fold_importance <- data.frame()

for(s in 0:(n_batches-1)) {
    # Calculate indices for current batch
    start_idx <- s * BATCH_SIZE + 1
    end_idx <- min((s + 1) * BATCH_SIZE, n_probes)
    sbetas <- betas_train[start_idx:end_idx, ]

    model <- randomForest(
        x = t(sbetas),
        y = as.factor(train_labels),
        ntree = N_TREES,
        importance = TRUE
    )

    fold_models[[s+1]] <- model
    importance <- data.frame(
        Model = s,
        Feature = rownames(model$importance),
        MeanDecreaseAccuracy = model$importance[,"MeanDecreaseAccuracy"],
        MeanDecreaseGini = model$importance[,"MeanDecreaseGini"]
    )
    fold_importance <- rbind(fold_importance, importance)
}

features <- fold_importance %>%
    arrange(desc(MeanDecreaseAccuracy))

saveRDS(fold_models, file.path(MODEL_DIR, "fold_models.rds"))
saveRDS(features, file.path(MODEL_DIR, "fold_importance.rds"))

# ========================
# MODEL DEVELOPMENT
# ========================
imp <- readRDS(file.path(MODEL_DIR, "fold_importance.rds"))
sel_probes3k <- rownames(imp %>% head(3000))
train_betas3k <- betas_train[sel_probes3k, ]

# TRAIN MODEL
model <- randomForest(
    x = t(train_betas3k),
    y = train_labels,
    mtry = 2,
    nodesize = 5,
    importance = TRUE
)

# OOB estimate of  error rate: 4.71%
# Confusion matrix:
#     BRAF V600E DICER1 Kinase Fusion Ras-like class.error
# BRAF V600E            21      0             0        1  0.04545455
# DICER1                 0      6             0        2  0.25000000
# Kinase Fusion          0      0            34        1  0.02857143
# Ras-like               0      0             0       20  0.00000000

saveRDS(model, file.path(MODEL_DIR, "final_driver_model.rds"))

# ========================
# PREDICT ON 3 INDETERMINATE CASES
# ========================
betas_test3k  <- betas_test[sel_probes3k, ]

pred <- predict(model, t(betas_test3k))
prob <- predict(model, t(betas_test3k), type = "prob")

pred_df <- data.frame(
    Sample_ID = colnames(betas_test3k),
    Predicted_Class = pred
)
prob_df <- as.data.frame(prob)
prob_df$Sample_ID <- rownames(prob_df)

results_df <- merge(pred_df, prob_df, by = "Sample_ID")

write.csv(results_df, file.path(MODEL_DIR, "indeterminate_predictions.csv"), row.names = FALSE)

print(results_df)
# Sample_ID Predicted_Class BRAF V600E DICER1 Kinase Fusion Ras-like
# 1 207686140028_R03C01   Kinase Fusion      0.154  0.030         0.736    0.080
# 2 207686140033_R02C01   Kinase Fusion      0.152  0.020         0.676    0.152
# 3 207700160043_R03C01        Ras-like      0.006  0.292         0.044    0.658

# ========================
# ONE V ALL BINARY
# ========================
unique_classes <- levels(train_labels)

for (class_name in unique_classes) {
    binary_target <- as.numeric(train_labels == class_name)
    train_data <- t(train_betas3k)

    binary_model <- randomForest(
        x = train_data,
        y = as.factor(binary_target),
        mtry = 2,
        nodesize = 5,
        importance = TRUE
    )
    clean <- clean_label(class_name)
    saveRDS(binary_model, file.path(MODEL_DIR, paste0(clean, "_shap_model.rds")))

    unified_model <- randomForest.unify(binary_model, train_data)
    shap <- treeshap(unified_model, train_data, interactions = TRUE)
    shp <- shapviz::shapviz(shap, X = train_data)

    saveRDS(shap, file.path(SHAP_DIR, paste0(clean, "_shap.rds")))
    saveRDS(shp, file.path(SHAP_DIR, paste0(clean, "_shpviz.rds")))

    # most important features
    most_imp <- sort(colMeans(abs(shap$shaps)), decreasing = TRUE)
    # saveRDS(most_imp, file.path(SHAP_DIR, paste0(clean, "_shap_importance.rds")))
    write.csv(data.frame(Feature = names(most_imp),
                         SHAP_Value = most_imp),
              file.path(SHAP_DIR, paste0(clean, "_shap_importance.csv")))
}

# ========================
# 4E-H. SHAP PLOTS
# ========================
for (class_name in unique_classes) {
    clean <- clean_label(class_name)

    # read in
    shp <- readRDS(file.path(SHAP_DIR, paste0(clean, "_shpviz.rds")))
    shap <- readRDS(file.path(SHAP_DIR, paste0(clean, "_shap.rds")))

    # beeswarm
    p <- sv_importance.shapviz(shp, kind = "beeswarm", max_display = 30)
    pdf(file.path(FIG_DIR, paste0(clean, "_ped85_beeswarm.pdf")), height = 6, width = 5)
    plot(p)
    dev.off()

    # contribution plot
    # p <- plot_contribution(shap, obs = 1, max_vars = 30)
    # pdf(file.path(fig_dir, paste0(clean, "_ped85_contribution.pdf")), height = 7, width = 10)
    # plot(p)
    # dev.off()
}

# ========================
# 4H. ENRICHMENT PLOT
# ========================
imp <- readRDS(file.path(MODEL_DIR, "fold_importance.rds"))
sel_probes3k <- rownames(imp %>% head(3000))

x <- testEnrichment(sel_probes3k, "TFBSconsensus", platform = "EPIC")
pdf(file.path(FIG_DIR, "fmodel_driver_enrichment_tfbs.pdf"),width = 2.8, height = 3.9, onefile=FALSE)
plotDotFDR(x, n_min = 30, n_max = 30)
dev.off()
