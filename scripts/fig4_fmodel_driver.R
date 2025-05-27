# ============================================================
#   FIG 4.
#   This script trains the final random forest classifier for
#   predicting driver group in thyroid cancer samples
#   using all 85 primary samples (excluding 3 indeterminate).
# ============================================================

# ========================
# CONFIGURATION
# ========================
set.seed(123)
N_WORKERS <- 20
DATE <- format(Sys.Date(), "%Y%m%d")

BASE_DIR <- here()
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

MODEL_DIR <- file.path(DATA_DIR, "final_driver_model")
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

print("Configuration complete.")

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

# Labels
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
    set.seed(123 + s)

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

saveRDS(fold_models, file.path(MODEL_DIR, paste0(DATE, "_fold_models.rds")))
saveRDS(features, file.path(MODEL_DIR, paste0(DATE, "_fold_importance.rds")))

# ========================
# MODEL DEVELOPMENT
# ========================
imp <- readRDS(file.path(MODEL_DIR, paste0(DATE, "_fold_importance.rds")))
sel_probes3k <- rownames(imp %>% head(3000))
train_betas3k <- betas_train[sel_probes3k, ]

# TRAIN MODEL
set.seed(123)
model <- randomForest(
    x = t(train_betas3k),
    y = train_labels,
    mtry = 2,
    nodesize = 5,
    importance = TRUE
)

# OOB estimate of  error rate: 5.88%
# Confusion matrix:
#               BRAF V600E DICER1 Kinase Fusion Ras-like class.error
# BRAF V600E            21      0             0        1  0.04545455
# DICER1                 0      7             0        1  0.12500000
# Kinase Fusion          1      0            33        1  0.05714286
# Ras-like               0      1             0       19  0.05000000

saveRDS(model, file.path(MODEL_DIR, paste0(DATE, "_final_driver_model.rds")))

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

write.csv(results_df, file.path(MODEL_DIR, paste0(DATE, "_indeterminate_predictions.csv")), row.names = FALSE)

# 207686140028_R03C01 207686140033_R02C01 207700160043_R03C01
# THY0138T            THY0170BT           THY0066AT
# Kinase Fusion       Kinase Fusion       Ras-like

