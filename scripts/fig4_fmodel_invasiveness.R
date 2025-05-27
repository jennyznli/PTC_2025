# ============================================================
#   FIG 4.
#   This script trains the final random forest classifier for
#   predicting invasiveness in thyroid cancer samples
#   using all 88 primary samples.
# ============================================================

# ========================
# CONFIGURATION
# ========================
set.seed(123)
N_WORKERS <- 20

BASE_DIR <- here()
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

MODEL_DIR <- file.path(DATA_DIR, "final_invasiveness_model")
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
# LOADING DATA
# ========================
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx")) %>%
    dplyr::filter(Lymph_Node == "F")

betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))

valid_ids <- intersect(ss$IDAT, colnames(betas))
ss <- ss[match(valid_ids, ss$IDAT), ]
betas <- betas[, valid_ids]

# Labels
labels <- as.factor(ss$Clinical_Invasiveness)
train_labels <- labels

# ========================
# FEATURE SELECTION
# ========================
n_probes <- nrow(betas)
n_batches <- ceiling(n_probes / BATCH_SIZE)
fold_models <- list()
fold_importance <- data.frame()

for (s in 0:(n_batches - 1)) {
    set.seed(123 + s)

    # Calculate indices for current batch
    start_idx <- s * BATCH_SIZE + 1
    end_idx <- min((s + 1) * BATCH_SIZE, n_probes)
    sbetas <- betas[start_idx:end_idx, ]

    model <- randomForest(
        x = t(sbetas),
        y = train_labels,
        ntree = N_TREES,
        importance = TRUE
    )

    fold_models[[s + 1]] <- model
    importance <- data.frame(
        Model = s,
        Feature = rownames(model$importance),
        MeanDecreaseAccuracy = model$importance[, "MeanDecreaseAccuracy"],
        MeanDecreaseGini = model$importance[, "MeanDecreaseGini"]
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
sel_probes3k <- imp %>% head(3000) %>% pull(Feature)
train_betas3k <- betas[sel_probes3k, ]

# TRAIN MODEL
model <- randomForest(
    x = t(train_betas3k),
    y = train_labels,
    mtry = 2,
    nodesize = 5,
    importance = TRUE
)

# OOB estimate of  error rate: 10.23%
# High Low class.error
# High   53   2  0.03636364
# Low     7  26  0.21212121

saveRDS(model, file.path(MODEL_DIR, "final_invasiveness_model.rds"))
