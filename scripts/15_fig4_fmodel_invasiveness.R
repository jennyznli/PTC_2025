# ============================================================
#   FIG 4.
#   This script trains the final random forest classifier for
#   predicting invasiveness in thyroid cancer samples
#   using all 88 primary samples, including:
#   S4D. SHAP beeswarm plot for final model.
#   4G. Enrichment plot of 3k most important probes.
# ============================================================

# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- here()
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

MODEL_DIR <- file.path(DATA_DIR, "final_invasiveness_model")
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
# LOADING DATA
# ========================
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx")) %>%
    dplyr::filter(Lymph_Node == "F")

betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))

valid_ids <- intersect(ss$IDAT, colnames(betas))
ss <- ss[match(valid_ids, ss$IDAT), ]
betas <- betas[, valid_ids]

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
print(model)

# OOB estimate of  error rate: 4.55%
# Confusion matrix:
#     High Low class.error
# High   56   1  0.01754386
# Low     3  28  0.09677419

saveRDS(model, file.path(MODEL_DIR, "final_invasiveness_model.rds"))

# ========================
# 4D. SHAP ANALYSIS
# ========================
train_data <- as.data.frame(t(train_betas3k))
train_labels <- factor(ifelse(as.factor(ss$Clinical_Invasiveness) == "Low", 0, 1),
                       levels = c(0, 1))

model_shp <- randomForest(
    x = train_data,
    y = train_labels,
    ntree = 500,
    mtry = 2,
    nodesize = 5,
    importance = TRUE
)
saveRDS(model_shp, file.path(SHAP_DIR, "invasiveness_shap_model.rds"))

unified_model <- randomForest.unify(model_shp, train_data)
shap <- treeshap(unified_model, train_data, interactions = TRUE)
shp <- shapviz::shapviz(shap, X = train_data)

saveRDS(shap, file.path(SHAP_DIR, "invasiveness_shap.rds"))
saveRDS(shp, file.path(SHAP_DIR, "invasiveness_shpviz.rds"))

# beeswarm plot
p <- sv_importance.shapviz(shp, kind = "beeswarm", max_display = 30)
pdf(file.path(FIG_DIR, "fmodel_invasiveness_beeswarm.pdf"), height = 6, width = 5)
plot(p)
dev.off()

# contribution plot
# p <- plot_contribution(shap, obs = 1, max_vars = 30)
# pdf(file.path(FIG_DIR, "fmodel_invasiveness_contribution.pdf"), height = 7, width = 10)
# plot(p)
# dev.off()

# ========================
# 4G. ENRICHMENT PLOT
# ========================
imp <- readRDS(file.path(MODEL_DIR, "fold_importance.rds"))
sel_probes3k <- imp %>% head(3000) %>% pull(Feature)

x <- testEnrichment(sel_probes3k, "TFBSconsensus", platform = "EPIC")
pdf(file.path(FIG_DIR, "fmodel_invasiveness_enrichment_tfbs.pdf"), width = 2.8, height = 3.9, onefile=FALSE)
plotDotFDR(x, n_min = 30, n_max = 30)
dev.off()

