# ============================================================
#   This script contains the analysis for Fig. 4, including:
#   t-SNE analysis of
#   2D Invasiveness results evaluation
#   2E Driver results evaluation
#   2F LN + primary t-SNE evaluation
# ============================================================
# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- "/home/lijz/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")
IDAT_DIR <- file.path(BASE_DIR, "IDAT")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# LOADING DATA
# ========================
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
betas <- readRDS(file.path(DATA_DIR, "ped_betas_imputed.rds"))

# ========================
# TSNE COORDS 136
# ========================
mtx = t(bSubMostVariable(betas, 3000))

set.seed(12345678)
pca_result <- prcomp(mtx, center = TRUE, scale. = TRUE)
var_explained <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
n_components <- which(var_explained >= 0.8)[1] #27
pca_data <- pca_result$x[, 1:n_components]

set.seed(123467)
tsne <- Rtsne(
    pca_data,
    dims = 2,
    perplexity = 12,
    max_iter = 2000,
    pca = FALSE
)

df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")

saveRDS(df, file.path(DATA_DIR, "ped98_tsne_coords.rds"))

# ========================
# 4F. t-SNE PLOT
# ========================
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx"))
df <- readRDS(here("data", "20250323_thyroid136_tsne_coords.rds"))

ss$tSNE1 <- df$tSNE1
ss$tSNE2 = df$tSNE2
ss$Actual_Age <- as.numeric(ss$Actual_Age)
ss$Confidence <- as.factor(ss$Confidence)

custom_shapes <- c("0" = 16,
                   "1" = 17)

PLOT_PATH <- file.path(here(), "figures")
PLOT_DATE <- "20250410"
PLOT_SIZE <- list(width = 5.5, height = 4.5)

common_theme <- theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray"),
    aspect.ratio = 1,
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "none"
)


