# ============================================================
#   This script contains the analysis for Fig. 3, including:
#   Joint t-SNE embedding of adult and pediatric cohorts
#
#   3A
#   3B
#   3C
#   3D

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
# TSNE COORDS
# ========================
betas = readRDS(here("data", "20250411_combined584_betas_processed.rds"))
ss <- read_xlsx(here("ss", "202504011_combined584.xlsx"))

mtx = t(bSubMostVariable(betas, 3000)) # 584 3000

set.seed(12345678)
pca_result <- prcomp(mtx, center = TRUE, scale. = TRUE)
var_explained <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
n_components <- which(var_explained >= 0.8)[1] #25
pca_data <- pca_result$x[, 1:n_components]

set.seed(123467)
tsne <- Rtsne(
    pca_data,
    dims = 2,
    perplexity = 24,
    max_iter = 2000,
    pca = FALSE
)

df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
ss$tSNE1 <- df$tSNE1
ss$tSNE2 = df$tSNE2

saveRDS(df, here("data", "20250411_combined584_tsne_coords.rds"))
write_xlsx(ss, here("ss", "202504011_combined584.xlsx"))

# ========================
# TSNE PLOTTING
# ========================
ss <- read_xlsx(here("ss", "202504011_combined584.xlsx"))
ss$Source <- ifelse(ss$Source == "THCA", "ADULT", "PED")
DATE = "20250411"

pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_invasiveness.pdf")), width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Sample_Group, shape = Source)) +
    scale_color_manual(values = cluster_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_drivers.pdf")), width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Driver_Group, shape = Source,
                              text = paste("Sample ID:", Sample_ID))) +
    scale_color_manual(values = driver_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

# ggplotly(p)

pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_sex.pdf")), width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Sex, shape = Source,
                              text = paste("Sample ID:", Sample_ID))) +
    scale_color_manual(values = sex_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()


pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_leuko.pdf")), width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Leukocyte_Fraction, shape = Source)) +
    scale_shape_manual(values = custom_shapes) +
    scale_color_gradientn(colors=parula(20)) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", color = NA),  # Added color = NA
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_N.pdf")), width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = N, shape = Source,
                              text = paste("Sample ID:", Sample_ID))) +
    scale_color_manual(values = n_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_M.pdf")), width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = M, shape = Source,
                              text = paste("Sample ID:", Sample_ID))) +
    scale_color_manual(values = m_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

# ========================
# NO LEGENDS
# ========================
# tSNE plot - Invasiveness
pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_invasiveness_nolegend.pdf")), width=5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Sample_Group, shape = Source)) +
    scale_color_manual(values = cluster_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(
        legend.position = "none",
        panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray")
    )
plot(p)
dev.off()

# tSNE plot - Drivers
pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_drivers_nolegend.pdf")), width=5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Driver_Group, shape = Source,
                              text = paste("Sample ID:", Sample_ID))) +
    scale_color_manual(values = driver_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(
        legend.position = "none",
        panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray")
    )
plot(p)
dev.off()

# tSNE plot - Sex
pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_sex_nolegend.pdf")), width=5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Sex, shape = Source,
                              text = paste("Sample ID:", Sample_ID))) +
    scale_color_manual(values = sex_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(
        legend.position = "none",
        panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray")
    )
plot(p)
dev.off()

# tSNE plot - Leukocyte Fraction
pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_leuko_nolegend.pdf")), width=5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Leukocyte_Fraction, shape = Source)) +
    scale_shape_manual(values = custom_shapes) +
    scale_color_gradientn(colors = parula(20)) +
    theme(
        legend.position = "none",
        panel.background = element_rect(fill = "#ffffff", color = NA),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray")
    )
plot(p)
dev.off()


pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_N_nolegend.pdf")), width=5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = N, shape = Source,
                              text = paste("Sample ID:", Sample_ID))) +
    scale_color_manual(values = n_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(
        legend.position = "none",
        panel.background = element_rect(fill = "#ffffff", color = NA),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray")
    )
plot(p)
dev.off()

pdf(file.path(here(), "figures", paste0(DATE, "_combined584_tsne_M_nolegend.pdf")), width=5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = M, shape = Source,
                              text = paste("Sample ID:", Sample_ID))) +
    scale_color_manual(values = m_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(
        legend.position = "none",
        panel.background = element_rect(fill = "#ffffff", color = NA),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray")
    )
plot(p)
dev.off()


