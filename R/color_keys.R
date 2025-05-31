# ========================
# COLOR & SHAPE KEYS
# ========================

methylation_colors <- c(
    "Hypo" = "#0c66bcff",
    "Hyper" = "#d61525ff"
)
cluster_colors <- c(
    "HI" = "#f8766dff",
    "HIL" = "#f7b456ff",
    "LI" = "#00bfc4ff",
    "DICER1" = "#b47cffff",
    "NA" = "lightgray"
)

invasiveness_colors <- c(
    "High" = "#f8766dff",
    "Low" = "#00bfc4ff",
    "NA" = "lightgray"
)

driver_colors <- c(
    "BRAF V600E" = "#ff6cc3ff",
    "Kinase Fusion" = "#20bb20ff",
    "Ras-like" = "#00b4f0ff",
    "DICER1" = "#b47cffff",
    "Indeterminate" = "gray",
    "NA" = "lightgray"
)

lymph_node_colors <- c(
    "T" = "#363636",
    "F" = "lightgray"
)


sex_colors <- c(
    "Female" = "#8c44f6ff",
    "Male" = "#00cd92ff"
)

t_colors <- c(
    "T0" = "lightgray",
    "T1" = "#ffe811ff",
    "T1a" = "#ffad29ff",
    "T1b" = "#f16500ff",
    "T2" = "#1cbf00ff",
    "T3" = "#00d3e7ff",
    "T3a" = "#2972ffff",
    "T3b" = "#9637ffff",
    "T4" = "#f263c0ff",
    "T4a" = "#c40cecff",
    "NA" = "lightgray"
)

n_colors <- c(
    "N0" = "#ffe93fff",
    "N1" = "#fc530a",
    "N1a" = "#ff8a04ff",
    "N1b" = "#fc0f0bff",
    "NX" = "lightgray",
    "NA" = "lightgray"
)

m_colors <- c(
    "M0" = "#fc71c5ff",
    "M1" = "#9808caff",
    "MX" = "lightgray",
    "NA" = "lightgray"
)

pair_colors <- c("#4f2288", "#538fcb", "#88CCEE", "#44AA99", "#117733", "#c9a90e", "#769b26",
                 "#FF7F0EFF", "#cf4a17", "#9e0b22", "#AA4499", "#6F63BBFF")

feature_colors <- c(
    "Common" = "blue",
    "Driver-specific" = "darkgreen",
    "Invasiveness-specific" = "purple"
)

custom_shapes <- c("ADULT" = 16,  # Circle
                   "PED" = 1)  # Empty circle

accuracy_shapes <- c("0" = 16,  # Circle
                   "1" = 17, # Triangle
                   "NA" = 16)   # Circle


