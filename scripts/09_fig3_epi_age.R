# ============================================================
#   FIG 3.
#   This script contains the script for epigenetic age analysis.
#   3E. Actual v predicted age.
#   3F. Age acceleration between invasiveness box plots
# ============================================================

# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- here()
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")
IDAT_DIR <- file.path(BASE_DIR, "IDAT")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# 3E. ACTUAL V. PREDICTED AGE
# ========================
ss_ped <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_adult <- read_excel(file.path(SS_DIR, "adult_master.xlsx"))

age_ped <- read.csv(file.path(DATA_DIR, "ped_age_horvath.csv"))
age_adult <- read.csv(file.path(DATA_DIR, "adult_age_horvath.csv"))

ss_ped$Age <- as.numeric(ss_ped$Age)
ss_ped$Predicted_Age <- as.numeric(age_ped[,2])
ss_adult$Age <- as.numeric(ss_adult$Age)
ss_adult$Predicted_Age <- as.numeric(age_adult[,2])

ss_ped <- ss_ped %>%
    select('Source',"Methylation_Clusters", 'Sample_ID', 'Age', 'Predicted_Age', 'Driver_Group', 'Clinical_Invasiveness', 'Sex')
ss_adult <- ss_adult %>%
    select('Source',"Methylation_Clusters", 'Sample_ID', 'Age', 'Predicted_Age', 'Driver_Group', 'Clinical_Invasiveness', 'Sex')

ss_ped$Fold_Change <- ss_ped$Predicted_Age / ss_ped$Age
ss_adult$Fold_Change <- ss_adult$Predicted_Age / ss_adult$Age

ss = rbind(ss_ped, ss_adult)

p <- ggplot(ss, aes(x = Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Methylation_Clusters, shape = Source), size = 1) +
    scale_color_manual(values = cluster_colors) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", size = 0.3) +
    scale_shape_manual(values = custom_shapes) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(file.path(FIG_DIR, "joint_actual_predicted_age_cluster.pdf"),
       plot = p, width = 6, height = 5)

# ========================
# 3F. BOX PLOT AGE ACCELERATION
# ========================
ss <- ss %>%
    mutate(
        Age_Acceleration = Predicted_Age - Age,
        Age_Acceleration_Residual = residuals(lm(Predicted_Age ~ Age))
    )

format_pval <- function(p) {
    if(p < 0.001) return("p < 0.001")
    return(sprintf("p = %.3f", p))
}

# remove the NAs....
ss_filter <- ss %>%
    filter(Clinical_Invasiveness != "NA")

stats_list <- lapply(unique(ss$Source), function(src) {
    data_subset <- ss_filter[ss_filter$Source == src, ]
    test <- wilcox.test(Age_Acceleration_Residual ~ Clinical_Invasiveness, data = data_subset)
    data.frame(
        Source = src,
        p.value = test$p.value,
        label = format_pval(test$p.value),
        y.position = max(data_subset$Age_Acceleration_Residual, na.rm = TRUE) + 2
    )
})
stats_df <- do.call(rbind, stats_list)

p <- ggplot(ss_filter, aes(x = Source, y = Age_Acceleration_Residual, fill = Clinical_Invasiveness)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                alpha = 0.3, size = 0.6) +
    scale_fill_manual(values = invasiveness_colors) +
    theme_minimal() +
    labs(x = "Cohort",
         y = "Age Acceleration (Years)") +
    theme(axis.text.x = element_text(hjust = 1)) +
    geom_text(data = stats_df,
              aes(x = Source, y = y.position, label = label),
              inherit.aes = FALSE,
              size = 3)
ggsave(file.path(FIG_DIR, "age_acceleration_invasiveness_cohort.pdf"),
       plot = p, width=5, height=4.5)

print("Wilcoxon test results by cohort:")
print(stats_df)

# Source   p.value     label y.position
# 1    PED 0.1023640 p = 0.102   32.27781
# 2  ADULT 0.5261311 p = 0.526   59.86178

# ORIGINAL VALUES...
# Source   p.value     label y.position
# 1    PED 0.0385872 p = 0.039   43.19281
# 2   THCA 0.2089089 p = 0.209   71.14655

# ========================
# FIG. 3E - JOINT
# ========================
ss_joint <- read_excel(file.path(SS_DIR, "joint_master.xlsx"))
age_joint <- read.csv(file.path(DATA_DIR, "joint_age_horvath.csv"))
indices <- match(ss_joint$Sample_ID, age_joint[,1])
age_joint <- age_joint[indices, , drop = FALSE]

ss_joint$Age <- as.numeric(ss_joint$Age)
ss_joint$Predicted_Age <- as.numeric(age_joint[,2])

ss_joint$Fold_Change <- ss_joint$Predicted_Age / ss_joint$Age

p <- ggplot(ss_joint, aes(x = Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Methylation_Clusters, shape = Source), size = 1) +
    scale_color_manual(values = cluster_colors) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", size = 0.3) +
    scale_shape_manual(values = custom_shapes) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(file.path(FIG_DIR, "joint_actual_predicted_age_cluster_2.pdf"),
       plot = p, width = 6, height = 5)

# ========================
# 3F. BOX PLOT AGE ACCELERATION
# ========================
ss_joint <- ss_joint %>%
    mutate(
        Age_Acceleration = Predicted_Age - Age,
        Age_Acceleration_Residual = residuals(lm(Predicted_Age ~ Age))
    )

format_pval <- function(p) {
    if(p < 0.001) return("p < 0.001")
    return(sprintf("p = %.3f", p))
}

# remove the NAs....
ss_filter <- ss_joint %>%
    filter(Clinical_Invasiveness != "NA")

stats_list <- lapply(unique(ss_filter$Source), function(src) {
    data_subset <- ss_filter[ss_filter$Source == src, ]
    test <- wilcox.test(Age_Acceleration_Residual ~ Clinical_Invasiveness, data = data_subset)
    data.frame(
        Source = src,
        p.value = test$p.value,
        label = format_pval(test$p.value),
        y.position = max(data_subset$Age_Acceleration_Residual, na.rm = TRUE) + 2
    )
})
stats_df <- do.call(rbind, stats_list)

p <- ggplot(ss_filter, aes(x = Source, y = Age_Acceleration_Residual, fill = Clinical_Invasiveness)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                alpha = 0.3, size = 0.6) +
    scale_fill_manual(values = invasiveness_colors) +
    theme_minimal() +
    labs(x = "Cohort",
         y = "Age Acceleration (Years)") +
    theme(axis.text.x = element_text(hjust = 1)) +
    geom_text(data = stats_df,
              aes(x = Source, y = y.position, label = label),
              inherit.aes = FALSE,
              size = 3)
ggsave(file.path(FIG_DIR, "age_acceleration_invasiveness_cohort2.pdf"),
       plot = p, width=5, height=4.5)

print("Wilcoxon test results by cohort:")
print(stats_df)

# Source   p.value     label y.position
# 1    PED 0.2052922 p = 0.205   26.98744
# 2  ADULT 0.4667626 p = 0.467   58.99673

