# ============================================================
#   This script calculates leukocyte fraction, mir200C,
#   epigenetic age using CytoMethIC and Horvath models.
# ============================================================

# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- "/home/lijz/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")
FIG_DIR <- file.path(BASE_DIR, "figures")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# LOADING DATA
# ========================
ss_ped <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_primary <- ss_ped %>% filter(Lymph_Node == "F")
ss_adult <- read_excel(file.path(SS_DIR, "adult_master.xlsx"))

betas_ped <- readRDS(file.path(DATA_DIR, "ped_betas_imputed.rds"))
betas_adult <- readRDS(file.path(DATA_DIR, "adult_betas_imputed.rds"))
betas_joint <- readRDS(file.path(DATA_DIR, "joint_betas_imputed.rds"))

# ========================
# LEUKOCYTE FRACTION
# ========================
lf_epic <- readRDS(file.path(DATA_DIR, "LeukoFrac_EPIC_20240614.rds"))
lf_hm450 <- readRDS(file.path(DATA_DIR, "LeukoFrac_HM450_20240614.rds"))

cat("Predicting leukocyte fractions...\n")

leuko_ped = cmi_predict(betas_ped, lf_epic)
cat("Pediatric predictions completed for", ncol(betas_ped), "samples\n")
leuko_adult = cmi_predict(betas_adult, lf_hm450)
cat("Adult predictions completed for", ncol(betas_adult), "samples\n")
leuko_joint = cmi_predict(betas_joint, lf_hm450)
cat("Joint predictions completed for", ncol(betas_joint), "samples\n")

# Helper function to format and save results
save_leuko_results <- function(leuko_pred, sample_names, filename) {
    df <- as.data.frame(sapply(leuko_pred, function(x) x))
    rownames(df) <- sample_names
    colnames(df) <- c("Leukocyte_Fraction")

    filepath <- file.path(DATA_DIR, filename)
    write.csv(df, filepath)
    cat("Saved", nrow(df), "predictions to", filename, "\n")

    return(df)
}

leuko_ped_df <- save_leuko_results(leuko_ped, colnames(betas_ped), "ped_leuko.csv")
leuko_adult_df <- save_leuko_results(leuko_adult, colnames(betas_adult), "adult_leuko.csv")
leuko_joint_df <- save_leuko_results(leuko_joint, colnames(betas_joint), "joint_leuko.csv")
cat("Completed leukocyte fraction.\n")

# ========================
# MIR200C ESTIMATE
# ========================
cat("Predicting mir200C...\n")

cmi <- readRDS(file.path(DATA_DIR,"MIR200C_EPIC_20240315.rds"))

# rows or columns?? maybe it's 2
mir <- apply(betas_ped, 1, function(x) cmi_predict(x, cmi, lift_over=TRUE))
cat("Pediatric MIR200C predictions completed for", length(mir), "samples\n")

df <- as.data.frame(do.call(rbind, lapply(mir, as.data.frame)))
rownames(df) <- colnames(betas_ped)
write.csv(df, file.path(DATA_DIR, "ped98_mir200.csv"))
# add the IDATs here??

cat("Completed mir200c analysis.\n")

# ========================
# EPIGENETIC AGE
# ========================
cat("Predicting epigenetic age...\n")

age <- readRDS(here("prev_data", "AgeClock_HumanHorvathN353_HM450.rds"))
betas_ped <- readRDS(file.path(DATA_DIR, "ped_betas_imputed.rds")) %>% mLiftOver("HM450")
cat("Pediatric EPICv2 -> HM450 lifted over.\n")
betas_adult <- readRDS(file.path(DATA_DIR, "adult_betas_imputed.rds"))
cat("Adult HM450 loaded.\n")
betas_joint <- readRDS(file.path(DATA_DIR, "joint_betas_imputed.rds")) %>% mLiftOver("HM450") # does this work??
cat("Joint ?? -> HM450 lifted over.\n")

df_ped <- apply(betas_ped, 2, function(column) predictAge(column, age))
rownames(df_ped) <- colnames(betas_ped)
cat("Pediatric calculated\n")
df_adult <- apply(betas_adult, 2, function(column) predictAge(column, age))
rownames(df_adult) <- colnames(betas_adult)
cat("Adult calculated.\n")
df_joint <- apply(betas_joint, 2, function(column) predictAge(column, age))
rownames(df_joint) <- colnames(betas_joint)
cat("Joint calculated.\n")

write.csv(df_ped, file.path(DATA_DIR, "ped98_age_horvath.csv"), row.names = TRUE)
write.csv(df_adult, file.path(DATA_DIR, "adult_age_horvath.csv"), row.names = TRUE)
write.csv(df_joint, file.path(DATA_DIR, "joint_age_horvath.csv"), row.names = TRUE) # check this??
cat("Completed epigenetic age.\n")

# age_predictions88 <- apply(betas88, 2, function(column) predictAge(column, model))
# write.csv(age_predictions88, here("data", "20240915_thyroid136_age_predictions_horvath.csv"), row.names = TRUE)
#
# age_predictions496 <- apply(betas496, 2, function(column) predictAge(column, model))
# write.csv(age_predictions496, here("data", "20240915_thca496_age_predictions_horvath.csv"), row.names = TRUE)
#


