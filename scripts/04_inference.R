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
betas_ped_uncollapsed <- readRDS(file.path(DATA_DIR, "ped_betas_imputed_uncollapsed.rds"))
betas_adult <- readRDS(file.path(DATA_DIR, "adult_betas_imputed.rds"))
betas_joint <- readRDS(file.path(DATA_DIR, "joint_betas_imputed.rds"))

# ========================
# EPIGENETIC AGE
# ========================
cat("Predicting epigenetic age...\n")

age_epicv2 <- readRDS(url("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/refs/heads/main/Anno/EPICv2/Clock_Horvath353.EPICv2.345.rds"))
age_hm450 <- readRDS(url("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/refs/heads/main/Anno/HM450/Clock_Horvath353.rds"))

df_ped <- as.matrix(apply(betas_ped_uncollapsed, 2, function(column) predictAge(column, age_epicv2)))
cat("Pediatric age predictions completed for", ncol(betas_ped_uncollapsed), "samples\n")
df_adult <- apply(betas_adult, 2, function(column) predictAge(column, age_hm450))
cat("Adult age predictions completed for", ncol(betas_adult), "samples\n")
df_joint <- apply(betas_joint, 2, function(column) predictAge(column, age_hm450))
cat("Joint age predictions completed for", ncol(betas_joint), "samples\n")

write.csv(df_ped, file.path(DATA_DIR, "ped_age_horvath.csv"), row.names = TRUE)
write.csv(df_adult, file.path(DATA_DIR, "adult_age_horvath.csv"), row.names = TRUE)
write.csv(df_joint, file.path(DATA_DIR, "joint_age_horvath.csv"), row.names = TRUE)

cat("Completed epigenetic age.\n")

# ========================
# LEUKOCYTE FRACTION
# ========================
cat("Predicting leukocyte fractions...\n")

lf_epic <- readRDS(url("https://github.com/zhou-lab/CytoMethIC_models/raw/refs/heads/main/models/LeukoFrac_EPIC_20240614.rds"))
lf_hm450 = readRDS(url("https://github.com/zhou-lab/CytoMethIC_models/raw/refs/heads/main/models/LeukoFrac_HM450_20240614.rds"))

cat("Lifting over pediatric EPICv2 betas to EPIC...\n")
betas_ped_lifted <- mLiftOver(betas_ped_uncollapsed, "EPIC")
# cat("Lifting over joint betas to HM450...\n")
# betas_joint_lifted <- mLiftOver(betas_joint, "EPIC")

leuko_ped <- cmi_predict(betas_ped_lifted, lf_epic)
cat("Pediatric leukocyte predictions completed for", ncol(betas_ped_lifted), "samples\n")
leuko_adult <- cmi_predict(betas_adult, lf_hm450)
cat("Adult leukocyte predictions completed for", ncol(betas_adult), "samples\n")
leuko_joint <- cmi_predict(betas_joint, lf_hm450)
cat("Joint leukocyte predictions completed for", ncol(betas_joint), "samples\n")

leuko_ped_df <- save_cytomethic_results(leuko_ped, "ped_leuko.csv", colnames(betas_ped_uncollapsed))
leuko_adult_df <- save_cytomethic_results(leuko_adult, "adult_leuko.csv", colnames(betas_adult))
leuko_joint_df <- save_cytomethic_results(leuko_joint, "joint_leuko.csv", colnames(betas_joint))

cat("Completed leukocyte fraction.\n")

# ========================
# MIR200C ESTIMATE
# ========================
cat("Predicting mir200C...\n")

cmi <- readRDS(url("https://github.com/zhou-lab/CytoMethIC_models/raw/refs/heads/main/models/MIR200C_EPIC_20240315.rds"))

# cat("Lifting over pediatric EPICv2 betas to EPIC...\n") - already done
# betas_ped_lifted <- mLiftOver(betas_ped_uncollapsed, "EPIC")
cat("Lifting over adult HM450 betas to EPIC...\n")
betas_adult_lifted <- mLiftOver(betas_adult, "EPIC")
cat("Lifting over joint betas to EPIC...\n")
betas_joint_lifted <- mLiftOver(betas_joint, "EPIC")

mir_ped <- cmi_predict(betas_ped_lifted, cmi)
cat("Pediatric MIR200C predictions completed for", ncol(betas_ped_lifted), "samples\n")
mir_adult <- cmi_predict(betas_adult_lifted, cmi)
cat("Adult MIR200C predictions completed for", ncol(betas_adult_lifted), "samples\n")
mir_joint <- cmi_predict(betas_joint_lifted, cmi)
cat("Joint MIR200C predictions completed for", ncol(betas_joint_lifted), "samples\n")

mir_ped_df <- save_cytomethic_results(mir_ped, "ped_mir200.csv", colnames(betas_ped))
mir_adult_df <- save_cytomethic_results(mir_adult, "adult_mir200.csv", colnames(betas_adult))
mir_joint_df <- save_cytomethic_results(mir_joint, "joint_mir200.csv", colnames(betas_joint))

cat("Completed mir200C analysis.\n")



