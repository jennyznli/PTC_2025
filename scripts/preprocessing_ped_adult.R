# ============================================================
#   This script performs preprocessing for
#   pediatric + adult betas.
# ============================================================

# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20
cat("Using", N_WORKERS, "parallel workers\n")

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
betas_ped <- readRDS(file.path(DATA_DIR, "ped_betas_processed.rds"))
cat("Original uncollapsed pediatric EPICv2 beta dimensions:", dim(betas_ped), "\n")
betas_ped <- betasCollapseToPfx(betas_ped, BPPARAM = BiocParallel::MulticoreParam(N_WORKERS))
cat("Collapsed EPICv2 beta dimensions:", dim(betas_ped), "\n")

betas_adult <- readRDS(file.path(DATA_DIR, "adult_betas_processed.rds"))
cat("Original HM450 beta dimensions:", dim(betas_adult), "\n")

ss_ped <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_adult <- read_excel(file.path(SS_DIR, "adult_master.xlsx"))

# ========================
# SEX CHROMOSOMES
# ========================
hm450_mft <- sesameAnno_readManifestTSV("HM450.hg38.manifest")
epicv2_mft <- sesameAnno_readManifestTSV("EPICv2.hg38.manifest")

# remove sex chromosome probes
hm450_sex <- hm450_mft %>%
    filter(CpG_chrm %in% c("chrX", "chrY"))
probes_hm450_sex <- hm450_sex$Probe_ID
cat("Found", length(probes_hm450_sex), "HM450 sex probes.\n")

epicv2_sex <- epicv2_mft %>%
    filter(CpG_chrm %in% c("chrX", "chrY"))
probes_epicv2_sex <- epicv2_sex$Probe_ID

# clean EPICv2 probe names
probes_epicv2_sex <- clean_probe_ids(probes_epicv2_sex)
cat("Found", length(probes_epicv2_sex), "EPICv2 sex probes.\n")

# ========================
# GET MASKS
# ========================
hm450_mask <- getMask(platform = "HM450", mask_names = "recommended")
cat("Found", length(hm450_mask), "HM450 probes to be masked.\n")

epicv2_mask <- getMask(platform = "EPICv2", mask_names = "recommended")
epicv2_mask <- clean_probe_ids(epicv2_mask)
cat("Found", length(epicv2_mask), "EPICv2 probes to be masked.\n")

remove_probes <- unique(c(hm450_mask, epicv2_mask,
                          probes_epicv2_sex, probes_hm450_sex))
cat("Found", length(remove_probes), "total EPICv2 & HM450 probes to be masked.\n")

# ========================
# COMBINING BETAS
# ========================
ped_probes <- rownames(betas_ped)
adult_probes <- rownames(betas_adult)
cat("Pediatric probe count:", length(ped_probes), "\n")
cat("Adult probe count:", length(adult_probes), "\n")

joint_probes <- intersect(ped_probes, adult_probes)
cat("Overlapping probes:", length(joint_probes), "\n")

diff_probes <- setdiff(joint_probes, remove_probes)
cat("Probes remaining after sex, masking:", length(diff_probes), "\n")

# probe types
probe_types <- table(substr(diff_probes, 1, 2))
cat("Probe types before filtering:\n")
print(probe_types)

# filter to only cg probes
cg_probes <- diff_probes[substr(diff_probes, 1, 2) == "cg"]
cat("CG only probes:", length(cg_probes), "\n")

joint_betas <- cbind(betas_ped[cg_probes,], betas_adult[cg_probes,])
cat("Joint beta matrix dimensions:", dim(joint_betas), "\n")

saveRDS(joint_betas, file.path(DATA_DIR, "joint_betas_processed.rds"))

# ========================
# IMPUTATION
# ========================
# clean and impute
joint_betas_clean <- cleanMatrixForClusterW(joint_betas)
joint_betas_imputed <- impute(joint_betas_clean, "HM450")

# save primary + tumor
cat("Joint beta matrix dimensions:", dim(joint_betas_imputed), "\n")
saveRDS(joint_betas_imputed, file.path(DATA_DIR, "joint_betas_imputed.rds"))
