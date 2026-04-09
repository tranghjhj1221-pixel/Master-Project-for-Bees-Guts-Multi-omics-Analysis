library(tibble)
library(tidyverse)
library(readxl)
library(reshape2)

# load file
dat <- read.csv(file.choose(), check.names = FALSE)
colnames(dat)

library(dplyr)

dat.sub <- dat %>%
  select(
    compound = `Compound_name`, #I manually create this column in excel to create a consistent nonmenclature names for lipid
    score = `preferred_annotation:score`,   # ← ADD THIS
    adduct = `preferred_annotation:adduct`,
    mz,
    rt,
    iin = `ion_identities:iin_id`,
    starts_with("datafile:Crop"),
    starts_with("datafile:Mid"),
    starts_with("datafile:Hind")
  )
View(dat.sub) #check out this file to see which annotation is the best
#choosing the best annotation because every ion in the same ion identity network should be assigned to the same compound name (ideally)

#pick the problematic ones (which has multiple compound names appear)
conflicts <- dat.sub %>%
  group_by(iin) %>%
  summarise(n_names = n_distinct(compound, na.rm = TRUE)) %>%
  filter(n_names > 1)

#arrange iin in a desc order interm of annotation score to inspec manually 
View(dat.sub %>%
  filter(iin %in% conflicts$iin) %>%
  arrange(iin, desc(score)))

#gonna collapse the feature in the same network as one, inspec manually showed some problematic feature, ask Ali about iid == 23


# ==============================
# SPLIT DATA
# ==============================
dat.safe <- dat.sub %>%
  filter(!iin %in% conflicts$iin)

dat.conflict <- dat.sub %>%
  filter(iin %in% conflicts$iin)

dat.no_iin <- dat.sub %>%
  filter(is.na(iin) | iin == "")

# ==============================
# COLLAPSE SAFE NETWORKS
# ==============================
area_cols <- grep(":area$", colnames(dat.sub), value = TRUE)

dat.safe <- dat.safe %>%
  mutate(total_area = rowSums(across(all_of(area_cols)), na.rm = TRUE))

dat.collapsed <- dat.safe %>%
  filter(!is.na(iin) & iin != "") %>%
  group_by(iin) %>%
  arrange(desc(score), desc(total_area), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# ==============================
# COMBINE EVERYTHING
# ==============================
dat.final <- bind_rows(dat.collapsed, dat.conflict, dat.no_iin)



# ==============================
# BUILD HEATMAP MATRIX
# ==============================
area_cols <- grep(":area$", colnames(dat.final), value = TRUE)

dat.hm <- dat.final[, area_cols]
dat.hm <- as.data.frame(dat.hm)

# ==============================
# USE COMPOUND + m/z + RT AS ROWNAMES
# ==============================
rownames(dat.hm) <- make.unique(paste0(
  ifelse(is.na(dat.final$compound), "Unknown", dat.final$compound),
  " | m/z=", round(dat.final$mz, 2),
  " | RT=", round(dat.final$rt, 2)
))

# ==============================
# ENSURE NUMERIC
# ==============================
dat.hm[] <- lapply(dat.hm, as.numeric)

# ==============================
# LOG2 TRANSFORM
# ==============================
dat.hm[dat.hm == 0] <- NA
dat.hm.log <- log2(dat.hm)

# ==============================
# HISTOGRAM (QC)
# ==============================
hist(unlist(dat.hm.log),
     breaks = 50,
     main = "Log2 Intensity Distribution",
     xlab = "log2(area)")
#histogram looks gucci

# ==============================
# BUILD HEATMAP MATRIX
# ==============================
area_cols <- grep(":area$", colnames(dat.final), value = TRUE)

dat.hm <- dat.final[, area_cols]
dat.hm <- as.data.frame(dat.hm)

# ==============================
# ENSURE NUMERIC
# ==============================
dat.hm[] <- lapply(dat.hm, as.numeric)

# ==============================
# LOG2 TRANSFORM
# ==============================
dat.hm[dat.hm == 0] <- NA
dat.hm.log <- log2(as.matrix(dat.hm))

# ==============================
# SELECT TISSUE COLUMNS
# ==============================
dat.hm.log <- dat.hm.log[, c(
  "datafile:Crop_P1-E-1_1_14911.d:area",
  "datafile:Mid_P1-E-2_1_14912.d:area",
  "datafile:Hind_P1-E-3_1_14913.d:area"
)]

colnames(dat.hm.log) <- c("Crop", "Midgut", "Hindgut")

# ==============================
# KEEP ROWS WITH AT LEAST 2 OBSERVED VALUES
# fixes clustering issues from rows that are all NA
# or have only 1 non-NA value
# ==============================
dat.hm.log <- dat.hm.log[rowSums(!is.na(dat.hm.log)) >= 2, , drop = FALSE]

# ==============================
# MANUAL ROW SCALING (IGNORE NA)
# ==============================
dat.hm.scaled <- t(apply(dat.hm.log, 1, function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - m) / s
}))

colnames(dat.hm.scaled) <- colnames(dat.hm.log)

# ==============================
# REMOVE ROWS THAT BECOME ALL NA AFTER SCALING
# ==============================
dat.hm.scaled <- dat.hm.scaled[rowSums(!is.na(dat.hm.scaled)) >= 2, , drop = FALSE]

# ==============================
# OPTIONAL SAFETY CHECK
# ==============================
if (nrow(dat.hm.scaled) < 2) {
  stop("Not enough rows with at least 2 non-NA values to cluster.")
}

# ==============================
# HEATMAP
# ==============================
library(pheatmap)
library(viridis)

pheatmap(
  dat.hm.scaled,
  scale = "none",
  color = viridis(100),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE,
  na_col = "grey90"
)



# ==============================
# REBUILD CLEAN LOG2 MATRIX
# IDENTIFY LIPIDS ENRICHED IN EACH REGION
# SUMMARIZE TOP DRIVERS
# ==============================

# 1) Extract the 3 tissue area columns
area_cols <- c(
  "datafile:Crop_P1-E-1_1_14911.d:area",
  "datafile:Mid_P1-E-2_1_14912.d:area",
  "datafile:Hind_P1-E-3_1_14913.d:area"
)

dat.hm <- dat.final[, area_cols, drop = FALSE]
dat.hm <- as.data.frame(dat.hm)

# 2) Make sure values are numeric
dat.hm[] <- lapply(dat.hm, function(x) as.numeric(as.character(x)))

# 3) Set row names from compound names
rownames(dat.hm) <- make.unique(
  ifelse(is.na(dat.final$compound) | dat.final$compound == "",
         "Unknown",
         dat.final$compound)
)

# 4) Replace zeros with NA, then log2 transform
dat.hm[dat.hm == 0] <- NA
dat.hm.log <- log2(as.matrix(dat.hm))

# 5) Rename columns for convenience
colnames(dat.hm.log) <- c("Crop", "Midgut", "Hindgut")

# 6) Keep rows with at least 1 observed value
x <- dat.hm.log[rowSums(!is.na(dat.hm.log)) > 0, , drop = FALSE]

# 7) For each lipid, assign the region with the highest log2 abundance
enriched_in <- apply(x, 1, function(v) {
  if (all(is.na(v))) return(NA_character_)
  colnames(x)[which.max(v)]
})

# 8) Use max log2 abundance as driver strength
driver_strength <- apply(x, 1, function(v) {
  if (all(is.na(v))) return(NA_real_)
  max(v, na.rm = TRUE)
})

# 9) Build summary table
driver_df <- data.frame(
  lipid = rownames(x),
  enriched_in = enriched_in,
  lipid_class = sub("\\(.*", "", rownames(x)),
  driver_strength = driver_strength,
  stringsAsFactors = FALSE
)

driver_df <- driver_df[!is.na(driver_df$enriched_in), , drop = FALSE]

# 10) Inspect counts by region
cat("\nCounts by region:\n")
print(table(driver_df$enriched_in, useNA = "ifany"))

cat("\nCounts by lipid class and region:\n")
print(table(driver_df$lipid_class, driver_df$enriched_in, useNA = "ifany"))

# 11) Top 10 strongest lipids per region
top_crop <- head(
  driver_df[driver_df$enriched_in == "Crop", ][
    order(-driver_df[driver_df$enriched_in == "Crop", ]$driver_strength),
  ],
  10
)

top_mid <- head(
  driver_df[driver_df$enriched_in == "Midgut", ][
    order(-driver_df[driver_df$enriched_in == "Midgut", ]$driver_strength),
  ],
  10
)

top_hind <- head(
  driver_df[driver_df$enriched_in == "Hindgut", ][
    order(-driver_df[driver_df$enriched_in == "Hindgut", ]$driver_strength),
  ],
  10
)

cat("\nTop Crop drivers:\n")
print(top_crop)

cat("\nTop Midgut drivers:\n")
print(top_mid)

cat("\nTop Hindgut drivers:\n")
print(top_hind)

# ==============================
# COMBINE TOP DRIVERS INTO ONE TABLE
# ==============================

top_crop$region <- "Crop"
top_mid$region  <- "Midgut"
top_hind$region <- "Hindgut"

top_all <- rbind(top_crop, top_mid, top_hind)

# reorder columns nicely
top_all <- top_all[, c("lipid", "lipid_class", "region", "driver_strength")]

# sort by region then strength
top_all <- top_all[order(top_all$region, -top_all$driver_strength), ]

# view table
View(top_all)
