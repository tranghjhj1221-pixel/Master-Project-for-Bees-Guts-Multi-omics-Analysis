library(tibble)
library(tidyverse)
library(readxl)
install.packages("reshape2")
library(reshape2)

# load file
dat <- read.csv(file.choose(), check.names = FALSE)

# remove rows already annotated as lipids
dat <- dat %>%
  filter(is.na(`preferred_annotation:annotation_method`) |
           `preferred_annotation:annotation_method` != "Lipid Annotation")

# clean names by removing trailing match scores like "(0.882)"
dat <- dat %>%
  mutate(
    name_clean = str_remove(
      `spectral_db_matches:spectral_db_matches`,
      "\\s*\\([0-9.]+\\)$"
    )
  )

# replace empty or missing names with mz/rt-based placeholders
dat <- dat %>%
  mutate(
    name_clean = ifelse(
      is.na(name_clean) | name_clean == "",
      paste0("unknown_mz", round(mz, 4), "_rt", round(rt, 2)),
      name_clean
    )
  )

# make names unique
dat$name_clean <- make.unique(dat$name_clean) #add surfix .1 .2 to similar names

# move cleaned names to row names
rownames(dat) <- NULL
library(tibble)

dat <- rownames_to_column(dat, var = "name_clean")

write.csv(dat,
          file = "clean_up_metabolites_run_from_R.csv",
          row.names = FALSE)
# keep only peak area columns
dat.area <- dat %>%
  select(contains("area"))

# log2 transform area values
# +1 avoids log2(0) = -Inf, but NA values remain NA <- NA is biologically meaningful 
dat.area.log <- log2(dat.area + 1) #log transform

# inspect histogram
dat.melt <- melt(as.matrix(dat.area.log))
hist(dat.melt$value[!is.na(dat.melt$value)], breaks = 100) #histogram looks good

# no pooled QC samples were included, so QC-CV filtering is skipped
# no detection-frequency filtering applied because each sample represents a different tissue
# and tissue-specific metabolites may only appear in one sample

dat.area.filt <- dat.area

# log-transform filtered area matrix
dat.area.log <- log2(dat.area.filt + 1)

View(dat.area.log)



library(pheatmap)
library(viridis)

# keep only the 3 tissue columns and rename them
dat.hm <- dat.area.log[, c("datafile:Crop_P1-A-1_562.d:area",
                           "datafile:Mid_P1-A-2_563.d:area",
                           "datafile:Hind_P1-A-3_564.d:area")]

colnames(dat.hm) <- c("Crop", "Midgut", "Hindgut")

# keep only annotated compounds
dat.annot <- dat.hm[!grepl("^unknown_", rownames(dat.hm)), ]

# manually scale each row while ignoring NA
scale_row_na <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  
  # if variance is zero or cannot be calculated, return 0s
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  
  (x - m) / s
}

dat.annot.scaled <- t(apply(dat.annot, 1, scale_row_na))
colnames(dat.annot.scaled) <- colnames(dat.annot)
rownames(dat.annot.scaled) <- rownames(dat.annot)

# later: keep only the most variable annotated compounds for a cleaner heatmap
rv <- apply(dat.annot.scaled, 1, var, na.rm = TRUE)
top_n <- min(50, nrow(dat.annot.scaled))
dat.annot.scaled <- dat.annot.scaled[order(rv, decreasing = TRUE)[1:top_n], ]

# fixed color scale
breaks_list <- seq(-2, 2, by = 0.1)

# draw heatmap
pheatmap(
  dat.annot.scaled,
  scale = "none", #already manually scal
  color = viridis(length(breaks_list) - 1, direction = -1),
  breaks = breaks_list,
  border_color = NA,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  na_col = "grey90", #keep missing values visible in grey color
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 9,
  fontsize_row = 7,
  fontsize_col = 10,
  main = "Annotated metabolites"
)


#eugh
#crop has dietary sugar, not filtered out
#midgut has amino acid, nuclotides, and NAD, and fiboflavin, digestion, absorption and metabolic processing
#seems like mid and hind share metabolic functioin

#Lets try a heat map with the lipid that got through in this metabolites samples
library(tibble)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(viridis)

# load file
dat <- read.csv(file.choose(), check.names = FALSE)

# keep only rows annotated as lipids
dat.lipid <- dat %>%
  filter(`preferred_annotation:annotation_method` == "Lipid Annotation")

# make a clean lipid name column
# first try spectral DB name, then preferred annotation compound name,
# then fall back to mz/rt if needed
dat.lipid <- dat.lipid %>%
  mutate(
    name_clean = `spectral_db_matches:spectral_db_matches`,
    name_clean = str_remove(name_clean, "\\s*\\([0-9.]+\\)$"),
    name_clean = ifelse(
      is.na(name_clean) | name_clean == "",
      `preferred_annotation:compound_name`,
      name_clean
    ),
    name_clean = ifelse(
      is.na(name_clean) | name_clean == "",
      paste0("unknown_lipid_mz", round(mz, 4), "_rt", round(rt, 2)),
      name_clean
    )
  )

# standardize a few common lipid names
dat.lipid <- dat.lipid %>%
  mutate(
    name_clean = case_when(
      name_clean %in% c("[1-MYRISTOYL-GLYCEROL-3-YL]PHOSPHONYLCHOLINE", "LPC 14:0") ~ "LPC(14:0)",
      name_clean %in% c("1-Oleoylglycerophosphocholine", "LPC 18:1") ~ "LPC(18:1)",
      name_clean %in% c("1-Palmitoylglycerophosphocholine", "LPC 16:0") ~ "LPC(16:0)",
      TRUE ~ name_clean
    )
  )

# make names unique
dat.lipid$name_clean <- make.unique(dat.lipid$name_clean)

# move names to rownames
rownames(dat.lipid) <- NULL
dat.lipid <- column_to_rownames(dat.lipid, var = "name_clean")

View(dat.lipid)

# keep only the 3 tissue area columns
dat.lipid.area <- dat.lipid[, c("datafile:Crop_P1-A-1_562.d:area",
                                "datafile:Mid_P1-A-2_563.d:area",
                                "datafile:Hind_P1-A-3_564.d:area")]

colnames(dat.lipid.area) <- c("Crop", "Midgut", "Hindgut")

# log transform
dat.lipid.log <- log2(dat.lipid.area + 1)

View(dat.lipid.log)

# manually scale each row while ignoring NA
scale_row_na <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

dat.lipid.scaled <- t(apply(dat.lipid.log, 1, scale_row_na))
colnames(dat.lipid.scaled) <- colnames(dat.lipid.log)
rownames(dat.lipid.scaled) <- rownames(dat.lipid.log)

# optional: keep top variable lipids for a cleaner heatmap
rv <- apply(dat.lipid.scaled, 1, var, na.rm = TRUE)
top_n <- min(50, nrow(dat.lipid.scaled))
dat.lipid.scaled <- dat.lipid.scaled[order(rv, decreasing = TRUE)[1:top_n], ]

# fixed color scale
breaks_list <- seq(-2, 2, by = 0.1)

# lipid heatmap
pheatmap(
  dat.lipid.scaled,
  scale = "none",
  color = viridis(length(breaks_list) - 1, direction = -1),
  breaks = breaks_list,
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  na_col = "grey90",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 9,
  fontsize_row = 7,
  fontsize_col = 10,
  main = "Lipid-annotated features across gut tissues"
)
#heatmap is expected because crop is showing less lipid, that makes sense

