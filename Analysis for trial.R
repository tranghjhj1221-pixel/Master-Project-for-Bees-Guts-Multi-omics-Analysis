library(readr)
library(dplyr)
library(stringr)

# 1) read the TSV file
pg <- read_tsv("report.pg_matrix.tsv")

# 2) inspect column names
colnames(pg)

# 3) rename the sample columns to something simple
colnames(pg) <- c(
  "Protein.Group",
  "Protein.Names",
  "Genes",
  "First.Protein.Description",
  "Blank",
  "Crop",
  "Hind",
  "Mid"
)

# 4) check the first few rows
head(pg)

#Creating metadata
metadata <- data.frame(
  sample = c("Blank", "Crop", "Hind", "Mid"),
  condition = c("Blank", "Crop", "Hind", "Mid")
)

metadata
#Check for contaiminants: 
sum(grepl("Cont_", pg$Protein.Group))

#Removing contaminants + keep useful columns
library(dplyr)

# remove contaminants
pg_clean <- pg %>%
  filter(!grepl("Cont_", Protein.Group))

# keep only what I need
pg_clean <- pg_clean %>%
  select(
    Protein.Group,
    Protein.Names,
    First.Protein.Description,
    Genes,
    Blank, Crop, Hind, Mid
  )
head(pg_clean)
dim(pg_clean)

#Sanity check for contaiminants: 
sum(grepl("Cont_", pg_clean$Protein.Group))

#Check for trypsin in the blank
pg %>%
  filter(
    grepl("trypsin", Protein.Names, ignore.case = TRUE) |
      grepl("trypsin", First.Protein.Description, ignore.case = TRUE) |
      grepl("trypsin", Genes, ignore.case = TRUE)
  ) %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description, Blank, Crop, Hind, Mid)

##Presence/absence testing (not gonna impute anything)
pa <- pg_clean %>%
  mutate(
    Blank_present = !is.na(Blank),
    Crop_present  = !is.na(Crop),
    Hind_present  = !is.na(Hind),
    Mid_present   = !is.na(Mid)
  )
#remove blank-associated 

view(pg_clean) #inspect everything manually

blank_hits <- pg_clean %>%
  filter(!is.na(Blank)) %>%
  select(
    Protein.Group,
    Protein.Names,
    Genes,
    Blank, Crop, Hind, Mid
  ) %>%
  arrange(desc(Blank))  # highest blank signal first

# View in spreadsheet format
View(blank_hits)
print(blank_hits$Protein.Names)
# OR print top 50 in console
head(blank_hits, 50)

##CREATE a classification for everything I see in blank, basically flagging everything that is suspicious in blank

pg_classified <- pg_clean %>%
  mutate(
    in_blank = !is.na(Blank),
    
    max_sample = pmax(Crop, Hind, Mid, na.rm = TRUE),
    
    sample_blank_ratio = max_sample / Blank,
    sample_blank_ratio = ifelse(
      is.infinite(sample_blank_ratio) | is.nan(sample_blank_ratio),
      NA,
      sample_blank_ratio
    ),
    
    class = case_when(
      !in_blank ~ "High_confidence_bee", #most trustworthy,not detected in blank
      in_blank & (is.na(sample_blank_ratio) | sample_blank_ratio < 2) ~ "Bee_in_blank_suspicious", #sample and blank at similar levels, low confidence
      in_blank & sample_blank_ratio >= 2 ~ "Bee_possible_carryover" #detected in blank but less prominent (sample signal is 3-10x times higher)
    )
  )

view(pg_classified)
view(pg_classified |>
       filter(class == "Bee_in_blank_suspicious"))

##I have decided to press on with Bee_possible carried over and High confidence bee
pg_final <- pg_classified %>%
  filter(class %in% c("High_confidence_bee", "Bee_possible_carryover"))

# Presence / absence
pg_pa <- pg_final %>%
  mutate(
    Crop_present = !is.na(Crop),
    Hind_present = !is.na(Hind),
    Mid_present  = !is.na(Mid)
  )

# Region-specific proteins
crop_only <- pg_pa %>%
  filter(Crop_present & !Hind_present & !Mid_present)

hind_only <- pg_pa %>%
  filter(Hind_present & !Crop_present & !Mid_present)

mid_only <- pg_pa %>%
  filter(Mid_present & !Crop_present & !Hind_present)

# Shared proteins (all 3 regions)
shared_all <- pg_pa %>%
  filter(Crop_present & Hind_present & Mid_present)

#Counts (sanity check)
cat("Crop only:", nrow(crop_only), "\n")
cat("Hind only:", nrow(hind_only), "\n")
cat("Mid only:", nrow(mid_only), "\n")
cat("Shared (all):", nrow(shared_all), "\n")

#Top proteins per region
top_crop <- crop_only %>%
  arrange(desc(Crop)) %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description, Crop, Hind, Mid, class) %>%
  head(20)

top_hind <- hind_only %>%
  arrange(desc(Hind)) %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description, Crop, Hind, Mid, class) %>%
  head(20)

top_mid <- mid_only %>%
  arrange(desc(Mid)) %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description, Crop, Hind, Mid, class) %>%
  head(20)

View(top_crop)
View(top_hind)
View(top_mid)


###Fold-change analysis, pair-wise run (kinda use crop as a baseline here) : 
library(dplyr)

# Compute log2 fold-changes
pg_fc <- pg_final %>%
  mutate(
    log2_Crop = log2(Crop + 1),
    log2_Hind = log2(Hind + 1),
    log2_Mid  = log2(Mid + 1),
    
    Mid_vs_Crop  = log2_Mid  - log2_Crop,
    Hind_vs_Crop = log2_Hind - log2_Crop,
    Hind_vs_Mid  = log2_Hind - log2_Mid
  )

#Top enriched in Mid vs Crop
top_mid_vs_crop <- pg_fc %>%
  arrange(desc(Mid_vs_Crop)) %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         Crop, Hind, Mid, Mid_vs_Crop) %>%
  head(20)

#Top enriched in Hind vs Crop
top_hind_vs_crop <- pg_fc %>%
  arrange(desc(Hind_vs_Crop)) %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         Crop, Hind, Mid, Hind_vs_Crop) %>%
  head(20)

# Top enriched in Hind vs Mid
top_hind_vs_mid <- pg_fc %>%
  arrange(desc(Hind_vs_Mid)) %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         Crop, Hind, Mid, Hind_vs_Mid) %>%
  head(20)


View(top_mid_vs_crop) #damn 400x fold

View(top_hind_vs_crop)
View(top_hind_vs_mid)


#Quick look at them show me dis:
#midgut is doing nutrient transport, detoxification (P450), active metabolism, epithelial structure protein
#hindgut = transport + detox + immune/microbiome interface (chitinase is so interesting)

##Comparision between multi class:
library(dplyr)

pg_multi <- pg_final %>%
  mutate(
    log2_Crop = log2(Crop + 1),
    log2_Hind = log2(Hind + 1),
    log2_Mid  = log2(Mid + 1),
    
    # differences
    Mid_vs_Crop  = log2_Mid  - log2_Crop,
    Mid_vs_Hind  = log2_Mid  - log2_Hind,
    
    Hind_vs_Crop = log2_Hind - log2_Crop,
    Hind_vs_Mid  = log2_Hind - log2_Mid,
    
    Crop_vs_Mid  = log2_Crop - log2_Mid,
    Crop_vs_Hind = log2_Crop - log2_Hind,
    
    # classification
    region_class = case_when(
      Mid_vs_Crop > 1 & Mid_vs_Hind > 1 ~ "Mid_specific",
      Hind_vs_Crop > 1 & Hind_vs_Mid > 1 ~ "Hind_specific",
      Crop_vs_Mid > 1 & Crop_vs_Hind > 1 ~ "Crop_specific",
      TRUE ~ "Shared_or_mixed"
    )
  )

# Check out for skewness?
library(dplyr)
library(reshape2)
library(ggplot2)

# 1) subset relevant columns
pg_subset <- pg_final %>%
  select(Genes, Crop, Hind, Mid)

# 2) melt (pivot longer)
pg_long <- melt(
  pg_subset,
  id.vars = "Genes",
  variable.name = "Sample",
  value.name = "Intensity"
)

# 3) remove NA
pg_long <- pg_long %>%
  filter(!is.na(Intensity))

# 4) log2 transform
pg_long <- pg_long %>%
  mutate(log2_Intensity = log2(Intensity + 1))

# 5) density plot
ggplot(pg_long, aes(x = log2_Intensity, fill = Sample)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Protein intensity distribution (log2)",
    x = "log2 Intensity",
    y = "Density"
  )
#Log2 Histogram looks great, cannot do qqplot here

###Heatmap Inspection
library(dplyr)
install.packages("pheatmap")
library(pheatmap)
library(tibble)

#creating matrix

mat <- pg_final %>%
  mutate(
    row_id = ifelse(is.na(Genes) | Genes == "", Protein.Group, Genes),
    row_id = make.unique(row_id)
  ) %>%
  select(row_id, Crop, Hind, Mid) %>%
  column_to_rownames("row_id") %>%
  as.matrix()


# 2) log2 transform
mat_log <- log2(mat + 1)

# 3) select top variable proteins
top_var <- apply(mat_log, 1, var, na.rm = TRUE) %>%
  sort(decreasing = TRUE) %>%
  head(100)

view(top_var)

mat_top <- mat_log[names(top_var), ]

# 4) simple imputation for plotting (low value)
mat_top[is.na(mat_top)] <- min(mat_top, na.rm = TRUE) - 1 


# 5) row scaling (important!)
mat_scaled <- t(scale(t(mat_top)))
#here, I'm comparing Crop, Mid, and Hind protein intensity for that protein, there is no benchmarking done yet
#not sure if I'll trust the z-score here, this is a rough estimate
#Let's try both heatmap


pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  main = "Top variable proteins"
)


pheatmap(
  mat_top,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  main = "Top variable proteins"
)

##Did not do a differential expression analysis here, will do with real samples

###PIck out top 20 proteins for GO Enrichment
library(dplyr)

# Start from proteins not detected in blank
pg_noblank <- pg_classified %>%
  filter(class == "High_confidence_bee") %>%
  mutate(
    log2_Crop = log2(Crop + 1),
    log2_Hind = log2(Hind + 1),
    log2_Mid  = log2(Mid + 1),
    
    Crop_vs_Mid  = log2_Crop - log2_Mid,
    Crop_vs_Hind = log2_Crop - log2_Hind,
    
    Mid_vs_Crop  = log2_Mid - log2_Crop,
    Mid_vs_Hind  = log2_Mid - log2_Hind,
    
    Hind_vs_Crop = log2_Hind - log2_Crop,
    Hind_vs_Mid  = log2_Hind - log2_Mid
  )

# Top 20 Crop-enriched proteins
top20_crop_enriched <- pg_noblank %>%
  filter(!is.na(Crop)) %>%
  mutate(crop_score = Crop_vs_Mid + Crop_vs_Hind) %>%
  arrange(desc(crop_score)) %>%
  select(
    Protein.Group, Protein.Names, Genes, First.Protein.Description,
    Crop, Hind, Mid, Crop_vs_Mid, Crop_vs_Hind, crop_score
  ) %>%
  head(20)

# Top 20 Mid-enriched proteins
top20_mid_enriched <- pg_noblank %>%
  filter(!is.na(Mid)) %>%
  mutate(mid_score = Mid_vs_Crop + Mid_vs_Hind) %>%
  arrange(desc(mid_score)) %>%
  select(
    Protein.Group, Protein.Names, Genes, First.Protein.Description,
    Crop, Hind, Mid, Mid_vs_Crop, Mid_vs_Hind, mid_score
  ) %>%
  head(20)

# Top 20 Hind-enriched proteins
top20_hind_enriched <- pg_noblank %>%
  filter(!is.na(Hind)) %>%
  mutate(hind_score = Hind_vs_Crop + Hind_vs_Mid) %>%
  arrange(desc(hind_score)) %>%
  select(
    Protein.Group, Protein.Names, Genes, First.Protein.Description,
    Crop, Hind, Mid, Hind_vs_Crop, Hind_vs_Mid, hind_score
  ) %>%
  head(20)

# View them
View(top20_crop_enriched)
View(top20_mid_enriched)
View(top20_hind_enriched)
##Save them as files
write.csv(top20_crop_enriched, "top20_crop_enriched.csv", row.names = FALSE)
write.csv(top20_mid_enriched,  "top20_mid_enriched.csv", row.names = FALSE)
write.csv(top20_hind_enriched, "top20_hind_enriched.csv", row.names = FALSE)



##LIMMA  has not been used here therefore the pvalues for GO term enrichment was FAKED by using the rank-based score 
  # I FAKE PVALUES HERE TO TEST OUT THE WORKFLOW
#MID
mid_erminej <- top20_mid_enriched %>%
  mutate(
    rank = row_number(desc(mid_score)),
    fake_p = rank^2 / (n() + 1)^2 ##basicially taking the rank devided by 21 (out of my top20 table) and square it so that it looks like pvalues
  ) %>%
  select(Genes, fake_p)

write.table(
  mid_erminej,
  file = "mid_erminej_input.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
#HIND
hind_erminej <- top20_hind_enriched %>%
  mutate(
    rank = row_number(desc(hind_score)),
    fake_p = rank^2 / (n() + 1)^2 
  ) %>%
  select(Genes, fake_p)

write.table(
  hind_erminej,
  file = "hind_erminej_input.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#CROP
crop_erminej <- top20_crop_enriched %>%
  mutate(
    rank = row_number(desc(crop_score)),
    fake_p = rank^2 / (n() + 1)^2
  ) %>%
  select(Genes, fake_p)

write.table(
  crop_erminej,
  file = "crop_erminej_input.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

View(crop_erminej)
head(pg_final$Genes)


## eDIT YOUR 
library(dplyr)
library(readr)
library(stringr)

anno <- read_tsv(file.choose())

erminej_anno <- anno %>%
  transmute(
    probe_id = Entry,
    gene_symbol = Entry,
    gene_name = Entry,
    go_ids = str_replace_all(`Gene Ontology IDs`, ";\\s*", "|")
  ) %>%
  filter(!is.na(go_ids), go_ids != "")

write.table(
  erminej_anno,
  file = "erminej_annotation_final.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
sum(erminej_anno$probe_id %in% pg_final$Protein.Group) #This should be >>0 otherwise it wouldnt run on ErmineJ
View(erminej_anno)
