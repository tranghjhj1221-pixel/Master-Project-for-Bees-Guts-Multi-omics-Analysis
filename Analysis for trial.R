library(readr)
library(dplyr)
library(stringr)

# 1) read the TSV file
pg <- read_tsv("C:/Trang/outputtrial2/report.pg_matrix.tsv")

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