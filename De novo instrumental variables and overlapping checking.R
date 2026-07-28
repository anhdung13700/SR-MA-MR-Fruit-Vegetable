library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(readr)
library(biomaRt)
library(phenoscanner)

Sys.setenv(OPENGWAS_JWT = "PUT YOUR OWN OPENGWAS KEYS")

exposure_ids <- c("ukb-b-16576", "ukb-b-3881", "ukb-b-8089", "ukb-b-1996", "ieu-b-4877","ukb-b-5779")
print("Extracting instruments from IEU OpenGWAS...")

exposure_dat <- extract_instruments(
  outcomes = exposure_ids,
  p1 = 5e-08,
  clump = TRUE,
  r2 = 0.001,
  kb = 10000
)
# Connect to Ensembl SNP database
snp_mart <- useEnsembl(
  biomart = "snp",
  dataset = "hsapiens_snp",
  mirror = "useast"   # or "asia"
)

# Retrieve gene annotation for SNPs
snp_info <- getBM(
  attributes = c(
    "refsnp_id",
    "chr_name",
    "chrom_start"
  ),
  filters = "snp_filter",
  values = exposure_dat$SNP,
  mart = snp_mart
)
gene_mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)
gene_annotation <- lapply(1:nrow(snp_info), function(i) {
  chr <- snp_info$chr_name[i]
  pos <- snp_info$chrom_start[i]
  
genes <- getBM(
    attributes = c("external_gene_name"),
    filters = c("chromosome_name", "start", "end"),
    values = list(chr, pos - 10000, pos + 10000),
    mart = gene_mart
  )
  
data.frame(
    SNP = snp_info$refsnp_id[i],
    Gene = paste(unique(genes$external_gene_name), collapse = ", ")

# Clean annotation table
gene_annotation <- lapply(1:nrow(snp_info), function(i) {
    chr <- snp_info$chr_name[i]
    pos <- snp_info$chrom_start[i]
    
genes <- getBM(
      attributes = c("external_gene_name"),
      filters = c("chromosome_name", "start", "end"),
      values = list(chr, pos - 10000, pos + 10000),
      mart = gene_mart
    )
    
    data.frame(
      SNP = snp_info$refsnp_id[i],
      Gene = paste(unique(genes$external_gene_name), collapse = ", ")
    )
  }) %>% bind_rows()

#Clean and Calculate R2 and F-statistics

#Map friendly names to the IDs
  id_map <- data.frame(
    id.exposure = c("ukb-b-16576", "ukb-b-3881", "ukb-b-8089", "ukb-b-1996","ieu-b-4877","ukb-b-5779"),
    Exposure_Name = c("Dried Fruit", "Fresh Fruit", "Cooked Vegetable", "Salad/Raw Vegetable","Smoking Initiation","Alcohol Intake Frequency"))
  
clean_snp_list <- exposure_dat %>%
    left_join(id_map, by = "id.exposure") %>%
    left_join(gene_annotation, by = "SNP", relationship = "many-to-many") %>%
    mutate(
      R2 = (beta.exposure^2) / (beta.exposure^2 + samplesize.exposure * se.exposure^2),
      F_statistic = (beta.exposure^2) / (se.exposure^2)
    ) %>%
    dplyr::select(   # <-- ADDED dplyr:: HERE TO PREVENT MASKING
      Exposure_Name,
      SNP,
      Gene,
      CHR = chr.exposure,
      POS = pos.exposure,
      Effect_Allele = effect_allele.exposure,
      Other_Allele = other_allele.exposure,
      EAF = eaf.exposure,
      Beta = beta.exposure,
      SE = se.exposure,
      P_value = pval.exposure,
      Samplesize = samplesize.exposure,
      R2,
      F_statistic
    ) %>%
    arrange(Exposure_Name, CHR, POS)

#Export to CSV
write_csv(clean_snp_list, "Supplementary_Table_Genetic_Instruments.csv")

#Summary Table for manuscript
summary_stats <- clean_snp_list %>%
  group_by(Exposure_Name) %>%
  summarise(
    Num_SNPs = n(),
    Mean_F_stat = mean(F_statistic),
    Total_R2 = sum(R2)
  )

print(summary_stats)

#Check overlapped between exposures
#dried vs fresh
#dried vs cooked
#dried vs salad
#fresh vs cooked
#fresh vs salad
#cooked vs salad
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(readr) # For exporting to CSV

# ---------------------------------------------------------
# 1. Define your traits and IDs
# ---------------------------------------------------------
traits <- list(
  "Fresh fruit" = "ukb-b-3881",
  "Dried fruit" = "ukb-b-16576",
  "Salad"       = "ukb-b-1996",
  "Cooked veg"  = "ukb-b-8089"
)

# ---------------------------------------------------------
# 2. Define the 6 specific pairs to compare
# ---------------------------------------------------------
comparisons <- list(
  c("Fresh fruit", "Salad"),
  c("Cooked veg", "Salad"),
  c("Cooked veg", "Dried fruit"),
  c("Cooked veg", "Fresh fruit"),
  c("Salad", "Cooked veg"),
  c("Salad", "Dried fruit"),
  c("Salad", "Fresh fruit")
)

# ---------------------------------------------------------
# 3. Loop through each pair, export immediately, clear memory
# ---------------------------------------------------------
for (pair in comparisons) {
  
  trait1_name <- pair[1]
  trait2_name <- pair[2]
  
  id1 <- traits[[trait1_name]]
  id2 <- traits[[trait2_name]]
  
  message(sprintf("Processing: %s vs %s...", trait1_name, trait2_name))
  
  # Step A: Extract significant SNPs for Trait 1
  exposure_dat <- extract_instruments(outcomes = id1, p1 = 5e-8,clump = TRUE,r2 = 0.001,kb = 10000)
  
  if (!is.null(exposure_dat) && nrow(exposure_dat) > 0) {
    
    # Step B: Look up those exact SNPs in Trait 2
    overlap_dat <- associations(variants = exposure_dat$SNP, id = id2)
    
    if (!is.null(overlap_dat) && nrow(overlap_dat) > 0) {
      
      # Step C: Add clear labels
      overlap_dat <- overlap_dat %>%
        mutate(
          Exposure_Trait = trait1_name,
          Exposure_ID = id1,
          Outcome_Trait = trait2_name,
          Outcome_ID = id2
        ) %>%
        relocate(Exposure_Trait, Outcome_Trait, Exposure_ID, Outcome_ID, .before = 1) 
      
      # Step D: Create a safe file name (replace spaces with underscores)
      safe_t1 <- gsub(" ", "_", trait1_name)
      safe_t2 <- gsub(" ", "_", trait2_name)
      file_name <- sprintf("overlap_%s_vs_%s.csv", safe_t1, safe_t2)
      
      # Step E: Export THIS specific result immediately
      write_csv(overlap_dat, file_name)
      message(sprintf(" -> Success! Saved to %s\n", file_name))
      
      # Clear the overlap data from memory
      rm(overlap_dat)
      
    } else {
      message(" -> No overlapping SNPs found in Outcome database.\n")
    }
  } else {
    message(" -> No significant SNPs found for Exposure.\n")
  }
  
  # Step F: Housekeeping
  # Clear the exposure data from memory
  if(exists("exposure_dat")) rm(exposure_dat)
  
  # Force R to clean up RAM
  gc() 
  
  # Pause for 3 seconds to let the API server breathe and prevent timeouts
  Sys.sleep(3) 
}

message("All comparisons finished safely!")