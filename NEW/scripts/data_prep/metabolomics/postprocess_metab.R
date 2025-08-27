## Postprocess Metabolomics Data: Building Analytical Dataset


library(tidyverse) ; library(data.table)
lapply(list.files("../../pantry/functions/", recursive = T, full.names = T), source)


######################################################### 
## Load & build metabolomics dataset
######################################################### 

metabs <- fread("../data/processed/no_batch_adj/plasma/merged_data/ln_merged_QCd_plasma_knowns.csv") %>%
  separate(sample_id, sep="_", into=c("id", "time")) %>%
  filter(!startsWith(id, "PREF"))  # 740 metabolites

analysis <- readRDS("../data/processed/premier_analysis.rda") %>%
  filter(!is.na(genotype))

# merge in premier genotypes & covariates
analysis2 <- left_join(
  analysis %>% select(id, age, sex, bmi, PC1z, PC2z, PC3z, genotype, meal_choice, genetic_cat_x_meal_choice), 
  metabs, by="id") %>%
  mutate_at("time", ~factor(ifelse(.=="BL", 0, .), levels=c(0,120,235,360)))

# load metabolomics info (built by PH)
met_info <- fread("../data/processed/met_info.csv") %>%
  mutate(Compound = gsub("_.*", "", Compound_Id), .after="Compound_Id")
metabolites <- names(analysis2 %>% select(starts_with("HMDB")))
metabclass <- fread("../data/raw/metabolomics/metabolites_subclasses.csv") %>%
  rename(HMDB=HMDBID)

## load external metabolomics data dictionary: RefMet
#devtools::install_github("metabolomicsworkbench/RefMet")
library(RefMet)
metab_dict <- refmet_map_df(gsub("_.*", "", metabolites))


###############################################################################
## 1) Check for duplicated metabolites (across data types) & select ONE with  
## stronger postprandial changes
###############################################################################

## Run LME for all metabolites on glucose 
lme_metab_mmtt_timeonly <- lapply(1:length(metabolites), function(m) {
  metab <- metabolites[m]
  Metabolite <- met_info$Name[met_info$HMDB_Id == metab]
  run_lme(exposure = "time", outcome=metab, outcome_label=Metabolite,
          covariates = "time", #coefficients_to_print = coefs_to_print, 
          data=analysis2 %>% filter(time %in% c(0,120,235)), digits=c(1,3)) %>%
    mutate(outcome=metab) 
}) %>% do.call(rbind.data.frame, .) ; lme_metab_mmtt_timeonly

#lme_metab_mmtt_timeonly %>% fwrite("../data/processed/tab_result_lme_metab_mmtt_timeonly.csv")
#lme_metab_mmtt_timeonly <- fread("../data/processed/tab_result_lme_metab_mmtt_timeonly.csv")

# ================================================================
## Select metabolites with duplicates across measurement types 
# ================================================================

# Count number of total metabolites 
length(gsub("_.*", "", metabolites)) # 740 metabolites

# Count number of duplicated metabolites (across methods: cn/cp/hn/hp)
metabolites_dup <- unique(gsub("_.*", "", metabolites[which(duplicated(gsub("_.*", "", metabolites)))])) ; length(metabolites_dup) # 62

# Among duplicates, which metabolite_method to keep (based on stronger postprandial changes/time)
metabolites_dup.keep <- lapply(metabolites_dup, function(metab) {
  (lme_metab_mmtt_timeonly %>% 
     filter(Exposure == "Time_(joint)") %>%
     filter(startsWith(outcome, metab)) %>%
     arrange(anovaP) %>% pull(outcome))[1]
}) %>% do.call(rbind, .)

## Select the *top* metabolites for duplicates (strongest ANOVA P-value for change/time ^^)
metabolites_unq <- analysis2 %>% select(-(starts_with(metabolites_dup) & !starts_with(metabolites_dup.keep))) %>% 
  select(starts_with("HMDB")) %>% names()

# Drop source flag from metabolite name
analysis2 <- analysis2 %>% select(!starts_with("HMDB") | starts_with(metabolites_unq)) %>%
  rename_with(., ~gsub("_.*", "", .), starts_with("HMDB"))

# Update metabolomics info to exclude duplicates
met_info <- met_info %>% filter(HMDB_Id %in% metabolites_unq) %>%
  mutate(HMDB = gsub("_.*", "", HMDB_Id))

metabolites <- analysis2 %>% select(starts_with("HMDB")) %>% names() ; length(metabolites) # 677

## Filter metabolite x time data to "kept" duplicates
lme_metab_mmtt_timeonly <- lme_metab_mmtt_timeonly %>% 
  filter(outcome %in% metabolites_unq) %>%
  mutate_at("outcome", ~gsub("_.*", "", .)) ; dim(lme_metab_mmtt_timeonly) #2031 15

# Select metabolites changing postprandially
#metabolites_pp <- lme_metab_mmtt_timeonly %>% filter(anovaP<0.05) %>% pull(outcome)
metabolites_pp <- lme_metab_mmtt_timeonly %>% filter(p<0.05) %>% pull(outcome) %>% unique() #396
length(metabolites_pp) #N=401

# Write .csv 
lme_metab_mmtt_timeonly %>% 
  fwrite("../data/processed/tab_prelim_lme_metab_mmtt_timeonly.csv")

## Save list of metabolites -------------
as.data.frame(metabolites) %>% fwrite("../data/processed/metbolites.txt", row.names = F)

## Metabolomics IDs
metabolomics_ids <- analysis2 %>% filter(!is.na(HMDB0034408)) %>% pull(id) %>% unique()

## Create indicator for who has metabolomics data
analysis <- analysis %>% mutate(has_metabolomics = ifelse(id %in% metabolomics_ids,1,0))
table(analysis$has_metabolomics) # 0=5; 1=17

analysis2 <- analysis2 %>% mutate(has_metabolomics = ifelse(id %in% metabolomics_ids,1,0))
saveRDS(analysis2, "../data/processed/premier_analysis2.rda")


