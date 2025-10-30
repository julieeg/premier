## Analyze metabolomics data

library(tidyverse) ; library(data.table)
lapply(list.files("../../pantry/functions/", recursive = T, full.names = T), source)

library(poolr)
library(ComplexHeatmap)


################################################################################ 
## Load & build metabolomics dataset 
################################################################################

# analysis dataset with metabolomics 
analysis2 <- readRDS("../data/processed/premier_analysis2.rda") %>%
  filter(!is.na(genotype))

metabolomics_ids <- analysis2 %>% filter(has_metabolomics==1) %>% pull(id) %>% unique()

# analysis dataset with phenotypes & glycemic traits
analysis <- readRDS("../data/processed/premier_analysis.rda") %>%
  filter(!is.na(genotype)) %>%
  mutate(has_metabolomics = ifelse(id %in% metabolomics_ids,1,0)) %>%
  mutate(tg_0=tg, tg_log_0=tg_log, hdl_0=hdl, ldl_0=ldl)


# load metabolomics info (built by PH)
met_info <- fread("../data/processed/met_info_processed.csv")
metabolites <- names(analysis2 %>% select(starts_with("HMDB"))) #817 TOTAL metabolites
metabolites_pp <- fread("../data/processed/metbolites_pp.txt")$metabolites_pp # 547 changing postprandially (p<0.05)
metabclass <- fread("../data/processed/metabolites_subclasses.csv") %>%
  rename(HMDB=HMDBID)

## load external metabolomics data dictionary: RefMet
#devtools::install_github("metabolomicsworkbench/RefMet")
#library(RefMet)
#metab_dict <- refmet_map_df(gsub("_.*", "", metabolites))
#metab_dict %>% fwrite("../data/processed/met_dictionary.csv")
metab_dict <- fread("../data/processed/met_dictionary.csv")


###############################################################################
## Descriptive statistics for metabolites
###############################################################################

## Lipid levels at each time point by genotype & meal ---------------

clinvars <- c("glucose", "insulin", "tg", "ldl", "hdl")
clinvars.lab <- c("Glucose", "Insulin", "TG", "LDL", "HDL")
names(clinvars.lab) <- clinvars

## Compare characteristics BY having metabolomics data
lapply(summary_table_vars.l, function(vars) {
  print_summary_table(
    analysis, vars_to_summarize = vars$vars, p_adjust = "none", digits = c(1,1,3),
    var_strata = "has_metabolomics", #var_strata_order = c("Yes"=1, "No"=0),
    p_types = "descriptive", p_smalln=T) %>%
    mutate_all(., ~gsub("[%]", "", .), .) }) %>% 
  do.call(rbind.data.frame, .) %>%
  rename("Has metabolomics"="1", "No metabolomics"="0") %>% 
  fwrite(paste0("../output/tab_descr_byMetabolomics.csv"), row.names = T)

table(analysis$genotype, analysis$has_metabolomics)
prop.table(table(analysis$genotype, analysis$has_metabolomics))
fisher.test(analysis$has_metabolomics, analysis$genotype)


# ==============================================================================
## Determine Number of Effective Tests (number of NON-correlated metabolites)  
## for multiple comparison correction = 30
# ==============================================================================

library(poolr)

# Known metabolites = 817
metabolites.df <- analysis2 %>% filter(time == 0) %>% select(all_of(metabolites_pp)) %>%
  filter(complete.cases(.)) # Remove NAs ; dim(metabolites.df): 22 x 547

# Effective tests = 30
corr_metab <- cor(metabolites.df) ## create correlation matrix
effective_tests <- meff(corr_metab, method="liji") ## 30 un-correlated metabolites 
effective_tests # 30


################################################################################
## Primary analysis: changes in metabolomics profile after standrdized MMTT
################################################################################

# =========================================================
## 2) Run LME for metab ~ genotype x time (categorical)
## Restricting to postprandial metabolites
# =========================================================

lme_metab_mmtt_geno_all <- lapply(1:length(metabolites), function(m) {
  metab <- metabolites[m]
  Metabolite <- met_info$Name[met_info$HMDB == metab][1] 
  run_lme(exposure = "genotype", outcome=metab, 
          covariates = "genotype*time", coefficients_to_print = c(Genotype="genotype", Time="time"), 
          data=analysis2 %>% filter(time %in% c(0,120,235)), digits=c(1,3)) %>%
    mutate(outcome=metab, Outcome=Metabolite) %>%
    select(Exposure, outcome, Outcome, beta, se, p, anovaP, Effect, lowCI, upCI) %>%
    filter(grepl(" x ", Exposure)) %>%
    mutate_at("Exposure", ~gsub(".* x", "HC Genotype x", .))
  }) %>% do.call(rbind.data.frame, .)

lme_metab_mmtt_geno_pp <- lme_metab_mmtt_geno %>% filter(outcome %in% metabolites_pp)

## Save as .csv
#lme_metab_mmtt_geno_all %>% fwrite("../output/tab_res_mmtt_metab_genoEffect.csv")
#lme_metab_mmtt_geno_pp %>% fwrite("../output/tab_res_mmtt_metab_genoEffect_pp.csv")
lme_metab_mmtt_geno <- fread("../output/tab_res_mmtt_metab_genoEffect.csv")


################################################################################
## Primary analysis stratified by GENOTYPE
################################################################################

genotypes <- c("HC genotype", "HF genotype")

lme_metab_mmtt_bygeno_all <- lapply(genotypes, function(geno) {
  lapply(1:length(metabolites), function(m) {
    metab <- metabolites[m]
    Metabolite <- met_info$Name[met_info$HMDB == metab][1]
    run_lme(exposure = "time", outcome=metab, 
            covariates = "time", coefficients_to_print = c(Time="time"), 
            data=analysis2 %>% filter(time %in% c(0,120,235) & genotype == geno), digits=c(1,3)) %>%
      mutate(outcome=metab, Outcome=Metabolite) %>%
      mutate(genotype = geno, .before=Exposure) %>%
      select(genotype, Exposure, outcome, Outcome, beta, se, p, anovaP, Effect, lowCI, upCI)
  }) %>% do.call(rbind.data.frame, .)
}) %>% do.call(rbind.data.frame, .)

## Save as .csv
#lme_metab_mmtt_bygeno_all %>% fwrite("../output/tab_res_mmtt_metab_timeEffect_byGeno.csv")
#lme_metab_mmtt_bygeno_pp %>% fwrite("../output/tab_res_mmtt_metab_timeEffect_byGeno_pp.csv")

lme_metab_mmtt_bygeno <- fread("../output/tab_res_mmtt_metab_timeEffect_byGeno.csv")
lme_metab_mmtt_bygeno_pp <- fread("../output/tab_res_mmtt_metab_timeEffect_byGeno_pp.csv")


# ==============================================================================
## Aggregate list of significant metabolites & fold changes, by genotype
# ===============================================================================

#dir.create("../data/processed/pathways")

##significant metabolites by time and genotype
genotypes=list(HC="HC genotype", HF="HF genotype")
times=c(m120="120 min", m235="235 min")

signif_metabs.l <- lapply(genotypes, function(g) {
  temp <-lapply(times, function(t) {
    lme_metab_mmtt_bygeno_pp %>% rename(time=Exposure) %>%
      mutate_at("time", ~paste0(gsub("Time_", "", .), " min")) %>%
      filter(genotype==g & time == t) %>%
      arrange(p) %>%
      ## Only keeping metabolites with P<0.05/30 significant & large changes for pathway enrichment ***
      filter(p<0.05/30 & abs(beta/log(2))>0.5) %>% 
      pull(outcome) %>% unique() }) ## NOTE: CONVERT BETA TO LOG2 FOLD CHANGE FOR THRESHOLD >0.5
  names(temp) <- names(times) ; temp
}) ; names(signif_metabs.l) <- names(genotypes)


### Extract metabolite list for common metabolites in ALL participants
#dir.create("../data/processed/pathways/metab_pEff_log2fc05")
#lme_metab_mmtt_timeonly <- fread("../data/processed/tab_prelim_lme_metab_mmtt_timeonly.csv")

lapply(c("120", "235"), function(time) {
  lme_metab_mmtt_timeonly %>% 
    filter(Exposure == paste0("Time_",time) & p<0.05/30 & abs(beta/log(2))>0.5) %>% 
    pull(outcome) %>%
    as.data.frame() %>% fwrite(paste0("../data/processed/pathways/all",time,"_pEff_log2fc05.txt"), col.names = F)
})


## Save list of ALL metabolites with LARGE & SIGNIF effects for each genotype and time point
signif_metabs.l$HC$m120 %>% as.data.frame() %>% fwrite("../data/processed/pathways/hc120_pEff_log2fc05.txt")
signif_metabs.l$HC$m235 %>% as.data.frame() %>% fwrite("../data/processed/pathways/hc235_pEff_log2fc05.txt")
signif_metabs.l$HF$m120 %>% as.data.frame() %>% fwrite("../data/processed/pathways/hf120_pEff_log2fc05.txt")
signif_metabs.l$HF$m235 %>% as.data.frame() %>% fwrite("../data/processed/pathways/hf235_pEff_log2fc05.txt")
union(signif_metabs.l$HC$m120, signif_metabs.l$HF$m120) %>% as.data.frame() %>% fwrite("../data/processed/pathways/m120_pEff_log2fc05.txt")
union(signif_metabs.l$HC$m235, signif_metabs.l$HF$m235) %>% as.data.frame() %>% fwrite("../data/processed/pathways/m235_pEff_log2fc05.txt")

# Union
intersect(signif_metabs.l$HC$m120, signif_metabs.l$HF$m120) %>% length() # 138
intersect(signif_metabs.l$HC$m235, signif_metabs.l$HF$m235) %>% length() # 131


## Look at the metabolite types reflected in each genotype/time
metab_dict %>% filter(Input.name %in% signif_metabs.l$HC$m120) %>% reframe(Super=n_pct(Sub.class))
metab_dict %>% filter(Input.name %in% c(intersect(signif_metabs.l$HC$m120, signif_metabs.l$HF$m120) )) %>% reframe(Super=n_pct(Sub.class))


## Merge full metabolite list with time effects for combined supplementary table
metab_suppl_info <- full_join(fread("../data/processed/met_info_processed.csv") %>% unique(.),
                              fread("../data/processed/tab_prelim_lme_metab_mmtt_timeonly.csv") %>% 
                                rename(HMDB=outcome), by="HMDB") %>%
  filter(Exposure != "Time_(joint)") %>%
  mutate(Beta.95CI=sprintf("%s (%s, %s)", round(beta,2), round(lowCI,2), round(upCI,2))) %>%
  select(Method, Compound, MZ, RT, Name, HMDB, Exposure, Beta.95CI, P) %>%
  pivot_wider(names_from = "Exposure", values_from=c(Beta.95CI, P)) %>%
  arrange(P_Time_120) %>%
  mutate(P_signif=ifelse(P_Time_120 < 0.05/30 | P_Time_235 < 0.05/30,1,0))
metab_suppl_info %>% fwrite("../output/tab_res_mmtt_metab_ppwithdescr.csv")


# ===============================================================================================
## After running Pathway Enrichment Analysis --> Find significant differences in BA metabolism
# ===============================================================================================

## Gather all bile acids
all_bile_acids <- c(metab_dict %>% filter(Main.class == "Bile acids") %>% 
  filter(HMDB_ID %in% metabolites_pp) %>% pull(HMDB_ID)) 

bile_acids <- metab_suppl_info %>% filter(P_signif == 1 & HMDB %in% all_bile_acids) %>%
  select(HMDB, Name, starts_with("P_")) %>% mutate(
    Abbrev = c("GLCA", "TDCA", "GDCA", "TCDCA", "GCDCA", "GCA", "TCA", "GUDCA", "DCA", "TUDCA", "Tα-MCA", "aHC", "CDCA", "CA")
  ) ; bile_acids 


bile_acids <- bile_acids %>% rowwise() %>%
  mutate(HC120_only = ifelse(HMDB %in% signif_metabs.l$HC$m120 & !HMDB %in% signif_metabs.l$HF$m120,1,0)) %>%
  mutate(HF120_only = ifelse(!HMDB %in% signif_metabs.l$HC$m120 & HMDB %in% signif_metabs.l$HF$m120,1,0)) %>%
  mutate(Both120_Geno = ifelse(HMDB %in% signif_metabs.l$HC$m120 & HMDB %in% signif_metabs.l$HF$m120,1,0)) %>%
  mutate(HC235_only = ifelse(HMDB %in% signif_metabs.l$HC$m235 & !HMDB %in% signif_metabs.l$HF$m235,1,0)) %>%
  mutate(HF235_only = ifelse(!HMDB %in% signif_metabs.l$HC$m235 & HMDB %in% signif_metabs.l$HF$m235,1,0)) %>%
  mutate(Both235_Geno = ifelse(HMDB %in% signif_metabs.l$HC$m235 & HMDB %in% signif_metabs.l$HF$m235,1,0))
bile_acids %>% fwrite("../data/processed/tab_descr_bileacids_mmtt.csv")


# ===========================================================
## Calculate 120 min fold change for all bile acids
# ===========================================================

bileacids_fc.df <- lapply(bile_acids$HMDB, function(ba) {
  analysis2 %>% select(id, genotype, time, all_of(ba)) %>%
    mutate_at("time", ~paste0("m", .)) %>%
    filter(time %in% c("m0", "m120")) %>%
    pivot_wider(names_from=time, values_from=ba) %>% 
    mutate(diff=m120-m0) %>%
    rename_with(., ~gsub("diff", ba, .)) %>%
    select(id, genotype, ba)
}) %>% reduce(full_join, by=c("id", "genotype"))


## Build dataset with significant metabolites, glucose, insulin and lipids
postprandial2 <- analysis %>% 
   #rename(tg_0=tg, hdl_0=hdl, ldl_0=ldl) %>%
  select(id, starts_with(c("glu", "insu", "tg", "ldl", "hdl")), -starts_with("tg_log"),
         genotype, meal_choice, genetic_cat_x_meal_choice,
         age, sex, bmi, PC1z, PC2z, PC3z, -contains("iAUC")) %>% 
  pivot_longer(starts_with(c("glu", "insu", "tg", "ldl", "hdl")), names_sep="_", names_to=c("metabolite", "time")) %>%
  filter(!is.na(time), !is.na(genotype)) %>%
  pivot_wider(names_from=metabolite) %>%
  mutate(time=factor(time, levels=c(0,30,60,120,180,235,270,300,360))) %>%
  arrange(time) ; head(postprandial2) 

## Add tg_log
postprandial_tg.log <- analysis %>% select(id, genotype, starts_with("tg_log")) %>%
  rename_with(., ~gsub("tg_log", "tg.log", .), starts_with("tg_log")) %>%
  rename(tg.log_0=tg.log) %>%
  pivot_longer(starts_with("tg.log_"), names_sep="_", names_to = c("metabolite", "time") ) %>%
  filter(!is.na(time), !is.na(genotype)) %>%
  pivot_wider(names_from=metabolite) %>%
  select(id, tg.log)

postprandial2 <- left_join(postprandial2, analysis2 %>% select(id, time, bile_acids$HMDB), by=c("id","time")) %>%
  left_join(postprandial_tg.log, by="id")


##################################################################################### 
# Post-hoc analyses of Bile Acid metabolites 
#####################################################################################

# Load in postprandial data from primary glycemic analyses
postprandial <- readRDS("../data/processed/postprandial_long.rda")

# ====================================================================================
## Prepare data with 120 fold change in BA, glucose/insulin and lipid metabolites
# ====================================================================================

ba <- left_join(
  # Bile Acids at 120
  postprandial2 %>% filter(time == 120) %>% select(id, age, sex, PC1z, PC2z, PC3z, bile_acids$HMDB,
                                                   "glucose", "insulin", "hdl", "ldl", "tg", "tg.log") %>%
    rename_at(c(bile_acids$HMDB, "glucose", "insulin", "hdl", "ldl", "tg", "tg.log"), ~paste0(., "_120")),
  # BA fold change from 0 to 120
  bileacids_fc.df %>% rename_at(bile_acids$HMDB, ~paste0(., "_120fc")), by = "id") %>% 
  # glucose and insulin 120 min iAUC
  #left_join(analysis %>% select(id, c(paste0(paste(rep(c("glucose", "insulin"), each=7), c(30, 60, 120, 180, 235, 270, 360), sep="_"), rep("iAUC.net",7)))), by="id") %>%
  # glucose and insulin at 30 min
  left_join(postprandial %>% filter(time == 30) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_30")), by="id") %>%
  # glucose and insulin at 60 min
  left_join(postprandial %>% filter(time == 60) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_60")), by="id") %>%
  # glucose and insulin at 60 min
  left_join(postprandial %>% filter(time == 180) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_180")), by="id") %>%
  # clinical metabolites at 235
  left_join(postprandial2 %>% filter(time == 235) %>% select(id, bile_acids$HMDB, "glucose", "insulin", "tg", "tg.log", "ldl", "hdl") %>% 
              rename_at(c(bile_acids$HMDB, "glucose", "insulin", "tg", "tg.log", "ldl", "hdl"), ~paste0(., "_235")), by="id")  %>%
  # bile acid metabolites at 0
  left_join(postprandial2 %>% filter(time == 0) %>% select(id, bile_acids$HMDB) %>% 
              rename_at(bile_acids$HMDB, ~paste0(., "_0")), by="id") %>%
  filter(complete.cases(genotype))


ba <- ba %>% 
  rowwise() %>%
  # Add total BA metabolites at 120 and 235
  mutate(totalBA_120 = rowSums(across(c(paste0(bile_acids$HMDB, "_120")))),
         totalBA_235 = rowSums(across(c(paste0(bile_acids$HMDB, "_235"))))) %>% 
  # Add total BA at 0
  left_join(postprandial2 %>% filter(time == 0) %>% select(id, bile_acids$HMDB) %>%
              mutate(totalBA_0 = rowSums(across(bile_acids$HMDB))) %>%
              select(id, totalBA_0), by="id") %>%
  mutate(totalBA_120fc = totalBA_120 - totalBA_0)


# =================================================================
# Compile all exposures and outcomes to text with LINEAR models
# =================================================================

## Exposures: 120 iAUC, 120 min levels, 120 min fold change
exposures <- paste0(bile_acids$HMDB, rep("_120fc",8))

## Outcomes:  120 iAUC, 120 min levels, 120 min fold change
outcomes <- c(paste0(rep(clinvars, each=2), rep(c("_120", "_235"), length(clinvars))),
              paste0(rep(c("glucose", "insulin"), each=4), rep(c("_60", "_180", "_120iAUC.net", "_60iAUC.net"))))


# ===================================================================
## Summary table of bile acid levels at 0, 120 and 235 by genotype
# ===================================================================

ba_formatted <- bile_acids$Abbrev ; names(ba_formatted) <- bile_acids$HMDB

lapply(c("0", "120", "235"), function(i) {
  vars <- ba_formatted ; names(vars) <- paste0(names(vars),"_",i)
  lapply(1:length(vars), function(v) {
    p = t.test(formula(paste(names(vars)[v],"~genotype")), data=ba)$p.value
    ba %>% group_by(genotype) %>% select(ba=names(vars)[v], genotype) %>%
      reframe(msd=mean_sd(ba)) %>% 
      mutate(BA=vars[[v]], ba=names(vars)[v]) %>%
      pivot_wider(names_from=genotype, values_from="msd") %>% as.data.frame() %>%
    mutate(p=round(p,3))
    }) %>% do.call(rbind.data.frame, .)
  }) %>% do.call(rbind.data.frame, .) %>%
  separate("ba", sep="_", into=c("ba", "time"))

## t-tests of 120 fild change in bile acids by genotype
lapply(names(bileacids.all), function(y) {
  summary(lm(formula(paste0(y, "_120fc~genotype")), data=ba))$coef
})


# =========================================================================
## Pearson correlations of 120m FC bile acids with glucose/insulin/TG
# =========================================================================

yvars <- c(paste(rep(c("glucose", "insulin"), each=11), c("30","60","120","180","235"), sep="_"),
           "tg_120", "tg_235")

ba_clinvars_corr <- lapply(exposures, function(x) {
  X <- met_info$Name[met_info$HMDB==gsub("_.*", "", x)][1]
  lapply(c(yvars), function(y) {
  cor <- cor.test(ba %>% pull(x), ba  %>%  pull(y))
  cbind.data.frame(BA=X, ba=x, y=y, cor=cor$estimate, p=cor$p.value, lowci=cor$conf.int[[1]], upci=cor$conf.int[[2]])
  }) %>% do.call(rbind.data.frame, .) 
})  %>% do.call(rbind.data.frame, .) %>% unique()
ba_clinvars_corr %>% arrange(p) %>% filter(endsWith(ba,"fc")) %>% filter(p<0.3)
#ba_clinvars_corr %>% fwrite("../output/tab_res_ba_clinvars_pearson.csv")


## Spearman corelations (given non-normality of distributions --> better suited for data)

ba_clinvars_corr_sp <- lapply(exposures, function(x) {
  X <- met_info$Name[met_info$HMDB==gsub("_.*", "", x)][1]
  lapply(c(yvars), function(y) {
    cor <- cor.test(ba %>% pull(x), ba %>% pull(y), method = "spearman")
    cbind.data.frame(BA=X, ba=x, y=y, cor=cor$estimate, p=cor$p.value)
  }) %>% do.call(rbind.data.frame, .) 
})  %>% do.call(rbind.data.frame, .) %>% unique()
ba_clinvars_corr_sp %>% arrange(p) %>% filter(endsWith(ba,"fc")) %>% filter(p<0.05)
#ba_clinvars_corr_sp %>% fwrite("../output/tab_res_ba_clinvars_spearman.csv")


## Test for genotype-specific correlations by stratifing by genotype
ba_clinvars_corr_bygeno <- lapply(genotypes, function(g) {
  lapply(exposures, function(x) {
    X <- met_info$Name[met_info$HMDB==gsub("_.*", "", x)][1]
    lapply(c(yvars), function(y) {
      cor <- cor.test(ba %>% filter(genotype ==g) %>% pull(x), ba %>% filter(genotype ==g) %>% pull(y))
      cbind.data.frame(Geno=g, BA=X, ba=x, y=y, cor=cor$estimate, p=cor$p.value, lowci=cor$conf.int[[1]], upci=cor$conf.int[[2]])
    }) %>% do.call(rbind.data.frame, .) 
  })  %>% do.call(rbind.data.frame, .) %>% unique()
}) %>% do.call(rbind.data.frame, .) ; ba_clinvars_corr_bygeno
ba_clinvars_corr_bygeno %>% arrange(p) %>% filter(endsWith(ba,"fc")) %>% filter(p<0.3)



#########################################################################################
## SELF-SELECTED MIXED MEALS: Run LM for change in metabolites to self-selected meals 
#########################################################################################

meals <- c("HC meal", "HF meal")
genotypes <- c("HC genotype", "HF genotype")


# ==============================================================================
## Exploratory metabolomics analysis stratified by meal type - Time only**
# ==============================================================================

lme_metab_ssmt_timeonly <- lapply(meals, function(meal) {
  lapply(1:length(metabolites), function(m) {
    metab <- metabolites[m]
    Metabolite <- met_info$Name[met_info$HMDB == metab][1]
    run_lme(exposure = "time", outcome=metab, 
            covariates = "time",coefficients_to_print = c(Time="time"), 
            data=analysis2 %>% filter(time %in% c(235,360) & meal_choice == meal), 
            digits=c(1,3)) %>%
      mutate(outcome=metab, Outcome=Metabolite) %>%
      mutate(meal_type = meal, genotype = "All", .before=Exposure) 
  }) %>% do.call(rbind.data.frame, .)
}) %>% do.call(rbind.data.frame, .) %>%
  dplyr::select(genotype, meal_type, Exposure, outcome, Outcome, beta, se, p, anovaP, Effect, lowCI, upCI)

## Save as .csv
#lme_metab_ssmt_timeonly %>% fwrite("../output/tab_res_ssmt_metab_timeonly_bymeal.csv")
lme_metab_ssmt_timeonly <- fread("../output/tab_res_ssmt_metab_timeonly_bymeal.csv")


# ==================================================================================================
## Calculate Kendall's Tao for agreement between first and second high(ER) carb meals
# Interpretting: https://blogs.sas.com/content/iml/2023/04/05/interpret-spearman-kendall-corr.html
# ==================================================================================================

# Tao for time effect in standardized vs HC self-selected MMTT
tao.dat <- inner_join(
  lme_metab_mmtt_timeonly %>%
    filter(Exposure == "Time_120") %>%
    arrange(beta) %>% mutate(mmtt = 1:nrow(.)) %>%
    dplyr::select(outcome, Outcome, mmtt), 
  lme_metab_ssmt_timeonly %>%
    filter(Exposure == "Time_360" & meal_type == "HC meal") %>%
    arrange(beta) %>% mutate(ssmt = 1:nrow(.)) %>%
    dplyr::select(outcome, Outcome, ssmt),
  by=c("Outcome", "outcome")) ; tao.dat

cor(tao.dat$mmtt, tao.dat$ssmt, method="kendall") # 0.457976 = *Moderate to strong* correlation !
cor.test(tao.dat$mmtt, tao.dat$ssmt, method="kendall") # 0.457976 = *Moderate to strong* correlation !
res <- cor.test(tao.dat$mmtt, tao.dat$ssmt, method = "kendall")
res$p.value


# ==============================================================================
## Time effect, stratified by meal type && genotype ********
# ==============================================================================

lme_metab_ssmt_bymealxgeno <- lapply(meals, function(meal) {
  lapply(genotypes, function(geno) {
    lapply(1:length(metabolites), function(m) {
      metab <- metabolites[m]
      Metabolite <- met_info$Name[met_info$HMDB == metab][1]
      run_lme(exposure = "time", outcome=metab, 
              covariates = "time", coefficients_to_print = c(Time="time"), 
              data=analysis2 %>% filter(time %in% c(235,360) & genotype == geno & meal_choice == meal), 
              digits=c(1,3)) %>%
        mutate(outcome=metab, Outcome=Metabolite) %>%
        mutate(meal_type = meal, genotype = geno, .before=Exposure) 
    }) %>% do.call(rbind.data.frame, .)
  }) %>% do.call(rbind.data.frame, .)
  }) %>% do.call(rbind.data.frame, .) %>%
  select(genotype, meal_type, Exposure, outcome, Outcome, beta, se, p, anovaP, Effect, lowCI, upCI)

## Save as .csv
#lme_metab_ssmt_bymealxgeno %>% fwrite("../output/tab_res_ssmt_metab_time_bymealxgeno.csv")
lme_metab_ssmt_bymealxgeno <- fread("../output/tab_res_ssmt_metab_time_bymealxgeno.csv")

lme_metab_ssmt_bymealxgeno %>% filter(p<0.05/23)


# ======================================================================
# Genotype effect, for each self-selected meal type
# ======================================================================

lme_metab_ssmt_bymeal_geno <- lapply(meals, function(meal) {
  lapply(1:length(metabolites), function(m) {
      metab <- metabolites[m]
      Metabolite <- met_info$Name[met_info$HMDB == metab][1]
      run_lme(exposure = "genotype*time", outcome=metab, 
              covariates = "genotype*time", coefficients_to_print = c(Time="time", Genotype="genotype"), 
              data=analysis2 %>% filter(time %in% c(235,360) & meal_choice == meal), digits=c(1,3)) %>%
        mutate(outcome=metab, Outcome=Metabolite) %>%
        mutate(meal_type = meal, .before=Exposure) 
    }) %>% do.call(rbind.data.frame, .)
  }) %>% do.call(rbind.data.frame, .) %>%
  filter(grepl(" x ", Exposure)) %>%
  mutate_at("Exposure", ~gsub(".* x", "HC Genotype x", .)) %>%
  select(meal_type, Exposure, outcome, Outcome, beta, se, p, anovaP, Effect, lowCI, upCI)

## Save as .csv
#lme_metab_ssmt_bymeal_geno %>% fwrite("../output/tab_res_ssmt_metab_bymeal_geno.csv")
lme_metab_ssmt_bymeal_geno <- fread("../output/tab_res_ssmt_metab_bymeal_geno.csv")


# ========================================================
## Pathway Enrichment Analysis for self-selected meals
# ========================================================

##significant metabolites by time and genotype
genotypes=list(HC="HC genotype", HF="HF genotype")
meals=c(HC="HC meal", HF="HF meal")

ssmt_metabs.l <- lapply(genotypes, function(g) {
    temp <-lapply(meals, function(m) {
      lme_metab_ssmt_bymealxgeno %>% rename(time=Exposure) %>%
        filter(time == "Time_360" & genotype==g & meal_type == m) %>%
        arrange(p) %>%
        filter(p<0.05/30 & abs(beta/log(2))>0.5) %>% pull(outcome) %>% unique() }) ## NOTE: CONVERT BETA TO LOG2 FOLD CHANGE FOR THRESHOLD >0.5 **TRYING >0.3**
    names(temp) <- names(meals) ; temp
  }) ; names(ssmt_metabs.l) <- names(genotypes)
ssmt_metabs.l


## List of metabolites changing for each genotype on each meal
ssmt_metabs.l$HC$HC %>% as.data.frame() %>% fwrite("../data/processed/pathways/hchc_pEff_log2fc05.txt")
ssmt_metabs.l$HC$HF %>% as.data.frame() %>% fwrite("../data/processed/pathways/hchf_pEff_log2fc05.txt")
ssmt_metabs.l$HF$HC %>% as.data.frame() %>% fwrite("../data/processed/pathways/hfhc_pEff_log2fc05.txt")
ssmt_metabs.l$HF$HF %>% as.data.frame() %>% fwrite("../data/processed/pathways/hfhf_pEff_log2fc05.txt")


# Are any BA changing in the ssmt?
lme_metab_ssmt_bymealxgeno %>% filter(outcome %in% names(bileacids.all)) %>%
  filter(p<0.05) %>%
  filter(p<0.05/30)

bileacids_fc_ssmt.df <- lapply(names(bileacids.all), function(x) {
  analysis2 %>% select(id, genotype, time, meal_choice, names(bileacids.all)) %>%
    filter(time %in% c(235, 360)) %>%
    pivot_longer(names(bileacids.all), values_to="value", names_to="bileacids") %>%
    unite("ba_time", bileacids, time, sep="_") %>%
    pivot_wider(names_from=ba_time, values_from="value")  %>%
    
    select(id, mt1=paste0(x, "_235"), mt2=paste0(x, "_360"), genotype, meal_choice) %>%
    mutate(mdiff=mt2 - mt1) %>%
    rename_with(., ~gsub("mdiff", paste0(x, "_360fc"), .)) %>%
    select(id, ends_with("fc"), genotype, meal_choice)
  }) %>% reduce(full_join,by=c("id","genotype", "meal_choice")) ; bileacids_fc_ssmt.df



# ==================================================
## Checking Bile acids in self-selection meal
# ==================================================

## Buld ba dataframe with long.form
ba.ssmt <- left_join(
  # Bile Acids at 235,360
  postprandial2 %>% filter(time == 235) %>% select(id, age, sex, PC1z, PC2z, PC3z, bile_acids$HMDB,
                                                   "glucose", "insulin", "hdl", "ldl", "tg") %>%
    rename_at(c(bile_acids$HMDB, "glucose", "insulin", "hdl", "ldl", "tg"), ~paste0(., "_235")),
  bileacids_fc_ssmt.df, by = "id") %>% 
  left_join(postprandial %>% filter(time == 270) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_270")), by="id") %>%
  left_join(postprandial %>% filter(time == 300) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_300")), by="id") %>%
  left_join(postprandial2 %>% filter(time == 360) %>% select(id, bile_acids$HMDB, "glucose", "insulin", "tg", "ldl", "hdl") %>% 
              rename_at(c(bile_acids$HMDB, "glucose", "insulin", "tg", "ldl", "hdl"), ~paste0(., "_360")), by="id")  %>%
  left_join(postprandial2 %>% filter(time == 0) %>% select(id, bile_acids$HMDB) %>% 
              rename_at(bile_acids$HMDB, ~paste0(., "_0")), by="id") %>%
  filter(complete.cases(genotype))

ba.ssmt <- ba.ssmt %>% 
  rowwise() %>%
  # Add total BA metabolites at 120 and 235
  mutate(totalBA_235 = rowSums(across(c(paste0(bile_acids$HMDB, "_235")))),
         totalBA_360 = rowSums(across(c(paste0(bile_acids$HMDB, "_360"))))) %>% 
  # Add total BA at 0
  left_join(postprandial2 %>% filter(time == 0) %>% select(id, bile_acids$HMDB) %>%
              mutate(totalBA_0 = rowSums(across(bile_acids$HMDB))) %>%
              select(id, totalBA_0), by="id") %>%
  mutate(totalBA_360fc = totalBA_360 - totalBA_235)


# Correlations of bile acids with insulin at 360?
yvars <- c(paste(rep(c("glucose", "insulin"), each=4), c("235", "270", "300", "360"), sep="_"), "tg_235", "tg_360")
ba_ssmt <- paste0(names(bileacids.all), "_360fc")

ba_clinvars_ssmt_sp <- lapply(ba_ssmt, function(x) {
  X <- met_info$Name[met_info$HMDB==gsub("_.*", "", x)][1]
  lapply(c(yvars), function(y) {
    cor <- cor.test(ba.ssmt %>% pull(x), ba.ssmt  %>%  pull(y))
    cbind.data.frame(BA=X, ba=x, y=y, cor=cor$estimate, p=cor$p.value, lowci=cor$conf.int[[1]], upci=cor$conf.int[[2]])
  }) %>% do.call(rbind.data.frame, .) 
})  %>% do.call(rbind.data.frame, .) %>% unique()
ba_clinvars_ssmt_sp %>% arrange(p) %>% filter(endsWith(ba,"fc")) %>% filter(p<0.05)
#ba_clinvars_corr %>% fwrite("../output/tab_res_ba_clinvars_pearson.csv")


## bile acids following SSMT
bileacids_fc_ssmt.df <- lapply(names(bileacids.all), function(ba) {
  analysis2 %>% select(id, genotype, meal_choice, time, ba) %>%
    mutate_at("time", ~paste0("m", .)) %>%
    filter(time %in% c("m235", "m360")) %>%
    pivot_wider(names_from=time, values_from=ba) %>% 
    mutate(diff=m360-m235) %>%
    rename_with(., ~gsub("diff", ba, .)) %>%
    select(id, genotype, meal_choice, ba)
}) %>% reduce(full_join, by=c("id", "genotype", "meal_choice"))


## Table of 120 min fold change by genotype & meal type
bileacids_fc_ssmt.df %>% 
  pivot_longer(bile_acids_to_plot) %>% 
  #pivot_longer(names(bileacids.all)) %>% 
  mutate_at("name", ~bileacids.all[name]) %>%
  mutate_at("genotype", ~relevel(., ref="HC genotype")) %>%
  group_by(genotype, meal_choice, name) %>% filter(!is.na(genotype)) %>%
  reframe(n=n(), mean=mean(exp(value)), sd=sd(exp(value))) %>% mutate(se=mean/sqrt(n)) %>%
  
  ggplot(aes(x=name, group=genotype, y=mean, ymin=mean-se, ymax=mean+se, color=genotype)) +
  facet_grid(rows = vars(meal_choice)) +
  geom_hline(yintercept = 1, linewidth=0.25, color="black") +
  geom_point(size=2, position=position_dodge(0.35)) + geom_errorbar(width=0.15, linewidth=0.45, position=position_dodge(0.35)) +
  scale_color_manual(values=parameters$geno$palette, name=parameters$geno$label) + 
  xlab("Bile Acid Metabolites") + ylab("120 min fold change")


## Table of 120 min fold change by genotype & meal type
bileacids_fc.df %>% 
  pivot_longer(names(bileacids.all)) %>% 
  mutate_at("name", ~bileacids.all[name]) %>%
  mutate_at("genotype", ~relevel(., ref="HC genotype")) %>%
  group_by(genotype, name) %>% filter(!is.na(genotype)) %>%
  reframe(n=n(), mean=mean(exp(value)), sd=sd(exp(value))) %>% mutate(se=mean/sqrt(n)) %>%
  
  ggplot(aes(x=name, group=genotype, y=mean, ymin=mean-se, ymax=mean+se, color=genotype)) +
  geom_hline(yintercept = 1, linewidth=0.25, color="black") +
  geom_point(size=2, position=position_dodge(0.35)) + geom_errorbar(width=0.15, linewidth=0.45, position=position_dodge(0.35)) +
  scale_color_manual(values=parameters$geno$palette, name=parameters$geno$label) 


## EOF

