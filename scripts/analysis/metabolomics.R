## Rscript for metabolomics data analysis

library(tidyverse)
library(data.table)
library(poolr)
library(ComplexHeatmap)

lapply(list.files("../../pantry/functions/", recursive = T, full.names = T), source)


################################################################################ 
## Load & build metabolomics dataset 
################################################################################

# analysis dataset with metabolomics 
analysis2 <- readRDS("../data/processed/premier_analysis2.rda") %>%
  filter(!is.na(genotype))

# confirm all 22 participants have metabolomics data
metabolomics_ids <- analysis2 %>% filter(has_metabolomics==1) %>% 
  pull(id) %>% unique() ; length(metabolomics_ids) # N = 22

# analysis dataset with phenotypes & glycemic traits
analysis <- readRDS("../data/processed/premier_analysis.rda") %>%
  filter(!is.na(genotype)) %>%
  mutate(has_metabolomics = ifelse(id %in% metabolomics_ids,1,0)) %>%
  # add lipid levels
  mutate(tg_0=tg, tg_log_0=tg_log, hdl_0=hdl, ldl_0=ldl)


# load metabolomics info file with descriptive metabolite names, compound IDs, weights
met_info <- fread("../data/processed/met_info_processed.csv")
metabolites <- names(analysis2 %>% select(starts_with("HMDB"))) #817 TOTAL metabolites
metabolites_pp <- fread("../data/processed/metbolites_pp.txt")$metabolites_pp # 547 changing postprandially (p<0.05)
metabclass <- fread("../data/processed/metabolites_subclasses.csv") %>%
  rename(HMDB=HMDBID)

## load external metabolomics data dictionary: RefMet
#devtools::install_github("metabolomicsworkbench/RefMet")
library(RefMet)
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

# ========================================================
## Check distributions of log2-transformed metabolites
# ========================================================

metab_pp_shapiro_time0 <- lapply(metabolites_pp, function(metab) {
  cbind.data.frame(HMDB=metab, p=shapiro.test(
    analysis2 %>% filter(time == 0) %>% pull(metab))$p.value)
  }) %>% do.call(rbind.data.frame, .) %>%
  mutate(p_level = case_when(
    p>=0.05 ~ "Normal", p<0.05 & p>=0.01 ~ "1",
    p <0.01 & p>= 0.001 ~ "2", p<0.001 ~ "3")) 

metab_pp_shapiro_time0 %>%
  ggplot(aes(x=1, fill=p_level)) + 
  geom_bar(stat="count")

prop.table(table(metab_pp_shapiro_time0$p_level))


################################################################################
## Primary analysis: changes in metabolomics profile after standardized MMTT
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

## Save as .csv
lme_metab_mmtt_geno_all %>% fwrite("../output/tab_res_mmtt_metab_genoEffect.csv")

# Restrict to postprandial metabolites
lme_metab_mmtt_geno_pp <- lme_metab_mmtt_geno %>% filter(outcome %in% metabolites_pp)
lme_metab_mmtt_geno_pp %>% fwrite("../output/tab_res_mmtt_metab_genoEffect_pp.csv")


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
lme_metab_mmtt_bygeno_all %>% fwrite("../output/tab_res_mmtt_metab_timeEffect_byGeno.csv")

# Restrict to postprandial metabolites
lme_metab_mmtt_bygeno_pp %>% fwrite("../output/tab_res_mmtt_metab_timeEffect_byGeno_pp.csv")


# ========================================================================
## Aggregate list of significant metabolites & fold changes, by genotype
# ========================================================================

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
      ## Filter to metabolites with P<0.05/30 & log2fc >0.5 for pathway enrichment ***
      filter(p<0.05/30 & abs(beta/log(2))>0.5) %>% 
      pull(outcome) %>% unique() }) 
  names(temp) <- names(times) ; temp
}) ; names(signif_metabs.l) <- names(genotypes)


### Extract metabolite list for common metabolites in ALL participants
#dir.create("../data/processed/pathways/metab_pEff_log2fc05")

lapply(c("120", "235"), function(time) {
  lme_metab_mmtt_timeonly %>% 
    filter(Exposure == paste0("Time_",time) & p<0.05/30 & abs(beta/log(2))>0.5) %>% 
    pull(outcome) %>%
    as.data.frame() %>% 
    fwrite(paste0("../data/processed/pathways/all",time,"_pEff_log2fc05.txt"), col.names = F)
})


## Save list of ALL metabolites with signif & large effects for each PGS & time
signif_metabs.l$HC$m120 %>% as.data.frame() %>% fwrite("../data/processed/pathways/hc120_pEff_log2fc05.txt")
signif_metabs.l$HC$m235 %>% as.data.frame() %>% fwrite("../data/processed/pathways/hc235_pEff_log2fc05.txt")
signif_metabs.l$HF$m120 %>% as.data.frame() %>% fwrite("../data/processed/pathways/hf120_pEff_log2fc05.txt")
signif_metabs.l$HF$m235 %>% as.data.frame() %>% fwrite("../data/processed/pathways/hf235_pEff_log2fc05.txt")
union(signif_metabs.l$HC$m120, signif_metabs.l$HF$m120) %>% as.data.frame() %>% fwrite("../data/processed/pathways/m120_pEff_log2fc05.txt")
union(signif_metabs.l$HC$m235, signif_metabs.l$HF$m235) %>% as.data.frame() %>% fwrite("../data/processed/pathways/m235_pEff_log2fc05.txt")

# Compare metabolites across PGS and time points
intersect(signif_metabs.l$HC$m120, signif_metabs.l$HF$m120) %>% length() # 138
intersect(signif_metabs.l$HC$m235, signif_metabs.l$HF$m235) %>% length() # 131
metab_dict %>% filter(Input.name %in% signif_metabs.l$HC$m120) %>% reframe(Super=n_pct(Sub.class))
metab_dict %>% filter(Input.name %in% c(intersect(signif_metabs.l$HC$m120, signif_metabs.l$HF$m120) )) %>% 
  reframe(Super=n_pct(Sub.class))

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


# ===============================================================
## Based on enrichment analysis, examine bile acid metabolites 
# ===============================================================

# ------------------------------------------------
## Filtering by p<0.05/30 & log2fc > 0.5
# ------------------------------------------------

## Gather ALL bile acids
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

# List only PRIMARY bile acids
bile_acids_primary <- bile_acids %>% filter(Abbrev %in% c("CA", "GCA", "TCA", "CDCA", "GCDCA", "TCDCA"))

# ===========================================================
## Calculate 120 min log2 fold change for all bile acids
# ===========================================================

bileacids_log2fc.df <- lapply(bile_acids$HMDB, function(ba) {
  # 120min fold change
  analysis2 %>% select(id, genotype, time, all_of(ba)) %>%
    mutate_at("time", ~paste0("m", .)) %>%
    filter(time %in% c("m0", "m120")) %>%
    pivot_wider(names_from=time, values_from=ba) %>% 
    mutate(diff_120fc=(m120-m0)/log(2)) %>% # log2 120 min fold change
    rename_with(., ~gsub("diff", ba, .)) %>%
    select(id, genotype, paste0(ba,"_120fc")) %>% left_join(
      #235 min fold change
  analysis2 %>% select(id, genotype, time, all_of(ba)) %>%
    mutate_at("time", ~paste0("m", .)) %>%
    filter(time %in% c("m0", "m235")) %>%
    pivot_wider(names_from=time, values_from=ba) %>% 
    mutate(diff_235fc=(m235-m0)/log(2)) %>% # 235 min fold change (from 0)
    rename_with(., ~gsub("diff", ba, .)) %>%
    select(id, genotype, paste0(ba,"_235fc")), by = c("id", "genotype"))
}) %>% reduce(full_join, by=c("id", "genotype"))


## Build dataset with metabolites, glucose, insulin and lipids
postprandial2 <- analysis %>% 
   #rename(tg_0=tg, hdl_0=hdl, ldl_0=ldl) %>%
  select(id, starts_with(c("glu", "insu", "tg", "ldl", "hdl")), -starts_with("tg_log"),
         genotype, meal_choice, genetic_cat_x_meal_choice,
         age, sex, bmi, PC1z, PC2z, PC3z, -contains("iAUC")) %>% 
  pivot_longer(starts_with(c("glu", "insu", "tg", "ldl", "hdl")), names_sep="_", names_to=c("metabolite", "time")) %>%
  filter(!is.na(time), !is.na(genotype)) %>%
  pivot_wider(names_from=metabolite) %>%
  mutate(time=factor(time, levels=c(0,30,60,120,180,235,270,300,360))) %>%
  mutate(tg.log = log(tg)) %>%
  arrange(time) ; head(postprandial2) 

postprandial2 <- left_join(postprandial2, analysis2 %>% select(id, time, bile_acids$HMDB), by=c("id","time"))
head(postprandial2)


################################################################################ 
# Post-hoc analyses of BA metabolites 
################################################################################

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
  bileacids_log2fc.df, by = "id") %>% 
  # glucose and insulin at 0 min
  left_join(postprandial %>% filter(time == 0) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_0")), by="id") %>%
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
  left_join(postprandial %>% filter(time == 0) %>% select(id, tg) %>% mutate(tg.log_0 = log(tg)), by = "id") %>% 
  filter(complete.cases(genotype)) %>%
  # Add changes for each clinical metabolites
  mutate(
    glucose_0to60 = glucose_60 - glucose_0, glucose_0to120 = glucose_120 - glucose_0, 
    glucose_0to180 = glucose_180 - glucose_0, glucose_0to235 = glucose_235 - glucose_0,
    insulin_0to60 = insulin_60 - insulin_0, insulin_0to120 = insulin_120 - insulin_0, 
    insulin_0to180 = insulin_180 - insulin_0, insulin_0to235 = insulin_235 - insulin_0,
    tg.log_0to120 = tg.log_120 - tg.log_0, tg.log_0to235 = tg.log_235 - tg.log_0 
  )


# =================================================================
# Compile all exposures and outcomes to text with LINEAR models
# =================================================================

## Exposures: 120 iAUC, 120 min levels, 120 min fold change
exposures <- paste0(bile_acids_primary$HMDB, rep("_120fc",14))

## Outcomes: 120 iAUC, 120 min levels, 120 min fold change
outcomes <- c(paste0(rep(clinvars, each=2), rep(c("_120", "_235"), length(clinvars))),
              paste0(rep(c("glucose", "insulin"), each=4), rep(c("_60", "_180", "_120iAUC.net", "_60iAUC.net"))))

# ===================================================================
## Summary table of bile acid levels at 0, 120 and 235 by genotype
# ===================================================================

ba_formatted <- bile_acids_primary$Abbrev ; names(ba_formatted) <- bile_acids_primary$HMDB

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

## t-tests of 120 fold change in bile acids by genotype
lapply(bile_acids_primary$HMDB, function(y) {
  summary(lm(formula(paste0(y, "_120fc~genotype+age+sex+PC1z+PC2z+PC3z")), data=ba))$coef[2,]
}) %>% do.call(rbind, .) %>% cbind(., bile_acids_primary$Name)


# =========================================================================
## Correlations of 120m log2 FC bile acids with glucose/insulin/TG
# =========================================================================

yvars <- c(paste(rep(c("glucose", "insulin"), each=11), c("30","60","120","180","235"), sep="_"),
           "tg.log_120", "tg.log_235")

yvars_change <- c(paste(rep(c("glucose", "insulin"), each=11), 
                        c("0to60","0to120","0to180","0to235"), sep="_"),
                  "tg.log_0to120", "tg.log_0to235")

# 120 min fold change -------
ba_clinvars_corr_sp_120fc <- lapply(exposures, function(x) {
  X <- met_info$Name[met_info$HMDB==gsub("_.*", "", x)][1]
  lapply(c(yvars), function(y) {
    cor <- cor.test(ba %>% pull(x), ba %>% pull(y), method = "spearman")
    cbind.data.frame(BA=X, ba=x, y=y, cor=cor$estimate, p=cor$p.value)
  }) %>% do.call(rbind.data.frame, .) 
})  %>% do.call(rbind.data.frame, .) %>% unique()
ba_clinvars_corr_sp_120fc %>% arrange(p) %>% filter(endsWith(ba,"fc")) %>% filter(p<0.05)
ba_clinvars_corr_sp_120fc %>% fwrite("../output/tab_res_ba_clinvars_xsect_spearman_120log2fc.csv")
#ba_clinvars_corr_sp_xsect_120fc <- fread("../output/tab_res_ba_clinvars_xsect_spearman_120log2fc.csv")


# 235 min fold change ---------------------------------
ba_clinvars_corr_sp_235fc <- lapply(gsub("120","235",exposures), function(x) {
  X <- met_info$Name[met_info$HMDB==gsub("_.*", "", x)][1]
  lapply(c(yvars), function(y) {
    cor <- cor.test(ba %>% pull(x), ba %>% pull(y), method = "spearman")
    cbind.data.frame(BA=X, ba=x, y=y, cor=cor$estimate, p=cor$p.value)
  }) %>% do.call(rbind.data.frame, .) 
})  %>% do.call(rbind.data.frame, .) %>% unique()
ba_clinvars_corr_sp_235fc %>% arrange(p) %>% filter(endsWith(ba,"fc")) %>% filter(p<0.05)


# Test for genotype*time (0,120) interactions for each bile acid?
lapply(1:length(bile_acids_primary$HMDB), function(y) {
  summary(lmerTest::lmer(formula(paste0(bile_acids_primary$HMDB[y],"~genotype*time+age+sex+PC1z+PC2z+PC3z+(1|id)")), 
                       data=analysis2 %>% filter(time %in% c(0,120))))$coef[9,]
}) %>% do.call(rbind, .)


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

# Tab for time effect in standardized vs HC self-selected MMTT
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

lme_metab_ssmt_bymealxgeno %>% filter(p<0.05/30)


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

# ==================================================
## Summarize log2 fc in BA by meal x genotype
# ==================================================

bileacids_log2fc_ssmt.df <- lapply(names(bileacids.all), function(x) {
  analysis2 %>% select(id, genotype, time, meal_choice, names(bileacids.all)) %>%
    filter(time %in% c(235, 360)) %>%
    pivot_longer(names(bileacids.all), values_to="value", names_to="bileacids") %>%
    unite("ba_time", bileacids, time, sep="_") %>%
    pivot_wider(names_from=ba_time, values_from="value")  %>%
    
    select(id, mt1=paste0(x, "_235"), mt2=paste0(x, "_360"), genotype, meal_choice) %>%
    mutate(mdiff= (mt2 - mt1)/log(2)) %>%
    rename_with(., ~gsub("mdiff", paste0(x, "_360fc"), .)) %>%
    select(id, ends_with("fc"), genotype, meal_choice)
}) %>% reduce(full_join,by=c("id","genotype", "meal_choice")) ; bileacids_log2fc_ssmt.df



## Build ba dataframe with long.form
ba.ssmt <- left_join(
  # Bile Acids at 235,360
  postprandial2 %>% filter(time == 235) %>% select(id, age, sex, PC1z, PC2z, PC3z, bile_acids$HMDB,
                                                   "glucose", "insulin", "hdl", "ldl", "tg") %>%
    rename_at(c(bile_acids$HMDB, "glucose", "insulin", "hdl", "ldl", "tg"), ~paste0(., "_235")),
  bileacids_log2fc_ssmt.df, by = "id") %>% 
  left_join(postprandial %>% filter(time == 270) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_270")), by="id") %>%
  left_join(postprandial %>% filter(time == 300) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_300")), by="id") %>%
  left_join(postprandial2 %>% filter(time == 360) %>% select(id, bile_acids$HMDB, "glucose", "insulin", "tg", "ldl", "hdl") %>% 
              rename_at(c(bile_acids$HMDB, "glucose", "insulin", "tg", "ldl", "hdl"), ~paste0(., "_360")), by="id")  %>%
  left_join(postprandial2 %>% filter(time == 0) %>% select(id, bile_acids$HMDB) %>% 
              rename_at(bile_acids$HMDB, ~paste0(., "_0")), by="id") %>%
  filter(complete.cases(genotype))

# Correlations of bile acids with insulin at 360?
yvars <- c(paste(rep(c("glucose", "insulin"), each=4), c("235", "270", "300", "360"), sep="_"), "tg_235", "tg_360")
ba_ssmt <- paste0(bile_acids_primary$HMDB, "_360fc")

ba_clinvars_ssmt_sp <- lapply(ba_ssmt, function(x) {
  X <- met_info$Name[met_info$HMDB==gsub("_.*", "", x)][1]
  lapply(c(yvars), function(y) {
    cor <- cor.test(ba.ssmt %>% pull(x), ba.ssmt  %>% pull(y), method="spearman")
    cbind.data.frame(BA=X, ba=x, y=y, cor=cor$estimate, p=cor$p.value)
  }) %>% do.call(rbind.data.frame, .) 
}) %>% do.call(rbind.data.frame, .) %>% unique()
ba_clinvars_ssmt_sp %>% arrange(p) %>% filter(endsWith(ba,"fc")) %>% filter(p<0.05)
#ba_clinvars_corr %>% fwrite("../output/tab_res_ba_clinvars_pearson.csv")

# Test for genotype*time (235,360) interactions for each bile acid?
lapply(c("HC meal", "HF meal"), function(m) {
  lapply(1:length(bile_acids_primary$HMDB), function(y) {
  summary(lmerTest::lmer(formula(paste0(bile_acids_primary$HMDB[y],"~genotype*time+age+sex+PC1z+PC2z+PC3z+(1|id)")), 
                         data=analysis2 %>% filter(time %in% c(235,360)) %>%
                           filter(meal_choice == m)))$coef[9,]
  }) %>% do.call(rbind, .) %>%
    as.data.frame() %>%
    mutate(bile_acid = bile_acids_primary$Abbrev, meal = m, .before=1)
  }) %>% 
  do.call(rbind.data.frame, .) %>%
  fwrite("../data/processed/tab_ssmt_lme_int_ba_primary_genoxtime.csv") # no significant genotype x time interactions 

## EOF

