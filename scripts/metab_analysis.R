## Analyze metabolomics data

library(tidyverse) ; library(data.table)
lapply(list.files("../../pantry/functions/", recursive = T, full.names = T), source)

library(poolr)
library(ComplexHeatmap)

######################################################### 
## Load & build metabolomics dataset
######################################################### 

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
met_info <- fread("../data/processed/met_info.csv") %>%
  mutate(Compound = gsub("_.*", "", Compound_Id), .after="Compound_Id") %>%
  mutate(HMDB = gsub("_.*","",HMDB_Id)) 
metabolites <- names(analysis2 %>% select(starts_with("HMDB")))
metabclass <- fread("../data/raw/metabolomics/metabolites_subclasses.csv") %>%
  rename(HMDB=HMDBID)

## load external metabolomics data dictionary: RefMet
#devtools::install_github("metabolomicsworkbench/RefMet")
library(RefMet)
metab_dict <- refmet_map_df(gsub("_.*", "", metabolites))


###############################################################################
## Describe lipid levels at each time point by genotype (& meal selection)
###############################################################################

## Set lists of variables & strata -----------------

# Demographics
demographic_vars <- c(
  age="Age, years", Sex="Sex", Race2="Race", Ethnicity="Ethnicity", bmi="BMI, kg/m2", 
  whr="Waist-Hip", sbp="SBP, mmHg", dbp="DBP, mmHg", hba1c="HbA1c")

# Fasting metabolites
metabolite_vars <- c(
  glucose="Glucose, mg/dL", insulin="Insulin", tg_log="log(Triglyceride)", tg="Triglyceride",
  chol="Total cholesterol", hdl="HDL cholesterol", ldl="LDL cholesterol", protein="Total protein", alb="Albumin", 
  globulin="Globulin", alb_globulin_ratio="Albumin:Globulin", tot_bilirubin_log="log(Total Bilirubin)", 
  dir_bilirubin_log="log(Direct Bilirubin)", ast="AST", alt_log="Log(ALT)")

# Demographics
survey_vars <- c(
  education="Education level", smoke_packyears="Smoking Status", alcohol_perweek_4lvl="Alcohol Use", 
  sleep_hrday="Sleep, hrs/day", physact_mvpa_perweek="Physical Activity, hr/week")

# Summary table descriptive vars
summary_table_vars.l <- list(
  list(label="demographics", vars=demographic_vars), 
  list(label="metabolites", vars=metabolite_vars),
  list(label="survey", vars=survey_vars)) ; summary_table_strata.l <- list(
  list(label="genecat", var="genotype",order=c("HC genotype", "HF genotype")),
  list(label="mealchoice", var="meal_choice", order=c("HC meal", "HF meal")))

# Clinical metabolites
clinical_metabolites.labs <- paste(
  rep(c("Glucose", "Insulin", "TG", "log(TG)", "HDL", "LDL"),each=4),
  c("at 0 min", "at 120 min", "at 235 min", "at 360 min")) ; clinical_metabolites <- paste0(
  rep(c("glucose", "insulin", "tg", "tg_log", "hdl", "ldl"),each=4),
  c("_0", "_120", "_235", "_360")) ; names(clinical_metabolites.labs) <- clinical_metabolites

clinvars <- c("glucose", "insulin", "tg", "ldl", "hdl")
clinvars.lab <- c("Glucose", "Insulin", "TG", "LDL", "HDL")
names(clinvars.lab) <- clinvars


# differences by genotype
print_summary_table(vars_to_summarize = clinical_metabolites.labs,
                    var_strata = "genotype", var_strata_order = c("HC genotype", "HF genotype"),
                    data=analysis %>% filter(has_metabolomics==1)) 

# differences by meal choice; Pvalue=t.test; Ptrend=lm~var+age+sex
print_summary_table(vars_to_summarize = clinical_metabolites.labs, 
                    var_strata = "meal_choice", var_strata_order = c("HC meal", "HF meal"),
                    digits = c(1,1,3), data=analysis) %>%
  fwrite("../output/tab_descr_glycemic&lipids_ssmt_bymeal.csv", row.names = T)

# post-doc tests adjusting for genotype
t.test(insulin_235~meal_choice, data=analysis)
summary(lm(insulin_235~meal_choice+genotype, data=analysis))
summary(lm(tg_log_235~meal_choice+genotype, data=analysis))

# differences by genotype x meal choice
print_summary_table(vars_to_summarize = clinical_metabolites.labs, var_strata = "genetic_cat_x_meal_choice", 
                    var_strata_order = c("HC genotype x HC meal", "HC genotype x HF meal", 
                                         "HF genotype x HC meal", "HF genotype x HF meal"),
                    data=analysis %>% filter(has_metabolomics==1)) 

## Descriptive statistics among participants with metabolomics data
lapply(summary_table_strata.l, function(strata){
  lapply(summary_table_vars.l, function(vars) {
    print_summary_table(
      analysis %>% filter(has_metabolomics==1), vars_to_summarize = vars$vars, p_adjust = "none", digits = c(1,1,3),
      var_strata = strata$var, var_strata_order = strata$order,
      p_types = "descriptive", p_smalln=T) #%>%
      #fwrite(paste0("../output/tab_descr_", vars$label, "_", strata$label, "_metabolomics.csv"), row.names = T)
  })
})

table(analysis$genotype, analysis$has_metabolomics)
prop.table(table(analysis$genotype, analysis$has_metabolomics))
fisher.test(analysis$has_metabolomics, analysis$genotype)

## Compare characteristics BY having metabolomics data
lapply(summary_table_vars.l, function(vars) {
  print_summary_table(
    analysis, vars_to_summarize = vars$vars, p_adjust = "none", digits = c(1,1,3),
    var_strata = "has_metabolomics", #var_strata_order = c("Yes"=1, "No"=0),
    p_types = "descriptive", p_smalln=T) %>%
      #fwrite(paste0("../output/tab_descr_", vars$label, "_", strata$label, "_metabolomics.csv"), row.names = T)
  mutate_all(., ~gsub("[%]", "", .), .)
  }) %>% do.call(rbind.data.frame, .) #%>%
  #fwrite("../output/tab_descr_all_metabolomics_YesNo.csv", row.names = T, col.names = T)

# ==============================================================================
## Determine Number of Effective Tests (number of NON-correlated metabolites)  
## for multiple comparison correction = 25
# ==============================================================================

library(poolr)

# Known metabolites = 366
metabolites.df <- analysis2 %>% filter(time == 0) %>% select(all_of(metabolites)) %>%
  filter(complete.cases(.)) # Remove NAs ; dim(metabolites.df): 17 x 677

# Effective tests = 25
corr_metab <- cor(metabolites.df) ## create correlation matrix
effective_tests <- meff(corr_metab, method="liji") ## 23 un-correlated metabolites 
effective_tests #23

## Heat Map of metabolite correlations 
#heatmap(corr_metab)


################################################################################
## Primary analysis: changes in metabolomics profile after standrdized MMTT
################################################################################

# =========================================================
## 2) Run LME for metab ~ genotype x time (categorical)
## Restricting to postprandial metabolites
# =========================================================

lme_metab_mmtt_unadj <- lapply(1:length(metabolites), function(m) {
  metab <- metabolites[m]
  Metabolite <- met_info$Name[met_info$HMDB == metab]
  run_lme(exposure = "genotype", outcome=metab, 
          covariates = "genotype*time", coefficients_to_print = c(Genotype="genotype", Time="time"), 
          data=analysis2 %>% filter(time %in% c(0,120,235)), digits=c(1,3)) %>%
    mutate(outcome=metab, Outcome=Metabolite) %>%
    select(Exposure, outcome, Outcome, beta, se, p, anovaP, Effect, lowCI, upCI) %>%
    filter(grepl(" x ", Exposure)) %>%
    mutate_at("Exposure", ~gsub(".* x", "HC Genotype x", .))
  }) %>% do.call(rbind.data.frame, .)

## Save as .csv
#lme_metab_mmtt_unadj %>% fwrite("../data/processed/tab_result_lme_metab_mmtt_unadj.csv")
lme_metab_mmtt_unadj <- fread("../data/processed/tab_result_lme_metab_mmtt_unadj.csv")


## Run additional sensitivity analyses in demographics & bmi-adjusted models
models <- list(primary="age+sex+PC1z+PC2z+PC3z+genotype*time",
               bmi="age+sex+PC1z+PC2z+PC3z+genotype*time+bmi")

lme_metab_postprand_mmtt_adj <- do.call(rbind.data.frame, lapply(1:2, function(mod) {
  lapply(1:length(metabolites), function(m) {
    metab <- metabolites[m]
    Metabolite <- met_info$Name[met_info$HMDB == metab]
    run_lme(exposure = "genotype", outcome=metab, 
            covariates = models[[mod]], coefficients_to_print = coefs_to_print, 
            data=analysis2 %>% filter(time %in% c(0,120,235)), digits=c(1,3)) %>%
      mutate(outcome=metab, outcome=Metabolite) %>%
      mutate(model=names(models[mod]))
    }) %>% do.call(rbind.data.frame, .) 
  })) %>% 
  arrange(outcome, model) ; lme_metab_postprand_mmtt_adj

# Save as .csv
#lme_metab_postprand_mmtt_adj %>% fwrite("../data/processed/tab_result_lme_metab_postpran_mmtt_adj.csv")
lme_metab_postprand_mmtt_adj <- fread("../data/processed/tab_result_lme_metab_postpran_mmtt_adj.csv")


################################################################################
## Primary analysis stratified by GENOTYPE
################################################################################

lme_metab_mmtt_bygeno <- lapply(genotypes, function(geno) {
  lapply(1:length(metabolites), function(m) {
    metab <- metabolites[m]
    Metabolite <- met_info$Name[met_info$HMDB == metab]
    run_lme(exposure = "time", outcome=metab, 
            covariates = "time", coefficients_to_print = c(Time="time"), 
            data=analysis2 %>% filter(time %in% c(0,120,235) & genotype == geno), digits=c(1,3)) %>%
      mutate(outcome=metab, Outcome=Metabolite) %>%
      mutate(genotype = geno, .before=Exposure) %>%
      select(genotype, Exposure, outcome, Outcome, beta, se, p, anovaP, Effect, lowCI, upCI)
  }) %>% do.call(rbind.data.frame, .)
}) %>% do.call(rbind.data.frame, .)

## Save as .csv
#lme_metab_mmtt_bygeno %>% fwrite("../data/processed/tab_result_lme_metab_mmtt_unadj_bygeno.csv")
lme_metab_mmtt_bygeno <- fread("../data/processed/tab_result_lme_metab_mmtt_unadj_bygeno.csv")


# ==============================================================================
## Aggregate list of significant metabolites & fold changes, by genotype
# ===============================================================================

#dir.create("../data/processed/pathways")

##significant metabolites by time and genotype
genotypes=list(HC="HC genotype", HF="HF genotype")
times=c(m120="120 min", m235="235 min")

signif_metabs.l <- lapply(genotypes, function(g) {
  temp <-lapply(times, function(t) {
    lme_metab_mmtt_bygeno %>% rename(time=Exposure) %>%
      mutate_at("time", ~paste0(gsub("Time_", "", .), " min")) %>%
      filter(genotype==g & time == t) %>%
      arrange(p) %>%
      filter(p<0.05/23 & abs(beta/log(2))>0.5) %>% pull(outcome) %>% unique() }) ## NOTE: CONVERT BETA TO LOG2 FOLD CHANGE FOR THRESHOLD >0.5
  names(temp) <- names(times) ; temp
}) ; names(signif_metabs.l) <- names(genotypes)


## Extract metabolite list for common metabolites in ALL participants
dir.create("../data/processed/enrichment/metab_pEff_log2fc05")
lapply(c("120", "235"), function(time) {
  lme_metab_mmtt_timeonly %>% filter(Exposure == paste0("Time_",time) & 
                                     p<0.05/23 & abs(beta/log(2))>0.5) %>% pull(outcome) %>% length()
  #as.data.frame() %>% fwrite(paste0("../data/processed/enrichment/metab_pEff_log2fc05/all",time,"_pEff_log2fc05.txt"), col.names = F)
})

signif_metabs.l$HC$m120 %>% as.data.frame() %>% fwrite("../data/processed/enrichment/metab_pEff_log2fc05/hc120_pEff_log2fc05.txt")
signif_metabs.l$HC$m235 %>% as.data.frame() %>% fwrite("../data/processed/enrichment/metab_pEff_log2fc05/hc235_pEff_log2fc05.txt")
signif_metabs.l$HF$m120 %>% as.data.frame() %>% fwrite("../data/processed/enrichment/metab_pEff_log2fc05/hf120_pEff_log2fc05.txt")
signif_metabs.l$HF$m235 %>% as.data.frame() %>% fwrite("../data/processed/enrichment/metab_pEff_log2fc05/hf235_pEff_log2fc05.txt")

# Union
intersect(signif_metabs.l$HC$m120, signif_metabs.l$HF$m120) %>% length() # 104
intersect(signif_metabs.l$HC$m235, signif_metabs.l$HF$m235) %>% length() # 91

## Look at the metabolite types reflected in each genotype/time
metab_dict %>% filter(Input.name %in% signif_metabs.l$HC$m120) %>% reframe(Super=n_pct(Sub.class))
metab_dict %>% filter(Input.name %in% c(intersect(signif_metabs.l$HC$m120, signif_metabs.l$HF$m120) )) %>% reframe(Super=n_pct(Sub.class))


# ===============================================================================================
## After running Pathway Enrichment Analysis --> Find significant differences in BA metabolism
# ===============================================================================================

# top metabolites_geno for HC at 120 and NOT for HF at 120
temp <- lme_metab_mmtt_bygeno %>% 
  filter(Exposure != "Time_(joint)") %>%
  rename(Time=Exposure) %>%
  mutate_at("Time", ~paste0(gsub("Time_", "", .), " min")) %>%
  mutate_at("Time", ~factor(., levels=c("120 min", "235 min"))) %>%
  filter(outcome %in% signif_metabs.l$HC$m120[!signif_metabs.l$HC$m120 %in% signif_metabs.l$HF$m120] 
                                           & genotype == "HC genotype" & Time == "120 min") %>%
  arrange(p) %>% select(outcome, Outcome)

# Bile Acids changing significantly & with large effects
bileacids <- c("HMDB0000637"="Glycochenodeoxycholic acid", "HMDB0000631"="Glycodeoxycholic acid", 
               "HMDB0000896"="Taurodeoxycholic acid", "HMDB0000951"="Taurochenodeoxycholic acid")
# note: GCDCA, TCDCA, GCA and TCA identified in MetaboAnalyst as Primary Bile Acid enrichment 

## All bileacids (not just those changing most for HC)
bileacids_all <- c(temp$Outcome[c(7,9,11,13,27,30,33,35,38)]) ; 
names(bileacids_all) <- c(temp$outcome[c(7,9,11,13,27,30,33,35,38)])
bileacids_all[c(8:9)] <- gsub(".*or ","",bileacids_all[c(8,9)])
bileacids_all_abbrev <- c("TCDCA", "TDCA", "GDCA", "GCDCA", "GCA", "TCA", "GUDCA", "TUDCA", "TbMCA") ; 
names(bileacids_all_abbrev) <- names(bileacids_all)

## ALL BILE ACID METABOLITES
#bileacids_all <- metab_dict %>% filter(Main.class == "Bile acids") %>% pull(Standardized.name)
#names(bileacids_all) <- metab_dict %>% filter(Main.class == "Bile acids") %>% pull(Input.name)

lme_metab_mmtt_bygeno %>% filter(outcome %in% names(bileacids_all)) %>% 
  filter(Exposure %in% c("Time_120", "Time_235")) %>%
  filter(p<0.05) %>%
  arrange(p)

# ===========================================================
## Calculate 120 min fold change for primary bile acid metabolites
# ===========================================================

bileacids_fc.df <- lapply(names(bileacids_all), function(ba) {
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
  select(id, starts_with(clinical_metabolites), -starts_with("tg_log"),
         genotype, meal_choice, genetic_cat_x_meal_choice,
         age, sex, bmi, PC1z, PC2z, PC3z, -contains("iAUC")) %>% 
  pivot_longer(starts_with(clinical_metabolites), names_sep="_", names_to=c("metabolite", "time")) %>%
  filter(!is.na(time), !is.na(genotype)) %>%
  pivot_wider(names_from=metabolite) %>%
  mutate(time=factor(time, levels=c(0,30,60,120,180,235,270,300,360))) %>%
  arrange(time) ; head(postprandial2) 

postprandial2 <- left_join(postprandial2, analysis2 %>% select(id, time, names(bileacids_all)), by=c("id","time")) 


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
  postprandial2 %>% filter(time == 120) %>% select(id, age, sex, PC1z, PC2z, PC3z, names(bileacids_all),
                                                   "glucose", "insulin", "hdl", "ldl", "tg") %>%
    rename_at(c(names(bileacids_all), "glucose", "insulin", "hdl", "ldl", "tg"), ~paste0(., "_120")),
  # BA fold change from 0 to 120
  bileacids_fc.df %>% rename_at(names(bileacids_all), ~paste0(., "_120fc")), by = "id") %>% 
  #left_join(analysis %>% select(id, starts_with(c("glucose", "insulin")) & ends_with("net")), by="id") %>%
  # glucose and insulin 120 min iAUC
  left_join(analysis %>% select(id, c(paste0(paste(rep(c("glucose", "insulin"), each=7), c(30, 60, 120, 180, 235, 270, 360), sep="_"), rep("iAUC.net",7)))), by="id") %>%
  # glucose and insulin at 30 min
  left_join(postprandial %>% filter(time == 30) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_30")), by="id") %>%
  # glucose and insulin at 60 min
  left_join(postprandial %>% filter(time == 60) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_60")), by="id") %>%
  # glucose and insulin at 60 min
  left_join(postprandial %>% filter(time == 180) %>% select(id, "glucose", "insulin") %>% rename_at(c("glucose", "insulin"), ~paste0(., "_180")), by="id") %>%
  # clinical metabolites at 235
  left_join(postprandial2 %>% filter(time == 235) %>% select(id, names(bileacids_all), "glucose", "insulin", "tg", "ldl", "hdl") %>% 
              rename_at(c(names(bileacids_all), "glucose", "insulin", "tg", "ldl", "hdl"), ~paste0(., "_235")), by="id")  %>%
  # bile acid metabolites at 0
  left_join(postprandial2 %>% filter(time == 0) %>% select(id, names(bileacids_all)) %>% 
              rename_at(names(bileacids_all), ~paste0(., "_0")), by="id") %>%
  filter(complete.cases(genotype))

ba <- ba %>% 
  rowwise() %>%
  # Add total BA metabolites at 120 and 235
  mutate(totalBA_120 = rowSums(across(c(paste0(names(bileacids_all), "_120")))),
         totalBA_235 = rowSums(across(c(paste0(names(bileacids_all), "_235"))))) %>% 
  # Add total BA at 0
  left_join(postprandial2 %>% filter(time == 0) %>% select(id, names(bileacids_all)) %>%
              mutate(totalBA_0 = rowSums(across(names(bileacids_all)))) %>%
              select(id, totalBA_0), 
            by="id") %>%
  mutate(totalBA_120fc = totalBA_120 - totalBA_0)


# =================================================================
# Compile all exposures and outcomes to text with LINEAR models
# =================================================================

## Exposures: 120 iAUC, 120 min levels, 120 min fold change
exposures <- paste0(names(bileacids_all), rep("_120fc",8))
## Outcomes:  120 iAUC, 120 min levels, 120 min fold change
outcomes <- c(paste0(rep(clinvars, each=2), rep(c("_120", "_235"), length(clinvars))),
              paste0(rep(c("glucose", "insulin"), each=4), rep(c("_60", "_180", "_120iAUC.net", "_60iAUC.net"))))

# ===================================================================
## Summary table of bile acid levels at 0, 120 and 235 by genotype
# ===================================================================

lapply(c("0", "120", "235"), function(i) {
  vars <- bileacids_all ; names(vars) <- paste0(names(vars),"_",i)
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


# =========================================================================
## Pearson correlations of 120m FC bile acids with glucose/insulin/TG
# =========================================================================

yvars <- c(paste(rep(c("glucose", "insulin"), each=11), c("30","60","120","180","235","30iAUC.net","60iAUC.net","120iAUC.net","180iAUC.net","235iAUC.net"), sep="_"),
           "tg_120", "tg_235")

ba_clinvars_corr <- lapply(exposures, function(x) {
  X <- met_info$Name[met_info$HMDB==gsub("_.*", "", x)][1]
  lapply(c(yvars), function(y) {
  cor <- cor.test(ba %>% pull(x), ba %>% pull(y))
  cbind.data.frame(BA=X, ba=x, y=y, cor=cor$estimate, p=cor$p.value, lowci=cor$conf.int[[1]], upci=cor$conf.int[[2]])
  }) %>% do.call(rbind.data.frame, .) 
})  %>% do.call(rbind.data.frame, .) %>% unique()
ba_clinvars_corr
ba_clinvars_corr %>% arrange(p) %>% filter(endsWith(ba,"fc")) %>% filter(p<0.3)

#ba_clinvars_corr %>% fwrite("../data/processed/tab_result_metab_clinvars_corr_pearson.csv")


# =============================================
## All exposures & outcomes for LM
# =============================================

ba_clinvars_lm <- lapply(exposures, function(x) {
  X <- met_info$Name[met_info$HMDB==gsub("_.*", "", x)][1]
  lapply(c(outcomes, "glucose_30", "glucose_60", "glucose_180", "insulin_30", "insulin_60", "insulin_180"), function(y) {
    print_lm(exposure = x, outcome = y, covariates = "genotype+age+sex+PC1z+PC2z+PC3z",
             label=X, data=ba)
  }) %>% do.call(rbind.data.frame, .) 
  })  %>% do.call(rbind.data.frame, .) %>% unique()

ba_clinvars_lm %>% arrange(p) %>% filter(endsWith(exposure,"fc")) %>% filter(p<0.3)

ba_clinvars_lm %>% fwrite("../data/processed/tab_metab_ba_clinvars_lm_all.csv")
ba_clinvars_lm %>% arrange(p) %>% filter(p<0.05)

## Isolate changes
ba_clinvars_lm %>% 
  filter(endsWith(outcome, "net")) %>%
  filter(endsWith(exposure, "fc")) %>%
  arrange(p)


# ==========================================================================================
## Mediation analysis: test whether BA explain any of the genotype --> glucose/insulin assoc.
# ==========================================================================================

## CAN'T test in a LMER since we don't have BA data at 30/60/180 min --> need to use standard lm
summary(lmerTest::lmer(glucose~genotype*time+age+sex+PC1z+PC2z+PC3z+(1|id)+HMDB0000951, data=postprandial2 %>% filter(time != 360)))

## Test mediation framework by calculating Percent Direct/Mediated Effects
library("mediation")
detach("package:mediation")

## test for mediation of 120fc total BA for genotype on glucose/insulin AT 60 min
mA <- lm(totalBA_120fc~genotype+age+sex+PC1z+PC2z+PC3z, data=ba)
mB <- lm(glucose_60~genotype+totalBA_120fc+age+sex+PC1z+PC2z+PC3z, data=ba)
med1 <- mediate(mA, mB, treat="genotype", treat.value = "HF genotype", control.value = "HC genotype", mediator="totalBA_120fc", robustSE=F, sims=1000)
summary(med1)

mA <- lm(totalBA_120fc~genotype+age+sex+PC1z+PC2z+PC3z, data=ba)
mB <- lm(insulin_235~genotype+totalBA_120fc+age+sex+PC1z+PC2z+PC3z, data=ba)
med2 <- mediate(mA, mB, treat="genotype", treat.value = "HF genotype", control.value = "HC genotype", mediator="totalBA_120fc", robustSE=F, sims=1000)
summary(med2)

## --> no significant, or tend, of mediated effect of total BA on glucose/insulin at 60


# test for mediation of 120fc total BA for genotype on glucose/insulin 120 min AUC
mA <- lm(totalBA_120fc~genotype+age+sex+PC1z+PC2z+PC3z, data=ba)
mB <- lm(glucose_120iAUC.net~genotype+totalBA_120fc+age+sex+PC1z+PC2z+PC3z, data=ba)
med3 <- mediate(mA, mB, treat="genotype", treat.value = "HF genotype", control.value = "HC genotype", mediator="totalBA_120fc", robustSE=F, sims=1000)
summary(med3)

mA <- lm(totalBA_120fc~genotype+age+sex+PC1z+PC2z+PC3z, data=ba)
mB <- lm(insulin_120iAUC.net~genotype+totalBA_120fc+age+sex+PC1z+PC2z+PC3z, data=ba)
med4 <- mediate(mA, mB, treat="genotype", treat.value = "HF genotype", control.value = "HC genotype", mediator="totalBA_120fc", robustSE=F, sims=1000)
summary(med4)




#########################################################################################
## SELF-SELECTED MIXED MEALS: Run LM for change in metabolites to self-selected meals 
#########################################################################################

meals <- c("HC meal", "HF meal")
genotypes <- c("HC genotype", "HF genotype")


# ==============================================================
## Exploratory metabolomics analysis stratified by meal type - Time only
# ==============================================================

lme_metab_ssmt_mealonly <- lapply(meals, function(meal) {
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
lme_metab_ssmt_bymeal %>% fwrite("../data/processed/tab_result_lme_metab_ssmt_unadj_bymeal.csv")
lme_metab_ssmt_bymeal <- fread("../data/processed/tab_result_lme_metab_ssmt_unadj_bymeal.csv")


# ==============================================================
## Exploratory metabolomics analysis stratified by meal type
# ==============================================================

lme_metab_ssmt_bymeal <- lapply(meals, function(meal) {
    lapply(1:length(metabolites), function(m) {
      metab <- metabolites[m]
      Metabolite <- met_info$Name[met_info$HMDB == metab][1]
      run_lme(exposure = "genotype", outcome=metab, 
              covariates = "genotype*time",coefficients_to_print = c(Genotype="genotype", Time="time"), 
              data=analysis2 %>% filter(time %in% c(235,360) & meal_choice == meal), 
              digits=c(1,3)) %>%
        mutate(outcome=metab, Outcome=Metabolite) %>%
        mutate(meal_type = meal, genotype = "All", .before=Exposure) 
    }) %>% do.call(rbind.data.frame, .)
  }) %>% do.call(rbind.data.frame, .) %>%
  dplyr::select(genotype, meal_type, Exposure, outcome, Outcome, beta, se, p, anovaP, Effect, lowCI, upCI)

## Save as .csv
#lme_metab_ssmt_bymeal %>% fwrite("../data/processed/tab_result_lme_metab_ssmt_unadj_bymeal.csv")
lme_metab_ssmt_bymeal <- fread("../data/processed/tab_result_lme_metab_ssmt_unadj_bymeal.csv")


# ==============================================================================
## Exploratory metabolomics analysis stratified by meal type & genotype ********
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
  dplyr::select(genotype, meal_type, Exposure, outcome, Outcome, beta, se, p, anovaP, Effect, lowCI, upCI)


## Save as .csv
#lme_metab_ssmt_bymealxgeno %>% fwrite("../data/processed/tab_result_lme_metab_ssmt_unadj_bymealxgeno_2.csv")
lme_metab_ssmt_bymealxgeno <- fread("../data/processed/tab_result_lme_metab_ssmt_unadj_bymealxgeno_2.csv")


# ============================
## Calculate Kendall's Tao for agreement between first and second high(ER) carb meals
# ============================

tao.dat <- inner_join(
  lme_metab_ssmt_bymeal %>%
    dplyr::filter(Exposure == "Genotype_HC Genotype_ x Time_360" & meal_type == "HC meal") %>%
    arrange(beta) %>% mutate(ssmt = 1:nrow(.)) %>%
    dplyr::select(outcome, Outcome, ssmt), 
  lme_metab_mmtt_unadj %>%
    filter(Exposure == "HC Genotype x Time_120") %>%
    arrange(beta) %>% mutate(mmtt = 1:nrow(.)) %>%
    dplyr::select(outcome, Outcome, mmtt),
  by=c("Outcome", "outcome")) ; tao.dat

cor(tao.dat$mmtt, tao.dat$ssmt, method="kendall") # 0.285


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
        filter(p<0.05/23 & abs(beta/log(2))>0.5) %>% pull(outcome) %>% unique() }) ## NOTE: CONVERT BETA TO LOG2 FOLD CHANGE FOR THRESHOLD >0.5 **TRYING >0.3**
    names(temp) <- names(meals) ; temp
  }) ; names(ssmt_metabs.l) <- names(genotypes)


## Extract metabolite list for common metabolites in ALL participants
dir.create("../data/processed/enrichment/ssmt/metab_pEff_log2fc025")
lapply(c("HC meal", "HF meal"), function(meal) {
  lme_metab_ssmt_mealonly %>% 
    filter(Exposure == "Time_360" & meal_type == meal & p<0.05/23 & abs(beta/log(2))>0.5) %>% pull(outcome) %>%
  as.data.frame() %>% fwrite(paste0("../data/processed/enrichment/ssmt/metab_pEff_log2fc05/all",meal,"_pEff_log2fc05.txt"), col.names = F)
})

ssmt_metabs.l$HC$HC %>% as.data.frame() %>% fwrite("../data/processed/enrichment/ssmt/metab_pEff_log2fc05/hc_hc_pEff_log2fc05.txt")
ssmt_metabs.l$HC$HF %>% as.data.frame() %>% fwrite("../data/processed/enrichment/ssmt/metab_p0Eff_log2fc05/hc_hf_pEff_log2fc05.txt")
ssmt_metabs.l$HF$HC %>% as.data.frame() %>% fwrite("../data/processed/enrichment/ssmt/metab_pEff_log2fc05/hf_hc_pEff_log2fc05.txt")
ssmt_metabs.l$HF$HF %>% as.data.frame() %>% fwrite("../data/processed/enrichment/ssmt/metab_pEff_log2fc05/hf_hf_pEff_log2fc05.txt")


## --> RUN ENRICHMENT ANALYSIS
## See selective bile acid enrichment. for HF on HC meal, and HC genotype on HF meal

# BA metabolites changing for HC/HF genotype that are bile acids: same for both genotypes x meals
hchf_ba <- ssmt_metabs.l$HC$HF[ssmt_metabs.l$HC$HF %in% names(bileacids_all)]
hfhc_ba <- ssmt_metabs.l$HF$HC[ssmt_metabs.l$HF$HC %in% names(bileacids_all)]

# BA metabolites identified in enrichment analysis for primary bile acid synthesis
metab_dict %>% filter(Input.name %in% hchf_ba)
bileacids_ssmt <- c("HMDB0000518" = "Chenodeoxycholic acid", "HMDB0000637" = "Glycochenodeoxycholic acid", 
                    "HMDB0000951" = "Taurochenodeoxycholic acid", "HMDB0000138" = "Glycocholic acid", 
                    "HMDB0000036" = "Taurocholic acid")

# All BA metabolites changing significantly for either genotype (regardles of KEGG pathway)
bileacids_ssmt_abbrev <- c("CDCA", "GCDCA", "TCDCA", "GCA", "TCA")
names(bileacids_ssmt_abbrev) <- names(bileacids_ssmt)
  
## Exposures: 120 iAUC, 120 min levels, 120 min fold change
## Outcomes:  120 iAUC, 120 min levels, 120 min fold change
exposures <- paste0(rep(c(names(bileacids_ssmt)), each=3), rep(c("_235", "_360fc", "_360"), length(bileacids_ssmt)))
outcomes <- c(paste0(rep(clinvars, each=2), rep(c("_235", "_360"), length(clinvars))))

bileacids_fc_ssmt.df <- lapply(c(names(bileacids_all),"HMDB0000518","HMDB0000124", "HMDB0000258"), function(x) {
  analysis2 %>% select(id, genotype, time, meal_choice, names(bileacids_all), "HMDB0000518","HMDB0000124", "HMDB0000258") %>%
    filter(time %in% c(235, 360)) %>%
    pivot_longer(c(names(bileacids_all),"HMDB0000518","HMDB0000124", "HMDB0000258"), values_to="value", names_to="bileacids") %>%
    unite("ba_time", bileacids, time, sep="_") %>%
    pivot_wider(names_from=ba_time, values_from="value")  %>%
    
    select(id, mt1=paste0(x, "_235"), mt2=paste0(x, "_360"), genotype, meal_choice) %>%
    mutate(mdiff=mt2 - mt1) %>%
    rename_with(., ~gsub("mdiff", paste0(x, "_360fc"), .)) %>%
    #rename_with(., ~gsub("mdiff", paste0(x, "_", 235, ".", 360, "_diff"), .)) %>%
    select(id, ends_with("fc"), genotype, meal_choice)
  }) %>% reduce(full_join,by=c("id","genotype", "meal_choice")) ; bileacids_fc_ssmt.df

bileacids_ssmt.df <- analysis2 %>% select(id, genotype, time, names(bileacids_ssmt)) %>%
  filter(time %in% c(235, 360)) %>%
  pivot_longer(names(bileacids_ssmt), values_to="value", names_to="bileacids") %>%
  unite("ba_time", bileacids, time, sep="_") %>%
  pivot_wider(names_from=ba_time, values_from="value")  ; bileacids_ssmt.df

analysis4 <- left_join(analysis %>% select(
  id, genotype, meal_choice, genetic_cat_x_meal_choice, age, sex, PC1z, PC2z, PC3z, bmi,
  all_of(outcomes), glucose_270, glucose_300, insulin_270, insulin_300),
  bileacids_fc_ssmt.df, by=c("id", "genotype")) %>%
  left_join(bileacids_ssmt.df, by=c("id", "genotype")) %>%
  mutate_at("meal_choice", ~factor(., levels=c("HF meal", "HC meal")))

bileacids_fc_ssmt.df <- bileacids_fc_ssmt.df %>% left_join(analysis %>% select(id,  meal_choice, genexmeal = genetic_cat_x_meal_choice), by="id")


# ==============================================================================================================
## 1. do changes in bile acids from 235 to 360 correlate with glucose, insulin and lipid levels at 235 and 360
# ==============================================================================================================

postprandial3 <- left_join(postprandial2, analysis2 %>% select(id, time, "HMDB0000518","HMDB0000124", "HMDB0000258"), by=c("id","time"))

# ====================================================================================
## Prepare data with 120 fold change in BA, glucose/insulin and lipid metabolites
# ====================================================================================

ba2 <- left_join(ba, 
  # Bile Acids at 235 and 360
  postprandial3 %>% filter(time == 360) %>% select(id, age, sex, PC1z, PC2z, PC3z, names(bileacids_all), "HMDB0000518", "HMDB0000124", "HMDB0000258",
                                                   "glucose", "insulin", "hdl", "ldl", "tg") %>%
    rename_at(c(names(bileacids_all), "HMDB0000518", "HMDB0000124", "HMDB0000258", "glucose", "insulin", "hdl", "ldl", "tg"), ~paste0(., "_360")), by="id") %>%
  left_join(postprandial3 %>% filter(time==235) %>% select(id, "HMDB0000518",  "HMDB0000124", "HMDB0000258") %>% 
              rename_at(c("HMDB0000518", "HMDB0000124", "HMDB0000258"), ~paste0(., "_235")), by="id") %>%
  #left_join(postprandial3 %>% filter(time==360) %>% select(id, "HMDB0000518",  "HMDB0000124", "HMDB0000258") %>% 
   #           rename_at(c("HMDB0000518", "HMDB0000124", "HMDB0000258"), ~paste0(., "_360")), by="id") %>%
  # BA fold change from 235 to 360
  left_join(bileacids_fc_ssmt.df %>% select("id", paste0(c(names(bileacids_all), "HMDB0000518",  "HMDB0000124", "HMDB0000258"), "_360fc")), by = "id") %>% 
  # glucose and insulin 360 min iAUC
  left_join(analysis %>% select(id, c(paste0(rep(c("glucose", "insulin"), each=4), c("_270", "_300")) )), by="id") %>%
  filter(complete.cases(genotype))

ba2 <- ba2 %>% 
  rowwise() %>%
  # Add total BA metabolites at 120 and 235
  mutate(totalBA_235 = rowSums(across(c(paste0(c(names(bileacids_all),"HMDB0000518"), "_235")))),
         totalBA_360 = rowSums(across(c(paste0(c(names(bileacids_all),"HMDB0000518"), "_360")))) ) %>%
  mutate(totalBA_360fc = totalBA_360 - totalBA_235)

yvars <- c(paste(rep(c("glucose", "insulin"), each=11), c("235","270","300","360", "270"), sep="_"),
           "tg_235", "tg_360")

## Run simple correlations 
exposures <- paste0(rep(c(names(bileacids_all),"HMDB0000518", "HMDB0000124", "HMDB0000258"), each=3), 
                    rep(c("_235", "_360","_360fc"), 12))

ba_clinvars_corr_2 <- lapply(exposures, function(x) {
  X <- met_info$Name[met_info$HMDB==gsub("_.*", "", x)][1]
  lapply(c(yvars), function(y) {
    cor <- cor.test(ba2 %>% pull(x), ba2 %>% pull(y))
    cbind.data.frame(BA=X, ba=x, y=y, cor=cor$estimate, p=cor$p.value, lowci=cor$conf.int[[1]], upci=cor$conf.int[[2]])
  }) %>% do.call(rbind.data.frame, .) 
})  %>% do.call(rbind.data.frame, .) %>% unique()
ba_clinvars_corr_2
ba_clinvars_corr_2 %>% arrange(p) #%>% filter(endsWith(ba,"fc")) %>% filter(p<0.3)

#ba_clinvars_corr %>% fwrite("../data/processed/tab_result_metab_clinvars_corr_pearson.csv")


## Check correlations of carbohydrte metabolites with glucose levels 
# HMDB0000258 = sucrose ; HMDB0000124 = F 6-P
cor.test(analysis2 %>% filter(time==360) %>% pull("HMDB0000258"), analysis$glucose_270iAUC.pos)



ba_clinvars_lm_ssmt <- lapply(exposures, function(x) {
    X <- bileacids_all[[which(names(bileacids_all) == gsub("_.*", "", x))]]
    lapply(c(outcomes,"glucose_270", "glucose_300", "insulin_270", "insulin_300"), function(y) {
      print_lm(exposure = x, outcome = y, covariates = "meal_choice+age+sex+PC1z+PC2z+PC3z",
               label=X, data=analysis4)
    }) %>% do.call(rbind.data.frame, .) 
}) %>% do.call(rbind.data.frame, .) ; ba_clinvars_lm_ssmt

ba_clinvars_lm_ssmt %>% arrange(p) %>% filter(endsWith(exposure, "fc")) %>% filter(p<0.2) 



## fold change outcome in bile acids with genotype * meal type interaction
lapply(names(bileacids_all), function(y) {
  Y <- bileacids_all[[which(names(bileacids_all) == gsub("_.*", "", y))]]
  y <- paste0(y, "_360fc")
  summary(lm(formula(paste0(y,"~genotype*meal_choice+age+sex+PC1z+PC2z+PC3z")), data=analysis4))$coef %>%
    as.data.frame() %>%
    filter(grepl(":",rownames(.))) %>%
    mutate(bileacid=y, BileAcid=Y)
  }) %>% do.call(rbind.data.frame, .)

# stratified by gentoype
lapply(genotypes, function(g) {
  lapply(names(bileacids_all), function(y) {
    Y <- bileacids_all[[which(names(bileacids_all) == gsub("_.*", "", y))]]
    y <- paste0(y, "_360fc")
    summary(lm(formula(paste0(y,"~meal_choice+age+sex+PC1z+PC2z+PC3z")), data=analysis4 %>% filter(genotype == g)))$coef %>%
      as.data.frame() %>%
      filter(grepl("meal_choice",rownames(.))) %>%
      mutate(bileacid=y, BileAcid=Y, genotyp=g)
  }) %>% do.call(rbind.data.frame, .)
}) %>% do.call(rbind.data.frame, .)


# =============================================================
## LME for metabolites at 235/360 - stratified by meal type
# =============================================================

lme_ba_ssmt_bymeal <- lapply(meals, function(meal) { 
  lapply(1:length(bileacids_all), function(m) {
    ba <- names(bileacids_all)[m]
    BA <- bileacids_all[[m]]
    run_lme(exposure = "genotype*time", outcome=ba, covariates = "time", 
            coefficients_to_print = c("genotype", "time"), 
            data=analysis2 %>% filter(time %in% c(235, 360) & meal_choice == meal), digits=c(1,3)) %>%
      mutate(outcome=ba, Outcome=BA, meal_type=meal) %>%
      filter(grepl(" x ", Exposure)) %>%
      select(Exposure, outcome, Outcome, meal_type, beta, se, p, anovaP, Effect, lowCI, upCI)
  }) %>% do.call(rbind.data.frame, .)
  }) %>% do.call(rbind.data.frame, .) ; lme_ba_ssmt_bymeal

lme_metab_ssmt_bymeal %>% arrange(p)

## Save as .csv
lme_metab_ssmt_bymeal %>% 
  mutate_at("Exposure", ~gsub(".* x", "HF meal x", .)) %>%
  fwrite("../data/processed/tab_result_lme_metab_ssmt_bymeal.csv")



## EOF

