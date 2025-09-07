# Statistical comparisons

library(tidyverse)
library(data.table)
library(lme4)
library(purrr)


units <- c("", "_mmol")
Units <- c("", "mmol/L")

unit<-units[1]
Unit<-Units[1]

lapply(list.files("../../pantry/functions/", recursive = T, full.names = T), source)

# Analysis dataset
analysis <- readRDS("../data/processed/premier_analysis.rda") %>%
  mutate(gene_x_meal = genetic_cat_x_meal_choice)


################################################################################
##  ---------------------  Build Variable Groups ---------------------------  ##
################################################################################

## Table 1 vars --------------------
 
# Demographics
demographic_vars <- c(
  age="Age, years", Sex="Sex", Race2="Race", Ethnicity="Ethnicity", bmi="BMI, kg/m2", 
  whr="Waist-Hip", sbp="SBP, mmHg", dbp="DBP, mmHg", 
  hba1c="HbA1c", homair="HOMA-IR", homab="HOMA-B")

# Fasting metabolites
metabolite_vars <- c(
  glucose="Glucose, mg/dL", insulin="Insulin", tg_log="log(Triglyceride)", tg="Triglyceride",
  chol="Total cholesterol", hdl="HDL cholesterol", ldl="LDL cholesterol", protein="Total protein")

# Demographics
survey_vars <- c(
  education="Education level", smoke_packyears="Smoking Status", alcohol_perweek_4lvl="Alcohol Use", 
  sleep_hrday="Sleep, hrs/day", physact_mvpa_perweek="Physical Activity, hr/week")


## Main outcomes --------------------

# Postprandial glycemic traits
postprandial_vars <- c(paste0(
  rep(c("Glucose ", "Insulin "), each=5), rep(c("0", "30", "60", "120", "180", "235"), 2), rep(" min", 10)))
names(postprandial_vars) <- c(paste0(
  rep(c("glucose", "insulin"), each=5), rep(c("", "_30", "_60", "_120", "_180", "_235"), 2))) 

# Second meal glycemic traits
postprandial_mmtt_vars <- c(paste0(
  rep(c("Glucose ", "Insulin "), each=3), rep(c("270", "300", "360"), 2), rep(" min", 6)))
names(postprandial_mmtt_vars) <- c(paste0(
  rep(c("glucose", "insulin"), each=3), rep(c("_270", "_300", "_360"), 2))) 

# iAUC glycemic traits
iAUC_vars <- c(paste0(
  rep(c("Glucose ", "Insulin "), each=8), rep(c("30", "60", "120", "180", "235", "270", "300", "360"), 2), rep("min iAUC (pos)",14), rep(paste0(" ", Unit), 10)))
names(iAUC_vars) <- c(paste0(
  rep(c("glucose_", "insulin_"), each=8), rep(c("30", "60", "120", "180", "235", "270", "300", "360"), 2),  rep("iAUC.pos",14)))



################################################################################
##  -------------------  Descriptive statistics ----------------------------  ##
################################################################################

## Set lists of variables & strata -----------------

summary_table_vars.l <- list(
  list(label="demographics", vars=demographic_vars), 
  list(label="metabolites", vars=metabolite_vars),
  list(label="survey", vars=survey_vars),
  list(label=paste0("postprandial",unit), vars=c(postprandial_vars, postprandial_mmtt_vars)),
  list(label=paste0("iAUC",unit), vars=iAUC_vars))

summary_table_strata.l <- list(
  list(label="genecat", var="genotype",order=c("HC genotype", "HF genotype")),
  list(label="mealchoice", var="meal_choice", order=c("HC meal", "HF meal")))


## Generate summary tables ----------------------

lapply(summary_table_strata.l, function(strata){
  lapply(summary_table_vars.l, function(vars) {
    print_summary_table(
      analysis, vars_to_summarize = vars$vars, p_adjust = "none", digits = c(1,1,3),
      var_strata = strata$var, var_strata_order = strata$order,
      p_types = "descriptive", p_smalln=T) %>%
      fwrite(paste0("../output/tab_descr_", vars$label, "_", strata$label, ".csv"), row.names = T)
  })
})

# Genotype effect on meal choice
chisq.test(analysis$genotype, analysis$meal_choice)
summary(glm(hc_meal~genotype+age+sex+PC1z+PC2z+PC3z, family=binomial("logit"), 
            data=analysis %>% mutate(hc_meal = ifelse(meal_choice=="HC meal",1,0))))

# Check meal choice effect effect on insulin levels at start of SSMT
t.test(insulin_235~meal_choice, data=analysis)
summary(lm(insulin_235~meal_choice+genotype, data=analysis))

summary(glm(hc_meal~insulin_235+genotype, family=binomial("logit"), 
            data=analysis %>% mutate(hc_meal = ifelse(meal_choice=="HC meal",1,0))))

## Mean differences in iAUC by genotype
summary(lm(glucose_120iAUC.net~genotype,data=analysis))
t.test(glucose_120iAUC.net~genotype,data=analysis)

## Additional comparisons across genotype x meal type 
print_summary_table(analysis, vars_to_summarize = summary_table_vars.l[[5]]$vars,
                    p_adjust = "none", digits = c(1,1,3), var_strata = "gene_x_meal", 
                    var_strata_order = c(unique(analysis$gene_x_meal)), p_print = F) %>%
  fwrite("../output/tab_descr_iAUC_genexmeal.csv", row.names = T)

anova(lm(glucose_30iAUC.pos~genotype, data=analysis))
t.test(glucose_60iAUC.pos~genotype,data=analysis)

################################################################################
##  ------------------ Postprandial glycemic responses  --------------------- ##
################################################################################

## Create a long dataframe, by time: postprandial
pp_vars <- c("0", "30", "60", "120", "180", "235", "270", "300", "360")
iAUC.pos_vars <- c("0iAUC.pos", "30iAUC.pos", "60iAUC.pos", "120iAUC.pos", 
                   "180iAUC.pos", "235iAUC.pos", "270iAUC.pos", "300iAUC.pos", "360iAUC.pos")

postprandial <- analysis %>% dplyr::select(id, starts_with("glucose"), starts_with("insulin"), 
                                    genotype, meal_choice, genetic_cat_x_meal_choice,
                                    age, sex, bmi, PC1z, PC2z, PC3z,
                                    smoke_packyears, alcohol_perweek_4lvl, physact_mvpa_perweek, sleep_hrday) %>%
  rename(glucose.mmtt_0=glucose, insulin.mmtt_0=insulin) %>%
  mutate(glucose.selected_0=glucose_235, insulin.selected_0=insulin_235) %>%
  pivot_longer(c(starts_with("glucose_"), starts_with("insulin_")),
               names_sep="_", names_to=c("metabolite", "time")) %>%
  pivot_wider(names_from=metabolite) %>%
  mutate(measure=ifelse(time %in% pp_vars, "pp", ifelse(time %in% iAUC.pos_vars, "iAUC.pos", NA))) %>%
  mutate(time=factor(time, levels=c(0,30,60,120,180,235,270,300,360))) %>%
  arrange(time) ; head(postprandial)

saveRDS(postprandial, paste0("../data/processed/postprandial_long",unit,".rda"))


############################################################################# 
##  Mixed-Meal Tolerance Tests (time 0-235)
############################################################################# 

# List metabolites & meal types
metabolites <- c("glucose", "insulin")
meal_types <- c("mmtt", "selected")
measures <- c("pp", "iAUC.pos")
coefs_to_print = c("Genotype"="genotype","Meal Choice"="meal_choice")

# Model covariates 
models_mmtt.l <- list(
  empty = list(main="genotype+time", int="genotype*time"), 
  primary=list(main="age+sex+PC1z+PC2z+PC3z+time", 
               int="age+sex+PC1z+PC2z+PC3z+time+genotype*time"),
  bmi=list(main="age+sex+PC1z+PC2z+PC3z+time+bmi", 
           int="age+sex+PC1z+PC2z+PC3z+time+genotype*time+bmi")
)

models_iAUC.pos.l <- list(
  empty="genotype",
  primary="age+sex+PC1z+PC2z+PC3z",
  bmi="age+sex+PC1z+PC2z+PC3z+bmi")


# =========================================================
## Levene's test for heterogeneity in response patterns
# =========================================================

library(car)
leveneTest(glucose ~ as.factor(time), data=postprandial %>% filter(time %in% c(0,30,60,120,180,235)))
leveneTest(insulin ~ as.factor(time), data=postprandial %>% filter(time %in% c(0,30,60,120,180,235)))


## Postprandial metabolites ----------------------------------

do.call(rbind, lapply(1:length(models_mmtt.l), function(m) {
  do.call(rbind, lapply(c("glucose", "insulin"), function(metab) {
    # subset data_long
    data_long <- postprandial %>% filter(time %in% c(0,30,60,120,180,235)) 
    # run lme models
    rbind(run_lme(exposure="genotype", outcome=metab, outcome_label = str_to_sentence(metab),
                  covariates=models_mmtt.l[[m]]$main, coefficients_to_print = coefs_to_print,
                  data_long=data_long, digits = c(1,3)),
          run_lme(exposure="genotype", outcome=metab, outcome_label = str_to_sentence(metab),
                  covariates=models_mmtt.l[[m]]$int, coefficients_to_print = coefs_to_print, 
                  data_long=data_long, digits = c(1,3)) %>%
            filter(grepl(" x ", Exposure))) 
  })) %>% mutate(Model=names(models_mmtt.l[m]), .before="Exposure")
  })) %>% 
  mutate(Meal.Type="MMTT", .before="Model") %>%
  select(Meal.Type, Model, Effect, Outcome, Exposure, beta, se, p, anovaF, anovaP, Beta.SE, P, P_signif, anovaP_signif) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  fwrite(paste0("../output/tab_res_mmtt_postprandial_mgdL.csv"))


## iAUCs ---------------------------------------------------

# primary model & adjusting for metabolite
do.call(rbind, lapply(1:length(models_iAUC.pos.l), function(m) {
  do.call(rbind, lapply(c("glucose", "insulin"), function(metab) {
    do.call(rbind, lapply(c(30,60,120,180,235), function(t) {
      y <- paste0(metab,"_",t,"iAUC.pos") ; 
      print_lm(exposure = "genotype", outcome=y, lm_trend = F,
               covariates = models_iAUC.pos.l[[m]], label=names(models_iAUC.pos.l[m]))
      }))
    }))
  })) %>% 
      mutate(Meal.Type="MMTT", Effect="Main") %>%
      mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.))) %>%
      mutate(lowCI=beta-1.96*se, upCI=beta+1.96*se) %>%
      mutate(Beta.SE=ifelse(is.na(beta)==T, "-", sprintf("%s (%s, %s)", round(beta,2), round(lowCI,2), round(upCI,2))), P=format_p(p)) %>%
      mutate(P_signif=format_p_star(p), anovaP_signif=format_p_star(f_p)) %>%
      select(Meal.Type, Model=model, Effect, Outcome=outcome, Exposure=exposure, N=n, beta, se, p, anovaF=f, anovaP=f_p, Beta.SE, P, P_signif, anovaP_signif) %>%
      mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
      mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
      mutate_all(., ~ifelse(is.na(.)==T | .=="NA (NA, NA)", "-", .)) %>%
      fwrite(paste0("../output/tab_res_mmtt_iAUCpos.csv"))



##############################################################################
## Self-Selected Meals (time 235-360)
##############################################################################

selected_strata.l <- list(meal=c("HC meal", "HF meal"), geno=c("HC genotype", "HF genotype"))
models_selected.l <- list(
  strat_meal=list(
    empty=list(main="genotype+time", int="genotype*time"),
    primary=list(main="age+sex+PC1z+PC2z+PC3z+genotype+time", 
                 int="age+sex+PC1z+PC2z+PC3z+genotype*time"),
    bmi=list(main="age+sex+PC1z+PC2z+PC3z+time+genotype+bmi", 
             int="age+sex+PC1z+PC2z+PC3z+genotype*time+bmi")),
  strat_geno=list(
    empty=list(main="meal_choice+time", int="meal_choice*time"),
    primary=list(main="age+sex+PC1z+PC2z+PC3z+meal_choice+time",
                 int="age+sex+PC1z+PC2z+PC3z+meal_choice*time"),
    bmi=list(main="age+sex+PC1z+PC2z+PC3z+meal_choice+time+bmi",
             int="age+sex+PC1z+PC2z+PC3z+meal_choice*time+bmi"))
)
 

## Overall meal effect on postprandial responses (all participants) ----------------------------
table(postprandial$meal_choice)

# meal choice x time
do.call(rbind, lapply(c("glucose", "insulin"), function(m) {
  # set model covariates & subset data_long
  models <- list(main=paste0("time"), 
                 int=paste0("meal_choice*time"))
  data_long <- postprandial %>% filter(time %in% c(235,270,300,360)) 
  # run lme models
  rbind(run_lme(exposure="meal_choice", outcome=m, outcome_label = str_to_sentence(m),
                covariates=paste0(models$main,"+age+sex"), coefficients_to_print = c(coefs_to_print, "Time"="time"), 
                data_long = data_long, digits = 1),
        run_lme(exposure="meal_choice", outcome=m, outcome_label = str_to_sentence(m),
                covariates=paste0(models$int,"+age+sex"), coefficients_to_print = c(coefs_to_print, "Time"="time"), 
                data_long = data_long,  digits = c(1,3)) %>%
          filter(grepl(" x ", Exposure))) 
  })) %>% mutate(Model="primary", .before="Exposure") %>% 
  mutate_all(., ~ifelse(is.na(.)==T, "-", .))  %>% 
  mutate(Meal.Type="Self-Selected", .before="Model") %>%
  select(Meal.Type, Model, Effect, Outcome, Exposure, beta, se, p, anovaF, anovaP, Beta.SE, P, P_signif, anovaP_signif) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  fwrite(paste0("../output/tab_res_ssmt_postprandial_mealEffect_main_mgdL.csv"))


## Stratified analyses: effect of genotype on postprandial responses by meal ----------------------------

# Set reference group**
postprandial$genotype <- relevel(postprandial$genotype, ref="HC genotype")

# stratified by meal type
do.call(rbind, lapply(1:length(models_selected.l$strat_meal), function(m) {
  do.call(rbind, lapply(c("glucose", "insulin"), function(metab) {
  data_long <- postprandial %>% filter(time %in% c(235,270,300,360)) 
  do.call(rbind, lapply(c("HC meal", "HF meal"), function(meal) {
    rbind(run_lme(exposure="genotype", outcome=metab, outcome_label = str_to_sentence(metab),
                  covariates=models_selected.l$strat_meal[[m]]$main, coefficients_to_print = coefs_to_print, 
                  data_long = data_long %>% filter(meal_choice == meal), digits=c(1,3)),
          run_lme(exposure="genotype", outcome=metab, outcome_label = str_to_sentence(metab),
                  covariates=models_selected.l$strat_meal[[m]]$int, coefficients_to_print = coefs_to_print, 
                  data_long = data_long %>% filter(meal_choice == meal), digits=c(1,3)) %>% 
            filter(grepl(" x ", Exposure))) %>%
      mutate(Level=meal) %>%
      mutate(Strata="Meal Choice", strata_level=meal, .before=Exposure)
      })) %>% mutate(Model=names(models_selected.l$strat_meal[m]), .before="Strata")
    }))
  })) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .))  %>% 
  mutate(Meal.Type="Self-Selected", .before="Model") %>%
  dplyr::select(Meal.Type, Model, Strata, Level, Effect, Outcome, Exposure, beta, se, p, anovaF, anovaP, Beta.SE, P, P_signif, anovaP_signif) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  fwrite(paste0("../output/tab_res_ssmt_postprandial_genoEffect_byMeal_refHF_mgdL.csv"))


# stratified by genotype
do.call(rbind, lapply(1:length(models_selected.l$strat_geno), function(m) {
  do.call(rbind, lapply(metabolites, function(metab) {
    data_long <- postprandial %>% filter(time %in% c(235,270,300,360)) 
    do.call(rbind, lapply(c("HC genotype", "HF genotype"), function(geno) {
      rbind(run_lme(exposure="meal_choice", outcome=metab, outcome_label = str_to_sentence(metab), digits=c(1,3),
                    covariates=models_selected.l$strat_geno[[m]]$main, coefficients_to_print = coefs_to_print, 
                    data_long = data_long %>% filter(genotype == geno)),
            run_lme(exposure="meal_choice", outcome=metab, outcome_label = str_to_sentence(metab), digits=c(1,3),
                    covariates=models_selected.l$strat_geno[[m]]$int, coefficients_to_print = coefs_to_print, 
                    data_long = data_long %>% filter(genotype == geno)) %>% 
              filter(grepl(" x ", Exposure))) %>%
        mutate(Level=geno) %>%
        mutate(Strata="Genotype", strata_level=geno, .before=Exposure)
    })) %>% mutate(Model=names(models_selected.l$strat_geno[m]), .before="Strata")
  }))
})) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .))  %>% 
  mutate(Meal.Type="Self-Selected", .before="Model") %>%
  select(Meal.Type, Model, Strata, Level, Effect, Outcome, Exposure, beta, se, p, anovaF, anovaP, Beta.SE, P, P_signif, anovaP_signif) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  fwrite(paste0("../output/tab_res_ssmt_postprandial_mealEffect_byGeno_refHF_mgdL.csv"))


## Effect of meal choice and genotype on iAUCs (stratified samples) ----------------------------------

print_summary_table(data=analysis, vars_to_summarize = c(
  glucose_270iAUC.pos="Glucose, 270 iAUC", glucose_300iAUC.pos="Glucose, 300 iAUC", glucose_360iAUC.pos="Glucose, 360 iAUC",
  insulin_270iAUC.pos="Insulin, 270 iAUC", insulin_300iAUC.pos="Insulin, 300 iAUC", insulin_360iAUC.pos="Insulin, 360 iAUC"), 
  var_strata="genetic_cat_x_meal_choice") %>%
  fwrite("../output/tab_descr_iAUC_genecat.csv")

# Glucose differences for SSMT
t.test(glucose_270iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal"))$p.value #P=0.947
t.test(glucose_270iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal"))$p.value #P=0.0421*
t.test(glucose_300iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal"))$p.value #P=0.721
t.test(glucose_300iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal"))$p.value #P=0.014*
t.test(glucose_360iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal"))$p.value #P=0.537 !
t.test(glucose_360iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal"))$p.value #P=0.03596* !

t.test(glucose_270iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HC genotype")) #P=0.04385
t.test(glucose_270iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HF genotype")) #P=0.1855
t.test(glucose_300iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HC genotype")) #P=0.0154
t.test(glucose_300iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HF genotype")) #P=0.1018

# Insulin differences for SSMT
t.test(insulin_270iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal")) #P=0.1791
t.test(insulin_270iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal")) #P=0.3778
t.test(insulin_300iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal")) #P=0.1445
t.test(insulin_300iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal")) #P=0.3387
t.test(insulin_360iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal")) #P=0.08116 !
t.test(insulin_360iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal")) #P=0.455 ! 

t.test(insulin_270iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HC genotype")) #P=0.3577
t.test(insulin_270iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HF genotype")) #P=0.1914
t.test(insulin_300iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HC genotype")) #P=0.2764
t.test(insulin_300iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HF genotype")) #P=0.1254

# 120-min Glucose & Insulin for SSMT stratified by Meal 
t.test(glucose_360iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HC genotype")) #P=0.003868**
t.test(glucose_360iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HF genotype")) #P=0.1224
t.test(insulin_360iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HC genotype")) #P=0.2653
t.test(insulin_360iAUC.pos~meal_choice, data=analysis %>% filter(genotype=="HF genotype")) #P=0.04498*

# Global interactions
anova(lm(glucose_360iAUC.pos~meal_choice*genotype, data=analysis)) 
anova(lm(insulin_360iAUC.pos~meal_choice*genotype, data=analysis)) 


## Stratified analyses: effect of genotype on iAUC by meal type ==========

do.call(rbind, lapply(metabolites, function(m) {
  do.call(rbind, lapply(c("HC meal", "HF meal"), function(meal) {
    do.call(rbind, lapply(c(270,300,360), function(t) {
      y <- paste0(m,"_",t,"iAUC.pos") ; 
      rbind.data.frame(
        # primary model
        print_lm(exposure = "genotype", outcome=y, lm_trend = F,
                 covariates = "age+sex+PC1z+PC2z+PC3z", label="Primary",
                 data = analysis %>% filter(meal_choice==meal)),
        # adjusting for metabolite at time0
        print_lm(exposure = "genotype", outcome=y, lm_trend = F,
                 covariates = paste0("age+sex+PC1z+PC2z+PC3z+",m,"_235"), label="AdjTime0",
                 data = analysis %>% filter(meal_choice==meal)),
        # adjusting bmi
        print_lm(exposure = "genotype", outcome=y, lm_trend = F,
                 covariates = paste0("age+sex+PC1z+PC2z+PC3z+bmi"), label="BMI",
        data = analysis %>% filter(meal_choice==meal))
        ) }))  %>% 
      mutate(strata="Meal Choice", strata_level=meal, .before=outcome) 
  }))
  })) %>%
  mutate(Meal.Type="Self-Selected", Effect="Main") %>%
  mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.))) %>%
  mutate(lowCI=beta-1.96*se, upCI=beta+1.96*se) %>%
  mutate(Beta.SE=ifelse(is.na(beta)==T, "-", sprintf("%s (%s, %s)", round(beta,1), round(lowCI,1), round(upCI,1))), P=format_p(p)) %>%
  mutate(P_signif=format_p_star(p), anovaP_signif=format_p_star(f_p)) %>%
  select(Strata=strata, Level=strata_level, Meal.Type, Model=model, Effect, Outcome=outcome, Exposure=exposure, N=n, beta, se, p, anovaF=f, anovaP=f_p, Beta.SE, P, P_signif, anovaP_signif) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  arrange(desc(Model), Level) %>% 
  fwrite(paste0("../output/tab_res_ssmt_iAUC_genoEffect_byMeal_mgdL.csv"))

## Stratified analyses: effect of meal type on iAUC by genotype ==========

do.call(rbind, lapply(metabolites, function(m) {
  do.call(rbind, lapply(c("HC genotype", "HF genotype"), function(geno) {
    do.call(rbind, lapply(c(270,300,360), function(t) {
      y <- paste0(m,"_",t,"iAUC.pos") ; 
      rbind.data.frame(
        # primary model
        print_lm(exposure = "meal_choice", outcome=y, lm_trend = F, digits=c(1,3),
                 covariates = "age+sex+PC1z+PC2z+PC3z", label="Primary",
                 data = analysis %>% filter(genotype==geno)),
        # adjusting for metabolite at time0
        print_lm(exposure = "meal_choice", outcome=y, lm_trend = F, digits=c(1,3),
                 covariates = paste0("age+sex+PC1z+PC2z+PC3z+",m,"_235"), label="AdjTime0",
                 data = analysis %>% filter(genotype==geno)),
        # adjusting bmi
        print_lm(exposure = "meal_choice", outcome=y, lm_trend = F, digits=c(1,3),
                 covariates = paste0("age+sex+PC1z+PC2z+PC3z+bmi"), label="BMI",
                 data = analysis %>% filter(genotype==geno))
      ) }))  %>% 
      mutate(strata="Genotype", strata_level=geno, .before=outcome) 
  }))
})) %>%
  mutate(Meal.Type="Self-Selected", Effect="Main") %>%
  mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.))) %>%
  mutate(lowCI=beta-1.96*se, upCI=beta+1.96*se) %>%
  mutate(Beta.SE=ifelse(is.na(beta)==T, "-", sprintf("%s (%s, %s)", round(beta,1), round(lowCI,1), round(upCI,1))), P=format_p(p)) %>%
  mutate(P_signif=format_p_star(p), anovaP_signif=format_p_star(f_p)) %>%
  select(Strata=strata, Level=strata_level, Meal.Type, Model=model, Effect, Outcome=outcome, Exposure=exposure, N=n, beta, se, p, anovaF=f, anovaP=f_p, Beta.SE, P, P_signif, anovaP_signif) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  arrange(desc(Model), Level) %>% 
  fwrite(paste0("../output/tab_res_ssmt_iAUC_mealEffect_byGeno_mgdL.csv"))


##################################################################
##  ---------- Sensitivity analyses  (Self-Selected) ---------- ##
##################################################################

## Additionally adjust for energy/carb/protein consumed 

# ==========
## MMTT
# ==========

# Check for differenes in energy/carb/fat/prot intake by genetic category

analysis %>% group_by(genotype) %>%
  reframe(Energy_MMTT=mean_sd(mmtt_energy_kcal),
          Carb_g_MMTT=mean_sd(mmtt_carb_g),
          Fat_g_MMTT=mean_sd(mmtt_fat_g),
          Prot_g_MMTT=mean_sd(mmtt_prot_g))

mmtt_meal_data <- c("mmtt_energy_kcal", "mmtt_carb_g", "mmtt_fat_g", "mmtt_prot_g")
lapply(mmtt_meal_data, function(y){
  summary(lm(formula(paste0(y, "~genotype++age+sex+PC1z+PC2z+PC3z")), data=analysis))
}) # No significant associations between genetic category with energy/carb/fat/protein in MMTTs


# ====================
## Self-Selected
# ====================

# Check for differenes in energy/carb/fat/prot intake by genetic category

analysis %>% group_by(genetic_cat_x_meal_choice) %>%
  reframe(Energy_selected=mean_sd(selected_energy_kcal),
          Carb_g_selected=mean_sd(selected_carb_g),
          Fat_g_selected=mean_sd(selected_fat_g),
          Prot_g_selected=mean_sd(selected_prot_g))

selected_meal_data <- c("selected_energy_kcal", "selected_carb_g", "selected_fat_g", "selected_prot_g")

mmtt_meal_vars <- c(mmtt_energy_kcal="Energy (kcal)", mmtt_carb_g="Carbohydrate (g)",
                    mmtt_fat_g="Fat (g)", mmtt_prot_g="Protein (g)")
print_summary_table(data=analysis, vars_to_summarize = mmtt_meal_vars, 
                    var_strata = "genotype", var_strata_order = c("HC genotype", "HF genotype"))

lapply(selected_meal_data, function(y){
  summary(lm(formula(paste0(y, "~genotype+meal_choice+age+sex+PC1z+PC2z+PC3z")), data=analysis))
}) # -significant ANOVA for genotype effects on selected_carb_g and selected_fat_g, after adjusting for meal_choice

lapply(selected_meal_data, function(y){
  anova(lm(formula(paste0(y, "~genotype*meal_choice+age+sex+bmi+PC1z+PC2z+PC3z")), data=analysis))
}) # -significant ANOVA for genotype*meal_choice interaction on selected_carb_g and selected_fat_g


# =================================================================
## Additionally adjust for demographic & lifestyle covariates
# =================================================================

m_lifestyle_ppMain <- paste0(m_primary_ppMain, "+smoke_packyears+alcohol_perweek_4lvl+physact_mvpa_perweek")
m_lifestyl_ppInt <- paste0(m_primary_ppInt, "+smoke_packyears+alcohol_perweek_4lvl+physact_mvpa_perweek")

do.call(rbind, lapply(metabolites, function(m) {
  data_long <- postprandial %>% filter(time %in% c(0,30,60,120,180)) 
  bind_rows(run_lme(exposure="genotype", outcome=m, outcome_label = str_to_sentence(m),
                    covariates=m_lifestyle_ppMain, coefficients_to_print = coefs_to_print, data_long),
            run_lme(exposure="genotype", outcome=m, outcome_label = str_to_sentence(m),
                    covariates=m_lifestyl_ppInt, coefficients_to_print = coefs_to_print, data_long) %>%
              filter(grepl(" x ", Exposure)))
})) #%>% fwrite("../output/tab_result_lme_glycemic_mmtt_primary_v2.csv")



## EOF 



