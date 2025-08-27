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
analysis <- readRDS("../data/processed/premier_analysis.rda")


####################################################
##  ----------  Build Variable Groups ----------  ##
####################################################

## Table 1 vars --------------------
 
# Demographics
demographic_vars <- c(
  age="Age, years", Sex="Sex", Race2="Race", Ethnicity="Ethnicity", bmi="BMI, kg/m2", 
  whr="Waist-Hip", sbp="SBP, mmHg", dbp="DBP, mmHg", 
  hba1c="HbA1c", homair="HOMA-IR", homab="HOMA-B")

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


## Main outcomes --------------------

# Postprandial glycemic traits
postprandial_vars <- c(paste0(
  rep(c("Glucose ", "Insulin "), each=5), rep(c("0", "30", "60", "120", "180"), 2), rep("min", 10), rep(paste0(" ",Unit), 10)))
names(postprandial_vars) <- c(paste0(
  rep(c("glucose", "insulin"), each=5), rep(c("", "_30", "_60", "_120", "_180"), 2))) 

# iAUC glycemic traits
iAUC_vars <- c(paste0(
  rep(c("Glucose ", "Insulin "), each=7), rep(c("30", "60", "120", "180", "235", "270", "300", "360"), 2), rep("min iAUC (pos)",14), rep(paste0(" ", Unit), 10)))
names(iAUC_vars) <- c(paste0(
  rep(c("glucose_", "insulin_"), each=7), rep(c("30", "60", "120", "180", "235", "270", "300", "360"), 2),  rep("iAUC.pos",14)))



#####################################################
##  ----------  Descriptive statistics ----------  ##
#####################################################

## Set lists of variables & strata -----------------

summary_table_vars.l <- list(
  list(label="demographics", vars=demographic_vars), 
  list(label="metabolites", vars=metabolite_vars),
  list(label="survey", vars=survey_vars),
  list(label=paste0("postprandial",unit), vars=postprandial_vars),
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
      p_types = "descriptive", p_smalln=T) #%>%
      #fwrite(paste0("../output/tab_descr_", vars$label, "_", strata$label, ".csv"), row.names = T)
  })
})
  


#############################################################
##  ---------- Postprandial glycemic responses  ---------- ##
#############################################################

#Aim 1: Determine the impact of genetic category on postprandial glycemic responses
#-linear mixed effects models of OGTT glucose
#-iAUC glucose & iAUC insulin
#-outcomes: OGTT glucose and insulin
#-covariates: age, sex, 5 genetic PCs, BMI, baseline glucose, time_0 levels

## MMTT
#main_effects_model=lmer(METAB ~ B0 + B1*genotype + B2*time + B3*{age, sex, bmi, PC1, PC2, PC3} + e)
#interaction_effects_model=lmer(METAB ~ B0 + B1*genotype + B2*time + B3*genotype*time + B4*{age, sex, bmi, PC1, PC2, PC3} + e)

## Self-Selected
#main_effects_model=lmer(METAB ~ B0 + B1*genotype + B2*time + B3*meal_choice + B4*{age, sex, bmi, PC1, PC2, PC3} + e)
#interaction_effects_model=lmer(METAB ~ B0 + B1*genotype + B2*time + B3*genotype*time + B4*meal_choice + B5*genotype*meal_choice + B6*{age, sex, bmi, PC1, PC2, PC3} + e)

## Create a long dataframe, by time
pp_vars <- c("0", "30", "60", "120", "180", "235", "270", "300", "360")
iAUC.pos_vars <- c("0iAUC.pos", "30iAUC.pos", "60iAUC.pos", "120iAUC.pos", "180iAUC.pos", "235iAUC.pos", "270iAUC.pos", "300iAUC.pos", "360iAUC.pos")
iAUC.tot_vars <- c("0iAUC.tot", "30iAUC.tot", "60iAUC.tot", "120iAUC.tot", "180iAUC.tot", "235iAUC.tot", "270iAUC.tot", "300iAUC.tot", "360iAUC.tot")

postprandial <- analysis %>% dplyr::select(id, starts_with("glucose"), starts_with("insulin"), 
                                    genotype, meal_choice, genetic_cat_x_meal_choice,
                                    age, sex, bmi, PC1z, PC2z, PC3z,
                                    smoke_packyears, alcohol_perweek_4lvl, physact_mvpa_perweek, sleep_hrday) %>%
  rename(glucose.mmtt_0=glucose, insulin.mmtt_0=insulin) %>%
  mutate(glucose.selected_0=glucose_235, insulin.selected_0=insulin_235) %>%
  pivot_longer(c(starts_with("glucose_"), starts_with("insulin_")),
               names_sep="_", names_to=c("metabolite", "time")) %>%
  pivot_wider(names_from=metabolite) %>%
  mutate(measure=ifelse(time %in% pp_vars, "pp", ifelse(time %in% iAUC.pos_vars, "iAUC.pos", ifelse(time %in% iAUC.tot_vars, "iAUC", NA)))) %>%
  mutate(time=factor(time, levels=c(0,30,60,120,180,235,270,300,360))) %>%
  arrange(time) ; head(postprandial)

#saveRDS(postprandial, paste0("../data/processed/postprandial_long",unit,".rda"))


################################################## 
## Mixed-Meal Tolerance Tests (time 0-235)
##################################################

# List metabolites & meal types
metabolites <- c("glucose", "insulin")
meal_types <- c("mmtt", "selected")
measures <- c("pp", "iAUC.pos", "iAUC")
coefs_to_print = c("Genotype"="genotype","Meal Choice"="meal_choice")

models_mmtt.l <- list(
  empty = list(main="genotype+time", int="genotype*time"), 
  primary=list(main="age+sex+PC1z+PC2z+PC3z+time", 
               int="age+sex+PC1z+PC2z+PC3z+time+genotype*time"),
  bmi=list(main="age+sex+PC1z+PC2z+PC3z+time+bmi", 
           int="age+sex+PC1z+PC2z+PC3z+time+genotype*time+bmi")
)

#                adjTime0=list(main=paste0("age+sex+bmi+PC1z+PC2z+PC3z+time+", m,".mmtt_0"), 
#                               int=paste0("age+sex+bmi+PC1z+PC2z+PC3z+time+genotype*time+",m,".mmtt_0")))

models_iAUC.pos.l <- list(
  empty="genotype",
  primary="age+sex+PC1z+PC2z+PC3z",
  bmi="age+sex+PC1z+PC2z+PC3z+bmi")
  #adjTime0=paste0("age+sex+bmi+PC1z+PC2z+PC3z+",m,"_0"))


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
  fwrite(paste0("../output/tab_result_mmtt_postprandial_mgdL.csv"))


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
      fwrite(paste0("../output/tab_result_mmtt_postprand_iAUCpos.csv"))


###############################################
## Self-Selected Meals (time 235-360)
###############################################

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
 

## Postprandial responses by meal type ----------------------------------

# meal choice x time
do.call(rbind, lapply(metabolites, function(m) {
  # set model covariates & subset data_long
  models <- list(main=paste0("time"), 
                 int=paste0("meal_choice*time"))
  data_long <- postprandial %>% filter(time %in% c(235,270,300,360)) 
  # run lme models
  rbind(run_lme(exposure="meal_choice", outcome=m, outcome_label = str_to_sentence(m),
                covariates=models$main, coefficients_to_print = c(coefs_to_print, "Time"="time"), 
                data_long = data_long, digits = 1),
        run_lme(exposure="meal_choice", outcome=m, outcome_label = str_to_sentence(m),
                covariates=models$int, coefficients_to_print = c(coefs_to_print, "Time"="time"), 
                data_long = data_long,  digits = c(1,3)) %>%
          filter(grepl(" x ", Exposure))) 
})) %>% mutate(Model="empty", .before="Exposure") %>% 
  mutate_all(., ~ifelse(is.na(.)==T, "-", .))  %>% 
  mutate(Meal.Type="Self-Selected", .before="Model") %>%
  select(Meal.Type, Model, Effect, Outcome, Exposure, beta, se, p, anovaF, anovaP, Beta.SE, P, P_signif, anovaP_signif) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  fwrite(paste0("../output/tab_result_postprandial_selected_mealonly_mgdL.csv"))


## iAUC by meal type (meal type x time)
do.call(rbind, lapply(metabolites, function(m) {
  do.call(rbind, lapply(c(270,300,360), function(t) {
    y <- paste0(m,"_",t,"iAUC.pos") ; 
    rbind.data.frame(
        print_lm(exposure = "meal_choice", outcome=y, lm_trend = F,
                 covariates = paste0("age+sex+bmi+",m,"_235"), label="AdjTime0",
                 data = analysis)) }))
  })) %>%
  mutate(Meal.Type="Self-Selected", Effect="Main") %>%
  mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.))) %>%
  mutate(lowCI=beta-1.96*se, upCI=beta+1.96*se) %>%
  mutate(Beta.SE=ifelse(is.na(beta)==T, "-", sprintf("%s (%s, %s)", round(beta,2), round(lowCI,2), round(upCI,2))), P=format_p(p)) %>%
  mutate(P_signif=format_p_star(p), anovaP_signif=format_p_star(f_p)) %>%
  select(Meal.Type, Model=model, Effect, Outcome=outcome, Exposure=exposure, N=n, beta, se, p, anovaF=f, anovaP=f_p, Beta.SE, P, P_signif, anovaP_signif) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  fwrite("../output/tab_result_lme_glycemic_selected_iauc_mealchoice_nogentics_v3.csv")



### Postprandial responses by genotype x meal type ----------------------------------
print_summary_table(data=analysis, vars_to_summarize = c(
  glucose_270iAUC.pos="Glucose, 270 iAUC", glucose_300iAUC.pos="Glucose, 300 iAUC", glucose_360iAUC.pos="Glucose, 360 iAUC",
  insulin_270iAUC.pos="Insulin, 270 iAUC", insulin_300iAUC.pos="Insulin, 300 iAUC", insulin_360iAUC.pos="Insulin, 360 iAUC"), 
                    var_strata="genetic_cat_x_meal_choice")

# Glucose at baseline for self-selected MMTT
t.test(glucose_270iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal")) #P=0.947
t.test(glucose_270iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal")) #P=0.0421
t.test(glucose_300iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal")) #P=0.721
t.test(glucose_300iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal")) #P=0.014
t.test(glucose_360iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal")) #P=0.537
t.test(glucose_360iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal")) #P=0.03596

# Insulin
t.test(insulin_270iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal")) #P=0.1791
t.test(insulin_270iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal")) #P=0.3778
t.test(insulin_300iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal")) #P=0.1445
t.test(insulin_300iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal")) #P=0.3387
t.test(insulin_360iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HC meal")) #P=0.08116
t.test(insulin_360iAUC.pos~genotype, data=analysis %>% filter(meal_choice=="HF meal")) #P=0.455


## Set refernece group**
postprandial$genotype <- relevel(postprandial$genotype, ref="HC genotype")
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
  fwrite(paste0("../output/tab_result_postprandial_selected_allmodels_bymeal_HFref_mgdL.csv"))


# stratified by Genetic Category
do.call(rbind, lapply(1:length(models_selected.l$strat_geno), function(m) {
  do.call(rbind, lapply(metabolites, function(metab) {
    data_long <- postprandial %>% filter(time %in% c(235,270,300,360)) 
    do.call(rbind, lapply(c("HC genotype", "HF genotype"), function(geno) {
      rbind(run_lme(exposure="meal_choice", outcome=metab, outcome_label = str_to_sentence(metab),
                    covariates=models_selected.l$strat_geno[[m]]$main, coefficients_to_print = coefs_to_print, 
                    data_long = data_long %>% filter(genotype == geno)),
            run_lme(exposure="meal_choice", outcome=metab, outcome_label = str_to_sentence(metab),
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
  fwrite(paste0("../output/tab_result_postprandial_selected_allmodels_bygeno_HFref_mgdL.csv"))



# iAUCs ---------------------------------------------------

# stratified by Meal Choice
tab_iAUC_selected_meal <- do.call(rbind, lapply(metabolites, function(m) {
  do.call(rbind, lapply(c("HC meal", "HF meal"), function(meal) {
    do.call(rbind, lapply(c(270,300,360), function(t) {
      y <- paste0(m,"_",t,"iAUC.net") ; 
      rbind.data.frame(
        # primary model
        print_lm(exposure = "genotype", outcome=y, lm_trend = F,
                 covariates = "age+sex+bmi+PC1z+PC2z+PC3z", label="Primary",
                 data = analysis %>% filter(meal_choice==meal)),
        # adjusting for metabolite at time0
        print_lm(exposure = "genotype", outcome=y, lm_trend = F,
                 covariates = paste0("age+sex+bmi+PC1z+PC2z+PC3z+",m,"_235"), label="AdjTime0",
        data = analysis %>% filter(meal_choice==meal))) }))  %>% 
      mutate(strata="Meal Choice", strata_level=meal, .before=outcome) 
  }))
  })) %>%
  mutate(Meal.Type="Self-Selected", Effect="Main") %>%
  mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.))) %>%
  mutate(lowCI=beta-1.96*se, upCI=beta+1.96*se) %>%
  mutate(Beta.SE=ifelse(is.na(beta)==T, "-", sprintf("%s (%s, %s)", round(beta,2), round(lowCI,2), round(upCI,2))), P=format_p(p)) %>%
  mutate(P_signif=format_p_star(p), anovaP_signif=format_p_star(f_p)) %>%
  select(Strata=strata, Level=strata_level, Meal.Type, Model=model, Effect, Outcome=outcome, Exposure=exposure, N=n, beta, se, p, anovaF=f, anovaP=f_p, Beta.SE, P, P_signif, anovaP_signif) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  arrange(desc(Model), Level)


# stratified by Genetic Category
tab_iAUC_selected_geno <- do.call(rbind, lapply(metabolites, function(m) {
  do.call(rbind, lapply(c("HC genotype", "HF genotype"), function(geno) {
    do.call(rbind, lapply(c(270,300,360), function(t) {
      y <- paste0(m,"_",t,"iAUC.pos") ; 
      rbind.data.frame(
        # primary model
        print_lm(exposure = "meal_choice", outcome=y, lm_trend = F,
                 covariates = "age+sex+bmi+PC1z+PC2z+PC3z", label="Primary",
                 data=analysis %>% filter(genotype==geno)),
        # adjusting for metabolite at time0
        print_lm(exposure = "meal_choice", outcome=y, lm_trend = F,
                 covariates = paste0("age+sex+bmi+PC1z+PC2z+PC3z+",m,"_235"), label="AdjTime0",
                 data=analysis %>% filter(genotype==geno) )) })) %>% 
    mutate(strata="Genetic Category", strata_level=geno, .before=outcome)
  })) })) %>%
  mutate(Meal.Type="Self-Selected", Effect="Main") %>%
  mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.))) %>%
  mutate(lowCI=beta-1.96*se, upCI=beta+1.96*se) %>%
  mutate(Beta.SE=ifelse(is.na(beta)==T, "-", sprintf("%s (%s, %s)", round(beta,2), round(lowCI,2), round(upCI,2))), P=format_p(p)) %>%
  mutate(P_signif=format_p_star(p), anovaP_signif=format_p_star(f_p)) %>%
  select(Strata=strata, Level=strata_level, Meal.Type, Model=model, Effect, Outcome=outcome, Exposure=exposure, N=n, beta, se, p, anovaF=f, anovaP=f_p, Beta.SE, P, P_signif, anovaP_signif) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  arrange(desc(Model), Level)


rbind.data.frame(tab_iAUC_selected_meal, tab_iAUC_selected_geno) %>% 
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  fwrite(paste0("../output/tab_result_lm_glycemic_selected_iAUC_allModels_v3",unit,".csv"))




#######################
## Beta cell Indices ##
#######################

# primary model
rbind.data.frame(
  do.call(rbind, lapply(c("mmtt", "selected"), function(x) {
    model <- "age+sex+bmi+PC1z+PC2z+PC3z"
    y <- paste0("matsuda_",x)
    print_lm(exposure = "genotype", outcome=y, lm_trend = F,
             covariates = "age+sex+bmi+PC1z+PC2z+PC3z", label="Total",
             data=analysis) })), 
  
  do.call(rbind, lapply(c("HC genotype", "HF genotype"), function(x) {
    model <- "age+sex+bmi+PC1z+PC2z+PC3z"
    print_lm(exposure = "meal_choice", outcome="matsuda_selected", lm_trend = F,
             covariates = "age+sex+bmi+PC1z+PC2z+PC3z", label=x,
             data=analysis %>% filter(genotype==x)) })), 
  
  do.call(rbind, lapply(c("HC meal", "HF meal"), function(x) {
    model <- "age+sex+bmi+PC1z+PC2z+PC3z"
    print_lm(exposure = "genotype", outcome="matsuda_selected", lm_trend = F,
             covariates = "age+sex+bmi+PC1z+PC2z+PC3z", label=x,
             data=analysis %>% filter(meal_choice==x)) }))
  ) %>% 
  mutate(across(c(beta,se,p,f,f_p), ~as.numeric(.))) %>%
  mutate(lowCI=beta-1.96*se, upCI=beta+1.96*se) %>%
  mutate(Beta.SE=ifelse(is.na(beta)==T, "-", sprintf("%s (%s, %s)", round(beta,2), round(lowCI,2), round(upCI,2))), P=format_p(p)) %>%
  mutate(P_signif=format_p_star(p)) %>%
  mutate(anovaP_signif=format_p_star(f_p)) %>%
  mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
  fwrite("../output/tab_result_lm_indices_selected_v3.csv")


## 



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
  summary(lm(formula(paste0(y, "~genotype++age+sex+bmi+PC1z+PC2z+PC3z")), data=analysis))
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
                    var_strata = "genotype", var_strata_order = c("Prefer Carb", "Prefer Fat")) %>%
  kable(caption="MMTT nutrient intakes by genetic category")


lapply(selected_meal_data, function(y){
  anova(lm(formula(paste0(y, "~genotype+meal_choice+age+sex+bmi+PC1z+PC2z+PC3z")), data=analysis))
}) # -significant ANOVA for genotype effects on selected_carb_g and selected_fat_g, after adjusting for meal_choice


lapply(selected_meal_data, function(y){
  anova(lm(formula(paste0(y, "~genotype*meal_choice+age+sex+bmi+PC1z+PC2z+PC3z")), data=analysis))
}) # -significant ANOVA for genotype*meal_choice interaction on selected_carb_g and selected_fat_g


# =============================================
## Self-Selected - Postprandial metbolites 
# =============================================

## 2x2 interactioin
#do.call(rbind, lapply(metabolites, function(m) {
#  times=c(235,270,300,360) ; main_effects_formatted="age+sex+bmi+PC1z+PC2z+PC3z+time+meal_choice+selected_carb_g+selected_fat_g"
#  interaction_effects_formatted="genotype*time+genotype+meal_choice"
#  do.call(rbind.data.frame, list(
#    pp=run_lme_postprandial(metabolite = m, times=times, measure="pp", meal_type = "selected",
#                            main_effects_formatted=main_effects_formatted, interaction_effects_formatted=interaction_effects_formatted),
#    iAUC=run_lme_postprandial(
#      metabolite = m, times=times[-1], measure="iAUC", meal_type = "selected",
#      main_effects_formatted=main_effects_formatted, interaction_effects_formatted=interaction_effects_formatted))
#  ) })) %>%
#  fwrite("../output/tab_result_lme_glycemic_selected_sensitivity_adjCarbFat.csv")


## 4-way interaction: Genetic Category x Meal Choice
#do.call(rbind, lapply(metabolites, function(m) {
#  times=c(235,270,300,360) ; main_formatted="age+sex+bmi+PC1z+PC2z+PC3z+time+selected_carb_g+selected_fat_g"
#  interaction_formatted="genetic_cat_x_meal_choice*time"
#  do.call(rbind.data.frame, list(
#    pp=run_lme_postprandial(exposure = "genetic_cat_x_meal_choice", metabolite = m, times=times, measure="pp", meal_type = "selected",
#                            main_effects_formatted=main_formatted, interaction_effects_formatted=interaction_formatted),
#    iAUC=run_lme_postprandial(exposure = "genetic_cat_x_meal_choice", metabolite = m, times=times[-1], measure="iAUC", meal_type = "selected",
#                              main_effects_formatted=main_formatted, interaction_effects_formatted=interaction_formatted))
#  ) })) %>%
#  fwrite("../output/tab_result_lme_glycemic_selected_4levelGxM_sensitivity_adjCarbFat.csv")


# =============================================
## Self-Selected - iAUC metbolites 
# =============================================

## 2x2 interaction
#do.call(rbind, lapply(iAUC_vars_selected, function(metabolite) {
#  bind_rows(print_lm(exposure = "genotype", outcome=metabolite, lm_trend = F,
#                     covariates = "age+sex+bmi+PC1z+PC2z+PC3z+meal_choice+selected_carb_g+selected_fat_g", label="Main") %>%
#              mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.))),
#            print_lm_interaction(exposure = "genotype", interaction = "meal_choice",
#                                 outcome=metabolite, 
#                                 model = "age+sex+bmi+PC1z+PC2z+PC3z+meal_choice+selected_carb_g+selected_fat_g",
#                                 label_model="Interaction", lm_trend = F, data=analysis,
#                                 label_interaction = "Genetic Category x Meal Choice") %>%
#              mutate(outcome=metabolite, model="Interaction", exposure=rownames(.), .before=n) %>%
#              filter(startsWith(exposure, "Genetic Category")) )
#})) %>% as.data.frame() %>% arrange(outcome, model) %>% 
#  mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.))) %>%
#  mutate(lowCI=beta-1.96*se, upCI=beta+1.96*se) %>%
#  mutate(Beta.SE=ifelse(is.na(beta)==T, "-", sprintf("%s (%s, %s)", round(beta,2), round(lowCI,2), round(upCI,2))), P=format_p(p)) %>%
#  mutate(P_signif=format_p_star(p)) %>%
#  fwrite("../output/tab_result_lm_glycemic_selected_sensitivity_adjCarbFat.csv", row.names = F)


# 4-level interaction: Self-Selected 
#do.call(rbind, lapply(iAUC_vars_selected, function(metabolite) {
#  print_lm(exposure = "genetic_cat_x_meal_choice", outcome=metabolite, lm_trend = F,
##           covariates = "age+sex+bmi+PC1z+PC2z+PC3z+selected_carb_g+selected_fat_g", label="Main") %>%
#    mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.)))
#})) %>% as.data.frame() %>% arrange(outcome, model) %>% 
#  mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.))) %>%
#  mutate(lowCI=beta-1.96*se, upCI=beta+1.96*se) %>%
#  mutate(Beta.SE=ifelse(is.na(beta)==T, "-", sprintf("%s (%s, %s)", round(beta,2), round(lowCI,2), round(upCI,2))), P=format_p(p)) %>%
#  mutate(P_signif=format_p_star(p)) %>%
#  fwrite("../output/tab_result_lm_glycemic_selected_4levelGxM_sensitivity_adjCarbFat.csv")


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



# ====================
## Sex-Stratified
# ====================

# postprandial (primary model)

do.call(rbind.data.frame, lapply(c("0", "1"), function(x){
    do.call(rbind, lapply(metabolites, function(m) {
      to_print<-coefs_to_print[c(1,3)]
      models <- list(main="age+bmi+PC1z+PC2z+PC3z+time+genotype*time",
                    int=paste0("age+bmi+PC1z+PC2z+PC3z+time+genotype*time+",m,".selected_0"))
      data_long <- postprandial %>% filter(time %in% c(0,30,60,120,180)) 
      # run lme models
      rbind(run_lme(exposure="genotype", outcome=m, outcome_label = str_to_sentence(m),
                    covariates=models$main, coefficients_to_print = to_print, 
                    data_long = data_long %>% filter(sex==x)),
            run_lme(exposure="genotype", outcome=m, outcome_label = str_to_sentence(m),
                    covariates=models$int, coefficients_to_print = to_print, 
                    data_long = data_long %>% filter(sex==x)) %>%
              filter(grepl(" x ", Exposure)))
    })) %>% mutate(Sex=ifelse(x==0,"Female","Male")) 
  })) %>% 
  mutate_all(., ~ifelse(is.na(.)==T | .=="NA (NA, NA)", "-", .)) %>%
  fwrite(paste0("../output/tab_result_lm_glycemic_selected_pp_sensitivity_sex_v3",unit,".csv"))

# iAUC (primary model)

# primary model & adjusting for metabolite at time0
  do.call(rbind, lapply(c("Female", "Male"), function(x) {
    do.call(rbind, lapply(metabolites, function(m) {
      do.call(rbind, lapply(c(30,60,120,180), function(t) {
        y <- paste0(m,"_",t,"iAUC.pos") 
        print_lm(exposure = "genotype", outcome=y, lm_trend = F,
                 covariates = "age+sex+bmi+PC1z+PC2z+PC3z", label="Primary",
                 data=analysis %>% filter(Sex==x)) })) })) %>%
      mutate(Sex=x) })) %>%
      mutate(across(c(n, beta, se, p, f, f_p), ~as.numeric(.))) %>%
      mutate(lowCI=beta-1.96*se, upCI=beta+1.96*se) %>%
      mutate(Beta.SE=ifelse(is.na(beta)==T, "-", sprintf("%s (%s, %s)", round(beta,2), round(lowCI,2), round(upCI,2))), P=format_p(p)) %>%
      mutate(P_signif=format_p_star(p), anovaP_signif=format_p_star(f_p)) %>%
      select(Model=model, Sex, Outcome=outcome, Exposure=exposure, N=n, beta, se, p, anovaF=f, anovaP=f_p, Beta.SE, P, P_signif, anovaP_signif) %>%
      mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
      arrange(desc(Model)) %>%
      mutate_all(., ~ifelse(is.na(.)==T, "-", .)) %>%
    mutate_all(., ~ifelse(is.na(.)==T | .=="NA (NA, NA)", "-", .)) %>%
    fwrite(paste0("../output/tab_result_lm_glycemic_selected_iAUC_sensitivity_sex_v3",unit,".csv"))




## EOF 



