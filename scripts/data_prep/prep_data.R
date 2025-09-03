# Prepare PREMIER data 

#Load libraries
library(Hmisc)
library(tidyverse)
library(data.table)
library(devtools)
library(ggpubr)
library(readxl)


# Load pantry "package"
lapply(list.files("../../pantry/functions/", recursive = T, full.names = T), source)


## CHOOSE UNITS FOR DATA PREPARATION (mg/dL or mmol/L)

units <- c("", "_mmol")
Units <- c("mg/dL", "mmol/L")

unit<-units[1] # ""
Unit<-Units[1] # mg/dL


##########################
## Build base dataframe ##
##########################

# Define descriptive variable levels
sex_labs <- c("Female"=0, "Male"=1)
yes_no_labs <- c("Yes"=1, "No"=0)
race_labs <- c("White"=1, "Black"=2, "Asian"=3, "More than 1 race"=4, "Unknown or not reported"=5, 
               "American Indian/Alskan native"=6, "Native Hawaiin or Other Pacific Islander"=7, "Other"=8)
asian_labs <- c("East Asian"=1, "South Asian"=2)
ethn_labs <- c("Hispanic"=1, "Not Hispanic"=2, "Unknown"=3)
language_labs <- c("English"=1, "Spanish"=2, "Other"=3)
recruitment_labs <- c("MGB Biobank"=1, "RODY"=2)
diab_tmt_labs <- c("Diet-controlled"=1, "Not on medication"=2, "Used to be on medication"=3)
complete_labs <- c("Incomplete"=0, "Unverified"=1, "Complete"=2)


# ===========================================
## Read In Participant & Recruitment Data 
# ===========================================

#Read In results data to get participant IDs
allbase.dat=fread('../data/raw/PREMIERResults_DATA_2024-12-18_1509.csv') %>%
  select(!starts_with("qi") & !starts_with("tf")) %>%
  filter(study_id != "test subject") %>%
  filter(study_id != "P26")

included_participants <- allbase.dat$study_id
length(included_participants) #N=24


#PremierNames data
names.dat <- fread("../data/raw/PREMIERNames_DATA_2025-01-07_1031.csv")
names_mrn.dat <- names.dat %>% 
  mutate(Recruitment=add_descr_labels(., "recruitment", recruitment_labs),
         Language=add_descr_labels(., "primary_language", language_labs)) %>%
  select(id=study_id, first_name, last_name, Recruitment, MRN=mrn, dob) %>%
  mutate_at("MRN", ~as.character(.)) %>%
  mutate(initials=paste0(substr(first_name, 1,1), substr(last_name, 1,1))) %>%
  filter(id %in% included_participants)
  

#Recruitment data from MGBB 
mgbb <- readxl::read_xlsx("../data/raw/PREMIER_StudyIDs_Extremes_BIOBNAK_NotRODY_EMPI_20250123.xlsx") %>%
  rename("MGH"="MGH...1") %>% select(-"MGH...5") %>%
  mutate(across(c("MGH", "EMPI", "BWH", "Subject.ID"), ~as.character(.)))


#Recruitment data from RODY
rody <- readxl::read_xlsx("../data/raw/Eligible.RODY.extremes_PREMIER.xlsx") %>%
  mutate_at(c("MGH", "EMPI", "Subject.ID", "MRN"), ~as.character(.)) %>%
  separate(MRN, into = c("MRN_1", "MRN_2", "MRN_3", "MRN_4", "MRN_5"), sep = ", ") %>%
  select(Subject.ID, study.group, starts_with("MRN")) %>%
  pivot_longer(paste0("MRN_", 1:5), values_to="MRN") %>%
  filter(!is.na(MRN))


#Re-calculated genetic category data
genetic_category <- fread("../data/raw/genetic_category_20250212.csv") %>%
  mutate_at("Subject.ID", ~as.character(.))


#Additional premier participant IDs pulled from mgbb 
addn_mgbb<-fread("../data/raw/addn_premier_subjectIDs_20250209.csv") %>% 
  select(Subject.ID, MRN) %>% 
  mutate(across(c("MRN", "Subject.ID"), ~as.character(.)))


#Meal selection data
meal_choice<-readxl::read_xlsx("../data/raw/PREMIER_Meal_Data.xlsx") %>% 
  select(id=`Participant  ID`, meal_choice=Choice, initials=Initials) %>%
  filter(!is.na(id))
meal_choice %>% print(n=30) # NOTE: two participant IDs missing; match using initials
meal_choice %>% filter(id=="Not provided") #Initials == BB & BG

#check which participants have missing meal_choice data
allbase.dat$study_id[which(!allbase.dat$study_id %in% meal_choice$id)] # P11 & P23
names_mrn.dat %>% filter(id %in% c("P11", "P23")) #P11=BB; P23=BG
meal_choice$id[meal_choice$id == "Not provided" & meal_choice$initials=="BB"] <- "P11"
meal_choice$id[meal_choice$id == "Not provided" & meal_choice$initials=="BG"] <- "P23"
meal_choice %>% print(n=30) # All values have complete participant IDs


#Meal consumption data
meal_data_names <- c(id="STUDY ID", energy_kcal="Energy", prot_g="Protein (g)", 
                     fat_g="Total lipid (fat) (g)", carb_g="Carbohydrate (g)", 
                     prot_pct="Protein (%)", fat_pct="Fat (%)", carb_pct="Carbohydrate (%)")
meal_data_mmtt=readxl::read_xlsx("../data/raw/PREMIER_Meal_Data_Breakfast_MMTT.xlsx") %>%
  select(all_of(meal_data_names)) %>% rename_with(~paste0("mmtt_", .), -id)
meal_data_selected=readxl::read_xlsx("../data/raw/PREMIER_Meal_Data_Lunch_Baked_Ziti.xlsx") %>%
  select(all_of(meal_data_names)) %>% rename_with(~paste0("selected_", .), -id)
meal_data <- inner_join(meal_data_mmtt, meal_data_selected, by = "id")



# ==================================================================
## Prepare additional MGBB phenotypes (genetic PCs & covariates)
# ==================================================================

#Genetic ancestry PCs
genetic_pcs<-fread("../data/raw/MGBB.64K.Genetic_Ancestry_premier.txt") %>%
  select(Subject.ID=`Biobank Subject ID`, paste0("PC", 1:5)) %>%
  mutate(Subject.ID=as.character(Subject.ID))


# MGBB survey data (demographic & lifestyle covariates)
survey_vars <- c(
  education="What is the highest grade in school that you finished?",
  smoke_gt100 = "Have you smoked at least 100 cigarettes in your lifetime?",
  smoke_agerange="During which age ranges did you smoke? (please check all ages that apply)",
  smoke_packyears="Smoke Exposure in Pack-years",
  alcohol_perweek="During the past year- how many alcoholic drinks (glass/bottle/can of beer; 4oz glass of wine; drink or shot of liquor) did you usually drink in a typical week?",
  physact_mvpa_perweek="Total moderate to high-intensity exercise (excludes walking/hiking) in hours per week",
  sleep_hrday="Total average sleep duration in hours per day"
)

survey_data<-fread("../data/raw/mgbb_phenos_premier.txt") %>%
  select(Subject.ID=id, all_of(survey_vars)) %>%
  mutate(Subject.ID=as.character(Subject.ID)) %>%
  mutate_all(., ~ifelse(.=="", NA, .)) %>%
  mutate_at("education", ~case_when(.=="Four year college"~"College degree", 
                                .== "Masters, doctoral, or professional degree"~"Graduate or Professional degree",
                                .=="Some college"~ "High School degree",
                                is.na(.) ~ as.character(NA))) %>%
  mutate_at("smoke_packyears", ~case_when(.=="Current Smokers"~"Current", .=="Past Smokers"~"Former",
                                          is.na(.) & smoke_gt100==0~"Never", 
                                          is.na(.) & smoke_gt100==1~"Former")) %>%
  mutate_at("alcohol_perweek", ~case_when(.=="None, or less than 1 per month"~"Less than 1 per month",
                                           .default=as.character(.))) %>%
  mutate("alcohol_perweek_4lvl" = case_when(alcohol_perweek %in% c("1-2 per day", "5-6 per week") ~ "5+ per week",
                                             alcohol_perweek %in% c("2-4 per week", "1 per week") ~ "1-4 per week",
                                             .default=as.character(alcohol_perweek))) %>%
  mutate(across(c("physact_mvpa_perweek", "sleep_hrday"), ~as.numeric(.))) %>%
  
  # Add meaningful factor levels
  mutate_at("education", ~factor(., levels=c("Graduate or Professional degree", "College degree", "High School degree"))) %>%
  mutate_at("smoke_packyears", ~factor(., levels=c("Never", "Former", "Current"))) %>%
  mutate_at("alcohol_perweek", ~factor(., levels=c("1-2 per day", "5-6 per week", "2-4 per week", "1 per week", "1-3 per month", "Less than 1 per month"))) %>%
  mutate_at("alcohol_perweek_4lvl", ~factor(., levels=c("5+ per week", "1-4 per week", "1-3 per month", "Less than 1 per month")))


# ========================================================
## Merge base + names/recruitment/genetic category data
# ========================================================

names_mrn.dat <- names_mrn.dat %>%
  left_join(rody %>% select(Subject.ID, MRN, study.group), by = "MRN") %>% 
  left_join(mgbb %>% select(MRN=MGH, study.group, Subject.ID), by = c("MRN"), relationship = "many-to-many") %>% 
  left_join(addn_mgbb, by = c("MRN")) %>%
  rowwise() %>%
  mutate(pps_cat = ifelse(!is.na(study.group.x), study.group.x, study.group.y)) %>%
  mutate(Subject.ID.z = ifelse(is.na(Subject.ID.x), Subject.ID.y, Subject.ID.x)) %>%
  mutate(Subject.ID = ifelse(is.na(Subject.ID.z), Subject.ID, Subject.ID.z)) %>%
  select(-ends_with(".x"), -ends_with(".y"), -ends_with(".z"))
  

## Load in base data 
base.dat <- allbase.dat %>% 
  mutate(Sex=add_descr_labels(., "sex", sex_labs),
         Race=add_descr_labels(., "race", race_labs),
         Asian=add_descr_labels(., "asian", asian_labs),
         Ethnicity=add_descr_labels(., "ethnicity", ethn_labs),
         Demographics_Complete=add_descr_labels(., "demographics_complete", complete_labs),
         Vitals_Complete=add_descr_labels(., "vitals_complete", complete_labs),
         Visit_Results_Complete=add_descr_labels(., "visit_results_complete", complete_labs)
  ) %>%
  select(id=study_id, age, sex, Sex, Race, Ethnicity, bmi, 
         whr=calculated_whr, sbp=systolic_bp, dbp=diastolic_bp, pulse=resting_pulse,
         starts_with("glucose"), starts_with("insulin"),
         hba1c, protein=total_protein, alb=albumin, globulin, alb_globulin_ratio=albumin_globulin_ratio,
         tot_bilirubin=total_bilirubin, dir_bilirubin=direct_bilirubin, alk_phosphatase=alkaline_phosphatase,
         ast, alt, creatinine, gfr, chol=total_cholesterol, tg=triglycerides, 
         hdl=hdl_cholesterol, ldl=ldl_cholesterol, vldl=vldl_cholesterol,
         wbc, rbc, hemoglobin, hematocrit, mcv, mch, mchc, rdw, platelet_count,
         mpv, neutrophils, lymphocytes, monocytes, eosinophils, basophils, 
         chol_120=total_chol_120m, chol_235=cool_235m, chol_360=chol_360, 
         hdl_120, hdl_235=hdl_235m, hdl_360=hdl_360,
         ldl_120, ldl_235=ldl_235m, ldl_360=ldl_360,
         tg_120=triglyc_120, tg_235=triglyc_235m, tg_360=triglycerides_360) %>%
  rename(glucose=glucose_baseline, insulin=insulin_baseline) %>%
  mutate_at("whr", ~as.numeric(.)) %>%
  mutate(glucose_0=glucose, insulin_0=insulin)
         

## Merge in genetic category (pps_cat) 
base.dat <- base.dat %>% left_join(
  names_mrn.dat %>% select(id, Subject.ID, MRN, pps_cat), by = "id")

## Merge in meal_choice
base.dat <- base.dat %>% left_join(meal_choice, by ="id")

## Merge in genetic ancestry PCs
base.dat <- base.dat %>% left_join(genetic_pcs, by="Subject.ID") %>%
  mutate(PC1z=zscore(PC1), PC2z=zscore(PC2), PC3z=zscore(PC3),
         PC4z=zscore(PC4), PC5z=zscore(PC5))

## Merge in meal consumption data
base.dat <- base.dat %>% left_join(meal_data)


## Merge in survey data
base.dat <- base.dat %>% left_join(survey_data, by = "Subject.ID")
saveRDS(base.dat, "../data/processed/premier_base_data.rda")


##########################
## Clean base dataframe ##
##########################

processed.dat <- base.dat

## Merge in re-calculated genetic category variable
processed.dat <- processed.dat %>% 
  left_join(genetic_category %>% mutate(Subject.ID=as.character(Subject.ID)), by = "Subject.ID") %>%
  rename(genetic_category=study.group.recalc) %>%
  mutate(genetic_category=as.factor(genetic_category)) %>%
  mutate(genetic_category.lab=ifelse(genetic_category==1, "Prefer Carb", ifelse(genetic_category==2, "Prefer Fat", NA))) %>%
  mutate(genetic_category.lab=factor(genetic_category.lab, levels=c("Prefer Carb", "Prefer Fat", NA)))
  

## Add interaction term for genetic category X meal choice
processed.dat <- processed.dat %>% mutate(
  genetic_cat_x_meal_choice = ifelse(
    !is.na(genetic_category.lab) & !is.na(meal_choice), 
    paste0(genetic_category.lab, " x ", meal_choice), NA)
)

# =======================================================================
## Check breakdown of categorical variables & recoded, as needed
# =======================================================================

cat_descr_vars <- c("Sex", "Race", "Ethnicity")
do.call(rbind.data.frame, lapply(cat_descr_vars, function(x) {
  processed.dat %>% select(var=x) %>%
    reframe(n_pct(var))
}))

# Recode Race --> only White/Black/Other
processed.dat <- processed.dat %>% 
  mutate(Race2=case_when(Race == "White" | Race == "Black" | 
                           Race == "Other" ~ Race,
                         TRUE ~ as.character(NA))) %>% 
  mutate(RaceEthn = paste0(Ethnicity, " ", Race2))


# =======================================================================
## Check distributions of continuous variables & apply transformations
# =======================================================================

# Descriptive variables
cont_descr_vars <- c("age", "bmi", "whr", "sbp", "dbp", "dbp", "hba1c")
plots_descr.l <- lapply(cont_descr_vars, function(var) {
  plot_continuous(var, data=processed.dat) + 
    geom_histogram(bins=15) })
ggarrange(plotlist = plots_descr.l)


# Metabolite variables
cont_metabolite_vars <- c("glucose", "insulin","tg", "chol", "hdl", "ldl", "protein", 
                          "alb", "globulin", "alb_globulin_ratio","tot_bilirubin", "dir_bilirubin", "ast", "alt")

plots_metab.l <- lapply(cont_metabolite_vars, function(var) plot_continuous(var, data=processed.dat) + 
                          geom_histogram(bins=15))
ggarrange(plotlist = plots_metab.l)


## Add log-transformed fasting metabolites
processed.dat <- processed.dat %>%
  mutate(tg_log = log(tg),
         tg_log_120 = log(tg_120),
         tg_log_235 = log(tg_235),
         tg_log_360 = log(tg_360)
         )

plots_metab_log.l <- lapply(paste0("tg_log",c("","_120","_235","_360")), function(var) plot_continuous(var, data=processed.dat))
ggarrange(plotlist = plots_metab_log.l)


# =======================================================================
## Add delta bm variables for all combinations (iAUC)
## iAUC: 0-30, 0-60, 0-120, 0-180
# Trapezoidal rule: ((height1+height2)/2 x width1[time])
# =======================================================================

## calculate iAUC using calcAUC
devtools::install_github("scrs-msu/auctime")
library(auctime)

make_iAUC <- function(metab, method, data) { 
  times.l <- list(x_30iAUC=c(0,30), 
                  x_60iAUC=c(0,30,60), 
                  x_120iAUC=c(0,30,60,120), 
                  x_180iAUC=c(0,30,60,120,180),
                  x_235iAUC=c(0,30,60,120,180,235),
                  x_270iAUC=c(235,270), 
                  x_300iAUC=c(235,270,300),
                  x_360iAUC=c(235,270,300,360)
                  )
  tag<-ifelse(method=="positive", ".pos", ifelse(method=="net", ".net", ""))
  AUC <- lapply(1:length(times.l), function(t) {
    calcAUC(data %>% select(subjects=id, c(paste0(metab, "_", times.l[[t]]))),
            biomarker=metab, method=method, subjects = T, interval = "minutes", 
            plot=T)$dataframe %>% as.data.frame() %>%
      rename_all(., ~gsub("x[_]", "", gsub("AUC", paste0(names(times.l)[t],tag), .)))
  }) %>% reduce(full_join,by="Subject") %>%
    rename(id=Subject) %>%
    mutate(across(contains("iAUC"), ~as.numeric(.)))
  AUC
}

processed.dat <- processed.dat %>% 
  left_join(make_iAUC("glucose", method="positive", data=.), by="id") %>%
  left_join(make_iAUC("glucose", method="net", data=.), by="id") %>%
  left_join(make_iAUC("insulin", method="positive", data=.), by="id") %>%
  left_join(make_iAUC("insulin", method="net", data=.), by="id")
 
# =======================================================================
## Calculate beta cell indices
# =======================================================================

## Calculate HOMA-IR & HOMA-B
processed.dat <- processed.dat %>% mutate(
  homair = (insulin*glucose)/405,
  homab = ifelse(glucose <63, NA, (360*insulin)/(glucose-63))
)


# =======================================================================
## Re-code genotype & meal type
# =======================================================================

# Rename as analysis
analysis <- processed.dat

analysis <- analysis %>% 
  mutate(genotype = case_when(genetic_category.lab=="Prefer Carb" ~ "HC genotype",
                              genetic_category.lab=="Prefer Fat" ~ "HF genotype")) %>%
  mutate_at("meal_choice", ~case_when(.=="High Carb" ~ "HC meal",
                                      .=="High Fat" ~ "HF meal")) %>%
  mutate_at("genetic_cat_x_meal_choice", ~
              gsub("Prefer Carb", "HC genotype", gsub("Prefer Fat", "HF genotype", gsub(
                "High Carb", "HC meal", gsub("High Fat", "HF meal", .))))) %>%
  mutate(genotype=factor(genotype, levels=c("HF genotype", "HC genotype"))) %>%
  filter(!is.na(genotype))


dim(analysis) #22 186
saveRDS(analysis, "../data/processed/premier_analysis.rda")


# =======================================================================
## Make dataframe in mmol/L & revise coding of genotype & meals
# =======================================================================

analysis_mmol<-analysis %>% mutate(across(starts_with(c("glucose")), ~./18))
#saveRDS(analysis_mmol, "../data/processed/premier_analysis_mmol.rda")


#EOF


