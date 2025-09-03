## script to match participants to genetic_category

# run after data_prep.R for required input files/functions

premier=readRDS("../data/processed/premier_analysis.rda")

# =======================================================
## Participants with missing genetic category
# =======================================================

missing_pps <- premier %>% filter(is.na(pps_cat))

missing_pps_info <- premier %>% 
  filter(id %in% missing_pps$id) %>% 
  select(id, dob, sex, race, ethnicity, visit_date, age, visit_date) %>%
  left_join(names.dat %>% select(id=study_id, mrn, first_name, last_name, phone_number, email, address, contacted_for_future_studi,
                                 recruitment, recruitment_source, studyidchange, v1_location, date), by = "id") %>%
  select(id, first_name, last_name, mrn, dob, sex, race, phone_number, email, address, visit_date, v1_location, recruitment, recruitment_source) 
missing_pps_info %>% write.csv("../data/processed/PREMIER_IDs_missing_GeneticCategory_20250124.csv")



# ===============================
## Lists of unmatched IDs 
# ===============================

## Matched
match_rody <- analysis %>% select(first_name, last_name, id, MRN, pps_cat.lab) %>% 
  inner_join(rody %>% select(MRN, Subject.ID, study.group), by = "MRN") %>%
  mutate(Subject.ID = as.character(Subject.ID))
match_mgbb <- analysis %>% select(first_name, last_name, id, MRN, pps_cat.lab) %>%
  inner_join(mgbb %>% select(MRN=MGH, Subject.ID, study.group), by = "MRN") %>%
  mutate(Subject.ID = as.character(Subject.ID))

## Unmatched (n=6)
unmatched <- analysis %>% select(first_name, last_name, id, MRN, pps_cat.lab) %>% 
  filter(is.na(pps_cat.lab)) %>%
  rowwise() %>% mutate(MRN_MGH=ifelse(MRN %in% mgbb$MGH | MRN %in% rody$MRN, 1, 0)) 


# =================================================
## Matching MRN/Subject.ID by hand (mgbb biobank)
# =================================================

#paste0(unmatched$MRN, collapse="|")
#mgbb="/humgen/florezlab/mgbb/releases/mgbb_phenotypes_v1/"
#'cat $mgbb/mgbb_dem_v1.csv | grep -E 4962468|2191969|5742712|2524029|4130989|5219380 '

#!# MISSING: Gary Duran #!#

add_from_mgbb <- as.data.frame(cbind(
  EMPI=c(111183304,104040480,110541591,107082588,102850311),
  Subject.ID=c(10029769,10087199,10080426,10010490,10096227),
  MRN=c("2524029","2191969","5742712","5219380","4130989"),
  sex=c("Female","Female","Female","Male","Male"),
  dob=IDateTime(c("1960-12-09","1975-04-27","1963-10-20","1990-01-10","1980-03-31"))
)) %>% select(-dob.itime)

add_from_mgbb %>% write.csv("../data/processed/addn_premier_subjectIDs_20250209.csv")

unmatched <- unmatched %>% inner_join(add_from_mgbb, by = "MRN")

subjectIDs <- c(match_rody$Subject.ID, match_mgbb$Subject.ID, unmatched$Subject.ID)
subjectIDs %>% write.table("./premier_subjectIDs.txt", row.names = F, col.names=F, quote=F)



# =================================================
## Matching based on re-calculated PGS
# =================================================

## 
mgbb %>% 
  pivot_longer(c(CHO.SCORE, FAT.SCORE), values_to="SCORE", names_to="MACRO") %>%
  ggplot(aes(x=SCORE, fill=as.factor(study.group), group=as.factor(study.group))) +
  facet_wrap(~study.group+MACRO, scale="free_y")+
  geom_histogram(alpha=0.75)+
  scale_fill_manual(values=c(palettes$NatMainAccent3$Level2))


mgbb %>% 
  pivot_longer(c(CHO.SCORE, FAT.SCORE), values_to="SCORE", names_to="MACRO") %>%
  mutate(PREMIER=ifelse(Subject.ID %in% subjectIDs, T, F)) %>%
  ggplot(aes(x=MACRO, y=SCORE, color=as.factor(study.group),
             shape=PREMIER, size=PREMIER, alpha=PREMIER)) +
  geom_point(position = position_jitter(0.25))+
  scale_color_manual(values=c(palettes$NatMainAccent3$Level2)) + 
  scale_shape_manual(values=c(21,19)) + scale_size_manual(values=c(1,3)) + 
  scale_alpha_manual(values=c(0.55,1))

mgbb %>% filter(Subject.ID %in% subjectIDs) %>%
  pivot_longer(c(CHO.SCORE, FAT.SCORE), values_to="SCORE", names_to="MACRO") %>%
  ggplot(aes(x=MACRO, y=SCORE, color=as.factor(study.group), group=Subject.ID)) + 
  geom_hline(yintercept = 0) +
  geom_point(size=3) + geom_line() + 
  scale_color_manual(values=c(palettes$NatMainAccent3$Level2))



# =================================================
## Matching based on re-calculated PGS
# =================================================
  
## load re-calculated scores (v2)
pgs <- fread("../data/processed/macro_allsnps_pgs.S1.sscore") %>%
  mutate(Subject.ID=gsub(".*[_]", "", gsub(".*[-]", "", IID))) %>% 
  select(Subject.ID, carb_score_center=SCORE1_AVG) %>% full_join(
    fread("../data/processed/macro_allsnps_pgs.S2.sscore") %>%
      mutate(Subject.ID=gsub(".*[_]", "", gsub(".*[-]", "", IID))) %>% 
      select(Subject.ID, fat_score_center=SCORE1_AVG),
    by="Subject.ID", relationship="many-to-many")


pgs <- pgs %>% 
  mutate(premier=ifelse(Subject.ID %in% subjectIDs, T, F)) %>%
  left_join(match_rody %>% select(Subject.ID, study.group) %>%
              bind_rows(match_mgbb %>% select(Subject.ID, study.group)),by="Subject.ID") %>%
  mutate(carb_score_center_quintiles=cut(carb_score_center, breaks=quantile(carb_score_center, probs=seq(0,1,0.2)), include.lowest = T, labels=1:5)) %>%
  mutate(fat_score_center_quintiles=cut(fat_score_center, breaks=quantile(fat_score_center, probs=seq(0,1,0.2)), include.lowest = T, labels=1:5)) %>%
  #Recalculate genetic category based on carb_quintile <=> fat_quintile
  mutate(study.group.recalc = ifelse(as.numeric(carb_score_center_quintiles) > as.numeric(fat_score_center_quintiles), 1, 2))

## Compare pgs determinations
table(pgs$study.group, pgs$study.group.recalc) # re-classified all 18 participants correctly!
pgs %>% filter(premier==T) %>% reframe(study.group.recals=table(study.group.recalc)) # 22 participants classified

p1<-pgs %>% #filter(Subject.ID %in% subjectIDs) %>%
  filter(premier==T) %>%
  pivot_longer(c(carb_score_center, fat_score_center), values_to="score", names_to="macro") %>%
  ggplot(aes(x=macro, y=score, color=as.factor(study.group.recalc), group=Subject.ID)) + 
  geom_hline(yintercept = 0) +
  geom_point(size=3) + geom_line() +
  scale_color_manual(values=c(palettes$NatMainAccent3$Level2, "black"), name="Study Group")
p1 %>% ggsave(filename="../output/fig_studygrouprecalc_20250212.pdf", height=5, width = 4.5)

genetic_category <- pgs %>% filter(premier==T) %>% select(Subject.ID, study.group.recalc) 
genetic_category %>% write.csv("../data/processed/genetic_category_20250212.csv")

p2 <- pgs %>% filter(premier==T) %>%
  ggplot(aes(x=carb_score_center, y=fat_score_center, color=as.factor(study.group.recalc))) + 
  geom_point(size=3) + 
  scale_color_manual(values=c(palettes$NatMainAccent3$Level2, "black"), name="Study Group") + 
  geom_hline(yintercept = 0, linewidth=0.55) + geom_vline(xintercept = 0, linewidth=0.55)
p2 %>% ggsave(filename="../output/fig_studygrouprecalc_quadrant_20250212.pdf", height=5, width = 4.5) 
  

  
