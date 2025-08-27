## Rscript to combine chromosome-stratified pgs


#args=commandArgs(trailingOnly=T)
#macro=args[1]

pgs_dir="../data/processed/pgs"


library(tidyverse);library(data.table)


macros.l=list()
scores.l=list()

macros=c("carb","fat","pro")


## For each macronutrient
for (macro in c("carb", "fat", "pro")) {

# Gather all chr.sscore files
for (i in 1:22) {
score=fread(paste0(pgs_dir,"/", macro,"_prs_chr",i,".sscore"))
names(score)=c("IID", paste0("SCORE_",i))
scores.l[[i]]=score
} 

# Combine into single macro_sum dataframe
m=which(macros == macro)

macros.l[[m]]=scores.l %>% reduce(inner_join,by="IID") %>% 
mutate(SUM=rowSums(across(starts_with("SCORE")))) %>% 
mutate(Subject.ID=gsub(".*[_]", "", gsub(".*[-]", "", IID))) %>% 
rename_at("SUM", ~gsub("SUM", paste0(macro,"_score_sum"), .)) %>%
select(Subject.ID, ends_with("sum")) 

cat(paste0("Done aggregating ", macro, " pgs sscore files.\n"))

}

## Reduce list of macronutrient dataframes into single dataframe
macros_scores=macros.l %>% reduce(inner_join, by = "Subject.ID")

macros_scores %>% fwrite(paste0(pgs_dir,"/macros_pgs_score.txt"), col.names=T)
head(macros_scores)

cat("Finished part 2: Combining macronutrient pgs scores into single macro_pgs_sum file!\n")

#EOF

