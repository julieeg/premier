## Model Building

analysis <- readRDS("../data/processed/premier_analysis.rda")
model_building <- analysis %>%
  filter(!is.na(genetic_category.lab))

# ==============================================================================
## Check relevance of genetic PC variables to glucose & insulin 120-min iAUC
# ==============================================================================

#model_building <- model_building %>%
 # mutate(bmi_cat = case_when(bmi<25~"lt25", bmi>=25 & bmi<30~"gt25_lt30", bmi>=30~"gt30"))

covariates <- list(age="Age, years", Sex="Sex, female", bmi="BMI, kg/m2", #bmi_cat="BMI Category", 
                whr="Waist-Hip", PC1="Genetic PC1", PC2="Genetic PC2", PC3="Genetic PC3", 
                PC4="Genetic PC4", PC5="Genetic PC5")


## Test for covariate associations with primary exposure

cov_exp <- do.call(rbind.data.frame, lapply(1:length(covariates), function(i) {
  dat <- model_building %>% select(cov=names(covariates)[[i]], genetic_category.lab)
  sum <- as.data.frame(matrix(NA, 1, 3, dimnames = list(NULL, c("covar", "estimate", "p"))))
  if(is.numeric(dat$cov)==T) {
    t_test <- t.test(cov~genetic_category.lab, data=dat)
    sum[,1:3] <- c(names(covariates)[[i]], t_test$statistic, p=t_test$p.value)
  } else {
    chisq_test <- fisher.test(dat$genetic_category.lab, dat$cov)
    sum[,1:3] <- c(names(covariates)[[i]], chisq_test$estimate, chisq_test$p.value)
  } ; return(sum)
})) ; View(cov_exp)
  
##--> No covariates significantly associated with primary exposure of interest

  

## Test for covariate associations with fasting metabolites

cov_out_fasting <- do.call(rbind.data.frame, lapply(c("glucose", "insulin"), function(m){
  do.call(rbind.data.frame, lapply(1:length(covariates), function(i) {
  mod <- summary(lm(formula(paste0(m, "~", names(covariates)[[i]])), data=model_building))
  coefs <-mod$coef[2,] ; r2 <- mod$r.squared ; f<-mod$fstatistic[1]
  sum<-as.data.frame(matrix(c(coefs, f, r2),1,6, dimnames=list(NULL, c("beta", "se", "t_val", "p_val", "f_val", "r2"))))
  sum$covar=names(covariates)[[i]] ; sum$metabolite=m
  return(sum) })) 
})) ; cov_out_fasting


## Test for covariate associations with primary outcomes 

cov_out_120iauc <- do.call(rbind.data.frame, lapply(c("glucose", "insulin"), function(m){
  do.call(rbind.data.frame, lapply(1:length(covariates), function(i) {
    mod <- summary(lm(formula(paste0(m, "_120iauc~", names(covariates)[[i]])), data=model_building))
    coefs <-mod$coef[2,] ; r2 <- mod$r.squared ; f<-mod$fstatistic[1]
    sum<-as.data.frame(matrix(c(coefs, f, r2),1,6, dimnames=list(NULL, c("beta", "se", "t_val", "p_val", "f_val", "r2"))))
    sum$covar=names(covariates)[[i]] ; sum$metabolite=m
    return(sum) })) 
})) ; cov_out_120iauc



# ================================
## Dig deeper into genetic PCs
# ================================
genetic_pcs <- c(paste0("PC",1:5))

# genetic PCs and genetic category
sum <- as.data.frame(matrix(NA, 5, 5, dimnames = list(NULL, c("covar", "mean_prefer_carb", "mean_prefer_fat", "t_stat", "p"))))
for(i in 1:5) {
  dat <- model_building %>% select(pc=paste0("PC",i), genetic_category.lab)
  t_test <- t.test(pc~genetic_category.lab, data=dat)
  sum[i,1:5] <- c(paste0("PC",i), t_test$estimate[1], t_test$estimate[2], t_test$statistic, p=t_test$p.value)
} ; sum


# Univariate model r2 for genetic PCs predicting fasting metabolites & iAUC values
ggarrange(
  cov_out_fasting %>% filter(covar %in% c(paste0("PC",1:5))) %>%
    ggplot(aes(x=covar, y=r2, fill=covar)) + 
    geom_bar(stat="identity") + 
    facet_grid(~metabolite)+
    scale_fill_manual(values=palettes$NatExt$Purples) + 
    ylab("Univeriate R-squared\nfasting metabolites") + xlab("Genetic PC") + 
    theme_bw(), 
  cov_out_120iauc %>% filter(covar %in% c(paste0("PC",1:5))) %>%
    ggplot(aes(x=covar, y=r2, fill=covar)) + 
    geom_bar(stat="identity") + 
    facet_grid(~metabolite)+
    scale_fill_manual(values=palettes$NatExt$Purples) + 
    ylab("Univariate R-squared\n120-min iAUC") + xlab("Genetic PC") + 
    theme_bw(),
  nrow=2, align="h", common.legend = T, legend = "none"
)


pcs_metabolites <- do.call(rbind.data.frame, lapply(c("fasting", "iauc"), function(type) {
  do.call(rbind.data.frame, lapply(c("glucose", "insulin"), function(m){
  metabolite_type <- ifelse(type == "fasting", m, paste0(m, "_120iauc"))
  #base model  
  mod_base <- summary(lm(formula(paste0(metabolite_type, "~genetic_category.lab")), data=model_building))
  coefs_base <- mod_base$coef[2,1:2] ; r2_base <- mod_base$r.squared ; p_base=mod_base$coef[2,4]
  sum_base <- as.data.frame(matrix(c(coefs_base, p_base, r2_base),1,4, dimnames=list(NULL, c("beta","se","p_val","r2")))) %>%
    mutate(model="base", metabolite=m, type=type, covar="No adjustment", .before=beta)
  #adj model
  sum_adj <- do.call(rbind.data.frame, lapply(1:length(genetic_pcs), function(i) {
    mod_adj <- summary(lm(formula(paste0(metabolite_type, "~genetic_category.lab+", genetic_pcs[i])), data=model_building))
    coefs_adj <- mod_adj$coef[2,1:2] ; r2_adj <- mod_adj$r.squared ; p_adj=mod_adj$coef[2,4]
    sum_adj<-as.data.frame(matrix(c(coefs_adj, p_adj, r2_adj),1,4, dimnames=list(NULL, c("beta", "se", "p_val", "r2")))) %>%
      mutate(model="adj", metabolite=m, type=type, covar=genetic_pcs[i], .before=beta) }))
  bind_rows(sum_base, sum_adj)
  }))
})) 


# Adjusted model r2 for genetic PCs predicting fasting metabolites & iAUC values
ggarrange(
  pcs_metabolites %>% 
    ggplot(aes(x=covar, y=r2, fill=covar)) + 
    facet_wrap(~type+metabolite, nrow = 2)+
    geom_bar(stat="identity") + 
    scale_fill_manual(values=palettes$NatExt$Purples) + 
    ylab("Univeriate R-squared\nfasting metabolites") + xlab("Genetic PC") + 
    theme_bw()
  )

## Adding genetic PCs to prior model
genetic_pcs_add <- c("PC1", "PC1+PC2", "PC1+PC2+PC3", "PC1+PC2+PC3+PC4", "PC1+PC2+PC3+PC4+PC5")

pcs_add_metabolites <- do.call(rbind.data.frame, lapply(c("fasting", "iauc"), function(type) {
  do.call(rbind.data.frame, lapply(c("glucose", "insulin"), function(m){
    metabolite_type <- ifelse(type == "fasting", m, paste0(m, "_120iauc"))
    #base model  
    mod_base <- summary(lm(formula(paste0(metabolite_type, "~genetic_category.lab")), data=model_building))
    coefs_base <- mod_base$coef[2,1:2] ; r2_base <- mod_base$r.squared ; p_base=mod_base$coef[2,4]
    sum_base <- as.data.frame(matrix(c(coefs_base, p_base, r2_base),1,4, dimnames=list(NULL, c("beta","se","p_val","r2")))) %>%
      mutate(model="base", metabolite=m, type=type, covar="No adjustment", .before=beta)
    #adj model
    sum_adj <- do.call(rbind.data.frame, lapply(1:length(genetic_pcs_add), function(i) {
      mod_adj <- summary(lm(formula(paste0(metabolite_type, "~genetic_category.lab+", genetic_pcs_add[i])), data=model_building))
      coefs_adj <- mod_adj$coef[2,1:2] ; r2_adj <- mod_adj$r.squared ; p_adj=mod_adj$coef[2,4]
      sum_adj<-as.data.frame(matrix(c(coefs_adj, p_adj, r2_adj),1,4, dimnames=list(NULL, c("beta", "se", "p_val", "r2")))) %>%
        mutate(model="adj", metabolite=m, type=type, covar=genetic_pcs_add[i], .before=beta) }))
    bind_rows(sum_base, sum_adj)
  }))
})) 


ggarrange(
  pcs_add_metabolites %>% 
    ggplot(aes(x=covar, y=r2, fill=covar)) + 
    facet_wrap(~type+metabolite, nrow = 2)+
    geom_bar(stat="identity") + 
    scale_fill_manual(values=palettes$NatExt$Purples) + 
    ylab("Univeriate R-squared\nfasting metabolites") + xlab("Genetic PC") + 
    theme_bw()
)










