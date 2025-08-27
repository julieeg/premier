## Descriptive chatacteristics

## Table 1b - Baseline metabolites
vars_to_summarise.metab <- c(glucose="Glucose, mg/dL", insulin="Insulin", 
                             tg_log="log(Triglyceride)", chol="Total cholesterol", hdl="HDL cholesterol", 
                             ldl="LDL cholesterol", protein="Total protein", alb="Albumin", 
                             globulin="Globulin", alb_globulin_ratio="Albumin:Globulin",
                             tot_bilirubin_log="log(Total Bilirubin)", dir_bilirubin_log="log(Direct Bilirubin)", 
                             ast="AST", alt_log="Log(ALT)")

metabolites_total <- print_summary_table(data=processed.dat, vars_to_summarize = vars_to_summarise.metab, p_adjust = "none") %>%
  bind_cols(
    Range=c("-", processed.dat %>% select(names(vars_to_summarise.metab)) %>% summarise_all(., ~range_formatted(.)) %>% t())
    ) %>% 
  rename(Summary=Total)
metabolites_strat <- print_summary_table(data=processed.dat, vars_to_summarize = vars_to_summarise.metab, p_adjust = "none", var_strata = "pps_cat", p_types = "descriptive")



## Figure 1 - Postprandial meal responses in total cohort (glucose & insulin)
pp1 <- c("", "_30", "_60", "_120", "_180")
vars_to_summarise.pp <-c(paste0(rep("glucose", 5), pp1), paste0(rep("insulin", 5), pp1)) 
names(vars_to_summarise.pp)<-vars_to_summarise.pp

print_summary_table(data=processed.dat, vars_to_summarize = vars_to_summarise.pp, p_adjust = "none")


plot_postprandialxtime <- function(metabolite, times = c(0, 30, 60, 120, 180), fill_color) {
  processed.dat %>% select(id, starts_with(metabolite)) %>%
    rename_with(~gsub(metabolite, "x", .), starts_with(metabolite)) %>%
    rename(x_0=x) %>%
    rename_with(., ~gsub("m", "", gsub("_log", "", .))) %>%
    select(id, paste0("x_", times)) %>%
    pivot_longer(-id, names_sep="_", names_to=c("metabolite", "time")) %>%
    mutate(across(c("time", "value"), ~as.numeric(.))) %>%
    ggplot(aes(x=time, y=value, group=id)) + 
    scale_x_continuous(breaks=times) +
    geom_line(alpha=0.55) + 
    stat_summary(aes(group = metabolite), geom = "point", fun = mean, size = 4, color=fill_color) + 
    stat_summary(aes(group = metabolite), geom = "line", fun = mean, linewidth = 2, color=fill_color) + 
    scale_color_manual(values=fill_color) + scale_fill_manual(values=fill_color) + 
    xlab("Time (minutes)") + ylab("Concentration") + 
    theme_bw() + 
    ggtitle(paste("Postprandial", metabolite, "at", paste0(times, collapse=", "), "min"))
}

plot_postprandialxtime("glucose", times=c(0,30,60,120,180), fill_color = palettes$NatureMainAccents_level2[1])
plot_postprandialxtime("insulin", times=c(0,30,60,120,180), fill_color = palettes$NatureMainAccents_level2[2])


plot_postprandialxtime 

