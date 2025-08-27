# Prepare metabolomics data following manning-lab-metabolomics-qc pipeline

# ====================================================
## 00_install_pkgs.R
# ====================================================

install.packages(c("dplyr","GGally","BiocManager","matrixStats","ggplot2","data.table","tibble","tidyverse","janitor","gridExtra","circlize","magick"), repos="https://cloud.r-project.org")
BiocManager::install(c("ComplexHeatmap","limma","swamp"),force=T)


# ====================================================
## 01_preprocess_mets.R
# ====================================================

invisible(sapply(c("data.table","dplyr","tidyverse"), library, character.only=T))

options(stringsAsFactors=FALSE)
res_path     <- "../data/processed/"
data_path    <- "../data/raw/metabolomics/"
dir.create(file.path(paste0(res_path,"/formatted_data/")), recursive=T)


# Functions ---------------------------------------------------------------
`%!in%` <- Negate(`%in%`)

rmDtRows <- function(DT, del.idxs) {
  keep.idxs <- setdiff(DT[, .I], del.idxs); # row indexes to keep
  cols = names(DT);
  DT.subset <- data.table(DT[[1]][keep.idxs]); # this is the subsetted table
  setnames(DT.subset, cols[1]);
  for (col in cols[2:length(cols)]) {
    DT.subset[, (col) := DT[[col]][keep.idxs]];
    DT[, (col) := NULL];  # delete
  }
  return(DT.subset);
}

keep_first <- function(dt, i, IDs) {
  IDs <- IDs[grep("HMDB0.*", IDs)]
  for (ID in IDs){
    inds <- which(dt[[i]][,hmdb_id_cols[i],with=F]==ID); dt[[i]][inds,]
    dt[[i]] <- rmDtRows(dt[[i]], inds[2]) # Just keep the first entry
  }
  return(dt[[i]])
}

selectMetInfoCols  <- function(dts, col_inds) {
  lapply(seq_along(dts), function(i) {
    if(is.na(col_inds[i])) { rep(NA,times=nrow(dts[[i]])-data_start_rows[i]+1) }
    else                   { dts[[i]][data_start_rows[i]:nrow(dts[[i]]), col_inds[i], with=F] }
  })
}

renameMethods <- function(method_vec) {
  ifelse(method_vec=="C8-pos", "cp",
         ifelse(method_vec=="C18-neg", "cn",
                ifelse(method_vec=="HIL-neg" | method_vec=="HILIC-neg", "hn",
                       ifelse(method_vec=="HIL-pos" | method_vec=="HILIC-pos", "hp", NA))))
}

addMethodSuffix <- function(compound_id_vec, method_vec) {
  mapply(compound_id_vec, method_vec, FUN=function(id,method) {
    if(is.na(id)) return(NA)
    if(is.na(method)) return(id)
    
    if(method=="C8-pos") return(paste0(id,"_cp"))
    if(method=="C18-neg") return(paste0(id,"_cn"))
    if(method=="HIL-pos") return(paste0(id,"_hp"))
    if(method=="HILIC-pos") return(paste0(id,"_hp"))
    if(method=="HILIC-neg") return(paste0(id,"_hn"))
    if(method=="HIL-neg") return(paste0(id,"_hn"))
  })
}

# Set up ------------------------------------------------------------------

# Read in data
data_paths <- list.files(data_path, pattern = ".csv", full.names = T)
dts <- lapply(data_paths, fread)
names(dts) <- gsub("\\..*", "", basename(data_paths))

# Data row start
data_start_rows <- unlist(lapply(dts, function(x) which(x[,1]!="")[1]+1))
data_start_cols <- unlist(lapply(seq_along(dts), function(x) which(as.numeric(dts[[x]][data_start_rows[x],])%%1==0 &
                                                                     as.numeric(dts[[x]][data_start_rows[x],]) > 2)[1]))

# Sample-related stuff all starts @ data_start_COL
sample_id_rows <- unlist(lapply(dts, function(x) which(x[,1]!="")[1]))
extr_date_rows <- rep(c(3), length(dts))

# Metabolite-related stuff all starts @ data_start_ROW
method_cols      <- unlist(lapply(dts, function(x) which(x == "Method"|x == "method", arr.ind = TRUE)[2]))
compound_id_cols <- unlist(lapply(dts, function(x) which(grepl("compound", x, ignore.case = T), arr.ind = TRUE)[1]))
mrm_cols <- unlist(lapply(dts, function(x) which(grepl("mrm", x, ignore.case = T), arr.ind = TRUE)[1]))
mz_cols  <- unlist(lapply(dts, function(x) which(grepl("mz", x, ignore.case = T), arr.ind = TRUE)[1]))
rt_cols  <- unlist(lapply(dts, function(x) which(grepl("rt", x, ignore.case = T), arr.ind = TRUE)[1]))
hmdb_id_cols  <- unlist(lapply(dts, function(x) which(grepl("hmdb", x, ignore.case = T), arr.ind = TRUE)[1]))
met_name_cols <- unlist(lapply(dts, function(x) which(grepl("metabolite", x, ignore.case = T), arr.ind = TRUE)[1]))

# Clean up data ------------------------------------------------------------------

for(dti in 1:length(dts)) {
  hmdbs_col <- unlist(dts[[dti]][,hmdb_id_cols[dti],with=F])
  hmdbs_col <- sub("\\*","",hmdbs_col)
  for(r in 1:nrow(dts[[dti]]))
    set(dts[[dti]], i=r, j=as.integer(hmdb_id_cols[dti]), hmdbs_col[r])
}

# Duplicated metabolites
ID_list <- lapply(seq_along(dts), function(i) {unname(unlist(unique(dts[[i]][,hmdb_id_cols[i], with=F])))})
dts <- lapply(seq_along(dts), function(i) {keep_first(dts, i, ID_list[[i]])})

dts_clean <- list()
sample_idss <- lapply(seq_along(dts), function(i) unlist(dts[[i]][ sample_id_rows[i], data_start_cols[i]:ncol(dts[[i]]) ]))

# Use compound id w/ method suffix, but for amines must use HMDB id.
met_idss <- lapply(seq_along(dts), function(i) {
  tmp <- if(is.na(compound_id_cols[i])) {unlist( dts[[i]][data_start_rows[i]:nrow(dts[[i]]), hmdb_id_cols[i]    , with=F] )}
  else                    {unlist( dts[[i]][data_start_rows[i]:nrow(dts[[i]]), compound_id_cols[i], with=F] )}
  tmp[duplicated(tmp)] <- paste0(tmp[duplicated(tmp)], "_dup") # If an ID is repeated
  return(tmp)
})

for(i in seq_along(dts)) {
  dts_clean[[i]] <- dts[[i]][data_start_rows[i]:nrow(dts[[i]]),
                             data_start_cols[i]:ncol(dts[[i]])]
  
  # Now data is numbers only, convert to numeric to save memory
  # (but not integer, b/c some counts are > MAX_INT and also amine measurements are floats)
  for(col in 1:ncol(dts_clean[[i]])) set(dts_clean[[i]], j=col, value=as.numeric(dts_clean[[i]][[col]]))
  
  # tranform, and naming
  dts_clean[[i]] <- t(dts_clean[[i]])
  
  colnames(dts_clean[[i]]) <- met_idss[[i]]
  dts_clean[[i]] <- data.table(sample_id=sample_idss[[i]], dts_clean[[i]])
  dts_clean[[i]] <- dts_clean[[i]][!is.na(sample_id)]
}

names(dts) <- gsub("\\..*", "", basename(data_paths))
names(dts_clean) <- names(dts)


# Write cleaned data ------------------------------------------------------

for(i in seq_along(dts_clean)) { write.table(dts_clean[[i]], paste0(res_path, "/formatted_data/", names(dts_clean)[i],"_formatted.txt"), row.names=F) }

# Get metabolite and sample info ------------------------------------------------------------------

extr_datess <- lapply(seq_along(dts), function(i) unlist(dts[[i]][ extr_date_rows[i], data_start_cols[i]:ncol(dts[[i]]) ]))
sample_idss <- lapply(seq_along(dts), function(i) unlist(dts[[i]][ sample_id_rows[i], data_start_cols[i]:ncol(dts[[i]]) ]))
unique_ids <- unique(unlist(sample_idss))

# sample id & cohort cols
sample_info <- data.table(matrix( "", nrow=length(unique_ids), ncol=length(dts)+1 ))
colnames(sample_info) <- c("sample_id", paste0(names(dts),"_batch"))
sample_info$sample_id <- unique_ids

for(i in seq_along(extr_datess)) {
  for(j in seq_along(extr_datess[[i]])) {
    sample_info[sample_id==sample_idss[[i]][j], i+1] <- extr_datess[[i]][j]
  }}

sample_info[, c(2:(length(dts)+1) ):= lapply(.SD, as.factor), .SDcols=2:(length(dts)+1)] # Convert batch columns to factors

sample_info[, is_control := !grepl("MCD",sample_id)]
sample_info <- sample_info[!is.na(sample_info$sample_id)]

met_info <- data.table(
  Compound_Id = unlist(selectMetInfoCols(dts, compound_id_cols)),
  HMDB_Id     = unlist(selectMetInfoCols(dts, hmdb_id_cols)),
  Name        = unlist(selectMetInfoCols(dts, met_name_cols)),
  MRM         = unlist(selectMetInfoCols(dts, mrm_cols)),
  MZ          = as.numeric(unlist(selectMetInfoCols(dts,mz_cols))),
  RT          = as.numeric(unlist(selectMetInfoCols(dts,rt_cols))),
  Method      = unlist(selectMetInfoCols(dts, method_cols))
)

met_info[met_info==""] <- NA
met_info[grepl("tandard",HMDB_Id), HMDB_Id:=NA]
met_info[met_info=="n/a"] <- NA
met_info <- met_info[!is.na(Compound_Id)]
met_info[nchar(HMDB_Id)!=11, HMDB_Id := sub("HMDB","HMDB00", HMDB_Id)]

# Merging knowns across cohorts where the matches exactly (both HMDB and method). 
met_info$Compound_Id <- addMethodSuffix(met_info$Compound_Id,met_info$Method )
met_info$HMDB_Id <- addMethodSuffix(met_info$HMDB_Id,met_info$Method )

# Change bulky method names in the Method column.
met_info$Method  <- renameMethods(met_info$Method)

# Write sample and metabolite info files -------------------------------------------------------------

write.csv(sample_info %>% as.data.frame(), paste0(res_path, "/sample_info.csv"),  row.names=F)
fwrite(met_info, paste0(res_path,"met_info.csv"))

# Messages ---------------------------------------------------------

print(paste0("Info written to:", res_path, ""))
print(paste0("Data written to:", res_path, "formatted_data/"))


# ====================================================
## 02_qc_mets.R
# ====================================================

invisible(sapply(c("GGally","matrixStats","ggplot2","data.table","tibble","tidyverse","janitor","gridExtra","swamp"), library, character.only=T))

options(stringsAsFactors=FALSE)
res_path     <- "../data/processed/no_batch_adj/plasma/"
dir.create(file.path(paste0(res_path,"/QCd_data/")), recursive=T)


# Functions ---------------------------------------------------------------
stat.mode <- function(x) { names(sort(-table(x)))[1] } # Statistical mode. Tabulate then sort then pick the first item.

t2 <- function(x) {
  rn <- rownames(x); cn <- colnames(x)
  x <- t(x)
  rownames(x) <- cn; colnames(x) <- rn
  return(x)
}

#Limited to categorical variables for now
PCAPlot <- function(df, color_var, n_PC, title="") {
  pca <- prcomp(df, center=T, scale.=T)
  
  var_explained <- scales::percent(summary(pca)$importance[2,], 0.1)
  for(i in 1:n_PC) colnames(pca$x)[i] <- paste0("PC",i," (",var_explained[i],")")
  
  ggpairs(data = data.frame(pca$x[,1:n_PC], check.names=F),
          lower = list(continuous = wrap("points", alpha=0.5, size=0.5, pch=1)),
          diag  = list(continuous = wrap("densityDiag", alpha=0.5)),
          #upper= list(continuous = wrap("cor"),
          axisLabels = "none",
          mapping = aes(color = color_var),
          legend  = c(2,1) ) + labs(color="legend_title") + ggtitle(title)
}
rm0VarRows <- function(df) {
  rows2keep <- rowVars(df, na.rm=T) > 0
  rows2keep[is.na(rows2keep)] <- FALSE
  print(paste(nrow(df)-sum(rows2keep),"/",nrow(df),"rows will be removed for having 0 variance"))
  df[rows2keep,]
}

rmHighMissingness <- function(df, thresh, verbosity=1) {
  rows2keep <- rowSums(is.na(df)) < thresh*ncol(df)
  if(verbosity>0) print(paste( sum(!rows2keep),"/",nrow(df),"rows will be removed for >",scales::percent(thresh),"missingness" ))
  if(verbosity>2 & sum(!rows2keep>0)) print(paste( "Rows removed:",names(df)[!rows2keep] ))
  df <- df[rows2keep,]
  
  cols2keep <- colSums(is.na(df)) < thresh*nrow(df)
  if(verbosity>0) print(paste( sum(!cols2keep),"/",ncol(df),"cols will be removed for >",scales::percent(thresh),"missingness" ))
  if(verbosity>1 & sum(!cols2keep>0)) print(paste( "Cols removed:",names(df)[!cols2keep] ))
  df <- df[,cols2keep]
}

winsorizeBySd <- function(df, n_sds) {
  t2(apply(df,1, function(r) {
    u <- mean(r) + sd(r)*n_sds
    l <- mean(r) - sd(r)*n_sds
    r <- sapply(r, function(x) {
      if(x > u) {x <- u}
      if(x < l) {x <- l}
      else      {x}})}))
}

batchAdj <- function(df, batch_factor) { # Median scaling. batch_factor should have same cols as df.
  if(length(levels(batch_factor)) < 2) { print("Need at least 2 batch levels to adjust for batch, returning df untouched."); return(df) }
  swamp::quickadjust.ref(df, batches=batch_factor, refbatch=stat.mode(batch_factor))$adjusted.data
}

zScoreRows <- function(df) {
  t2(apply(df,1, function(r) {
    m <- mean(r)
    s <- sd(r)
    r <- sapply(r, function(x) (x-m)/s)
  }))}

plotSkewKurt <- function(df, title) {
  n <- ncol(df)
  ms <- rowMeans(df)
  sds <- rowSds(df)
  skews <- sapply(1:nrow(df), function(r) (sum((df[r,]-ms[r])^3)/sds[r]^3)/n   )
  kurts <- sapply(1:nrow(df), function(r) (sum((df[r,]-ms[r])^4)/sds[r]^4)/n-3 )
  sk <- data.frame(skews=skews, kurts=kurts)
  
  ggplot(sk, aes(x=skews, y=kurts)) +
    geom_point(size=0.5) + ggtitle(title) +
    geom_vline(xintercept = -0.5, linetype="dotted", color="blue", linewidth=1) +
    geom_vline(xintercept =  0.5, linetype="dotted", color="blue", linewidth=1) + 
    geom_hline(yintercept = -2.0, linetype="dotted", color="blue", linewidth=1) +
    geom_hline(yintercept =  2.0, linetype="dotted", color="blue", linewidth=1)
}

# Read in data ------------------------------------------------------------
data_files <- list.files(path=paste0("../data/processed/formatted_data/"), pattern=".txt", full.names = T)
dict_files <- list.files(path="../data/processed/", pattern="info.csv", full.names = T)
sample_info <- fread(dict_files[2])

# Format ------------------------------------------------------------------

dfs <- lapply(data_files, function(filename) {
  df <- as.matrix(fread(filename))
  rownames(df) <- df[,"sample_id"]; df <- df[,-1]
  mode(df) <- "numeric"
  df <- t2(df)
  non_control <- colnames(df) %in% sample_info$sample_id[!sample_info$is_control]
  if(sum(non_control)>0) df <- df[,non_control] # only non-control samples
  return(df)
})


batch_info_cols <- grep("batch", colnames(sample_info), value = T)
batch_list <- lapply(1:length((batch_info_cols)), function(i) {
  df <- dfs[[i]]
  batch <- sample_info[sample_id %in% colnames(df), sample_id, eval(batch_info_cols[i])]
  batch <- batch[complete.cases(batch),]
  batch <- batch[match(colnames(df),sample_id),] # Reorder batch to match the order of samples
  batch <- factor(unlist(batch[,1])) # Need only the dates column
})

df_labels <- gsub("_batch", "", batch_info_cols)
names(dfs) <- df_labels
names(batch_list) <- names(dfs)

# QC -----------------------------------------------------------------------

tmp <- mapply(dfs,names(dfs),batch_list, SIMPLIFY=F, FUN=function(df,nm,batch_factor) {
  print(nm)
  df <- rm0VarRows(df) # 1
  df <- rmHighMissingness(df, 0.25) # 2
  df <- t2(apply(df,1, function(r) { r[is.na(r)] <- min(r,na.rm=T)/2; r })) # 3
  
  df_l2 <- log2(df+1) # 5a
  df_l2_z <- zScoreRows(log2(df+1)) # 5b
  df_inv_norm <- t2(apply(df,1, function(r) qnorm( (rank(r)-0.5)/length(r) ) )) # 5c
  df_ln <- log(df+1) # 5d
  df_ln_z <- zScoreRows(log(df+1))
  
  dfs <- list(default=df, l2=df_l2, l2_z=df_l2_z, inv_norm=df_inv_norm, ln=df_ln, ln_z=df_ln_z)
  
  dfs <- lapply(dfs, winsorizeBySd, 5) # 4
  dfs <- lapply(dfs, batchAdj, batch_factor) # Comment out for non-adjusted QC'd data for PCA plots
})


# Transform ---------------------------------------------------------------

# tmp is a list of QC transforms (l2, l2_z, etc.) per cohort. Rearrange to List of cohorts per transform. 
dfss <- list(default  = lapply(tmp, function(cohort) cohort[["default" ]]),
             l2       = lapply(tmp, function(cohort) cohort[[   "l2"   ]]),
             l2_z     = lapply(tmp, function(cohort) cohort[[  "l2_z"  ]]),
             inv_norm = lapply(tmp, function(cohort) cohort[["inv_norm"]]),
             ln       = lapply(tmp, function(cohort) cohort[[   "ln"   ]]),
             ln_z     = lapply(tmp, function(cohort) cohort[[  "ln_z"  ]]))

# Figures -----------------------------------------------------------------
outfile <- paste0(res_path, "/QCd_data/QC_figures.pdf")
pdf(outfile)
for(i in seq_along(dfss[1])) {
  plot_data_column = function (col) {
    dfss[[col]][[i]] %>% as.data.frame() %>% mutate(mean=rowMeans(., na.rm = TRUE)) %>% 
      ggplot(aes(x=mean)) + 
      geom_histogram(color="black", fill="white", bins = 30) +
      ggtitle(paste(names(dfss[[col]])[i],names(dfss)[col],"distribution")) +
      xlab(names(dfss)[col])
  }
  myplots <- lapply(seq_along(dfss), plot_data_column)
  do.call("grid.arrange", c(myplots, ncol=2))
}

for(i in seq_along(dfss[1])) {
  plot_data_column = function (col) {
    plotSkewKurt(dfss[[col]][[i]], 
                 title=paste(names(dfss[[col]])[i],names(dfss)[col]))
  }
  myplots <- lapply(seq_along(dfss), plot_data_column)
  do.call("grid.arrange", c(myplots, ncol=2))
}

for(i in seq_along(dfss)) {
  for(j in seq_along(dfss[[i]])) {
    df <- as.data.frame(dfss[[i]][[j]])
    df_t <- t(df) %>% as.data.frame
    data <- rownames_to_column(df_t, var = "Sample_Id") %>% 
      mutate(group = ifelse(str_detect(Sample_Id,'QC'),1,0))
    
    plot_data_column = function (data, column) {
      ggplot(data, aes(x = 1:nrow(data), y = get(column), col=group))+
        geom_point() + ylab("sample") + xlab(column) + 
        ggtitle(paste(names(dfss[[i]])[j],names(dfss)[i], ":", column)) + 
        theme_bw() + theme(legend.position="none")
    }
    myplots <- lapply(colnames(data)[sample(c(2:length(colnames(data))),4)], plot_data_column, data = data)
    do.call("grid.arrange", c(myplots, ncol=2))
  }}
garbage <- dev.off()

# Write files -------------------------------------------------------------

for(type in names(dfss)) {dir.create(file.path(paste0(res_path, "/QCd_data/",type)), showWarnings=F)}

for(i in seq_along(dfss)) {
  for(j in seq_along(dfss[[i]])) {
    write.csv(dfss[[i]][[j]], paste0(res_path, "/QCd_data/", names(dfss)[i],"/", names(dfss[[i]])[j],"_QCd_",names(dfss)[i],".csv"))
  }}

# Message -----------------------------------------------------------------
print(paste("QC written to:", paste0(res_path, "/QCd_data/")))



# ====================================================
## 03_merge_mets.R 
## NOTE: JEG changed "Read data over" to grab ln_z metabolomics files**
# ====================================================

library(dplyr)
library(tidyr)

##Read data over
cn <- read.csv("../data/processed/no_batch_adj/plasma/QCd_data/ln/23_0525_MER_PREMIER_C18-neg_QCd_ln.csv")
cp <- read.csv("../data/processed/no_batch_adj/plasma/QCd_data/ln/23_0525_MER_PREMIER_C18-pos_QCd_ln.csv")
hn <- read.csv("../data/processed/no_batch_adj/plasma/QCd_data/ln/23_0525_MER_PREMIER_HILIC-pos_QCd_ln.csv")
hp <- read.csv("../data/processed/no_batch_adj/plasma/QCd_data/ln/23_0525_MER_PREMIER_C18-neg_QCd_ln.csv")
output_fnm <-  "../data/processed/no_batch_adj/plasma/merged_data/ln_merged_QCd_plasma_knowns.csv"
met_info <- fread("../data/processed/met_info.csv",na.strings="")
dir.create(dirname(output_fnm))

names(cn)[1] <- "Compound_Id"
names(cp)[1] <- "Compound_Id"
names(hn)[1] <- "Compound_Id"
names(hp)[1] <- "Compound_Id"

# Deal with muscle "_repeat" samples in muscle. Remove duplicates (preferring the samples with _repeat), then rename the ids to remove "_repeat".
fix_repeats <- function(df) {
  id_no_repeat <- sub("_repeat","",colnames(df))
  cols2keep <- grepl("_repeat",colnames(df))   |   !( duplicated(id_no_repeat) | duplicated(id_no_repeat,fromLast=T) )
  df <- df[,cols2keep]
  colnames(df) <- sub("_repeat","",colnames(df))
  df
}

cn <- cn %>%
  mutate(Compound_Id = paste(Compound_Id, "cn", sep = "_")) %>% fix_repeats

cp <- cp %>%
  mutate(Compound_Id = paste(Compound_Id, "cp", sep = "_")) %>% fix_repeats

hn <- hn %>%
  mutate(Compound_Id = paste(Compound_Id, "hn", sep = "_")) %>% fix_repeats

hp <- hp %>%
  mutate(Compound_Id = paste(Compound_Id, "hp", sep = "_")) %>% fix_repeats


## Merging the four methods
# Assuming df1 and df2 are your data frames
merged_df <- bind_rows(cn, cp, hn, hp)

#Selecting only relevant columns for merging
info <- met_info[, c("Compound_Id","Name", "HMDB_Id")]

## Merging with metab info
data <- merge(info, merged_df, by = "Compound_Id" )

#Write the files
#data %>% fwrite("../data/processed/no_batch_adj/plasma/merged_data/inv_norm_merged_QCd_plasma.csv")

## Selecting only known metabolites
selected_rows <- data[!is.na(HMDB_Id) & HMDB_Id!="NA",]

## Remove columns and keeping only HMDB_Id column
df <- subset(selected_rows, select = -c(Compound_Id,Name))

#Remove the controls
df <- df[!grepl("PREF",HMDB_Id)]

# Removing metabolite duplicates
df <- distinct(df, HMDB_Id, .keep_all = TRUE)

# Transposing data samples as rows, metbaolites as columns
transposed_df <- pivot_longer(df, cols = -HMDB_Id, names_to = "sample_id", values_to = "Value")
transposed_df <- pivot_wider(transposed_df, names_from = HMDB_Id, values_from = Value)

fwrite(transposed_df, output_fnm)

## END OF PIPELINE



