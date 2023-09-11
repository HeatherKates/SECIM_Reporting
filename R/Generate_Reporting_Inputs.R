Generate_Report_Inputs <-function(client,samples_to_drop=NULL,mzmine_version,ReferenceLevel=NULL,Input,contrast_var,
                                  num_meta,anova_formula=NULL,lm_model=NULL,test_type,subset=NULL,metid_DB_file){
  #client="Dudeja.Serum"
  #samples_to_drop=c("")
  #mzmine_version=2
  #ReferenceLevel="WT"
  #Input <- "Client_Input_Sheets/Dudeja-Metabolomics_Serum.xlsx" 
  #contrast_var="Class"
  
  #For stats subfunction
  #num_meta=1
  #subset=NULL,
  #anova_formula=as.formula(paste("id ~", "Class","+","Error(ID)")),lm_model=NULL,test_type="t.test",
  #metid_DB_file="kegg_ms1_database0.0.3.rda"

library(dplyr)
library(qs)
library(lme4)
library(impute)
library(data.table)
library(foreach)
library(parallel)
library(emmeans)
library(broom)
library(ggplot2)
library(gridExtra)
library(grid)
library(stringr)
library(nlme)
library(readxl)
library(omu)
library(metid)
library(MetaboAnalystR)

source("SECIM_Reporting/R/SECIM_Metabolomics.R")
source("SECIM_Reporting/R/Norm_Plots.R")
source("SECIM_Reporting/R/annotate_metabolites_mass_dataset.R")
source("SECIM_Reporting/R/mzIdentify_mass_dataset.R")

##Mode Neg

##RULES FOR INPUT DATA: 
#No special characters in names or metadata except for "_" and ".". 
#Names of samples are in columns in the peaktable and the string before the "p" or "n" must match the sample names in the metadata 
metadata <- read_excel(Input,sheet="Sample.data")
peakdata <- read_excel(Input,sheet="Peaktable.neg")

#################################################
######Standardize peak table columns#############
#################################################
if(mzmine_version==3){
  colnames(peakdata) <- c("id","rt","mz","internal_ID","spectral_ID",colnames(peakdata[6:ncol(peakdata)]))
  } else
    {
    colnames(peakdata) <- c("id","mz","rt","compound",colnames(peakdata[5:ncol(peakdata)]))
    peakdata <- peakdata %>% dplyr::relocate("rt",.before="mz")
    }
#################################################
######Rules to process two columns of IDs########
#################################################
#This assumes that the first one is named "compounds_db_identity:compound_name" and the second is "spectral_ID"
#NOTE: THE COLUMNS MAY NEED TO BE CONVERTED TO "TEXT" IN EXCEL TO AVOID BEING CONVERTED TO `NA` WHEN READ IN
if("internal_ID" %in% colnames(peakdata)){
#Create a column "confidence" that is empty
  peakdata$Confidence <- ""
  # Given a data frame with columns "A", "B", and "C"
  # Assuming the data frame is named "df"

  # 1) If column "B" is empty and column "A" is not empty, enter value "1" in column "C"
  peakdata$Confidence <- ifelse(is.na(peakdata$`spectral_ID`) & !is.na(peakdata$`internal_ID`), 1, peakdata$Confidence)

  # 1B) If column B == column A, make column B NA
  #peakdata$`spectral_ID` <- ifelse(peakdata$`spectral_ID` == peakdata$`internal_ID`, NA, peakdata$`spectral_ID`)

  # 2) If column "B" is not empty and column "A" is empty, enter value "2" in column "C" and convert the text in column "B" to uppercase
  peakdata$Confidence <- ifelse(!is.na(peakdata$`spectral_ID`) & is.na(peakdata$`internal_ID`), 2, peakdata$Confidence)
  peakdata$`spectral_ID` <- ifelse(!is.na(peakdata$`spectral_ID`) & is.na(peakdata$`internal_ID`), toupper(peakdata$`spectral_ID`), peakdata$`spectral_ID`)

  # 3) If column "B" is not empty and column "A" is not empty, replace the text in column "B" with ""
  peakdata$`spectral_ID` <- ifelse(!is.na(peakdata$`spectral_ID`) & !is.na(peakdata$`internal_ID`), "", peakdata$`spectral_ID`)

  # 4) Combine columns "A" and "B" into a single column named "row identify (main ID)" by populating that column with whichever column "A" or "B" is not empty for each row.
  peakdata$compound_name <- ifelse(!is.na(peakdata$`internal_ID`), peakdata$`internal_ID`, peakdata$`spectral_ID`)

  # Remove the original columns "A" and "B"
  peakdata <- peakdata[, !(names(peakdata) %in% c("internal_ID", "spectral_ID"))]

  #Move the compound_name to the position that is normal
  peakdata <- peakdata %>% relocate(compound_name,.after=mz)

  #Save the confidence information to be added back after stats
  Neg_pre_stats_conf_levels <- peakdata %>% select(`id`,Confidence)

  #Remove the confidence column
  peakdata <- peakdata[, -which(names(peakdata) == "Confidence")]
  }

#################################################
######Steps to process column names #############
#################################################

#If there are duplicate colnames, make them unique and print a warning
# Check for duplicate column names

if (any(duplicated(colnames(peakdata)))) {
  # Get the duplicate column colnames
  duplicate_cols <- colnames(peakdata)[duplicated(colnames(peakdata))]
  
  # Make duplicate column colnames unique
  colnames(peakdata)[5:ncol(peakdata)] <- make.names(colnames(peakdata)[5:ncol(peakdata)], unique = TRUE)

  # Print a warning message
  warning(paste("Duplicate column colnames found and made unique:",
                paste(duplicate_cols, collapse = ", ")))
}
#v10 Replace "-" with "_" in sample names because metaboanalyst will convert one and mess it up
colnames(peakdata) <- gsub("-","_",colnames(peakdata))
metadata$Sample.Name <- gsub("-","_",metadata$Sample.Name)


#Remove any rows with "adduct" in the fourth column
if(mzmine_version==2){
  peakdata <- filter(peakdata,!grepl("adduct",compound))
  #Remove any rows with "Complex" in the fourth column
  peakdata <- filter(peakdata,!grepl("Complex",compound))
  }
#################################################
#################END#############################
#################################################

#################################################
######Steps to process metadata #################
#################################################

#Custom: if needed, add an A to the beginning of reference groups - should phase this out by updating v10 to set reference group/custom contrasts even for pairwise=TRUE
if(!is.null(ReferenceLevel)){
    metadata$Class <- gsub(ReferenceLevel,paste0("Zref_",ReferenceLevel),metadata$Class)
}

if (length(samples_to_drop)>0){
  metadata <- metadata %>% filter(!Sample.Name %in% samples_to_drop)
}

#################################################
#################END#############################
#################################################

#################################################
######Create main data object####################
#################################################
#Create the data object by rbinding the metadata and peak data
#assume there are four cols before samples
# Indices of columns from peakdata that match any Sample.Name from metadata
matching_indices <- unlist(lapply(metadata$Sample.Name, function(x) {
  grep(x, colnames(peakdata))
}))

# Keep only columns 1:4 and the matched columns
peakdata <- peakdata[, c(1:4, matching_indices)]

# Rename the columns after the 4th to only have the Sample.Name
colnames(peakdata)[5:ncol(peakdata)] <- metadata$Sample.Name


peaks <- peakdata[,metadata$Sample.Name]
data <- data.frame(t(rbind(c(colnames(metadata),peakdata$id),merge(metadata,data.frame(t(peaks)),by.x="Sample.Name",by.y=0))))
colnames(data) <- data[1,];data <- data[-1,];data <- data %>% dplyr::rename("id"="Sample.Name");rownames(data) <- 1:nrow(data)
dataset <- data

#################################################
#################END#############################
#################################################


#################################################
######Run the function that does the stats#######
#################################################
#data=data;peakdata=peakdata;num_meta=1;original_data=data;contrast_var="Class";anova_formula=Metabolite~Class+Error(ID);
#test_type="t.test";subset=list(list("GI_MB_TBI","GII_MB_Sham"),list("GI_CX_TBI","GII_CX_Sham"));lm_model=Metabolite~Class+(1|ID);emmeans_var=~Class;mode="Pos";metid_DB_file="kegg_ms1_database0.0.3.rda";client="Hendrich-Wang"
#options for metid_DB are kegg_ms1_database0.0.3.rda,"hmdb_database0.0.3.rda,"bloodexposome_database1.0.rda"
#Example of anova_formula=as.formula(paste("Metabolite ~", "Class","+","Error(ID)"))
#debug(SECIM_Metabolomics)
neg.output <- SECIM_Metabolomics(dataset=data,peakdata=peakdata,num_meta=num_meta,original_data=data,contrast_var=contrast_var,
  subset=subset,
  anova_formula=anova_formula,lm_model=lm_model,test_type=test_type,emmeans_var=contrast_var,mode="Neg",
  metid_DB_file=metid_DB_file,client=client,metadata=metadata)

#Positive mode

metadata <- read_excel(Input,sheet="Sample.data")

peakdata <- read_excel(Input,sheet="Peaktable.pos")

#################################################
######Standardize peak table columns#############
#################################################
if(mzmine_version==3){
  colnames(peakdata) <- c("id","rt","mz","internal_ID","spectral_ID",colnames(peakdata[6:ncol(peakdata)]))
}else{
  colnames(peakdata) <- c("id","mz","rt","compound",colnames(peakdata[5:ncol(peakdata)]))
  peakdata <- peakdata %>% dplyr::relocate("rt",.before="mz")
}
#################################################
######Rules to process two columns of IDs########
#################################################
#This assumes that the first one is named "compounds_db_identity:compound_name" and the second is "spectral_ID"
#NOTE: THE COLUMNS MAY NEED TO BE CONVERTED TO "TEXT" IN EXCEL TO AVOID BEING CONVERTED TO `NA` WHEN READ IN
if("internal_ID" %in% colnames(peakdata)){
  #Create a column "confidence" that is empty
  peakdata$Confidence <- ""
  # Given a data frame with columns "A", "B", and "C"
  # Assuming the data frame is named "df"
  
  # 1) If column "B" is empty and column "A" is not empty, enter value "1" in column "C"
  peakdata$Confidence <- ifelse(is.na(peakdata$`spectral_ID`) & !is.na(peakdata$`internal_ID`), 1, peakdata$Confidence)
  
  # 1B) If column B == column A, make column B NA
  #peakdata$`spectral_ID` <- ifelse(peakdata$`spectral_ID` == peakdata$`internal_ID`, NA, peakdata$`spectral_ID`)
  
  # 2) If column "B" is not empty and column "A" is empty, enter value "2" in column "C" and convert the text in column "B" to uppercase
  peakdata$Confidence <- ifelse(!is.na(peakdata$`spectral_ID`) & is.na(peakdata$`internal_ID`), 2, peakdata$Confidence)
  peakdata$`spectral_ID` <- ifelse(!is.na(peakdata$`spectral_ID`) & is.na(peakdata$`internal_ID`), toupper(peakdata$`spectral_ID`), peakdata$`spectral_ID`)
  
  # 3) If column "B" is not empty and column "A" is not empty, replace the text in column "B" with ""
  peakdata$`spectral_ID` <- ifelse(!is.na(peakdata$`spectral_ID`) & !is.na(peakdata$`internal_ID`), "", peakdata$`spectral_ID`)
  
  # 4) Combine columns "A" and "B" into a single column named "row identify (main ID)" by populating that column with whichever column "A" or "B" is not empty for each row.
  peakdata$compound_name <- ifelse(!is.na(peakdata$`internal_ID`), peakdata$`internal_ID`, peakdata$`spectral_ID`)
  
  # Remove the original columns "A" and "B"
  peakdata <- peakdata[, !(names(peakdata) %in% c("internal_ID", "spectral_ID"))]
  
  #Move the compound_name to the position that is normal
  peakdata <- peakdata %>% relocate(compound_name,.after=mz)
  
  #Save the confidence information to be added back after stats
  Pos_pre_stats_conf_levels <- peakdata %>% select(`id`,Confidence)
  
  #Remove the confidence column
  peakdata <- peakdata[, -which(names(peakdata) == "Confidence")]
}

#################################################
######Steps to process column names #############
#################################################

#If there are duplicate colnames, make them unique and print a warning
# Check for duplicate column names

if (any(duplicated(colnames(peakdata)))) {
  # Get the duplicate column colnames
  duplicate_cols <- colnames(peakdata)[duplicated(colnames(peakdata))]
  
  # Make duplicate column colnames unique
  colnames(peakdata)[5:ncol(peakdata)] <- make.names(colnames(peakdata)[5:ncol(peakdata)], unique = TRUE)
  
  # Print a warning message
  warning(paste("Duplicate column colnames found and made unique:",
                paste(duplicate_cols, collapse = ", ")))
}

#v10 Replace "-" with "_" in sample names because metaboanalyst will convert one and mess it up
colnames(peakdata) <- gsub("-","_",colnames(peakdata))
metadata$Sample.Name <- gsub("-", "_",metadata$Sample.Name)

#Remove any rows with "adduct" in the fourth column
if(mzmine_version==2){
  peakdata <- filter(peakdata,!grepl("adduct",compound))
  #Remove any rows with "Complex" in the fourth column
  peakdata <- filter(peakdata,!grepl("Complex",compound))
  }
#################################################
#################END#############################
#################################################

#################################################
######Steps to process metadata #################
#################################################

#Custom: if needed, add an A to the beginning of reference groups - should phase this out by updating v10 to set reference group/custom contrasts even for pairwise=TRUE
if(!is.null(ReferenceLevel)){
  metadata$Class <- gsub(ReferenceLevel,paste0("Zref_",ReferenceLevel),metadata$Class)
  }

if (length(samples_to_drop)>0){
  metadata <- metadata %>% filter(!Sample.Name %in% samples_to_drop)
  }

#################################################
#################END#############################
#################################################

#################################################
######Create main data object####################
#################################################
#Create the data object by rbinding the metadata and peak data
#assume there are four cols before samples

# Indices of columns from peakdata that match any Sample.Name from metadata
matching_indices <- unlist(lapply(metadata$Sample.Name, function(x) {
  grep(x, colnames(peakdata))
}))

# Keep only columns 1:4 and the matched columns
peakdata <- peakdata[, c(1:4, matching_indices)]

# Rename the columns after the 4th to only have the Sample.Name
colnames(peakdata)[5:ncol(peakdata)] <- metadata$Sample.Name

peaks <- peakdata[,metadata$Sample.Name]
data <- data.frame(t(rbind(c(colnames(metadata),peakdata$id),merge(metadata,data.frame(t(peaks)),by.x="Sample.Name",by.y=0))))
colnames(data) <- data[1,];data <- data[-1,];data <- data %>% dplyr::rename("id"="Sample.Name");rownames(data) <- 1:nrow(data)
dataset <- data

#################################################
#################END#############################
#################################################


#################################################
######Run the function that does the stats#######
#################################################
#data=data;peakdata=peakdata;num_meta=1;original_data=data;contrast_var="Class";anova_formula=Metabolite~Class+Error(ID);
#test_type="t.test";subset=list(list("GI_MB_TBI","GII_MB_Sham"),list("GI_CX_TBI","GII_CX_Sham"));lm_model=Metabolite~Class+(1|ID);emmeans_var=~Class;mode="Pos";metid_DB_file="kegg_ms1_database0.0.3.rda";client="Hendrich-Wang"
#options for metid_DB are kegg_ms1_database0.0.3.rda,"hmdb_database0.0.3.rda,"bloodexposome_database1.0.rda"
#Example of anova_formula=as.formula(paste("Metabolite ~", "Class","+","Error(ID)"))

pos.output <- SECIM_Metabolomics(
  dataset=data,peakdata=peakdata,num_meta=num_meta,original_data=data,contrast_var=contrast_var,
  subset=subset,
  anova_formula=anova_formula,lm_model=lm_model,test_type=test_type,emmeans_var=contrast_var,mode="Pos",
  metid_DB_file=metid_DB_file,client=client,metadata=metadata)

#save.image(paste(client,"functionoutput.RDATA",sep="_"))
#################################################
#################END#############################
#################################################


#Reduce the pos and neg results combined to representative peaks per compound
combined_results <- rbind(neg.output[[1]],pos.output[[1]])
combined_results$temp.lc.Metabolite <- tolower(combined_results$compound)

#Reducing for duplicate KEGG IDs
HP.KEGG.list <- combined_results %>% filter(Level=="1") %>% dplyr::select(KEGG) %>% unlist() %>% na.omit() #KEGG IDs that are in the data with high confidence
combined_results <- combined_results %>% mutate(status=case_when(KEGG %in% HP.KEGG.list~"HP"))#A status of "HP" indicates that this KEGG has been identified with confidence level "1" for some peak in the table
orig.combined_results <- combined_results #save for anova processing

if(test_type=="t.test"){ #t test
  #For cases when 1+ of the duplicate KEGG IDs has a confidence level "1", pick the one with conf level 1 (then sort by p-value)
  temp.HP.Keggs <- combined_results %>% filter(status=="HP") %>% arrange(KEGG,Level,p.value)%>% distinct(KEGG, .keep_all = TRUE)
  #For cases when 1+ of the duplicates does not have a confidence level "1", pick the best match and then the best p-value
  temp.notHP.Keggs <- combined_results %>% filter(is.na(status)) %>% filter(!is.na(KEGG)) %>% arrange(KEGG,plyr::desc(mz.match.score),p.value) %>% distinct(KEGG, .keep_all = TRUE)

  #For peaks with Compound.names but no KEGG IDs
  #For cases when 1+ of the duplicates has a confidence level "1"
  HP.names.list <- combined_results %>% filter(Level=="1") %>% dplyr::select(compound) %>% unlist() %>% na.omit() #compounds that are in the data with high confidence
  combined_results <- combined_results %>% mutate(status=case_when(compound %in% HP.names.list~"HP"))#A status of "HP" indicates that this KEGG has been identified with confidence level "1" for some peak in the table
  temp.HP.names <- combined_results %>% filter(status=="HP") %>% filter(is.na(KEGG)) %>% arrange(temp.lc.Metabolite,Level,p.value)%>% distinct(temp.lc.Metabolite, .keep_all = TRUE)
  #For cases when 1+ of the duplicates does not have a confidence level "1"
  temp.notHP.names <- combined_results %>% filter(is.na(status)) %>% arrange(temp.lc.Metabolite,plyr::desc(mz.match.score),p.value) %>% distinct(temp.lc.Metabolite, .keep_all = TRUE)
  }
if (test_type=="nostats") {
    #For cases when 1+ of the duplicate KEGG IDs has a confidence level "1", pick the one with conf level 1 (then sort by p-value)
    temp.HP.Keggs <- combined_results %>% filter(status=="HP") %>% arrange(KEGG,Level,mz.match.score)%>% distinct(KEGG, .keep_all = TRUE)
    #For cases when 1+ of the duplicates does not have a confidence level "1", pick the best match and then the highest log2FC
    temp.notHP.Keggs <- combined_results %>% filter(is.na(status)) %>% filter(!is.na(KEGG)) %>% arrange(KEGG,plyr::desc(mz.match.score),desc(log2FC)) %>% distinct(KEGG, .keep_all = TRUE)
    
    #For peaks with Compound.names but no KEGG IDs
    #For cases when 1+ of the duplicates has a confidence level "1"
    HP.names.list <- combined_results %>% filter(Level=="1") %>% dplyr::select(compound) %>% unlist() %>% na.omit() #compounds that are in the data with high confidence
    combined_results <- combined_results %>% mutate(status=case_when(compound %in% HP.names.list~"HP"))#A status of "HP" indicates that this KEGG has been identified with confidence level "1" for some peak in the table
    temp.HP.names <- combined_results %>% filter(status=="HP") %>% filter(is.na(KEGG)) %>% arrange(temp.lc.Metabolite,Level,desc(log2FC))%>% distinct(temp.lc.Metabolite, .keep_all = TRUE)
    #For cases when 1+ of the duplicates does not have a confidence level "1"
    temp.notHP.names <- combined_results %>% filter(is.na(status)) %>% arrange(temp.lc.Metabolite,plyr::desc(mz.match.score),desc(log2FC)) %>% distinct(temp.lc.Metabolite, .keep_all = TRUE)
    } 
if (test_type %in% c("anova","lm","lme")){
    #For cases when 1+ of the duplicate KEGG IDs has a confidence level "1", pick the one with conf level 1 (then sort by p-value)
    temp.HP.Keggs <- combined_results %>% filter(status=="HP") %>% arrange(KEGG,Level,adj.p.value) %>% distinct(KEGG, .keep_all = TRUE)
    #For cases when 1+ of the duplicates does not have a confidence level "1", pick the best match and then the best p-value
    temp.notHP.Keggs <- combined_results %>% filter(is.na(status)) %>% filter(!is.na(KEGG)) %>% arrange(KEGG,plyr::desc(mz.match.score),adj.p.value) %>% distinct(KEGG, .keep_all = TRUE)
  
    #For peaks with Compound.names but no KEGG IDs
    #For cases when 1+ of the duplicates has a confidence level "1"
    temp.HP.names <- combined_results %>% filter(status=="HP") %>% filter(is.na(KEGG)) %>% arrange(temp.lc.Metabolite,Level,adj.p.value)%>% distinct(temp.lc.Metabolite, .keep_all = TRUE)
    #For cases when 1+ of the duplicates does not have a confidence level "1"
    temp.notHP.names <- combined_results %>% filter(is.na(status)) %>% arrange(temp.lc.Metabolite,plyr::desc(mz.match.score),adj.p.value) %>% distinct(temp.lc.Metabolite, .keep_all = TRUE)
    }
#Put them back together
combined_results <- rbind(temp.HP.Keggs,temp.notHP.Keggs,temp.HP.names,temp.notHP.names) %>% dplyr::select(-status)%>% dplyr::select(-temp.lc.Metabolite) %>% distinct()
combined_results <- combined_results %>% mutate(Confidence=case_when(Level==3~"LOW CONFIDENCE ID",.default=""))

combined_results$row.mode <- paste(combined_results$id,combined_results$mode,sep="_")

#Make a key list of values row.id and mode that are saved in the final combined results to extract from the model fits if needed
orig.combined_results$row.mode <- paste(orig.combined_results$id,orig.combined_results$mode,sep="_")

retained_peaks <- combined_results$row.mode
#Does this really need to be done or should the full results of the lm/anova be reported? Not sure
orig.combined_results <- orig.combined_results %>% filter(row.mode %in% retained_peaks)
combined_results <- orig.combined_results
#Merge samples' raw intensities with main contrast results 
combined_proc <- rbind(neg.output[[3]],pos.output[[3]]) %>% dplyr::select(c(id,mode,metadata$Sample.Name))
combined_proc$row.mode <- paste(combined_proc$id,combined_proc$mode,sep="_")
combined_results <- merge(combined_results,combined_proc,by="row.mode")

#Merge samples' normalized intensities with main contrast results
combined_norm <- rbind(neg.output[[4]],pos.output[[4]]) %>% dplyr::select(c(id,mode,metadata$Sample.Name))

combined_norm$row.mode <- paste(combined_norm$id,combined_norm$mode,sep="_")

colnames(combined_norm) <- gsub(paste0("(", paste(metadata$Sample.Name, collapse="|"), ")"), 
                      paste0("norm.", "\\1"), 
                      colnames(combined_norm))

combined_results <- merge(combined_results,combined_norm,by="row.mode")

Client_Data_Download <- list()
Client_Data_Download <- c(list(combined_results),pos.output,neg.output)
names(Client_Data_Download)[1] <- "report_results"
#For t-test studies, remove the empty placeholder elements
Client_Data_Download <- Client_Data_Download[!names(Client_Data_Download) %in% c("PosEmpty","NegEmpty")]
#Remove unneeded and redundant columns
for (i in 1:length(Client_Data_Download)){
  tryCatch({
  if("mode.x" %in% colnames(Client_Data_Download[[i]])){
    Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-c("mode.x","mode.y"))
  }
    if("Row.names.y" %in% colnames(Client_Data_Download[[i]])){
      Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-c("Row.names.x","Row.names.y"))
    }
    if("status" %in% colnames(Client_Data_Download[[i]])){
      Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-c("status"))
    }
    if(c("id.y") %in% colnames(Client_Data_Download[[i]])){
      Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-c("id.x","id.y"))
    }
  if(c("temp.lc.Metabolite") %in% colnames(Client_Data_Download[[i]])){
    Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-c("temp.lc.Metabolite"))
  }
    if ("row.ID" %in% colnames(Client_Data_Download[[i]]) & "id" %in% colnames(Client_Data_Download[[i]])){
      Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-c("id"))
    } else if ("id" %in% colnames(Client_Data_Download[[i]])) {
        Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::rename("row.ID"="id")}


Client_Data_Download[["metadata"]]["Sample.name"] <-  rownames(Client_Data_Download[["metadata"]]) 
Client_Data_Download[["metadata"]] <- Client_Data_Download[["metadata"]] %>% relocate("Sample.name")

empty_columns <- colSums(is.na(Client_Data_Download[[i]]) | Client_Data_Download[[i]] == "") == nrow(Client_Data_Download[[i]])
Client_Data_Download[[i]] <- Client_Data_Download[[i]][,!empty_columns]
}, error = function(e) {
  message(paste0("Plot data not processed, OK", i, ": ", conditionMessage(e)))
})
}
###########################################################################
###############Adjust confidence levels for pre-stats compounds############
###########################################################################
#Run this if there is another annotation in the mzmine output
if(exists("Pos_pre_stats_conf_levels")){
  Pos_pre_stats_conf_levels$Mode <- "Pos"
  Neg_pre_stats_conf_levels$Mode <- "Neg"
  Pre_stats_conf <- rbind(Pos_pre_stats_conf_levels,Neg_pre_stats_conf_levels)
  Pre_stats_conf$id = as.character(Pre_stats_conf$id)
  # Loop through each row in df1
  for(n in 1:length(Client_Data_Download)){
    for (i in 1:nrow(Pre_stats_conf)) {
      # Get the row and mode values from Pre_stats_conf
      row_value <- Pre_stats_conf$id[i]
      mode_value <- Pre_stats_conf$Mode[i]
  
      # Find the matching row in Client_Data_Download[[i]]
      matching_row <- Client_Data_Download[[n]]$id == row_value & Client_Data_Download[[n]]$mode == mode_value
  
      # Check if there is a match in Client_Data_Download[[i]] and Pre_stats_conf's confidence value is not an empty string
      if (any(matching_row) && Pre_stats_conf$Confidence[i] != "") {
        # Update the 'Level' value in Client_Data_Download[[i]] with the corresponding 'Confidence' value from Pre_stats_conf
        Client_Data_Download[[n]]$Level[matching_row] <- Pre_stats_conf$Confidence[i]
      }
  }
}
}
  return(Client_Data_Download)
}
#save(file=paste(client,".ReportingInput",Sys.Date(),".RDATA",sep=""),Client_Data_Download)

