###################################
###Begin normalization and stats###
###################################
#Changes v5 to v6 just to deal with slightly different Generating_Reporting_Inputs.R to accommodate different way of naming input dataset rows/cols
#Need to make an updated v6 for the next t-test that doesn't deal with "Metabolite" but with "id"
#v10 - added argument to specify whether emmeans should be pairwise or given a list of custom contrasts
#test_type = "lm","lme","anova","t.test"
#contrast_var = "Class"
#anova_formula= Metabolite ~ Class
#lm_model
#emmeans_var=~Class
#emmeans_contrasts=ALL or like: `list(c("Level1", "Level2"), c("Level3", "Level4"))`
#I just need to find out why there is an "X" before the sample names in the dataset 4/20 4:32
SECIM_Metabolomics <-function(dataset,peakdata,num_meta,original_data,contrast_var,anova_formula,lm_model,
                              test_type,subset,SECIM_column,emmeans_var,mode,metid_DB_file,client,
                              metadata,paired=FALSE,batch_correct=NULL,
                              rowNorm="SumNorm",transNorm="LogNorm",scaleNorm="ParetoNorm",samples_to_drop_post_norm){
  

  # Store the names of objects in the global environment before loading the file
  before <- ls()
  print(ls)
  load(metid_DB_file)
  metid_DB <- setdiff(ls(), before)[2]
  print(metid_DB)
  
  dataset[is.na(dataset)] <- 0
  #makes a metadataset frame out of the header rows of dataset
  md <- data.frame(t(data.frame(dataset[1:num_meta,])))
  colnames(md) <- md[1,]
  md <- md %>% slice(-1)
  rownames(md) <- gsub("^X","",rownames(md))
  
  if (num_meta==1){
  dataset[dataset==0] <- NA
  write.csv(file=paste(client,mode,"metab.in.csv"),dataset,row.names = FALSE)
  mSet<-InitDataObjects("pktable", "stat", FALSE)
  mSet<-Read.TextData(mSet, paste(client,mode,"metab.in.csv"), "colu", "disc")
  
  mSet<-SanityCheckDataHRK(mSet)
  #Get the metadataset from mSet so the order matches throughout
  md <- data.frame(mSet[["dataSet"]][["meta.info"]])
  rownames(md) <- mSet[["dataSet"]][["url.smp.nms"]]
  } else {
  dataset[dataset==0] <- NA
  dataset <- dataset[-c(1:num_meta),]
  write.csv(file=paste(client,mode,"metab.in.csv"),dataset,row.names = FALSE)
  write.csv(file=paste(client,mode,"metab.meta.in.csv"),md,row.names=TRUE)
  mSet<-InitDataObjects("pktable", "mf", FALSE)
  mSet<-SetDesignType(mSet, "multi")
  mSet<-Read.TextDataTs(mSet, paste(client,mode,"metab.in.csv"), "colmf")
  mSet<-ReadMetaData(mSet,paste(client,mode,"metab.meta.in.csv"))
  mSet<-SanityCheckDataHRK(mSet)
  #Get the metadataset from mSet so the order matches throughout
  md <- data.frame(mSet[["dataSet"]][["meta.info"]])
  #rownames(md) <- mSet[["dataSet"]][["url.smp.nms"]]
}
  
  #Use metaboanalystR instead of custom functions
  
  mSet<-RemoveMissingPercent(mSet, percent=0.5)
  mSet<-ImputeMissingVar(mSet, method="knn_var")
  SanityCheckDataHRK(mSet)
  #Blank-feature filtering
  mSet<-FilterVariable(mSet, var.filter="iqr", "F", 10) #changed this, I guess there was a metaboanalystR update when reinstalled but that may lead to other issues
  #Normalization
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, rowNorm=rowNorm, transNorm=transNorm, scaleNorm=scaleNorm, ratio=FALSE, ratioNum=20)
  proc.data <- qread("data_proc.qs") #for plotting
  norm.data <- qread("complete_norm.qs") #for plotting
  
  #Write the pre and post normalization plots
  plots <- Norm_Plots(proc.data=proc.data,norm.data=norm.data)
  
  #Add metadataset for client download
  data.final <- data.frame(t(norm.data))
  colnames(data.final) <- gsub("^X","",colnames(data.final))
  
  data.final$id <- rownames(data.final) #Quadruple check the rownames in data final match the id names from mzmine
  data.final <- data.final %>% relocate(id)
  if(num_meta==1){
  data.final <- rbind(dataset[1:num_meta,],data.final)
  
  data.proc <- data.frame(t(proc.data))
  colnames(data.proc) <- gsub("^X","",colnames(data.proc))
  
  data.proc$id<- rownames(data.proc)
  data.proc <- data.proc %>% relocate(id)
  # Get the column names from dataset[1:num_meta,]
  meta_cols <- colnames(dataset[1:num_meta,])
  # Order the columns of data.proc based on meta_cols
  data.proc <- data.proc[, order(colnames(data.proc) %in% meta_cols)]
  
  data.proc <- rbind(dataset[1:num_meta,],data.proc)
  }else{
    #md <- md[match(colnames(data.final)[2:ncol(dataset)], rownames(md)), ]
    #tempmd <- data.frame(t(md));tempmd$id <- NA;tempmd <- tempmd %>% relocate(id);colnames(tempmd) <- gsub("^X","",colnames(tempmd))
    tmd <- data.frame(t(md))
    tmd$id <- ""
    tmd <- tmd %>% dplyr::relocate(id)
    colnames(tmd) <- gsub("^X","",colnames(tmd))
    data.final <- rbind(tmd,data.final)
    
    data.proc <- data.frame(t(proc.data))
    colnames(data.proc) <- gsub("^X","",colnames(data.proc))
    data.proc$id <- rownames(data.proc)
    data.proc <- data.proc %>% relocate(id)
    #colnames(data.proc) <- gsub("X","",colnames(data.proc))#v6 This was needed but problematic. See if it's needed again
    # Get the column names from dataset[1:num_meta,]
    #meta_cols <- colnames(dataset[1:num_meta,])
    # Order the columns of data.proc based on meta_cols
    #data.proc <- data.proc[, order(colnames(data.proc) %in% meta_cols)]
    #browser()
    #data.proc <- data.proc[, colnames(tmd)]
    data.proc <- rbind(tmd,data.proc)  
  }
  

  ###############
  #####STATS#####
  ###############
  ##Ref:https://ucdavis-bioinformatics-training.github.io/2018-September-Bioinformatics-Prerequisites/friday/linear_models.html
  ###############
  
  #Batch correction with 
  if (batch_correct=="ComBat"){
    # Extracting batch information
    batch_info <- as.numeric(data.final[num_meta, -1])  # Assuming the last row of metadata contains batch information and excluding the 'id' column
    
    # Extracting and preparing the matrix for ComBat
    data_matrix <- as.matrix(data.final[-(1:num_meta), -1])  # Excludes the first three rows (metadata) and the first column ('id')
    rownames(data_matrix) <- data.final[-(1:num_meta), 1]  # Sets row names to the 'id' of features
    data_matrix <- apply(data_matrix, 2, as.numeric)  # Ensure numeric data
    # Assuming data_matrix is correctly formatted as numeric matrix
    original_feature_names <- rownames(data.final[-(1:num_meta), -1])   # Features, e.g., compounds or genes
    original_sample_names <- colnames(data_matrix)   # Sample identifiers
    
    # Now, your data is ready for batch correction. `data_matrix` is your input matrix, and `batch_info` is your batch vector.
    library(sva)
    
    # Applying ComBat for batch correction
    data_corrected <- ComBat(dat = data_matrix, batch = batch_info, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
    rownames(data_corrected) <- original_feature_names
    colnames(data_corrected) <- original_sample_names
    data_final_corrected <- as.data.frame(data_corrected)
    data_final_corrected$id <- rownames(data.final)[c((num_meta+1):nrow(data.final))]
    rownames(data_final_corrected) <- data_final_corrected$id
    data_final_corrected <- data_final_corrected %>% dplyr::relocate(id)  
    #add metadata columns
    data_final_corrected <- rbind(data.final[c(1:num_meta),],data_final_corrected)
    data.final <- data_final_corrected
  }
  if(batch_correct=="limma"){
    library(limma)
    batch_info <- as.numeric(data.final[num_meta, -1])  # Assuming the third row contains batch information and excluding the 'id' column
    
    # Extracting and preparing the matrix for ComBat
    data_matrix <- as.matrix(data.final[-(1:num_meta), -1])  # Excludes the first three rows (metadata) and the first column ('id')
    rownames(data_matrix) <- data.final[-(1:num_meta), 1]
    # Ensure data_matrix is numeric
    data_matrix <- apply(data_matrix, 2, as.numeric)
    # Assuming data_matrix is correctly formatted as numeric matrix
    original_feature_names <- rownames(data.final[-(1:num_meta), -1])  # Features, e.g., compounds or genes
    original_sample_names <- colnames(data_matrix)   # Sample identifiers
    
    # Apply removeBatchEffect
    data_matrix_corrected <- removeBatchEffect(data_matrix, batch=batch_info)
    
    # Restore row names (features/metabolites) if they were removed
    rownames(data_matrix_corrected) <- original_feature_names
    # Convert the corrected matrix back to a data frame, transposing it so that features become columns
    data_corrected_df <- as.data.frame(data_matrix_corrected)
    
    # Add the sample identifiers as the first column
    data_corrected_df$id <- original_feature_names
    
    # Use dplyr::relocate to move the 'id' column to the first position
    data_corrected_df <- dplyr::relocate(data_corrected_df, id)
    
    # Re-add the metadata rows at the top of the dataframe
    # Assuming you've stored your metadata similarly to how you extracted batch_info
    data_final_corrected <- rbind(data.final[c(1:num_meta),],data_corrected_df)
    data.final <- data_final_corrected
    
  }
  #After data.final has been created, create a data.final.subset that is just for the contrast of interest IF samples_to_drop_post_norm is not NULL
  if(!is.null(samples_to_drop_post_norm)){
  samples_to_drop_post_norm <- gsub("-","_",samples_to_drop_post_norm)
  #save the total normalized data for PCA and download
  data.final.total <- data.final
  #Remove the samples to be dropped
  data.final.subset <- data.final %>% select(!samples_to_drop_post_norm)
  #Rename the subset data for compatibility with downstream code
  data.final <- data.final.subset
  }
  
  if (test_type=="t.test"){
    print(contrast_var)
    #Step7: Stats, T-test
    ttest=list()
    ttest <- foreach (i = (num_meta+1):nrow(data.final),.packages=c("dplyr","stats"))%do% {
      temp <- data.frame(t(data.final[c(1:(num_meta),i),]))
      if(num_meta==1){
        colnames(temp) <- c(contrast_var,"Metabolite")
      }else{
      colnames(temp) <- c(colnames(temp[1,][1:num_meta]),"Metabolite")
      }
      temp <- temp[-1,]
      temp$Metabolite <- as.numeric(temp$Metabolite)
      #v11 Adds row ID to temp
      temp$rowID <- data.final$id[i][length(data.final$id[i])]
      exp1 <- expr(Metabolite ~ !!ensym(contrast_var)) #changed from !!ensym to ensym
      if(is.null(subset)){
        if(paired==TRUE){
          # Ensure the data is ordered consistently for both groups
          temp <- temp[order(temp$rowID, temp$Subject), ]
          
          # Extract the Metabolite values for each group
          group1 <- temp$Metabolite[temp$Class == levels(as.factor(temp[[contrast_var]]))[1]]
          group2 <- temp$Metabolite[temp$Class == levels(as.factor(temp[[contrast_var]]))[2]]
          
          # Check if both groups have the same length
          if(length(group1) == length(group2)) {
            # Run a paired t-test
            ttest.res <- tidy(t.test(group1, group2, paired = TRUE))
          } else {
            stop("The groups do not have the same number of observations for a paired t-test.")
          }
          
          # The result is stored in ttest.res
          ttest.res$contrast <- paste(levels(as.factor(temp[[contrast_var]]))[1],"-",levels(as.factor(temp[[contrast_var]]))[2])
        }else {
        ttest.res <- tidy(t.test(formula = eval(exp1), data = temp))
        ttest.res$contrast <- paste(levels(as.factor(temp[[contrast_var]]))[1],"-",levels(as.factor(temp[[contrast_var]]))[2])
        }
        ttest.res
        }else{
          if(paired==TRUE){
          #subset=list(list("GI_MB_TBI","GII_MB_Sham"),list("GI_CX_TBI","GII_CX_Sham"))
          ttest.res=list()
          for(n in 1:length(subset)){
           
            temp <- temp[order(temp$rowID, temp$Subject), ]
            temp <- temp %>% filter(Class %in% subset[[n]])
            
            # Extract the Metabolite values for each group
            group1 <- temp$Metabolite[temp$Class == levels(as.factor(temp[[contrast_var]]))[1]]
            group2 <- temp$Metabolite[temp$Class == levels(as.factor(temp[[contrast_var]]))[2]]
            
            # Check if both groups have the same length
            if(length(group1) == length(group2)) {
              # Run a paired t-test
              ttest.res <- tidy(t.test(group1, group2, paired = TRUE))
            } else {
              stop("The groups do not have the same number of observations for a paired t-test.")
            }
            
            
          ttest.res[[n]]$contrast <- paste(subset[[n]][1],"-",subset[[n]][2])
          }
            ttest.res <- do.call("rbind",ttest.res)
            ttest.res$id <- temp$rowID[[1]]
          } else{
            ttest.res=list()
            for(n in 1:length(subset)){
              ttest.res[[n]] <- tidy(t.test(formula = eval(exp1), data = temp %>% filter(Class %in% subset[[n]])))
              ttest.res[[n]]$contrast <- paste(subset[[n]][1],"-",subset[[n]][2])
            }
            ttest.res <- do.call("rbind",ttest.res)
            ttest.res$id <- temp$rowID[[1]]
          }
            ttest.res 
          }
        }
    ttest.results <- do.call("rbind", ttest)
    if(is.null(subset)){
    ttest.results$id <- data.final$id[(num_meta+1):nrow(data.final)]
    }
    ttest.results <- relocate(ttest.results,id)
    } 
  if (test_type %in% c("anova","lm","lme")){
    #Step 7: Stats, anova or lm 7/12, I think this section needs "Metabolite" changed to "id" and it should also be changed in the function argument for model
    fit=list()
    fit.results=list()
    emmeans=list()

    if(test_type=="anova"){ #I re-wrote this on 5/11/23 when there were >1 metadataset, but this might now not be right for 1 metadata. Need to go back and decide how to deal with that - two functions probably.
        fit_emmeans <- foreach (i = (num_meta+1):nrow(data.final), .packages = c("dplyr","emmeans","stats")) %do% {
          tempdf <- data.frame(t(data.final[c(1:(num_meta), i), ]))
          if(num_meta==1){
          colnames(tempdf) <- c((tempdf)[,1][1:num_meta], "id") #12.2 vs 12.1 changed "Metabolite" to "id"
          }else{
            colnames(tempdf) <- c(colnames(tempdf)[1:num_meta], "id")#12.2 vs 12.1 changed "Metabolite" to "id"
          }
          tempdf <- tempdf[-1, ]
          tempdf$id <- as.numeric(tempdf$id)
          tempdf$Class <- as.character(tempdf$Class)
          tempdf$Class <- as.factor(tempdf$Class)
          if("ID" %in% colnames(tempdf)){ #this is for repeated measures
            tempdf$ID <- as.factor(tempdf$ID)}
          #fit <- aov(as.formula(anova_formula), data = tempdf)
          fit <- do.call(aov, args = list(anova_formula, tempdf))
          #if("ALL" %in% emmeans_contrasts){
          emmeans_obj <- tidy(pairs(emmeans(fit,emmeans_var,data=tempdf)))
          list(fit = fit, emmeans = emmeans_obj)
        }
        # Extract fit and emmeans from the list
        fit <- lapply(fit_emmeans, function(x) x$fit)
        emmeans <- lapply(fit_emmeans, function(x) x$emmeans)
        
      }
    if (test_type=="lm"){
        fit_emmeans <- foreach (i = (num_meta+1):nrow(data.final),.packages=c("dplyr","stats","emmeans"))%do% {
        temp <- data.frame(t(data.final[c(1:(num_meta),i),]))
        colnames(temp) <- c(temp[1,][1:num_meta],"Metabolite"); temp <- temp[-1,]
        temp$Metabolite <- as.numeric(temp$Metabolite)
        fit <- lm(lm_model, data = temp)
        
        # emmeans analysis
        emmeans_obj <- emmeans(fit, specs = emmeans_var)
        pairwise_emmeans_obj <- tidy(pairs(emmeans_obj))
        
        # Return both fit and emmeans results
        list(fit = fit, emmeans = pairwise_emmeans_obj)
        }
        # Extract model fits and emmeans summary tables from the list
        fit <- lapply(fit_emmeans, function(x) x$fit)
        emmeans <- lapply(fit_emmeans, function(x) x$emmeans)
      }
    if (test_type=="lme"){
      fit_emmeans <- foreach(i = (num_meta+1):nrow(data.final), .packages = c("dplyr", "stats", "lmerTest", "emmeans", "broom")) %do% {
        temp <- data.frame(t(data.final[c(1:(num_meta), i), ]))
        colnames(temp) <- c(colnames(temp)[1:num_meta], "id")
        temp <- temp[-1, ]
        temp$id <- as.numeric(temp$id)
        #Convert variables to factors
        temp$Class <- as.factor(temp$Class)
        temp$Subject <- as.factor(temp$Subject)
        temp$Batch <- as.factor(temp$Batch)
        
        # Fit the model
        fit <- lmerTest::lmer(lm_model, data = temp)
        
        # emmeans analysis
        emmeans_obj <- emmeans(fit, specs = emmeans_var)
        pairwise_emmeans_obj <- tidy(pairs(emmeans_obj))
        
        # Return both fit and emmeans results
        list(fit = fit, emmeans = pairwise_emmeans_obj)
      }
      
      # Extract model fits and emmeans summary tables from the list
      fit <- lapply(fit_emmeans, function(x) x$fit)
      emmeans <- lapply(fit_emmeans, function(x) x$emmeans)
      
    } 
    if(test_type %in% c("anova","lm")){
      fit.results <- lapply(fit,tidy)
      }
    if (test_type == "lme") {
      library(broom.mixed)
      # Using broom.mixed to tidy up the lmerTest model object
      fit.results <- lapply(fit, function(x) tidy(x))
    }
    
    names(fit.results) <- data.final$id[(num_meta+1):nrow(data.final)]
    fit.results <- do.call("rbind", fit.results)
    fit.results$id <- rownames(fit.results)
    fit.results$id  <- gsub("\\.[0-9\\+]","",fit.results$id) #added the rowID key
  
    names(emmeans) <- data.final$id[(num_meta+1):nrow(data.final)]
    emmeans.results <- do.call("rbind", emmeans)
    emmeans.results$id <- rownames(emmeans.results)
    
    emmeans.results$id <- gsub("\\.[0-9\\+]","",emmeans.results$id) #added the rowID key
  }

  ######################
  #####GROUP-MEANS######
  ######################
  #Step 1: Define the contrasts, groups, and samples
  group_samples <- list()
  if (test_type=="t.test"){
    contrast_vec <- gsub(" ","",levels(as.factor(ttest.results$contrast)))  
    contrast_vec <- sapply(contrast_vec, function(x) gsub("[()]", "", x))
  } 
  if (test_type %in% c("anova","lm","lme")) {
    contrast_vec <- gsub(" ","",levels(as.factor(emmeans.results$contrast)))
    contrast_vec <- sapply(contrast_vec, function(x) gsub("[()]", "", x)) #emmeans will introduce "(" into the contrast names to deal with special chars in variables
  }
  if (test_type == "nostats"){
    contrast_vec <- combn(levels(as.factor(metadata$Class)), 2, FUN = function(x) paste0(x[1], "-", x[2]), simplify = FALSE)
  }
  
  if (test_type == "nostats"){
    group_vec <- unique(unlist(sapply(as.list(contrast_vec),function(x) str_split(x,"-"))))} else {
  group_vec <- unique(unlist(sapply(as.list(contrast_vec),function(x) str_split(x,"-"))))
    }

  
  for (i in 1:length(group_vec)){
    group_samples[[i]] <- rownames(md %>% filter(!!as.symbol(contrast_var) == group_vec[[i]]))
  }
  
  #Step 2: Calculate the group means from the processed peak intensity data
  group_means <- list()
  for (i in 1:length(group_vec)){
    group_means[[i]] <- data.frame(t(proc.data)) %>% dplyr::rowwise() %>% dplyr::mutate(!!group_vec[[i]] := mean(c_across(matches(group_samples[[i]]))))
  }
  #Means is formatted with one row per peak and one column per group
  means <- data.frame(matrix(ncol=0,nrow=nrow(data.frame(t(proc.data)))))
  for (i in 1:length(group_vec)){
    means <- cbind(means,group_means[[i]][,ncol(group_means[[i]])])
  }
  rownames(means) <- colnames(proc.data)
  means$`id` <- as.numeric(rownames(means))
  
  #Step 3: Calculate the per-contrast fold-changes
  #######################
  #CONTRAST FOLD-CHANGES#
  #######################
  contrast_fold_changes <- list()
  for (i in 1:length(contrast_vec)){
    contrast_fold_changes[[i]] <- means %>% dplyr::rowwise() %>% 
      dplyr::mutate(!!contrast_vec[[i]]:=!!as.symbol(str_split(contrast_vec[[i]],"-")[[1]][1])/!!as.symbol(str_split(contrast_vec[[i]],"-")[[1]][2]))%>%
      dplyr::select(c(!!contrast_vec[[i]],id))
    log2FC <- log2(contrast_fold_changes[[i]][contrast_vec[[i]]])
    contrast_fold_changes[[i]] <- cbind(contrast_fold_changes[[i]],log2FC)
    colnames(contrast_fold_changes[[i]]) <- c("FC","id","log2FC")
    contrast_fold_changes[[i]]["contrast"] <- contrast_vec[[i]]
    contrast_fold_changes[[i]] <- contrast_fold_changes[[i]][,c(2,4,1,3)]
  }
  names(contrast_fold_changes) <- contrast_vec
  all_fold_changes <- do.call("rbind",contrast_fold_changes)
  means_FC <- merge(means,all_fold_changes,by="id") #The means and FC dataset is now keyed by variable "rowID" NOT rowname
  means_FC$id <- as.character(means_FC$id)
  means_FC$contrast <- gsub("-"," - ",means_FC$contrast)
  
  #rownames(means_FC) <- colnames(proc.data)
  
  ###########################
  ########Peak Annotation####
  ###########################
  
  #Step 11: metID of the metabolites in the peaktable
  #Make a new object that has the old mzmine style colnames because that is what the function expects
  peakdataformetid <- peakdata
  peakdataformetid <- peakdataformetid[, c(1, 3, 2,4:ncol(peakdataformetid))]
  peakdataformetid[-c(1, 4)] <- lapply(peakdataformetid[-c(1, 4)], function(x) as.numeric(as.character(x)))
  
  colnames(peakdataformetid) <- c("row ID","row m/z","row retention time",colnames(peakdata)[4:ncol(peakdata)])
  metid = convet_mzmine2mass_dataset(x = peakdataformetid %>% dplyr::select(!compound) ,rt_unit = "minute")
  if(mode=="Pos"){
  metid <-
    annotate_metabolites_mass_dataset(object = metid, 
                                      ms1.match.ppm = 5, 
                                      rt.match.tol = 10001, 
                                      polarity = "positive",
                                      database = get(metid_DB),
                                      column = "rp_custom",
                                      #column = "rp",
                                      threads=4,
                                      candidate.num=1)
  } else if(mode=="Neg"){
    metid <-
      annotate_metabolites_mass_dataset(object = metid, 
                                        ms1.match.ppm = 10, 
                                        rt.match.tol = 10001, 
                                        polarity = "negative",
                                        database = get(metid_DB),
                                        column = "rp_custom",
                                        threads=4,
                                        candidate.num=1)
  }
  metid.result <- merge(
    metid@variable_info,metid@annotation_table,by="variable_id",all.x=TRUE)
  
  #Convert rt back to minute and fix to two decimal places
  metid.result<- mutate(metid.result, 
                        rt = rt/60)
  metid.result <- mutate(metid.result, 
                         rt = format(round(rt,digits=2),nsmall=2))
  #Fix the mass to four decimal places
  metid.result <- mutate(metid.result, 
                         mz = format(round(mz,digits=4),nsmall=4))
  #Get KEGG hierarchy for peaks named using metid
  metid.result <- metid.result %>% dplyr::rename(KEGG=KEGG.ID)
  metid.result <- as.data.frame(assign_hierarchy(count_data = metid.result, keep_unknowns = TRUE, identifier = "KEGG"))
  ##metid.result$variable_id == peakdata$`row ID` so can be used for downstream combining
  
  #Step 12: KEGG ID of metabolites in the peaktable
  # Step 12: KEGG ID of metabolites in the peaktable
  
  # Check both mode and SECIM_column values
  if (mode == "Pos") {
    if (SECIM_column == "Evospere-PFP") {
      KEGG.compound <- read.csv("/home/hkates/blue_garrett/SECIM_Reporting/InputFiles/library_positive1_MK_QE2.csv") %>%
        distinct() %>%
        mutate(KEGG = na_if(KEGG, "null"))%>%
                 mutate(KEGG = na_if(KEGG, "None"))%>%
                          mutate(KEGG = na_if(KEGG, ""))
    } else if (SECIM_column == "ACE-PFP") {
      KEGG.compound <- read.csv("/home/hkates/blue_garrett/SECIM_Reporting/InputFiles/Positive_ACE-PFPs.csv") %>%
        distinct() %>%
        mutate(KEGG = na_if(KEGG, "null"))%>%
        mutate(KEGG = na_if(KEGG, "None"))%>%
        mutate(KEGG = na_if(KEGG, ""))
    }
  } else if (mode == "Neg") {
    if (SECIM_column == "Evospere-PFP") {
      KEGG.compound <- read.csv("/home/hkates/blue_garrett/SECIM_Reporting/InputFiles/library_negative1_MK02.csv") %>%
        mutate(KEGG = na_if(KEGG, "null"))%>%
        mutate(KEGG = na_if(KEGG, "None"))%>%
        mutate(KEGG = na_if(KEGG, ""))
    } else if (SECIM_column == "ACE-PFP") {
      KEGG.compound <- read.csv("/home/hkates/blue_garrett/SECIM_Reporting/InputFiles/Negative_ACE-PFP.csv") %>%
        mutate(KEGG = na_if(KEGG, "null"))%>%
        mutate(KEGG = na_if(KEGG, "None"))%>%
        mutate(KEGG = na_if(KEGG, ""))
    }
  }
  
  #KEGG.compound <- KEGG.compound %>% dplyr::select(c("name","KEGG"))
  KEGG.compound <- KEGG.compound %>% mutate_at(c('KEGG'), ~na_if(., "")) 
  #If a compound name has >1 KEGG IDs, save the first one (I don't know a better way to pick atm)
  KEGG.compound <- KEGG.compound %>% group_by(name) %>% slice(1)
  peakdata.KEGG <- merge(peakdata,KEGG.compound,by.x="compound",by.y="name",all.x=TRUE)
  #Deal with cases where the same name has two KEGG IDs
  peakdata.KEGG <- as.data.frame(assign_hierarchy(count_data = peakdata.KEGG, keep_unknowns = TRUE, identifier = "KEGG")) 
  peakdata.KEGG <- peakdata.KEGG[,setdiff(colnames(peakdata.KEGG), metadata$Sample.Name)]
  
  #Reduce the metid.result to only variable_ids that is.na("row.identity..main.ID.") and peakdata.KEGG to id that are NOT NA
  rowid.SECIM.na <- dplyr::filter(data.frame(peakdata.KEGG),is.na(peakdata.KEGG$compound)) %>% dplyr::select(id) %>% unlist()#this is the row IDs of peaks that are NA 
  metid.result <- data.frame(metid.result) %>% filter(variable_id %in% rowid.SECIM.na) #Get only the metid peaks for which SECIM ID was NA
  peakdata.KEGG <- data.frame(peakdata.KEGG) %>% filter(!id %in% rowid.SECIM.na)
  #Prepare the two datasetframes for rbind
  peakdata.KEGG$Level=1 #If the KEGG ID came from SECIM, confidence is a "1"
  metid.result <- metid.result %>% dplyr::rename("id"="variable_id")
  metid.result$id <- as.double(metid.result$id)
  metid.result$mz <- as.double(metid.result$mz)
  metid.result$rt <- as.double(metid.result$rt)
  #
  peakdata.KEGG$id <- as.double(peakdata.KEGG$id)
  peakdata.KEGG$mz <- as.double(peakdata.KEGG$mz)
  peakdata.KEGG$rt <- as.double(peakdata.KEGG$rt)
  metid.result <- metid.result %>% dplyr::rename("compound"="Compound.name")
  peak_annotation <- dplyr::bind_rows(peakdata.KEGG, metid.result)
  #If there is no compound from SECIM or metID, use the mz_rt
  peak_annotation <- peak_annotation %>% mutate(compound = case_when(is.na(compound) ~ paste(mz,rt,sep="_"),.default = compound)) %>%
  mutate(KEGG = case_when(KEGG=="" ~ NA,.default = KEGG)) 
  peak_annotation$mode <- mode
                            
  outputs_list <- list()
  # Merge results with metabolite annotation and name the dataset to be saved with mode
  if (test_type == "t.test") {
    outputs_list[[1]] <- merge(ttest.results, peak_annotation, by = "id", all.x = TRUE)
    outputs_list[[1]]$contrast <- gsub("[()]", "", outputs_list[[1]]$contrast)
    outputs_list[[1]] <- outputs_list[[1]] %>% inner_join(means_FC, by = c('id', 'contrast'))
    outputs_list[[1]] <- outputs_list[[1]] %>% relocate(contrast, compound)
    outputs_list[[1]]$adj.p.value <- p.adjust(outputs_list[[1]]$p.value, method = "fdr")
    outputs_list[[1]] <- outputs_list[[1]] %>% relocate(adj.p.value, .after = p.value)
    outputs_list[[2]] <- print("Empty for t test")
  } 
  if (test_type %in% c("lm", "anova", "lme")) {
    outputs_list[[1]] <- merge(emmeans.results, peak_annotation, by = "id", all.x = TRUE)
    outputs_list[[1]]$contrast <- gsub("[()]", "", outputs_list[[1]]$contrast)
    
    if (length(unique(emmeans.results$contrast)) > 1) {
      outputs_list[[1]] <- outputs_list[[1]] %>% inner_join(means_FC, by = c('id', 'contrast'))
      outputs_list[[1]] <- outputs_list[[1]] %>% relocate(contrast, compound)
    } else {
      outputs_list[[1]]$adj.p.value <- p.adjust(outputs_list[[1]]$p.value, method = "fdr")
      outputs_list[[1]] <- outputs_list[[1]] %>% inner_join(means_FC, by = c('id', 'contrast'))
      outputs_list[[1]] <- outputs_list[[1]] %>% relocate(contrast, compound)
    }
    
    outputs_list[[2]] <- merge(fit.results, peak_annotation, by = "id", all.x = TRUE)
    outputs_list[[2]] <- outputs_list[[2]] %>% relocate(compound)
  }
  
  if (test_type == "nostats") {
    outputs_list[[1]] <- merge(means_FC, peak_annotation, by = "id", all.x = TRUE)
    outputs_list[[1]]$contrast <- gsub("[()]", "", outputs_list[[1]]$contrast)
    outputs_list[[1]] <- outputs_list[[1]] %>% relocate(contrast, compound)
    outputs_list[[2]] <- print("Empty for nostats")
  }
  
  outputs_list[[3]] <- merge(data.proc, peak_annotation, by = "id", all.x = TRUE)
  outputs_list[[3]] <- outputs_list[[3]][c(which(outputs_list[[3]]$id == "Class"), setdiff(seq_len(nrow(outputs_list[[3]])), which(outputs_list[[3]]$id == "Class"))),]
  outputs_list[[3]] <- outputs_list[[3]][c(which(is.na(outputs_list[[3]]$id)), setdiff(seq_len(nrow(outputs_list[[3]])), which(is.na(outputs_list[[3]]$id)))),]
  outputs_list[[3]] <- outputs_list[[3]] %>% relocate(compound)
  
  outputs_list[[4]] <- merge(data.final, peak_annotation, by = "id", all.x = TRUE)
  outputs_list[[4]] <- outputs_list[[4]][c(which(outputs_list[[4]]$id == "Class"), setdiff(seq_len(nrow(outputs_list[[4]])), which(outputs_list[[4]]$id == "Class"))),]
  outputs_list[[4]] <- outputs_list[[4]][c(which(is.na(outputs_list[[4]]$id)), setdiff(seq_len(nrow(outputs_list[[4]])), which(is.na(outputs_list[[4]]$id)))),]
  outputs_list[[4]] <- outputs_list[[4]] %>% relocate(compound)
  if (!is.null(samples_to_drop_post_norm)) {
    outputs_list[[if (mode == "Neg") 8 else 7]] <- merge(data.final.total, peak_annotation, by = "id", all.x = TRUE)
    outputs_list[[if (mode == "Neg") 8 else 7]] <- outputs_list[[if (mode == "Neg") 8 else 7]][c(which(outputs_list[[if (mode == "Neg") 8 else 7]]$id == "Class"), setdiff(seq_len(nrow(outputs_list[[if (mode == "Neg") 8 else 7]])), which(outputs_list[[if (mode == "Neg") 8 else 7]]$id == "Class"))),]
    outputs_list[[if (mode == "Neg") 8 else 7]] <- outputs_list[[if (mode == "Neg") 8 else 7]][c(which(is.na(outputs_list[[if (mode == "Neg") 8 else 7]]$id)), setdiff(seq_len(nrow(outputs_list[[if (mode == "Neg") 8 else 7]])), which(is.na(outputs_list[[if (mode == "Neg") 8 else 7]]$id)))),]
    outputs_list[[if (mode == "Neg") 8 else 7]] <- outputs_list[[if (mode == "Neg") 8 else 7]] %>% relocate(compound)
  }
  
  outputs_list[[5]] <- arrangeGrob(plots[["plot_env"]][["p1"]], plots[["plot_env"]][["p3"]], plots[["plot_env"]][["p2"]], plots[["plot_env"]][["p4"]], ncol = 2, top = textGrob("Feature View"))
  outputs_list[[6]] <- arrangeGrob(plots[["plot_env"]][["p5"]], plots[["plot_env"]][["p7"]], plots[["plot_env"]][["p6"]], plots[["plot_env"]][["p8"]], ncol = 2, top = textGrob("Sample View"))
  
  if (mode == "Neg") {
      if (!is.null(samples_to_drop_post_norm)) {
        outputs_list[[7]] <- md %>% dplyr::filter(!rownames(.) %in% samples_to_drop_post_norm)
      } else {
        outputs_list[[7]] <- md
      }
    
    if (test_type == "t.test") {
      names(outputs_list) <- c(
        paste0(mode, ".ttest.metab"),
        paste0(mode, "Empty"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView"),
        "metadata"
      )

      if (!is.null(samples_to_drop_post_norm)) {
        names(outputs_list) <- c(na.omit(names(outputs_list)), paste0(mode, ".total.normalized.data"))
      }
    } 
    if (test_type %in% c("anova", "lm", "lme")) {
      names(outputs_list) <- c(
        paste0(mode, ".emmeans.results.metab"),
        paste0(mode, ".fit.results.metab"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView"),
        "metadata"
      )
      if (!is.null(samples_to_drop_post_norm)) {
        names(outputs_list) <- c(na.omit(names(outputs_list)), paste0(mode, ".total.normalized.data"))
      }
    }
    if (test_type == "nostats") {
      names(outputs_list) <- c(
        paste0(mode, ".FCanlaysis.metab"),
        paste0(mode, "Empty"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView"),
        "metadata"
      )
      if (!is.null(samples_to_drop_post_norm)) {
        names(outputs_list) <- c(na.omit(names(outputs_list)), paste0(mode, ".total.normalized.data"))
      }
    }
  } else {
    if (test_type == "t.test") {
      names(outputs_list) <- c(
        paste0(mode, ".ttest.metab"),
        paste0(mode, "Empty"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView")
      )
      if (!is.null(samples_to_drop_post_norm)) {
        names(outputs_list) <- c(na.omit(names(outputs_list)), paste0(mode, ".total.normalized.data"))
      }
    }
    if (test_type %in% c("anova", "lm", "lme")) {
      names(outputs_list) <- c(
        paste0(mode, ".emmeans.results.metab"),
        paste0(mode, ".fit.results.metab"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView")
      )
      if (!is.null(samples_to_drop_post_norm)) {
        names(outputs_list) <- c(na.omit(names(outputs_list)), paste0(mode, ".total.normalized.data"))
      }
    }
    if (test_type == "nostats") {
      names(outputs_list) <- c(
        paste0(mode, ".FCanlaysis.metab"),
        paste0(mode, "Empty"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView")
      )
      if (!is.null(samples_to_drop_post_norm)) {
        names(outputs_list) <- c(na.omit(names(outputs_list)), paste0(mode, ".total.normalized.data"))
      }
    }
  }
  
  return(outputs_list)
}

# List of files to remove
files_to_remove <- c("/blue/timgarrett/hkates/complete_norm.qs", "/blue/timgarrett/hkates/data_orig.qs",
                     "/blue/timgarrett/hkates/data_prefilter_iqr.csv", 
                     "/blue/timgarrett/hkates/data_proc.qs", "/blue/timgarrett/hkates/prenorm.qs",
                     "/blue/timgarrett/hkates/preproc.qs", "/blue/timgarrett/hkates/row_norm.qs")

# Remove specific files
#file.remove(files_to_remove)

# For patterns, use list.files with a pattern and then remove those files
neg_metab_files <- list.files(pattern = "/blue/timgarrett/hkates/.* Neg metab\\.in\\.csv$")
pos_metab_files <- list.files(pattern = "/blue/timgarrett/hkates/.* Pos metab\\.in\\.csv$")

# Remove files matching patterns
#file.remove(c(neg_metab_files, pos_metab_files))
