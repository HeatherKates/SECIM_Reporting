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
                              test_type,subset,emmeans_var,mode,metid_DB_file,client,metadata){
  

  # Store the names of objects in the global environment before loading the file
  before <- ls()
  print(ls)
  load(metid_DB_file)
  metid_DB <- setdiff(ls(), before)[2]
  print(metid_DB)
  
  dataset[is.na(dataset)] <- 0
  #makes a metadataset frame out of the header rows of dataset
  md <- data.frame(t(data.frame(dataset[1:num_meta,]))); colnames(md) <- md[1,]; md <- md %>% slice(-1)
  
  if (num_meta==1){
  dataset[dataset==0] <- NA
  write.csv(file=paste(client,mode,"metab.in.csv"),dataset,row.names = FALSE)
  mSet<-InitDataObjects("pktable", "stat", FALSE)
  mSet<-Read.TextData(mSet, paste(client,mode,"metab.in.csv"), "colu", "disc")
  mSet<-SanityCheckData(mSet)
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
  mSet<-SanityCheckData(mSet)
  #Get the metadataset from mSet so the order matches throughout
  md <- data.frame(mSet[["dataSet"]][["meta.info"]])
  rownames(md) <- mSet[["dataSet"]][["url.smp.nms"]]
}
  
  #Use metaboanalystR instead of custom functions
  
  mSet<-RemoveMissingPercent(mSet, percent=0.5)
  mSet<-ImputeMissingVar(mSet, method="knn_var")
  SanityCheckData(mSet)
  #Blank-feature filtering
  mSet<-FilterVariable(mSet, filter="iqr", "F", 10)
  #Normalization
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "SumNorm", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
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
    md <- md[match(colnames(data.final)[2:ncol(dataset)], rownames(md)), ]
    tempmd <- data.frame(t(md));tempmd$id <- NA;tempmd <- tempmd %>% relocate(id)
    data.final <- rbind(tempmd,data.final)
    
    data.proc <- data.frame(t(proc.data))
    colnames(data.proc) <- gsub("^X","",colnames(data.proc))
    data.proc$id <- rownames(data.proc)
    data.proc <- data.proc %>% relocate(id)
    colnames(data.proc) <- gsub("X","",colnames(data.proc))#v6
    # Get the column names from dataset[1:num_meta,]
    meta_cols <- colnames(dataset[1:num_meta,])
    # Order the columns of data.proc based on meta_cols
    data.proc <- data.proc[, order(colnames(data.proc) %in% meta_cols)]
    data.proc <- rbind(tempmd,data.proc)  
  }
  

  ###############
  #####STATS#####
  ###############
  ##Ref:https://ucdavis-bioinformatics-training.github.io/2018-September-Bioinformatics-Prerequisites/friday/linear_models.html
  ###############
  
  if (test_type=="t.test"){
    #Step7: Stats, T-test
    ttest=list()
    ttest <- foreach (i = (num_meta+1):nrow(data.final),.packages=c("dplyr","stats"))%do% {
      temp <- data.frame(t(data.final[c(1:(num_meta),i),]))
      colnames(temp) <- c(temp[1,][1:num_meta],"Metabolite"); temp <- temp[-1,]
      temp$Metabolite <- as.numeric(temp$Metabolite)
      #v11 Adds row ID to temp
      temp$rowID <- data.final$id[i][length(data.final$id[i])]
      exp1 <- expr(Metabolite ~ !!ensym(contrast_var)) #changed from !!ensym to ensym
      if(is.null(subset)){
        ttest.res <- tidy(t.test(formula = eval(exp1), data = temp))
        ttest.res$contrast <- paste(levels(as.factor(temp[[contrast_var]]))[1],"-",levels(as.factor(temp[[contrast_var]]))[2])
        ttest.res
        }else{
          #subset=list(list("GI_MB_TBI","GII_MB_Sham"),list("GI_CX_TBI","GII_CX_Sham"))
          ttest.res=list()
          for(n in 1:length(subset)){
          ttest.res[[n]] <- tidy(t.test(formula = eval(exp1), data = temp %>% filter(Class %in% subset[[n]])))
          ttest.res[[n]]$contrast <- paste(subset[[n]][1],"-",subset[[n]][2])
          }
            ttest.res <- do.call("rbind",ttest.res)
            ttest.res$id <- temp$rowID[[1]]
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
        fit <- foreach (i = (num_meta+1):nrow(data.final),.packages=c("dplyr","stats"))%do% {
        temp <- data.frame(t(data.final[c(1:(num_meta),i),]))
        colnames(temp) <- c(temp[1,][1:num_meta],"Metabolite"); temp <- temp[-1,]
        temp$Metabolite <- as.numeric(temp$Metabolite)
        lm(lm_model, data = temp)
        }
      }
    if (test_type=="lme"){
      fit <- foreach (i = (num_meta+1):nrow(data.final),.packages=c("dplyr","stats"))%do% {
        temp <- data.frame(t(data.final[c(1:(num_meta),i),]))
        colnames(temp) <- c(temp[1,][1:num_meta],"id"); temp <- temp[-1,]
        temp$id <- as.numeric(temp$id)
        fit <- lmerTest::lmer(lm_model,data=temp) #Where was I: for mixed-models, emmeans needs to run on this, so figure out a way to write output suitable for emmeans and for saving the dataset frame results
        fit.results <- data.frame(summary(fit)[[10]])
        fit.results$Coefficient <- rownames(fit.results)
        rownames(fit.results) <- NULL
        return(list(fit,fit.results))
      }
    } 
    if(test_type %in% c("anova","lm")){
      fit.results <- lapply(fit,tidy)
      }
    if (test_type=="lme"){
      fit.results <- lapply(fit,'[[', 2)
      fit <-  lapply(fit,'[[', 1)
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
  
  group_vec <- unique(unlist(sapply(as.list(contrast_vec),function(x) str_split(x,"-"))))

  
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
  colnames(peakdataformetid) <- c("row ID","row m/z","row retention time",colnames(peakdata)[4:ncol(peakdata)])
  metid = convet_mzmine2mass_dataset(x = peakdataformetid %>% dplyr::select(!compound) ,rt_unit = "minute")
  if(mode=="Pos"){
  metid <-
    annotate_metabolites_mass_dataset(object = metid, 
                                      ms1.match.ppm = 5, 
                                      rt.match.tol = 10001, 
                                      polarity = "positive",
                                      database = get(metid_DB),
                                      #column = "rp_custom",
                                      column = "rp",
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
  if(mode=="Pos"){
  KEGG.compound <- read.csv("/blue/timgarrett/hkates/Garrett/Reporting/Positive_Garrett_MetaboliteStd_Library_RP_edited2022-2-4TJG_KEGG-ALL caps.csv")%>% distinct()
  }else if (mode=="Neg"){
  KEGG.compound <- read.csv("/blue/timgarrett/hkates/Garrett/Reporting/Negative_Garrett_MetaboliteStd_Library_RP_edited211001JGC_KEGG.csv") %>% distinct()
  }
  #KEGG.compound$name <- gsub(" $","",KEGG.compound$name)
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
  metid.result <- metid.result %>% dplyr::rename("compound"="Compound.name")
  peak_annotation <- dplyr::bind_rows(peakdata.KEGG, metid.result)
  #If there is no compound from SECIM or metID, use the mz_rt
  peak_annotation <- peak_annotation %>% mutate(compound = case_when(is.na(compound) ~ paste(mz,rt,sep="_"),.default = compound)) %>%
  mutate(KEGG = case_when(KEGG=="" ~ NA,.default = KEGG)) 
  peak_annotation$mode <- mode
                            
  outputs_list <- list()
  #Merge results with metabolite annotation and name the dataset to be saved with mode
  if(test_type=="t.test"){
    outputs_list[[1]] <- merge(ttest.results,peak_annotation,by="id",all.x=TRUE)
    outputs_list[[1]]$contrast <- gsub("[()]", "",outputs_list[[1]]$contrast)
    outputs_list[[1]] <- outputs_list[[1]] %>% inner_join(means_FC, 
                                                          by=c('id', 'contrast'))
    outputs_list[[1]] <- outputs_list[[1]] %>% relocate(contrast,compound)
    outputs_list[[1]]$adj.p.value <- p.adjust(outputs_list[[1]]$p.value,method="fdr")
    outputs_list[[1]] <- outputs_list[[1]] %>% relocate(adj.p.value, .after = p.value)
    outputs_list[[2]] <- print("Empty for t test")
    
  } 
  if (test_type %in% c("lm","anova","lme")){
    outputs_list[[1]] <- merge(emmeans.results,peak_annotation,by="id",all.x=TRUE)
    #In case emmeans adds "()" due to special chars
    outputs_list[[1]]$contrast <- gsub("[()]", "",outputs_list[[1]]$contrast)
    outputs_list[[1]] <- outputs_list[[1]] %>% inner_join(means_FC, 
                           by=c('id', 'contrast'))
    outputs_list[[1]] <- outputs_list[[1]] %>% relocate(contrast,compound)
    outputs_list[[2]] <- merge(fit.results,peak_annotation,by="id",all.x=TRUE)
    outputs_list[[2]] <- outputs_list[[2]] %>% relocate(compound)
  }
  if (test_type == "nostats"){
  outputs_list[[1]] <- merge(means_FC,peak_annotation,by="id",all.x=TRUE)
  outputs_list[[1]]$contrast <- gsub("[()]", "",outputs_list[[1]]$contrast)
  outputs_list[[1]] <- outputs_list[[1]] %>% relocate(contrast,compound)
  outputs_list[[2]] <- print("Empty for nostats")
  }

  #Merge filtered dataset with metabolite annotation
  #outputs_list[[3]] <- merge(data.frame(t(proc.data)),peak_annotation,by.x=0,by.y="id",all.x=TRUE)
  outputs_list[[3]] <- merge(data.proc,peak_annotation,by="id",all.x=TRUE)
  outputs_list[[3]] <- outputs_list[[3]][c(which(outputs_list[[3]]$id  == "Class"), setdiff(seq_len(nrow(outputs_list[[3]])), which(outputs_list[[3]]$id == "Class"))),]
  #Added in v 12.2 because for >1 metadata category, the above does not move the metadata rows to the top.
  outputs_list[[3]] <- outputs_list[[3]][c(which(is.na(outputs_list[[3]]$id)), setdiff(seq_len(nrow(outputs_list[[3]])),
                                                                                      which(is.na(outputs_list[[3]]$id)))),]
  outputs_list[[3]] <- outputs_list[[3]] %>% relocate(compound)
  
  outputs_list[[4]] <- merge(data.final,peak_annotation,by="id",all.x=TRUE)
  outputs_list[[4]] <- outputs_list[[4]][c(which(outputs_list[[4]]$id  == "Class"), setdiff(seq_len(nrow(outputs_list[[4]])), which(outputs_list[[4]]$id == "Class"))),]
  #Added in v 12.2 because for >1 metadata category, the above does not move the metadata rows to the top.
  outputs_list[[4]] <- outputs_list[[4]][c(which(is.na(outputs_list[[4]]$id)), setdiff(seq_len(nrow(outputs_list[[4]])),
                                                                                       which(is.na(outputs_list[[4]]$id)))),]
  outputs_list[[4]] <- outputs_list[[4]] %>% relocate(compound)
  
  #return plots
  outputs_list[[5]] <- arrangeGrob(plots[["plot_env"]][["p1"]],plots[["plot_env"]][["p3"]],plots[["plot_env"]][["p2"]],plots[["plot_env"]][["p4"]],ncol=2,top=textGrob("Feature View"))
  outputs_list[[6]] <- arrangeGrob(plots[["plot_env"]][["p5"]],plots[["plot_env"]][["p7"]],plots[["plot_env"]][["p6"]],plots[["plot_env"]][["p8"]],ncol=2,top=textGrob("Sample View"))
  if(mode=="Neg"){
  outputs_list[[7]] <- md
  #name the lists
  if(test_type=="t.test"){
    names(outputs_list) <- c(paste0(mode,".ttest.metab"),paste0(mode,"Empty"),paste0(mode,".processed.data"),
                             paste0(mode,".normalized.data"),paste0(mode,".FeatureView"),paste0(mode,".SampleView"),"metadata")
  } else{
    names(outputs_list) <- c(paste0(mode,".fit.results.metab"),paste0(mode,".emmeans.results.metab"),
                             paste0(mode,".processed.data"),
                             paste0(mode,".normalized.data"),paste0(mode,".FeatureView"),paste0(mode,".SampleView"),"metadata")
  }
  }else{
    #outputs_list[[7]] <- md
    #name the lists
    if(test_type=="t.test"){
      names(outputs_list) <- c(paste0(mode,".ttest.metab"),paste0(mode,"Empty"),paste0(mode,".processed.data"),
                               paste0(mode,".normalized.data"),paste0(mode,".FeatureView"),paste0(mode,".SampleView"))
    } else{
      names(outputs_list) <- c(paste0(mode,".emmeans.results.metab"),paste0(mode,".fit.results.metab"),
                               paste0(mode,".processed.data"),
                               paste0(mode,".normalized.data"),paste0(mode,".FeatureView"),paste0(mode,".SampleView"))
    }
    
  }
  return(outputs_list)
}
#pos.output <- outputs_list
#saveRDS(object = pos.output,file="Sumners.posoutput.5182023.Rdata")
