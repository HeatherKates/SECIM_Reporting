#'Sanity Check Data
#'@description SanityCheckData is used for data processing, and performs a basic sanity 
#'check of the uploaded content, ensuring that the data is suitable for further analysis. 
#'The function will return a message if the data has successfully passed the check
#'and is deemed suitable for further analysis. If it fails, the function will return a 0.
#'The function will perform the check directly onto the mSet$dataSet object, and must 
#'be performed immediately after reading in data. 
#'The sanity check function evaluates the accuracy of sample and class labels, data structure, 
#'deals with non-numeric values, removes columns that are constant across all samples (variance = 0), 
#'and by default replaces missing values with half of the original minimal positive value in your dataset.
#'@usage SanityCheckData(mSetObj=NA)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@import qs
#'@export
#'
SanityCheckDataHRK <- function(mSetObj=NA){
  mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj);
  if(file.exists("data_orig.qs")){  
    orig.data <- qs::qread("data_orig.qs");
  } else {
    return(0);
  }  
  msg <- NULL;
  cls <- mSetObj$dataSet$orig.cls;
  mSetObj$dataSet$small.smpl.size <- 0;
  
  # check class info only for one factor data
  # For "mf", there is a dedicated page/step "SanityCheckMeta" for this
  
  if(mSetObj$dataSet$cls.type == "disc"){
    
    # added mSetObj$dataSet$pair.checked to allow edit group function names not overwritten by original files
    if(mSetObj$dataSet$paired & !(mSetObj$dataSet$pair.checked)){ 
      msg<-c(msg,"Samples are paired.");
      # need to first set up pair information if not csv file
      if(!(mSetObj$dataSet$type=="conc" | mSetObj$dataSet$type=="specbin" | mSetObj$dataSet$type=="pktable" )){
        pairs <- ReadPairFile();
        # check if they are of the right length
        if(length(pairs)!=length(mSetObj$dataSet$url.smp.nms)){
          AddErrMsg("Error: the total paired names are not equal to sample names.");
          return(0);
        }else{
          # matching the names of the files
          inx<-match(rownames(orig.data), names(pairs));
          #check if all matched exactly
          if(sum(is.na(inx))>0){
            AddErrMsg("Error: some paired names not match the sample names.");
            return(0);
          }else{
            mSetObj$dataSet$pairs <- pairs[inx];
          }
        }
      }
      
      pairs <- mSetObj$dataSet$pairs;
      
      # check if QC samples are present
      qc.hits <- tolower(as.character(cls)) %in% "qc";
      if(sum(qc.hits) > 0){
        AddErrMsg("<font color='red'>Error: QC samples not supported in paired analysis mode.</font>");
        AddErrMsg("You can perform QC filtering using regular two-group labels.");
        AddErrMsg("Then re-upload your data (without QC samples) for paired analysis.");
        return(0);
      }else{
        pairs <- as.numeric(pairs);
      }  
      
      label <- as.numeric(pairs);
      cls <- as.factor(ifelse(label>0,1,0));
      mSetObj$dataSet$pairs <- label;       
      lev <- unique(pairs);
      uni.cl <- length(lev);
      uni.cl.abs <- uni.cl/2;             
      sorted.pairs <- sort(pairs,index=TRUE);
      
      if(!all(sorted.pairs$x==c(-uni.cl.abs:-1,1:uni.cl.abs))){
        AddErrMsg("There are some problems in paired sample labels! ");
        if(uni.cl.abs != round(uni.cl.abs)){
          duplicates <- pairs[duplicated(pairs)]
          dup.msg <- paste0("Duplicated labels:", duplicates)
          AddErrMsg(paste("The total samples must be of even number!", dup.msg));
        }else{
          AddErrMsg(paste("And class labels between ",-uni.cl.abs,
                          " and 1, and between 1 and ",uni.cl.abs,".",sep=""));
        }
        return(0);
      } 
      
      msg <- c(msg,"The labels of paired samples passed sanity check.");
      msg <- c(msg, paste("A total of", uni.cl.abs, "pairs were detected."));
      # make sure paired samples are sorted 1:n/2 and -1:-n/2
      
      x<-sorted.pairs$ix[(uni.cl.abs+1):uni.cl]
      y<-sorted.pairs$ix[uni.cl.abs:1]
      index<-as.vector(cbind(x,y));
      cls<-cls[index];
      pairs <- pairs[index];
      orig.data<- orig.data[index,];
      
      mSetObj$dataSet$pairs <- pairs;
      mSetObj$dataSet$orig.cls <- cls;
      
      #add sync for paired names
      mSetObj$dataSet$url.smp.nms <- mSetObj$dataSet$url.smp.nms[index];
      
      mSetObj$dataSet$pair.checked <- TRUE;
      #qs::qsave(orig.data, file="data_orig.qs");
      
    } else {
      
      # check for class labels at least two replicates per class but QC and BLANK
      
      cls.lbl <- mSetObj$dataSet$orig.cls;
      qb.inx <- tolower(cls.lbl) %in% c("qc", "blank");
      if(sum(qb.inx) > 0){
        cls.Clean <- as.factor(as.character(cls.lbl[!qb.inx])); # make sure drop level
      } else {
        cls.Clean <- cls.lbl;
      }
      
      # allow it pass to sanity check and correct there
      if (FALSE) {
      if(anal.type != "network" & anal.type != "mf"){ # add exception for DSPC correlation network 
        if(min(table(cls.Clean)) < 2 | length(levels(cls.Clean)) < 2){
          AddErrMsg(paste ("A total of", length(levels(cls.Clean)), "groups found with", length(cls.Clean), "samples."));
          AddErrMsg("<font color='red'>At least <b>two</b> groups and <b>two replicates</b> per group are required for analysis</font>!");
          if(length(levels(cls.Clean)) > 10){
            AddErrMsg("<font color='red'>It seems the number of groups is big. Make sure to specify the correct format (i.e. samples in <b>columns</b> or <b>rows</b>) in the Data Upload page</font>");
            return(-2);
          }else{
            AddErrMsg("You can click the <b>Edit Groups</b> button below to see the group labels for each sample and make corrections.");
            return(-1);
          }
        }
      }
    } else if(anal.type == "mf"){
        if(min(table(cls.Clean)) < 3 | length(levels(cls.Clean)) < 2){
          msg <- c(msg, paste ("A total of", length(levels(cls.Clean)), "groups found with", length(cls.Clean), "samples."));
          msg <- c(msg, "The primary factor is highly possible a continuous variable.")
        }
      }
      
      if("NMDR_id" %in% names(mSetObj$dataSet)){
        msg <- c(msg, paste("Study", mSetObj$dataSet$NMDR_id, "was successfully downloaded from the Metabolomics Workbench!"))
      }
      if(!mSetObj$dataSet$paired){
        msg <- c(msg,"Samples are not paired.");
      }else{
        msg <- c(msg,"Samples are paired.");
      }
    }
    
    # checking if too many groups but a few samples in each group
    cls.lbl <- mSetObj$dataSet$orig.cls;
    # need to exclude QC or blank
    qb.inx <- tolower(cls.lbl) %in% c("qc", "blank");
    if(sum(qb.inx) > 0){
      cls.lbl <- as.factor(as.character(cls.lbl[!qb.inx])); # make sure drop level
    }
    min.grp.size <- min(table(cls.lbl));
    cls.num <- length(levels(cls.lbl));
    if((cls.num/min.grp.size > 3) & (anal.type != "mf")){
      mSetObj$dataSet$small.smpl.size <- 1;
      msg <- c(msg, "<font color='red'>Too many groups with very small number of replicates!</font>");
      msg <- c(msg, "<font color='red'>Only a subset of methods will be available for analysis!</font>");
    }
    
    
    msg <- c(msg, paste(cls.num, "groups were detected in samples."));
    
    
    if("NMDR_id" %in% names(mSetObj$dataSet)){
      msg <- c(msg, paste("Study", mSetObj$dataSet$NMDR_id, "group labels:", paste0(unique(cls.lbl), collapse = ", ")))
    }
    
    mSetObj$dataSet$cls.num <- cls.num;
    mSetObj$dataSet$min.grp.size <- min.grp.size;
    
    
    ord.inx <- order(mSetObj$dataSet$orig.cls);
    mSetObj$dataSet$orig.cls <- cls[ord.inx];
    mSetObj$dataSet$url.smp.nms <- mSetObj$dataSet$url.smp.nms[ord.inx];
    orig.data <- orig.data[ord.inx, , drop=FALSE];
    qs::qsave(orig.data, file="data_orig.qs");
    if(mSetObj$dataSet$paired){
      mSetObj$dataSet$pairs <- mSetObj$dataSet$pairs[ord.inx];
    }
    
  }
  msg<-c(msg,"Only English letters, numbers, underscore, hyphen and forward slash (/) are allowed.");
  msg<-c(msg,"<font color=\"orange\">Other special characters or punctuations (if any) will be stripped off.</font>");
  
  int.mat <- orig.data;
  
  if(ncol(int.mat)==1){
    if(anal.type=="roc"){
      mSetObj$dataSet$roc_cols <- 1;
    } else {
      AddErrMsg("<font color='red'>One-column data is only supported for biomarker analysis.</font>");
      return(0);
    }
  } else {
    mSetObj$dataSet$roc_cols <- 2;
  }
  
  # check numerical matrix
  rowNms <- rownames(int.mat);
  colNms <- colnames(int.mat);
  naNms <- sum(is.na(int.mat));
  
  for(c in 1:ncol(int.mat)) {
    if(class(int.mat[,c]) == "integer64"){
      int.mat[,c] <- as.double(int.mat[,c]);
    }
  }
  
  num.mat <- apply(int.mat, 2, as.numeric)
  
  if(sum(is.na(num.mat)) > naNms){
    # try to remove "," in thousand seperator if it is the cause
    num.mat <- apply(int.mat,2,function(x) as.numeric(gsub(",", "", x)));
    if(sum(is.na(num.mat)) > naNms){
      msg<-c(msg,"<font color=\"red\">Non-numeric values were found and replaced by NA.</font>");
    }else{
      msg<-c(msg,"All data values are numeric.");
    }
  }else{
    msg<-c(msg,"All data values are numeric.");
  }
  
  int.mat <- num.mat;
  rownames(int.mat) <- rowNms;
  colnames(int.mat)<- colNms;
  
  # check for columns with all constant (var =0)
  varCol <- apply(int.mat, 2, var, na.rm=T);
  
  constCol <- (varCol == 0 | is.na(varCol));
  constNum <- sum(constCol, na.rm=T);
  if(constNum > 0){
    msg<-c(msg, paste("<font color=\"red\">", constNum, "features with a constant or single value across samples were found and deleted.</font>"));
    int.mat <- int.mat[,!constCol, drop=FALSE];
  }
  
  # check zero, NA values
  totalCount <- nrow(int.mat)*ncol(int.mat);
  naCount <- sum(is.na(int.mat));
  naPercent <- round(100*naCount/totalCount,1)
  #  print(naCount)
  mSetObj$dataSet$missingCount <- naCount;
  
  msg<-c(msg, paste("A total of ", naCount, " (", naPercent, "%) missing values were detected.", sep=""));
  msg<-c(msg, "<u>By default, missing values will be replaced by 1/5 of min positive values of their corresponding variables</u>",
         "Click the <b>Proceed</b> button if you accept the default practice;",
         "Or click the <b>Missing Values</b> button to use other methods.");
  
  mSetObj$dataSet$proc.cls <- mSetObj$dataSet$cls <- mSetObj$dataSet$orig.cls;
  
  if(is.null(mSetObj$dataSet$meta.info)){
    mSetObj$dataSet$meta.info <- data.frame(mSetObj$dataSet$cls);
    colnames(mSetObj$dataSet$meta.info) = "Class";
  }
  
  # make sure the meta.info is synchronized with data
  ##if(substring(mSetObj$dataSet$format,4,5)=="mf"){
  ##  my.sync <- .sync.data.metadata(int.mat, mSetObj$dataSet$meta.info);
  ##  int.mat <- my.sync$data;
  ##  mSetObj$dataSet$meta.info <- my.sync$metadata;
  ##}
  
  qs::qsave(as.data.frame(int.mat), "preproc.qs");
  
  mSetObj$msgSet$check.msg <- c(mSetObj$msgSet$read.msg, msg);
  

    print(c("Successfully passed sanity check!", msg))
  return(MetaboAnalystR:::.set.mSet(mSetObj));
}

