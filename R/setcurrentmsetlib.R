SetCurrentMsetLib <- function(mSetObj=NA, libname, excludeNum=0){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(libname=="self"){
    ms.list <- mSetObj$dataSet$user.mset;
    ms.list <- lapply(ms.list, function(x) unique(unlist(strsplit(x, "; ", fixed=TRUE))));
    current.msetlib <- vector("list", 3)
    names(current.msetlib) <- c("name", "member", "reference")
    mSetObj$analSet$msetlibname <- libname;
  } else {
    if(!.on.public.web & grepl("kegg", libname)){ # api only for KEGG msets
      mSetObj$api$libname <- libname
      mSetObj$api$excludeNum = excludeNum
      mSetObj$analSet$msetlibname <- libname
      return(.set.mSet(mSetObj));
    }
    
    # create a named list, use the ids for list names
    # https://github.com/xia-lab/MetaboAnalystR/issues/172
    ms.list <- iconv(current.msetlib[1], from = 'utf8', to = 'utf8');
    ms.list <- lapply(ms.list, function(x) unique(unlist(strsplit(x, "; ", fixed=TRUE))));
    names(ms.list) <- current.msetlib[1];
  }
  
  # total uniq cmpds in the mset lib
  mSetObj$dataSet$uniq.count <- length(unique(unlist(ms.list, use.names = FALSE)));
  
  # update current.mset and push to global env
  current.msetlib$member <- ms.list;
  
  current.msetlib <<- current.msetlib;
  return(.set.mSet(mSetObj));
}

