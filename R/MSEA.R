library(MetaboAnalystR)
source("/blue/timgarrett/hkates/SECIM_Reporting/R/calculatehyperscore.R")
source("/blue/timgarrett/hkates/SECIM_Reporting/R/general_data_utils.R")
source("/blue/timgarrett/hkates/SECIM_Reporting/R/setcurrentmsetlib.R")
.on.public.web=FALSE

## When input is a list

# Create vector consisting of compounds for enrichment analysis 
rm(mSet)
my.vec <- Client_Data_Download[["report_results"]] %>% filter(contrast==contrasts[[i]]) %>% filter(!!sym(p_type)<0.05) %>% filter(ID_confidence=="High") %>% dplyr::select(compound) %>% unlist()
# Remove anything after "+" and anything within "()"
#my.vec <-  gsub("\\+.*", "", my.vec)
#my.vec <- gsub("\\(.*?\\)", "", my.vec)



#my.vec <- readRDS("test_my_vec.RDATA")
# Create mSetObj
mSet<-InitDataObjects("list", "msetora", FALSE)
## [1] "MetaboAnalyst R objects initialized ..."
#Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, my.vec);

# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name");
# Example compound name map
# Create the mapping results table
mSet<-CreateMappingResultTable(mSet)

# Set the metabolite filter
mSet<-SetMetabolomeFilter(mSet, F);

# Select metabolite set library, refer to 
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 0)
#mSetTest <- mSet
# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)
## [1] "Loaded files from MetaboAnalyst web-server."
# Plot the ORA, bar-graph
mSet<-PlotORA(mSet, paste0("/blue/timgarrett/hkates/SECIM_Reporting/",client,"MSEA.", "png"), 300, width=NA)
#code for R markdown ![MSEA results for significantly changed compounds](paste(client,"MSEA.", "png"))
# List of files to remove
files_to_remove <- c("msea_ora_result.csv", "name_map.csv", "syn_nms.qs", 
                     "tosend.rds", "compound_db.qs")

# Remove specific files
file.remove(files_to_remove)