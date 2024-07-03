# Create vector consisting of compounds for enrichment analysis 
rm(mSet)
mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec <- Client_Data_Download[["report_results"]] %>% filter(contrast==contrasts[[i]]) %>% filter(!!sym(p_type)<0.05) %>% filter(ID_confidence=="High") %>% dplyr::select(compound) %>% unlist()
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)
## [1] "Loaded files from MetaboAnalyst web-server."
# Plot the ORA, bar-graph
mSet<-PlotORA(mSet, paste0("/blue/timgarrett/hkates/SECIM_Reporting/",client,"MSEA.", "png"), dpi=300, width=NA)
#code for R markdown ![MSEA results for significantly changed compounds](paste(client,"MSEA.", "png"))
result_tbl <- read.csv("msea_ora_result.csv")
colnames(result_tbl) <- c("SMPDB Pathway",colnames(result_tbl)[2:ncol(result_tbl)])
# List of files to remove
files_to_remove <- c("msea_ora_result.csv", "name_map.csv", "syn_nms.qs", 
                     "tosend.rds", "compound_db.qs")

# Remove specific files
file.remove(files_to_remove)
