# Create vector consisting of compounds for enrichment analysis 
rm(mSet)
mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec <- Client_Data_Download[["report_results"]] %>% filter(contrast==contrasts[[i]]) %>% filter(!!sym(p_type)<0.05) %>% filter(ID_confidence=="High") %>% dplyr::select(compound) %>% unlist()
#remove adducts
# Function to remove common adducts and clean compound names
# Define a comprehensive list of adducts, ordered from longest to shortest
adducts <- c("_M\\+HCOO", "_M\\+CH3COO", "_M\\+NH4", "_M\\+Na", "_M-H", "_M\\+K", "_M\\+H", "\\+HCOO", "\\+CH3COO", "\\+NH4", "\\+Na", "\\+K", "\\+H", "-H", "\\+Cl", "-Cl")

# Function to remove adducts and preceding symbols
remove_adducts <- function(compound, adduct_list) {
  for (adduct in adduct_list) {
    # Create a pattern to match the adduct and the preceding "+" or "-" symbol
    pattern <- paste0("([+_-])?", adduct, "$")
    
    # Remove the adduct and preceding symbol if found
    if (grepl(pattern, compound)) {
      compound <- sub(pattern, "", compound)
      break  # Exit loop after first match to prevent further removal
    }
  }
  return(compound)
}

# Apply the function to the vector
cleaned_cmpd.vec <- sapply(cmpd.vec, remove_adducts, adduct_list = adducts)

# View the cleaned vector
cleaned_cmpd.vec

# View the cleaned vector
cmpd.vec <- cleaned_cmpd.vec
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
