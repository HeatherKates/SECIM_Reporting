# Load the pathview library
library(pathview)
library(KEGGREST)

# Function to get pathway categories
get_pathway_categories <- function(pathways) {
  pathway_info <- keggGet(pathways)
  categories <- sapply(pathway_info, function(info) {
    if ("CLASS" %in% names(info)) {
      return(info$CLASS)
    } else {
      return(NA)
    }
  })
  names(categories) <- pathways
  return(categories)
}

# Example compound data
up_compounds=c("CE1554","C00185","C02710","C00815")
down_compounds=c("C01718","C00711","C00303","C00328")
#Add +1 or -1 to each compound
cpd_data <- c(setNames(rep(1, length(up_compounds)), up_compounds),
              setNames(rep(-1, length(down_compounds)), down_compounds))

# Function to get pathways for a given compound
get_pathways <- function(compound_id) {
  pathway_list <- keggLink("pathway", paste0("cpd:", compound_id))
  return(gsub("path:", "", pathway_list))
}

# Get pathways for each compound
compound_pathways <- unique(unlist(lapply(names(cpd_data), get_pathways)))

# Get pathway categories
pathway_categories <- get_pathway_categories(compound_pathways)

# Filter out global pathways
global_categories <- c("Metabolism", "Genetic Information Processing", "Environmental Information Processing",
                       "Cellular Processes", "Organismal Systems", "Human Diseases", "Drug Development",
                       "Metabolic pathways", "Biosynthesis of secondary metabolites",
                       "Microbial metabolism in diverse environments", "Carbon metabolism", "Biosynthesis of antibiotics")
filtered_pathways <- names(pathway_categories)[!sapply(pathway_categories, function(cat) {
  any(global_categories %in% cat)
})]

# Filter out global pathways
filtered_pathways <- names(pathway_categories)[!sapply(pathway_categories, function(cat) {
  any(global_categories %in% cat)
})]

# Add pathways without categories to filtered pathways if needed
# Here we exclude them, but you can include them if you decide so
# Identify pathways without categories
pathways_without_categories <- names(pathway_categories)[is.na(pathway_categories)]
pathways_with_categories <- names(pathway_categories)[!is.na(pathway_categories)]
filtered_pathways <- filtered_pathways[!filtered_pathways %in% pathways_without_categories]


# Convert pathways to "hsa" format
filtered_pathways <- gsub("map","hsa",filtered_pathways)

# Function to visualize pathways
visualize_pathway <- function(pathway_id, cpd_data) {
  tryCatch({
    pathview(gene.data = NULL, cpd.data = cpd_data, pathway.id = pathway_id, species = "hsa")
  }, error = function(e) {
    message(sprintf("Error in pathway %s: %s", pathway_id, e$message))
  })
}

# Loop through each filtered pathway and visualize
for (pathway_id in filtered_pathways) {
  visualize_pathway(pathway_id, cpd_data)
}

