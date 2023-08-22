library(rmarkdown)
source("R/Generate_Reporting_Inputs.R")

#Create the Reporting Input
ReportInput <- Generate_Report_Inputs(client="Dudeja.Serum",samples_to_drop=NULL,mzmine_version=2,
                                               ReferenceLevel="Serum_WT",Input="InputFiles/Dudeja-Metabolomics_Serum.xlsx",
                                               contrast_var="Class",num_meta=1,anova_formula=NULL,lm_model=NULL,
                                               test_type="t.test",subset=NULL,metid_DB_file="/blue/timgarrett/hkates/Garrett/Reporting/kegg_ms1_database0.0.3.rda")
saveRDS(ReportInput,file="Dudeja.Serum.ReportingInput2023-08-17.RDATA")

#Parameters for Report Generator

#ReportInput <- "Dudeja.Serum.ReportingInput2023-08-17.RDATA"
Grouping_Variable="Class"
filter_method="IQR"
norm_method="Sum"
ind_var="Class"
num_of_metadata <- 2 #includes sample names
num_samples=10
num_groups=2
contrast_var="Class"
boxplot_var=~Class
#class_order <- levels(as.factor(Client_Data_Download[["metadata"]]$Class))
test_type="t.test" #t.test, anova, lmm,"repeated_measures_anova"
class_order <- c("Serum_WT","Serum_KO")
drop_compounds <- c("Sodium bicarbonate")

PI <- "Dudeja"
Institution <- ""
Department <- ""
StudyContact <-"" 
Project <- ""
StudyTitle <- ""
Hypothesis <- ""
StudySummary <- ""

#Run the Report Generator
render("R/REPORT_GENERATOR.Rmd", output_file = "Dudeja.Serum.Report.html")
#When finished, save this file AS using client name for reproducibility (must be named Report_params.R at runtime)