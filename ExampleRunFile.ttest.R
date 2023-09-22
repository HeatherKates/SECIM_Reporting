library(rmarkdown)
source("SECIM_Reporting/R/Generate_Reporting_Inputs.R")

client="Example"
#anova_formula= Metabolite ~ Class
#Create the Reporting Input
ReportInput <- Generate_Report_Inputs(client=client,samples_to_drop=NULL,mzmine_version=2,
                                               ReferenceLevel="WT",
                                      Input="SECIM_Reporting/InputFiles/Example.ttest.xlsx",
                                               contrast_var="Class",num_meta=1,
                                      anova_formula= id ~ Class,lm_model=NULL,
                                               test_type="t.test",subset=NULL,
                                      metid_DB_file="SECIM_Reporting/InputFiles/kegg_ms1_database0.0.3.rda")

#Parameters for Report Generator

Grouping_Variable="Class"
filter_method="IQR"
norm_method="Sum"
ind_var="Class"
num_of_metadata <- 2 #includes sample names
num_samples=6
num_groups=2
contrast_var="Class"
boxplot_var=~Class
#class_order <- levels(as.factor(Client_Data_Download[["metadata"]]$Class))
test_type="t.test" #"t.test", "anova", "lmm","repeated_measures_anova","nostats"
#class_order <- c("Serum_WT","Serum_KO")
drop_compounds <- c("Sodium bicarbonate")

PI <- "Example PI name"
Institution <- "Example Institution Name"
Department <- "Example Department"
StudyContact <-"Example Study Contact" 
Project <- "Example Project"
StudyTitle <- "Example Study Title"
Hypothesis <- "Example Hypothesis"
StudySummary <- "Example Study Summary"
SampleType <- "Example Sample Type"

#Run the Report Generator
render("SECIM_Reporting/R/REPORT_GENERATOR.Rmd", output_file = paste0(client,".Report.html"))
