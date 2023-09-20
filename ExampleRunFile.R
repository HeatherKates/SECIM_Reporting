library(rmarkdown)
source("SECIM_Reporting/R/Generate_Reporting_Inputs.R")

#test_type = "t.test","anova","lm","lme","nostats"
client="Dudeja.Serum"
#anova_formula= Metabolite ~ Class
#Create the Reporting Input
ReportInput <- Generate_Report_Inputs(client=client,samples_to_drop=NULL,mzmine_version=2,
                                               ReferenceLevel="Serum_WT_pre",
                                      Input="SECIM_Reporting/InputFiles/Dudeja-Metabolomics_Serum.mini.anova.xlsx",
                                               contrast_var="Class",num_meta=1,
                                      anova_formula= id ~ Class,lm_model=NULL,
                                               test_type="anova",subset=NULL,
                                      metid_DB_file="SECIM_Reporting/InputFiles/kegg_ms1_database0.0.3.rda")
#saveRDS(ReportInput,file=paste0(client,".ReportingInput.",Sys.Date(),".RDATA"))

#Parameters for Report Generator

#ReportInput <- readRDS("Dudeja.Serum.ReportingInput.2023-08-25.RDATA")
Grouping_Variable="Class"
filter_method="IQR"
norm_method="Sum"
ind_var="Class"
num_of_metadata <- 2 #includes sample names
num_samples=9
num_groups=3
contrast_var="Class"
boxplot_var=~Class
#class_order <- levels(as.factor(Client_Data_Download[["metadata"]]$Class))
test_type="anova" #"t.test", "anova", "lmm","repeated_measures_anova","nostats"
#class_order <- c("Serum_WT","Serum_KO")
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
render("SECIM_Reporting/R/REPORT_GENERATOR.Rmd", output_file = paste0(client,".Report.html"))
