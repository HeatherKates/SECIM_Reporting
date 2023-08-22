#Update this file
#library(rmarkdown)
#render("REPORT_GENERATOR.Rmd", output_file = "Dudeja.Fecal.Report.html")
#When finished, save this file AS using client name for reproducibility (must be named Report_params.R at runtime)

ReportInput <- "Dudeja.Serum.ReportingInput2023-08-17.RDATA"

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

PI <- "Dudeja"
Institution <- ""
Department <- ""
StudyContact <-"" 
Project <- ""
StudyTitle <- ""
Hypothesis <- ""
StudySummary <- ""

drop_compounds <- c("Sodium bicarbonate")