Grouping_Variable="Class"
filter_method="IQR"
norm_method="Sum"
Input <- "Client_Input_Sheets/Dudeja-Metabolomics_Fecal.xlsx"
ind_var="Class"
num_of_metadata <- 2 #includes sample names
num_samples=10
num_groups=2
contrast_var="Class"
boxplot_var=~Class
load("Dudeja.Fecal.ReportingInput2023-08-17.RDATA")
class_order <- levels(as.factor(Client_Data_Download[["metadata"]]$Class))
test_type="t.test" #t.test, anova, lmm,"repeated_measures_anova"