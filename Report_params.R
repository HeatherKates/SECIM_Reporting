ReportInput <- "Dudeja.Fecal.ReportingInput2023-08-17.RDATA"

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
class_order <- c("Fecal_WT","Fecal_KO")

PI <- ""
Institution <- ""
Department <- ""
StudyContact <-"" 
Project <- ""
StudyTitle <- ""
Hypothesis <- ""
StudySummary <- ""

drop_compounds <- c("Acetone", "Sodium cyanide", "Vinyl Chloride", "sodium fluoride", "Magnesium hydroxide", "Acetone cynohydrin","Acetaldehyde","Sulfate","Cyanamide","EDTA","Mercury(2+)","Magnesium oxide","Chloric acid","Sinapate","Hydrochloric acid","Formate","Benzene","Potassium dichromate","Calcium formate","Radium-228","CO","Chloroform","Propan-2-ol","H2O","Acetate","Acetylene","Dichloromethane")