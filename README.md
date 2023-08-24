# SECIM_Reporting
These files (future R package, unnamed) are used to analyze and report metabolomics data created by Southeast Center fo Integrated Metabolomics.

The example run file ExampleRunFile.R can be run using `Rscript ExampleRunFile.R` to generate an html report including downloadable data and statistical test results

# The NoStats Branch 

The goal of this branch is to add required options to allow for generating a report when no statistical analysis is performed. In the SECIM_Metabolomics.R function, fold-changes should be calculated for all pairwise samples' comparisons, and the list of top metabolites for the report should be those with the greatest fold-changes for each pairwise samples' contrast- but the central reporting etc. will be most like ANOVA+emmeans, because we want to see all samples in the plots etc.