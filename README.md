# noobomics
A very quick and unsophisticated approach to exploring RNAseq and Proteomic data to distill out hypothesis-driven candidates.

To install this package simply use devtools and install.github() as such:

library(devtools)
install_github("KyleStiers/noobomics/Desktop/noobomics/noobomics")

Prepare your data in excel or whatever you prefer. This package requires the frames be already sorted (Sort function in Excel) by ID - as long as they are treated the same the package will work. As of now IDs must match exactly to count, in the future I'll likely add some substring searching.

The 1st data frame should look something like this in the csv:
ID, rna_WT, rna_Mutant1, rna_Mutant2
abcd,10,20,30

The second frame should be exactly the same format but with unique headers:
ID, prot_WT, prot_Mutant1, prot_Mutant2
abcd, 2, 5, 7

Then once you have your two dataframes (df1, df2) loaded in simply run:

noobomics(df1, df2)

This will plot the log2() transformed data as density distributions, scatterplots with Pearson correlation coefficients, a heatmap (with default clustering), and finally the biplot of the Principal Component Analysis.

I have many features planned for the future. Please feel free to also suggest things, give notes.
