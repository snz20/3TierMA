# Three-Tiered Meta Analysis of Gene Expression Profiles of Co-morbid Diseases
[Massachusetts Institute of Technology](http://web.mit.edu)

[Computer Science and Artificial Intelligence Laboratory](www.csail.mit.edu)

[Theory of Computation Group](http://theory.lcs.mit.edu)

[Berger Lab for Computation and Biology](people.csail.mit.edu/bab/index.html)

Corresponding author for study: bab@mit.edu

Queries regarding usage of three-tiered meta analysis pipeline: nazeen@mit.edu

---
## Synopsis
We introduce a three-tiered meta analysis approach for studying the shared genetics of co-ocurring disease conditions in patients from their gene expression profiles. Due to heterogeneity of the nature of microarray studies performed on different patient populations it is difficult to study them through a single lens. We employ Fisher's classic method of combined probabilities to overcome this difficulty. We employ our pipeline to study the shared biology of Autism Spectrum Disorder (ASD) and its highly prevalent co-morbidities and find a shared innate immunity component among them in the following study:

**Sumaiya Nazeen†, Nathan P. Palmer†, Bonnie Berger* and Isaac S. Kohane**. **_Integrative analysis of genetic datasets reveals a shared innate immune component in autism spectrum disorder and its co-morbidities_**. (accepted for publication in Genome Biology, 2016)

†: equal contribution; *: correspondence

---
## Instructions
1. Use the R code in DiffExprAnalysis-Indiv.R to perform differential expression analysis on individual microarray studies selected under each disease.
2. Use the pyhton code in combineOutputOfDiffExp.py to combine the output files (step 1) per disease.
3. Open the output of step 2 in MS Excel and use Excel tools for calculating Fisher's combined p-values per gene. (We have provided the excel file with the formula zipped under the name fisherComboP_disease_gene.xls.zip.
4. Use the R code in fdrAdjust.R to adjust Fisher's combine p-values calculated in step 3. The input file should consist of two tab-separated columns . Follow the instructions in the source file to generate the list of significant genes for each disease.
5. Use the python code in calcOverlap.py to find the overlap between disease gene set and pathway gene sets. The shell script enrichment.sh contains the shell commands to calculate the overlaps using calcOverlap.py.
6. Use the R code in EnrichmentAnalysis.R on the output from step 5 to perform hypergeometric enrichment analysis of the disease gene sets in the pathway gene sets.
7. Use the pyhton code in combineOutputOfEnrichment.py to combine the output files (step 6) per pathway dataset.
8. Open the output of step 7 in MS Excel and use Excel tools for calculating Fisher's combined p-values per pathway. Filter out the non-significant pathways in ASD and order the pathways in the ascending order of p-values. (We have provided the excel file with the formula under the name fisherComboP_disease_pathway.xls".)

---
## Data
Pathway gene sets used in our study can be downloaded from [here](http://groups.csail.mit.edu/cb/3TierMA/data/Pathways.zip). The expression data used in our study can be downloaded directly from GEO using the "GEOquery" package in R. The accession numbers for the GEO studies are listed [here](http://groups.csail.mit.edu/cb/3TierMA/data/Accessions.docx).

---
## License
This work is published under MIT License. Permission is hereby granted to whoever obtains a copy of the files associated with our pipeline to use, copy, modify, merge, publish, and distribute free of charge. For details see the license file
