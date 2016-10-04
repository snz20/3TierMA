# Three-Tiered Meta Analysis of Gene Expression Profiles of Co-morbid Diseases
[Massachusetts Institute of Technology](http://web.mit.edu)</br>
[Computer Science and Artificial Intelligence Laboratory](www.csail.mit.edu)</br>
[Theory of Computation Group](http://theory.lcs.mit.edu)</br>
[Berger Lab for Computation and Biology](people.csail.mit.edu/bab/index.html)</br>

Corresponding author for study: bab@mit.edu</br>
Queries regarding usage of three-tiered meta analysis pipeline: nazeen@mit.edu</br>

---
## Synopsis
We introduce a three-tiered meta analysis approach for studying the shared genetics of co-ocurring disease conditions in patients from their gene expression profiles. Due to heterogeneity of the nature of microarray studies performed on different patient populations it is difficult to study them through a single lens. We employ Fisher's classic method of combined probabilities to overcome this difficulty. We employ our pipeline to study the shared biology of Autism Spectrum Disorder (ASD) and its highly prevalent co-morbidities and find a shared innate immunity component among them in the following study:

**Sumaiya Nazeen†, Nathan P. Palmer†, Bonnie Berger* and Isaac S. Kohane**. **_Integrative analysis of genetic datasets reveals a shared innate immune component in autism spectrum disorder and its co-morbidities_**. (accepted for publication in Genome Biology, 2016)

†: equal contribution; *: correspondence

---
## Instructions
1. Use the R code in _DiffExprAnalysis-Indiv.R_ to perform differential expression analysis on individual microarray studies selected under each disease.
2. Use the pyhton code in _combineOutputOfDiffExp.py_ to combine the output files (step 1) per disease.
3. Open the output of step 2 in MS Excel and use Excel tools for calculating Fisher's combined p-values per gene. We have provided the excel file with the formula zipped under the name _fisherComboP_disease_gene.xls.zip_.
4. Use the R code in _fdrAdjust.R_ to adjust Fisher's combine p-values calculated in step 3. The input file should consist of two tab-separated columns. Follow the instructions in the source file to generate the list of significant genes for each disease.
5. Use the python code in _calcOverlap.py_ to find the overlap between disease gene set and pathway gene sets. The shell script _enrichment.sh_ contains the shell commands to calculate the overlaps using _calcOverlap.py_.
6. Use the R code in _EnrichmentAnalysis.R_ on the output from step 5 to perform hypergeometric enrichment analysis of the disease gene sets in the pathway gene sets.
7. Use the pyhton code in _combineOutputOfEnrichment.py_ to combine the output files (step 6) per pathway dataset.
8. Open the output of step 7 in MS Excel and use Excel tools for calculating Fisher's combined p-values per pathway. Filter out the non-significant pathways in ASD and order the pathways in the ascending order of p-values. We have provided the excel file with the formula under the name _fisherComboP_disease_pathway.xls_.
9. For the combined significant pathways (adjusted p-value < 0.05), calculate minimum Bayes factor and minimum posterior probability of null hypothesis individually in each disease as well as in the combined scenario as described in the Methods section. We have provided the exel file with the formulae under the name _BFandPosteriorProbabilityCalculation.xlsx_.

---
## Data
Pathway gene sets used in our study are provided in the _Pathways.zip_ file. The expression data used in our study can be downloaded directly from GEO using the "GEOquery" package in R. The accession numbers for the GEO studies are listed in the _Accessions.docx_ file.

---
## License
This work is published under MIT License. Permission is hereby granted to whoever obtains a copy of the files associated with our pipeline to use, copy, modify, merge, publish, and distribute free of charge. For details see the license file.
