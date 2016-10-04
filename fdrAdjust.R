## Code to perform FDR adjustment on Fisher's combined p-values of disease genes
## Written by Sumaiya Nazeen, nazeen@mit.edu
## R version 2.15.1
## Platform: 3.13.0-35-generic #62-Ubuntu (64-bit)

## Input file is two-column file of the format <gene, fisher_p>. The columns must have headers.
## Output file is a three-column file of the format <gene, fisher_p, BY_adjusted_p>. 
## Note that, the output file will have quote marks around the column headers and gene symbols.
## Shell command to remove quotes from the output file: sed -i 's/"//g' FILENAME  
## Select the significant genes using the cutoff: BY_adjusted_p < 0.05
## Shell command to select significant genes using the cutoff: 
##### awk -F'\t' 'NR>1{OFS="\t";if($3<0.05) print $1,$3}' all_genes_adjp.tab > sig_genes_by.tab

## Define Utility Function

fdrAdjust <- function(inFile, outFile){
        x = read.table(file=inFile, header=TRUE)
        colnames(x) = c("gene","fisher_p")
        BY_adjusted_p = p.adjust(x[,"fisher_p"],"BY") #for other adjustments replace "BY" with one of {"bonferroni","BH","none"}
        z = cbind(x, BY_adjusted_p)
        write.table(z, file=outFile, sep='\t', row.names=FALSE, col.names=TRUE)
}

