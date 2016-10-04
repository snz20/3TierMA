## Code to run differential expression analysis on individual microarray study selected under each disease
## Written by Sumaiya Nazeen, nazeen@mit.edu
## R version 2.15.1
## Platform: 3.13.0-35-generic #62-Ubuntu (64-bit)

## This code will allow one to reproduce the differential expression analyses in Nazeen et al., 2015

##### Goes through the following steps:
## For each of the microarray studies selected under a disease
##### 1) Get the datasets from GEO
##### 2) Obtain updated annotations
##### 3) Perform differential expression analysis with limma
##### 4) Annotate probes by gene symbols and output

## Obtain or load the necessary libraries

## Uncomment these lines if you do not have these packages
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("Biobase", "GEOquery"))

library("Biobase")
library("GEOquery")
library("limma")

## Define Utility Functions

readGEO <- function(seriesName, annotName){
	gset <- getGEO(seriesName, GSEMatrix =TRUE)
	if (length(gset) > 1) # check if multiple series
		idx <- grep(annotName, attr(gset, "names")) 
	else 
		idx <- 1
	gset <- gset[[idx]]

	# make proper column names to match toptable 
	fvarLabels(gset) <- make.names(fvarLabels(gset))

	# log2 transform - must do for limma
	ex <- exprs(gset)
	qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)) # remove NAs and count samples in quantiles
	LogC <- (qx[5] > 100) ||
			  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
			  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
	if (LogC) { ex[which(ex <= 0)] <- NaN
	  exprs(gset) <- log2(ex) }

	return(gset)
}

labelSamples <- function(gset, colName, gTerm0, gTerm1){
	pdat <- pData(gset)
	eset <- exprs(gset)
	sml <- rep(0,length(colnames(eset)))
	sml[grep(gTerm0,t(pdat[colName]))] <- "G0"
	sml[grep(gTerm1,t(pdat[colName]))] <- "G1"
	label <- as.factor(sml)
	
	return(label)
}

## This method is used for performing the differential expression analysis.
## The parameter "method" can be one of the following: {"bonferroni","BY","BH","none"}

analyseSeparate <- function(gset, label, method){
	# set up the data and proceed with analysis
	gset$description <- label # assigning labels to the samples
	design <- model.matrix(~ description + 0, gset) # generate design matrix
	colnames(design) <- levels(label)
	fit <- lmFit(gset, design)
	cont.matrix <- makeContrasts(G0-G1, levels=design) # construct contrast matrix
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2, 0.01)
	tT <- topTable(fit2, adjust=method, sort.by="B", number=dim(fit2)[1])
	return(tT)
}

# Annotate with NCBI supplied annotations
annotateAndOutput <- function(gset, tT, outFile){
	# load NCBI platform annotation
	gpl <- annotation(gset)
	platf <- getGEO(gpl, AnnotGPL=TRUE)
	ncbifd <- data.frame(attr(dataTable(platf), "table"))

	# replace original platform annotation
	tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
	tT <- merge(tT, ncbifd, by="ID")
	tT <- tT[order(tT$P.Value), ]  # restore correct order

	tT <- subset(tT, select=c("ID","Gene.symbol","adj.P.Val","P.Value","t","B","logFC","Gene.title"))
	tT <- tT[tT[,2] != "",];
	write.table(tT, file=outFile, row.names=F, sep="\t")
}

## Get the dataset for each study under disease and perform differential expression analysis.
## The following code uses "BY" adjustment. We have also performed the analyses using Bonferroni, BH adjustment and no adjustment.

## Asthma

a1 <- readGEO("GSE19187","GPL6244")
pa1 <- pData(a1)
filInd <- c(grep("asthma",t(pa1["title"])),grep("Healthy",t(pa1["title"]))) # filter asthma and control samples
fa1 <- a1
exprs(fa1) <- exprs(a1)[,filInd]
pData(fa1) <- pData(a1)[t(filInd),]
la1 <- labelSamples(fa1,"title","asthma","Healthy")
ta1 <- analyseSeparate(fa1, la1, "BY") # to identify differentially expressed genes with other FDR adjustments, replace "BY" with appropriate test 
annotateAndOutput(fa1, ta1, "diff_a1.tab")

a2 <- readGEO("GSE27011","GPL6244")
la2 <- labelSamples(a2,"title","A_","ctrl_")
ta2 <- analyseSeparate(a2, la2, "BY")
annotateAndOutput(a2, ta2, "diff_a2.tab")

a3 <- readGEO("GSE45251","GPL4133")
la3 <- labelSamples(a3,"title","TNFa","control")
ta3 <- analyseSeparate(a3, la3, "BY")
annotateAndOutput(a3, ta3, "diff_a3.tab")

a4 <- readGEO("GSE470","GPL8300")
la4 <- labelSamples(a4,"description","Astmatic","Normal")
ta4 <- analyseSeparate(a4, la4, "BY")
annotateAndOutput(a4, ta4, "diff_a4.tab")

a5 <- readGEO("GSE8190","GPL1708")
la5 <- labelSamples(a5,"title","[34]:","[12]:")
ta5 <- analyseSeparate(a5, la5, "BY")
annotateAndOutput(a5, ta5, "diff_a5.tab")

# Bacterial & Viral Infection
i1 <- readGEO("GSE40396","GPL10558")
li1 <- labelSamples(i1,"characteristics_ch1.5","group: [ABEHR]","Control")
ti1 <- analyseSeparate(i1, li1, "BY")
annotateAndOutput(i1, ti1, "diff_i1.tab")

i2 <- readGEO("GSE42026","GPL6947")
li2 <- labelSamples(i2,"title","WB-[HRb]","control")
ti2 <- analyseSeparate(i2, li2, "BY")
annotateAndOutput(i2, ti2, "diff_i2.tab")

i3 <- readGEO("GSE47172","GPL96")
li3 <- labelSamples(i3,"title","case","control")
ti3 <- analyseSeparate(i3, li3, "BY")
annotateAndOutput(i3, ti3, "diff_i3.tab")

i4 <- readGEO("GSE34205","GPL570")
li4 <- labelSamples(i4,"description","infection","healthy")
ti4 <- analyseSeparate(i4, li4, "BY")
annotateAndOutput(i4, ti4, "diff_i4.tab")

## Chronic Kidney Disease
ckd1 <- readGEO("GSE43484","GPL571")
lckd1 <- labelSamples(ckd1,"title","Uremic","Healthy")
tckd1 <- analyseSeparate(ckd1, lckd1, "BY")
annotateAndOutput(ckd1, tckd1, "diff_ckd1.tab")

ckd2 <- readGEO("GSE41030","GPL15950")
lckd2 <- labelSamples(ckd2,"title","uremic","control")
tckd2 <- analyseSeparate(ckd2, lckd2, "BY")
tckd2 <- subset(tckd2, select=c("ID","GENE_SYMBOL","adj.P.Val","P.Value","t","B","logFC","GENE_NAME"))
tckd2 <- tckd2[tckd2[,2] != "",];
write.table(tckd2, file="diff_ckd2.tab", row.names=F, sep="\t")

ckd3 <- readGEO("GSE38117","GPL4134")
lckd3 <- labelSamples(ckd3,"title","left","right")
tckd3 <- analyseSeparate(ckd3, lckd3, "BY")
annotateAndOutput(ckd3, tckd3, "diff_ckd3.tab")

ckd4 <- readGEO("GSE48041","GPL10787")
lckd4 <- labelSamples(ckd4,"description","hydronephrotic","control")
tckd4 <- analyseSeparate(ckd4, lckd4, "BY")
tckd4 <- subset(tckd4, select=c("ID","GENE_SYMBOL","adj.P.Val","P.Value","t","B","logFC","GENE_NAME"))
tckd4 <- tckd4[tckd4[,2] != "",];
write.table(tckd4, file="diff_ckd4.tab", row.names=F, sep="\t")

ckd5 <- readGEO("GSE15072","GPL96")
lckd5 <- labelSamples(ckd5,"title","Patient","Healthy")
tckd5 <- analyseSeparate(ckd5, lckd5, "BY")
annotateAndOutput(ckd5, tckd5, "diff_ckd5.tab")

## Cerebral Palsy
cp1 <- readGEO("GSE16447","GPL570")
lcp1 <- as.factor(c("G0","G1","G0","G0","G0","G0","G1","G1","G1"))
tcp1 <- analyseSeparate(cp1, lcp1, "BY")
annotateAndOutput(cp1, tcp1, "diff_cp1.tab")

cp2 <- readGEO("GSE31243","GPL571")
lcp2 <- labelSamples(cp2,"title","Palsy","Control")
tcp2 <- analyseSeparate(cp2, lcp2, "BY")
annotateAndOutput(cp2, tcp2, "diff_cp2.tab")

## Dilated Cardiomyopathy
dc1 <- readGEO("GSE29819","GPL570")
ldc1 <- labelSamples(dc1,"description"," failing "," control ")
tdc1 <- analyseSeparate(dc1, ldc1, "BY")
annotateAndOutput(dc1, tdc1, "diff_dc1.tab")

dc2 <- readGEO("GSE42955","GPL6244")
ldc2 <- labelSamples(dc2,"title","cardiomyopathy","Normal")
tdc2 <- analyseSeparate(dc2, ldc2, "BY")
annotateAndOutput(dc2, tdc2, "diff_dc2.tab")

## Ear Infection / Otitis Media
ei <- readGEO("GSE49128","GPL1261")
lei <- labelSamples(ei,"title","Hflu","control")
tei <- analyseSeparate(ei, lei, "BY")
annotateAndOutput(ei, tei, "diff_ei.tab")

## Inflammatory Bowel Disease
ibd1 <- readGEO("GSE11223","GPL1708")
libd1 <- labelSamples(ibd1,"title","UC","Normal")
tibd1 <- analyseSeparate(ibd1, libd1, "BY")
annotateAndOutput(ibd1, tibd1, "diff_ibd1.tab")

ibd2 <- readGEO("GSE3365","GPL96")
libd2 <- labelSamples(ibd2,"title","[UC]","Normal")
tibd2 <- analyseSeparate(ibd2, libd2, "BY")
annotateAndOutput(ibd2, tibd2, "diff_ibd2.tab")

ibd3 <- readGEO("GSE38713","GPL570")
libd3 <- labelSamples(ibd3,"description","UC","control")
tibd3 <- analyseSeparate(ibd3, libd3, "BY")
annotateAndOutput(ibd3, tibd3, "diff_ibd3.tab")

ibd4 <- readGEO("GSE9452","GPL570")
libd4 <- labelSamples(ibd4,"characteristics_ch1","patient","Control")
tibd4 <- analyseSeparate(ibd4, libd4, "BY")
annotateAndOutput(ibd4, tibd4, "diff_ibd4.tab")

## Upper Respiratory Infection
uri1 <- readGEO("GSE24132","GPL570")
luri1 <- labelSamples(uri1,"title","RSV","mock")
turi1 <- analyseSeparate(uri1, luri1, "BY")
annotateAndOutput(uri1, turi1, "diff_uri1.tab")
 
uri2 <- readGEO("GSE35940","GPL7202")
luri2 <- labelSamples(uri2,"title","lung_A","mock")
turi2 <- analyseSeparate(uri2, luri2, "BY")
annotateAndOutput(uri2, turi2, "diff_uri2.tab")

## Autism Spectrum Disorder
asd1 <- readGEO("GSE25507","GPL570")
lasd1 <- labelSamples(asd1,"characteristics_ch1.1","autism","control")
tasd1 <- analyseSeparate(asd1, lasd1, "BY")
annotateAndOutput(asd1, tasd1, "diff_asd1.tab")

asd2 <- readGEO("GSE7329","GPL1708")
lasd2 <- labelSamples(asd2,"title","autism","control")
tasd2 <- analyseSeparate(asd2, lasd2, "BY")
annotateAndOutput(asd2, tasd2, "diff_asd2.tab")

asd3 <- readGEO("GSE28521","GPL6883")
lasd3 <- labelSamples(asd3,"characteristics_ch1","autism","controls")
tasd3 <- analyseSeparate(asd3, lasd3, "BY")
annotateAndOutput(asd3, tasd3, "diff_asd3.tab")

asd4 <- readGEO("GSE26415","GPL6480")
lasd4 <- labelSamples(asd4,"characteristics_ch1","Autism","nonaustistic")
tasd4 <- analyseSeparate(asd4, lasd4, "BY")
annotateAndOutput(asd4, tasd4, "diff_asd4.tab")

asd5 <- readGEO("GSE6575","GPL570")
pasd5 <- pData(asd5)
filInd <- c(grep("Autism",t(pasd5["characteristics_ch1"])),grep("General",t(pasd5["characteristics_ch1"]))) # filter Autism and control samples
fasd5 <- asd5
exprs(fasd5) <- exprs(asd5)[,filInd]
pData(fasd5) <- pData(asd5)[t(filInd),]
lasd5 <- labelSamples(fasd5,"characteristics_ch1","Autism","General")
tasd5 <- analyseSeparate(fasd5, lasd5, "BY")
annotateAndOutput(fasd5, tasd5, "diff_asd5.tab")

asd6 <- readGEO("GSE18123","GPL570")
lasd6 <- labelSamples(asd6,"characteristics_ch1","[AP]","CONTROL")
tasd6 <- analyseSeparate(asd6, lasd6, "BY")
annotateAndOutput(asd6, tasd6, "diff_asd6.tab")

asd7 <- readGEO("GSE18123","GPL6244")
lasd7 <- labelSamples(asd7,"characteristics_ch1","[AP]","CONTROL")
tasd7 <- analyseSeparate(asd7, lasd7, "BY")
annotateAndOutput(asd7, tasd7, "diff_asd7.tab")

## Epilepsy
ep1 <- readGEO("GSE32534","GPL570")
lep1 <- labelSamples(ep1,"characteristics_ch1.2",": epilepsy",": no")
tep1 <- analyseSeparate(ep1, lep1, "BY")
annotateAndOutput(ep1, tep1, "diff_ep1.tab")

ep2 <- readGEO("GSE6834","GPL4757")
pep2 <- pData(ep2)
filInd <- grep("mTLE",t(pep2["title"])) # filter epileptic and control samples
fep2 <- ep2
exprs(fep2) <- exprs(ep2)[,filInd]
pData(fep2) <- pData(ep2)[t(filInd),]
lep2 <- labelSamples(fep2,"title","neocortex mTLE","Control")
tep2 <- analyseSeparate(fep2, lep2, "BY")
annotateAndOutput(fep2, tep2, "diff_ep2.tab")

ep3 <- readGEO("GSE6614","GPL81")
lep3 <- labelSamples(ep3,"characteristics_ch1.3","seizures","untreated")
tep3 <- analyseSeparate(ep3, lep3, "BY")
annotateAndOutput(ep3, tep3, "diff_ep3.tab")

ep4 <- readGEO("GSE47516","GPL1261")
lep4 <- labelSamples(ep4,"title","knockout","wild")
tep4 <- analyseSeparate(ep4, lep4, "BY")
annotateAndOutput(ep4, tep4, "diff_ep4.tab")

ep5 <- readGEO("GSE20977","GPL6947")
lep5 <- labelSamples(ep5,"title","patient","Control")
tep5 <- analyseSeparate(ep5, lep5, "BY")
annotateAndOutput(ep5, tep5, "diff_ep5.tab")

ep6 <- readGEO("GSE22225","GPL570")
lep6 <- labelSamples(ep6,"title","Patient","Control")
tep6 <- analyseSeparate(ep6, lep6, "BY")
annotateAndOutput(ep6, tep6, "diff_ep6.tab")

ep7 <- readGEO("GSE16969","GPL570")
lep7 <- labelSamples(ep7,"characteristics_ch1","surgical","autopsy")
tep7 <- analyseSeparate(ep7, lep7, "BY")
annotateAndOutput(ep7, tep7, "diff_ep7.tab")

## Schizophrenia
s1 <- readGEO("GSE53987","GPL570")
ps1 <- pData(s1)
filInd <- c(grep("schiz_",t(ps1["title"])),grep("control_",t(ps1["title"]))) # filter schiz and control samples
fs1 <- s1
exprs(fs1) <- exprs(s1)[,filInd]
pData(fs1) <- pData(s1)[t(filInd),]
ls1 <- labelSamples(fs1,"title","schiz","control")
ts1 <- analyseSeparate(fs1, ls1, "BY")
annotateAndOutput(fs1, ts1, "diff_s1.tab")

s2 <- readGEO("GSE27383","GPL570")
ls2 <- labelSamples(s2,"title","patient","control")
ts2 <- analyseSeparate(s2, ls2, "BY")
annotateAndOutput(s2, ts2, "diff_s2.tab")

s3 <- readGEO("GSE48072","GPL10558")
ls3 <- labelSamples(s3,"title","case","control")
ts3 <- analyseSeparate(s3, ls3, "BY")
annotateAndOutput(s3, ts3, "diff_s3.tab")

s4 <- readGEO("GSE46509","GPL1352")
ls4 <- labelSamples(s4,"title","Schizophrenia","Control")
ts4 <- analyseSeparate(s4, ls4, "BY")
annotateAndOutput(s4, ts4, "diff_s4.tab")

s5 <- readGEO("GSE37981","GPL1352")
ls5 <- labelSamples(s5,"title","Schizophrenia","Control")
ts5 <- analyseSeparate(s5, ls5, "BY")
annotateAndOutput(s5, ts5, "diff_s5.tab")

s6 <- readGEO("GSE21935","GPL570")
ls6 <- labelSamples(s6,"title","Scz","Control")
ts6 <- analyseSeparate(s6, ls6, "BY")
annotateAndOutput(s6, ts6, "diff_s6.tab")

s7 <- readGEO("GSE25673","GPL6244")
ls7 <- labelSamples(s7,"title","Patient","Control")
ts7 <- analyseSeparate(s7, ls7, "BY")
annotateAndOutput(s7, ts7, "diff_s7.tab")

s8 <- readGEO("GSE21138","GPL570")
ls8 <- labelSamples(s8,"title","Scz","Control")
ts8 <- analyseSeparate(s8, ls8, "BY")
annotateAndOutput(s8, ts8, "diff_s8.tab")

s9 <- readGEO("GSE17612","GPL570")
ls9 <- labelSamples(s9,"title","Scz","Control")
ts9 <- analyseSeparate(s9, ls9, "BY")
annotateAndOutput(s9, ts9, "diff_s9.tab")

s10 <- readGEO("GSE12654","GPL3345")
ps10 <- pData(s10)
filInd <- c(grep("schizophrenia",t(ps10["title"])),grep("control",t(ps10["title"]))) # filter schiz and control samples
fs10 <- s10
exprs(fs10) <- exprs(s10)[,filInd]
pData(fs10) <- pData(s10)[t(filInd),]
ls10 <- labelSamples(fs10,"title","schizophrenia","control")
ts10 <- analyseSeparate(fs10, ls10, "BY")
annotateAndOutput(fs10, ts10, "diff_s10.tab")

s11 <- readGEO("GSE12649","GPL96")
ps11 <- pData(s11)
filInd <- c(grep("schizophrenia",t(ps11["title"])),grep("control",t(ps11["title"]))) # filter schiz and control samples
fs11 <- s11
exprs(fs11) <- exprs(s11)[,filInd]
pData(fs11) <- pData(s11)[t(filInd),]
ls11 <- labelSamples(fs11,"title","schizophrenia","control")
ts11 <- analyseSeparate(fs11, ls11, "BY")
annotateAndOutput(fs11, ts11, "diff_s11.tab")

## Muscular Distrophy
md1 <- readGEO("GSE42806","GPL570")
lmd1 <- labelSamples(md1,"title","TMD","Control")
tmd1 <- analyseSeparate(md1, lmd1, "BY")
annotateAndOutput(md1, tmd1, "diff_md1.tab")

md2 <- readGEO("GSE36398","GPL6244")
lmd2 <- labelSamples(md2,"title","FSHD","control")
tmd2 <- analyseSeparate(md2, lmd2, "BY")
annotateAndOutput(md2, tmd2, "diff_md2.tab")

md3 <- readGEO("GSE9397","GPL96")
lmd3 <- labelSamples(md3,"title","FSHD","NHM")
tmd3 <- analyseSeparate(md3, lmd3, "BY")
annotateAndOutput(md3, tmd3, "diff_md3.tab")

md4 <- readGEO("GSE6011","GPL96")
lmd4 <- labelSamples(md4,"title","DMD","normal")
tmd4 <- analyseSeparate(md4, lmd4, "BY")
annotateAndOutput(md4, tmd4, "diff_md4.tab")