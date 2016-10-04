## Code to run hypergeometric enrichment analysis on individual microarray study selected under each disease
## Written by Sumaiya Nazeen, nazeen@mit.edu
## R version 2.15.1
## Platform: 3.13.0-35-generic #62-Ubuntu (64-bit)

## This code will allow one to reproduce the hypergeometric enrichment analyses in Nazeen et al., 2015

## Define Utility Function
penrich <- function(inFile, outFile, tot_genes, d_genes){
        x = read.table(file=inFile, header=FALSE)
		colnames(x) = c("name","ov","tot","geneList")
		p = apply(x,1,function(x) g = phyper(strtoi(x["ov"])-1,strtoi(x["tot"]),tot_genes-strtoi(x["tot"]),d_genes,lower.tail=FALSE))
		ad_p = p.adjust(p,"bonferroni")
		z = cbind(x,p,-log(p), ad_p)
		write.table(z, file=outFile, sep='\t', row.names=FALSE, col.names=TRUE)
}

## Make directories named {"all", "cp", "kegg", "bio", "reactome", "pid"} before running the following commands in R

## All C2 Gene Sets
penrich("eall_asd.tab","all/phyp_asd.tab",25222,1258)
penrich("eall_asthma.tab","all/phyp_asthma.tab",23047,858)
penrich("eall_inf.tab","all/phyp_inf.tab",24334,3630)
penrich("eall_ckd.tab","all/phyp_ckd.tab",48968,416)
penrich("eall_cp.tab","all/phyp_cp.tab",21246,220)
penrich("eall_dc.tab","all/phyp_dc.tab",24529,349)
penrich("eall_ei.tab","all/phyp_ei.tab",39636,7206)
penrich("eall_ibd.tab","all/phyp_ibd.tab",22291,2547)
penrich("eall_md.tab","all/phyp_md.tab",24529,517)
penrich("eall_s.tab","all/phyp_s.tab",25500,149)
penrich("eall_uri.tab","all/phyp_uri.tab",42544,59)
penrich("eall_ep.tab","all/phyp_ep.tab",45189,4)

## All Canonical Pathway Gene Sets
penrich("cp_asd.tab","cp/phyp_asd.tab",25222,1258) 
penrich("cp_asthma.tab","cp/phyp_asthma.tab",23047,858)
penrich("cp_inf.tab","cp/phyp_inf.tab",24334,3630)
penrich("cp_cd.tab","cp/phyp_cd.tab",96,14)
penrich("cp_ckd.tab","cp/phyp_ckd.tab",48968,416)
penrich("cp_cp.tab","cp/phyp_cp.tab",21246,220)
penrich("cp_dc.tab","cp/phyp_dc.tab",24529,349)
penrich("cp_ei.tab","cp/phyp_ei.tab",39636,7206)
penrich("cp_ibd.tab","cp/phyp_ibd.tab",22291,2547)
penrich("cp_md.tab","cp/phyp_md.tab",24529,517)
penrich("cp_s.tab","cp/phyp_s.tab",25500,149)
penrich("cp_uri.tab","cp/phyp_uri.tab",42544,59)
penrich("cp_ep.tab","cp/phyp_ep.tab",45189,4)

## KEGG Pathway Gene Sets
penrich("k_asd.tab","kegg/phyp_asd.tab",25222,1258)   
penrich("k_asthma.tab","kegg/phyp_asthma.tab",23047,852)   
penrich("k_inf.tab","kegg/phyp_inf.tab",24334,3630)   
penrich("k_cd.tab","kegg/phyp_cd.tab",96,14)   
penrich("k_ckd.tab","kegg/phyp_ckd.tab",48968,416)   
penrich("k_cp.tab","kegg/phyp_cp.tab",21246,220)   
penrich("k_dc.tab","kegg/phyp_dc.tab",24529,349)   
penrich("k_ei.tab","kegg/phyp_ei.tab",39636,7206)     
penrich("k_ibd.tab","kegg/phyp_ibd.tab",22291,2547)   
penrich("k_md.tab","kegg/phyp_md.tab",24529,517)      
penrich("k_s.tab","kegg/phyp_s.tab",25500,149)   
penrich("k_uri.tab","kegg/phyp_uri.tab",42544,59)   
penrich("k_ep.tab","kegg/phyp_ep.tab",45189,4)

## Biocarta Gene Sets
penrich("b_asd.tab","bio/phyp_asd.tab",25222,1258)
penrich("b_asthma.tab","bio/phyp_asthma.tab",23047,858)
penrich("b_inf.tab","bio/phyp_inf.tab",24334,3630)
penrich("b_cd.tab","bio/phyp_cd.tab",96,14)
penrich("b_ckd.tab","bio/phyp_ckd.tab",48968,416)
penrich("b_cp.tab","bio/phyp_cp.tab",21246,220)
penrich("b_dc.tab","bio/phyp_dc.tab",24529,349)
penrich("b_ei.tab","bio/phyp_ei.tab",39636,7206)
penrich("b_ibd.tab","bio/phyp_ibd.tab",22291,2547)
penrich("b_md.tab","bio/phyp_md.tab",24529,517)
penrich("b_s.tab","bio/phyp_s.tab",25500,149)
penrich("b_uri.tab","bio/phyp_uri.tab",42544,59)
penrich("b_ep.tab","bio/phyp_ep.tab",45189,4)

## Reactome Gene Sets
penrich("r_asd.tab","reactome/phyp_asd.tab",25222,1258)
penrich("r_asthma.tab","reactome/phyp_asthma.tab",23047,858)
penrich("r_inf.tab","reactome/phyp_inf.tab",24334,3630)
penrich("r_cd.tab","reactome/phyp_cd.tab",96,14)
penrich("r_ckd.tab","reactome/phyp_ckd.tab",48968,416)
penrich("r_cp.tab","reactome/phyp_cp.tab",21246,220)
penrich("r_dc.tab","reactome/phyp_dc.tab",24529,349)
penrich("r_ei.tab","reactome/phyp_ei.tab",39636,7206)
penrich("r_ibd.tab","reactome/phyp_ibd.tab",22291,2547)
penrich("r_md.tab","reactome/phyp_md.tab",24529,517)
penrich("r_s.tab","reactome/phyp_s.tab",25500,149)
penrich("r_uri.tab","reactome/phyp_uri.tab",42544,59)
penrich("r_ep.tab","reactome/phyp_ep.tab",45189,4)

## PID Gene Sets    
penrich("pid_asd.tab","pid/phyp_asd.tab",25222,1258)
penrich("pid_asthma.tab","pid/phyp_asthma.tab",23047,858)
penrich("pid_inf.tab","pid/phyp_inf.tab",24334,3630)
penrich("pid_cd.tab","pid/phyp_cd.tab",96,14)
penrich("pid_ckd.tab","pid/phyp_ckd.tab",48968,416)
penrich("pid_cp.tab","pid/phyp_cp.tab",21246,220)
penrich("pid_dc.tab","pid/phyp_dc.tab",24529,349)
penrich("pid_ei.tab","pid/phyp_ei.tab",39636,7206)
penrich("pid_ibd.tab","pid/phyp_ibd.tab",22291,2547)
penrich("pid_md.tab","pid/phyp_md.tab",24529,517)
penrich("pid_s.tab","pid/phyp_s.tab",25500,149)
penrich("pid_uri.tab","pid/phyp_uri.tab",42544,59) 
penrich("pid_ep.tab","pid/phyp_ep.tab",45189,4)

