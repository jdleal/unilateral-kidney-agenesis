# unilateral-kidney-agenesis
module load R
R
# load SAIGE
.libPaths(c(".../R/x86_64-pc-linux-gnu-library/3.6", .libPaths()))
library(SAIGEgds)
# load other packages
.libPaths(c(".../R/x86_64-pc-linux-gnu-library/3.6", .libPaths()))
library(SeqArray)
library(SNPRelate)
library(GWASTools)
setwd(".../missingKidney/")
# open the GDS file for genetic relationship matrix (GRM)
# save PLINK files in .../missingKidney/data/
grm_fn2 <- snpgdsBED2GDS("./data/round8_unpruned.bed", "./data/round8_unpruned.fam", "./data/round8_unpruned.bim", "./data/round8_unpruned_2.gds")
grm_gds2 <- SNPRelate::snpgdsOpen(grm_fn2)
#seqClose(grm_fn2) # won't work until it is a seqarray file 
closefn.gds(grm_gds2) # works for general gds files
gds.fn3 <- "./data/round8_unpruned_2.gds"
seqSNP2GDS(gds.fn3, "out.gds", storage.option="LZMA_RA", major.ref=TRUE,
           ds.type=c("packedreal16", "float", "double"), optimize=TRUE, digest=TRUE,
           verbose=TRUE)
seqSummary("out.gds", varname=NULL, check=c("default", "none", "full"),
           verbose=TRUE)
pheno <- read.table("./data/unilateralKidneyAgenesisPhenotype.txt", header=TRUE, as.is=TRUE)
head(pheno2)
# perform the association analysis, model 0
glmm0 <- seqFitNullGLMM_SPA(missing ~ 1, pheno, "out.gds", trait.type="binary", num.thread=15)
assoc0 <- seqAssocGLMM_SPA("out.gds", glmm0, mac=10)
head(assoc0)
write.csv(assoc0,'./output/assoc0.csv')
seqClose("out.gds") # won't work until it is a seqarray file 
