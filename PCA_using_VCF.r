library(vcfR)
library('adegenet')
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap) 
setwd("C:/Users/Clarke/Desktop/stacks results/spe80")	
getwd() 
vcf <- read.vcfR("populations.snps.vcf")

# add real SNP.names
aa.genlight <- vcfR2genlight(vcf, n.cores=1)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")
spemap=read.table("spemap.tsv",sep = "\t")
pop(aa.genlight)<-spemap[,2]
aa.genlight 
indNames(aa.genlight) # check individual names 
as.matrix(aa.genlight)[1:16,1:10] # see tiny bit of the data 
 
###plot total AFS of the dataset
mySum <- glSum(aa.genlight, alleleAsUnit = TRUE)
barplot(table(mySum), col="blue", space=0, xlab="Allele counts",
        main="Distribution of ALT allele counts in total dataset") 
pca.1 <- glPca(aa.genlight, nf=300, n.cores=1) 
pdf ("PCA_all_SNPs_ax12.pdf", width=14, height=7)
col=c("brown1","aquamarine2","khaki1","lightslateblue","gray38")
g1 <- s.class(pca.1$scores, pop(aa.genlight), xax=1, yax=2,
              col=transp(col,.6),
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T,
              pgrid.draw =F, plot = FALSE)
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels =
                 list(box = list(draw = FALSE),
                      
                      optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off() 
