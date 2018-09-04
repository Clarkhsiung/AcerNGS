#QCtest
library(vcfR)
library('adegenet')
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
setwd("C:/Users/Clarke/Desktop/stacks results/spe70")	
getwd() 	

##Extract sample deapth
vcf <- read.vcfR("populations.snps.vcf")
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE) 

pdf("DP_RAD_data.pdf", width = 10, height=5) # boxplot
par(mar=c(8,4,1,1)) 
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)",
        las=2, cex=0.4, cex.axis=0.5)
dev.off()

#zoom to smaller values
pdf("DP_RAD_data_zoom.pdf", width = 10, height=5) # boxplot
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)",
        las=2, cex=0.4, cex.axis=0.5, ylim=c(0,50))
abline(h=8, col="red")
dev.off() 

# add real SNP.names
aa.genlight <- vcfR2genlight(vcf, n.cores=1)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3) 
aa.genlight 
indNames(aa.genlight) # check individual names 
as.matrix(aa.genlight)[1:16,1:10] # see tiny bit of the data 
glPlot (aa.genlight) # takes some time 
x <- summary(t(as.matrix(aa.genlight)))
write.table(x[7,], file = "missing.persample.txt", sep = "\t")

##Caculate missing rate
vcf=read.vcfR("populations.snps.vcf")
gq <- extract.gt(vcf, element="GQ", as.numeric=TRUE)
dp <- extract.gt(vcf, element="DP", as.numeric=TRUE)
myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- myMiss/nrow(vcf)
miss=as.data.frame(myMiss)
miss$sample=rownames(miss)
miss$pop=substr(miss$sample,1,3)
miss$spe=substr(miss$pop,1,2)
miss$loc=substr(miss$pop,3,3)
ar.miss=subset(miss,spe=="AR")
ak.miss=subset(miss,spe=="AK")
library(ggplot2)
colo=c("brown1","aquamarine2","khaki1","lightslateblue","gray38")
theme_set(theme_classic())
g=ggplot(miss,aes(x=loc,y=myMiss))+facet_grid(spe~.)+geom_hline(yintercept = 0.7,color="red")
g+geom_boxplot()+ylab("Missing rate")+xlab("Sample location")+ggtitle("Miss rate per population")
g+geom_point()
null=dev.off()


#Creat a list to remove sample
remov=subset(miss,myMiss>=0.8)
remove.sample=remov$sample
remove.sample=as.vector(remove.sample)
genind=vcfR2genind(vcf)
#Remove sample
for(i in remove.sample){
  genind=genind[indNames(genind)!= i]
}
pop(genind)=substr(indNames(genind),1,2)
genpop=genind2genpop(genind)
genpop


##plot missing heatmap
library(grur)
library(radiator)
strata=individuals2strata(
  data = "c:\\Users\\Clarke\\Desktop\\pop.txt",1,2,
  filename = "strata.tsv"
)
imiss=missing_visualization(
  data="populations.snps.vcf",
  strata = "strata.tsv"
)
imiss$heatmap
imiss$ibm.plots
ggplot2::ggsave(
  filename = "heatmap.missing_spe.png",
  plot = imiss$heatmap)

ggplot2::ggsave(
  filename = "pct.missing.plot.strata.POP_ID.pdf",
  plot =imiss$pct.missing.total$pct.missing.plot.strata.POP_ID)
