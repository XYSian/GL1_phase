#!/usr/bin/env Rscript 
#parse parameter
options()$repos 

options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror

# # https://bioconductor.org/packages/release/bioc/html/GEOquery.html
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("KEGG.db",ask = F,update = F)
# BiocManager::install(c("GSEABase","GSVA","clusterProfiler" ),ask = F,update = F)
# BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
# BiocManager::install(c("org.Hs.eg.db","hgu133plus2.db" ),ask = F,update = F)
if(!requireNamespace("BioManager",quietly=TRUE))
  insatll.packages("BioManager")
BiocManager::install("subread")
BiocManager::install("Rsubread")
BiocManager::install("limma")
BiocManager::install("edgeR")
library(argparser,quietly = TRUE)
#Creat a parser
p<- arg_parser("run featureCounts and calculate FPKM/TPM")

#Add command line argumnets

p<-add_argument(p,"--bam",help="input: bam file",type="character")

p<-add_argument(p,"--gtf",help="input: gtf file",type="character")

p<-add_argument(p,"--output",help="out prefix",type="character")

#Parse the command line arguments

argv<-parse_args(p)
library(Rsubread)
library(limma)
library(edgeR)

bamFile<- argv$bam
gtfFile<- argv$gtf
nthreads<- 1
outFilePref<- argv$output

outStatsFilePath<- paste(outFilePref, '.log',sep='');
outCountsFilePath<- paste(outFilePref,'.count',sep='');

fCountsList=featureCounts(bamFile,annot.ext=gtfFile,isGTFAnnotationFile=TRUE,nthreads=nthreads,isPairedEnd=TRUE)
dgeList=DGEList(counts=fCountsList$counts,genes=fCountsList$annotation)
fpkm=rpkm(dgeList,dgeList$genes$Length)
tpm=exp(log(fpkm)-log(sum(fpkm))+log(1e6))

write.table(fCountsList$stat,outStatsFilePath,sep="\t",col.names = FALSE,row.names = FALSE,quote=FALSE)
featureCounts=cbind(fCountsList$annotation[,1],fCountsList$counts,fpkm,tpm)
colnames(featureCounts)=c('gene_id','counts','fpkm','tpm')

write.table(featureCounts,outCountsFilePath,sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
