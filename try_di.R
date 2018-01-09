#Try GCSA approach with Di's data from SLIC

require(MASS)
require(stats)

library(ASGSCA)

options(scipen = 999) #turn off scientific notation
didat<-read.table('SLICmergedPROcounttest.txt',header=TRUE,stringsAsFactors=FALSE)
#Gene 1 - rs10255943 to rs2396765 (13 SNPs) foxp2
#Gene 2 - rs10271363 to rs17170742, 23 SNPs cntnap2
#Gene 3 - rs16957277 to rs16957385, 7 SNPs (edited)
#Gene 4 - rs11149652 to rs4782962, 18 SNPs atp2c2
didat<-didat[,3:length(didat)] #strip off first 2 columns
didat[is.na(didat)]<-1 #just to get this to run I have replaced NA with 1!

#Now try model with all 4 genes together
ObservedVar=c(colnames(didat))
LatentVar=c("Gene1","Gene2","Gene3","Gene4","Language")
#W0 is I x L matrix where rows are genotypes and traits, columns are genes and latent phenotype
W0=matrix(rep(0,5*length(ObservedVar)),nrow=length(ObservedVar),ncol=5, dimnames=list(ObservedVar,LatentVar))
W0[1:13,1]=1
W0[14:36,2]=1
W0[37:43,3]=1
W0[44:61,4]=1
W0[62:64,5]=1

#B0 is L x L matrix with 0s and 1s, where 1 indicates arrow from latent geno in row to latent phenotype in column
B0=matrix(rep(0,5*5),nrow=5,ncol=5, dimnames=list(LatentVar,LatentVar))
B0[1:4,5]=1
mynperm=100 #probably need more than this but to test prog use 100 for speed
myfit<-GSCA(didat,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=TRUE,path=NULL,nperm=mynperm)


