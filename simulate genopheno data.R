#simulate genotype-phenotype data
#started 6th Jan 2018

#simulate data for analysis with 20 SNPs with some linkage disequibrium, and 3 correlated phenotyes
#simulate 3 situations
# 1. no relationship between geno and pheno
# 2. strong relationship (r = .75) for 2 SNPs with pheno
# 3. weak relationship (r = .1) for 10 SNPs with pheno

#Aim: start by using GSCA approach and see if it captures the associations

require(doBy)
require(tidyverse)
require(MASS)
require(stats)

options(scipen = 999) #turn off scientific notation
nrun <- 1 #N runs to simulate
nsnp<-24 #N SNPs to simulate
snpnames<-paste0('snp', 1:nsnp)
maf<-runif(nsnp,.25,.5) #minor allele freq set to value from .25 to .5
ncases<-10000
mydata<-data.frame(matrix(nrow=ncases,ncol=(3+sum(nsnp))))
#We will simulate large datasets (ncases) for 3 phenotypes and 24 SNPs 
#We can then sample from the saved data when running simulations.

#In each dataset we have one phenotype with 3 correlated indicators (last 3 cols)
#Dataset name indicates whether any corelation between SNPs, and whether any phenogeno correl
#DatasetRand_N: No correlation between SNPs, no assoc with phenotype
#Dataset3Block_N: SNPs in blocks of 8 with correlation within block dependent on distance
#DatasetRand_2GP_5: No correlation between SNPs, 2 SNPs correl with phenotype = .5
datafilename <- c('DatasetRand_N.csv','Dataset3Block_N.csv','DatasetRand_2GP_5')
thisfile<-1 #Indicates which type of correlation structure in datafile
mydatafile<- datafilename[thisfile]

#SNP values are 0, 1 or 2.
colnames(mydata)[(nsnp+1):(nsnp+3)]<-c('NwdRepPheno','LangPheno','NeurodevPheno')

#-----------------------------------------------------------------
#Simulate the 3 phenotypes as correlated zscores, correlation is mycorr
mymean<-rep(0,3)
mycorr<-.75
mycov<-matrix(mycorr,nrow=3,ncol=3)
diag(mycov)<-rep(1,3)
mydata[,(nsnp+1):(nsnp+3)]=mvrnorm(n = ncases, mymean, mycov)
#-----------------------------------------------------------------
#Now for each gene, simulate the SNPS as correlated zscores, correlation is mycorr
#We will assume correlation depends on how close SNPs are, ranging from .9 for adjacent and stepping 
#down by .1 for each additional distance (where sequence order is proxy for distance)

startcol<-1
endcol<-nsnp
  h <-nsnp
  mymean<-rep(0,h)
  mycorr<-0 #default is uncorrelated
  mycov<-matrix(mycorr,nrow=h,ncol=h)
  diag(mycov)<-rep(1,h)
  if(thisfile>1){
  for (i in 1:h){
    for (j in 1:h){
      k<- abs(i-j)
      mycov[i,j]<-(1-5*k/100)
      mycov[j,i]<-mycov[i,j] #correlation determined by distance between SNPs!
      
      
      if(k==0){mycov[i,j]<-1}
    }
  }
  }
  mydata[,startcol:endcol]=mvrnorm(n = ncases, mymean, mycov)
  colnames(mydata)[startcol:endcol]<-snpnames

#-------------------------------------------------------------------------
#Convert gene scores to integers: 0-2 for autosomal

firstcol<-1
lastcol<-nsnp
p<-c(0,0,0) #initialise a vector to hold p values for different N alleles
for (i in 1:nsnp){
    p[1]<-(1-maf[i])^2
    p[2]<-2*(1-maf[i])*maf[i]
    p[3]<-maf[i]^2

  #now use p-values to assign z-score cutoffs that convert to 0,1,2 or 3 minor alleles
  temp<-mydata[,i]
  w<-which(temp<qnorm(p[1]))
  mydata[w,i]<-0
  w<-which(temp>qnorm(p[1]))
  mydata[w,i]<-1
  w<-which(temp>qnorm(p[2]+p[1]))
  mydata[w,i]<-2

}
write.table(mydata, "dummydata.txt", sep=",",row.names=FALSE) 

prun<-vector() #initialise empty vector to hold p-values
for (myn in 1:nrun){
  #need to read in data for required sample size
#use ASCSCA to analyse

# library(ASGSCA)
ObservedVar=colnames(mydata)[1:(nsnp+3)]
LatentVar=c("CNTNAP2","Neurodev")

#W0 is I x L matrix where rows are genotypes and traits, columns are genes and 'clinical pathways'
W0=matrix(rep(0,2*(3+nsnp)),nrow=3+nsnp,ncol=2, dimnames=list(ObservedVar,LatentVar))
W0[1:nsnp,1]=W0[(nsnp+1):(nsnp+3),2]=1

#B0 is L x L matrix with 0s and 1s, where 1 indicates arrow from latent variable in row to latent variable in column
B0=matrix(rep(0,2*2),nrow=2,ncol=2, dimnames=list(LatentVar,LatentVar))
B0[1,2]=1

 GSCA(mydata,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=FALSE,path=NULL,nperm=1000)

#myfit<-GSCA(mydata,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=TRUE,path=NULL,nperm=1000)

prun<-c(prun,myfit$pvalues[1,2])
}
#This version of program gave 6 sig pvalues (< .05) in 110 runs!