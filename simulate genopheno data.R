#simulate genotype-phenotype data
#started 6th Jan 2018

#simulate data for analysis with 20 SNPs with some linkage disequibrium, and 3 correlated phenotyes
#simulate 3 situations
# 1. no relationship between geno and pheno
# 2. strong relationship (r = .5) for small N SNPs with pheno
# 3. weak relationship (r = .1) for large N SNPs with pheno

#Aim: start by using GSCA approach and see if it captures the associations

require(doBy)
require(tidyverse)
require(MASS)
require(stats)

library(ASGSCA)

options(scipen = 999) #turn off scientific notation

#We will simulate large datasets for 3 correlated phenotypes and 24 SNPs 
#We can then sample from the saved data when running simulations.
nsnp<-24 #N SNPs to simulate - for correlated case divided into 3 blocks of 8, correlation within block
snpnames<-paste0('snp', 1:nsnp)
maf<-runif(nsnp,.25,.5) #minor allele freq set to value from .25 to .5
ncases<-100000
mydata<-data.frame(matrix(nrow=ncases,ncol=(3+sum(nsnp))))

#In each dataset we have one phenotype with 3 correlated indicators (last 3 cols)
#Dataset name indicates whether any corelation between SNPs, and whether any phenogeno correl
#e.g. DatasetRand_N: No correlation between SNPs, no assoc with phenotype
#DatasetUncorr_2SNP_5: No correlation between SNPs, 2 SNPs correl with phenotype = .5
#Dataset3Block_N: SNPs in blocks of 8 with correlation within block dependent on distance

#----------------------------------------------------------------------------
# User specifies how many SNPs have effect on phenotype and how big an effect here
# These values are ignored if thisfile is 1 or 3 (no effect)
nsnpeff<-3
gpcov<-.3
n2sim<-3 #N phenotypes

#----------------------------------------------------------------------------
myfilenames <- c('DatasetRand_N',
                 paste0('DatasetUncorr_',nsnpeff,'SNP_',10*gpcorr),
                 'Dataset3Block_N',
                 paste0('Dataset3Block_',nsnpeff,'SNP_',10*gpcorr))
thisfile<-4 #User specified: Indicates which type of correlation structure in datafile (could use loop here)
mydatafile<- myfilenames[thisfile]

#SNP values are 0, 1 or 2, but we start with random normal numbers

#-----------------------------------------------------------------
#simulate the SNPS as correlated zscores, correlation is mycorr
#-----------------------------------------------------------------
#Where there is correlation, assume correlation depends on how close SNPs are, 
#achieved by making size of correlation proportional to distance
# (where sequence order is proxy for distance)

h <-nsnp
mymean<-rep(0,h) #vector of meeans to use in mvrnorm function
mycorr<-0 #default is uncorrelated
mycov2<-matrix(mycorr,nrow=h,ncol=h) #cov matrix for mvrnorm
diag(mycov2)<-rep(1,h) #ones on diagonal of cov matrix
if(thisfile>2){ #only add correlation for conditions where thisfile is 3 or 4
  for (i in 1:h){
    irange<-1+as.integer((i-1)/8) #irange specifies haplotype block for i
    for (j in 1:h){
      jrange<- 1+as.integer((j-1)/8) #jrange specifies haplotype block for j
      if(irange==jrange){
       k<- abs(i-j)
      thisr<-1-8*k/100 #tweaked so magnitude of correlation declines with distance between SNPs
      if(thisr<0){thisr<-0}
      mycov2[i,j]<-thisr
      mycov2[j,i]<-mycov2[i,j] #correlation determined by distance between SNPs!
      if(k==0){mycov2[i,j]<-1}
      }
    }
  }
}
mydata[,1:nsnp]=mvrnorm(n = ncases, mymean, mycov2)
colnames(mydata)[1:nsnp]<-snpnames

#-----------------------------------------------------------------
#Simulate the 3 phenotypes as correlated zscores, correlation is mycorr

#If any SNPs correlated with pheno, they will be included here as well
#SNPs with effect always at the end of the snp columns

if(thisfile==2 || thisfile==4){
  n2sim<-n2sim+nsnpeff}
mymean<-rep(0,n2sim)
mycorr<-.75
mycov<-matrix(mycorr,nrow=n2sim,ncol=n2sim)
if (nsnpeff>0){
#adjust covs for the SNP columns
#set correlation between SNPs as specified in previous block of code
mycov[1:nsnpeff,1:nsnpeff]<-mycov2[1:nsnpeff,1:nsnpeff]
mycov[1:nsnpeff,(nsnpeff+1):n2sim] <- gpcov #set cov value to gpcov for the SNPs in the matrix
mycov[(nsnpeff+1):n2sim,1:nsnpeff] <- gpcov
}
#then set diagonal values to 1 for all
diag(mycov)<-rep(1,n2sim)
mydata[,(nsnp-nsnpeff+1):(nsnp+3)]=mvrnorm(n = ncases, mymean, mycov)

colnames(mydata)[(nsnp+1):(nsnp+3)]<-c('NwdRepPheno','LangPheno','NeurodevPheno')

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
thisdatafile<-paste0(mydatafile,'.csv')
write.table(mydata, thisdatafile, sep=",",row.names=FALSE) 
myr<-cor(mydata) #to check you have desired correlation structure, View(myr)
#----------------------------------------------------------------------
# For each of nrun runs take a sample and analyse
#----------------------------------------------------------------------
nrun <- 100 #N runs to simulate
mybigdata<-read.table(thisdatafile, header=TRUE,sep=",",stringsAsFactors=FALSE) 
myfitsummary<-data.frame(matrix(,nrow=nrun,ncol=nsnp+5)) #initialise dataframe to hold results
colnames(myfitsummary)[1:2]<-c('Path_GenoPheno','pvalue_GenoPheno')
colnames(myfitsummary)[3:(nsnp+5)]<-colnames(mydata) #columns for weights


nsub<-120 #N subjects - set to resemble our study
for (myn in 1:nrun){
  if(myn%%10==0){print(myn)} #show progress on screen
  #read in data for required sample size
  myrows<-sample(ncases,nsub) 
  mysample<-mybigdata[myrows,]
#use ASCSCA to analyse

ObservedVar=colnames(mysample)[1:(nsnp+3)]
LatentVar=c("CNTNAP2","Neurodev")

#W0 is I x L matrix where rows are genotypes and traits, columns are genes and latent phenotype
W0=matrix(rep(0,2*(3+nsnp)),nrow=3+nsnp,ncol=2, dimnames=list(ObservedVar,LatentVar))
W0[1:nsnp,1]=W0[(nsnp+1):(nsnp+3),2]=1

#B0 is L x L matrix with 0s and 1s, where 1 indicates arrow from latent geno in row to latent phenotype in column
B0=matrix(rep(0,2*2),nrow=2,ncol=2, dimnames=list(LatentVar,LatentVar))
B0[1,2]=1

# GSCA(mysample,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=FALSE,path=NULL,nperm=1000)
# for quick scrutiny of weights use this -but for pvalues need slow version using path.test

mynperm=100 #probably need more than this but to test prog use 100 for speed
myfit<-GSCA(mysample,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=TRUE,path=NULL,nperm=mynperm)
myfitsummary[myn,1]<-myfit$Path[1,2]
myfitsummary[myn,2]<-myfit$pvalues[1,2]
myfitsummary[myn,3:(nsnp+2)]<-myfit$Weight[1:nsnp,1]
myfitsummary[myn,(nsnp+3):(nsnp+5)]<-myfit$Weight[(nsnp+1):(nsnp+3),2]
myfitsummary[myn,2]<-myfit$pvalues[1,2]
}

write.table(myfitsummary, paste0(mydatafile,'_results.csv'), sep=",",row.names=FALSE) 
w<-which(myfitsummary$pvalue_GenoPheno<.05)
print(paste('N SNPs with true effect = ',nsnpeff,'; effect size = ',gpcov))
runtype<-'Uncorrelated SNPS'
if (thisfile>2){runtype<-'SNPs in 3 correlated blocks (see myr for correlations)'}
print(runtype)
print(paste0('Proportion of runs with p < .05 = ',length(w),' out of ',nrun))
