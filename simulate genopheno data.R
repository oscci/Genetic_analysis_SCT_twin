#simulate genotype-phenotype data
#started 6th Jan 2018

#simulate data for analysis with 20 SNPs with some linkage disequibrium, and 3 correlated phenotyes
#simulate 3 situations
# 1. no relationship between geno and pheno
# 2. strong relationship (r = .5) for 2 SNPs with pheno
# 3. weak relationship (r = .1) for 10 SNPs with pheno

#Aim: start by using GSCA approach and see if it captures the associations

require(doBy)
require(tidyverse)
require(MASS)
require(stats)

options(scipen = 999) #turn off scientific notation

nsnp<-24 #N SNPs to simulate
snpnames<-paste0('snp', 1:nsnp)
maf<-runif(nsnp,.25,.5) #minor allele freq set to value from .25 to .5
ncases<-100000
mypsummary<-data.frame(matrix(c(0,0,0),nrow=1)) #to hold pvalues for each run
mydata<-data.frame(matrix(nrow=ncases,ncol=(3+sum(nsnp))))
#We will simulate large datasets (ncases) for 3 phenotypes and 24 SNPs 
#We can then sample from the saved data when running simulations.

#In each dataset we have one phenotype with 3 correlated indicators (last 3 cols)
#Dataset name indicates whether any corelation between SNPs, and whether any phenogeno correl
#DatasetRand_N: No correlation between SNPs, no assoc with phenotype
#DatasetUncorr_2SNP_5: No correlation between SNPs, 2 SNPs correl with phenotype = .5
#Dataset3Block_N: SNPs in blocks of 8 with correlation within block dependent on distance
colnames(mypsummary) <- c('DatasetRand_N','DatasetUncorr_2SNP_5','Dataset3Block_N','Dataset3Block_2SNP_5')
thisfile<-3 #Indicates which type of correlation structure in datafile
mydatafile<- colnames(mypsummary)[thisfile]

#SNP values are 0, 1 or 2.

#-----------------------------------------------------------------
#simulate the SNPS as correlated zscores, correlation is mycorr
#-----------------------------------------------------------------
#We will assume correlation depends on how close SNPs are, ranging from .9 for adjacent and stepping 
#down by .1 for each additional distance (where sequence order is proxy for distance)

startcol<-1
endcol<-nsnp
h <-nsnp
mymean<-rep(0,h)
mycorr<-0 #default is uncorrelated
mycov2<-matrix(mycorr,nrow=h,ncol=h)
diag(mycov2)<-rep(1,h)
if(thisfile>2){ #only add correlation for conditions where thisfile>2
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
mydata[,startcol:endcol]=mvrnorm(n = ncases, mymean, mycov2)
colnames(mydata)[1:nsnp]<-snpnames

#-----------------------------------------------------------------
#Simulate the 3 phenotypes as correlated zscores, correlation is mycorr

#If any SNPs correlated with pheno, they will be simulated here as well
n2sim<-3
nsnpeffect<-0 #N SNPs with true effect on phenotype

if(thisfile==2 || thisfile==4){
  nsnpeffect<-2 #you need to vary this number to specify N SNPs with effect
  n2sim<-n2sim+nsnpeffect
  gpcov<-.5} #covariance for geno-pheno - this is under user control
mymean<-rep(0,n2sim)
mycorr<-.75
mycov<-matrix(mycorr,nrow=n2sim,ncol=n2sim)
if (nsnpeffect>0){
#adjust covs for the SNP columns
#set correlation between SNPs as specified previously
mycov[1:nsnpeffect,1:nsnpeffect]<-mycov2[1:nsnpeffect,1:nsnpeffect]
mycov[1:nsnpeffect,(nsnpeffect+1):n2sim] <- gpcov #set cov value to gpcov for the SNPs in the matrix
mycov[(nsnpeffect+1):n2sim,1:nsnpeffect] <- gpcov
}
#then set diagonal values to 1 for all
diag(mycov)<-rep(1,n2sim)
mydata[,(nsnp-nsnpeffect+1):(nsnp+3)]=mvrnorm(n = ncases, mymean, mycov)

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

#----------------------------------------------------------------------
# For each of nrun runs take a sample and analyse
#----------------------------------------------------------------------
nrun <- 100 #N runs to simulate
mybigdata<-read.table(thisdatafile, header=TRUE,sep=",",stringsAsFactors=FALSE) 
myfitsummary<-data.frame(matrix(,nrow=nrun,ncol=nsnp+5))
colnames(myfitsummary)[1:2]<-c('Path_GenoPheno','pvalue_GenoPheno')
colnames(myfitsummary)[3:(nsnp+5)]<-colnames(mydata)


nsub<-120
for (myn in 1:nrun){
  if(myn%%20==0){print(myn)}
  #read in data for required sample size
  myrows<-sample(ncases,nsub)
  mysample<-mybigdata[myrows,]
#use ASCSCA to analyse

# library(ASGSCA)
ObservedVar=colnames(mysample)[1:(nsnp+3)]
LatentVar=c("CNTNAP2","Neurodev")

#W0 is I x L matrix where rows are genotypes and traits, columns are genes and 'clinical pathways'
W0=matrix(rep(0,2*(3+nsnp)),nrow=3+nsnp,ncol=2, dimnames=list(ObservedVar,LatentVar))
W0[1:nsnp,1]=W0[(nsnp+1):(nsnp+3),2]=1

#B0 is L x L matrix with 0s and 1s, where 1 indicates arrow from latent variable in row to latent variable in column
B0=matrix(rep(0,2*2),nrow=2,ncol=2, dimnames=list(LatentVar,LatentVar))
B0[1,2]=1

# GSCA(mysample,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=FALSE,path=NULL,nperm=1000)
mynperm=100 #really need more than this but to test prog using 100 for speed
myfit<-GSCA(mysample,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=TRUE,path=NULL,nperm=mynperm)
myfitsummary[myn,1]<-myfit$Path[1,2]
myfitsummary[myn,2]<-myfit$pvalues[1,2]
myfitsummary[myn,3:(nsnp+2)]<-myfit$Weight[1:nsnp,1]
myfitsummary[myn,(nsnp+3):(nsnp+5)]<-myfit$Weight[(nsnp+1):(nsnp+3),2]
myfitsummary[myn,2]<-myfit$pvalues[1,2]
}

write.table(myfitsummary, paste0(mydatafile,'_results.csv'), sep=",",row.names=FALSE) 
