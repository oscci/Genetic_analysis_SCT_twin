# Checking impact of ascertainment bias on regression data

# Our SCT sample is mixture of prenatal and postnatal identified cases
# We know the latter will have more impairment
# Want to check the impact on power etc of having this mixture
# The problem is similar to the one about using cutoffs for groups that was discussed in SQING

# Program based on SQING Simulation of genotype-phenotype relations by DVM Bishop
# See script 'simulating_genopheno_cutoffs.R' on https://osf.io/nxspw
# Note, however, that in SQING we looked at sampling bias that led to restriction of range
# in phenotypes. Here we look at oversampling of impaired cases that does not affect the range of the combined sample.

# This version started 18th October 2017 by D. V. M.Bishop

# Conclusion from simulation: inclusion of selected cases does reduce power to detect 
# regression of phenotype on genotype, but the effect is small enough not to give concern.
#-------------------------------------------------------------------------
# From sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Sierra 10.12.6

library(MASS) #for mvrnorm function to make multivariate normal distributed vars
library(tidyverse)
options(scipen=999) #disable scientific notation.
# Use real data to determine % with ascertainment bias, and how big the impact
dir<-"/Users/dorothybishop/Dropbox/ERCadvanced/project SCT analysis/SCT_ASD_analysis/Project_files/Data/"
mydata <- read.csv(paste0(dir,"SCTData_DATA_2017-11-01_1815.csv"))
mydata<-filter(mydata,trisomy<9) #remove case of isochromosome
#deal with missing data
for (mycol in 68:228){
  mymiss=which(mydata[,mycol]>900)
  mydata[mymiss,mycol]=NA
}
require(psych)
mydata$subgp<-1
w<-c(which(mydata$why_tested==2),which(mydata$why_tested==3))
mydata$subgp[w]<-2
describeBy(mydata$sent_rep_ss,mydata$subgp)
#Shows that difference in mean for groups by asc bias is around .9
#(Consistent with Paul's analysis of language factor)
#-------------------------------------------------------------------------
# Specify parameters to create correlated variables
#-------------------------------------------------------------------------
#NB actual data 50% selected bcs behav/neuro problems; with effect size d = .9
nVar<-2 #number of simulated variables 
myM<-0 # Mean score for simulated variables
myVar<-1 #Variance for simulated variables
myN<-130 #set sample size per group (You can vary this to see the effect)

nSims<-5000 #arbitary N simulations
mycutoff<- -.25 #final group will be a mix of general population and below cutoff
#This value selected to give mean difference (Cohen's d) between selected/unselected cases around .9
#(Can check at end by comparing mean.selected and mean.unselected)
mylo.p<-c(.01,.25,.5,.99) #proportion of cases selected as below cutoff (we will loop through these values)
#1st value allows estimate of totally unselected; last value for totally selected

corrlist<-c(.29) #actual correlation between vars; corresponds to effect size for impact of genotype, d, of .5^2 
#NB; correl with genotype is lower than with z because genotype is noncontinuous
#value of .29 here gives correl with genotype of around .25
summarytable<-data.frame(matrix(rep(NA,12*length(mylo.p)*nSims),ncol=12)) #initialise table to hold results
colnames(summarytable)<-c('Nsub','truer','cutoff','proplo','p.r','p.chi','N_aa','N_Aa','N_AA',
                          'mean_aa','mean_Aa','mean_AA')
for (myCorr in corrlist){ #correlation between variables; can loop through various values
  #----------------------------------------------------------------------------------------
  # Generate a sample from a multivariate normal distribution with the specified correlation matrix
  #----------------------------------------------------------------------------------------
  
  myCov<-matrix(rep(myCorr,nVar*nVar),nrow=nVar) #rep(x,y) will generate y values of x
  diag(myCov)<-rep(myVar,nVar)  #variance on diagonal
  mydata<-data.frame(mvrnorm(n = myN*100, rep(myM,nVar), myCov)) #make big population to select from
  colnames(mydata)<-c('myz','pheno') #correlated geno and pheno, where myz is geno
  
  #---------------------------------------------------------------------------------------
  # convert the random normal deviates to genotype values of 0, 1 or 2, depending on MAF
  MAF<-.5 # minor allele frequency
  Naa<-MAF*MAF
  NAA<-(1-MAF)^2
  NaA<-1-Naa-NAA
  mydata$mygeno<-1 #initialise mygeno with genotype Aa
  critA<--qnorm(NAA)
  mydata[mydata$myz>critA,3]<-2 #cases with 2 copies of maj allele, AA
  crita<-qnorm(Naa)
  mydata[mydata$myz<crita,3]<-0  #cases with no copies of maj allele, aa
  #---------------------------------------------------------------------------------------
  #test significant of correlation of genotype/phenotype for full population: sanity check
  rtest<-cor.test(mydata$pheno,mydata$mygeno)
  myr<-rtest$estimate
  
  regtest<-summary(lm(pheno~mygeno,data=mydata))
  myp<-regtest$coefficients[8] #pvalue for linear regression of pheno on geno
  
  #---------------------------------------------------------------------------------------
  #Check if distribution of genotypes is in line with expectation from MAF
  myevalue<-c(as.integer(myN*NAA),as.integer(myN*NaA),as.integer(myN*Naa)) #expected distribution of genotypes
  myovalue<-table(mydata$mygeno[1:myN]) #observed distribution of genotypes
  mychi<-chisq.test(rbind(myevalue,myovalue))
  #---------------------------------------------------------------------------------------
  #Start simulation loop here  
  i <-0
  for (k in 1:nSims){
    
    for (j in mylo.p){ #range of cutoffs to loop through
      
      i<-i+1 #increase index for summary table row
      proplo<-j#proportion drawn from group with low scores
      
      #record basic parameters for this run
      summarytable[i,1]<-myN
      summarytable[i,2]<-myr
      summarytable[i,3]<- mycutoff 
      summarytable[i,4]<-proplo #proportion of low ('late ascertained') cases this run
      
      lobit<-filter(mydata,pheno<mycutoff) #Cases that are below cutoff
      n.sel<-as.integer(myN*j) #N cases selected because below cutoff
      s<-sample(nrow(lobit),n.sel,replace=FALSE) #sample this N cases without replacement from those below cutoff
      selsample<-lobit[s,]
      
      n.unsel<-myN-n.sel #Remainder of cases are unselected sample from big population
      u<-sample(nrow(mydata),n.unsel,replace=FALSE) #select without replacement
      usample<-mydata[u,] #unbiased sample from all of mydata
      
      usample$source<-0 #keep a record of origin of subject as 0 or 1
      selsample$source<-1
      
      totsample<-rbind(selsample,usample) #create single sample with selected + unselected
      
      rtest<-cor.test(totsample$pheno,totsample$mygeno) #now check the correlation in the new sample
      myp<-rtest$p.value #alternative way of getting p for regression
      if (myp< .001){myp <- .001} #set floor for low p-values to help plotting
      
      #check chi square for distribution of genotypes
      Nsample<-nrow(totsample)
      myevalue[2]<-as.integer(Nsample*NaA)
      myevalue[1]<-as.integer(Nsample*Naa)
      myevalue[3]<-as.integer(Nsample*NAA)
      myovalue<-table(totsample$mygeno)
      mychi<-chisq.test(rbind(myevalue,myovalue))
      myp2<-mychi$p.value
      if (myp2< .001){myp2<- .001} #set floor for low p-values to help plotting
      
      #in summarytable, for this run, record p-values for regression and chi square analyses, 
      # as well as means for 3 genotypes
      summarytable[i,5]<-myp
      summarytable[i,6]<-myp2
      summarytable[i,7:9]<-myovalue
      mymeans<-aggregate(totsample$pheno,by=list(totsample$mygeno),FUN=mean)
      summarytable[i,10:12]<-mymeans[,2]
    }
  }
}
#----------------------------------------------------------------------------------------
summarytable$log.pr<-log10(summarytable[,5]) #logs of pvalues easier for some plots etc
summarytable$log.pchi<-log10(summarytable[,6])

myrange=1:length(mylo.p)
rdsname<-paste0('rds_N',myN,'_r_',myCorr,'.rds')
saveRDS(summarytable,file=rdsname)#save file in R format
#can read RDS back in with readRDS(file=rdsname)
#-------------------------------------------------------------
# aggregate data to get means by proportion who are late ascertained

aggdata<-aggregate(summarytable[,5:14], by=list(summarytable$proplo),
                   FUN=mean, na.rm=TRUE)
colnames(aggdata)[1]<- 'proplo'

myhead<-'Log p-values for \n regression of phenotype on genotype'
mysub<-paste('Population r =',myCorr,': Sample size =',myN)
plot(aggdata$proplo,aggdata$log.pr,main=myhead,sub=mysub,
     xlab='% late-ascertained',ylab='log10 p-value')
axis(side=1,at=proplo)
lines(aggdata$proplo,aggdata$log.pr,col='blue',type='o',pch=18)
abline(h=log10(.05),col='red',lty=2)


#rdsname<-'rds_N120_r_0.25.rds' #can substitute name of previous saved file here
rdsname<-paste0('rds_N',myN,'_r_',myCorr,'.rds')
allpower1<-allpower2<-0 #initialise
readRDS(file=rdsname)
thisloop<-0
for (j in mylo.p){ #range of proportions to loop through
  thisloop<-thisloop+1
  mybit<-filter(summarytable, proplo==j)
  thisN<-nrow(mybit)
  #for power for 1 tailed test, find runs where p.r < .1
  if (thisloop==1)
  {
    allpower2<-as.integer(100*length(which(mybit$p.r<.05))/thisN)
    allpower1<-as.integer(100*length(which(mybit$p.r<.1))/thisN)}
  if (thisloop>1){
    allpower1<-c(allpower1,as.integer(100*length(which(mybit$p.r<.1))/thisN))
    allpower2<-c(allpower2,as.integer(100*length(which(mybit$p.r<.05))/thisN))
  }
}
aggdata$power1<-allpower1
aggdata$power2<-allpower2

View(aggdata)

#Check mean difference for selected and unselected
mean.unselected<-(aggdata[1,4]*aggdata[1,7]+aggdata[1,5]*aggdata[1,8]+aggdata[1,6]*aggdata[1,9])/
  sum(aggdata[1,4:6])
mean.selected<-(aggdata[4,4]*aggdata[4,7]+aggdata[4,5]*aggdata[4,8]+aggdata[4,6]*aggdata[4,9])/
  sum(aggdata[4,4:6])

