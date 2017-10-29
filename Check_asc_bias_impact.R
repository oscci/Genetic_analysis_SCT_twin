# Checking impact of ascertainment bias on regression data

# Our SCT sample is mixture of prenatal and postnatal identified cases
# We know the latter will have more impairment
# Want to check the impact on power etc of having this mixture
# The problem is similar to the one about using cutoffs for groups that was discussed in SQING

# Program based on SQING Simulation of genotype-phenotype relations by DVM Bishop
# See script 'simulating_genopheno_cutoffs.R' on https://osf.io/nxspw
# Note, however, that in SQING we looked at sampling bias that led to restriction of range
# in phenotypes. Here we look at oversampling of impaired cases that does not affect the range.

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
#-------------------------------------------------------------------------
# Specify parameters to create correlated variables
#-------------------------------------------------------------------------
nVar<-2 #number of simulated variables 
myM<-0 # Mean score for simulated variables
myVar<-1 #Variance for simulated variables
myN<-125 #set sample size per group (You can vary this to see the effect)
i <-0
nSims<-5000 #arbitary N simulations
mycutoff<- -1 #final group will be a mix of general population and below cutoff
mylo.p<-c(.01,.2,.4,.6,.8) #proportion of cases selected as below cutoff (we will loop through these values)
corrlist<-c(.25) #actual correlation between vars
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
  colnames(mydata)<-c('myz','pheno')
  
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
  for (k in 1:nSims){
    
    for (j in mylo.p){ #range of cutoffs to loop through
      
      i<-i+1 #increase index for summary table row
      proplo<-j#proportion drawn from group with low scores
      
      #record basic parameters for this run
      summarytable[i,1]<-myN
      summarytable[i,2]<-myr
      summarytable[i,3]<- mycutoff 
      summarytable[i,4]<-proplo #proportion of low ('late ascertained') cases this run
      
      lobit<-filter(mydata,pheno<j) #Cases that are below cutoff
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



colnames(aggdata)[5:7]<-c('aa','aA','AA')

#rdsname<-'rds_N120_r_0.25.rds' #can substitute name of previous saved file here
rdsname<-paste0('rds_N',myN,'_r_',myCorr,'.rds')
allpower<-0 #initialise
readRDS(file=rdsname)
thisloop<-0
for (j in mylo.p){ #range of proportions to loop through
  thisloop<-thisloop+1
  mybit<-filter(summarytable, proplo==j)
  thisN<-nrow(mybit)
  #for power for 1 tailed test, find runs where p.r < .1
  if (thisloop==1)
  {allpower<-as.integer(100*length(which(mybit$p.r<.1))/thisN)}
  if (thisloop>1){
    allpower<-c(allpower,as.integer(100*length(which(mybit$p.r<.1))/thisN))
  }
}
aggdata$power<-allpower


