# Checking impact of ascertainment bias on regression data

# Our SCT sample is mixture of prenatal and postnatal identified cases
# We know the latter will have more impairment
# Want to check the impact on power etc of having this mixture
# The problem is similar to the one about using cutoffs for groups that was discussed in SQING

# Program based on SQING Simulation of genotype-phenotype relations by DVM Bishop
# This version started 18th October 2017

# Conclusion from simulation: inclusion of selected cases does reduce p-values for 
# regression of phenotype on genotype, but the effect is small enough not to give concern.
#-------------------------------------------------------------------------

library(MASS) #for mvrnorm function to make multivariate normal distributed vars
library(tidyverse)
options(scipen=999) #disable scientific notation.
#-------------------------------------------------------------------------
# Specify parameters to create correlated variables
#-------------------------------------------------------------------------
nVar<-2 #number of simulated variables 
myM<-0 # Mean score for simulated variables
myVar<-1 #Variance for simulated variables
myN<-140 #set sample size per group (You can vary this to see the effect)
i <-0
nSims<-5000 #arbitary N simulations
mycutoff<- -1 #differs from SQING bcs final group is a mix of cutoff and noncutoff
mylo.p<-c(.01,.2,.4,.6,.8) #proportion of cases selected as below cutoff
corrlist<-c(.25)
summarytable<-data.frame(matrix(rep(NA,12*length(mylo.p)*nSims),ncol=12)) #initialise table to hold results
colnames(summarytable)<-c('Nsub','truer','cutoff','proplo','p.r','p.chi','N_aa','N_Aa','N_AA',
'mean_aa','mean_Aa','mean_AA')
for (myCorr in corrlist){ #correlation between variables; can loop through various values
  #----------------------------------------------------------------------------------------
  # Generate a sample from a multivariate normal distribution with the specified correlation matrix
  #----------------------------------------------------------------------------------------
  
  myCov<-matrix(rep(myCorr,nVar*nVar),nrow=nVar) #rep(x,y) means generate y values of x
  diag(myCov)<-rep(myVar,nVar)  
  mydata<-data.frame(mvrnorm(n = myN*100, rep(myM,nVar), myCov)) #make big population to select from
  colnames(mydata)<-c('myz','pheno')
  MAF<-.5 # minor allele frequency
  # convert the random normal deviates to genotype values of 0, 1 or 2, depending on MAF
  Naa<-MAF*MAF
  NAA<-(1-MAF)^2
  NaA<-1-Naa-NAA
  mydata$mygeno<-1
  critA<--qnorm(NAA)
  mydata[mydata$myz>critA,3]<-2
  crita<-qnorm(Naa)
  mydata[mydata$myz<crita,3]<-0
  
  #test significant of correlation of genotype/phenotype for full population
  rtest<-cor.test(mydata$pheno,mydata$mygeno)
  myr<-rtest$estimate
  
  regtest<-summary(lm(pheno~mygeno,data=mydata))
  myp<-regtest$coefficients[8] #pvalue for linear regression of pheno on geno
  myevalue<-c(as.integer(myN/4),as.integer(myN/2),as.integer(myN/4)) #expected distribution of genotypes
  myovalue<-table(mydata$mygeno[1:myN]) #observed distribution of genotypes
  mychi<-chisq.test(rbind(myevalue,myovalue))
  
for (k in 1:nSims){

for (j in mylo.p){ #range of cutoffs to loop through

i<-i+1 #increase index for summary table row
proplo<-j#proportion drawn from group with low scores
summarytable[i,1]<-myN
summarytable[i,2]<-myr
summarytable[i,3]<- mycutoff 
summarytable[i,4]<-proplo #proportion of low ('late ascertained') cases this run

lobit<-filter(mydata,pheno<j)
n.sel<-as.integer(myN*j)
s<-sample(nrow(lobit),n.sel,replace=FALSE)
selsample<-lobit[s,]
n.unsel<-myN-n.sel
u<-sample(nrow(mydata),n.unsel,replace=FALSE)
usample<-mydata[u,] #unbiased sample from all of mydata
usample$source<-0
selsample$source<-1
totsample<-rbind(selsample,usample)
rtest<-cor.test(totsample$pheno,totsample$mygeno)
myp<-rtest$p.value #alternative way of getting p for regression
if (myp< .001){myp <- .001} #set floor for low p-values to help plotting
Nsample<-nrow(totsample)
myevalue[2]<-as.integer(Nsample*NaA)
myevalue[1]<-as.integer(Nsample*Naa)
myevalue[3]<-as.integer(Nsample*NAA)
myovalue<-table(totsample$mygeno)
mychi<-chisq.test(rbind(myevalue,myovalue))
myp2<-mychi$p.value
if (myp2< .001){myp2<- .001} #set floor for low p-values to help plotting
summarytable[i,5]<-myp
summarytable[i,6]<-myp2
summarytable[i,7:9]<-myovalue
mymeans<-aggregate(totsample$pheno,by=list(totsample$mygeno),FUN=mean)
summarytable[i,10:12]<-mymeans[,2]
}
}
}
#----------------------------------------------------------------------------------------
summarytable$log.pr<-log10(summarytable[,5])
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


#Illustrate how means and N change with cutoffs at 3 levels
# par(mfrow=c(1,3)) #set up plot frame with 3 plots in 1 row
# 
# barplot(as.matrix(aggdata[1,5:7]),main='No cutoff',ylim=c(-.2,1.25),ylab='zscore')
# text(.6,.04,paste('N = ',as.integer(aggdata[1,2])),cex=.8)
# text(1.8,.08,paste('N = ',as.integer(aggdata[1,3])),cex=.8)
# text(3.1,.35,paste('N = ',as.integer(aggdata[1,4])),cex=.8)
# 
# barplot(as.matrix(aggdata[3,5:7]),main='Cutoff = 0',ylim=c(-.2,1.25))
# text(.6,.8,paste('N = ',as.integer(aggdata[3,2])),cex=.8)
# text(1.8,.85,paste('N = ',as.integer(aggdata[3,3])),cex=.8)
# text(3.1,.95,paste('N = ',as.integer(aggdata[3,4])),cex=.8)
# 
# barplot(as.matrix(aggdata[4,5:7]),main='Cutoff = .5',ylim=c(-.2,1.25))
# text(.6,1.15,paste('N = ',as.integer(aggdata[4,2])),cex=.8)
# text(1.8,1.18,paste('N = ',as.integer(aggdata[4,3])),cex=.8)
# text(3.1,1.235,paste('N = ',as.integer(aggdata[4,4])),cex=.8)


readRDS(file=rdsname)
for (j in mylo.p){ #range of proportions to loop through
mybit<-filter(summarytable, proplo==j)
thisN<-nrow(mybit)
if (j==.01)
{allpower<-as.integer(100*length(which(mybit$p.r<.1))/thisN)}
if (j>.01){
allpower<-c(allpower,as.integer(100*length(which(mybit$p.r<.1))/thisN))
}
}
aggdata$power<-allpower


