#Aim is to depict contrast between Double Hit model and Increased CNV burden model for sex chromosome trisomies.
#Simulate distribution of CNVs for both models, with probabilistic rule linking CNV size to presence of phenotype
#Currently 10 binary phenotypes which are summed with added error to give continuous score.

#NB CNV size currently modeled as normally distributed; would be better if modelled to 
#reflect fact that the freq goes down with size, as discussed here:
# Conrad et al: 2006 Nature genetics  38(1):75-81.: A high-resolution survey of deletion polymorphism in the human genome.

#See also this review: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4472309/
#database where structural variation is cataloged (the Database of Genomic Variants or DGV,
#http://projects.tcag.ca/variation/).

require(MASS)
myvars<-c('sct','sizecnvA','sizecnvB','phenoA1','phenoA2','phenoA3','phenoA4','phenoA5',
          'phenoA6','phenoA7','phenoA8','phenoA9','phenoA10',
          'phenoAsum','phenoB1','phenoB2','phenoB3','phenoB4','phenoB5',
          'phenoB6','phenoB7','phenoB8','phenoB9','phenoB10','phenoBsum')

#set up dataframe to save simulation results
ncases<-120 #number of cases each for sct and control
mydata<-data.frame(matrix(NA,nrow=ncases*2,ncol=length(myvars)))
colnames(mydata)<-myvars
range1<-1:ncases
range2<-(ncases+1):(ncases*2)
mydata[range1,1]<-0
mydata[range2,1]<-1 #half control and half with trisomy
mymean<-0; mysd<-1; #mean and SD for size of cnv - arbitrary
mydata[,2]<-rnorm(ncases*2) #random number ZSCORE represents size cnv
#model B - higher rate of cnv in trisomy cf control - add .8 sd to cnv size
mydata[range1,3]<-mydata[range1,2]
mydata[range2,3]<-mydata[range2,2]+mysd*.8

#model A: same rate of CNV but bigger effect with SCT - Double Hit
#relation between sizecnv and neurodev is probabilistic; use lower cutoff for disorder for sct

  pvect<-pnorm(mydata$sizecnvA) #sizecnvA at this point is a zscore, so can just convert to pvalue
  for (i in 1:10){ #phenotypes model A (have 10, each independently determined)
    rvect<-runif(ncases) #random number from 0 to 1
    dvect<-pvect-rvect #if pvect value is high, then v likely to be +ve and vice versa
    mydata[,(i+3)]<-0 #whether or not phenotype present for this variable
    mycutoff<-.7 #the higher this value, the lower the correl between cnv/pheno
    w<-which(dvect>mycutoff) #if dvect is above cutoff phenotype present
    mydata[w,(i+3)]<-1
    sctcutoff<-.3 #can try changing this and mycutoff to see how relationship affected
    w<-which(dvect[range2]> sctcutoff) #different cutoff for sct cases
    mydata[(w+ncases),(i+3)]<-1 #nb add ncases, because reference is just to part of vector in formula
  }
  mydata$phenoAsum<-rowSums(mydata[,4:13])+rnorm(ncases)+3 #rnorm + 3 adds random error

#add 3 to the cnv variables, so they could plausibly be Mb (ie megabase units, million base pairs)
mydata[,2:3]  <-(3+mydata[,2:3])/3-.3


#png(filename = "CNV_models.png", width = 500, height = 300)
par(mfrow=c(1,2))
  
plot(mydata$sizecnvA,mydata$phenoAsum,col=(mydata$sct+1),
     xlab='CNV burden (arbitrary units)',ylab='Phenotype severity',main='A. Double Hit model',
     ylim=c(0,15),xlim=c(0,2))
legend(-.05, 15, legend=c("Comparison", "SCT"),
       col=c("black", "red"), pch=1,cex=0.75)

conphenomean<-mean(mydata$phenoAsum[range1])
sctphenomean<-mean(mydata$phenoAsum[range2])
abline(h=conphenomean, col='black')
abline(h=sctphenomean, col='red')
concnvmeana<-mean(mydata$sizecnvA[range1])
sctcnvmeana<-mean(mydata$sizecnvA[range2])
abline(v=concnvmeana, col='black')
abline(v=sctcnvmeana, col='red')


#model B: higher rate of CNV but same effect with SCT - added burden model
#otherwise identical : relation between sizecnv and neurodev is probabilistic; use lower cutoff for disorder for sct

pvect<-pnorm(mydata$sizecnvB)
for (i in 1:10){ #phenotypes model B
  rvect<-runif(ncases)
  dvect<-pvect-rvect
  mydata[,(i+14)]<-0
  w<-which(dvect>mycutoff)
  mydata[w,(i+14)]<-1
}
mydata$phenoBsum<-rowSums(mydata[,15:24])+rnorm(ncases)+3


plot(mydata$sizecnvB,mydata$phenoBsum,col=(mydata$sct+1),pch=1,
     xlab='CNV burden (arbitrary units)',ylab='Phenotype severity',
     main='B. Increased CNV Burden model', ylim=c(0,15),xlim=c(0,2))
legend(-.05, 15, legend=c("Comparison", "SCT"),
       col=c("black", "red"), pch=1,cex=0.75)
conphenomean<-mean(mydata$phenoBsum[range1])
sctphenomean<-mean(mydata$phenoBsum[range2])
abline(h=conphenomean, col='black')
abline(h=sctphenomean, col='red')
concnvmeanb<-mean(mydata$sizecnvB[range1])
sctcnvmeanb<-mean(mydata$sizecnvB[range2])
abline(v=concnvmeanb, col='black')
abline(v=sctcnvmeanb, col='red')
#dev.off()

# #consider impact of testing only those with phen sev>4
# w1<-which(mydata$phenoAsum[range1]>5)
# w2<- ncases+which(mydata$phenoAsum[range2]>5)
# concnvmeanAx<-mean(mydata$sizecnvA[w1])
# sctcnvmeanAx<-mean(mydata$sizecnvA[w2])
# concnvmeanBx<-mean(mydata$sizecnvB[w1])
# sctcnvmeanBx<-mean(mydata$sizecnvB[w2])
# 
# myrA<-cor.test(mydata$sizecnvA,mydata$phenoAsum)
# myrB<-cor.test(mydata$sizecnvB,mydata$phenoBsum)

#--------------------------------------------------------------
# Power analysis
#--------------------------------------------------------------
#Assume beta for group 1 between CNV and phenotype sev is .1
#Step through values of beta for group 2 to find smallest value
#at which effect of group is seen in regression
par(mfrow=c(1,1))
options(scipen = 999) #turn off scientific notation
nsim=10
bigmat<-data.frame(matrix(NA,nrow=nsim*8,ncol=5))
names(bigmat)<-c('condition','b1','b2','p','sig')
mydat <-data.frame(matrix(NA,nrow=240,ncol=3))
colnames(mydat)<-c('group','cnv','pheno')
mydat[1:120,1]<-0 #assign group
mydat[121:240,1]<-1
mydat[,2]<-rnorm(240) #random N for cnv - needs to be positive
myrow=0
for (i in 2:2){
  for (n in 1:nsim){
    myrow<-myrow+1
    myi<-i/10
mydat[1:120,3]<-1+.1*mydat[1:120,2]+rnorm(120)/5
mydat[121:240,3]<-1+myi*mydat[121:240,2]+rnorm(120)/5
r1<-cor(mydat[1:120,2:3])[1,2]
r2<-cor(mydat[121:240,2:3])[1,2]
myreg<-summary(lm(pheno~cnv*group,data=mydat))
myp<-myreg$coefficients[4,4] #pvalue for interaction term
bigmat[myrow,1]<-i
bigmat[myrow,2]<-r1
bigmat[myrow,3]<-r2
bigmat[myrow,4]<-round(myp,3)
bigmat[myrow,5]<-0
if(myp<.05){
  bigmat[myrow,5]<-1
  }
}
}
plot(mydat$cnv,mydat$pheno,col=mydat$group+1)