### R code from vignette source 'ASGSCA.Rnw'

###################################################
### code chunk number 1: ASGSCA.Rnw:200-227
###################################################
library(ASGSCA)
data("QCAHS")

#Names of all the observed variables: the SNPs then the traits
colnames(QCAHS)

#Extract the variables of interest
QCAHS1=data.frame(QCAHS$TaqIB,QCAHS$HindIII,QCAHS$G1302A,QCAHS$G1564A,QCAHS$G308A,QCAHS$G238A,
                  QCAHS$HDL,QCAHS$LDL,QCAHS$APOB,QCAHS$TG,QCAHS$Glucose,QCAHS$Insulin)

#Names of the observed variables used in this example
ObservedVar=c("TaqIB","HindIII","G1302A","G1564A","G308A","G238A","HDL","LDL","APOB",
              "TG","Glucose","Insulin")
colnames(QCAHS1)=ObservedVar

#Define the vector of the latent variables names
LatentVar=c("CETP","LPL","PGC","TNFa","Lipid metabolism","Energy metabolism")

#Construction of the matrices W0 and B0 describing the model illustrated in Figure 2.
#W0 is I x L matrix where rows are genotypes and traits, columns are genes and 'clinical pathways'
W0=matrix(rep(0,12*6),nrow=12,ncol=6, dimnames=list(ObservedVar,LatentVar))
W0[1,1]=W0[2,2]=W0[3:4,3]=W0[5:6,4]=W0[7:10,5]=W0[8:12,6]=1

#B0 is L x L matrix with 0s and 1s, where 1 indicates arrow from latent variable in row to latent variable in column
B0=matrix(rep(0,6*6),nrow=6,ncol=6, dimnames=list(LatentVar,LatentVar))
B0[1:3,5]=B0[3:4,6]=1

W0
B0


###################################################
### code chunk number 2: ASGSCA.Rnw:235-236
###################################################
GSCA(QCAHS1,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=FALSE,path=NULL,nperm=1000)


###################################################
### code chunk number 3: ASGSCA.Rnw:256-258
###################################################
set.seed(2)
GSCA(QCAHS1,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=TRUE,path=NULL,nperm=1000)


###################################################
### code chunk number 4: ASGSCA.Rnw:269-271
###################################################
set.seed(2)
GSCA(QCAHS1,W0, B0,latent.names=LatentVar,estim=FALSE,path.test=TRUE,path=NULL,nperm=1000)


###################################################
### code chunk number 5: ASGSCA.Rnw:280-286
###################################################
set.seed(2)

path0=matrix(c(2,3,5,6),ncol=2)
path0
GSCA(QCAHS1,W0, B0,latent.names=LatentVar, estim=FALSE,path.test=TRUE,path=path0,
     nperm=1000)


###################################################
### code chunk number 6: ASGSCA.Rnw:297-309
###################################################
ObservedVar=colnames(QCAHS)
ObservedVar

#Define the vector of the latent variables names
LatentVar=c("CETP","APOC3","ABCA1","FABP-2","APOA1","APOE","HL","LPL","MTP","PON1","PON2","PCSK9",
            "PGC","ADIPO","PPARg2","TNFa","eNOS","a23AR","b1AR","b2AR","b3AR","ACE","AGT","AGTR1","LEPR",
            "Lipid metabolism", "Energy metabolism","BP control")

#The matrices W0 and B0 describing the model illustrated in Figure 2.
data(W0); data(B0)
dim(W0)
dim(B0)


###################################################
### code chunk number 7: ASGSCA.Rnw:326-336
###################################################
#set.seed(4)
#ResQCAHS=GSCA(QCAHS,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=TRUE,path=NULL,nperm=1000)
data("ResQCAHS")
indices <- which(ResQCAHS$pvalues<0.05, arr.ind=TRUE)
ind.row=indices[,1] ; ind.col=indices[,2]
Significant<- matrix(rep(0,nrow(indices)*3),ncol=3);colnames(Significant)=c("Gene", "Pathway", "pval")
Significant[,1] <- rownames(ResQCAHS$pvalues)[ind.row]
Significant[,2] <- colnames(ResQCAHS$pvalues)[ind.col]
Significant[,3]<-ResQCAHS$pvalues[indices]
Significant


