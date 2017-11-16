#Identifying genes for study
#Program by DVM Bishop, started 30/09/17
#Updated 3/10/17 to include scatterplot of candidate genes
#Updated 8/10/17 to include search based on PAR genes
#Updated 15/11/17 to include additional genes from Dianne but to exclude those on X chromosome
#Additional genes from Di from protein-protein interaction database
# DLG4 (aka PSD95), NLGN3, PRIM1, PRIM2, UBB, DAP2 (aka SAPAP2, DLGAP2), NRXN1, BRK1
# NLGN3 excluded because on X

#Initial gene list by Dianne and Nuala selected on basis of association with language.
#This is 149 genes.xlsx
# From sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Sierra 10.12.6
readdir<-"~/Dropbox/ERCadvanced/project SCT analysis/SCT genetic analysis/biblio_data/"
writedir<-readdir

#search term in Scopus
# TITLE-ABS-KEY(ATP2C2 OR ROBO1 OR DCDC2 OR C2ORF3 OR MRPL19 OR AUTS2 OR BDNF OR 
# GRIN2B OR FOXP2 OR CNTNAP2 OR KIAA0319 OR DYX1C1 OR CMIP OR DRD2 OR APOE OR 
# ATP13A4 OR ASPM OR AP4E1 OR ARID1B OR S100B OR CCDC136 OR FLNC OR ERC1 OR NRXN1 
# OR BCL11A OR FOXP1 OR SETBP1 OR DISC1 OR GRM3 OR GRIN2A OR COMT OR GNPTG OR 
# NAGPA OR RBFOX2 OR GNPTAB OR DRD4 OR ELP4 OR NOP9 OR DLG4 
# OR PSD95 OR PRIM1 OR PRIM2 OR UBB OR DAP2 OR SAPAP2 OR DLGAP2 OR BRK1) 
# AND KEY(synap* OR neurexin OR neuroligin OR autism OR language OR dyslexia 
#         OR reading OR SLI) AND KEY(human) AND DOCTYPE(ar) 
# 
# 
# On 15/11/2017 generated 1903 results.


bibfile<-'gene_synapse_language_20171115.bib'
#bibfile<-'scopus-par.bib' #Use this version for PAR (pseudoautosomal region) genes

require(bibliometrix)
#https://cran.r-project.org/web/packages/bibliometrix/vignettes/bibliometrix-vignette.html
require(tidyverse)
require(ggrepel) #avoids overlap in scatterplot labels
mybib<-paste0(readdir,bibfile)
D <- readFiles(mybib)

mydf <- convert2df(D, dbsource = "scopus", format = "bibtex")
mydf<-select(mydf,AU,TI,SO,DE,ID,AB,PY) #retain these 7 columns from mydf
#the ID field has the Scopus keywords and DE has author's keywords
nrecords<-nrow(mydf)
# NetMatrix <- biblioNetwork(mydf, analysis = "co-occurrences", network = "keywords", sep = ";")
# 
# # Plot the network
# net=networkPlot(NetMatrix, n = 20, Title = "Keyword Co-occurrences", type = "kamada", size=T)

genelist<-c("ATP2C2", "ROBO1", "DCDC2", "C2ORF3", "MRPL19", "AUTS2", "BDNF", "GRIN2B", "FOXP2", 
            "CNTNAP2", "KIAA0319", "DYX1C1", "CMIP", "DRD2", "APOE", "ATP13A4", "ASPM", "AP4E1", 
            "ARID1B", "S100B", "CCDC136", "FLNC", "ERC1", "NRXN1", "BCL11A", "FOXP1", "SETBP1", 
            "DISC1", "GRM3", "GRIN2A", "COMT", "GNPTG", "NAGPA", "RBFOX2", "GNPTAB", "DRD4", 
            "ELP4", "NOP9", "DLG4 ", "PSD95", "PRIM1", "PRIM2", "UBB", "DAP2", "SAPAP2", 
            "DLGAP2", "BRK1")

# Alternative genelist used with PAR genes
# genelist<-c('AKAP17A', 'ASMT', 'ASMTL', 'ASMTL-AS1', 'CD99', 'CD99P1', 'CRLF2', 'CSF2RA', 'DHRSX', 'DHRSX-IT1', 
#             'FABP5P13', 'GTPBP6', 'IL3RA', 'LINC00102', 'LINC00106', 'LINC00685', 'MIR3690', 'MIR6089', 'PLCXD1', 
#             'PPP2R3B', 'P2RY8', 'SHOX', 'SLC25A6', 'XG', 'ZBED1', 'AMD1P2', 'DDX11L16', 'DPH3P2', 'IL9R', 'SPRY3',
#             'ELOCP24', 'TRPC6P', 'VAMP7', 'WASH6P', 'WASIR1')

keywordlist<-c('SYNAP','NEUREXIN','NEUROLIGIN','AUTIS','LANGUAGE','READING','SPEECH','DYSLEXIA','SPECIFIC LANGUAGE IMPAIRMENT')

ngene<-length(genelist)
nkey<-length(keywordlist)

#add columns to mydf to indicate which genes are included
addbit<-matrix(data=0,ncol=ngene,nrow=nrecords)
colnames(addbit)<-genelist
#colnames(addbit)[18]<-'AR' #shorten for column headings, though use full term for search
addbit2<-matrix(data=0,ncol=nkey,nrow=nrecords)
colnames(addbit2)<-keywordlist
colnames(addbit2)[9]<-'SLI'#shorten for column headings, though use full term for search
mydf<-cbind(mydf,addbit,addbit2)

#initialise table to hold summary results
mytab<-data.frame(matrix(data=NA,ncol=nkey+2,nrow=ngene))
mytab[,2:(nkey+2)]<-0
for (i in 1:ngene){
  mytab[i,1]<-genelist[i]
}
colnames(mytab)[2:(nkey+1)]<-keywordlist
colnames(mytab)[10]<-'SLI' #shortened name
colnames(mytab)[1]<-'gene'
colnames(mytab)[nkey+2]<-'papers'


#search each row for gene and, if found, for keywords
#NB for author keywords field DE, for Scopus keywords, field ID
for (j in 1:nrecords){
  mytext<-c(mydf$TI[j],mydf$AB[j],mydf$ID[j],mydf$DE[j]) #all relevant text from title/abs/keywords for this record
  for (i in 1:ngene){
    g<-genelist[i]
    if(length(grep(g,mytext)) >0)
    {
      mytab[i,11]<-mytab[i,11]+1 #add counter for total N papers in mytab
      mydf[j,(i+7)]<-mydf[j,(i+7)]+1 #mark record in main data frame to show this gene mentioned
      for (k in 1:nkey){
        thiskey<-keywordlist[k]
        if(length(grep(thiskey,mytext)) >0)
        {mytab[i,(k+1)]<-mytab[i,(k+1)]+1
        mydf[j,(k+7+ngene)]<-mydf[j,(k+7+ngene)]+1}
      }
    }
  }
}

#deal with alternative gene names
w<-which(mydf$DAP2 >0)
mydf$DLGAP2[w]<-1
w<-which(mydf$SAPAP2 >0)
mydf$DLGAP2[w]<-1

w<-which(mydf$DLG4 >0)
mydf$PSD95[w]<-1

#put cols with duplicate names to zero so they won't get counted
mydf$DAP2<-0
mydf$SAPAP2<-0
mydf$DLG4<-0



writebit<-paste0(writedir,"gene_synapse_summary_20171115.csv")
write.table(mytab, writebit, sep=",",row.names=FALSE) 

mytab<-filter(mytab,papers>2) #only consider if at least 3 papers
myrows<-nrow(mytab)
for (i in 1:myrows){

mytab$p_neuro[i]<-max(mytab$SYNAP[i],mytab$NEUREXIN[i],mytab$NEUROLIGIN[i])/mytab$papers[i]
mytab$p_lang[i]<-max(mytab$AUTIS[i],mytab$language[i],mytab$reading[i],mytab$dyslexia[i],mytab$SLI[i],mytab$SPEECH[i])/mytab$papers[i]
}
mytab2<-filter(mytab,p_neuro>0)
mytab2<-filter(mytab2,p_lang>0)
ggplot(mytab2, aes(x= p_lang, y= p_neuro, label=gene))+
  xlab('Proportion with ASD/language/reading')+
  ylab('Proportion with synapse/neurexin/neuroligin')+
  ggtitle('                Papers from SCOPUS search')+
  geom_point() +geom_text_repel(aes(label=gene),cex=3)

#To inspect articles for a specific gene, you can proceed as follows (example with NRXN1)

myNRXN1<-filter(mydf,NRXN1==1) #just those articles featuring this gene
myNRXN1<-select(myNRXN1,AU,TI,AB) #focus on author/abstract only
View(myNRXN1)
write.csv(myNRXN1,'NRXN1papers.csv')

#I then annotated these in xlsx. Can read back to categorise types of study
#---------------------------------------------------------------------------
# Plotting details of the papers
mygene<-'NRXN1'
myfilename<-paste0(mygene,'papers.csv')
mygenedata<-read.csv(myfilename,stringsAsFactors=FALSE)
#---------------------------------------------------------------------------
#codes are: C: cell
#           A: animal model
#           S: case study
#           M: group study of mutation/rare variant/CNV
#           V: group study of common variant
#           X: methods
#           RC: review with functional/cellular focus
#           RG: review with focus on phenotypes
#           O: other
#           Z: not relevant to this gene

#convert rare categories to 'Other'
w<-cbind(which(mygenedata$type=='R'),which(mygenedata$type=='RC'),
         which(mygenedata$type=='RG'),which(mygenedata$type=='X'),
         which(mygenedata$type=='O') ,which(mygenedata$type=='Z'))
mygenedata$type[w]<-'Other'

w<-which(mygenedata$type=='A')
mygenedata$type[w]<-'Animal model'
w<-which(mygenedata$type=='C')
mygenedata$type[w]<-'Cells'
w<-which(mygenedata$type=='S')
mygenedata$type[w]<-'Case study'
w<-which(mygenedata$type=='M')
mygenedata$type[w]<-'Rare variants/CNV'
w<-which(mygenedata$type=='V')
mygenedata$type[w]<-'Common variants'

mytab<-data.frame(table(mygenedata$type))

colnames(mytab)[1]<-'Type'
#reorder factor levels
mytab$Type = factor(mytab$Type,
                    levels(mytab$Type)[c(3,1,2,6,4,5)])
#also need to reorder the rows of mytab
mytab<-mytab[c(3,1,2,6,4,5),]
# Add addition columns, needed for drawing with geom_rect.
mytab$fraction = mytab$Freq / sum(mytab$Freq)
#mytab = mytab[order(mytab$fraction), ]
mytab$ymax = cumsum(mytab$Freq)
mytab$ymin = c(0, head(mytab$ymax, n=-1))
library(ggplot2)
p2 = ggplot(mytab, aes(fill=Type, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(colour="grey30") +
  coord_polar(theta="y") +
  xlim(c(0, 4)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  labs(title="NRXN1")
p2

#Now focus just on the studies of rare and common variants
mygenepheno<-rbind(filter(mygenedata,type=='Rare variants/CNV'),filter(mygenedata,type=='Common variants'))

mypheno=data.frame(c(sum(mygenepheno$lang.etc),sum(mygenepheno$asd),sum(mygenepheno$id),
          sum(mygenepheno$sz),sum(mygenepheno$epilepsy),sum(mygenepheno$other)))
colnames(mypheno)<-'Freq'
mypheno$Phenotype<-as.factor(c('Language/speech','ASD','Intellectual disability','Schizophrenia','Epilepsy','Other'))
mypheno$Phenotype = factor(mypheno$Phenotype ,
                    levels(mypheno$Phenotype)[c(4,1,3,6,2,5)])

mypheno$ymax = cumsum(mypheno$Freq)
mypheno$ymin = c(0, head(mypheno$ymax, n=-1))
p3 = ggplot(mypheno, aes(fill=Phenotype, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  scale_fill_brewer(palette="RdBu")+
  geom_rect(colour="grey30") +
  coord_polar(theta="y") +
  xlim(c(0, 4)) +
  theme_bw() +
  geom_text( aes(label = 'NRXN1',x=0, y=0), size=5, fontface="bold")+
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(axis.title=element_blank())+
  theme(panel.grid=element_blank())+
  theme(panel.border=element_blank())+
  labs(title="")
p3


