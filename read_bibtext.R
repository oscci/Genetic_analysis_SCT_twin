#Identifying genes for study
#Program by DVM Bishop, started 30/09/17

#Initial gene list by Dianne and Nuala selected on basis of association with language/literacy/laterality.
#This is 149 genes.xlsx

readdir<-"~/Dropbox/ERCadvanced/project SCT analysis/SCT genetic analysis/biblio_data/"
writedir<-readdir

#search term in Scopus
#TITLE-ABS-KEY(ATP2C2 OR ROBO1 OR DCDC2 OR C2ORF3 OR MRPL19 OR AUTS2 OR BDNF OR GRIN2B OR FOXP2 OR CNTNAP2 
#OR KIAA0319 OR DYX1C1 OR CMIP OR DRD2 OR APOE OR ATP13A4 OR ASPM OR androgen receptor OR NLGN4X OR NLGN4Y 
#OR CCKAR OR AP4E1 OR ARID1B OR S100B OR CCDC136 OR FLNC OR ERC1 OR NRXN1 OR DIP2A OR FMR1 OR BCL11A OR FOXP1
#OR SETBP1 OR DISC1 OR GRM3 OR GRIN2A OR COMT OR PCSK6 OR SRPX2 OR LRRTM1 OR GNPTG OR NAGPA OR RBFOX2 OR 
#GNPTAB OR DRD4 OR ELP4 OR NOP9) 
#AND KEY(synap* OR neurexin OR neuroligin OR autism OR language OR dyslexia OR reading OR SLI) 
#AND KEY(human) AND NOT TITLE(review)

#On 1st Oct 2017 generated 1203 results.

bibfile<-'gene_synapse_lang.bib'

require(bibliometrix)
#https://cran.r-project.org/web/packages/bibliometrix/vignettes/bibliometrix-vignette.html
require(tidyverse)

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

genelist<-c('ATP2C2', 'ROBO1', 'DCDC2', 'C2ORF3', 'MRPL19', 'AUTS2', 'BDNF', 'GRIN2B', 'FOXP2',
            'CNTNAP2', 'KIAA0319', 'DYX1C1', 'CMIP', 'DRD2', 'APOE', 'ATP13A4', 'ASPM', 'ANDROGEN RECEPTOR', 
            'NLGN4X', 'NLGN4Y', 'CCKAR', 'AP4E1', 'ARID1B', 'S100B', 'CCDC136', 'FLNC', 'ERC1', 
            'NRXN1', 'DIP2A', 'FMR1', 'BCL11A', 'FOXP1', 'SETBP1', 'DISC1', 'GRM3', 'GRIN2A', 'COMT',
            'PCSK6', 'SRPX2', 'LRRTM1', 'GNPTG', 'NAGPA', 'RBFOX2', 'GNPTAB', 'DRD4', 'ELP4', 'NOP9')
#NB AR gene is problematic because AR is embedded in many words, so spelt out here

keywordlist<-c('SYNAP','NEUREXIN','NEUROLIGIN','AUTIS','LANGUAGE','READING','SPEECH','DYSLEXIA','SPECIFIC LANGUAGE IMPAIRMENT')

ngene<-length(genelist)
nkey<-length(keywordlist)

#add columns to mydf to indicate which genes are included
addbit<-matrix(data=0,ncol=ngene,nrow=nrecords)
colnames(addbit)<-genelist
addbit2<-matrix(data=0,ncol=nkey,nrow=nrecords)
colnames(addbit2)<-keywordlist
mydf<-cbind(mydf,addbit,addbit2)

#initialise table to hold summary results
mytab<-data.frame(matrix(data=NA,ncol=nkey+2,nrow=ngene))
mytab[,2:(nkey+2)]<-0
for (i in 1:ngene){
  mytab[i,1]<-genelist[i]
}
colnames(mytab)[2:(nkey+1)]<-keywordlist
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
        if(length(grep(k,mytext)) >0)
        {mytab[i,(k+1)]<-mytab[i,(k+1)]+1
        mydf[j,(k+54)]<-mydf[j,(k+54)]+1}
      }
    }
  }
}
  
writebit<-paste0(writedir,"gene_summary.csv")
write.table(mytab, writebit, sep=",",row.names=FALSE) 
  
  