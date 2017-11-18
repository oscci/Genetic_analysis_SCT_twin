#Reading in _papers files to summarise coded information re papers
library(ggplot2)
library(tidyverse)
#Can't get this to print all on one pdf - not sure why

readdir<-"~/Dropbox/ERCadvanced/project SCT analysis/SCT genetic analysis/biblio_data/gene_specific_csv/"
writedir<-readdir
targetlist<-c('NRXN1','DLGAP2','GRIN2B','CNTNAP2','FOXP2','BCL11A','ATP2C2')
ntargets<-length(targetlist)

#mypng<-"geneplots.png"
#png(mypng,width=10,height=14)
#quartz()
#par(mfrow=c(7,2))

for (i in 1:ntargets)
{
  myname<-targetlist[i]
  
  genefile<-paste0(readdir,targetlist[i],'_papers_done.csv')
  mygenedata<-read.csv(genefile,stringsAsFactors=FALSE)
  #---------------------------------------------------------------------------
  #codes are: C: cell
  #           A: animal model
  #           S: case study
  #           M: group study of mutation/rare variant/CNV
  #           V: group study of common variant
  #           O: other (includes reviews/not relevant/etc)
  
  w<-which(mygenedata$Code=='A')
  mygenedata$Code[w]<-'Animal model'
  w<-which(mygenedata$Code=='C')
  mygenedata$Code[w]<-'Cells'
  w<-which(mygenedata$Code=='S')
  mygenedata$Code[w]<-'Case study'
  w<-which(mygenedata$Code=='M')
  mygenedata$Code[w]<-'Rare variants/CNV'
  w<-which(mygenedata$Code=='V')
  mygenedata$Code[w]<-'Common variants'
  w<-which(mygenedata$Code=='O')
  mygenedata$Code[w]<-'Other'
  
  mytab<-data.frame(table(mygenedata$Code))
  
  colnames(mytab)[1]<-'Type'
  #reorder factor levels
  mytab$Type = factor(mytab$Type,
                      levels(mytab$Type)[c(3,1,2,6,4,5)])
  #also need to reorder the rows of mytab
  mytab<-mytab[c(3,1,2,6,4,5),]
  # Add additional columns
  mytab<-filter(mytab,Freq>0)
  mytab$ymax = cumsum(mytab$Freq)
  mytab$ymin = c(0, head(mytab$ymax, n=-1))
  mynpapers<-sum(mytab$Freq,na.rm=TRUE)
  mymidlab1<-paste0(mynpapers,' papers')
  
  p2=ggplot(mytab, aes(fill=Type, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
    geom_rect(colour="grey30") +
    coord_polar(theta="y") +
    xlim(c(0, 4)) +
    theme_bw() +
    geom_text( aes(label = mymidlab1,x=0, y=0))+
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    theme(axis.title=element_blank())+
    theme(panel.grid=element_blank())+
    theme(panel.border=element_blank())+
    labs(title=myname)
  p2
  #Now focus just on the studies of rare and common variants
  #------------------------------------------------------------------------
  colnames(mygenedata)[12:15]<-c('ID','EPILEPSY','SZ','other') #to ensure compatible names
  mygenepheno<-rbind(filter(mygenedata,Code=='Rare variants/CNV'),filter(mygenedata,Code=='Common variants'))
  mynpapers<-nrow(mygenepheno)
  if (mynpapers>1){
    mypheno=data.frame(c(sum(mygenepheno$LANGUAGE),sum(mygenepheno$AUTIS),sum(mygenepheno$ID),
                         sum(mygenepheno$SZ),sum(mygenepheno$EPILEPSY),sum(mygenepheno$other)))
    colnames(mypheno)<-'Freq'
    mypheno$Phenotype<-as.factor(c('Language','ASD','Intellectual disability','Schizophrenia','Epilepsy','Other'))
    mypheno$Phenotype = factor(mypheno$Phenotype ,
                               levels(mypheno$Phenotype)[c(4,1,3,6,2,5)])
    mypheno<-filter(mypheno,Freq>0)
    mypheno$ymax = cumsum(mypheno$Freq)
    mypheno$ymin = c(0, head(mypheno$ymax, n=-1))
    mynmentions<-sum(mypheno$Freq,na.rm=TRUE)
    mymidlabel<-paste0(mynmentions,'\nphenotypes in \n',mynpapers,' papers')
    
    p3=ggplot(mypheno, aes(fill=Phenotype, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
      scale_fill_brewer(palette="RdBu")+
      geom_rect(colour="grey30") +
      coord_polar(theta="y") +
      xlim(c(0, 4)) +
      theme_bw() +
      geom_text( aes(label = mymidlabel,x=0, y=0))+
      theme(panel.grid=element_blank()) +
      theme(axis.text=element_blank()) +
      theme(axis.ticks=element_blank()) +
      theme(axis.title=element_blank())+
      theme(panel.grid=element_blank())+
      theme(panel.border=element_blank())+
      labs(title='')
    p3
  }
  #plot_grid(p2,p3, labels = c('All', 'Human group studies'))
  
}
#dev.off()