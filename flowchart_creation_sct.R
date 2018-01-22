#http://rich-iannone.github.io/DiagrammeR/graphviz_and_mermaid.html
# DVM Bishop 27th October 2017
# modified by Paul Thompson on 19th Jan 2018

# From sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Sierra 10.12.6

# Based on previous versions done by Paul Thompson/Alex Wilson
# This one does flowchart for genetics study with both SCT and twins

library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(svglite)
library(rsvg)
library(png)
library(tidyverse)
library(stringr)
#Function to detect if OS is PC or Macintosh and set path accordingly
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
os1=get_os()
dir<-"/Users/dorothybishop/Dropbox/ERCadvanced/project SCT analysis/Data from redcap/"
if(os1=="windows"){
  dir<-"C:\\Users\\pthompson\\Dropbox\\project SCT analysis\\Data from redcap\\"}

#Sorry had to add this bit as the link didn't work for me. i.e. manual down load and redirect link.
#dir.PT<-"c:/Users/pthompson/Desktop/"

main.data <- read.csv(paste0(dir,"SCTData_DATA_2018-01-22_1351.csv"))
main.data<-filter(main.data,trisomy<9) #remove one isochromosome case
names(main.data)[1]<-"record_id"

#Referral source - revised for genetics paper: NHS referral not relevant here
# nhscases<-c(70, 212, 213, 214, 215, 218, 221, 222, 223, 225, 229, 231,
#             233, 234, 237, 239, 240, 241, 247, 251, 252, 254, 255, 257,
#             258, 264, 265, 268, 270, 273, 274, 276, 283, 285, 297, 304,
#             306, 307, 310, 318, 319, 321, 322, 323, 327, 328, 329, 331,
#             334, 335, 336, 337, 340, 342, 346, 348, 351, 352, 355, 357,
#             358, 360, 362, 365, 367, 368, 369, 370)
# nhsset<-filter(main.data,record_id %in% nhscases)
# othset<-setdiff(main.data,nhsset)

#for ease of computation, record dna_ok values of 9 to zero
#(for 9 Di did not have a sample ;for 0 there was sample but poor quality)
w<-which(main.data$dna_ok==9)
main.data$dna_ok[w]<-0

    ascbias<-which(main.data$why_tested==3)
    group2<-main.data[ascbias,] #these are divided
    #according to whether tested because of neurodev concerns
    group1<-main.data[-ascbias,] #tested bcs other med concerns
    label1<-'Medical'

 
  n.A<-dim(group2)[1]
  n.B<-dim(group1)[1]
  
  y1=50;y2=34;y3=58
  y1<-length(which(group1$pre_postnatal_diag==0))
  y2<-length(which(group1$pre_postnatal_diag==1))
  y3<-nrow(group2)
  no_dna.B<-length(which(group2$dna_ok==0))
  no_dna.A<-length(which(group1$dna_ok==0))
  #subset by karyotype
  xxx2<-subset(group2,trisomy==1 & dna_ok==1)
  xxy2<-subset(group2,trisomy==2 & dna_ok==1)
  xyy2<-subset(group2,trisomy==3 & dna_ok==1)
  xxx1<-subset(group1,trisomy==1 & dna_ok==1)
  xxy1<-subset(group1,trisomy==2 & dna_ok==1)
  xyy1<-subset(group1,trisomy==3 & dna_ok==1)
  
  n.C<-dim(xxx1)[1]
  n.D<-dim(xxy1)[1]
  n.E<-dim(xyy1)[1]
  n.F<-dim(xxx2)[1]
  n.G<-dim(xxy2)[1]
  n.H<-dim(xyy2)[1]
  

  #now create flow chart ; TB denotes top to bottom
  #Need to add labels along the side: top row 'Reason for testing or Time of testing'
  #Then trisomy, then with DAWBA data
  
  myflow1<-grViz("
              digraph a_nice_graph {
              
              # node definitions with substituted label text
              node [shape = plaintext, fontname = Helvetica]
              X[label='@@1']
              Y[label= 'Excluded: DNA quality']
              Z[label = 'Trisomy']
              
              node [shape=square]
              A[label='@@2',fillcolor=lightBlue]
              B[label='@@3']
              
              C[label='@@4']
              D[label='@@5']
              E[label='@@6']
              F[label='@@7']
              G[label='@@8']
              H[label='@@9']
              I[label='@@10']
              J[label='@@11']
     
 
              
              # edge definitions with the node IDs
              A -> C
              B -> D
              C -> {E,F,G}
              D -> {H,I,J}
             
              
              
              X -> Y [alpha=0,color='white']
              Y -> Z [alpha=0,color='white']
              }
              
              [1]: paste0('Reason for testing\\nTime of diagnosis')
              [2]: paste0('Medical\\n', 'Prenatal: N = ',y1,'\\n', 'Postnatal: N = ',y2)
              [3]: paste0('Neurodevelopmental\\ndisorder\\nPostnatal: N = ',y3 )
              [4]: paste0('N = ',no_dna.A)
              [5]: paste0('N = ',no_dna.B)
              [6]: paste0('XXX','\\n', 'N = ',n.C)
              [7]: paste0('XXY',' \\n', 'N = ',n.D)
              [8]: paste0('XYY',' \\n', 'N = ',n.E)
              [9]: paste0('XXX',' \\n', 'N = ',n.F)
              [10]: paste0('XXY',' \\n', 'N = ',n.G)
              [11]: paste0('XYY',' \\n', 'N = ',n.H)
       
              ")
  
  

myflow1 %>% export_svg %>% charToRaw %>% rsvg %>% png::writePNG("SCT_flow_20180121.png")

#Now do plot for twin data

dir<-"/Users/dorothybishop/Dropbox/ERCadvanced/project SCT analysis/Data from redcap/"
if(os1=="windows"){
  dir<-"C:\\Users\\pthompson\\Dropbox\\project twin kids\\Project_files\\Data\\"}

#dir.PT<-"c:/Users/pthompson/Desktop/"

twin.data <- read.csv(paste0(dir,"TwinsData_DATA_2018-01-22_0744.csv"))
w<-which(twin.data$dna_ok==9)
twin.data$dna_ok[w]<-0
twin.data$zygo2<-twin.data$zygosity
fam_ids<-unique(twin.data$fam_id)
w<-which(twin.data$zygo2==3)
twin.data$zygo2[w]<-2 #zygo2 is zygosity just as MZ and DZ
#now categorise according to whether one or both twins has parental concern
#relevant variable is splang_conc, coded 1-3 for varying degrees of lang severity and 4 for reading
#any values 1-4 are regarded as concern
twin.data$parconcern<-twin.data$splang_conc
w<-which(twin.data$splang_conc>0)
twin.data$parconcern[w]<-1
#twin.data$parent_rep_12 captures whether 0, 1 or 2 children in pair raise concern
for (i in 1:length(fam_ids)){
twin.data$parent_rep_12[which(twin.data$fam_id==fam_ids[i])]<-sum(twin.data$parconcern[which(twin.data$fam_id==fam_ids[i])])
}
twin.data$parent_rep_pair<-twin.data$parent_rep_12
#now recode so as just to show whether any concern in either twin
w<-which(twin.data$parent_rep_12>0)
twin.data$parent_rep_pair[w]<-1
excl1<-nrow(filter(twin.data,parent_rep_pair==0,dna_ok==0))
excl2<-nrow(filter(twin.data,parent_rep_pair==1,dna_ok==0))
noconcern<-length(which(twin.data$parent_rep_pair==0))/2 #i.e. both twins have zero concern
concern<-length(which(twin.data$parent_rep_pair>0))/2 #i.e. one or both twins have concern

mytwins<-filter(twin.data,dna_ok>0)#both members per pair, exclude those with no DNA

mygirls<-filter(mytwins,female==1)
myt<-table(mygirls$zygo2,mygirls$parent_rep_pair)
colnames(myt)<-c('none','concern')
rownames(myt)<-c('MZ','DZ')
n.MZn.g<-myt[1,1]
n.DZn.g<-myt[2,1]
n.MZc.g<-myt[1,2]
n.DZc.g<-myt[2,2]
myboys<-filter(mytwins,female==0)
myt1<-table(myboys$zygo2,myboys$parent_rep_pair)
colnames(myt1)<-c('none','concern')
rownames(myt1)<-c('MZ','DZ')
n.MZn.b<-myt1[1,1]
n.DZn.b<-myt1[2,1]
n.MZc.b<-myt1[1,2]
n.DZc.b<-myt1[2,2]
#now create flow chart ; TB denotes top to bottom

label3<-'Parental concern re language'

print(grViz("
            digraph twinflow {
            
            # node definitions with substituted label text
            node [shape = plaintext, fontname = Helvetica]
            X[label='@@1']
             Y[label= 'Exclude: No DNA']
            Z[label= 'Total N twin children']
            
            
            node [shape=square]
            A[label='@@2',fillcolor=lightBlue]
            B[label='@@3']
            
            C[label='@@4']
            D[label='@@5']
            E[label='@@6']
            F[label='@@7']
            G[label='@@8']
            H[label='@@9']
            
            
            # edge definitions with the node IDs
            A ->{G}
            B ->{H}
            G -> {C D}
            H -> {E F}
            
            
            
            X -> Y [alpha=0,color='white']
            Y -> Z [alpha=0,color='white']
            }
            
            [1]: paste0(label3, ':\\n ',' ')
            [2]: paste0('Neither twin:\\n N = ',noconcern,' pairs')
            [3]: paste0('One or both\\n twins:\\n N = ',concern, ' pairs')
            [4]: paste0('MZ:\\nN = ',n.MZn.b,' boys\\n',n.MZn.g,' girls')
            [5]: paste0('DZ:\\nN = ',n.DZn.b,' boys\\n',n.DZn.g,' girls')
            [6]: paste0('MZ:\\nN = ',n.MZc.b,' boys\\n',n.MZc.g,' girls')
            [7]:  paste0('DZ:\\nN = ',n.DZc.b,' boys\\n',n.DZc.g,' girls')
            [8]: paste0('N children = ',excl1)
            [9]: paste0('N children = ',excl2)
            
            "))

#flow %>% export_svg %>% charToRaw %>% rsvg %>% png::writePNG("twins_flow_180115.png")

#Counts of ASD, ID, and hearing.
library(gmodels)

CrossTable(twin.data$neurodev_diagnosis,twin.data$lang_probs,prop.r = F,prop.c = F,prop.t = F,chisq = F,expected = F,prop.chisq = F)


#Save information on random twin inclusion for those with zygostiy data
dibt<-select(twin.data,record_id,zygo2,randomtwininclude)
dibt<-filter(dibt,zygo2<9)
write_csv(dibt,'twinatrandom.csv')