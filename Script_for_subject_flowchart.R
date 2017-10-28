#http://rich-iannone.github.io/DiagrammeR/graphviz_and_mermaid.html
# DVM Bishop 27th October 2017
# Based on previous versions done by Paul Thompson/Alex Wilson
# This one does flowchart for genetics study with both SCT and twins

library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(svglite)
library(rsvg)
library(png)
library(xlsx)
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
  dir<-"C:\\Users\\wilsona\\Dropbox\\project SCT analysis\\Data from redcap\\"}

#Sorry had to add this bit as the link didn't work for me. i.e. manual down load and redirect link.
#dir.PT<-"c:/Users/pthompson/Desktop/"

main.data <- read.csv(paste0(dir,"SCTData_DATA_2017-10-27_1035.csv"))

names(main.data)[1]<-"record_id"

#Referral source
nhscases<-c(70, 212, 213, 214, 215, 218, 221, 222, 223, 225, 229, 231,
            233, 234, 237, 239, 240, 241, 247, 251, 252, 254, 255, 257,
            258, 264, 265, 268, 270, 273, 274, 276, 283, 285, 297, 304,
            306, 307, 310, 318, 319, 321, 322, 323, 327, 328, 329, 331,
            334, 335, 336, 337, 340, 342, 346, 348, 351, 352, 355, 357,
            358, 360, 362, 365, 367, 368, 369, 370)
nhsset<-filter(main.data,record_id %in% nhscases)
othset<-setdiff(main.data,nhsset)

y<-table(nhsset$pre_postnatal_diag)
z<-table(othset$pre_postnatal_diag)
y1=y[1]
y2=y[2]
z1=z[1]
z2=z[2]



for (i in 1:1){ #change to 2 to see rates by reason for testing
  #Count for diagnosis pre or postnatal
  prenatals<-subset(main.data,pre_postnatal_diag==0)
  postnatals<-subset(main.data,pre_postnatal_diag==1)
  label1<-'Prenatal'
  label2<-'Postnatal'
  label3<-'When SCT identified'
  if (i==2){
    ascbias<-c(which(main.data$why_tested==2),which(main.data$why_tested==3))
    postnatals<-main.data[ascbias,] #these are no longer pre-post but divided
    #according to whether tested because of behav/neuro concerns
    prenatals<-main.data[-ascbias,] #tested bcs other med concerns
    label1<-'Medical'
    label2<-'Neurodev.'
    label3<-'Reason for testing'
  }
  n.A<-dim(postnatals)[1]
  n.B<-dim(prenatals)[1]
  
  #subset by karyotype
  xxx2<-subset(postnatals,trisomy==1)
  xxy2<-subset(postnatals,trisomy==2)
  xyy2<-subset(postnatals,trisomy==3)
  xxx1<-subset(prenatals,trisomy==1)
  xxy1<-subset(prenatals,trisomy==2)
  xyy1<-subset(prenatals,trisomy==3)
  
  n.C<-dim(xxx1)[1]
  n.D<-dim(xxy1)[1]
  n.E<-dim(xyy1)[1]
  n.F<-dim(xxx2)[1]
  n.G<-dim(xxy2)[1]
  n.H<-dim(xyy2)[1]
  
  #now check DNA status
  n.I=length(which(xxx1$dna_ok>0))
  n.J=length(which(xxy1$dna_ok>0))
  n.K=length(which(xyy1$dna_ok>0))
  n.L=length(which(xxx2$dna_ok>0))
  n.M=length(which(xxy2$dna_ok>0))
  n.N=length(which(xyy2$dna_ok>0))
  
  
  #now create flow chart ; TB denotes top to bottom
  #Need to add labels along the side: top row 'Reason for testing or Time of testing'
  #Then trisomy, then with DAWBA data
  
  print(grViz("
              digraph a_nice_graph {
              
              # node definitions with substituted label text
              node [shape = plaintext, fontname = Helvetica]
              X[label='@@1']
              Y[label= 'Trisomy']
              Z[label = 'With DNA']
              
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
              K[label='@@12']
              L[label='@@13']
              M[label='@@14']
              N[label='@@15']
              
              # edge definitions with the node IDs
              A -> {C D E}
              B -> {F G H}
              C -> I
              D -> J
              E -> K
              F -> L
              G -> M
              H -> N
              
              
              X -> Y [alpha=0,color='white']
              Y -> Z [alpha=0,color='white']
              }
              
              [1]: paste0(label3, ':\\n ',' Recruited from')
              [2]: paste0(label1,':\\n', 'NHS: N = ',y1,':\\n', 'Other: N = ',y2)
              [3]: paste0(label2,':\\n', 'NHS: N = ',z1,':\\n' ,'Other: N = ',z2)
              [4]: paste0('XXX',':\\n', 'N = ',n.C)
              [5]: paste0('XXY',' :\\n', 'N = ',n.D)
              [6]: paste0('XYY',' :\\n', 'N = ',n.E)
              [7]: paste0('XXX',' :\\n', 'N = ',n.F)
              [8]: paste0('XXY',' :\\n', 'N = ',n.G)
              [9]: paste0('XYY',' :\\n', 'N = ',n.H)
              [10]: paste0('N = ',n.I)
              [11]: paste0('N = ',n.J)
              [12]: paste0('N = ',n.K)
              [13]: paste0('N = ',n.L)
              [14]: paste0('N = ',n.M)
              [15]: paste0('N = ',n.N)
              "))
  
  
}
flow %>% export_svg %>% charToRaw %>% rsvg %>% png::writePNG("SCT_flow.png")

#Now do plot for twin data

dir<-"/Users/dorothybishop/Dropbox/ERCadvanced/project twin kids/Project_files/Data/"
if(os1=="windows"){
  dir<-"C:\\Users\\wilsona\\Dropbox\\project twin kids\\Project_files\\Data\\"}

#dir.PT<-"c:/Users/pthompson/Desktop/"

twin.data <- read.csv(paste0(dir,"TwinsData_DATA_2017-10-28_0620.csv"))
twin.data$zygo2<-twin.data$zygosity
fam_ids<-unique(twin.data$fam_id)
w<-which(twin.data$zygo2==3)
twin.data$zygo2[w]<-2
for (i in 1:length(fam_ids)){
twin.data$parent_rep_pair[which(twin.data$fam_id==fam_ids[i])]<-sum(twin.data$parental_report_dld[which(twin.data$fam_id==fam_ids[i])])
}

w<-which(twin.data$parent_rep_pair>0)
twin.data$parent_rep_pair[w]<-1
excl1<-nrow(filter(twin.data,parent_rep_pair==0,zygosity_di==9,randomtwininclude==1))
excl2<-nrow(filter(twin.data,parent_rep_pair==1,zygosity_di==9,randomtwininclude==1))
noconcern<-length(which(twin.data$parent_rep_pair==0))/2 #i.e. both twins have zero concern
concern<-length(which(twin.data$parent_rep_pair>0))/2 #i.e. one or both twins have concern

mytwins<-filter(twin.data,randomtwininclude==1,zygosity_di<9)#one per pair, exclude those with no DNA

mygirls<-filter(mytwins,female==1,zygo2<9)
myt<-table(mygirls$zygo2,mygirls$parent_rep_pair)
colnames(myt)<-c('none','concern')
rownames(myt)<-c('MZ','DZ')
n.MZn.g<-myt[1,1]
n.DZn.g<-myt[2,1]
n.MZc.g<-myt[1,2]
n.DZc.g<-myt[2,2]
myboys<-filter(mytwins,female==0,zygo2<9)
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
            Z[label= 'One twin per pair at random']
            
            
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
            [8]: paste0('N = ',excl1)
            [9]: paste0('N = ',excl2)
            
            "))

flow %>% export_svg %>% charToRaw %>% rsvg %>% png::writePNG("twins_flow.png")
