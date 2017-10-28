#Produces flow chart showing group membership and exclusions
#Original version September 2017 for PeerJ submission of twin laterality paper.

#Updated 28th Oct 2017 to incorporate updated twin file with random twin selection

#install.packages('DiagrammeR')
library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(svglite)
library(rsvg)
library(png)
library(xlsx)
library(tidyverse)
library(stringr)

#dir<-"C:\\Users\\Alex\\Dropbox\\project twin kids\\"
dir<-"/Users/dorothybishop/Dropbox/ERCadvanced/project twin kids/Project_files/Data/"

#main.data <- read.csv(paste0(dir,"TwinsData_DATA_2017-08-19_1139.csv"))
main.data <- read.csv(paste0(dir,"TwinsData_DATA_2017-10-28_0620.csv")) #NB latest version
#this version has additional columns, so need to make column references dynamic later on

my.ncol<-ncol(main.data)

#record exclusions for reporting in the paper
n.ASD<-sum(main.data$lang_probs==2)
n.ID<-sum(main.data$lang_probs==3)
n.hearing<-sum(main.data$lang_probs==4)
n.other<-sum(main.data$include==0)-(n.ASD+n.ID+n.hearing)


#main.data codes whether parent reports their child as having language problems (1) or not (0). Now need to code whether twin pair is concordant
#for DLD (2), discordant for DLD (1), or neither child has DLD (0)
fam_ids<-unique(main.data$fam_id)
main.data<-cbind(main.data,matrix(0,ncol=2,nrow=length(main.data$record_id)))#add 2 new columns
colnames(main.data)[(my.ncol+1):(my.ncol+2)]<-c("parent_rep_pair","DLD_pair")
for(i in 1:length(fam_ids)){
  mytwins_DLDstatus<-main.data$lang_probs[which(main.data$fam_id==fam_ids[i])]
  if(mytwins_DLDstatus[1]==0 & mytwins_DLDstatus[2]==0){
    main.data$DLD_pair[which(main.data$fam_id==fam_ids[i])]<-0
  }
  if(sum(mytwins_DLDstatus==1)){
    main.data$DLD_pair[which(main.data$fam_id==fam_ids[i])]<-1
  }
  if(mytwins_DLDstatus[1]==1 & mytwins_DLDstatus[2]==1){
    main.data$DLD_pair[which(main.data$fam_id==fam_ids[i])]<-2
  }
  if(sum(main.data$include[which(main.data$fam_id==fam_ids[i])])==1){
    includedtwin<-which(main.data$include[which(main.data$fam_id==fam_ids[i])]==1)
    includedindex<-which(main.data$fam_id==fam_ids[i])[includedtwin]
    excludedtwin<-which(main.data$include[which(main.data$fam_id==fam_ids[i])]==0)
    excludedindex<-which(main.data$fam_id==fam_ids[i])[excludedtwin]
    main.data$DLD_pair[includedindex]<-9
    main.data$DLD_pair[excludedindex]<-99
  }
  if(sum(main.data$include[which(main.data$fam_id==fam_ids[i])])==0){
    main.data$DLD_pair[which(main.data$fam_id==fam_ids[i])]<-999
  }
  main.data$parent_rep_pair[which(main.data$fam_id==fam_ids[i])]<-sum(main.data$parental_report_dld[which(main.data$fam_id==fam_ids[i])])
  #main.data$DLD_pair[which(main.data$fam_id==fam_ids[i])]<-sum(main.data$lang_probs[which(main.data$fam_id==fam_ids[i])])
}
table(main.data$parent_rep_pair)
#subset by parental report of DLD or no language problems reported
main.bothDLD<-subset(main.data,parent_rep_pair==2)
main.oneDLD<-subset(main.data,parent_rep_pair==1)
main.noDLD<-subset(main.data,parent_rep_pair==0)

#first level - parental report
N1<-dim(main.bothDLD)[[1]] #report DLD both
N2<-dim(main.oneDLD)[[1]] #report DLD one
N3<-dim(main.noDLD)[[1]] #no report of DLD

#second level - language problems identified on testing
#(1) in pairs where both twins are reported as having current/historic difficulties
N4a<-dim(subset(main.bothDLD,DLD_pair==2 & include==1))[[1]]
n4b<-dim(subset(main.bothDLD,DLD_pair==1 & include==1))[[1]]
n4c<-dim(subset(main.bothDLD,DLD_pair==0 & include==1))[[1]]
n4d<-dim(subset(main.bothDLD,DLD_pair==999 & include==0))[[1]] #indicates both excluded
n4d1<-dim(subset(main.bothDLD,DLD_pair==99 & include==0))[[1]] #indicates one excluded
n4b1<-dim(subset(main.bothDLD,DLD_pair==9 & include==1 & lang_probs == 1))[[1]] #indicates DLD where fellow twin excluded
n4c1<-dim(subset(main.bothDLD,DLD_pair==9 & include==1 & lang_probs == 0))[[1]] #indicates non-DLD where fellow twin excluded

N4b<-n4b+n4b1; N4c<-n4c+n4c1; N4d<-n4d+n4d1

#(2) in pairs where one twin is reported as having current/historic difficulties
N5a<-dim(subset(main.oneDLD,DLD_pair==2 & include==1))[[1]]
n5b<-dim(subset(main.oneDLD,DLD_pair==1 & include==1))[[1]]
n5c<-dim(subset(main.oneDLD,DLD_pair==0 & include==1))[[1]]
n5d<-dim(subset(main.oneDLD,DLD_pair==999 & include==0))[[1]] #indicates both excluded
n5d1<-dim(subset(main.oneDLD,DLD_pair==99 & include==0))[[1]] #indicates one excluded
n5b1<-dim(subset(main.oneDLD,DLD_pair==9 & include==1 & lang_probs == 1))[[1]] #indicates DLD where fellow twin excluded
n5c1<-dim(subset(main.oneDLD,DLD_pair==9 & include==1 & lang_probs == 0))[[1]] #indicates non-DLD where fellow twin excluded

N5b<-n5b+n5b1; N5c<-n5c+n5c1; N5d<-n5d+n5d1

#(2) in pairs where neither twin is reported as having current/historic difficulties
N6a<-dim(subset(main.noDLD,DLD_pair==2 & include==1))[[1]]
n6b<-dim(subset(main.noDLD,DLD_pair==1 & include==1))[[1]]
n6c<-dim(subset(main.noDLD,DLD_pair==0 & include==1))[[1]]
n6d<-dim(subset(main.noDLD,DLD_pair==999 & include==0))[[1]] #indicates both excluded
n6d1<-dim(subset(main.noDLD,DLD_pair==99 & include==0))[[1]] #indicates one excluded 
n6b1<-dim(subset(main.noDLD,DLD_pair==9 & include==1 & lang_probs == 1))[[1]] #indicates DLD where fellow twin excluded
n6c1<-dim(subset(main.noDLD,DLD_pair==9 & include==1 & lang_probs == 0))[[1]] #indicates non-DLD where fellow twin excluded

N6b<-n6b+n6b1; N6c<-n6c+n6c1; N6d<-n6d+n6d1

#third level - diagnosis
N7<-dim(subset(main.data,lang_probs==1 & include==1))[[1]]
N8<-dim(subset(main.data,lang_probs==0 & include==1))[[1]]
N9<-dim(subset(main.data,include==0))[[1]]

#filter for for useable fTCD data
main.useableTCD<-main.data[which(!is.na(main.data$laterality_index)),]
n.extreme.li<-length(which(abs(main.useableTCD$laterality_index)>10))
main.useableTCD<-filter(main.useableTCD,abs(main.useableTCD$laterality_index)<10 & n_trials> 11)

#fourth level - usable doppler data
N7a<-dim(subset(main.useableTCD,lang_probs==1 & include==1))[[1]]
N8a<-dim(subset(main.useableTCD,lang_probs==0 & include==1))[[1]]

#captured<-capture.output(
flow<-grViz("
  digraph a_nice_graph {
            
            # node definitions with substituted label text
            node [shape = plaintext, fontname = Helvetica, fontsize = 30]
            W[label='Parental Report']
            X[label='Cognitive Testing\nand Exclusions']
            Y[label='Diagnosis of DLD']
            Z[label='Useable fTCD data']
            
            node [shape=sqaure,margin=0.5,fontsize=30]
            A[label='@@1']
            B[label='@@2']
            C[label='@@3']
            
            D[label='@@4']
            E[label='@@5']
            F[label='@@6']
            G[label='@@7']
            
            H[label='@@8']
            I[label='@@9']

            J[label='@@10']
            K[label='@@11']
            
            L[label='@@12']
            M[label='@@13']
            
            # edge definitions with the node IDs
            edge [minlen=3,headport=n]
            A -> {D E} [arrowhead = vee]
            B -> {F G} [arrowhead = vee]
            C -> {H I} [arrowhead = vee]
            D -> {J K} [arrowhead = vee]
            F -> {J K} [arrowhead = vee]
            H -> {J K} [arrowhead = vee]

            J -> L
            K -> M
            
            W -> X [alpha=0,color='white']
            X -> Y [alpha=0,color='white']
            Y -> Z [alpha=0,color='white']
  }

  [1]: paste0('Language problems\\nBoth twins, N = ',N1)
  [2]: paste0('Language problems\\nOne twin, N = ',N2)
  [3]: paste0('No problems  \\nreported, N = ',N3)
  [4]: paste0('Both DLD, N = ',N4a,'\\nOne DLD, N = ',N4b,'\\nTD, N = ',N4c)
  [5]: paste0('Exclusions,\\nN = ',N4d)
  [6]: paste0('Both DLD, N = ',N5a,'\\nOne DLD, N = ',N5b,'\\nTD, N = ',N5c)
  [7]: paste0('Exclusions,\\nN = ',N5d)
  [8]: paste0('Both DLD, N = ',N6a,'\\nOne DLD, N = ',N6b,'\\nTD, N = ',N6c)
  [9]: paste0('Exclusions,\\nN = ',N6d)
  [10]: paste0('DLD,\\nN = ',N7)
  [11]: paste0('TD,\\nN = ',N8)
  [13]: paste0('N = ',N7a)
  [14]: paste0('N = ',N8a)
  ")
#)
flow %>% export_svg %>% charToRaw %>% rsvg %>% png::writePNG("flow.png")
