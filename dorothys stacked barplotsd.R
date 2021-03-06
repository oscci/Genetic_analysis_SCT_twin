#modified by DB 3rd Oct 2017
library(tidyverse)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(likert)

myuser='db' #select mac options

#mydir<-"C:/Users/pthompson/Dropbox/project SCT analysis/Data from redcap/"
mydir<-"~/Dropbox/ERCadvanced/project SCT analysis/Data from redcap/"
writedir<-mydir
sct.data<-read.csv(paste0(mydir,"SCTData_DATA_2017-10-03_1246.csv"))

short.sct2<-select(sct.data,record_id,age_at_test,trisomy,pre_postnatal_diag,wasi_matrices_ss,wasi_block_design_ss,wasi_vocab_ss,
                  wdck_jhsn_ss,sent_rep_ss,nonword_rep_ss,oromotor_ss,nara_acc_ss,nara_comp_ss,nara_rate_ss,
                  towre_words_ss,towre_nonwords_ss,phab_pic_ss,phab_digit_ss)

###recode data to categories
for(i in 5:length(short.sct2))
{
  short.sct2[,i]<-car::recode(short.sct2[,i],"996=NA;997=NA;998=NA;999=NA")
}


short.sct2$wasi_mat_cat<-cut(short.sct2$wasi_matrices_ss,breaks= c(0,30,40,74), labels=c(1:3))
short.sct2$wasi_block_cat<-cut(short.sct2$wasi_block_design_ss,breaks= c(0,30,40,max(short.sct2$wasi_block_design_ss,na.rm=T)), labels=c(1:3))
short.sct2$wasi_vocab_cat<-cut(short.sct2$wasi_vocab_ss,breaks= c(0,30,40,max(short.sct2$wasi_vocab_ss,na.rm=T)), labels=c(1:3))

short.sct2$WJ_cat<-cut(short.sct2$wdck_jhsn_ss,breaks= c(0,70,85,max(short.sct2$wdck_jhsn_ss,na.rm=T)), labels=c(1:3))

short.sct2$sent_rep_cat<-cut(short.sct2$sent_rep_ss,breaks= c(0,4,7,max(short.sct2$sent_rep_ss,na.rm=T)), labels=c(1:3))
short.sct2$nonword_rep_cat<-cut(short.sct2$nonword_rep_ss,breaks= c(0,4,7,max(short.sct2$nonword_rep_ss,na.rm=T)), labels=c(1:3))
short.sct2$oromotor_cat<-cut(short.sct2$oromotor_ss,breaks= c(0,4,7,max(short.sct2$oromotor_ss,na.rm=T)), labels=c(1:3))

short.sct2$nara_acc_cat<-cut(short.sct2$nara_acc_ss,breaks= c(0,70,85,max(short.sct2$nara_acc_ss,na.rm=T)), labels=c(1:3))
short.sct2$nara_comp_cat<-cut(short.sct2$nara_comp_ss,breaks= c(0,70,85,max(short.sct2$nara_comp_ss,na.rm=T)), labels=c(1:3))
short.sct2$nara_rate_cat<-cut(short.sct2$nara_rate_ss,breaks= c(0,70,85,max(short.sct2$nara_rate_ss,na.rm=T)), labels=c(1:3))

short.sct2$towre_words_cat<-cut(short.sct2$towre_words_ss,breaks= c(0,70,85,max(short.sct2$towre_words_ss,na.rm=T)), labels=c(1:3))
short.sct2$towre_nonwords_cat<-cut(short.sct2$towre_nonwords_ss,breaks= c(0,70,85,max(short.sct2$towre_nonwords_ss,na.rm=T)), labels=c(1:3))

short.sct2$phab_pic_cat<-cut(short.sct2$phab_pic_ss,breaks= c(0,70,85,max(short.sct2$phab_pic_ss,na.rm=T)), labels=c(1:3))
short.sct2$phab_digit_cat<-cut(short.sct2$phab_digit_ss,breaks= c(0,70,85,max(short.sct2$phab_digit_ss,na.rm=T)), labels=c(1:3))
names(short.sct2)[1]<-"record_id"

short.sct3<-select(short.sct2,record_id,age_at_test,trisomy,pre_postnatal_diag,wasi_mat_cat, wasi_block_cat, wasi_vocab_cat, WJ_cat, sent_rep_cat, nonword_rep_cat,     
oromotor_cat, nara_acc_cat, nara_comp_cat, nara_rate_cat, towre_words_cat, towre_nonwords_cat,
phab_pic_cat, phab_digit_cat) 

colnames(short.sct3)<-c('ID','age','trisomy','pre_postnatal_diag','Matrices','Blocks','Vocab','Sent. comp',
                        'Sent Rep','Nonword rep','Oromotor','Read acc',
                        'Read comp','Read rate','Word read','Nonword read','Pic naming','Digit naming')
####plot likert bar
short.sct3$trisomy[short.sct3$trisomy==9]<-NA
short.sct3$trisomy<-as.factor(short.sct3$trisomy)
short.sct3$pre_postnatal_diag<-as.factor(short.sct3$pre_postnatal_diag)

#all

my.lik<-likert(short.sct3[5:18])
plot(my.lik)



#Pre-natal
# windows(width=22,height=8)
# if(myuser=='db'){
#   quartz(width=22,height=8)
#}
  pdf("prenatal_plot.pdf",width=16,height=8)


items_XXX<-short.sct3[short.sct3$trisomy==1&short.sct3$pre_postnatal_diag==0,5:18]
items_XXY<-short.sct3[short.sct3$trisomy==2&short.sct3$pre_postnatal_diag==0,5:18]
items_XYY<-short.sct3[short.sct3$trisomy==3&short.sct3$pre_postnatal_diag==0,5:18]



my.lik_XXX<-likert(items_XXX)
p1<-plot(my.lik_XXX,wrap=30,centered=F,group.order=my.lik$results[order(my.lik$results[,4],decreasing=T),1])+ ggtitle('XXX') + theme(legend.position="none") + scale_fill_manual(values = c("red", "purple", "blue"))
#
my.lik_XXY<-likert(items_XXY)
p2<-plot(my.lik_XXY,wrap=30,centered=F,group.order=my.lik$results[order(my.lik$results[,4],decreasing=T),1])+ ggtitle('XXY') + theme(legend.position="none")+ scale_fill_manual(values = c("red", "purple", "blue"))
#
my.lik_XYY<-likert(items_XYY)
p3<-plot(my.lik_XYY,wrap=30,centered=F,group.order=my.lik$results[order(my.lik$results[,4],decreasing=T),1])+ ggtitle('XYY') + theme(legend.position="none")+ scale_fill_manual(values = c("red", "purple", "blue"))



grid.arrange(p1, p2, p3, ncol=3, top = "Pre-natal")
dev.off()
##########################################
pdf("postnatal_plot.pdf",width=16,height=8)
#Post-natal
items_XXX<-short.sct3[short.sct3$trisomy==1&short.sct3$pre_postnatal_diag==1,5:18]
items_XXY<-short.sct3[short.sct3$trisomy==2&short.sct3$pre_postnatal_diag==1,5:18]
items_XYY<-short.sct3[short.sct3$trisomy==3&short.sct3$pre_postnatal_diag==1,5:18]

#
my.lik_XXX<-likert(items_XXX)
p4<-plot(my.lik_XXX,wrap=30,centered=F,group.order=my.lik$results[order(my.lik$results[,4],decreasing=T),1])+ ggtitle('XXX') + theme(legend.position="none")+ scale_fill_manual(values = c("red", "purple", "blue"))
#
my.lik_XXY<-likert(items_XXY)
p5<-plot(my.lik_XXY,wrap=30,centered=F,group.order=my.lik$results[order(my.lik$results[,4],decreasing=T),1])+ ggtitle('XXY') + theme(legend.position="none")+ scale_fill_manual(values = c("red", "purple", "blue"))
#
my.lik_XYY<-likert(items_XYY)
p6<-plot(my.lik_XYY,wrap=30,centered=F,group.order=my.lik$results[order(my.lik$results[,4],decreasing=T),1])+ ggtitle('XYY') + theme(legend.position="none")+ scale_fill_manual(values = c("red", "purple", "blue"))


grid.arrange(p4, p5, p6, ncol=3,top="Post-natal")

dev.off()

pdf("prenatal_lang_plot.pdf",width=11,height=4)

items_XXX<-short.sct3[short.sct3$trisomy==1&short.sct3$pre_postnatal_diag==0,5:11]
items_XXY<-short.sct3[short.sct3$trisomy==2&short.sct3$pre_postnatal_diag==0,5:11]
items_XYY<-short.sct3[short.sct3$trisomy==3&short.sct3$pre_postnatal_diag==0,5:11]

my.lik_XXX<-likert(items_XXX)
p1<-plot(my.lik_XXX,wrap=30,centered=F,
         plot.percent.low=FALSE,plot.percent.neutral=FALSE,plot.percent.high=FALSE ,
         group.order=my.lik$results[order(my.lik$results[1:7,4],decreasing=T),1])+
  ggtitle('XXX') + theme(legend.position="none") + 
  scale_fill_manual(values = c("red", "purple", "blue"))
                                                       
#
my.lik_XXY<-likert(items_XXY)
p2<-plot(my.lik_XXY,wrap=30,centered=F,
         plot.percent.low=FALSE,plot.percent.neutral=FALSE,plot.percent.high=FALSE ,
         group.order=my.lik$results[order(my.lik$results[1:7,4],decreasing=T),1])+
  ggtitle('XXY') + theme(legend.position="none") + 
  scale_fill_manual(values = c("red", "purple", "blue"))
#
my.lik_XYY<-likert(items_XYY)
p3<-plot(my.lik_XYY,wrap=30,centered=F,
         plot.percent.low=FALSE,plot.percent.neutral=FALSE,plot.percent.high=FALSE ,
         group.order=my.lik$results[order(my.lik$results[1:7,4],decreasing=T),1])+
  ggtitle('XYY') + theme(legend.position="none") + 
  scale_fill_manual(values = c("red", "purple", "blue"))
#
grid.arrange(p1, p2, p3, ncol=3, top = "Pre-natal")
dev.off()
##########################################
pdf("postnatal_lang_plot.pdf",width=11,height=4)
#Post-natal
items_XXX<-short.sct3[short.sct3$trisomy==1&short.sct3$pre_postnatal_diag==1,5:11]
items_XXY<-short.sct3[short.sct3$trisomy==2&short.sct3$pre_postnatal_diag==1,5:11]
items_XYY<-short.sct3[short.sct3$trisomy==3&short.sct3$pre_postnatal_diag==1,5:11]

#
my.lik_XXX<-likert(items_XXX)
p4<-plot(my.lik_XXX,wrap=30,centered=F,group.order=my.lik$results[order(my.lik$results[1:7,4],decreasing=T),1])+ ggtitle('XXX') + theme(legend.position="none")+ scale_fill_manual(values = c("red", "purple", "blue"))
#
my.lik_XXY<-likert(items_XXY)
p5<-plot(my.lik_XXY,wrap=30,centered=F,group.order=my.lik$results[order(my.lik$results[1:7,4],decreasing=T),1])+ ggtitle('XXY') + theme(legend.position="none")+ scale_fill_manual(values = c("red", "purple", "blue"))
#
my.lik_XYY<-likert(items_XYY)
p6<-plot(my.lik_XYY,wrap=30,centered=F,group.order=my.lik$results[order(my.lik$results[1:7,4],decreasing=T),1])+ ggtitle('XYY') + theme(legend.position="none")+ scale_fill_manual(values = c("red", "purple", "blue"))


grid.arrange(p4, p5, p6, ncol=3,top="Post-natal")

dev.off()