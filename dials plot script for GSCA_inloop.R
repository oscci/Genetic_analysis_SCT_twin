#Script modified by DVMB 21/9/18
#Created a loop for all plots
#Plots saved as RDS objects so can be reloaded and assembled into one big plot
#Increased font size for power
#Changed power labels to proportions rather than %
#Changed layout of big plots and removed 1 gene versions from the clustered plot (so it is 2 x 2)
#Modified some other labels

# library
library(tidyverse)
library(viridis)
library(ggplot2)
library(grid)
library(gridExtra)


mydir<-"~/Dropbox/ERCadvanced/project SCT analysis/GSCA validation/simulations_results2/"
# 1 gene, cor=0.4
plot.data<-read.csv(paste0(mydir,"PT_simulation_results_MC_adj_formatted.csv"))
thisplot<-0
for (c in 0:1){
  for(e in c(4,1)){
   for (g in c(1,2,4)){
     
      thisplot<-thisplot+1
      plothead<-paste0(LETTERS[thisplot],': ',g,' gene; r = ',e/10)
      data<-plot.data[plot.data$ngenes==g &plot.data$eff_size==(e/10)&plot.data$cluster==c,]
      data<-tidyr::gather(data,"gene_ind","nsig",nsig1:nsig4)
      
      data<-data.frame(individual.lab=data$Npats,group=data$nseff,observation=data$nsnp, value=data$nsig,gene_ind=data$gene_ind)
      
      data$individual<-paste0(data$individual.lab,":",data$group,":",data$gene_ind,":",data$observation)
      
      w<-which(is.na(data$value))
      if(length(w)>0){
      data<-data[-w,]
      }
      # data<-data.frame(individual.lab=data$Npats,group=data$nseff,observation=data$nsnp, value=data$nsig1)
      # data$individual<-paste0(data$individual.lab,":",data$group,":",data$observation)
      # 

# Transform data in a tidy format (long format)
#data = data %>% gather(key = "observation", value="value", -c(1,2))

data$group<-as.factor(data$group)
data$observation<-as.factor(data$observation)
data$individual<-as.factor(data$individual)
data$group<-as.factor(data$group)

# Set a number of 'empty bar' to add at the end of each group
empty_bar=2
nObsType=nlevels(as.factor(data$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar*nObsType )
data=rbind(data, to_add)
data=data %>% arrange(group, individual)
data$id=rep( seq(1, nrow(data)/nObsType) , each=nObsType)
data$power<-data$value/100
data$power[data$value<10]<-NA #use proportions rather than % for bar labels and remove labels for low numbers
#Plot easier to understand if scale for N participants, N Snps and power are not all integers!

# Get the name and the y position of each label
label_data= data %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
#label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$hjust <-1
label_data$angle<-ifelse(angle < -90, angle+180, angle)

label_data$individual <- substring(label_data$individual,1,3)
label_data$tot<-label_data$tot/100
label_data$tot[label_data$tot<.10]<-NA



# prepare a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
p<-ggplot(data) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
  geom_text(aes(x=as.factor(id), y=value, fill=observation,label=power),size = 4, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete=TRUE) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 300, xend = start, yend = 300), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),4), y = c(0, 100, 200, 300), label = c("0", "100", "200", "300") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-160,330) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=330, label=individual, hjust=hjust), color="black", 
            fontface="bold",alpha=0.6, size=6, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", 
               alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0), colour = "black", 
            alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE) +
  annotate("text", x = 0, y = -160, label = 'atop(bold("\n  N SNPs\nwith effect"))', parse = TRUE,size=5)+
  annotate("text", x = 0, y = 320, label = plothead, size=6)+
  guides(fill=guide_legend(title="N SNPs analysed"))+theme(text = element_text(size=25))
p #need to explicitly plot - ggsave will save most recent plot
plotname<-paste0('p',thisplot,'.png')
ggsave(plotname) #save as png
savename<-paste0('p',thisplot)
saveRDS(p,savename) #save as RDS object - can then read it in again
    }
  }
}

p1<-readRDS('p1')
p2<-readRDS('p2')
p3<-readRDS('p3')
p4<-readRDS('p4')
p5<-readRDS('p5')
p6<-readRDS('p6')

png(filename="dials_nocluster_MC.png",height=20,width=20,units="in",res=100)
#windows()
big_p<-ggpubr::ggarrange(p1,p4,p2,p5,p3,p6,ncol=2,nrow=3,
                         common.legend = TRUE, legend = "bottom",hjust=-0.6,font.label = list(size = 25, face = "bold"))


library(grid)
print(big_p, vp=viewport(angle=0))
dev.off()

#file.show("dials_nocluster_MC.png")

#NB files p7 and p10 are identical to cluster condition because just one gene
p1<-readRDS('p8')
p2<-readRDS('p9')
p3<-readRDS('p11')
p4<-readRDS('p12')

png(filename="dials_cluster_MC.png",height=13.4,width=20,units="in",res=100)
#windows()
big_p<-ggpubr::ggarrange(p1,p3,p2,p4,ncol=2,nrow=2,
                         common.legend = TRUE, legend = "bottom",hjust=-0.6,font.label = list(size = 25, face = "bold"))


library(grid)
print(big_p, vp=viewport(angle=0))
dev.off()

#
file.show("dials_nocluster_MC.png")