
####
#### Final Figure 2 plots
####

###Also includes code for supplementary figures S1, S2, S3, and S4

###Panels 1A-1B: PCA
PCA_samples<-allsamples_PCAinfo %>% filter(!is.na(PC1),total>=50000,!is.na(plot),isMom) %>%
  mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2","LM")),
         Year=factor(Year,levels=c(2019,2021,2022))) %>% slice_sample(n=nrow(.),replace = F)

PC1v2_legend<-PCA_samples %>% 
  ggplot(aes(x=PC1,y=PC2,color=Year,shape=stream)) + geom_point() + theme_bw() + 
  xlab(paste0("PC1 (",round(allnorthmoms_eigen$values[1]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  ylab(paste0("PC2 (",round(allnorthmoms_eigen$values[2]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  scale_color_manual(values=colors_yr[c(2,3,4)]) + 
  scale_shape_manual(values=c(1,2,0),labels=c("CAC_S1","CAC_S2","LM"))
PC1v2<-PCA_samples %>%
  ggplot(aes(x=PC1,y=PC2,color=Year,shape=stream)) + geom_point() + theme_bw() + 
  xlab(paste0("PC1 (",round(allnorthmoms_eigen$values[1]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  ylab(paste0("PC2 (",round(allnorthmoms_eigen$values[2]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  scale_color_manual(values=colors_yr[c(2,3,4)],guide="none") + 
  scale_shape_manual(values=c(1,2,0),guide="none")
PC1v2_legend
PC1v2

PC3v4_legend<-PCA_samples %>% filter(!(PCAcluster %in% c("PC1_nasutus","PC2_LM"))) %>%
  ggplot(aes(x=PC3,y=PC4,color=Year,shape=stream)) + geom_point() + theme_bw() + 
  xlab(paste0("PC3 (",round(allnorthmoms_eigen$values[3]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  ylab(paste0("PC4 (",round(allnorthmoms_eigen$values[4]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  scale_color_manual(values=colors_yr[c(2,3,4)]) + 
  scale_shape_manual(values=c(1,2,0),labels=c("CAC_S1","CAC_S2","LM"))

PC3v4<-PCA_samples %>% filter(!(PCAcluster %in% c("PC1_nasutus","PC2_LM"))) %>%
  ggplot(aes(x=PC3,y=PC4,color=Year,shape=stream)) + geom_point() + theme_bw() + 
  xlab(paste0("PC3 (",round(allnorthmoms_eigen$values[3]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  ylab(paste0("PC4 (",round(allnorthmoms_eigen$values[4]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  scale_color_manual(values=colors_yr[c(2,3,4)],guide="none") + 
  scale_shape_manual(values=c(1,2,0),guide="none")
PC3v4_legend
PC3v4

bothPCs<-plot_grid(PC1v2,PC3v4,nrow=2,align="hv",axis="tblr")
bothPCs
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/SEP_versions/PCs_1to4.png",bothPCs,device = 'png',height = 6,width=4,dpi=600,units='in')

###Panel 1C: PC1 vs. ancestry
ancestry_samples<-allsamples_PCAinfo %>% filter(!is.na(K1_prop),!is.na(plot),total>=50000)

PC1_ancestry_lm<-lm(PC1~nas_prop,data=ancestry_samples)

PC1_vs_ancestry<-ggplot(ancestry_samples,aes(x=nas_prop,y=PC1)) + geom_point() + 
  theme_bw() + xlab("Hybrid index") + 
  ylab(paste0("PC1 (",round(allnorthmoms_eigen$values[1]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  geom_abline(slope=PC1_ancestry_lm$coefficients[2],intercept=PC1_ancestry_lm$coefficients[1],color="red") + 
  geom_text(x=0.75,y=0,label=paste0("R^2=",round(summary(PC1_ancestry_lm)$r.squared,3)))

PC1_vs_ancestry
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/May_versions/PC1_vs_ancestry.png",PC1_vs_ancestry,device = 'png',height = 3,width=3,dpi=600,units='in')

###Panel 1D: Histogram of ancestry
histogram_samples<-allsamples_PCAinfo %>% filter(isMom,total>=50000,!is.na(plot),isSOOK %in% c("no","F1")) %>%
  mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2","LM")),
         Year=factor(Year,levels=c(2012,2019,2021,2022)))

streamlabels<-c("CAC_S1","CAC_S2","LM")
names(streamlabels)<-c("CAC_stream1","CAC_stream2","LM")

histograms_legend<-ggplot(histogram_samples,aes(x=nas_prop,fill=Year)) + 
  facet_grid(stream~.,scales='free_y',space='free_y',labeller=labeller(stream=streamlabels)) + 
  geom_histogram(binwidth=0.05,boundary=1,position="stack",color="black",linewidth=0.2) + 
  #geom_histogram(binwidth=0.05,boundary=1,position="identity",fill=NA,color="black") + 
  theme_bw() + 
  scale_fill_manual(values=colors_yr) + 
  scale_y_continuous(breaks=c(0:6)*10) + 
  xlab("Hybrid index") + ylab("Count of maternal samples")
histograms_legend

histograms<-ggplot(histogram_samples,aes(x=nas_prop,fill=Year)) + 
  facet_grid(stream~.,scales='free_y',space='free_y',labeller=labeller(stream=streamlabels)) + 
  geom_histogram(binwidth=0.05,boundary=1,position="stack",color="black",linewidth=0.2) + 
  #geom_histogram(binwidth=0.05,boundary=1,position="identity",fill=NA,color="black") + 
  theme_bw() + 
  scale_fill_manual(values=colors_yr,guide="none") + 
  scale_y_continuous(breaks=c(0:6)*10) + 
  xlab("Hybrid index") + ylab("Count of maternal samples")
histograms

###Putting it all together
plot_grid(PC1v2,PC1_vs_structure,PC3v4,histograms,nrow=2,align="tblr",axis="hv")
Rcomplete_Fig2<-plot_grid(
  plot_grid(PC1v2,PC3v4,nrow=2,align="tblr",axis="hv",rel_heights=c(1,1)),
  plot_grid(PC1_vs_structure,histograms,nrow=2,align="v",axis='lr',rel_heights=c(1,2)),
  nrow=1
)
Rcomplete_Fig2_legend<-plot_grid(
  plot_grid(PC1v2_legend,PC3v4,nrow=2,align="tblr",axis="hv",rel_heights=c(1,1)),
  plot_grid(PC1_vs_structure,histograms_legend,nrow=2,align="v",axis='lr',rel_heights=c(1,2)),
  nrow=1
)
Rcomplete_Fig2_legend
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/SEP_versions/Rcomplete_Fig2.pdf",Rcomplete_Fig2,device = 'pdf',height = 6,width=6,dpi=600,units='in')
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/SEP_versions/Rcomplete_Fig2_legend.pdf",Rcomplete_Fig2_legend,device = 'pdf',height = 6,width=6,dpi=600,units='in')


### Figure S1: ancestry vs. structure
structure_ancestry_lm<-lm(K1_prop~nas_prop,data=ancestry_samples)

ancestry_vs_structure<-ggplot(ancestry_samples,aes(x=nas_prop,y=K1_prop)) + geom_point() + 
  theme_bw() + xlab("Hybrid index (ancestry_hmm)") + ylab("Admixture proportion (NGSAdmix, K=2)") + 
  geom_abline(slope=structure_ancestry_lm$coefficients[2],intercept=structure_ancestry_lm$coefficients[1],color="red") + 
  geom_text(x=0.75,y=0.25,label=paste0("R^2=",round(summary(structure_ancestry_lm)$r.squared,3)))

ancestry_vs_structure
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/SEP_versions/ancestry_vs_structure.pdf",ancestry_vs_structure,device = 'pdf',height = 3.1,width=3,dpi=600,units='in')

### Figure S2: additional PCA analysis on only CAC guttatus+admixed individuals
infolder<-"~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Outputs/angsd_Dec2024"

CACadmixedonlyPCAmoms_cov<-read.table(file.path(infolder,"CACadmixedonlyPCA.IM62_v3.CACpanel100.maf20.nonmiss60.cov"),header=F)
CACadmixedonlyPCAmoms_eigen<-eigen(CACadmixedonlyPCAmoms_cov)
CACadmixedonlyPCAmoms_PCAvariance<-CACadmixedonlyPCAmoms_eigen$values/sum(CACadmixedonlyPCAmoms_eigen$values)

CACadmixedonly_newangsd<-read.table(file.path(infolder,"moms_to_PCA_CACadmixedonly_bamlist.txt"),header=F,col.names="bampath") %>%
  bind_cols(read.table(file.path(infolder,"CACadmixedonlyPCA.IM62_v3.panel38_thinned1kb.K2.qopt"),header=F,col.names=c("Structure.Q1","Structure.Q2"))) %>%
  bind_cols(as.data.frame(CACadmixedonlyPCAmoms_eigen$vectors[,1:5])) %>%
  left_join(select(allsamples_sookinfo,c(sampleID,bampath,nas_prop,heterozygosity,Year,stream,plot,plot_renamed,isSOOK)),by="bampath")

library(RColorBrewer)
library(pals)
RdYlBu<-brewer.pal(11,name="RdYlBu")
PRGn<-brewer.pal(11,name="PRGn")
RdGy<-brewer.pal(11,name="RdGy")
colors_plotpalette=c(RdGy[1],RdYlBu[c(1,3,5)],PRGn[c(5,3,2)],PRGn[c(11,9,7)],RdYlBu[c(11,9,8)])

doesItCorrelate<-CACadmixedonly_newangsd %>% mutate(Year=factor(Year)) %>%
  ggplot(aes(x=V1,y=nas_prop,color=plot_renamed,shape=stream)) + 
  geom_point() + theme_bw() + 
  ylab(paste0("Hybrid index")) + 
  xlab(paste0("PC1 (",round(CACadmixedonlyPCAmoms_PCAvariance[1],4)*100,"%)")) + 
  scale_color_manual(values=colors_plotpalette,name="Plot",
                     guide=guide_legend(override.aes=list(shape=c(rep(16,7),rep(17,6))))) + 
  scale_shape_manual(values=c(16,17),guide="none")

colored_newPCA<-CACadmixedonly_newangsd %>% mutate(Year=factor(Year)) %>%
  ggplot(aes(x=V1,y=V2,color=plot_renamed,shape=stream)) + 
  geom_point() + theme_bw() + 
  xlab(paste0("PC1 (",round(CACadmixedonlyPCAmoms_PCAvariance[1],4)*100,"%)")) + 
  ylab(paste0("PC2 (",round(CACadmixedonlyPCAmoms_PCAvariance[2],4)*100,"%)")) + 
  scale_color_manual(values=colors_plotpalette,name="Plot",
                     guide=guide_legend(override.aes=list(shape=c(rep(16,7),rep(17,6))))) + 
  scale_shape_manual(values=c(16,17),guide="none")

colored_newPCA_v34<-CACadmixedonly_newangsd %>% mutate(Year=factor(Year)) %>%
  ggplot(aes(x=V3,y=V4,color=plot_renamed,shape=stream)) + 
  geom_point() + theme_bw() + 
  xlab(paste0("PC3 (",round(CACadmixedonlyPCAmoms_PCAvariance[3],3)*100,"%)")) + 
  ylab(paste0("PC4 (",round(CACadmixedonlyPCAmoms_PCAvariance[4],4)*100,"%)")) + 
  scale_color_manual(values=colors_plotpalette,name="Plot",
                     guide=guide_legend(override.aes=list(shape=c(rep(16,7),rep(17,6))))) + 
  scale_shape_manual(values=c(16,17,0),guide="none")

newSuppPCA<-doesItCorrelate+colored_newPCA+colored_newPCA_v34+plot_layout(design="A\nB\nC",guides = "collect")
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/DEC_revision/CAConlyPCA_byplot.pdf",newSuppPCA,device="pdf",units = "in",dpi=600,width=3,height=6)

### Figure S3: triangle plots
triangle_samples<-allsamples_PCAinfo %>% filter(isMom,total>=50000,!is.na(plot),Year!=2012) %>% 
  mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2","LM")),
         Year=factor(Year,levels=c(2019,2021,2022))) %>% slice_sample(n=nrow(.),replace=F)

triangle<-ggplot(triangle_samples,aes(x=nas_prop,y=heterozygosity,color=Year,shape=stream)) + 
  geom_point() + theme_bw() + ylim(0,1) + 
  scale_color_manual(values=colors_yr[c(2,3,4)]) + 
  scale_shape_manual(values=c(16,17,15),breaks=c("CAC_stream1","CAC_stream2","LM"),labels=c("CAC_S1","CAC_S2","LM")) + 
  geom_segment(x=0,y=0,xend=0.5,yend=1,color="black") + geom_segment(x=0.5,y=1,xend=1,yend=0,color="black") + 
  xlab("Hybrid index (ancestry_hmm)") + ylab("Ancestry heterozygosity")
triangle

triangle_offspring<-allsamples_PCAinfo %>% filter(!isMom,total>=50000,!is.na(plot),Year!=2012) %>% 
  mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2","LM")),
         Year=factor(Year,levels=c(2019,2021,2022))) %>% filter(!is.na(stream)) %>% slice_sample(n=nrow(.),replace=F)
triangle_offs<-ggplot(triangle_offspring,aes(x=nas_prop,y=heterozygosity,color=Year,shape=stream)) + 
  geom_point() + theme_bw() + ylim(0,1) + 
  scale_color_manual(values=colors_yr[c(2,3,4)]) + 
  geom_segment(x=0,y=0,xend=0.5,yend=1,color="black") + geom_segment(x=0.5,y=1,xend=1,yend=0,color="black") + 
  scale_shape_manual(values=c(16,17,15),breaks=c("CAC_stream1","CAC_stream2","LM"),labels=c("CAC_S1","CAC_S2","LM")) + 
  xlab("Hybrid index (ancestry_hmm)") + ylab("Ancestry heterozygosity")
triangle_offs

two_triangles<-plot_grid(triangle,triangle_offs,align="hv",axis='tblr')
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/SEP_versions/two_triangles.pdf",two_triangles,device = 'pdf',height = 4,width=10,dpi=600,units='in')

###Figure S4: PC5 sookensis
PC5<-PCA_samples %>%
  ggplot(aes(x=PC1,y=PC5,color=Year,shape=stream)) + geom_point() + theme_bw() + 
  xlab(paste0("PC1 (",round(allnorthmoms_eigen$values[1]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  ylab(paste0("PC5 (",round(allnorthmoms_eigen$values[5]/sum(allnorthmoms_eigen$values),4)*100,"%)")) + 
  scale_color_manual(values=colors_yr[c(2,3,4)]) + 
  scale_shape_manual(values=c(16,17,15),breaks=c("CAC_stream1","CAC_stream2","LM"),labels=c("CAC_S1","CAC_S2","LM"))
PC5  
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/SEP_versions/PC5_sookID.pdf",PC5,device = 'pdf',height = 3,width=4,dpi=600,units='in')
