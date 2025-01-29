
####
#### Final Figure 4 plots
####

### Also includes code for Supplemental figures S7 and S8

library(gam)
library(mgcv)

### Panels 4A, 4B, and 4C: plot histograms for three streams
plot_hists_stream1only<-allsamples_PCAinfo %>% filter(isMom,total>=50000,!is.na(plot),isSOOK %in% c("no","F1"),stream %in% c("CAC_stream1","CAC_stream2","LM")) %>%
  mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2","LM")),
         Year=factor(Year,levels=c(2012,2019,2021,2022)),
         plot_renamed=factor(plot_renamed,levels=paste0("S1_",rev(1:8)))) %>% 
  filter(stream=="CAC_stream1") %>%
  ggplot(aes(x=nas_prop,fill=Year)) + geom_histogram(binwidth=0.05,boundary=1,position="stack") + 
  facet_grid(plot_renamed~"CAC_S1",scales="free_y") + theme_bw() + 
  scale_fill_manual(values=colors_yr,guide="none") + xlab("") + ylab("Maternal samples per plot (years are stacked)") + 
  scale_x_continuous(breaks=c(0:5)/5,name="Hybrid index") + scale_y_continuous(breaks=c(0:2)*20,limits=c(0,40))

plot_hists_stream2only<-allsamples_PCAinfo %>% filter(isMom,total>=50000,!is.na(plot),isSOOK %in% c("no","F1"),stream %in% c("CAC_stream1","CAC_stream2","LM")) %>%
  mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2","LM")),
         Year=factor(Year,levels=c(2012,2019,2021,2022)),
         plot_renamed=factor(plot_renamed,levels=paste0("S2_",rev(1:7)))) %>% 
  filter(stream=="CAC_stream2") %>%
  ggplot(aes(x=nas_prop,fill=Year)) + geom_histogram(binwidth=0.05,boundary=1,position="stack") + 
  facet_grid(plot_renamed~"CAC_S2",scales="free_y") + theme_bw() + 
  scale_fill_manual(values=colors_yr[c(1,3,4)],guide="none") + 
  scale_x_continuous(breaks=c(0:5)/5,name="Hybrid index") + scale_y_continuous(breaks=c(0:2)*8,limits=c(0,16),name="")

plot_hists_LMonly<-allsamples_PCAinfo %>% filter(isMom,total>=50000,!is.na(plot),isSOOK %in% c("no","F1"),stream %in% c("CAC_stream1","CAC_stream2","LM")) %>%
  mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2","LM")),
         Year=factor(Year,levels=c(2012,2019,2021,2022)),
         plot_renamed=factor(str_replace(plot_renamed,"LM","LM_"),levels=paste0("LM_",rev(c(1:8))))) %>% 
  filter(stream=="LM") %>%
  ggplot(aes(x=nas_prop,fill=Year)) + geom_histogram(binwidth=0.05,boundary=1,position="stack") + 
  facet_grid(plot_renamed~"LM",scales="free_y") + theme_bw() + 
  scale_fill_manual(values=colors_yr[3],guide="none") + 
  scale_x_continuous(breaks=c(0:5)/5,limits=c(0,1),name="Hybrid index") + scale_y_continuous(breaks=c(0:2)*2,limits=c(0,4),name="")

plot_hists_threestreams<-plot_grid(plot_hists_stream1only,plot_hists_stream2only,plot_hists_LMonly,ncol = 3)


### Figure 4D: plotwise phenology: prepare phenology data
phenology_ready<-phenology %>% mutate(open_flowers=as.integer(open_flowers),plot_renamed=sapply(plot,get_univ_plot,key=plot_key,zone="renamed")) %>% 
  filter(open_flowers>=0,!is.na(open_flowers)) %>%
  group_by(Year,plot_renamed,Days_start_Apr1) %>% summarize(open_flowers=sum(open_flowers,na.rm=T)) %>% 
  mutate(flowers_proportion_of_max=open_flowers/max(open_flowers,na.rm=T),Year=factor(Year)) %>% ungroup() %>% 
  mutate(plot_renamed=factor(plot_renamed),
         stream=ifelse(str_starts(plot_renamed,"S1"),"CAC_stream1",
                       ifelse(str_starts(plot_renamed,"S2"),"CAC_stream2","LM")))

phenology_quartiles<-phenology_ready %>% group_by(Year,plot_renamed) %>% filter(open_flowers>0) %>%
  summarize(min_flowerdate=min(Days_start_Apr1),
            Q1_flowerdate=wtd.quantile(Days_start_Apr1,weights=open_flowers,probs=0.25),
            median_flowerdate=wtd.quantile(Days_start_Apr1,weights=open_flowers,probs=0.5),
            Q3_flowerdate=wtd.quantile(Days_start_Apr1,weights=open_flowers,probs=0.75),
            max_flowerdate=max(Days_start_Apr1))

HI_plot_quartiles<-allsamples_PCAinfo %>% filter(isMom,total>=50000,!is.na(plot),isSOOK %in% c("no","F1"),stream %in% c("CAC_stream1","CAC_stream2")) %>%
  mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2")),
         #Year=factor(Year,levels=c(2022,2021,2019,2012))) %>%
         Year=factor(Year,levels=c(2012,2019,2021,2022))) %>%
  group_by(Year,stream,plot_renamed) %>% summarize(min_HI=min(nas_prop),
                                                     Q1_HI=quantile(nas_prop,0.25),
                                                     median_HI=median(nas_prop),
                                                     Q3_HI=quantile(nas_prop,0.75),
                                                     max_HI=max(nas_prop)) %>%
  inner_join(phenology_quartiles,by=c("Year","plot_renamed"))

###Figure 4D plot
plot_quartiles_minmax<-HI_plot_quartiles %>% ungroup() %>%
   ggplot(aes(x=median_flowerdate,y=median_HI,color=Year,shape=stream)) +
  geom_point() + theme_bw() + 
  facet_grid(Year~.) + 
  geom_segment(aes(x=min_flowerdate,xend=max_flowerdate,y=median_HI,yend=median_HI,color=Year),linewidth=0.2) + 
  scale_color_manual(values=colors_yr[c(1,2,4)],guide="none") + scale_shape_discrete(guide="none") + 
  scale_x_continuous(breaks=c(1,31,62,92),limits=c(1,92),
                     minor_breaks=NULL,name="Flowering date range",
                     labels=c("April 1","May 1","June 1","July 1")) + 
  ylab("Median hybrid index") + theme(axis.text.x=element_text(angle=45,hjust=1))

plot_quartiles_minmax


### Figure 4E: individual fruit timings
# Plot betareg confidence intervals using gam
markedfruits_strict_modelled<-markedfruits_strict %>% mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2")),
                                                             Year=factor(Year,levels=c(2019,2022))) %>%
  filter(mom_isSOOK == "no") %>%
  mutate(mom_nas_prop=ifelse(mom_nas_prop==1,0.999,mom_nas_prop))

mygam = gam(mom_nas_prop~Days_start_Apr1*Year*stream,data=markedfruits_strict_modelled,family=betar(link="logit"))
min <- min(markedfruits_strict_modelled$Days_start_Apr1)
max <- max(markedfruits_strict_modelled$Days_start_Apr1)
new.x.2019 <- data.frame(Days_start_Apr1 = seq(min, max, length.out = 1000),Year="2019",stream="CAC_stream1")
new.x.2022.stream1 <- data.frame(Days_start_Apr1 = seq(min, max, length.out = 1000),Year="2022",stream="CAC_stream1")
new.x.2022.stream2 <- data.frame(Days_start_Apr1 = seq(min, max, length.out = 1000),Year="2022",stream="CAC_stream2")
new.x<-rbind(new.x.2019,new.x.2022.stream1,new.x.2022.stream2)
new.y <- predict(mygam, newdata = new.x, se.fit = TRUE, type="response")
new.y <- data.frame(new.y)
addThese <- data.frame(new.x, new.y)
addThese <- rename(addThese, mom_nas_prop = fit, SE = se.fit)
addThese <- mutate(addThese, lwr = mom_nas_prop - 1.96 * SE, upr = mom_nas_prop + 1.96 * SE) # calculating the 95% confidence interval

fruit_timing_modelled_se<-markedfruits_strict_modelled %>% 
  ggplot(aes(x=Days_start_Apr1,y=mom_nas_prop,color=Year,shape=stream,fill=Year,linetype=stream)) + 
  geom_smooth(data = addThese, aes(ymin = lwr, ymax = upr), stat = 'identity',alpha=0.3) + 
  geom_point() + theme_bw() + 
  scale_color_manual(values=colors_yr[c(2,4)],guide="none",aesthetics = c("color","fill")) + 
  scale_shape_discrete(guide="none") + 
  ylab("Hybrid index") + 
  scale_x_continuous(breaks=c(1,31,62,92),limits=c(1,92),
                     minor_breaks=NULL,
                     labels=c("April 1","May 1","June 1","July 1"),
                     name="Date of marked open flower") + 
  geom_hline(yintercept=0.15,color="grey",linetype="dashed") + geom_hline(yintercept=0.95,color="grey",linetype="dashed") + 
  geom_vline(xintercept=44,color="grey",linetype="dashed") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  scale_linetype_discrete(guide="none")
fruit_timing_modelled_se

fruit_timing_modelled_se_legend<-markedfruits_strict_modelled %>% 
  ggplot(aes(x=Days_start_Apr1,y=mom_nas_prop,color=Year,shape=stream,fill=Year,linetype=stream)) + 
  geom_point() + theme_bw() + 
  scale_color_manual(values=colors_yr[c(2,4)],aesthetics=c("color","fill")) + 
  scale_shape_discrete(name="Stream",labels=c("CAC_S1","CAC_S2"),guide=guide_legend(override.aes=list("color"="black"))) + 
  geom_smooth(data = addThese, aes(ymin = lwr, ymax = upr), stat = 'identity') + 
  ylab("Maternal hybrid index") + xlab("Date of marked open flower") + 
  scale_x_continuous(breaks=c(1,31,62,92),limits=c(1,92),
                     minor_breaks=NULL,name=NULL,
                     labels=c("April 1","May 1","June 1","July 1")) + 
  geom_hline(yintercept=0.15,color="grey",linetype="dashed") + geom_hline(yintercept=0.95,color="grey",linetype="dashed") + 
  geom_vline(xintercept=44,color="grey",linetype="dashed") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  scale_linetype_discrete(name="Stream",labels=c("CAC_S1","CAC_S2"),guide=guide_legend(override.aes=list("color"="black")))
fruit_timing_modelled_se_legend
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/DEC_revision/fruit_timing_byancestry_modelled_withse.pdf",fruit_timing_modelled_se_legend,device = 'pdf',height = 3,width=5,dpi=600,units='in')

### Putting it all together
Fig4_newlayout<-plot_grid(
  plot_grid(plot_hists_stream1only,plot_hists_stream2only,plot_hists_LMonly,ncol=3),
  plot_grid(plot_quartiles_minmax,fruit_timing_modelled_se,ncol=2),
  nrow=2
)

ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/DEC_revision/Fig3_newlayout.pdf",Fig3_newlayout,device = 'pdf',height = 8,width=8,dpi=600,units='in')

### Figure S7 phenology and flowering maxes
###S7A: stream 1 plots
plot_hists_stream1only_noname<-allsamples_PCAinfo %>% filter(isMom,total>=50000,!is.na(plot),isSOOK %in% c("no","F1"),stream %in% c("CAC_stream1","CAC_stream2","LM")) %>%
  mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2","LM")),
         Year=factor(Year,levels=c(2012,2019,2021,2022)),
         plot_renamed=factor(plot_renamed,levels=paste0("S1_",1:8)),
         plot_quickname=factor(as.numeric(plot_renamed),levels=rev(c(1:8)))) %>% 
  filter(stream=="CAC_stream1") %>%
  ggplot(aes(x=nas_prop,fill=Year)) + geom_histogram(binwidth=0.05,boundary=1,position="stack") + 
  facet_grid(plot_quickname~.,scales="free_y") + theme_bw() + 
  scale_fill_manual(values=colors_yr,guide="none") + xlab("") + ylab("Maternal samples per plot (years are stacked)") + 
  scale_x_continuous(breaks=c(0:5)/5,name="Hybrid index") + scale_y_continuous(breaks=c(0:2)*20,limits=c(0,40)) + 
  theme(strip.text.y=element_text(angle=0))

###S7B: stream 1 phenology
filter_list<-c("1_10","1_05","1_01","1_03","1_04","302A","303","305","312")
plotwise_flowerdata<-phenology %>% mutate(open_flowers=as.integer(open_flowers),plot_renamed=sapply(plot,get_univ_plot,key=plot_key,zone="renamed")) %>% 
  filter(open_flowers>=0,!is.na(open_flowers)) %>%
  group_by(Year,plot_renamed,plot,Days_start_Apr1) %>% summarize(open_flowers=sum(open_flowers,na.rm=T)) %>% 
  mutate(flowers_proportion_of_max=open_flowers/max(open_flowers,na.rm=T),Year=factor(Year)) %>% ungroup() %>% 
  mutate(stream=ifelse(str_starts(plot_renamed,"S1"),"CAC_stream1",ifelse(str_starts(plot_renamed,"S2"),"CAC_stream2","LM"))) %>%
  filter(stream=="CAC_stream1") %>% 
  mutate(plot_quickname=factor(as.numeric(factor(plot_renamed)),levels=rev(c(1:8))))
plotwise_maxflowers<-plotwise_flowerdata %>% 
  group_by(Year,plot_quickname,plot) %>% summarize(max_open=max(open_flowers,na.rm=T)) %>%
  summarize(median_max_open=median(max_open))
plotwise_phenology<-plotwise_flowerdata %>%
  filter(!(plot %in% filter_list)) %>%
  ggplot(aes(x=Days_start_Apr1,y=flowers_proportion_of_max,color=Year,group=plot)) + facet_grid(plot_quickname~.) + 
  geom_line() + 
  scale_y_continuous(breaks=c(0,0.5,1)) + scale_x_continuous(breaks=c(30,60,90)) + 
  scale_color_manual(values=colors_yr[c(1,2,4)],guide="none") + 
  theme_bw() + ylab("Open flowers (Proportion of max)") + 
  scale_x_continuous(breaks=c(1,31,62,92),limits=c(1,92),
                     minor_breaks=NULL,name="Flower date",
                     labels=c("April 1","May 1","June 1","July 1")) + 
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.y=element_text(angle=0))

###S7C: stream 1 flower maxes

plotwise_maxflowers_bars<-plotwise_maxflowers %>%
  ggplot(aes(y=Year,x=median_max_open,fill=Year)) + 
  facet_grid(plot_quickname~.) + 
  geom_col(position=position_dodge(preserve = "single"),width=1,just = 0.5) + 
  scale_fill_manual(values=colors_yr[c(1,2,4)],guide="none") + 
  theme_bw() + xlab("Max open flowers") + 
  theme(strip.text.y=element_text(angle=0))

### Figure S8 phenology and flowering maxes for stream 2
###S8A: stream 2 plots
plot_hists_stream2only_noname<-allsamples_PCAinfo %>% filter(isMom,total>=50000,!is.na(plot),isSOOK %in% c("no","F1"),stream %in% c("CAC_stream1","CAC_stream2","LM")) %>%
  mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2","LM")),
         Year=factor(Year,levels=c(2012,2019,2021,2022)),
         plot_renamed=factor(plot_renamed,levels=paste0("S2_",1:7)),
         plot_quickname=factor(as.numeric(plot_renamed),levels=rev(c(1:7)))) %>% 
  filter(stream=="CAC_stream2") %>%
  ggplot(aes(x=nas_prop,fill=Year)) + geom_histogram(binwidth=0.05,boundary=1,position="stack") + 
  facet_grid(plot_quickname~.,scales="free_y") + theme_bw() + 
  scale_fill_manual(values=colors_yr[c(1,3,4)],guide="none") + xlab("") + ylab("Maternal samples per plot (years are stacked)") + 
  scale_x_continuous(breaks=c(0:5)/5,name="Hybrid index") + scale_y_continuous(breaks=c(0:2)*20,limits=c(0,40)) + 
  theme(strip.text.y=element_text(angle=0))

###S8B: stream 2 phenology
filter_list_S2<-c("2_09","2_05","2_06","2_07","2_01","2_02")
plotwise_flowerdata_stream2<-phenology %>% mutate(open_flowers=as.integer(open_flowers),plot_renamed=sapply(plot,get_univ_plot,key=plot_key,zone="renamed")) %>% 
  filter(open_flowers>=0,!is.na(open_flowers)) %>%
  group_by(Year,plot_renamed,plot,Days_start_Apr1) %>% summarize(open_flowers=sum(open_flowers,na.rm=T)) %>% 
  mutate(flowers_proportion_of_max=open_flowers/max(open_flowers,na.rm=T),Year=factor(Year)) %>% ungroup() %>% 
  mutate(stream=ifelse(str_starts(plot_renamed,"S1"),"CAC_stream1",ifelse(str_starts(plot_renamed,"S2"),"CAC_stream2","LM"))) %>%
  filter(stream=="CAC_stream2") %>% 
  mutate(plot_quickname=factor(as.numeric(factor(plot_renamed)),levels=rev(c(1:7))))
plotwise_maxflowers_stream2<-plotwise_flowerdata_stream2 %>% 
  group_by(Year,plot_quickname,plot) %>% summarize(max_open=max(open_flowers,na.rm=T)) %>%
  summarize(median_max_open=median(max_open))
plotwise_phenology_stream2<-plotwise_flowerdata_stream2 %>%
  filter(!(plot %in% filter_list_S2)) %>%
  ggplot(aes(x=Days_start_Apr1,y=flowers_proportion_of_max,color=Year,group=plot)) + facet_grid(plot_quickname~.) + 
  geom_line() + 
  scale_y_continuous(breaks=c(0,0.5,1)) + scale_x_continuous(breaks=c(30,60,90)) + 
  scale_color_manual(values=colors_yr[c(1,4)],guide="none") + 
  theme_bw() + ylab("Open flowers (Proportion of max)") + 
  scale_x_continuous(breaks=c(1,31,62,92),limits=c(1,92),
                     minor_breaks=NULL,name="Flower date",
                     labels=c("April 1","May 1","June 1","July 1")) + 
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.y=element_text(angle=0))

###S8C: stream 2 flower maxes

plotwise_maxflowers_bars_stream2<-plotwise_maxflowers_stream2 %>%
  ggplot(aes(y=Year,x=median_max_open,fill=Year)) + 
  facet_grid(plot_quickname~.) + 
  geom_col(position=position_dodge(preserve = "single"),width=1,just = 0.5) + 
  scale_fill_manual(values=colors_yr[c(1,4)],guide="none") + 
  theme_bw() + xlab("Max open flowers") + 
  theme(strip.text.y=element_text(angle=0))


### Putting it together: Figures S7 and S8
Rcomplete_FigS7<-plot_hists_stream1only_noname+plotwise_phenology+plotwise_maxflowers_bars+plot_layout(design="ABC")
Rcomplete_FigS8_stream2<-plot_hists_stream2only_noname+plotwise_phenology_stream2+plotwise_maxflowers_bars_stream2+plot_layout(design="ABC")

ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/DEC_revision/Stream1_all_FigS7_top.pdf",Rcomplete_FigS7,device = 'pdf',height = 6,width=8,dpi=600,units='in')
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/DEC_revision/Stream2_all_FigS8_top.pdf",Rcomplete_FigS8_stream2,device = 'pdf',height = 6,width=8,dpi=600,units='in')

###NOT USED: old version of Figure 4D phenology vs. HI
plot_quartiles<-HI_plot_quartiles %>% ungroup() %>%
  #add_row(Year="2021","Q1_flowerdate"=0,"Q3_flowerdate"=0,"median_HI"=0,"median_flowerdate"=0,"Q1_HI"=0,"Q3_HI"=0,"min_HI"=0,"max_HI"=0,"min_flowerdate"=0,"max_flowerdate"=0,"stream"="CAC_stream1") %>%
  #filter(stream=="CAC_stream1") %>% 
  ggplot(aes(x=median_flowerdate,y=median_HI,color=Year,shape=stream)) +
  geom_point() + theme_bw() + 
  geom_segment(aes(x=Q1_flowerdate,xend=Q3_flowerdate,y=median_HI,yend=median_HI,color=Year),linewidth=0.2) + 
  geom_segment(aes(x=median_flowerdate,xend=median_flowerdate,y=Q1_HI,yend=Q3_HI,color=Year),linewidth=0.2) + 
  #geom_segment(aes(x=min_flowerdate,xend=max_flowerdate,y=median_HI,yend=median_HI,color=Year),linetype="dotted",linewidth=0.2) + 
  #geom_segment(aes(x=median_flowerdate,xend=median_flowerdate,y=min_HI,yend=max_HI,color=Year),linetype="dotted",linewidth=0.2) + 
  #geom_point(aes(x=max_flowerdate,y=median_HI,color=Year,shape=NULL),shape=4,size=0.5) + 
  #geom_point(aes(x=min_flowerdate,y=median_HI,color=Year,shape=NULL),shape=4,size=0.5) + 
  scale_color_manual(values=colors_yr[c(1,2,4)],guide="none") + scale_shape_discrete(guide="none") + 
  scale_x_continuous(breaks=c(1,31,62,92),limits=c(1,92),
                     minor_breaks=NULL,name=NULL,
                     labels=c("April 1","May 1","June 1","July 1")) + 
  ylab("Hybrid index by plot") + theme(axis.text.x=element_text(angle=45,hjust=1))
plot_quartiles_legend<-HI_plot_quartiles %>% ungroup() %>%
  #add_row(Year="2021","Q1_flowerdate"=0,"Q3_flowerdate"=0,"median_HI"=0,"median_flowerdate"=0,"Q1_HI"=0,"Q3_HI"=0,"min_HI"=0,"max_HI"=0,"min_flowerdate"=0,"max_flowerdate"=0,"stream"="CAC_stream1") %>%
  #filter(stream=="CAC_stream1") %>% 
  ggplot(aes(x=median_flowerdate,y=median_HI,color=Year,shape=stream)) +
  geom_point() + theme_bw() + 
  geom_segment(aes(x=Q1_flowerdate,xend=Q3_flowerdate,y=median_HI,yend=median_HI,color=Year),linewidth=0.2) + 
  geom_segment(aes(x=median_flowerdate,xend=median_flowerdate,y=Q1_HI,yend=Q3_HI,color=Year),linewidth=0.2) + 
  geom_segment(aes(x=min_flowerdate,xend=max_flowerdate,y=median_HI,yend=median_HI,color=Year),linetype="dotted",linewidth=0.2) + 
  geom_segment(aes(x=median_flowerdate,xend=median_flowerdate,y=min_HI,yend=max_HI,color=Year),linetype="dotted",linewidth=0.2) + 
  #geom_point(aes(x=max_flowerdate,y=median_HI,color=Year,shape=NULL),shape=4,size=0.5) + 
  #geom_point(aes(x=min_flowerdate,y=median_HI,color=Year,shape=NULL),shape=4,size=0.5) + 
  scale_color_manual(values=colors_yr[c(1,2,4)]) + scale_shape_discrete(labels=c("CAC_S1","CAC_S2")) + 
  scale_x_continuous(breaks=c(1,31,62,92),limits=c(1,92),
                     minor_breaks=NULL,name=NULL,
                     labels=c("April 1","May 1","June 1","July 1")) + 
  ylab("Hybrid index by plot") + theme(axis.text.x=element_text(angle=45,hjust=1))
plot_quartiles

