
###Calculations of RI

phenology_forRI_allplots<-phenology_ready %>% filter(stream=="CAC_stream1") %>% mutate(plot_group=ifelse(plot_renamed %in% c("S1_5","S1_8"),"guttatus",
                                                                                           ifelse(plot_renamed %in% c("S1_2","S1_4"),"hybrid",
                                                                                                  ifelse(plot_renamed %in% c("S1_3","S1_6"),"nasutus","drop")))) %>%
  filter(plot_group!="drop") %>% group_by(Year,plot_group,Days_start_Apr1) %>% summarize(total_flowers=sum(open_flowers,na.rm=T)) %>%
  pivot_wider(names_from=plot_group,values_from=total_flowers) %>% 
  mutate(guttatus=ifelse(is.na(guttatus),0,guttatus),
         hybrid=ifelse(is.na(hybrid),0,hybrid),
         nasutus=ifelse(is.na(nasutus),0,nasutus)) %>% ungroup() %>% group_by(Year) %>%
  mutate(count_GN=(guttatus*nasutus)/(guttatus+nasutus),
         count_GH=(guttatus*hybrid)/(guttatus+hybrid),
         count_NH=(nasutus*hybrid)/(nasutus+hybrid)) %>%
  summarize(GN=sum(count_GN,na.rm=T),GH=sum(count_GH,na.rm=T),NH=sum(count_NH,na.rm=T),
            Gtotal=sum(guttatus,na.rm=T),Htotal=sum(hybrid,na.rm=T),Ntotal=sum(nasutus,na.rm=T)) %>%
  mutate(prob_G_N=GN/Gtotal,prob_N_G=GN/Ntotal,
         prob_G_H=GH/Gtotal,prob_H_G=GH/Htotal,
         prob_N_H=NH/Ntotal,prob_H_N=NH/Htotal)

phenology_forRI<-phenology_ready %>% filter(stream=="CAC_stream1") %>% mutate(plot_group=ifelse(plot_universal %in% c("603"),"guttatus",
                                                                                                         ifelse(plot_universal %in% c("605"),"hybrid",
                                                                                                                ifelse(plot_universal %in% c("602"),"nasutus","drop")))) %>%
  filter(plot_group!="drop") %>% group_by(Year,plot_group,Days_start_Apr1) %>% summarize(total_flowers=sum(open_flowers,na.rm=T)) %>%
  pivot_wider(names_from=plot_group,values_from=total_flowers) %>% 
  mutate(guttatus=ifelse(is.na(guttatus),0,guttatus),
         hybrid=ifelse(is.na(hybrid),0,hybrid),
         nasutus=ifelse(is.na(nasutus),0,nasutus)) %>% ungroup() %>% group_by(Year) %>%
  mutate(count_GN=(guttatus*nasutus)/(guttatus+nasutus),
         count_GH=(guttatus*hybrid)/(guttatus+hybrid),
         count_NH=(nasutus*hybrid)/(nasutus+hybrid)) %>%
  summarize(GN=sum(count_GN,na.rm=T),GH=sum(count_GH,na.rm=T),NH=sum(count_NH,na.rm=T),
            Gtotal=sum(guttatus,na.rm=T),Htotal=sum(hybrid,na.rm=T),Ntotal=sum(nasutus,na.rm=T)) %>%
  mutate(prob_G_N=GN/Gtotal,prob_N_G=GN/Ntotal,
         prob_G_H=GH/Gtotal,prob_H_G=GH/Htotal,
         prob_N_H=NH/Ntotal,prob_H_N=NH/Htotal)

phenology_RI_summarytable<-phenology_forRI %>% select(Year,starts_with("prob")) %>% 
  pivot_longer(cols = starts_with("prob"),names_to="Comparison",values_to = "Probability") %>%
  mutate(RI=1-(2*Probability))

phenology_RI_summarytable_allplots<-phenology_forRI_allplots %>% select(Year,starts_with("prob")) %>% 
  pivot_longer(cols = starts_with("prob"),names_to="Comparison",values_to = "Probability") %>%
  mutate(RI=1-(2*Probability))

#〖RI〗_Combined=〖RI〗_M+((1-〖RI〗_M )*〖RI〗_P )
#RIcombined=RI_M+((1-RI_M)*RI_P)

###Phenological RI for stream 2
phenology_forRI_stream2<-phenology_ready %>% filter(stream=="CAC_stream2",Year=="2022") %>% mutate(plot_group=ifelse(plot_renamed %in% c("S2_1","S2_2"),"nasutus", "hybrid")) %>%
  filter(plot_group!="drop") %>% group_by(Year,plot_group,Days_start_Apr1) %>% summarize(total_flowers=sum(open_flowers,na.rm=T)) %>%
  pivot_wider(names_from=plot_group,values_from=total_flowers) %>% 
  mutate(hybrid=ifelse(is.na(hybrid),0,hybrid),
         nasutus=ifelse(is.na(nasutus),0,nasutus)) %>% ungroup() %>% group_by(Year) %>%
  mutate(count_NH=(nasutus*hybrid)/(nasutus+hybrid)) %>%
  summarize(NH=sum(count_NH,na.rm=T),
            Htotal=sum(hybrid,na.rm=T),Ntotal=sum(nasutus,na.rm=T)) %>%
  mutate(prob_N_H=NH/Ntotal,prob_H_N=NH/Htotal)
phenology_RI_summarytable_stream2<-phenology_forRI_stream2 %>% select(Year,starts_with("prob")) %>% 
  pivot_longer(cols = starts_with("prob"),names_to="Comparison",values_to = "Probability") %>%
  mutate(RI=1-(2*Probability))
