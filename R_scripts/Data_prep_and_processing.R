
###Required libraries
library(Hmisc)
library(readxl)
library(raster)
library(zoo)
library(cowplot)
library(ggbeeswarm)
library(ggmap)
library(ggrepel)
library(tmap)
library(ggpattern)
library(lme4)


library(tidyverse)

### Load and prep all data and metadata
datafolder<-"" #fill in location of data
outputfolder<-"" #fill in destination location

allsamples<-read_csv("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/Combined_allsamples_metadata_ancestry.csv")
samples2012<-read_excel("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Data/Datasheets/Official_sample_metadata/CAC.STRUCT.Results.May2014.xlsx",sheet="K2.indivs.plots") %>%
  mutate(Plot=str_remove(str_remove(str_replace(Plot,"\\.","_"),"near_"),"by_"),group="CAC2012",momID=indiv,fruitID=NA,offspringID=NA,isMom=T,total=Inf,heterozygosity=NA,Year=2012,Zone=NA) %>% 
  rename("sampleID"="indiv","plot"="Plot") %>%
  dplyr::select(-`M. guttatus`,-`M. nasutus`,-K2.order) %>%
  right_join(read_excel("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/CAC2012_hapmix.50K.ancestryprops.xlsx"),by="sampleID")
bamfiles_data<-read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/allbamfiles_50000calls_correct_withancestry.txt",col.names=c("bampath","group","isMom","plot","nas_prop")) %>%
  separate(bampath,into=c("V1","V2","V3","V4","V5","bamfile"),sep="/",remove=F) %>% 
  separate(bamfile,into=c("prefix"),sep="\\.",extra = "drop",remove=T) %>%
  mutate(sampleID=str_remove(prefix,"-combined")) %>%
  select(bampath,sampleID)
phenology<-read_excel("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/all_phenology_formatted.xlsx")
phenology_plus2324<-read_excel("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/all_phenology_formatted.xlsx",sheet="add23_24")
flowerdata2022<-read_excel("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/2022field-data.xlsx",sheet="Name-corrected")
flowerdata2019<-read_excel("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Data/Datasheets/Official_sample_metadata/CAC2019_fruit_dates.xlsx")
PRISM<-read_xlsx("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/PRISM_temp_precip_2010-2022_fmt.xlsx",sheet = "Data") %>%
  separate(Date,into=c("Year","Month"),sep="-",remove = F) %>% 
  mutate(rollsum_6mo_precip_mm=rollsum(ppt_mm,k=6,align = "right",fill=NA),
         trange_degreesC=tmax_degreesC-tmin_degreesC,
         month.index=c(1:nrow(.)))
PRISM_daily<-read_xlsx("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/PRISM_temp_precip_2010-2022_fmt.xlsx",sheet = "Daily") %>%
  #separate(Date,into=c("Month","Day","Year2"),sep="/",remove = F) %>% 
  mutate(rollsum_14day_precip_mm=rollsum(ppt_mm,k=14,align = "right",fill=NA),
         rollmean_14day_tmean_degreesC=rollmean(tmean_degreesC,k=14,align = "right",fill=NA))
panel100<-read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/allyears_CACpanel100_bamlist.txt",col.names="bampath") %>% 
  separate(bampath,into=c("V1","V2","V3","V4","V5","bamfile"),sep="/",remove=F) %>% 
  separate(bamfile,into=c("prefix"),sep="\\.",extra = "drop",remove=T) %>%
  mutate(sampleID=str_remove(prefix,"-combined")) %>% select(sampleID) %>%
  left_join(allsamples_PCAinfo)
  

plot_key<-read_csv("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/Plot_key.csv",col_types = "cc")

colors_sp<-c("#993929","#CC892F","#FFDA35") #G,A,N
#colors_yr<-c("#E228E2","red","skyblue","#005200")
colors_yr<-c("#B039E3","#FF7800","#86CEEB","#015100")

stream1_plots<-c("601","601B","602","603","604","605","606")
stream2_plots<-c("607","608","609","610","610B","611","612")
LM_plots<-paste0("LM",c(1:8))

###Processing functions
get_univ_plot<-function(x,key,zone=TRUE) {
  if (zone=="renamed") {
    if (is.na(x) || str_sub(x,1,2)=="LM") {
      return(x)
    } else if (x %in% key$Previous) {
      return(key[key$Previous==x,]$Renamed)
    } else {
      return(NA)
    }
  }
  if (zone) {
    return(key$Zone[match(x,key$Plot)])
  }
  if (is.na(x) || str_sub(x,1,1)=="6" || str_sub(x,1,2)=="LM") {
    return(x)
  } else if (x %in% key$Previous) {
    return(key[key$Previous==x,]$CAC2022)
  } else {
    return(NA)
  }
}

###Processing all samples

allsamples_with2012<-bind_rows(allsamples,samples2012)
allsamples_info<-allsamples_with2012 %>% 
  mutate(momfullID=paste0(group,"_",momID),
         plot_universal=sapply(plot,get_univ_plot,key=plot_key,zone=F),
         plot_renamed=sapply(plot,get_univ_plot,key=plot_key,zone="renamed"),
         stream=ifelse(plot_universal %in% stream1_plots,"CAC_stream1",
                       ifelse(plot_universal %in% stream2_plots,"CAC_stream2",
                              ifelse(group=="LM2021","LM",
                                     ifelse(group=="CAC-E2021","CAC-E",NA)))),
         streamyear=paste0(stream,"_",Year),
         nasprop_cohort=ifelse(nas_prop>0.8,"nas",
                               ifelse(nas_prop>0.44 & nas_prop<0.55 & !is.na(heterozygosity) & heterozygosity>0.80,"sook_or_F1",
                                      ifelse(nas_prop>0.5,"admix_high",
                                             ifelse(nas_prop>0.15,"admix","gut"))))) %>%
  filter(sampleID!="KK047") %>% #exclude KK047 mom because two sequencing batches disagreed on ancestry
  left_join(bamfiles_data,by="sampleID")

allsamples_info_50000<-allsamples_info %>% filter(total>=50000)

###identifying and flagging sookensis samples
moms_sook_potential<-allsamples_info %>% filter(isMom,nasprop_cohort=="sook_or_F1")
offspring_sook_potential<-allsamples_info %>% filter(!isMom,nasprop_cohort=="sook_or_F1")
moms_sook_ids<-moms_sook_potential %>% mutate(isSOOK=ifelse(momfullID %in% offspring_sook_potential$momfullID,"confirmed",
                                                            ifelse(plot %in% c("611","612"),"likely","F1"))) %>%
  select(sampleID,isSOOK)
offspring_sook_ids<-offspring_sook_potential %>% mutate(isSOOK=ifelse(momfullID %in% moms_sook_potential$momfullID,"confirmed",
                                                                      ifelse(plot %in% c("611","612"),"likely","F1"))) %>% 
  select(sampleID,isSOOK)
all_sook_ids<-rbind(moms_sook_ids,offspring_sook_ids)

allsamples_sookinfo<-allsamples_info %>% left_join(all_sook_ids,by="sampleID") %>%
  mutate(isSOOK=ifelse(!is.na(isSOOK),isSOOK,"no"))

### PCA and Structure data
allnorthmoms_structure<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/angsd_Apr2024/allnorthmoms_nofilter.IM62_v3.CACpanel100.maf20.nonmiss60.combined.qopt.samples",
                                            header=F,col.names=c("bampath","K1_prop","K2_prop"))

allnorthmoms_PCAraw<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/angsd_Apr2024/allnorthmoms_nofilter.IM62_v3.CACpanel100.maf20.nonmiss60.combined.cov")
allnorthmoms_eigen<-eigen(allnorthmoms_PCAraw)
allnorthmoms_PCs<-cbind(allnorthmoms_structure,as.data.frame(allnorthmoms_eigen$vectors[,1:6]))

#define clusters based on PCA
allsamples_PCAinfo<-allsamples_sookinfo %>%
  left_join(allnorthmoms_PCs,by="bampath") %>%
  rename("PC1"="V1","PC2"="V2","PC3"="V3","PC4"="V4","PC5"="V5","PC6"="V6") %>%
  mutate(PCAcluster=ifelse(PC1>0.075,"PC1_nasutus",
                           ifelse(PC2>0.04,"PC2_LM",
                                  ifelse(PC3>0.025,"PC3_lower_east",
                                         ifelse(PC3<(-0.025) & PC4>0.02,"PC4_mid_east",
                                                ifelse(PC5>0.2,"PC5_sookensis","PC3_PC4_westeast"))))))

###Filtering based on unlikely mom-offspring combos
###See "~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Code/borice_inputprep.R"
strict_filter_list<-rbind(
  read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2019_strict_borice_families.txt",col.names=c("bampath","BORICE_FAM","paroff")),
  read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2021_strict_borice_families.txt",col.names=c("bampath","BORICE_FAM","paroff")),
  read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2022_strict_borice_families.txt",col.names=c("bampath","BORICE_FAM","paroff"))
)
allsamples_strictfilter<-allsamples_PCAinfo %>% filter(bampath %in% strict_filter_list$bampath)

### Flower date info for moms and offspring
momdates2019<-flowerdata2019 %>% dplyr::select(momID,Days_start_Apr1) %>% 
  group_by(momID) %>% summarize(first_marked=min(Days_start_Apr1),last_marked=max(Days_start_Apr1)) %>% mutate(group="CAC2019")

offspringdates2019<-allsamples_PCAinfo %>% filter(group=="CAC2019",!isMom) %>% select(momID,fruitID) %>% unique() %>%
  mutate(fruitIDtemp=str_remove(fruitID,"fr-"),
         Fruit_color=str_to_lower(str_remove(fruitIDtemp,"[AB]$"))) %>% 
  left_join(flowerdata2019,by=c("momID","Fruit_color")) %>% select(momID,fruitID,Days_start_Apr1) %>% mutate(group="CAC2019")

momdates2022<-flowerdata2022 %>% select(momID,Days_start_Apr1) %>% 
  group_by(momID) %>% summarize(first_marked=min(Days_start_Apr1),last_marked=max(Days_start_Apr1)) %>% mutate(group="CAC2022")

offspringdates2022<-flowerdata2022 %>% filter(!is.na(GH_progeny)) %>% 
  dplyr::select(momID,Color_marked,Days_start_Apr1) %>% 
  dplyr::rename("fruitID"="Color_marked") %>% mutate(group="CAC2022")


momdates_tocombine<-allsamples_PCAinfo %>% filter(total>=50000,group %in% c("CAC2019","CAC2022"),!(str_detect(sampleID,"Census")),isMom) %>%
  left_join(rbind(momdates2019,momdates2022),by=c("momID","group")) %>% select(-fruitID,-offspringID,-isMom) %>%
  rename("mom_nas_prop"="nas_prop","mom_het"="heterozygosity","mom_sampleID"="sampleID",
         "mom_total"="total","mom_isSOOK"="isSOOK","mom_nasprop_cohort"="nasprop_cohort","mom_bampath"="bampath")
markedoffspring<-allsamples_PCAinfo %>% filter(total>=50000,group %in% c("CAC2019","CAC2022"),!(str_detect(sampleID,"Census")),!isMom) %>%
  left_join(rbind(offspringdates2019,offspringdates2022),by=c("group","momID","fruitID"))
markedfruits<-markedoffspring %>% group_by(group,momID,fruitID) %>%
  summarize(fruitmean_nas_prop=mean(nas_prop),Days_start_Apr1=mean(Days_start_Apr1),fruitmean_het=mean(heterozygosity,na.rm=T)) %>% 
  left_join(momdates_tocombine,by=c("group","momID")) %>% 
  mutate(momfruit_sway=fruitmean_nas_prop-mom_nas_prop,momfruit_hetsway=fruitmean_het-mom_het) %>%
  filter(!is.na(mom_nas_prop),!is.na(Days_start_Apr1))

#155 marked fruits without filtering (except total>=50000)

### Flower date info for moms and offspring
momdates_tocombine_strict<-allsamples_strictfilter %>% filter(group %in% c("CAC2019","CAC2022"),!(str_detect(sampleID,"Census")),isMom) %>%
  left_join(rbind(momdates2019,momdates2022),by=c("momID","group")) %>% select(-fruitID,-offspringID,-isMom) %>%
  rename("mom_nas_prop"="nas_prop","mom_het"="heterozygosity","mom_sampleID"="sampleID",
         "mom_total"="total","mom_isSOOK"="isSOOK","mom_nasprop_cohort"="nasprop_cohort","mom_bampath"="bampath")
markedoffspring_strict<-allsamples_strictfilter %>% filter(group %in% c("CAC2019","CAC2022"),!(str_detect(sampleID,"Census")),!isMom) %>%
  left_join(rbind(offspringdates2019,offspringdates2022),by=c("group","momID","fruitID"))
markedfruits_strict<-markedoffspring_strict %>% group_by(group,momID,fruitID) %>%
  summarize(fruitmean_nas_prop=mean(nas_prop),Days_start_Apr1=mean(Days_start_Apr1),fruitmean_het=mean(heterozygosity,na.rm=T)) %>% 
  left_join(momdates_tocombine_strict,by=c("group","momID")) %>% 
  mutate(momfruit_sway=fruitmean_nas_prop-mom_nas_prop,momfruit_hetsway=fruitmean_het-mom_het) %>%
  filter(!is.na(mom_nas_prop),!is.na(Days_start_Apr1))

#143 marked fruits with strict filtering

###Fruit means before strict filter
allfruits_moms_lax<-allsamples_PCAinfo %>% filter(!(str_detect(sampleID,"Census")),isMom,total>=50000,isSOOK %in% c("no","F1"),Year!="2012") %>%
  rename("mom_nas_prop"="nas_prop","mom_het"="heterozygosity","mom_sampleID"="sampleID",
         "mom_total"="total","mom_isSOOK"="isSOOK","mom_nasprop_cohort"="nasprop_cohort","mom_bampath"="bampath") %>%
  select(-fruitID)
allfruits_offspring_lax<-allsamples_PCAinfo %>% filter(!isMom,total>=50000,isSOOK %in% c("no","F1"),Year!="2012") %>% 
  select(group,momID,fruitID,sampleID,bampath,total,nas_prop,heterozygosity,isSOOK) %>%
  left_join(allfruits_moms_lax,by=c("group","momID")) %>% mutate(offspring_sway=nas_prop-mom_nas_prop,offspring_hetsway=heterozygosity-mom_het) %>%
  filter(!is.na(mom_total))
allfruits_fruits_lax<-allsamples_PCAinfo %>% filter(!isMom,total>=50000,isSOOK %in% c("no","F1"),Year!="2012") %>% 
  select(group,momID,fruitID,sampleID,bampath,total,nas_prop,heterozygosity,isSOOK) %>%
  group_by(group,momID,fruitID) %>%
  summarize(fruitmean_nas_prop=mean(nas_prop,na.rm=T),fruitmean_het=mean(heterozygosity,na.rm=T)) %>%
  left_join(allfruits_moms_lax,by=c("group","momID")) %>% mutate(momfruit_sway=fruitmean_nas_prop-mom_nas_prop,momfruit_hetsway=fruitmean_het-mom_het) %>%
  filter(!is.na(mom_total))

###Fruit means with looser filter
loose_filter_list<-rbind(
  read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2019_loose_borice_families.txt",col.names=c("bampath","BORICE_FAM","paroff")),
  read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2021_loose_borice_families.txt",col.names=c("bampath","BORICE_FAM","paroff")),
  read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2022_loose_borice_families.txt",col.names=c("bampath","BORICE_FAM","paroff"))
)
allsamples_loosefilter<-allsamples_PCAinfo %>% filter(bampath %in% loose_filter_list$bampath)

allfruits_moms_loose<-allsamples_loosefilter %>% filter(!(str_detect(sampleID,"Census")),isMom,total>=50000,isSOOK %in% c("no","F1"),Year!="2012") %>%
  rename("mom_nas_prop"="nas_prop","mom_het"="heterozygosity","mom_sampleID"="sampleID",
         "mom_total"="total","mom_isSOOK"="isSOOK","mom_nasprop_cohort"="nasprop_cohort","mom_bampath"="bampath") %>%
  select(-fruitID)
allfruits_offspring_loose<-allsamples_loosefilter %>% filter(!isMom,total>=50000,isSOOK %in% c("no","F1"),Year!="2012") %>% 
  select(group,momID,fruitID,sampleID,bampath,total,nas_prop,heterozygosity,isSOOK) %>%
  left_join(allfruits_moms_lax,by=c("group","momID")) %>% mutate(offspring_sway=nas_prop-mom_nas_prop,offspring_hetsway=heterozygosity-mom_het) %>%
  filter(!is.na(mom_total))
allfruits_fruits_loose<-allsamples_loosefilter %>% filter(!isMom,total>=50000,isSOOK %in% c("no","F1"),Year!="2012") %>% 
  select(group,momID,fruitID,sampleID,bampath,total,nas_prop,heterozygosity,isSOOK) %>%
  group_by(group,momID,fruitID) %>%
  summarize(fruitmean_nas_prop=mean(nas_prop,na.rm=T),fruitmean_het=mean(heterozygosity,na.rm=T)) %>%
  left_join(allfruits_moms_lax,by=c("group","momID")) %>% mutate(momfruit_sway=fruitmean_nas_prop-mom_nas_prop,momfruit_hetsway=fruitmean_het-mom_het) %>%
  filter(!is.na(mom_total))


### Fruit means including 2021
allfruits_moms<-allsamples_strictfilter %>% filter(!(str_detect(sampleID,"Census")),isMom) %>%
  rename("mom_nas_prop"="nas_prop","mom_het"="heterozygosity","mom_sampleID"="sampleID",
         "mom_total"="total","mom_isSOOK"="isSOOK","mom_nasprop_cohort"="nasprop_cohort","mom_bampath"="bampath") %>%
  select(-fruitID)
allfruits_offspring<-allsamples_strictfilter %>% filter(!isMom) %>% 
  select(group,momID,fruitID,sampleID,bampath,total,nas_prop,heterozygosity,isSOOK) %>%
  left_join(allfruits_moms,by=c("group","momID")) %>% mutate(offspring_sway=nas_prop-mom_nas_prop,offspring_hetsway=heterozygosity-mom_het)
allfruits_fruits<-allsamples_strictfilter %>% filter(!isMom) %>% 
  select(group,momID,fruitID,sampleID,bampath,total,nas_prop,heterozygosity,isSOOK) %>%
  group_by(group,momID,fruitID) %>%
  summarize(fruitmean_nas_prop=mean(nas_prop,na.rm=T),fruitmean_het=mean(heterozygosity,na.rm=T)) %>%
  left_join(allfruits_moms,by=c("group","momID")) %>% mutate(momfruit_sway=fruitmean_nas_prop-mom_nas_prop,momfruit_hetsway=fruitmean_het-mom_het)


###add selfing data
#Outself_2019<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/CAC2019_strict/strict_1.boriceready.CAC2019_strict.IM62_v3.CACpanel100.maf20.nonmiss60.cleaned.OutSelf.pp.txt",header=T) %>%
#  rowwise() %>% mutate(prob.self=self/(self+out)) %>% ungroup()
#Outself_2019<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/CAC2019_strict2/run1.boriceready.CAC2019_strict.IM62_v3.CACpanel100.maf20.nonmiss50.thinned10kb.cleaned.OutSelf.pp.txt",header=T) %>%
#  rowwise() %>% mutate(prob.self=self/(self+out)) %>% ungroup()
#Outself_2019<-read.table("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Outputs/hwe/borice_hwe/strict_hwe1.boriceready.CAC2019_strict.IM62_v3.hwe_based_SNPpanel.cleaned.OutSelf.pp.txt",header=T) %>%
#  rowwise() %>% mutate(prob.self=self/(self+out)) %>% ungroup()
Outself_2019<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/Aims_thinned10kb/AIMs_thinned10kb_1.boriceready.CAC2019_strict.IM62_v3.AIMs_thinned10kb_cleaned.OutSelf.pp.txt",header=T) %>%
   rowwise() %>% mutate(prob.self=self/(self+out)) %>% ungroup()
boricemom_info_2019<-read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2019_strict_borice_families_key.txt",header=T)
Outself_probs_mominfo_2019<-Outself_2019 %>% left_join(boricemom_info_2019,by=c("family"="BORICE_FAM")) %>% 
  mutate(selfcall=ifelse(prob.self<=0.1,"out",ifelse(prob.self>=0.9,"self",NA)),momfullID=paste0(group,"_",momID)) %>%
  group_by(momfullID,family) %>% summarize(nas_prop=mean(nas_prop),count_self=sum(selfcall=="self",na.rm=T),count_out=sum(selfcall=="out",na.rm=T),total_called=count_self+count_out,fraction_selfed=count_self/(total_called))
Outself_probs_binned_2019<-Outself_probs_mominfo_2019 %>% mutate(nas_prop_bin=factor(floor((nas_prop-0.0001)*10)/10,levels=c(0:9)/10)) %>%
  group_by(nas_prop_bin) %>% summarize(count_self=sum(count_self,na.rm=T),count_out=sum(count_out,na.rm=T)) %>%
  pivot_longer(cols = c(count_self,count_out),names_to="category",values_to="count")

#Outself_2021CAC<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/CAC2021_strictCAC/strictCAC1.boriceready.CAC2021_strictCAC.IM62_v3.CACpanel100.maf20.nonmiss60.cleaned.OutSelf.pp.txt",header=T) %>%
#  rowwise() %>% mutate(prob.self=self/(self+out)) %>% ungroup()
Outself_2021CAC<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/Aims_thinned10kb/strictCAC_AIMs_thinned10kb_1.boriceready.CAC2021_strictCAC.IM62_v3.AIMs_thinned10kb_cleaned.OutSelf.pp.txt",header=T) %>%
  rowwise() %>% mutate(prob.self=self/(self+out)) %>% ungroup()
  
boricemom_info_2021CAC<-read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2021_strictCAC_borice_families_key.txt",header=T)
Outself_probs_mominfo_2021CAC<-Outself_2021CAC %>% left_join(boricemom_info_2021CAC,by=c("family"="BORICE_FAM")) %>% 
  mutate(selfcall=ifelse(prob.self<=0.1,"out",ifelse(prob.self>=0.9,"self",NA)),momfullID=paste0(group,"_",momID)) %>%
  group_by(momfullID,family) %>% summarize(nas_prop=mean(nas_prop),count_self=sum(selfcall=="self",na.rm=T),count_out=sum(selfcall=="out",na.rm=T),total_called=count_self+count_out,fraction_selfed=count_self/(total_called))
Outself_probs_binned_2021CAC<-Outself_probs_mominfo_2021CAC %>% mutate(nas_prop_bin=factor(floor((nas_prop-0.0001)*10)/10,levels=c(0:9)/10)) %>%
  group_by(nas_prop_bin) %>% summarize(count_self=sum(count_self,na.rm=T),count_out=sum(count_out,na.rm=T)) %>%
  pivot_longer(cols = c(count_self,count_out),names_to="category",values_to="count")

#Outself_2022<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/CAC2022_strict/strict_1.boriceready.CAC2022_strict.IM62_v3.CACpanel100.maf20.nonmiss60.cleaned.OutSelf.pp.txt",header=T) %>%
#  rowwise() %>% mutate(prob.self=self/(self+out)) %>% ungroup()
#Outself_2022<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/CAC2022_strict2/run1.boriceready.CAC2022_strict.IM62_v3.CACpanel100.maf20.nonmiss50.thinned10kb.cleaned.OutSelf.pp.txt",header=T) %>%
#  rowwise() %>% mutate(prob.self=self/(self+out)) %>% ungroup()
#Outself_2022<-read.table("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Outputs/hwe/borice_hwe/strict_hwe1.boriceready.CAC2022_strict.IM62_v3.hwe_based_SNPpanel.cleaned.OutSelf.pp.txt",header=T) %>%
#  rowwise() %>% mutate(prob.self=self/(self+out)) %>% ungroup()
Outself_2022<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/Aims_thinned10kb/AIMs_thinned10kb_1.boriceready.CAC2022_strict.IM62_v3.AIMs_thinned10kb_cleaned.OutSelf.pp.txt",header=T) %>%
  rowwise() %>% mutate(prob.self=self/(self+out)) %>% ungroup()
boricemom_info_2022<-read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2022_strict_borice_families_key.txt",header=T)
Outself_probs_mominfo_2022<-Outself_2022 %>% left_join(boricemom_info_2022,by=c("family"="BORICE_FAM")) %>% 
  mutate(selfcall=ifelse(prob.self<=0.1,"out",ifelse(prob.self>=0.9,"self",NA)),momfullID=paste0(group,"_",momID)) %>%
  group_by(momfullID,family) %>% summarize(nas_prop=mean(nas_prop),count_self=sum(selfcall=="self",na.rm=T),count_out=sum(selfcall=="out",na.rm=T),total_called=count_self+count_out,fraction_selfed=count_self/(total_called))
Outself_probs_binned_2022<-Outself_probs_mominfo_2022 %>% mutate(nas_prop_bin=factor(floor((nas_prop-0.0001)*10)/10,levels=c(0:9)/10)) %>%
  group_by(nas_prop_bin) %>% summarize(count_self=sum(count_self,na.rm=T),count_out=sum(count_out,na.rm=T)) %>%
  pivot_longer(cols = c(count_self,count_out),names_to="category",values_to="count")

Outself_moms_combine<-rbind(Outself_probs_mominfo_2019,Outself_probs_mominfo_2021CAC,Outself_probs_mominfo_2022) %>%
  left_join(allfruits_moms,by="momfullID")
Outself_moms_combine_binned<-Outself_moms_combine %>% mutate(nas_prop_bin=factor(floor((nas_prop-0.00001)*10)/10,levels=c(0:9)/10)) %>%
  group_by(Year,nas_prop_bin) %>% summarize(count_self=sum(count_self,na.rm=T),count_out=sum(count_out,na.rm=T)) %>%
  pivot_longer(cols = c(count_self,count_out),names_to="category",values_to="count")

Outself_offspring_2019<-read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2019_strict_borice_families.txt",header=F) %>%
  filter(V3=="off") %>% select(V1) %>% separate(V1,into=c(NA,NA,NA,NA,NA,"filename"),sep="/") %>% 
  mutate(sampleID=str_remove(filename,".IM62_v3.fds.bam")) %>% select(sampleID) %>% bind_cols(Outself_2019) %>%
  mutate(selfcall=ifelse(prob.self<=0.1,"out",ifelse(prob.self>=0.9,"self",NA))) %>% 
  left_join(allsamples_PCAinfo,by="sampleID") %>% 
  select(sampleID,out,self,prob.self,selfcall,fruitID,total,nas_prop,heterozygosity,momfullID,nasprop_cohort,bampath,isSOOK) %>%
  left_join(allfruits_moms,by=c("momfullID")) %>%
  mutate(offspring_sway=nas_prop-mom_nas_prop,offspring_hetsway=heterozygosity-mom_het,offspring_hetprop=ifelse(mom_het==0,-1,heterozygosity/mom_het))

Outself_offspring_2021CAC<-read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2021_strictCAC_borice_families.txt",header=F) %>%
  filter(V3=="off") %>% select(V1) %>% separate(V1,into=c(NA,NA,NA,NA,NA,"filename"),sep="/") %>% 
  mutate(sampleID=str_remove(filename,".IM62_v3.fds.bam")) %>% select(sampleID) %>% bind_cols(Outself_2021CAC) %>%
  mutate(selfcall=ifelse(prob.self<=0.1,"out",ifelse(prob.self>=0.9,"self",NA))) %>% 
  left_join(allsamples_PCAinfo,by="sampleID") %>% 
  select(sampleID,out,self,prob.self,selfcall,fruitID,total,nas_prop,heterozygosity,momfullID,nasprop_cohort,bampath,isSOOK) %>%
  left_join(allfruits_moms,by=c("momfullID")) %>%
  mutate(offspring_sway=nas_prop-mom_nas_prop,offspring_hetsway=heterozygosity-mom_het,offspring_hetprop=ifelse(mom_het==0,-1,heterozygosity/mom_het))

Outself_offspring_2022<-read.table("~/Documents/GitHub/Tn5_sequencing/borice_hq/newinputs/CAC2022_strict_borice_families.txt",header=F) %>%
  filter(V3=="off") %>% select(V1) %>% separate(V1,into=c(NA,NA,NA,NA,NA,"filename"),sep="/") %>% 
  mutate(sampleID=str_remove(filename,".IM62_v3.fds.bam")) %>% select(sampleID) %>% bind_cols(Outself_2022) %>%
  mutate(selfcall=ifelse(prob.self<=0.1,"out",ifelse(prob.self>=0.9,"self",NA))) %>% 
  left_join(allsamples_PCAinfo,by="sampleID") %>% 
  select(sampleID,out,self,prob.self,selfcall,fruitID,total,nas_prop,heterozygosity,momfullID,nasprop_cohort,bampath,isSOOK) %>%
  left_join(allfruits_moms,by=c("momfullID")) %>%
  mutate(offspring_sway=nas_prop-mom_nas_prop,offspring_hetsway=heterozygosity-mom_het,offspring_hetprop=ifelse(mom_het==0,-1,heterozygosity/mom_het))

Outself_offspring_combine<-rbind(Outself_offspring_2019,Outself_offspring_2021CAC,Outself_offspring_2022)

Outself_fruits_combine<-Outself_offspring_combine %>% group_by(Year,momID,fruitID) %>% 
  summarize(count_self=sum(selfcall=="self",na.rm=T),count_out=sum(selfcall=="out",na.rm=T),total_called=count_self+count_out,fraction_selfed=count_self/(total_called)) %>%
  left_join(markedfruits_strict,by=c("Year","momID","fruitID"))

###sibships
sibships_2019<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/Aims_thinned10kb/AIMs_thinned10kb_1.boriceready.CAC2019_strict.IM62_v3.AIMs_thinned10kb_cleaned.sibships.pp.txt",header=T) %>%
  filter((full+half)>=36) %>%
  mutate(sibcall=ifelse(full>half,"full",ifelse(half>full,"half",NA)))
sibships_2021<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/Aims_thinned10kb/strictCAC_AIMs_thinned10kb_1.boriceready.CAC2021_strictCAC.IM62_v3.AIMs_thinned10kb_cleaned.sibships.pp.txt",header=T) %>%
  filter((full+half)>=36) %>%
  mutate(sibcall=ifelse(full>half,"full",ifelse(half>full,"half",NA)))
sibships_2022<-read.table("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Results/Bioinformatic_outputs/borice/Aims_thinned10kb/AIMs_thinned10kb_1.boriceready.CAC2022_strict.IM62_v3.AIMs_thinned10kb_cleaned.sibships.pp.txt",header=T) %>%
  filter((full+half)>=36) %>%
  mutate(sibcall=ifelse(full>half,"full",ifelse(half>full,"half",NA)))

sibships_info_2019<-sibships_2019 %>% group_by(family) %>% summarize(fullsibs=sum(sibcall=="full",na.rm=T),halfsibs=sum(sibcall=="half",na.rm=T)) %>% 
  mutate(fraction_halfsibs=halfsibs/(halfsibs+fullsibs)) %>%
  left_join(boricemom_info_2019,by=c("family"="BORICE_FAM")) %>%
  mutate(nas_prop_bin=factor(floor((nas_prop-0.00001)*10)/10,levels=c(0:9)/10)) %>%
  mutate(counted_sibships=fullsibs+halfsibs,plot=as.character(plot),Year="2019")
sibships_info_2021<-sibships_2021 %>% group_by(family) %>% summarize(fullsibs=sum(sibcall=="full",na.rm=T),halfsibs=sum(sibcall=="half",na.rm=T)) %>% 
  mutate(fraction_halfsibs=halfsibs/(halfsibs+fullsibs)) %>%
  left_join(boricemom_info_2021CAC,by=c("family"="BORICE_FAM")) %>%
  mutate(nas_prop_bin=factor(floor((nas_prop-0.00001)*10)/10,levels=c(0:9)/10)) %>%
  mutate(counted_sibships=fullsibs+halfsibs,plot=as.character(plot),Year="2021")
sibships_info_2022<-sibships_2022 %>% group_by(family) %>% summarize(fullsibs=sum(sibcall=="full",na.rm=T),halfsibs=sum(sibcall=="half",na.rm=T)) %>% 
  mutate(fraction_halfsibs=halfsibs/(halfsibs+fullsibs)) %>%
  left_join(boricemom_info_2022,by=c("family"="BORICE_FAM")) %>%
  mutate(nas_prop_bin=factor(floor((nas_prop-0.00001)*10)/10,levels=c(0:9)/10)) %>%
  mutate(counted_sibships=fullsibs+halfsibs,plot=as.character(plot),Year="2022")

allsibships<-bind_rows(sibships_info_2019,sibships_info_2021,sibships_info_2022) %>% 
  mutate(counted_sibships=factor(counted_sibships),
         nas_prop_cohort=ifelse(nas_prop>0.8,"nas",
                                ifelse(nas_prop>0.5,"admix_high",
                                       ifelse(nas_prop>0.15,"admix","gut"))))
# allsibships %>% ggplot(aes(x=counted_sibships,y=fraction_halfsibs)) + 
#   geom_boxplot(outlier.shape=NA) + geom_beeswarm(aes(color=Year)) + theme_bw() + 
#   scale_color_manual(values=colors_yr[c(2,3,4)])
# 
# allsibships %>% ungroup() %>% group_by(Year,nas_prop_cohort) %>% 
#   summarize(n=n(),contains_halfsibs=sum(fraction_halfsibs>0),
#             count_halfsibs=sum(halfsibs),
#             count_sibships=sum(halfsibs+fullsibs)) %>%
#   mutate(fraction_halfsibs=count_halfsibs/count_sibships) %>% View()
# allsibships %>% filter(fraction_halfsibs>0) %>% 
#   group_by(Year,nas_prop_bin) %>% 
#   summarize(n=n())
# allsibships %>% filter(fraction_halfsibs>0) %>% 
#   group_by(Year,nas_prop_bin) %>% 
#   summarize(n=n())

#readcount info
readcount_info<-read.table("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/Everything_readcount_info.txt",header=T)
readcount_info_fullbatch<-readcount_info %>% filter(str_starts(Sample,"KK|CAC2019|CAC2021|Moms2022|Census2022|Prog2022|LM2021"))
readcount_info_fullbatch %>% summarize(total=sum(Nreadpairs),
                                       mean=mean(Nreadpairs),
                                       median=median(Nreadpairs),
                                       standdev=sqrt(var(Nreadpairs)),
                                       min=min(Nreadpairs),
                                       max=max(Nreadpairs))
readcount_info_fullbatch %>% mutate(Sample=str_remove(Sample,"-combined")) %>% inner_join(allsamples_PCAinfo,by=c("Sample"="sampleID")) %>% 
  filter(total>=50000,!is.na(plot)) %>%
  summarize(nsamples=n(),
            total=sum(Nreadpairs),
            mean=mean(Nreadpairs),
            median=median(Nreadpairs),
            standdev=sqrt(var(Nreadpairs)),
            min=min(Nreadpairs),
            max=max(Nreadpairs))
bamcoverage<-read.table("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Outputs/coverage_info/summary_bam_coverage_CAC2019_2021.txt",
                        header=F,col.names = c("sampleID","nreads_mapped","mean_coverage","std_coverage",
                                               "coverage1X","coverage2X","coverage3X","coverage5X","coverage10X")) %>%
  bind_rows(read.table("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Outputs/coverage_info/summarize_bam_coverage_CAC2022.txt",
                       header=F,col.names = c("sampleID","nreads_mapped","mean_coverage","std_coverage",
                                              "coverage1X","coverage2X","coverage3X","coverage5X","coverage10X")))
readcount_info_with_bamcoverage<-right_join(bamcoverage,readcount_info_fullbatch,by=c("sampleID"="Sample")) %>%
  mutate(sampleID=str_remove(sampleID,"-combined")) %>% filter(!duplicated(sampleID)) %>%
  left_join(allsamples_PCAinfo,by="sampleID")
#write.table(readcount_info_with_bamcoverage,"~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/Allsamples_readcount_bamcoverage_info.txt",quote = F,row.names=F,sep = '\t')

readcount_info_with_bamcoverage %>% ggplot(aes(x=Nreadpairs,nreads_mapped)) + theme_bw() + geom_point() + 
  geom_abline(slope=1,intercept=0,color="red")
readcount_info_with_bamcoverage %>% summarize(mapped=sum(nreads_mapped,na.rm=T),raw=sum(Nreadpairs,na.rm=T)) %>% 
  mutate(mappability=mapped/raw)
readcount_info_with_bamcoverage %>% summarize(n=sum(!is.na(mean_coverage)),mean_coverage=mean(mean_coverage,na.rm=T))
readcount_info_with_bamcoverage %>% mutate(coverage1X=as.numeric(str_remove(coverage1X,"%")),
                                           coverage2X=as.numeric(str_remove(coverage2X,"%"))) %>% 
  summarize(mean_coverage1X=mean(coverage1X,na.rm=T),
            mean_coverage2X=mean(coverage2X,na.rm=T))
readcount_info_with_bamcoverage_useful<-readcount_info_with_bamcoverage %>% 
  filter(sampleID %in% allsamples_PCAinfo$sampleID,!is.na(stream),!is.na(total))
readcount_info_with_bamcoverage_finalset<-readcount_info_with_bamcoverage %>% 
  filter(sampleID %in% allsamples_PCAinfo$sampleID,total>=50000,!is.na(stream))
readcount_info_with_bamcoverage_finalset %>% summarize(mapped=sum(nreads_mapped,na.rm=T),raw=sum(Nreadpairs,na.rm=T)) %>% 
  mutate(mappability=mapped/raw)
readcount_info_with_bamcoverage_finalset %>% summarize(n=sum(!is.na(mean_coverage)),mean_coverage=mean(mean_coverage,na.rm=T))
readcount_info_with_bamcoverage_finalset %>% mutate(coverage1X=as.numeric(str_remove(coverage1X,"%")),
                                           coverage2X=as.numeric(str_remove(coverage2X,"%"))) %>% 
  summarize(mean_coverage1X=mean(coverage1X,na.rm=T),
            mean_coverage2X=mean(coverage2X,na.rm=T))

readcount_info_nasutus_samples<-readcount_info %>% mutate(Sample=str_remove(Sample,"-combined")) %>% 
  inner_join(allsamples_PCAinfo,by=c("Sample"="sampleID")) %>% 
  filter(nas_prop>0.85) %>%
  mutate(rawcoverage_perbase=Nreadpairs*300/430000000) %>% arrange(desc(rawcoverage_perbase))
readcount_info_nasutus_familypool<-readcount_info_nasutus_samples %>% group_by(Year,stream,plot,momfullID) %>% 
  summarize(mom_present=(sum(isMom)>0),n_offspring=sum(!isMom),total_readpairs=sum(Nreadpairs),total_rawcoverage_perbase=sum(rawcoverage_perbase)) %>%
  arrange(desc(total_rawcoverage_perbase))

#write.table(readcount_info_nasutus_samples,"~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Other/Nasutus_popgen/nasutus_samples_info.txt",sep="\t",quote=F,row.names=F)
#write.table(readcount_info_nasutus_familypool,"~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Other/Nasutus_popgen/nasutus_familypools_info.txt",sep="\t",quote=F,row.names=F)

DPRnasutus_samples_list<-DPRfullfilter %>% rename("nas_prop"="prop_p2","total"="count_total") %>% filter(nas_prop>0.85) %>% mutate(proportion_called=total/(208560*2)) %>% arrange(proportion_called) 
readcount_info_DPRnasutus<-readcount_info %>% mutate(Sample=str_remove(Sample,"-combined")) %>% 
  inner_join(DPRnasutus_samples_list,by=c("Sample"="sampleID")) %>% 
  filter(nas_prop>0.85) %>%
  mutate(rawcoverage_perbase=Nreadpairs*300/430000000) %>% arrange(desc(rawcoverage_perbase))
#write.table(readcount_info_DPRnasutus,"~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Other/Nasutus_popgen/DPRnasutus_samples_info.txt",sep="\t",quote=F,row.names=F)

CACnasutus_simple<-readcount_info_nasutus_samples %>% mutate(area="CACarea",sitecode=ifelse(stream=="LM","LM","CAC")) %>% select(Year,area,sitecode,stream,plot,plot_renamed,momfullID,Sample,nas_prop,heterozygosity,total,Nreadpairs,rawcoverage_perbase) %>% arrange(rawcoverage_perbase)
DPRnasutus_simple<-readcount_info_DPRnasutus %>% mutate(Year=2021,area="DPRarea",stream=NA,plot=NA,plot_renamed=NA,momfullID=str_remove(Sample,"_[0-9]+$")) %>% select(Year,area,sitecode,stream,plot,plot_renamed,momfullID,Sample,nas_prop,heterozygosity,total,Nreadpairs,rawcoverage_perbase) %>% arrange(rawcoverage_perbase)
all_nasutus_simple<-bind_rows(CACnasutus_simple,DPRnasutus_simple)
#write.table(all_nasutus_simple,"~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Other/Nasutus_popgen/all_nasutus_info.txt",sep="\t",quote=F,row.names=F)

###get bam lists
moms_to_PCA<-allsamples_PCAinfo %>% 
  filter(isMom,isSOOK %in% c("no","F1"),Year!="2012",total>=50000,!is.na(plot)) %>%
  select(bampath)
moms_to_PCA_CACadmixedonly<-allsamples_PCAinfo %>% 
  filter(isMom,isSOOK %in% c("no","F1"),Year!="2012",total>=50000,!is.na(plot),
         stream!="LM",nas_prop<0.85) %>%
  select(bampath)

#write.table(moms_to_PCA,"~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Outputs/moms_to_PCA_bamlist.txt",sep="\t",quote=F,row.names=F)
#write.table(moms_to_PCA_CACadmixedonly,"~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Outputs/moms_to_PCA_CACadmixedonly_bamlist.txt",sep="\t",quote=F,row.names=F)

##panel of 100 representative
#representative100<-allsamples_PCAinfo %>% 
#  filter(isMom,isSOOK %in% c("no","F1"),Year!="2012",total>=200000,!is.na(plot)) %>% 
#  slice_sample(n=100)
representative100 %>% mutate(nas_prop_bin=factor(ceiling(nas_prop*20)/20)) %>% 
  group_by(nas_prop_bin) %>% summarize(n())
representative100_bamlist<-representative100 %>% select(bampath)

#write.table(representative100_bamlist,"~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Outputs/representative100_bamlist.txt",sep="\t",quote=F,row.names=F)

