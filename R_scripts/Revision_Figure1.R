
####
#### Final Figure 1 plots
####

#GPS data
site_locations<-read_excel("~/OneDrive - University of Georgia/Sweigart_Lab/Ch2_Natural_hybrids/Data/Datasheets/CAC-site-locations.xlsx")

allplots_GPS<-read_excel("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Datasheets/plot_GPS.xlsx") %>%
  mutate(plot_renamed=sapply(plot,get_univ_plot,key=plot_key,zone="renamed")) %>%
  separate(plot_renamed,into=c("str","plot_renamed_simple"),sep="_|M",remove = F)

LM_GPS<-allplots_GPS %>% filter(stream=="LM")
CAC_GPS<-allplots_GPS %>% filter(stream!="LM")

#libraries
library(geodata)
library(sf)
library(tidyterra)
library(ggmap)
library(ggrepel)
library(ggspatial)
library(patchwork)

### PRISM precipitation and temperature data (panel 1B and 1C)
####1B: Temp and Precipitation Data

PRISM_daily_precip<-PRISM_daily %>% mutate(Year=factor(as.character(Year),levels=as.character(c(2010:2022)))) %>%
  filter(Year %in% as.character(c(2012:2022))) %>%
  ggplot(aes(x=Julian_date,y=rollsum_14day_precip_mm,group=Year,color=Year)) +
  geom_line() +
  scale_color_manual(values=c(colors_yr[1],rep("grey",6),colors_yr[2],"grey",colors_yr[3:4]),guide="none") +
  #scale_linetype_manual(values=c("solid",rep("dotted",6),"solid","dotted",rep("solid",2)),guide="none") +
  scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),labels=c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  ylab("Precipitation (mm)") + xlab("Month") +
  geom_vline(xintercept=91,linetype="dashed",color="black") +
  geom_vline(xintercept=182,linetype="dashed",color="black") +
  theme_bw() +
  theme(text=element_text(size=unit(8,units="pt")),
        axis.title=element_text(face="bold"))

PRISM_daily_temp<-PRISM_daily %>% mutate(Year=factor(as.character(Year),levels=as.character(c(2010:2022)))) %>%
  filter(Year %in% as.character(c(2012:2022))) %>%
  ggplot(aes(x=Julian_date,y=rollmean_14day_tmean_degreesC,group=Year,color=Year)) +
  geom_line() +
  scale_color_manual(values=c(colors_yr[1],rep("grey",6),colors_yr[2],"grey",colors_yr[3:4]),guide="none") +
  #scale_linetype_manual(values=c("solid",rep("dotted",6),"solid","dotted",rep("solid",2)),guide="none") +
  scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),labels=c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  ylab("Temperature (C)") + xlab("Month") +
  geom_vline(xintercept=91,linetype="dashed",color="black") +
  geom_vline(xintercept=182,linetype="dashed",color="black") +
  theme_bw() +
  theme(text=element_text(size=unit(8,units="pt")),
        axis.title=element_text(face="bold"))

###State and province outlines for panel 1A
us <- gadm(country="USA",level=1,path = "~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch3Real_genomescans/USA_outline")
canada <- gadm(country="CAN",level=1,path = "~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch3Real_genomescans/CAN_outline")

us.states<-us[us$NAME_1 %in% c("California","Oregon","Washington","Idaho","Nevada","Arizona","Montana","Utah","Wyoming"),]
can.states<-canada[canada$NAME_1 %in% c("British Columbia","Alberta"),]

###Top half panel 1A
pnw<-ggplot(us.states) +
  geom_spatvector(data=us.states.plus,fill=NA) +
  geom_spatvector(data=can.states,fill=NA) + theme_bw() +
  coord_sf(xlim=c(-125,-115),ylim=c(42,50),crs=4326) +
  #geom_point(data=cac_location,fill=colors_sp[3],shape=15,aes(group=NULL),size=5) +
  geom_rect(xmin=-121.405,xmax=-121.355,ymin=45.69,ymax=45.72,color="black",fill="transparent") +
  #geom_point(aes(x=longitude,y=latitude),data=site_locations,group=NA,size=1) +
  scale_y_continuous(breaks=c(42,44,46,48,50),name="Latitude",expand = c(0,0)) +
  scale_x_continuous(breaks=c(-124,-120,-116),name="Longitude",expand = c(0,0)) +
  theme(text=element_text(size=unit(8,units="pt")),
        axis.title=element_text(face="bold"),
        axis.text.x=element_text(angle=45,hjust=1))
pnw


###transparent map outlines with site info
columbia_transparent<-ggplot(sites_to_draw,aes(x=longitude,y=latitude)) +
  coord_sf(crs=4326,default_crs = sf::st_crs(4326)) +
  annotation_scale(text_cex=0) +
  geom_point(fill=colors_sp[3],color="black",shape=21,size=2) +
  scale_y_continuous(limits=c(45.69,45.72),breaks=c(45.69,45.7,45.71,45.72),name="Latitude",expand = c(0,0)) +
  scale_x_continuous(limits=c(-121.405,-121.355),breaks=c(-121.40,-121.38,-121.36),name="Longitude",expand = c(0,0)) +
  theme(text=element_text(size=unit(8,units="pt")),
        axis.text.x=element_text(angle=45,hjust=1),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.title=element_text(face="bold"))
CAC_transparent<-ggplot(CAC_GPS,aes(x=longitude,y=latitude,label=plot_renamed_simple)) +
  coord_sf(crs=4326,default_crs = sf::st_crs(4326)) +
  annotation_scale() +
  scale_y_continuous(limits=c(45.710,45.713),breaks=c(45.710,45.711,45.712,45.713),name="Latitude",expand = c(0,0)) +
  scale_x_continuous(limits=c(-121.36775,-121.36425),breaks=c(-121.367,-121.366,-121.365),name="Longitude",expand = c(0,0)) +
  geom_label(size=4) +
  theme(text=element_text(size=unit(8,units="pt")),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.title=element_text(face="bold"))
LM_transparent<-ggplot(LM_GPS,aes(x=longitude,y=latitude,label=plot_renamed_simple)) +
  coord_sf(crs=4326,default_crs = sf::st_crs(4326)) +
  annotation_scale() +
  scale_y_continuous(limits=c(45.702,45.705),breaks=c(45.702,45.703,45.704,45.705),name="Latitude",expand = c(0,0)) +
  scale_x_continuous(limits=c(-121.39625,-121.39275),breaks=c(-121.396,-121.395,-121.394,-121.393),name="Longitude",expand = c(0,0)) +
  geom_label(size=4) +
  theme(text=element_text(size=unit(8,units="pt")),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.title=element_text(face="bold"))

fullmap_transparent<-pnw+PRISM_daily_precip+
  columbia_transparent+PRISM_daily_temp+
  LM_transparent+CAC_transparent+
  plot_layout(design="AABBBB\nCCDDDD\nEEEFFF",heights=c(1,1,2))

ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/JAN_fixes/fullmaps_transparent.png",fullmap_transparent,device = "png",height=8,width=8,dpi=600,units="in",bg = "transparent")

###OLD version using Google maps

CACgooglemap<-get_map(location = c(left=-121.405,bottom=45.690,right=-121.355,top=45.720),maptype="satellite",source="google",api_key=my_key)
sites_to_draw<-site_locations %>% filter(Site %in% c("LM","CAC_stream1","CAC_stream2"))

columbia<-ggmap(CACgooglemap) +
  coord_sf(crs=4326,default_crs = sf::st_crs(4326)) +
  #annotation_scale(text_cex=0,bar_cols = c("white","white"),height = unit(0.05,"in"),pad_x = unit(0.1,"npc"),pad_y=unit(0.1,"npc"),) +
  annotation_scale(text_cex=0,bar_cols = c("white","white")) +
  #scale_x_continuous(expand = c(0,0),name="Longitude") +
  geom_point(data=sites_to_draw,aes(x=longitude,y=latitude),color=colors_sp[3],size=1) +
  #geom_label_repel(data=site_locations,aes(x=longitude,y=latitude,label=Site)) +
  #geom_rect(xmin=-121.368,xmax=-121.364,ymin=45.709,ymax=45.713,color="white",fill="transparent") +
  scale_y_continuous(breaks=c(45.69,45.7,45.71,45.72),name="Latitude",expand = c(0,0)) +
  scale_x_continuous(breaks=c(-121.40,-121.38,-121.36),name="Longitude",expand = c(0,0)) +
  theme(axis.text.x=element_text(angle=45,hjust=1))
columbia


dualmap<-pnw + columbia + plot_layout(design="
A
B
")

map_plus_PRISMs<-pnw+PRISM_daily_precip+columbia+PRISM_daily_temp+
  plot_layout(design="AB\nCD",widths = c(1,3))
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/DEC_revision/pnw_map_plus_PRISMs.pdf",map_plus_PRISMs,device = "pdf",height=4,width=8,dpi=600,units="in")

###GOOGLE maps of each stream
CACsitegooglemap<-get_map(location = c(left=-121.368,bottom=45.708,right=-121.364,top=45.716),maptype="satellite",source="google",api_key=my_key)
CACsitegooglemap2<-get_map(location = c(left=-121.368,bottom=45.709,right=-121.364,top=45.713),maptype="satellite",source="google",api_key=my_key)

LMsitegooglemap1<-get_map(location = c(left=-121.396,bottom=45.700,right=-121.393,top=45.708),maptype="satellite",source="google",api_key=my_key)
LMsitegooglemap2<-get_map(location = c(left=-121.396,bottom=45.701,right=-121.393,top=45.705),maptype="satellite",source="google",api_key=my_key)

CACstreams<-ggmap(CACsitegooglemap) + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(labels = ~ .x,breaks = c(-121.367,-121.366,-121.365),expand = c(0,0),name="Longitude") +
  scale_y_continuous(labels = ~ .x,expand = c(0,0),name="Latitude") +
  #geom_label(data=CAC_GPS,aes(x=longitude,y=latitude,label=plot_renamed_simple))
  geom_label(data=filter(CAC_GPS,latitude>=45.712),aes(x=longitude,y=latitude,label=plot_renamed_simple),size=4)
CACstreams2<-ggmap(CACsitegooglemap2) + ylab(NULL) + xlab(NULL) +
  coord_sf(crs=4326,default_crs = sf::st_crs(4326)) + annotation_scale(pad_y=unit(0.4,"in"),text_cex=0) +
  scale_x_continuous(labels = ~ .x,breaks = c(-121.367,-121.366,-121.365),expand = c(0,0),name="Longitude") +
  scale_y_continuous(labels = ~ .x,expand = c(0,0),name="Latitude") +
  #geom_label(data=CAC_GPS,aes(x=longitude,y=latitude,label=plot_renamed_simple))
  geom_label(data=filter(CAC_GPS,latitude<45.712),aes(x=longitude,y=latitude,label=plot_renamed_simple),size=4)
CACstreams
CACstreams2

LMstream_map1<-ggmap(LMsitegooglemap1) + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(labels = ~ .x,breaks = c(-121.396,-121.395,-121.394,-121.393),expand = c(0,0),name="Longitude") +
  scale_y_continuous(labels = ~ .x,expand = c(0,0),name="Latitude") +
  #geom_label(data=LM_GPS,aes(x=longitude,y=latitude,label=plot_renamed_simple))
  geom_label(data=filter(LM_GPS,latitude>=45.703),aes(x=longitude,y=latitude,label=plot_renamed_simple),size=4)
LMstream_map2<-ggmap(LMsitegooglemap2) + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(labels = ~ .x,breaks = c(-121.396,-121.395,-121.394,-121.393),expand = c(0,0),name="Longitude") +
  scale_y_continuous(labels = ~ .x,expand = c(0,0),name="Latitude") +
  #geom_label(data=LM_GPS,aes(x=longitude,y=latitude,label=plot_renamed_simple)) +
  geom_label(data=filter(LM_GPS,latitude<=45.703),aes(x=longitude,y=latitude,label=plot_renamed_simple),size=4) +
  coord_sf(crs=4326,default_crs = sf::st_crs(4326)) + annotation_scale(pad_y=unit(0.4,"in"),text_cex=0)

fourmaps<-plot_grid(LMstream_map1,CACstreams,LMstream_map2,CACstreams2,nrow=2)
fourmaps
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/DEC_revision/streams_fourmaps.pdf",fourmaps,device = "pdf",height=8,width=8,dpi=600,units="in")


###NOT USED: Natural earth data + open street map, open source map features

library(rnaturalearthdata)
library(rnaturalearth)
library(osmdata)
get_ne_rivers<-ne_download(category="physical",type="rivers_north_america",scale=10)
get_ne_roads<-ne_download(category="cultural",type="roads_north_america",scale=10)
get_ne_parks<-ne_download(category="cultural",type="parks_and_protected_lands_point",scale=10)
get_osm_river<-opq(bbox=c(-121.405,45.69,-121.355,45.72)) %>%
  add_osm_feature(key="natural",value="water") %>%
  osmdata_sf()
get_osm_streams<-opq(bbox=c(-121.405,45.69,-121.355,45.72)) %>%
  add_osm_feature(key="waterway",value="stream") %>%
  osmdata_sf()
get_osm_reserve<-opq(bbox=c(-121.405,45.69,-121.355,45.72)) %>%
  add_osm_feature(key="leisure",value="nature_reserve") %>%
  osmdata_sf()
get_osm_grassland<-opq(bbox=c(-121.405,45.69,-121.355,45.72)) %>%
  add_osm_feature(key="natural",value="grassland") %>%
  osmdata_sf()
get_osm_roads<-opq(bbox=c(-121.405,45.69,-121.355,45.72)) %>%
  add_osm_features(features = c ("\"highway\"=\"motorway\"",
                                 "\"highway\"=\"trunk\"",
                                 "\"highway\"=\"primary\"",
                                 "\"highway\"=\"secondary\"",
                                 "\"highway\"=\"tertiary\"",
                                 "\"highway\"=\"motorway_link\"",
                                 "\"highway\"=\"trunk_link\"",
                                 "\"highway\"=\"primary_link\"",
                                 "\"highway\"=\"secondary_link\"",
                                 "\"highway\"=\"tertiary_link\""
  )) %>%
  osmdata_sf()
get_osm_rails<-opq(bbox=c(-121.405,45.69,-121.355,45.72)) %>%
  add_osm_feature(key="railway",value="rail") %>%
  osmdata_sf()
get_osm_treecover<-opq(bbox=c(-121.405,45.69,-121.355,45.72)) %>%
  add_osm_features(features = c ("\"natural\"=\"tree\"",
                                 "\"natural\"=\"wood\"",
                                 "\"natural\"=\"scrub\"")) %>%
  osmdata_sf()

osm_basemap<-ggplot(get_osm_river$osm_multipolygons) +
  coord_sf(crs=4326,default_crs = sf::st_crs(4326)) +
  scale_y_continuous(limits=c(45.69,45.72),breaks=c(45.69,45.7,45.71,45.72),name="Latitude",expand = c(0,0)) +
  scale_x_continuous(limits=c(-121.405,-121.355),breaks=c(-121.40,-121.38,-121.36),name="Longitude",expand = c(0,0)) +
  #geom_spatvector(data=get_ne_rivers,color="blue",fill="blue") +
  #geom_spatvector(data=get_ne_roads,color="red") +
  geom_rect(xmin=-121.405,xmax=-121.355,ymin=45.69,ymax=45.72,fill="tan") +
  geom_sf(data=get_osm_reserve$osm_multipolygons,fill="lightgreen",color="black") +
  geom_sf(data=get_osm_reserve$osm_polygons,fill="lightgreen",color="black") +
  geom_sf(data=get_osm_roads$osm_lines,fill="darkgrey",color="black") +
  geom_sf(data=get_osm_roads$osm_multilines,fill="darkgrey",color="black") +
  geom_sf(data=get_osm_rails$osm_lines,fill="lightgrey",color="black") +
  geom_sf(data=get_osm_rails$osm_multilines,fill="lightgrey",color="black") +
  geom_sf(data=get_osm_treecover$osm_multipolygons,fill="darkgreen",color="black") +
  geom_sf(data=get_osm_treecover$osm_polygons,fill="darkgreen",color="black") +
  geom_sf(data=get_osm_river$multipolygons,fill="blue",color="black") +
  geom_sf(data=get_osm_river$osm_polygons,fill="blue",color="black") +
  geom_sf(data=get_osm_streams$osm_lines,fill="blue",color="blue") +
  geom_sf(data=get_osm_streams$osm_multilines,fill="blue",color="blue") +
  theme(text=element_text(size=unit(8,units="pt")))
threestreams<-osm_basemap +
  geom_point(data=sites_to_draw,aes(x=longitude,y=latitude),shape=21,fill=colors_sp[3],color="black",size=3) +
  annotation_scale(text_cex=0) +
  theme(text=element_text(size=unit(8,units="pt")),
        axis.text=element_text(face="bold"))

CAC_osm<-osm_basemap +
  annotation_scale() +
  scale_y_continuous(limits=c(45.710,45.713),breaks=c(45.710,45.711,45.712,45.713),name="Latitude",expand = c(0,0)) +
  scale_x_continuous(limits=c(-121.368,-121.364),breaks=c(-121.367,-121.366,-121.365),name="Longitude",expand = c(0,0)) +
  geom_label(data=CAC_GPS,aes(x=longitude,y=latitude,label=plot_renamed_simple),size=4) +
  theme(text=element_text(size=unit(8,units="pt")),
        axis.text=element_text(face="bold"))

LM_osm<-osm_basemap +
  annotation_scale() +
  scale_y_continuous(limits=c(45.702,45.705),breaks=c(45.702,45.703,45.704,45.705),name="Latitude",expand = c(0,0)) +
  scale_x_continuous(limits=c(-121.3965,-121.3925),breaks=c(-121.396,-121.395,-121.394,-121.393),name="Longitude",expand = c(0,0)) +
  geom_label(data=LM_GPS,aes(x=longitude,y=latitude,label=plot_renamed_simple),size=4) +
  theme(text=element_text(size=unit(8,units="pt")),
        axis.text=element_text(face="bold"))

osmfullfig<-pnw+PRISM_daily_precip+
  threestreams+PRISM_daily_temp+
  LM_osm+CAC_osm+
  plot_layout(design="AABBBB\nCCDDDD\nEEEFFF",heights=c(1,1,2))
