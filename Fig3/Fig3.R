rm(list=ls())

#time to make Fig3A

library(ggplot2)
#install.packages('ggtree')
library(ggtree)
library(ape)
library(ggnewscale)
library(tidyverse)

homewd= "/Users/SophiaHorigan/Documents/GitHub/Mada-Bat-Astro/"

setwd(paste0(homewd, "/Fig3"))

#load the fig3a tree
treeA <-  read.tree(file = paste0(homewd, "Fig3/mam_fg_500bs.newick"))

#root it


rooted.tree.A <- root(treeA, which(treeA$tip.label == "NC_002470"))
#take a quick look in base R
plot(rooted.tree.A)

#load tree data prepared from elsewhere
dat <- read.csv(file=paste0(homewd,"Fig3/astro_meta_full.csv"), header = T, stringsAsFactors = F)
head(dat)
#dat <- dat[-c(10:17)]
#check subgroup names
unique(dat$Genus)

colz = c("Mamastrovirus" = "royalblue", "Avastrovirus" = "tomato")

#pick order for the labels
dat$Genus <- factor(dat$Genus, levels = c("Mamastrovirus", "Avastrovirus"))   

#and add a "novel" category
dat$novel = 0
dat$novel[dat$Geo_Location=="Madagascar"] <- 1
dat$novel <- as.factor(dat$novel)

#rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)

#take a glance
p <- ggtree(rooted.tree.A) %<+% dat + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = c(.2,.85), legend.title = element_blank())
p #why are some of the tips grey?

#now get new tip labels
dat$old_tip_label <- dat$Accession
dat$new_label <- NA
dat$new_label[!is.na(dat$Strain)] <- paste(dat$Accession[!is.na(dat$Strain)], " | ", 
                                           dat$Strain[!is.na(dat$Strain)], " | ", 
                                           dat$Host[!is.na(dat$Strain)], " | ",
                                           dat$Geo_Location[!is.na(dat$Strain)], " | ",
                                           dat$Year[!is.na(dat$Strain)])

#dat$new_label[is.na(dat$strain)] <- paste(dat$accession_num[is.na(dat$strain)], " | ", 
#                                         dat$host[is.na(dat$strain)], " | ",
#                                          dat$country[is.na(dat$strain)], " | ",
#                                          dat$collection_year[is.na(dat$strain)])

#after Gwen checks these, can later manually edit any that are "NA"

#make sure to sort in order
tree.dat <- data.frame(old_tip_label=rooted.tree.A$tip.label, num =1:length(rooted.tree.A$tip.label))
head(tree.dat)
head(dat)

#for some reason there are quotations around our new samples, need to remove
#SEE IF YOU CAN FIX THIS IN THE EXPORT FROM GENEIOUS
## NEED TO SEE WHERE THE NAMES END UP IN TREE.DAT
tree.dat[42,1] = "F_MIZ141_RR034B_198_NODE_4_length_6593_cov_808.543744"
tree.dat[43,1] = "F_MIZ141_RR034B_198_NODE_5_length_6456_cov_49.595532"

tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = T, sort = F)

names(tree.dat)

tree.dat$tip_label <- tree.dat$new_label
tree.dat <- dplyr::select(tree.dat, tip_label, Accession, Strain, Host, Bat_host, Geo_Location, Year, Genus, novel, old_tip_label, Family, Animal)

rooted.tree.A$tip.label <- tree.dat$tip_label

tree.dat$Bat_host[tree.dat$Bat_host==0] <- "non-bat host"
tree.dat$Bat_host[tree.dat$Bat_host==1] <- "bat host"
tree.dat$Bat_host <- as.factor(tree.dat$Bat_host)
shapez = c("bat host" =  24, "non-bat host" = 21)
colz2 = c('1' =  "yellow", '0' = "white")

## Tree 1: Colored by genus, shape by bat/non-bat
p1 <- ggtree(rooted.tree.A) %<+% tree.dat + geom_tippoint(aes(color=Genus, fill=Genus, shape=Bat_host)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  scale_color_manual(values=colz) + 
  scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.85), legend.title = element_blank()) +
  xlim(c(0,5.6))
p1
########FAMILY#############
## color tip by family name
p1_family <- ggtree(rooted.tree.A) %<+% tree.dat + geom_tippoint(aes(color=Family, fill=Family, shape=Bat_host)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.6), legend.title = element_blank()) +
  xlim(c(0,6))
p1_family

#collapsing branches by family where possible
##need to match colors to color bar, or just remove color bar and label on the figure
node1 <- MRCA(p1_family, 'NC_024297  |  Bovine astrovirus  |  Bos taurus  |  China  |  2014', 'NC_023629  |  Bovine astrovirus B76/HK  |  Bos taurus  |  Hong Kong  |  NA')
c1 <- scaleClade(p1_family, node1, .3) %>% collapse(node1, 'max', fill='red')

node2 <- MRCA(p1_family, 'NC_018702  |  Murine astrovirus  |  Mouse  |  USA  |  2011', 'NC_036583  |  Rodent astrovirus  |  Rattus norvegicus  |  China  |  NA')
c2 <- collapse(c1, node2, 'max', fill='green')

node3 <- MRCA(p1_family, 'NC_023675  |  Porcine astrovirus 4  |  Sus scrofa  |  USA  |  2010','NC_016896  |  Astrovirus wild boar/WBAstV-1/2011/HUN  |  Sus scrofa  |  Hungary  |  2011')
c3 <- collapse(c2, node3, 'max', fill='purple')

node4 <- MRCA(p1_family, 'NC_019027  |  Astrovirus VA4  |  Homo sapiens  |  Nepal  |  2008', 'NC_024472  |  Burkina Faso astrovirus  |  Homo sapiens  |  Burkina Faso  |  2010')
c4 <- collapse(c3, node4, 'max', fill='darkgreen')

scaleClade(c4, node4, .3)

##try highlighting and shrinking instead
h1 <- p1_family + geom_highlight(node1, 'red')

## color tip by common name
p1_common <- ggtree(rooted.tree.A, size = 1) %<+% tree.dat + geom_tippoint(aes(color=Animal, fill=Animal, shape=Bat_host, size=1)) +
  geom_nodelab(size=1.8,nudge_x = -0.11, nudge_y = .4) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3,size=3, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.05,.6), legend.title = element_blank()) +
  xlim(c(0,10.5))
p1_common


ggsave(file = paste0(homewd, "/final-figures/Fig3_poster_4.png"),
       plot = p1_family,
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=4)#, 






#########################
##### And Fig 3 B #####
#########################

#load the fig3a tree
treeB <-  read.tree(file = paste0(homewd, "Fig3/mam_rdrp_100bs.newick"))

#remove quotes
treeB$tip.label <- gsub("'", '', treeB$tip.label)

#root it
rooted.tree.B <- root(treeB, which(treeB$tip.label == "NC_002470 - RdRp"))

#take a quick look in base R
plot(rooted.tree.B)

#load tree data prepared from elsewhere
datB <- read.csv(file=paste0(homewd,"Fig3/astro_meta_rdrp.csv"), header = T, stringsAsFactors = F)
head(datB)

#check subgroup names
unique(datB$Genus)

colz = c("Mamastrovirus" = "royalblue", "Avastrovirus" = "tomato")

#pick order for the labels
datB$Genus <- factor(datB$Genus, levels = c("Mamastrovirus", "Avastrovirus"))   

#and add a "novel" category
datB$novel = 0
datB$novel[datB$Geo_Location=="Madagascar"] <- 1
datB$novel <- as.factor(datB$novel)

#rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)

#take a glance
p2 <- ggtree(rooted.tree.B) %<+% datB + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = c(.2,.85), legend.title = element_blank())
p2 #why are some of the tips grey?

#now get new tip labels
datB$old_tip_label <- datB$Accession
datB$new_label <- NA
datB$new_label[!is.na(datB$Strain)] <- paste(datB$Accession[!is.na(datB$Strain)], " | ", 
                                           datB$Strain[!is.na(datB$Strain)], " | ", 
                                           datB$Host[!is.na(datB$Strain)], " | ",
                                           datB$Geo_Location[!is.na(datB$Strain)], " | ",
                                           datB$Year[!is.na(datB$Strain)])

#dat$new_label[is.na(dat$strain)] <- paste(dat$accession_num[is.na(dat$strain)], " | ", 
#                                         dat$host[is.na(dat$strain)], " | ",
#                                          dat$country[is.na(dat$strain)], " | ",
#                                          dat$collection_year[is.na(dat$strain)])

#after Gwen checks these, can later manually edit any that are "NA"

#make sure to sort in order
tree.datB <- data.frame(old_tip_label=rooted.tree.B$tip.label, num =1:length(rooted.tree.B$tip.label))
head(tree.datB)
head(datB)

tree.datB <- merge(tree.datB, datB, by = "old_tip_label", all.x = T, sort = F)

names(tree.datB)

tree.datB$tip_label <- tree.datB$new_label
tree.datB <- dplyr::select(tree.datB, tip_label, Accession, Strain, Host, Bat_host, Geo_Location, Year, Genus, novel, old_tip_label, Family, Animal)

rooted.tree.B$tip.label <- tree.datB$tip_label

tree.datB$Bat_host[tree.datB$Bat_host==0] <- "non-bat host"
tree.datB$Bat_host[tree.datB$Bat_host==1] <- "bat host"
tree.datB$Bat_host <- as.factor(tree.datB$Bat_host)
shapez = c("bat host" =  24, "non-bat host" = 21)
colz2 = c('1' =  "yellow", '0' = "white")

## Tree 1: Colored by genus, shape by bat/non-bat
pb <- ggtree(rooted.tree.B) %<+% tree.datB + geom_tippoint(aes(color=Genus, fill=Genus, shape=Bat_host)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  scale_color_manual(values=colz) + 
  scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.85), legend.title = element_blank()) +
  xlim(c(0,5.6))
pb

## color tip by family name
pB_family <- ggtree(rooted.tree.B) %<+% tree.datB + geom_tippoint(aes(color=Family, fill=Family, shape=Bat_host)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.6), legend.title = element_blank()) +
  xlim(c(0,5.6))
pB_family

## color tip by common name
pB_common <- ggtree(rooted.tree.B, size = 1) %<+% tree.datB + geom_tippoint(aes(color=Animal, fill=Animal, shape=Bat_host, size=1)) +
  geom_nodelab(size=1.8,nudge_x = -0.11, nudge_y = .4) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3,size=3, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.05,.6), legend.title = element_blank()) +
  xlim(c(0,10.5))
pB_common


ggsave(file = paste0(homewd, "/final-figures/Fig3B_poster_3.png"),
       plot = pB_common,
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=4)#, 



#####now combine the two together somehow


###wrking on p1

p1 <- ggtree(rooted.tree.A) %<+% tree.dat + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), size=3, show.legend = F) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = .9) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3,  show.legend=F, size=3, hjust = -.08) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  geom_treescale(fontsize=4, x=1,y=124, linesize = .5) + 
  xlim(c(0,4))
p1

#and flip some clades
node_flip_Embeco_Merbeco1 = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.A$tip.label == "NC_006213  |  HCoV_OC43  |  Homo_sapiens  |  USA  |  1960"  ))
node_flip_Merbeco_Sarbeco1 = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.A$tip.label == "MZ081380  |  SARSr_CoV  |  Rhinolophus_stheno  |  China  |  2020"))
node_flip_Embeco_Sarbeco1 = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "MZ081380  |  SARSr_CoV  |  Rhinolophus_stheno  |  China  |  2020" ),which(rooted.tree.A$tip.label == "NC_006213  |  HCoV_OC43  |  Homo_sapiens  |  USA  |  1960"  ))

p1.2 <- p1 %>% ggtree::rotate(node = node_flip_Embeco_Sarbeco1 )


#p1.2 <- p1 %>% ggtree::rotate(node = node_flip_Merbeco_Sarbeco1)
#p1.3 <- p1.2 %>% ggtree::rotate(node = node_flip_Embeco_Sarbeco1)


#collapse the alpha clade (all bat CoVs)
#alpha_node = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "NC_048211 | Suncus_murinus | China | 2015" ),which(rooted.tree.A$tip.label == "NC_018871 | bat | China | 2021" ))

p1.2.leg <- ggtree(rooted.tree.A) %<+% tree.dat + 
  geom_tippoint(aes(color=sub_group, shape=bat_host), size=3) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = .9) +
  scale_fill_manual(values=colz) + 
  scale_color_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3,  show.legend=F, size=3, hjust = -.08) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size=12)) +
  geom_treescale(fontsize=4, x=1,y=124, linesize = .5) + 
  xlim(c(0,4))
p1.2.leg

#separate legend
leg.all <- cowplot::get_legend(p1.2.leg)


#new p2
p2.1 <- ggtree(rooted.tree.B) %<+% datB + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), show.legend = F, size=5) +
  geom_nodelab(size=3,nudge_x = -.03, nudge_y = .4) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.1) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "right", legend.title = element_blank()) +
  geom_treescale(fontsize=4, x=.3,y=50, linesize = .5) + 
  xlim(c(0,1.5))
p2.1 

#add lineage clade labels bars

#nodebase
clade.a <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "EF065514  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005" ),which(rooted.tree.B$tip.label == "HM211098  |  HKU9  |  Rhinolophus_sinicus  |  China  |  2005" ))
clade.b <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "OK020089  |  Rousettus_madagascariensis  |  Madagascar  |  2018" ),which(rooted.tree.B$tip.label == "MG693172  |  Eidolon_helvum  |  Cameroon  |  2013" ))
clade.c <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_030886  |  GCCDC1  |  Rousettus_leschenaulti  |  China  |  2014" ),which(rooted.tree.B$tip.label == "MT350598  |  GCCDC1  |  Eonycteris_spelaea  |  Singapore  |  2016" ))
clade.d <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "MK211379  |  GX2018  |  Rhinolophus_affinis  |  China  |  2016" ),which(rooted.tree.B$tip.label == "MK492263  |  BatCoV92  |  Cynopteris_brachyotis  |  Singapore  |  2015" ))
clade.e <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "OK020087  |  Pteropus_rufus  |  Madagascar  |  2018" ),which(rooted.tree.B$tip.label == "OK067319  |  Pteropus_rufus  |  Madagascar  |  2018" ))




p2.1 <- ggtree(rooted.tree.B) %<+% datB + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), show.legend = F, size=5) +
  geom_nodelab(size=3,nudge_x = -.03, nudge_y = .4) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.1) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "right", legend.title = element_blank()) +
  geom_treescale(fontsize=4, x=.3,y=50, linesize = .5) + 
  geom_cladelabel(node = clade.a, label = "HKU9", offset = 1, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.b, label = "atop(African,italic(Eidolon))", offset = .8, fontsize = 6.5, color="tomato", parse=T) +
  geom_cladelabel(node = clade.c, label = "GCCDC1", offset = 1.05, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.d, label = "BtCoV92 /\nGX2018", offset = 1.03, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.e, label = "atop(Madagascar,italic(Pteropus))", offset = .7, fontsize = 6.5, color="tomato" , parse = T) +
  xlim(c(0,1.8))
p2.1 


#great, now need to flip some of the clases to match plot on the left

node_flip_Embeco_Merbeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.B$tip.label == "NC_006213  |  HCoV_OC43  |  Homo_sapiens  |  USA  |  1960"  ))
node_flip_Sarbeco_Hibeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_025217  |  Hipposideros_pratti  |  China  |  2013" ),which(rooted.tree.B$tip.label == "NC_004718  |  SARS_CoV  |  Homo_sapiens  |  Canada  |  2003" ))
node_flip_Embeco_Nobeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.B$tip.label == "EF065516  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005"   ))


p2.2 <- p2.1 %>% ggtree::rotate(node = node_flip_Embeco_Merbeco)
p2.3 <- p2.2 %>% ggtree::rotate(node = node_flip_Sarbeco_Hibeco)
p2.4 <- p2.3 %>% ggtree::rotate(node = node_flip_Embeco_Nobeco)


Fig3 <- cowplot::plot_grid(p1.2,p2.4, ncol=2, nrow=1, labels = c("(A)", "(B)"), label_size = 22, label_x = .03, label_y = .98)

Fig3all <- cowplot::plot_grid(Fig3,leg.all, ncol=1, nrow=2, rel_heights = c(1,.1))


#and save to the final figures

ggsave(file = paste0(homewd, "/final-figures/Fig3.png"),
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=4)#, 
