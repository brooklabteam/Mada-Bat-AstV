rm(list=ls())

#time to make Fig3A
#install.packages('ggplot2')
library(ggplot2)
#install.packages('ggtree')
library(ggtree)
library(ape)
#install.packages('ape')
#nstall.packages('ggnewscale')
library(ggnewscale)
#install.packages('tidyverse')
library(tidyverse)
library(dplyr)
#install.packages('dplyr')
library("ggsci")


homewd= "/Users/sophiahorigan/Documents/GitHub/Mada-Bat-Astro/"

setwd(paste0(homewd, "/ML-trees"))


#########################
##### ML Full Genome #####
#########################

#load the fig3a tree
treeA <-  read.tree(file = paste0(homewd, "ML-trees/full-ML-tree-FPB"))

#root it

rooted.tree.A <- root(treeA, which(treeA$tip.label == "NC_002470"))
#take a quick look in base R
plot(rooted.tree.A)

#load tree data prepared from elsewhere
dat <- read.csv(file=paste0(homewd,"ML-trees/astro_meta_full.csv"), header = T, stringsAsFactors = F)
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
#dat$novel <- as.factor(dat$novel)

#rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)

#take a glance
p <- ggtree(rooted.tree.A) %<+% dat + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = c(.2,.85), legend.title = element_blank())
p #why are some of the tips grey?

#now get new tip labels
dat$old_tip_label <- dat$Accession
dat$new_label <- NA
#dat$new_label[!is.na(dat$Strain)] <- paste(dat$Accession[!is.na(dat$Strain)], " | ", 
                                           #dat$Strain[!is.na(dat$Strain)], " | ", 
                                          # dat$Host[!is.na(dat$Strain)], " | ",
                                          # dat$Geo_Location[!is.na(dat$Strain)], " | ",
                                          # dat$Year[!is.na(dat$Strain)])

dat$new_label[!is.na(dat$Strain)] <- paste(dat$Accession[!is.na(dat$Strain)], " | ", 
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
#tree.dat[5,1] = "F_MIZ141_1"
#tree.dat[6,1] = "F_MIZ141_2"
#tree.dat[3,1] = "NC_002470"

tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = T, sort = F)

names(tree.dat)

tree.dat$tip_label <- tree.dat$new_label
tree.dat <- dplyr::select(tree.dat, tip_label, Accession, Strain, Host, Bat_host, Geo_Location, Year, Genus, novel, old_tip_label, Family, Animal, Order)

rooted.tree.A$tip.label <- tree.dat$tip_label

#shape by novelty
tree.dat$novel[tree.dat$novel==0] <- "previously published"
tree.dat$novel[tree.dat$novel==1] <- "novel"
tree.dat$novel <- as.factor(tree.dat$novel)
shapez = c("novel" =  24, "previously published" = 21)
colz2 = c('novel' =  "yellow", 'previously published' = "white")

#shape by bat host
tree.dat$Bat_host[tree.dat$Bat_host==0] <- "non-bat host"
tree.dat$Bat_host[tree.dat$Bat_host==1] <- "bat host"
tree.dat$Bat_host <- as.factor(tree.dat$Bat_host)
shapez = c("bat host" =  24, "non-bat host" = 21)

tree.dat$Clade <- tree.dat$Strain
#need to check these against literature.
tree.dat$Clade[tree.dat$Host == "Sus scrofa" ] <- "Porcine AstV" 
tree.dat$Clade[tree.dat$Host == "Bos taurus" ] <- "Bovine AstV"
tree.dat$Clade[tree.dat$Host == "Homo sapiens" ] <- "Human AstV" 
tree.dat$Clade[tree.dat$Host == "Mus musculus" ] <- "Murine AstV" 
tree.dat$Clade[tree.dat$Host == "Felis catus" ] <- "Feline AstV" 
tree.dat$Clade[tree.dat$Host == "Marmota himalayana" ] <- "Marmot AstV" 
tree.dat$Clade[tree.dat$Host == "Ovis aries" ] <- "Bovine AstV" 
tree.dat$Clade[tree.dat$Host == "Mouse" ] <- "Murine AstV" 
tree.dat$Clade[tree.dat$Host == "Oryctolagus cuniculus" ] <- "Leporine AstV" 
tree.dat$Clade[tree.dat$Host == "Budorcas taxicolor tibetana" ] <- "Bovine AstV"
tree.dat$Clade[tree.dat$Bat_host=="bat host"] <- "Bat AstV"
tree.dat$Clade[tree.dat$Host == "Canis lupus familiaris" ] <- "Canine AstV"
tree.dat$Clade[tree.dat$Host == "Camelus dromedarius" ] <- "Camel AstV"
tree.dat$Clade[tree.dat$Host == "Neogale vison" ] <- "Mink AstV"
tree.dat$Clade[tree.dat$Host == "Rattus norvegicus" ] <- "Murine AstV"
tree.dat$Clade[tree.dat$Host == "Meleagris gallopavo" ] <- "AvastV"

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
  theme(legend.position = c(.9,.6), legend.title = element_blank()) +
  xlim(c(0,6))
p1_family

## color tip by order, shape by chiroptera, color by novelty
## CURRENT
p1_order <- ggtree(rooted.tree.A, size=.8) %<+% tree.dat +
  geom_tippoint(aes(color=Clade, fill=Clade), size=3.8) +
  scale_color_manual(values=c(AvastV = "azure4", 'Bat AstV' = "#0CB702", 'Bovine AstV' = "#00BFC4", 'Camel AstV' = "#00A9FF", 'Canine AstV' = "#8494FF", 'Feline AstV' = "#C77CFF", 'Human AstV' = "#ED68ED", 'Leporine AstV' = "#FF68A1", 'Marmot AstV' = "#E68613", 'Mink AstV' = "#CD9600", 'Murine AstV' = "#ABA300", 'Porcine AstV' = "#7CAE00")) +
  scale_fill_manual(values=c(AvastV = "azure4", 'Bat AstV' = "#0CB702", 'Bovine AstV' = "#00BFC4", 'Camel AstV' = "#00A9FF", 'Canine AstV' = "#8494FF", 'Feline AstV' = "#C77CFF", 'Human AstV' = "#ED68ED", 'Leporine AstV' = "#FF68A1", 'Marmot AstV' = "#E68613", 'Mink AstV' = "#CD9600", 'Murine AstV' = "#ABA300", 'Porcine AstV' = "#7CAE00")) +
  #scale_fill_npg() +
  #scale_color_npg() +
  #geom_nodelab(size=3,nudge_x = -.064, nudge_y = .6) +
  geom_treescale(fontsize=3) + 
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill=novel), geom = "label", label.size = 0, alpha=.3, size=5, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.8,.6), legend.title = element_blank()) +
  xlim(c(0,6))

#quartz()
p1_order

p12.dat <- p1_order$data
p12.dat$Bootstrap <- NA
Bootstrap<-p12.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p12.dat$label)] <- as.numeric(p12.dat$label[(length(tree.dat$tip_label)+1):length(p12.dat$label)])#fill with label

ppp <- p1_order %<+% p12.dat +
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = F), shape=21, stroke=0, size=2.5)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))

ppp






p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=2, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,14))

p1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree)
p0.dat <- p0$data
p0.dat$Bootstrap <- NA
Bootstrap<-p0.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p0.dat$label)] <- as.numeric(p0.dat$label[(length(tree.dat$tip_label)+1):length(p0.dat$label)])#fill with label

#add bootstrap values to original plot
p1.1 <- p1  %<+% p0.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"))
p1.1




##nodes to collapose
node1 <- MRCA(p1_order, 'NC_037655  |  Sichuan takin astrovirus  |  Budorcas taxicolor tibetana  |  China  |  2013', 'NC_023674  |  Porcine astrovirus 2  |  Sus scrofa  |  USA  |  2010')
node2 <- MRCA(p1_order, 'NC_019026  |  Astrovirus VA3  |  Homo sapiens  |  India  |  2005', 'NC_024472  |  Burkina Faso astrovirus  |  Homo sapiens  |  Burkina Faso  |  2010')

p_order2 <- p1_order %>% collapse(node=node1) +
  geom_point2(aes(subset=(node==node1)), shape=23, size=5, fill='#F8766D')
#p_order3 <- collapse(p_order2, node=node2) +
  #geom_point2(aes(subset=(node==node2)), shape=23, size=5, fill="#619CFF")
p_order2

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


ggsave(file = paste0(homewd, "final-figures/Fig2.png"),
       plot = p1_order,
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=4)#,






#####################
## SWIO RdRp
#####################

##Madabat tree
#load the fig3a tree
treeC <-  read.tree(file = paste0(homewd, "ML-trees/RdRp-SWIO-ML-tree-FPB"))

#remove quotes
treeC$tip.label <- gsub("'", '', treeC$tip.label)

#root it
#double rooting bc of bat avastrovirus
#1. get the node number so we can root by node
rooted.tree.C <- root(treeC, which(treeC$tip.label == "NC_002470"))

#take a quick look in base R
plot(rooted.tree.C)

#load tree data prepared from elsewhere
datC <- read.csv(file=paste0(homewd,"ML-trees/astro_meta_rdrp_SWIO.csv"), header = T, stringsAsFactors = F)
head(datC)

#manually change names to prevent errors
#datB[1,1] <- "F_MIZ141_RR034B_198_NODE_4_length_6593_cov_808.543744"
#datB[2,1] <- "F_MIZ141_RR034B_198_NODE_5_length_6456_cov_49.595532"
#check subgroup names
unique(datC$Genus)

colz = c("Mamastrovirus" = "royalblue", "Avastrovirus" = "tomato")

#pick order for the labels
datC$Genus <- factor(datC$Genus, levels = c("Mamastrovirus", "Avastrovirus"))   

#and add a "novel" category
#datB$novel = 0
#datB$novel[datB$Geo_Location=="Madagascar"] <- 1
#datB$novel <- as.factor(datB$novel)

#rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)

#take a glance
p3 <- ggtree(rooted.tree.C) %<+% datC + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = c(.2,.85), legend.title = element_blank())
p3 #why are some of the tips grey?

#now get new tip labels
datC$old_tip_label <- datC$Accession
datC$new_label <- NA
datC$new_label[!is.na(datC$Strain)] <- paste(datC$Accession[!is.na(datC$Strain)])


#datC$new_label[!is.na(datC$Strain)] <- paste(datC$Accession[!is.na(datC$Strain)], " | ", 
#                                             datC$Family[!is.na(datC$Strain)], "|" ,
#                                             datC$Host[!is.na(datC$Strain)])

#dat$new_label[is.na(dat$strain)] <- paste(dat$accession_num[is.na(dat$strain)], " | ", 
#                                         dat$host[is.na(dat$strain)], " | ",
#                                          dat$country[is.na(dat$strain)], " | ",
#                                          dat$collection_year[is.na(dat$strain)])

#after Gwen checks these, can later manually edit any that are "NA"

#make sure to sort in order
tree.datC <- data.frame(old_tip_label=rooted.tree.C$tip.label, num =1:length(rooted.tree.C$tip.label))
head(tree.datC)
head(datC)

tree.datC <- merge(tree.datC, datC, by = "old_tip_label", all.x = T, sort = F)

names(tree.datC)

tree.datC$tip_label <- tree.datC$new_label
tree.datC <- dplyr::select(tree.datC, tip_label, Accession, Strain, Host, Geo_Location, Year, Genus, novel, old_tip_label, Family, Animal, suborder)

rooted.tree.C$tip.label <- tree.datC$tip_label

tree.datC$novel[tree.datC$novel==0] <- "previously published"
tree.datC$novel[tree.datC$novel==1] <- "novel"
tree.datC$novel <- as.factor(tree.datC$novel)
shapez = c("novel" =  24, "previously published" = 21)
colz2 = c("novel" =  "yellow", "previously published" = "white")

tree.datC$Geo_Location <- as.factor(tree.datC$Geo_Location)
shape_geo = c("Madagascar" = 17, "Mozambique" = 16, "Reunion Island" = 15, "Unknown" = 4)

tree.datC$suborder <- as.factor(tree.datC$suborder)
shape_suborder = c("Yangochiroptera" = 24, "Yinpterochiroptera" = 21)

## color tip by family
pC_family <- ggtree(rooted.tree.C) %<+% tree.datC + geom_tippoint(aes(color=Family, fill=Family, shape=novel)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  #scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.6), legend.title = element_blank()) +
  xlim(c(0,5.6))
pC_family

## color tip by geo_location
pC_geo <- ggtree(rooted.tree.C) %<+% tree.datC + geom_tippoint(aes(color=Geo_Location, fill=Geo_Location, shape=novel)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  #scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.6), legend.title = element_blank()) +
  xlim(c(0,5.6))
pC_geo

## color by sub-order, shape by location
##CURRENT
pC_suborder <- ggtree(rooted.tree.C) %<+% tree.datC + geom_tippoint(aes(color=suborder, fill=suborder, shape=Geo_Location)) +
  #geom_hilight(node = 220, fill = "darkgoldenrod1", alpha = 0.6) +
  geom_nodelab(aes(subset = as.numeric(label) > 30), size=1.5, nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shape_geo) + 
  new_scale_fill() +
  geom_tiplab(aes(fill=novel, color=suborder), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  scale_color_manual(values=c(Yangochiroptera = "darkslategray4", Yinpterochiroptera = "orangered3")) +
  #scale_color_manual(values=c(Hipposideridae = "#F8766D", Rhinonycteridae = "#E68613", Pteropodidae = "#CD9600", Vespertilionidae = "#00BFC4", Nycteridae = "#00B8E7", Molossidae = "#00A9FF", Miniopteridae = "#00C19A", Phasianidae = "azure4")) +
  theme(legend.position = c(.1,.6), legend.title = element_blank()) +
  geom_treescale(x= 2.15, y= .05, fontsize=1.5) + 
  xlim(c(0,5.6))
pC_suborder


## rotate clades
MRCA(pC_suborder, 'KY575655  |  Vespertilionidae | Myotis goudoti', 'MH013985  |  Rhinonycteridae | Triaenops afer')
pR <- ggtree::rotate(pC_suborder, 176)
pR2 <- ggtree::rotate(pR, 206)
pR3 <- ggtree::rotate(pR2, 207)
pR4 <- ggtree::rotate(pR3, 208)
pR4


#collapse nodes
node1 = MRCA(pC_suborder, 'KY575657', 'KY575658')
p1 <- pC_suborder %>% collapse(node=node1)
p1




## color by location, shape by suborder
pC_suborder2 <- ggtree(rooted.tree.C) %<+% tree.datC + geom_tippoint(aes(color=Geo_Location, fill=Geo_Location, shape=suborder)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shape_suborder) + 
  new_scale_fill() +
  geom_tiplab(geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  #scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.6), legend.title = element_blank()) +
  xlim(c(0,5.6))
pC_suborder2

ggsave(file = paste0(homewd, "/final-figures/ML-Mada-bat-RdRp_sept1.png"),
       plot = pB_family,
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=4)#, 



#####################
## Africa RdRp
#####################

##Madabat tree
#load the fig3a tree
treeD <-  read.tree(file = paste0(homewd, "Phylogenies/RdRp-Africa-ML-tree"))

#remove quotes
treeD$tip.label <- gsub("'", '', treeD$tip.label)

#root it
#double rooting bc of bat avastrovirus
#1. get the node number so we can root by node
rooted.tree.D <- root(treeD, which(treeD$tip.label == "NC_002470"))

#take a quick look in base R
plot(rooted.tree.D)

#load tree data prepared from elsewhere
datD <- read.csv(file=paste0(homewd,"Phylogenies/astro_meta_rdrp_Africa.csv"), header = T, stringsAsFactors = F)
head(datD)

#manually change names to prevent errors
#datB[1,1] <- "F_MIZ141_RR034B_198_NODE_4_length_6593_cov_808.543744"
#datB[2,1] <- "F_MIZ141_RR034B_198_NODE_5_length_6456_cov_49.595532"
#check subgroup names
unique(datD$Genus)

colz = c("Mamastrovirus" = "royalblue", "Avastrovirus" = "tomato")

#pick order for the labels
datD$Genus <- factor(datD$Genus, levels = c("Mamastrovirus", "Avastrovirus"))   

#and add a "novel" category
#datB$novel = 0
#datB$novel[datB$Geo_Location=="Madagascar"] <- 1
#datB$novel <- as.factor(datB$novel)

#rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)

#take a glance
p4 <- ggtree(rooted.tree.D) %<+% datD + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = c(.2,.85), legend.title = element_blank())
p4 #why are some of the tips grey?

#now get new tip labels
datD$old_tip_label <- datD$Accession
datD$new_label <- NA
datD$new_label[!is.na(datD$Strain)] <- paste(datD$Accession[!is.na(datD$Strain)], " | ", 
                                             datD$Family[!is.na(datD$Strain)], "|" ,
                                             datD$Host[!is.na(datD$Strain)])

#dat$new_label[is.na(dat$strain)] <- paste(dat$accession_num[is.na(dat$strain)], " | ", 
#                                         dat$host[is.na(dat$strain)], " | ",
#                                          dat$country[is.na(dat$strain)], " | ",
#                                          dat$collection_year[is.na(dat$strain)])

#after Gwen checks these, can later manually edit any that are "NA"

#make sure to sort in order
tree.datD <- data.frame(old_tip_label=rooted.tree.D$tip.label, num =1:length(rooted.tree.D$tip.label))
head(tree.datD)
head(datD)

tree.datD <- merge(tree.datD, datD, by = "old_tip_label", all.x = T, sort = F)

names(tree.datD)

tree.datD$tip_label <- tree.datD$new_label
tree.datD <- dplyr::select(tree.datD, tip_label, Accession, Strain, Host, Geo_Location, Year, Genus, novel, old_tip_label, Family, Animal)

rooted.tree.D$tip.label <- tree.datD$tip_label

tree.datD$novel[tree.datD$novel==0] <- "previously published"
tree.datD$novel[tree.datD$novel==1] <- "novel"
tree.datD$novel <- as.factor(tree.datD$novel)
shapez = c("novel" =  24, "previously published" = 21)
colz2 = c('1' =  "yellow", '0' = "white")

#tree.datD$Geo_Location <- as.factor(tree.datD$Geo_Location)
#shape_geo = c("Madagascar" = 24, "Mozambique" = 21, "Reunion Island" = 15, "Unknown" = 4, "Egypt")

tree.datD$suborder <- as.factor(tree.datD$suborder)
shape_suborder = c("Yangochiroptera" = 24, "Yinpterochiroptera" = 21)

## color tip by family
pD_family <- ggtree(rooted.tree.C) %<+% tree.datD + geom_tippoint(aes(color=Family, fill=Family, shape=novel)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  #scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.6), legend.title = element_blank()) +
  xlim(c(0,5.6))
pD_family

## color tip by geo_location
pD_geo <- ggtree(rooted.tree.D) %<+% tree.datD + geom_tippoint(aes(color=Geo_Location, fill=Geo_Location, shape=novel)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  #scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.6), legend.title = element_blank()) +
  xlim(c(0,5.6))
pD_geo

## color by sub-order, shape by location
pD_suborder <- ggtree(rooted.tree.C) %<+% tree.datD + geom_tippoint(aes(color=suborder, fill=suborder, shape=Geo_Location)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shape_geo) + 
  new_scale_fill() +
  geom_tiplab(geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  #scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.6), legend.title = element_blank()) +
  xlim(c(0,5.6))
pD_suborder

## color by location, shape by suborder
pD_suborder2 <- ggtree(rooted.tree.C) %<+% tree.datD + geom_tippoint(aes(color=Geo_Location, fill=Geo_Location, shape=suborder)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shape_suborder) + 
  new_scale_fill() +
  geom_tiplab(geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  #scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.1,.6), legend.title = element_blank()) +
  xlim(c(0,5.6))
pD_suborder2

ggsave(file = paste0(homewd, "/final-figures/ML-Mada-bat-RdRp_sept1.png"),
       plot = pB_family,
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=4)#, 










#####################
## Combining figures
#####################

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










##graveyard

#########################
##### ML RdRp Mada Bat #####
#########################
#**this tree is being absorbed into the SWIO tree

##Madabat tree
#load the fig3a tree
treeB <-  read.tree(file = paste0(homewd, "Phylogenies/RdRp-Mada-ML-tree"))

#remove quotes
treeB$tip.label <- gsub("'", '', treeB$tip.label)

#root it
rooted.tree.B <- root(treeB, which(treeB$tip.label == "NC_002470___RdRp"))

#take a quick look in base R
plot(rooted.tree.B)

#load tree data prepared from elsewhere
datB <- read.csv(file=paste0(homewd,"Phylogenies/astro_meta_rdrp_madabat.csv"), header = T, stringsAsFactors = F)
head(datB)

#manually change names to prevent errors
#datB[1,1] <- "F_MIZ141_RR034B_198_NODE_4_length_6593_cov_808.543744"
#datB[2,1] <- "F_MIZ141_RR034B_198_NODE_5_length_6456_cov_49.595532"
#check subgroup names
unique(datB$Genus)

colz = c("Mamastrovirus" = "royalblue", "Avastrovirus" = "tomato")

#pick order for the labels
datB$Genus <- factor(datB$Genus, levels = c("Mamastrovirus", "Avastrovirus"))   

#and add a "novel" category
#datB$novel = 0
#datB$novel[datB$Geo_Location=="Madagascar"] <- 1
#datB$novel <- as.factor(datB$novel)

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
                                             datB$Host[!is.na(datB$Strain)])

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
tree.datB <- dplyr::select(tree.datB, tip_label, Accession, Strain, Host, Geo_Location, Year, Genus, novel, old_tip_label, Family, Animal)

rooted.tree.B$tip.label <- tree.datB$tip_label

tree.datB$novel[tree.datB$novel==0] <- "previously published"
tree.datB$novel[tree.datB$novel==1] <- "novel"
tree.datB$novel <- as.factor(tree.datB$novel)
shapez = c("novel" =  24, "previously published" = 21)
colz2 = c('1' =  "yellow", '0' = "white")

## Tree 1: Colored by family, shaped by novel or not novel
# doesn't work
pb <- ggtree(rooted.tree.B) %<+% tree.datB + geom_tippoint(aes(color=Family, fill=Family, shape=novel)) +
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
pB_family <- ggtree(rooted.tree.B) %<+% tree.datB + geom_tippoint(aes(color=Family, fill=Family, shape=novel)) +
  geom_nodelab(size=1.5,nudge_x = -.05, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  #scale_fill_manual(values=colz2) + 
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


ggsave(file = paste0(homewd, "/final-figures/ML-Mada-bat-RdRp_sept1.png"),
       plot = pB_family,
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=4)#, 


