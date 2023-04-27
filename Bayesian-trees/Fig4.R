rm(list=ls())


library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(lubridate)
#library(rBt)
library(treeio)


#make Bayesian timetree from Astrovirus strict molecular clock model

#first, read in the tree

homewd= "/Users/sophiahorigan/Documents/GitHub/Mada-Bat-Astro/"

setwd(paste0(homewd, "/Bayesian-trees"))


tree <-  read.beast(file = paste0(homewd, "Bayesian-trees/astro-full-beast-MCC"))
#tree <-  read.beast(file = paste0(homewd, "Bayesian-trees/TESTT"))


#test <- root(tree, which(tree@phylo$tip.label == "NC_002470_1985-07-15"))


treedat <- cbind.data.frame(tip_name = tree@phylo$tip.label)
treedat$beast_name <-treedat$tip_name
#tree <- read.annot.beast(file = paste0(homewd, "/Fig4/beast-out/AllNobeco/NobecoStrict/AvgNobecoStrictNexus.trees"))

#tree$node.label <- round(tree$posterior,2)
#treedat <- cbind.data.frame(tip_name = tree$tip.label)
treedat$accession_num <- sapply(strsplit(treedat$tip_name, "_"), function(x) x[[1]])
treedat$accession_num[treedat$accession_num=="NC"] <- c("NC_001943", "NC_002469", "NC_004579", "NC_011400", "NC_013060", "NC_013443", "NC_015935", "NC_016155", "NC_016896", "NC_018702", "NC_019026", "NC_019027", "NC_019028", "NC_019494", "NC_022249", "NC_023629", "NC_023630", "NC_023631", "NC_023632", "NC_023636", "NC_023674", "NC_023675", "NC_024297", "NC_024472", "NC_024498", "NC_024701", "NC_025346", "NC_025379", "NC_026814", "NC_027711", "NC_030922", "NC_033792", "NC_033821", "NC_034974", "NC_036583", "NC_037655")
treedat$accession_num[treedat$accession_num=="F"] <- c("OQ606244")

#names(treedat)[names(treedat)=="tip_name"] <- "beast_name"


#and load data of corresponding tree

dat <- read.csv(file = "astro_full_beast_meta.csv", header = T, stringsAsFactors = F)
dat$collection_date <- as.Date(dat$collection_date)


#test 

mrsd.dat <- max(dat$collection_date)
p1 <- ggtree(tree, mrsd=mrsd.dat)  + theme_tree2()  +geom_nodelab()

tree.dat <- p1$data
node.sub <- dplyr::select(tree.dat, node, x)
names(node.sub) <-  c("node", "nodetime")

#and 
head(dat)

dat$clade <- dat$strain
#need to check these against literature.
dat$clade[dat$host == "Sus scrofa" ] <- "Porcine AstV" 
dat$clade[dat$host == "Bos taurus" ] <- "Bovine AstV"
dat$clade[dat$host == "Homo sapiens" ] <- "Human AstV" 
dat$clade[dat$host == "Mus musculus" ] <- "Murine AstV" 
dat$clade[dat$host == "Felis catus" ] <- "Feline AstV" 
dat$clade[dat$host == "Marmota himalayana" ] <- "Marmot AstV" 
dat$clade[dat$host == "Ovis aries" ] <- "Bovine AstV" 
dat$clade[dat$host == "Mouse" ] <- "Murine AstV" 
dat$clade[dat$host == "Oryctolagus cuniculus" ] <- "Leporine AstV" 
dat$clade[dat$host == "Budorcas taxicolor tibetana" ] <- "Bovine AstV"
dat$clade[dat$bat_host==1] <- "Bat AstV"
dat$clade[dat$host == "Canis lupus familiaris" ] <- "Canine AstV"
dat$clade[dat$host == "Camelus dromedarius" ] <- "Camel AstV"
dat$clade[dat$host == "Neogale vison" ] <- "Mink AstV"
dat$clade[dat$host == "Rattus norvegicus" ] <- "Murine AstV"
dat$clade[dat$host == "Meleagris gallopavo" ] <- "AvastV"




dat.plot <- merge(treedat, dat, by="accession_num", all.x = T, sort=F)

head(dat.plot)
dat.plot$new_label = NA
dat.plot$new_label[!is.na(dat.plot$strain)] <- paste(dat.plot$accession_num[!is.na(dat.plot$strain)], " | ", 
                                                     dat.plot$strain[!is.na(dat.plot$strain)], " | ", 
                                                     dat.plot$host[!is.na(dat.plot$strain)], " | ",
                                                     dat.plot$country[!is.na(dat.plot$strain)], " | ",
                                                     dat.plot$collection_year[!is.na(dat.plot$strain)])

dat.plot$new_label[is.na(dat.plot$strain)] <- paste(dat.plot$accession_num[is.na(dat.plot$strain)], " | ", 
                                                    dat.plot$host[is.na(dat.plot$strain)], " | ",
                                                    dat.plot$country[is.na(dat.plot$strain)], " | ",
                                                    dat.plot$collection_year[is.na(dat.plot$strain)])


tree@phylo$tip.label <- dat.plot$new_label

dat.sub <- dplyr::select(dat.plot, new_label, collection_date, country, clade)
head(dat.sub)
dat.sub$clade <- as.factor(dat.sub$clade)

dat.sub$novel = "no"
dat.sub$novel[dat.sub$country=="Madagascar"] <- "yes"

colz2 = c('yes' =  "yellow", 'no' = "white")

p3 <-ggtree(tree, mrsd=mrsd.dat, size=.8) %<+% dat.sub +
  geom_tippoint(aes(color=clade), size=4) +
  scale_color_manual(values=c('Bat AstV' = "#0CB702", 'Bovine AstV' = "#00BFC4", 'Camel AstV' = "#00A9FF", 'Canine AstV' = "#8494FF", 'Feline AstV' = "#C77CFF", 'Human AstV' = "#ED68ED", 'Leporine AstV' = "#FF68A1", 'Marmot AstV' = "#E68613", 'Mink AstV' = "#CD9600", 'Murine AstV' = "#ABA300", 'Porcine AstV' = "#7CAE00")) +
  scale_fill_manual(values=c('Bat AstV' = "#0CB702", 'Bovine AstV' = "#00BFC4", 'Camel AstV' = "#00A9FF", 'Canine AstV' = "#8494FF", 'Feline AstV' = "#C77CFF", 'Human AstV' = "#ED68ED", 'Leporine AstV' = "#FF68A1", 'Marmot AstV' = "#E68613", 'Mink AstV' = "#CD9600", 'Murine AstV' = "#ABA300", 'Porcine AstV' = "#7CAE00")) +
  #geom_nodelab(size=2.5,nudge_x = -21, nudge_y = .7) +
 # geom_nodelab(aes(label=new.nodel.lab), size=4,nudge_x = -55, nudge_y = -1,  color="firebrick", fontface=2, geom="label", fill="white") +
  #geom_nodelab(size=1.8,nudge_x = -.05, nudge_y = .7) +
  #geom_treescale(fontsize=2.5) + 
  theme_tree2() +
  theme(legend.position = c(.01,.7), plot.margin = unit(c(.2,24,3,3), "lines"), legend.title = element_blank()) +
  coord_cartesian(clip = "off") + 
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=3, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") +
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1)) +
  scale_x_continuous(breaks=c(1620, 1720, 1820, 1920, 2020),
                     labels=c(400, 300, 200, 100, 0)) +
  xlab("years to MRCA") +
ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3,  show.legend=F, size=3) + scale_fill_manual(values=colz2) 
  #xlim(c(0,6))

p3


ggsave(file = paste0(homewd, "/final-figures/Fig4.png"),
       units="mm",  
       width=90, 
       height=60, 
       #limitsize = F,
       scale=3)#, 




#calculate how long ago Rousettus split from Eidolon
orig.date <- round(node.sub$nodetime[33],0)

nodeRous <- MRCA(tree, which(tree@phylo$tip.label =="OQ606244  |  Bat astrovirus  |  Rousettus madagascariensis  |  Madagascar  |  2018"), which(tree@phylo$tip.label =="MG693176  |  Bat astrovirus  |  Eidolon helvum  |  Cameroon  |  2013"))
Rous.date <- round(node.sub$nodetime[nodeRous],0)

Rous.date <- Rous.date - 2018



