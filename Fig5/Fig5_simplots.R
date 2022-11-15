rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

homewd="/Users/sophiahorigan/Documents/Github/Mada-Bat-Astro/"
setwd(paste0(homewd, "/Fig5/"))


#########################################
##  All bats ##
###################################


simplot2 <- read.csv(file = "all_bat.csv", header = T, stringsAsFactors = F)
head(simplot2)

#move to long
long.sim2 <- melt(simplot2, id.vars = c("pointer"), measure.vars = c("MG693176", "MN832787", "MT734809", "MZ218053", "MZ218054"))

head(long.sim2)
unique(long.sim2$variable)
long.sim2$variable <- as.character(long.sim2$variable)
#long.sim2$variable[long.sim2$variable=="Alternate"] <- long.sim2$Alternate_ID[long.sim2$variable=="Alternate"] 
#unique(long.sim2$variable)
names(long.sim2)[names(long.sim2)=="variable"] <- "strain"

long.sim2$strain[long.sim2$strain=="MG693176"] <- "Eidolon helvum"
long.sim2$strain[long.sim2$strain=="MN832787"] <- "Myotis daubentonii 1"
long.sim2$strain[long.sim2$strain=="MT734809"] <- "Chiroptera spp."
long.sim2$strain[long.sim2$strain=="MZ218053"] <- "Myotis daubentonii 2"
long.sim2$strain[long.sim2$strain=="MZ218054"] <- "Myotis daubentonii 3"

long.sim2$strain <- factor(long.sim2$strain, levels = c("Eidolon helvum", "Myotis daubentonii 1", "Myotis daubentonii 2", "Myotis daubentonii 3", "Chiroptera spp."))

#and plot

long.sim2$value[long.sim2$value<0] <- 0
long.sim2$value <- long.sim2$value/100
#long.sim2$Query[long.sim2$Query=="Pteropus_rufus"] <- "Pteropus rufus"
#long.sim2$Query[long.sim2$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"

genome.df <- data.frame(position = c(1, 2600,
                                     2601, 4150,
                                     4151, 6700), 
                        gene = rep(c("ORF1a", "ORF1b", "ORF2"), each=2))


genome.df$gene <- factor(genome.df$gene, levels = unique(genome.df$gene))


colz2= c("Eidolon helvum"="firebrick3", "Myotis daubentonii 1" = "royalblue", "Myotis daubentonii 2" = "lightblue","Myotis daubentonii 3" = "blue","Chiroptera spp." = "goldenrod")



allbat <- ggplot(long.sim2) + geom_line(aes(x=pointer, y=value, color=strain), size=1) +
  #geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05,  fill=gene), color="black") +
  theme_bw() + xlab("genome position") + ylab("nucleotide similarity") + ylim(0,1) +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  scale_color_manual(values=colz2) + #coord_cartesian(ylim=c(-.1,1)) +
  scale_x_continuous(breaks=c(0,2000,4000,6000), labels = c(0,2000, 4000,6000))

allbat





#####################
## EXTRA 
#####################








#add in the bootscan
bootplot <- read.csv(file = "all_bootscan.csv", header = T, stringsAsFactors = F)
head(bootplot)

#move to long
long.boot <- melt(bootplot, id.vars = c("CenterPos", "ID_Alternative", "Query"), measure.vars = c("HKU9", "Eidolon_helvum", "Alternative"))

head(long.boot)
unique(long.boot$variable)
long.boot$variable <- as.character(long.boot$variable)
long.boot$variable[long.boot$variable=="Alternative"] <- long.boot$ID_Alternative[long.boot$variable=="Alternative"] 
unique(long.boot$variable)
names(long.boot)[names(long.boot)=="variable"] <- "strain"

long.boot$strain_label <- NA

#ylab(bquote("r"^"*"~", virus growth")) 

long.boot$strain[long.boot$strain=="Eidolon_helvum"] <- "E. helvum bat coronavirus"
long.boot$strain[long.boot$strain=="Pteropus_rufus"] <- "P. rufus Nobecovirus"
long.boot$strain[long.boot$strain=="Rousettus_madagascariensis"] <-  "R. madagascariensis Nobecovirus"

long.boot$strain <- factor(long.boot$strain, levels = c("HKU9", "E. helvum bat coronavirus", "P. rufus Nobecovirus", "R. madagascariensis Nobecovirus"))

long.boot$value[long.boot$value<0] <- 0
long.boot$Query[long.boot$Query=="Pteropus_rufus"] <- "Pteropus rufus"
long.boot$Query[long.boot$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"



p3 <- ggplot(long.boot) + geom_line(aes(x=CenterPos, y=value, color=strain), show.legend = F, size=.9) +
  geom_ribbon(data=genome.df2, aes(x=position, ymin=-10, ymax=-5, fill=gene), color="black", show.legend = F) +
  facet_grid(~Query) + theme_bw() + xlab("genome position") + ylab("% of Permuted Trees") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin = unit(c(.1, .1,.1, .3), "cm"),
        legend.text = element_text(face="italic"),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  scale_color_manual(values=colz)

p3







#and together
Fig5top <- cowplot::plot_grid(p1,p2,p3, nrow=3, ncol = 1, labels= c("(A)", "(B)", "(C)"), label_size = 16, label_x = -.01, rel_heights = c(1,.9,1.2))

#Fig5 <- cowplot::plot_grid(Fig5top, leg1, nrow = 1, ncol = 2, rel_widths = c(1,.2))


Fig5 <- cowplot::plot_grid(Fig5top, leg1, nrow = 2, ncol = 1, rel_heights = c(1,.1))

ggsave(file = paste0(homewd, "/final-figures/Fig5.png"),
       plot=Fig5,
       units="mm",  
       width=95, 
       height=70, 
       #limitsize = F,
       scale=4)#, 



#######################################
##  MIZ141 vs human ##
#####################################


#id.plot <- read.csv(file = paste0(homewd, "/Fig5/AA_identity_P_ruf_R_mad.csv"), header = T, stringsAsFactors = F)

dat1 <- read.csv(file = paste0(homewd, "Fig5/all_bat.csv"), header = T, stringsAsFactors = F)

#move to long
id.plot <- melt(dat1, id.vars = c("pointer"), measure.vars = c("F_MIZ141_1"))
id.plot$variable <- as.character(id.plot$variable)
#id.plot$variable[id.plot$variable=="Alternate"] <- id.plot$Alternate_ID[id.plot$variable=="Alternate"]
head(id.plot)


names(id.plot)[names(id.plot)=="variable"] <- "strain"

#id.plot$strain[id.plot$strain=="Eidolon_helvum"] <- "E. helvum bat coronavirus"
#id.plot$strain[id.plot$strain=="Pteropus_rufus"] <- "P. rufus Nobecovirus"
#id.plot$strain[id.plot$strain=="Rousettus_madagascariensis"] <- "R. madagascariensis Nobecovirus"

#id.plot$strain <- factor(id.plot$strain, levels = c("HKU9", "GCCDC1", "GX2018.BtCoV92", "E. helvum bat coronavirus", "P. rufus Nobecovirus", "R. madagascariensis Nobecovirus"))
#id.plot$strain <- factor(id.plot$strain, levels = c("HKU9", "E. helvum bat coronavirus", "P. rufus Nobecovirus", "R. madagascariensis Nobecovirus"))


#and plot

id.plot$value[id.plot$value<0] <- 0
# id.plot$Query[id.plot$Query=="Pteropus_rufus"] <- "Pteropus rufus"
# id.plot$Query[id.plot$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"
# id.plot$Query <- factor(id.plot$Query, labels = c("Pteropus rufus", "Rousettus madagascariensis"))
S

genome.df <- data.frame(position = c(1, 2600,
                                     2601, 4150,
                                     4151, 6700), 
                        gene = rep(c("ORF1a", "ORF1b", "ORF2"), each=2))

id.plot$value <- id.plot$value/100
genome.df$gene <- factor(genome.df$gene, levels = unique(genome.df$gene))

colz= c("F_MIZ141_1"="goldenrod", "F_MIZ141_2" = "forestgreen")
#colz= c("HKU9"="firebrick3", "E. helvum bat coronavirus" = "royalblue", "P. rufus Nobecovirus" = "goldenrod", "R. madagascariensis Nobecovirus" = "forestgreen")

bat_human <- ggplot(id.plot) + geom_line(aes(x=pointer, y=value, color=strain), size=1) +
  geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05,  fill=gene), color="black") +
  theme_bw() + xlab("genome position") + ylab("nucleotide similarity") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  scale_color_manual(values=colz) + #coord_cartesian(ylim=c(-.1,1)) +
  scale_x_continuous(breaks=c(0,2000,4000,6000), labels = c(0,2000, 4000,6000))

bat_human

leg1 <- cowplot::get_legend(bat_human)


## amino acid

# p1 <- ggplot(id.plot) + geom_line(aes(x=pointer, y=value, color=strain), size=1, show.legend = F) +
#   geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05, fill=gene), color="black", show.legend = F) +
#   facet_grid(~Query) + theme_bw() + xlab("genome position") + ylab("amino acid similarity") +
#   theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
#         strip.background = element_rect(fill="white"),
#         plot.margin = unit(c(.1, .1,0, .2), "cm"),
#         legend.text = element_text(face="italic"),
#         axis.text.y = element_text(size=12), axis.title.y = element_text(size=14),
#         axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
#   scale_color_manual(values=colz)
# p1
