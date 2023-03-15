rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gggenes)
library(wesanderson)
#library(LaCroixColoR)
#install.packages("LaCroixColoR")
library(devtools)

## first, let's start out with amino acid similarity

homewd="/Users/shorigan/Documents/Github/Mada-Bat-Astro/"
setwd(paste0(homewd, "/Fig5/"))

# load in symplot csv output file 

dat1 <- read.csv(file = paste0(homewd, "orf.csv"), header = T, stringsAsFactors = F)

# move to long 
# this will change orientation of csv file for easier visualization

id.plot <- melt(dat1, id.vars = c("pointer"), measure.vars = c("NC_001918",  "KJ934637"))

# since our variable name is our accession number, change to character
# then head to make sure structure is fine 

id.plot$variable <- as.character(id.plot$variable)
head(id.plot)

# let's then change the name of our variable to strain
# this is essential for visualization at the end

names(id.plot)[names(id.plot)=="variable"] <- "strain"

#id.plot$strain[id.plot$strain=="Eidolon_helvum"] <- "E. helvum Kobuvirus"

#id.plot$strain <- factor(id.plot$strain, levels = c("E. helvum Kobuvirus",  "R. madagascariensis Kobuvirus"))
head(id.plot)


# ignore this

#id.plot$value[id.plot$value<0] <- 0
#id.plot$Query[id.plot$Query=="Eidolon_helvum"] <- "Eidolon helvum"
#id.plot$Query[id.plot$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"
#id.plot$Query <- factor(id.plot$Query, labels = c("Pteropus rufus", "Rousettus madagascariensis"))

# 
# genome.df <- data.frame(position = c(1, 4316,
#                                      4316, 7002,
#                                      7002, 8324,
#                                      8324, 8580,
#                                      8580, 8659,
#                                      8659, 8881,
#                                      8881, 9378,
#                                      9378, 9820), 
#                         gene = rep(c("ORF1a", "ORF1b", "S", "NS3", "E", "M", "N", "NS7"), each=2))

# define genome positions 
# you'll have to have previously characterized the genome to know these positions 

genome.df <- data.frame(position = c(1, 185,
                                     186, 556,
                                     557, 779,
                                     780, 1021,
                                     1022, 1157,
                                     1158, 1322,
                                     1323, 1657, 
                                     1658, 1750, 
                                     1751, 1776, 
                                     1777, 1966, 
                                     1967, 2472), 
                        Peptide = rep(c("L", "VP0", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"), each=2))

# normalize similarity 
id.plot$value <- id.plot$value/100
genome.df$Peptide <- factor(genome.df$Peptide, levels = unique(genome.df$Peptide))

# give the sequences that are being plotted a color

colz= c("NC_001918" = "firebrick3", "KJ934637" = "forestgreen")

# plot 

 p3 <- ggplot(id.plot) + geom_line(aes(x=pointer, y=value, color=strain), size=1) + 
    geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05,  fill=Peptide), color="black") + 
    facet_grid() + theme_bw() + xlab("Genome position") + ylab("Amino acid similarity") +
    theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
          strip.background = element_rect(fill="white"), 
          legend.position = "bottom", legend.direction = "horizontal",# legend.box = "vertical",
          legend.text = element_text(face="italic", size = 12),
          axis.text = element_text(size=12), axis.title = element_text(size=14)) + 
   ggtitle("Amino acid similarity to OP287812") + 
   scale_fill_manual(values = wes_palette("Darjeeling2", 11, type = "continuous"))
 
 
 

## next, let's do our full genome nucleotide similarity

homewd="/Users/flg9/Desktop/Developer/brook_lab/symplot/nucleotide/"
setwd(paste0(homewd))

dat2 <- read.csv(file = paste0(homewd, "nt_alignment.csv"), header = T, stringsAsFactors = F)

# move to long
# this will change orientation of csv file for easier visualization

id.plot2 <- melt(dat2, id.vars = c("pointer"), measure.vars = c("KJ934637",  "NC_001918"))

# since our variable name is our accession number, change to character
# then head to make sure structure is fine 

id.plot2$variable <- as.character(id.plot2$variable)
head(id.plot2)


# let's then change the name of our variable to strain
# this is essential for visualization at the end

names(id.plot2)[names(id.plot2)=="variable"] <- "strain"

#id.plot2$strain[id.plot2$strain=="Eidolon_helvum"] <- "E. helvum Kobuvirus"
#id.plot2$strain[id.plot2$strain=="Rousettus_madagascariensis"] <- "R. madagascariensis Kobuvirus"

#id.plot$strain <- factor(id.plot$strain, levels = c("HKU9", "GCCDC1", "GX2018.BtCoV92", "E. helvum bat coronavirus", "P. rufus Nobecovirus", "R. madagascariensis Nobecovirus"))
#id.plot$strain <- factor(id.plot$strain, levels = c("E. helvum Kobuvirus",  "R. madagascariensis Kobuvirus"))
#head(id.plot)


#and plot

#id.plot$value[id.plot$value<0] <- 0
#id.plot$Query[id.plot$Query=="Eidolon_helvum"] <- "Eidolon helvum"
#id.plot$Query[id.plot$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"
#id.plot$Query <- factor(id.plot$Query, labels = c("Pteropus rufus", "Rousettus madagascariensis"))

# 
# genome.df <- data.frame(position = c(1, 4316,
#                                      4316, 7002,
#                                      7002, 8324,
#                                      8324, 8580,
#                                      8580, 8659,
#                                      8659, 8881,
#                                      8881, 9378,
#                                      9378, 9820), 
#                         gene = rep(c("ORF1a", "ORF1b", "S", "NS3", "E", "M", "N", "NS7"), each=2))

# define genome positions 
genome.df2 <- data.frame(position = c(1, 679,
                                     680, 1234,
                                     1235, 2347,
                                     2348, 3016,
                                     3017, 3742,
                                     3743, 4150,
                                     4151, 4645, 
                                     4646, 5650, 
                                     5651, 5929, 
                                     5930, 6007, 
                                     6008, 6577,
                                     6578, 7984,
                                     7985, 8505), 
                        Peptide = rep(c("5UTR", "L", "VP0", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D", "3UTR"), each=2))


# normalize similarity 
id.plot2$value <- id.plot2$value/100
genome.df2$Peptide <- factor(genome.df2$Peptide, levels = unique(genome.df2$Peptide))

#colz= c("HKU9"="firebrick3", "GCCDC1"="magenta", "GX2018.BtCoV92" = "purple", "E. helvum bat coronavirus" = "royalblue", "P. rufus Nobecovirus" = "goldenrod", "R. madagascariensis Nobecovirus" = "forestgreen")
#colz= c("HKU9"="firebrick3", "E. helvum bat coronavirus" = "royalblue", "P. rufus Nobecovirus" = "goldenrod", "R. madagascariensis Nobecovirus" = "forestgreen")

#colz= c("E. helvum Kobuvirus" = "firebrick3", "R. madagascariensis Kobuvirus" = "forestgreen")

#colz= c("NC_001918" = "firebrick3", "KJ934637" = "forestgreen")


# and plot 
p2 <- ggplot(id.plot2) + geom_line(aes(x=pointer, y=value, color=strain), size=1) + 
  scale_y_continuous(limits = c(-1, 1)) + 
  geom_ribbon(data=genome.df2, aes(x=position, ymin=-.1, ymax= -.05,  fill=Peptide), color="black") + 
  facet_grid() + theme_bw() + xlab("Genome position") + ylab("Nucleotide similarity") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal", #legend.box = "vertical",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) + 
  ggtitle("Nucleotide similarity to OP287812") + ylim(-.1,1) +
  scale_fill_manual(values = wes_palette("Darjeeling2", 13, type = "continuous"))


## now let's plot coverage 

homewd="/Users/flg9/Desktop/Developer/brook_lab/symplot/coverage/"
setwd(paste0(homewd))

datcovg <- read.csv(file = paste0(homewd, "Coverage.csv"), header = T, stringsAsFactors = F)
position <- read.csv(file = paste0(homewd, "position.csv"), header = T, stringsAsFactors = F)

#datcovg$Coverage <- as.character(datcovg$Coverage)
names(datcovg)[names(datcovg)=="Coverage"] <- "Coverage"

# this will give raw coverage

datcovg$Coverage <- datcovg$Coverage/100

# this will give rpm

datcovg$Coverage <- datcovg$Coverage/31.282
#genome.df2$Peptide <- factor(genome.df2$Peptide, levels = unique(genome.df2$Peptide))
colnames(genome.df2) <- c("Position", "Peptide")

datcovg2 <- right_join(position, datcovg, by = "Position")

#define genome positions
genome.df3 <- data.frame(Position = c(1, 679,
                                      680, 1234,
                                      1235, 2347,
                                      2348, 3016,
                                      3017, 3742,
                                      3743, 4150,
                                      4151, 4645, 
                                      4646, 5650, 
                                      5651, 5929, 
                                      5930, 6007, 
                                      6008, 6577,
                                      6578, 7984,
                                      7985, 8667), 
                         Peptide = rep(c("5UTR", "L", "VP0", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D", "3UTR"), each=2))


p1 <- ggplot(datcovg2) + geom_area(aes(x=Position, y=Coverage, fill = Peptide), size=1) +
  geom_ribbon(data=genome.df3, aes(x = Position, ymin=0, ymax=0,fill = Peptide), color="black") + 
  facet_grid() + theme_bw() + xlab("Genome position") + ylab("Coverage (rpm)") + 
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) + 
  ggtitle("OP287812 Coverage") + scale_x_continuous(breaks=c(0,2000/1,4000/1,6000/1, 8000/1), 
                                                    labels = c(0,2000, 4000, 6000, 8000)) + 
  scale_fill_manual(values = wes_palette("Darjeeling2", 13, type = "continuous"))
  #lacroix_palette("PeachPear", n = 13, type = "continuous")

#Darjeeling2
#Moonrise2

ggplot(datcovg2) + geom_area(aes(x=Position, y=Coverage, fill = Peptide), size=1) +
  geom_ribbon(data=genome.df3, aes(x = Position, ymin=0, ymax=0,fill = Peptide), color="black") + 
  facet_grid() + theme_bw() + xlab("Genome position") + ylab("Coverage (rpm)") + 
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) + 
  ggtitle("OP287812 Coverage") + scale_x_continuous(breaks=c(0,2000/1,4000/1,6000/1, 8000/1), 
                                                    labels = c(0,2000, 4000, 6000, 8000)) + 
  scale_fill_manual(values = wes_palette("Darjeeling2", 13, type = "continuous"))
 
ggplot(datcovg2) + geom_area(aes(x=Position, y=Coverage, fill = Peptide), size=1) +
  geom_ribbon(data=genome.df3, aes(x = Position, ymin=0, ymax=0,fill = Peptide), color="black") + 
  facet_grid() + theme_bw() + xlab("Genome position") + ylab("Coverage (rpm)") + 
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) + 
  ggtitle("OP287812 Coverage") + scale_x_continuous(breaks=c(0,2000/1,4000/1,6000/1, 8000/1), 
                                                    labels = c(0,2000, 4000, 6000, 8000)) + 
  scale_fill_manual(values = wes_palette("Moonrise2", 13, type = "continuous"))

ggplot(datcovg2) + geom_area(aes(x=Position, y=Coverage, fill = Peptide), size=1) +
  geom_ribbon(data=genome.df3, aes(x = Position, ymin=0, ymax=0,fill = Peptide), color="black") + 
  facet_grid() + theme_bw() + xlab("Genome position") + ylab("Coverage (rpm)") + 
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) + 
  ggtitle("OP287812 Coverage") + scale_x_continuous(breaks=c(0,2000/1,4000/1,6000/1, 8000/1), 
                                                    labels = c(0,2000, 4000, 6000, 8000)) + 
  scale_fill_manual(values = wes_palette("GrandBudapest1", 13, type = "continuous"))

  
  
fig5 <- cowplot::plot_grid(p1,p2,p3, labels = "AUTO", 
                           label_size = 16, nrow = 3, ncol = 1, 
                           label_colour = "#D32D41")
fig5




#fig5 <- cowplot::plot_grid(p1,p2,p3, nrow=3, ncol = 1, 
    #                       labels= c("(A)", "(B)", "(C)"), 
    #                       label_size = 16, label_x = -.01, 
     #                      rel_heights = c(1,.9,1.2))


#move to long



# cara's code on final fig
#and together
#Fig5top <- cowplot::plot_grid(p1,p2,p3, nrow=3, ncol = 1, labels= c("(A)", "(B)", "(C)"), label_size = 16, label_x = -.01, rel_heights = c(1,.9,1.2))

#Fig5 <- cowplot::plot_grid(Fig5top, leg1, nrow = 1, ncol = 2, rel_widths = c(1,.2))


# Fig5 <- cowplot::plot_grid(Fig5top, leg1, nrow = 2, ncol = 1, rel_heights = c(1,.1))

ggsave(file = paste0(homewd, "/final-figures/Fig5.png"),
       plot=Fig5,
       units="mm",  
       width=95, 
       height=70, 
       #limitsize = F,
       scale=4)#, 