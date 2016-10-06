require(ggplot2)
require(tidyr)
require(gridExtra)

setwd("~/Documents/PhD/Maela carriage length/")

# Read in data and merge. For each snp, the p-value, R^2 with lead snps and consequences from VEP
consequences <- read.table("vep.all.csq", header = T, stringsAsFactors = F)

pvals <- read.table("plink_baps.plot", header=F, stringsAsFactors = F)
pvals <- pvals[,c(3,4)]
colnames(pvals) <- c("Location", "Pval")

merged_data <- merge(pvals, consequences, all.x=T)
merged_data[is.na(merged_data[,3]),3] <- "intergenic_variant"

ld_scores <- read.table("lead_snps_ld.txt", header=F)
colnames(ld_scores) <- c("LeadSNP", "Location", "R2")
ld_scores = spread(ld_scores, LeadSNP, R2) # long to wide format
colnames(ld_scores) = c("Location", "SNP1", "SNP2", "SNP3", "SNP4")

secondary_ld_scores <- read.table("second_lead_snp_ld.txt", header=F)
colnames(secondary_ld_scores) <- c("Pos", "Location", "SNP1_second")
all_ld_scores <- merge(ld_scores, secondary_ld_scores, all = T)

merged_data <- merge(merged_data, all_ld_scores, all.x = T)
non_syn <- merged_data$Consequence != "intergenic_variant" & merged_data$Consequence != "synonymous_variant"
names(non_syn) <- "NonSyn"
merged_data <- cbind(merged_data, non_syn)

genes <- read.table("gene_plots/gene_locations.txt", stringsAsFactors = F)
ypos <- (rep(c(-1,1), nrow(genes)/2))
genes <- data.frame(genes, ypos=ypos)
genes[genes$V1 < genes$V2, 3] = 1
genes[genes$V1 >= genes$V2, 3] = -1

# Manhattan plot for first locus
manhattan = ggplot(merged_data) + 
  geom_point(aes(x=Location, y=Pval, 
                 colour=cut(SNP1,seq(0,1,0.2), include.lowest = T), 
                 shape=factor(non_syn))) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  scale_shape_manual(values=c(16, 3)) +
  scale_color_manual(values=c("dark blue", "light blue", "green", "gold", "red"), drop=F) +
  geom_hline(yintercept = -log(1.4*10^-8,10)) +
  ylab("-log10(P-value)") +
  xlim(789000,820000) + ylim(0,10) + theme(legend.position="none")

# Genes for the locus
genes_plot = ggplot(genes) + 
  geom_segment(aes(x=V1, xend=V2, y=ypos, yend=ypos), 
               arrow = arrow(angle = 30, length = unit(0.10, "npc"), 
                              type="closed"), colour="blue", size=2) + 
  theme_bw(base_size = 16, base_family = "Helvetica") +
  xlim(789000,820000) + ylim(-2,2) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

png(filename = "gene_plots/SNP1_locus.png", width = 1000, height = 700)

gA <- ggplotGrob(manhattan)
gB <- ggplotGrob(genes_plot)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=1, heights=c(4, 1))

dev.off()

# Manhattan plot for first locus
manhattan = ggplot(merged_data) + 
  geom_point(aes(x=Location, y=Pval, 
                 colour=cut(SNP1_second,seq(0,1,0.2), include.lowest = T), 
                 shape=factor(non_syn))) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  scale_shape_manual(values=c(16, 3)) +
  scale_color_manual(values=c("dark blue", "light blue", "green", "gold", "red"), drop=F) +
  geom_hline(yintercept = -log(1.4*10^-8,10)) +
  ylab("-log10(P-value)") +
  xlim(789000,820000) + ylim(0,10) + theme(legend.position="none")

# Genes for the locus
genes_plot = ggplot(genes) + 
  geom_segment(aes(x=V1, xend=V2, y=ypos, yend=ypos), 
               arrow = arrow(angle = 30, length = unit(0.10, "npc"), 
                             type="closed"), colour="blue", size=2) + 
  theme_bw(base_size = 16, base_family = "Helvetica") +
  xlim(789000,820000) + ylim(-2,2) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

png(filename = "gene_plots/SNP1_second_locus.png", width = 1000, height = 700)

gA <- ggplotGrob(manhattan)
gB <- ggplotGrob(genes_plot)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=1, heights=c(4, 1))

dev.off()

# Manhattan plot for second locus
manhattan = ggplot(merged_data) + 
  geom_point(aes(x=Location, y=Pval, 
                 colour=cut(SNP2,seq(0,1,0.2), include.lowest = T), 
                 shape=factor(non_syn))) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  scale_shape_manual(values=c(16, 3)) +
  scale_color_manual(values=c("dark blue", "light blue", "green", "gold", "red"), drop=F) +
  geom_hline(yintercept = -log(1.4*10^-8,10)) +
  ylab("-log10(P-value)") +
  xlim(1152500,1160600) + theme(legend.position="none")

# Genes for the locus
genes_plot = ggplot(genes) + 
  geom_segment(aes(x=V1, xend=V2, y=ypos, yend=ypos), 
               arrow = arrow(angle = 30, length = unit(0.10, "npc"), 
                             type="closed"), colour="blue", size=2) + 
  theme_bw(base_size = 16, base_family = "Helvetica") +
  xlim(1152500,1160600) + ylim(-2,2) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

png(filename = "gene_plots/SNP2_locus.png", width = 1000, height = 700)

gA <- ggplotGrob(manhattan)
gB <- ggplotGrob(genes_plot)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=1, heights=c(4, 1))

dev.off()

# Manhattan plot for third locus
manhattan = ggplot(merged_data) + 
  geom_point(aes(x=Location, y=Pval, 
                 colour=cut(SNP3,seq(0,1,0.2), include.lowest = T), 
                 shape=factor(non_syn))) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  scale_shape_manual(values=c(16, 3)) +
  scale_color_manual(values=c("dark blue", "light blue", "green", "gold", "red"), drop=F) +
  geom_hline(yintercept = -log(1.4*10^-8,10)) +
  ylab("-log10(P-value)") +
  xlim(2034000,2074500) + ylim(0,12.5) + theme(legend.position="none")

# Genes for the locus
genes_plot = ggplot(genes) + 
  geom_segment(aes(x=V1, xend=V2, y=ypos, yend=ypos), 
               arrow = arrow(angle = 30, length = unit(0.10, "npc"), 
                             type="closed"), colour="blue", size=2) + 
  theme_bw(base_size = 16, base_family = "Helvetica") +
  xlim(2034000,2074500) + ylim(-2,2) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

png(filename = "gene_plots/SNP3_locus.png", width = 1000, height = 700)

gA <- ggplotGrob(manhattan)
gB <- ggplotGrob(genes_plot)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=1, heights=c(4, 1))

dev.off()

# Manhattan plot for fourth locus
manhattan = ggplot(merged_data) + 
  geom_point(aes(x=Location, y=Pval, 
                 colour=cut(SNP4,seq(0,1,0.2), include.lowest = T), 
                 shape=factor(non_syn))) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  scale_shape_manual(values=c(16, 3)) +
  scale_color_manual(values=c("dark blue", "light blue", "green", "gold", "red"), drop=F) +
  geom_hline(yintercept = -log(1.4*10^-8,10)) +
  ylab("-log10(P-value)") +
  xlim(2185000,2192000) + ylim(0,11.5) + theme(legend.position="none")

# Genes for the locus
genes_plot = ggplot(genes) + 
  geom_segment(aes(x=V1, xend=V2, y=ypos, yend=ypos), 
               arrow = arrow(angle = 30, length = unit(0.10, "npc"), 
                             type="closed"), colour="blue", size=2) + 
  theme_bw(base_size = 16, base_family = "Helvetica") +
  xlim(2185000,2192000) + ylim(-2,2) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

png(filename = "gene_plots/SNP4_locus.png", width = 1000, height = 700)

gA <- ggplotGrob(manhattan)
gB <- ggplotGrob(genes_plot)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=1, heights=c(4, 1))

dev.off()

# Overall manhattan plot
png(filename = "gene_plots/whole_manhattan.png", width = 2000, height = 1000)
ggplot(merged_data) + 
  geom_point(aes(x=Location, y=Pval)) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  geom_hline(yintercept = -log(1.4*10^-8,10)) +
  ylab("-log10(P-value)") +
  theme(legend.position="none")
dev.off()