#### init ####

library(phyloseq)
library(ggplot2)

source("/Users/mdavi/Documents/R_functions/Rscripts/Rscripts/rename_taxa.R")
source("/Users/mdavi/Documents/R_functions/Rscripts/Rscripts/calculate_adonis2.R")
source("/Users/mdavi/Documents/R_functions/Rscripts/Rscripts/taxa_facet_barplot_asv.R")

colours <-  c(
  "#F0A3FF",
  "#0075DC",
  "#993F00",
  "#4C005C",
  "#2BCE48",
  "#FFCC99",
  "#808080",
  "#94FFB5",
  "#8F7C00",
  "#9DCC00",
  "#C20088",
  "#003380",
  "#FFA405",
  "#FFA8BB",
  "#426600",
  "#FF0010",
  "#5EF1F2",
  "#00998F",
  "#740AFF",
  "#990000",
  "#FFFF00"
)

#### true asvs (full) ####

ps <- readRDS("all_ITS/ps.all_ITS.2018-11-01.RDS")
ps <- rename_taxa(ps, add_ASV_2tax = F)
ps.dsmz <- prune_samples(ps@sam_data$Project %in% c("DSMZ"),ps)
ps.dsmz <- prune_samples(ps.dsmz@sam_data$primer_set %in% c("FULL"),ps.dsmz)
ps.dsmz <- prune_samples(ps.dsmz@sam_data$Library=="ITS_0012",ps.dsmz)
ps.dsmz <- prune_taxa(taxa_sums(ps.dsmz)!=0, ps.dsmz)

# write.csv(ps.dsmz@otu_table, "mock_members_otu_table.csv")
# write.csv(ps.dsmz@tax_table, "mock_members_tax_table.csv")

ps.dsmz.t <- transform_sample_counts(ps.dsmz, function(x) x/sum(x))
ps.dsmz.t

p <- taxa_facet_barplot_asv(ps = ps.dsmz.t, Group = "Sample_Name", facet1 = "Sample_Name", facet2 = "Library", rank = "ASV", lumpNA = T)
p + labs(title="Pure amplicons") + theme(panel.spacing = unit(0.4, "lines")) + theme(text = element_text(size = 10),
                                                                                          axis.text.x = element_text(angle = 90, hjust = 1))

x <- c()
asvs <- list()
for(i in 1:12){x[i] <- sum(unname(ps.dsmz.t@otu_table[i,order(ps.dsmz.t@otu_table[i,], decreasing = T)[1:10]])>0.01)
asvs[[i]] <- taxa_names(ps.dsmz.t)[order(ps.dsmz.t@otu_table[i,], decreasing = T)[1:x[i]]]
}

true_asvs <- unique(unlist(asvs))
ps@tax_table[true_asvs,1:7]
true_asvs <- true_asvs[true_asvs!="ASV_10094"] # primer dimer / FP

# Torulaspora delbrueckii is missing from the pure amplicon libraries
# distill these from the pool

ps.Td <- prune_samples(ps@sam_data$Sample_Name %in% grep("full", grep("1-1-poo", ps@sam_data$Sample_Name, value = T), value = T), ps)
ps.Td <- prune_taxa(taxa_sums(ps.Td)!=0, ps.Td)

ps.Td@otu_table[,ps.Td@tax_table[,7] %in% "s__delbrueckii"]
true_asvs <- unique(c(true_asvs,"ASV_26504"))

cbind(ps@tax_table[true_asvs,1:7],nchar(ps@tax_table[true_asvs,8]),ps@tax_table[true_asvs,8])

#### pool ITS_0009 ####

ps <- readRDS("all_ITS/ps.all_ITS.2018-11-01.RDS")
ps <- rename_taxa(ps, add_ASV_2tax = F)
ps.pool_9 <- prune_samples(ps@sam_data$Library=="ITS_0009",ps)
ps.pool_9 <- prune_samples(ps.pool_9@sam_data$Sample_Name %in% grep("pool",ps.pool_9@sam_data$Sample_Name, value = T), ps.pool_9)
ps.pool_9 <- prune_taxa(taxa_sums(ps.pool_9)!=0, ps.pool_9)
ps.pool_9@sam_data$Dilution <- -1000
ps.pool_9@sam_data$Dilution[grep("-2",ps.pool_9@sam_data$Sample_Name)] <- -999
ps.pool_9@sam_data$Dilution[grep("0-1",ps.pool_9@sam_data$Sample_Name)] <- 10
ps.pool_9@sam_data$Dilution[grep("1-0",ps.pool_9@sam_data$Sample_Name)] <- -10
ps.pool_9@sam_data$Dilution[grep("1000-1",ps.pool_9@sam_data$Sample_Name)] <- -3
ps.pool_9@sam_data$Dilution[grep("10-1",ps.pool_9@sam_data$Sample_Name)] <- -1
ps.pool_9@sam_data$Dilution[grep("1-1",ps.pool_9@sam_data$Sample_Name)] <- 0
ps.pool_9@sam_data$Dilution[grep("1-10",ps.pool_9@sam_data$Sample_Name)] <- 1
ps.pool_9@sam_data$Dilution <- as.numeric(ps.pool_9@sam_data$Dilution )
ps.pool_9@sam_data$Dilution <- as.factor(ps.pool_9@sam_data$Dilution)

p <- taxa_facet_barplot_asv(ps = ps.pool_9, Group = "Polymerase", facet1 = "Dilution", facet2 = "primer_set", rank = "Genus", lumpNA = T)
p$data$facet1 <- as.factor(as.numeric(p$data$facet1))
levels(p$data$facet1)[levels(p$data$facet1)=="-1000"] <- "NC1"
levels(p$data$facet1)[levels(p$data$facet1)=="-999"] <- "NC2"
levels(p$data$facet1)[levels(p$data$facet1)=="-10"] <- "Feces"
levels(p$data$facet1)[levels(p$data$facet1)=="10"] <- "Mock"
p + labs(title="Genus level composition mock samples") + theme(panel.spacing = unit(0.4, "lines")) + theme(text = element_text(size = 10),
                                                                                                           axis.text.x = element_text(angle = 90, hjust = 1))
ps.pool_9.bits <- prune_samples(ps.pool_9@sam_data$primer_set=="BITS", ps.pool_9)
ps.pool_9.full <- prune_samples(ps.pool_9@sam_data$primer_set=="FULL", ps.pool_9)  

p <- taxa_facet_barplot_asv(ps = ps.pool_9.bits, Group = "Polymerase", facet1 = "Dilution", facet2 = "primer_set", rank = "ASV", lumpNA = T)
p$data$facet1 <- as.factor(as.numeric(p$data$facet1))
p + theme(text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 1))

levels(p$data$facet1)[levels(p$data$facet1)=="-1000"] <- "NC1"
levels(p$data$facet1)[levels(p$data$facet1)=="-999"] <- "NC2"
levels(p$data$facet1)[levels(p$data$facet1)=="-10"] <- "Feces"
levels(p$data$facet1)[levels(p$data$facet1)=="10"] <- "Mock"
p + labs(title="Genus level composition mock samples") + theme(panel.spacing = unit(0.4, "lines")) + theme(text = element_text(size = 10),
                                                                                                           axis.text.x = element_text(angle = 90, hjust = 1))

p <- taxa_facet_barplot_asv(ps = ps.pool_9.full, Group = "Polymerase", facet1 = "Dilution", facet2 = "primer_set", rank = "ASV", lumpNA = T)
p$data$facet1 <- as.factor(as.numeric(p$data$facet1))
levels(p$data$facet1)[levels(p$data$facet1)=="-1000"] <- "NC1"
levels(p$data$facet1)[levels(p$data$facet1)=="-999"] <- "NC2"
levels(p$data$facet1)[levels(p$data$facet1)=="-10"] <- "Feces"
levels(p$data$facet1)[levels(p$data$facet1)=="10"] <- "Mock"
p + labs(title="Genus level composition mock samples") + theme(panel.spacing = unit(0.4, "lines")) + theme(text = element_text(size = 10),
                                                                                                           axis.text.x = element_text(angle = 90, hjust = 1))

ps.pool_9.full.t <- transform_sample_counts(ps.pool_9.full, function(x) x/sum(x))
ps.pool_9.bits.t <- transform_sample_counts(ps.pool_9.bits, function(x) x/sum(x))

ps.pool_9.full.t.tv <- prune_taxa(taxa_names(ps.pool_9.full.t) %in% true_asvs, ps.pool_9.full.t)

psm <- psmelt(ps.pool_9.full.t.tv)
psm$OTU2 <- paste0(psm$OTU, psm$Genus, psm$Species)
ggplot(psm, aes(x=Polymerase, y=Abundance, color=OTU2)) + geom_point() +geom_line(aes(group=OTU2)) + 
 facet_grid(~Dilution) +
 scale_y_log10() + 
 #scale_color_manual(values=colours) + 
  geom_label(aes(label=OTU))

ord <- ordinate(rarefy_even_depth(ps.pool_9, rngseed = 211202), "PCoA", "bray")
plot_ordination(physeq = ps.pool_9, ordination = ord, color = "Polymerase", label="Sample_Name" ) + facet_wrap(~primer_set)


ord <- ordinate(ps.pool_9.full.t, "PCoA", "bray")
plot_ordination(ps.pool_9.full.t, ordination = ord, color = "Polymerase", label="Sample_Name" ) + facet_wrap(~primer_set)

ord <- ordinate(ps.pool_9.bits.t, "PCoA", "bray")
plot_ordination(ps.pool_9.bits.t, ordination = ord, color = "Polymerase", label="Sample_Name" ) + facet_wrap(~primer_set)


#### theoretical composition ####

qPCR <- read.csv("dsmz_qPCR.csv", header = T, sep=";")
qPCR

ps.pool_9_genus_mocks <- tax_glom(prune_samples(as.numeric(as.character(ps.pool_9@sam_data$Dilution))>-4, ps.pool_9), taxrank = "Genus", NArm = F)
psm <- psmelt(transform_sample_counts(ps.pool_9_genus_mocks, function(x) x/sum(x)))

colnames(psm)
psm2 <- psm[,c("Abundance","Sample_Name","primer_set","Polymerase","Dilution","Genus")]
df <- data.frame(Abundance=qPCR$N0_.indiv_eff./sum(qPCR$N0_.indiv_eff.),
                 Sample_Name="Theo",
                 primer_set="FULL",
                 Polymerase="NONE",
                 Dilution="Theoretical",
                 Genus=paste0("g__",qPCR$Genus ))

psm3 <- data.frame(Abundance=c(psm2$Abundance, df$Abundance),
                   Sample_Name=c(as.character(psm2$Sample_Name), as.character(df$Sample_Name)),
                   primer_set=c(as.character(psm2$primer_set), as.character(df$primer_set)),
                   Polymerase=c(as.character(psm2$Polymerase), as.character(df$Polymerase)),
                   Dilution=c(as.character(psm2$Dilution), as.character(df$Dilution)),
                   Genus=c(as.character(psm2$Genus), as.character(df$Genus)))

psm3$Genus[!psm3$Genus %in% c(levels(df$Genus),"g__Vishniacozyma","g__Vishniacozyma","g__Torulaspora")] <- NA
psm3$Genus <- droplevels(psm3$Genus)
levels(psm3$Genus)[levels(psm3$Genus)=="g__Cryptococcus"] <- "g__Vishniacozyma"

ggplot(psm3, aes(x=Polymerase, y=Abundance, fill=Genus)) + geom_bar(stat="identity") + 
  facet_grid(primer_set~Dilution, scales="free_x", space="free_x") + 
  scale_fill_manual(values=colours) + 
  NULL


psm4 <- psm3[!is.na(psm3$Genus),]
psm4 <- psm4[psm4$Sample_Name!="Theo",]

rownames(qPCR) <- paste0("g__",qPCR$Genus)
rownames(qPCR)[rownames(qPCR)=="g__Cryptococcus"] <- "g__Vishniacozyma"

psm4$Theoretical <- qPCR[as.character(psm4$Genus),]$N0_.indiv_eff.
psm4$Theoretical[is.na(psm4$Theoretical)] <- 0

ggplot(psm4, aes(x=Theoretical+0.0001, y=Abundance+0.0001, group=Sample_Name, color=Genus)) + 
  geom_point(size=3) + 
  facet_wrap(~Sample_Name) + 
  scale_y_log10() + 
  scale_x_log10() + 
  scale_color_manual(values=colours) + 
  geom_smooth(method='lm',  formula = y~x) + 
  NULL

psm5 <- psm4
psm5 <- psm5[!psm5$Genus %in% c("g__Torulaspora","g__Aureobasidium","g__Rhodotorula"),] # wrong qPCR or no data
psm5 <- psm5[!(psm5$Genus=="g__Vishniacozyma" & psm5$Abundance<0.01),] # these are cryptococcus classified
psm5$f1  <- paste0(psm5$Polymerase, "_", psm5$primer_set)

ggplot(psm5, aes(x=Theoretical, y=Abundance, group=Sample_Name, color=Genus)) + 
  geom_point(size=3) + 
  #  facet_wrap(~Sample_Name) + 
  facet_grid(Dilution~f1) + 
#  scale_y_log10() + 
#  scale_x_log10() + 
  scale_color_manual(values=colours) + 
  geom_smooth(method='lm',  formula = y~x) + 
  #  xlim(c(0.00001,2)) + 
  #  ylim(c(0.00001,2)) + 
  #ggrepel::geom_label_repel(aes(label=Genus)) + 
  NULL


ggplot(psm5, aes(x=Theoretical, y=Abundance, group=Sample_Name, color=Genus)) + 
  geom_point(size=3) + 
#  facet_wrap(~Sample_Name) + 
  facet_grid(Dilution~f1) + 
  scale_y_log10() + 
  scale_x_log10() + 
  scale_color_manual(values=colours) + 
  geom_smooth(method='lm',  formula = y~x) + 
  #ggrepel::geom_label_repel(aes(label=Genus)) + 
  NULL

psm6 <- psm5[psm5$Sample_Name=="0-1-pool-Phu-BITS",]

# qPCR$nchar_ITS1 <- "X"
# qPCR[qPCR$DSMZ=="1333",]$nchar_ITS1 <- 550
# qPCR[qPCR$DSMZ=="27152",]$nchar_ITS1 <- 390
# qPCR[qPCR$DSMZ=="3847",]$nchar_ITS1 <- 380 #~ 302 amplicon
# qPCR[qPCR$DSMZ=="1079",]$nchar_ITS1 <- 380 #~ 302 amplicon
# qPCR[qPCR$DSMZ=="27485",]$nchar_ITS1 <- 370 #~ 302 amplicon
# qPCR[qPCR$DSMZ=="70403",]$nchar_ITS1 <- 340 #~ 270 amplicon
# qPCR[qPCR$DSMZ=="27774",]$nchar_ITS1 <- 340 #~ 251 amplicon
# qPCR[qPCR$DSMZ=="3433",]$nchar_ITS1 <- 300 #~ 221 amplicon
# 
# 
# qPCR[qPCR$DSMZ=="1386",]$nchar_ITS1 <- 257+80
# qPCR[qPCR$DSMZ=="1623",]$nchar_ITS1 <- 302+80
# qPCR[qPCR$DSMZ=="6170",]$nchar_ITS1 <- 302+80
# qPCR[qPCR$DSMZ=="6381",]$nchar_ITS1 <- 286+80

ggplot(psm6, aes(x=Theoretical, y=Abundance, group=Sample_Name, color=Genus)) + 
  geom_point(size=3) + 
  #  facet_wrap(~Sample_Name) + 
  facet_grid(Dilution~f1) + 
  scale_y_log10() + 
  scale_x_log10() + 
  scale_color_manual(values=colours) + 
  geom_smooth(method='lm',  formula = y~x) + 
  #ggrepel::geom_label_repel(aes(label=Genus)) + 
  NULL

# #### mod ps ####
# 
# ps <- readRDS("all_ITS/ps.all_ITS.2018-11-01.RDS")
# 
# ps.dsmz <- prune_samples(ps@sam_data$Project %in% c("Pool","DSMZ"),ps)
# ps.dsmz <- prune_samples(ps@sam_data$primer_set %in% c("FULL"),ps.dsmz)
# #ps.dsmz <- prune_samples(ps.dsmz@sam_data$Library=="ITS_0012",ps.dsmz)
# ps.dsmz <- prune_taxa(taxa_sums(ps.dsmz)>0,ps.dsmz)
# ps.dsmz@sam_data$Sample_Name <- make.unique(ps.dsmz@sam_data$Sample_Name)
# 
# 
# ps.dsmz.t <- transform_sample_counts(ps.dsmz, function(x) x/sum(x))
# ps.dsmz.t <- rename_taxa(ps.dsmz.t, add_ASV_2tax = T)
# 
# top21 <- names(sort(taxa_sums(ps.dsmz.t), decreasing = T)[1:40])
# 
# 
# psm <- psmelt(prune_taxa(taxa_names(ps.dsmz.t) %in% top21, ps.dsmz.t))
# psm$ASV2 <- psm$ASV
# levels(psm$ASV2) <- make.unique(substr(levels(psm$ASV2),1,10))
# 
# psm$ASV3 <- paste0(psm$Genus, psm$Species)
# 
# ggplot(psm, aes(x=Sample_Name, y=Abundance, fill=ASV)) + geom_bar(stat="identity") + facet_wrap(~primer_set)
# ggplot(psm, aes(x=Sample_Name, y=Abundance, fill=ASV)) + geom_bar(stat="identity") + facet_wrap(~Library, scales="free_x")
# 
# ggplot(psm, aes(x=Sample_Name, y=Abundance, fill=ASV2)) + geom_bar(stat="identity") + facet_wrap(~Library)
# ggplot(psm, aes(x=Sample_Name, y=Abundance, fill=Species)) + geom_bar(stat="identity") + facet_wrap(~Library)
# ggplot(psm, aes(x=Sample_Name, y=Abundance, fill=Genus)) + geom_bar(stat="identity") + facet_wrap(~Library)
# 
# 
# #### FULL TRUE Set ####
# 
# ps <- readRDS("all_ITS/ps.all_ITS.2018-11-01.RDS")
# ps.dsmz <- prune_samples(ps@sam_data$Project %in% c("DSMZ"),ps)
# ps.dsmz <- prune_samples(ps.dsmz@sam_data$primer_set %in% c("FULL"),ps.dsmz)
# ps.dsmz <- prune_samples(ps.dsmz@sam_data$Library=="ITS_0012",ps.dsmz)
# ps.dsmz <- prune_taxa(taxa_sums(ps.dsmz)!=0, ps.dsmz)
# 
# write.csv(ps.dsmz@otu_table, "mock_members_otu_table.csv")
# write.csv(ps.dsmz@tax_table, "mock_members_tax_table.csv")
# 
# 
# ps.dsmz.t <- transform_sample_counts(ps.dsmz, function(x) x/sum(x))
# ps.dsmz.t
# 
# x <- c()
# asvs <- list()
# for(i in 1:12){x[i] <- sum(unname(ps.dsmz.t@otu_table[i,order(ps.dsmz.t@otu_table[i,], decreasing = T)[1:10]])>0.05)
# asvs[[i]] <- taxa_names(ps.dsmz.t)[order(ps.dsmz.t@otu_table[i,], decreasing = T)[1:x[i]]]
# }
# 
# true_asvs <- unique(unlist(asvs))
# true_asvs <- c("CTTGGTCATTTAGAGGAACTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTAGAGAAATCTATATGAATGAAGTTAGAGGACGTCTAAAGATACTGTAAGAGAGGATCAGGTTCAAGACCAGCGCTTAATTGCGCGATACGTCTTGTGCGTGCTTCCCAGAGGTGACAAACACAAACAACTTTTTATTATTATAAACCAGTCAAAACCAATTTCGTTATGAAATTAAAAATATTTAAAACTTTCAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGC",
#                true_asvs)
# 
# #### ps.dsmz.pool
# 
# ps <- readRDS("all_ITS/ps.all_ITS.2018-11-01.RDS")
# 
# ps.dsmz <- prune_samples(ps@sam_data$Project %in% c("Pool","DSMZ"),ps)
# ps.dsmz <- prune_samples(ps.dsmz@sam_data$primer_set %in% c("FULL"),ps.dsmz)
# ps.dsmz <- prune_samples(ps.dsmz@sam_data$Library!="ITS_0010",ps.dsmz)
# ps.dsmz <- prune_samples(ps.dsmz@sam_data$Library!="ITS_0008",ps.dsmz)
# ps.dsmz <- prune_taxa(taxa_sums(ps.dsmz)>0,ps.dsmz)
# ps.dsmz@sam_data$Sample_Name <- make.unique(ps.dsmz@sam_data$Sample_Name)
# 
# ps.dsmz.t <- transform_sample_counts(ps.dsmz, function(x) x/sum(x))
# #ps.dsmz.t <- rename_taxa(ps.dsmz.t, add_ASV_2tax = T)
# 
# ps.dsmz.t <- prune_taxa(taxa_sums(ps.dsmz.t)>0.001,ps.dsmz.t)
# ps.dsmz.t
# 
# ps.dsmz.t <- transform_sample_counts(ps.dsmz, function(x) x/sum(x))
# ps.dsmz.t.tv <- prune_taxa(taxa_names(ps.dsmz.t) %in% true_asvs, ps.dsmz.t)
# 
# psm <- psmelt(prune_taxa(taxa_names(ps.dsmz.t) %in% top21, ps.dsmz.t))
# psm$ASV2 <- psm$ASV
# levels(psm$ASV2) <- make.unique(substr(levels(psm$ASV2),1,10))
# 
# psm$ASV3 <- paste0(psm$Genus, psm$Species, psm$ASV2)
# 
# pdf("ordination_dsmz_pool.pdf")
# 
# ps.dsmz.t <- transform_sample_counts(ps.dsmz, function(x) x/sum(x))
# ps.dsmz.t.tv <- prune_taxa(taxa_names(ps.dsmz.t) %in% true_asvs, ps.dsmz.t)
# ps.dsmz.t.fv <- prune_taxa(!taxa_names(ps.dsmz.t) %in% true_asvs, ps.dsmz.t)
# 
# psm <- psmelt(ps.dsmz.t.tv)
# psm$ASV2 <- psm$ASV
# levels(psm$ASV2) <- make.unique(substr(levels(psm$ASV2),1,10))
# psm$ASV3 <- paste0(psm$Genus, psm$Species, psm$ASV2)
# 
# psm$TP <- psm$ASV
# levels(psm$TP) <- levels(psm$TP) %in% true_asvs
# 
# taxa_facet_barplot_asv(rename_taxa(ps.dsmz.t.tv), Group = "Sample_Name",facet2 = "Project", facet1 = "Library", rank = "ASV", lumpNA = T)
# 
# ggplot(psm, aes(x=Sample_Name, y=Abundance, fill=ASV3)) + 
#   geom_bar(stat="identity") + 
#   facet_wrap(Project~Library, scales="free_x") + 
#   theme(text = element_text(size = 10),axis.text.x = element_text(angle = 90, hjust = 1))
# 
# ps.dsmz.t <- transform_sample_counts(ps.dsmz, function(x) x/sum(x))
# ps.dsmz.t@sam_data$Dilution <- gsub("water","",gsub("_.water.*","",gsub("pool","",gsub("-pool.*","",gsub("DSMZ.*","",ps.dsmz.t@sam_data$Sample_Name)))))
# ps.dsmz.t@sam_data$Dilution <- gsub("\\..*","",ps.dsmz.t@sam_data$Dilution)
# 
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution==""] <- NA
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1000-1"] <- 3
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1-1"] <- 0
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="10-1"] <- 1
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1-10"] <- -1
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="0-1"] <- -10
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1-0"] <- 10
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1000"] <- 3
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="10000"] <- 4
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="100000"] <- 5
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1000000"] <- 6
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="10000000"] <- 7
# 
# 
# # di-do 15:30 of 17:00
# # zo 7:45
# 
# 
# ps.dsmz.t@sam_data$Dilution <- as.numeric(ps.dsmz.t@sam_data$Dilution)
# ord <- ordinate(ps.dsmz.t, "PCoA", "bray")
# plot_ordination(ps.dsmz.t, ord, label="Sample_Name", color="Dilution") + scale_color_gradient(low="red", high="green")
# dev.off()
# 
# 
# length <- aggregate(t(ps.dsmz@otu_table), by=list(nchar(taxa_names(ps.dsmz))), FUN=sum)
# 
# rownames(length) <- length[,1]
# length <- length[,-1]
# length[,1]
# 
# df <- data.frame(ps.dsmz@sam_data)
# 
# df$length <- barplot(apply(length, 2, function(x) sum(x*as.numeric(rownames(length)))), las=2)
# 
# ggplot(df, aes(x=Sample_Name, y=length)) + geom_point() + facet_wrap(~Polymerase, scales="free_x")
# 
# 
# 
# #### observations ####
# 
# # loss of small amplicons ITS0012 pool1000; possibly due to ampure purification
# # difference in composition between Taq/FULL
# 
# ps <- readRDS("all_ITS/ps.all_ITS.2018-11-01.RDS")
# 
# ps.dsmz <- prune_samples(ps@sam_data$Project %in% c("Pool","DSMZ"),ps)
# ps.dsmz <- prune_samples(ps.dsmz@sam_data$primer_set %in% c("FULL"),ps.dsmz)
# ps.dsmz <- prune_samples(ps.dsmz@sam_data$Library!="ITS_0010",ps.dsmz)
# ps.dsmz.ITS_0009 <- prune_samples(ps.dsmz@sam_data$Library=="ITS_0009",ps.dsmz)
# 
# ps.dsmz.ITS_0009@sam_data$Sample_Name <- make.unique(ps.dsmz.ITS_0009@sam_data$Sample_Name)
# ps.dsmz.t <- transform_sample_counts(ps.dsmz.ITS_0009, function(x) x/sum(x))
# ps.dsmz.t <- transform_sample_counts(ps.dsmz, function(x) x/sum(x))
# ps.dsmz.t@sam_data$Dilution <- gsub("water","",gsub("_.water.*","",gsub("pool","",gsub("-pool.*","",gsub("DSMZ.*","",ps.dsmz.t@sam_data$Sample_Name)))))
# ps.dsmz.t@sam_data$Dilution <- gsub("\\..*","",ps.dsmz.t@sam_data$Dilution)
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution==""] <- NA
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1000-1"] <- 3
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1-1"] <- 0
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="10-1"] <- 1
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1-10"] <- -1
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="0-1"] <- -10
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1-0"] <- 10
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1000"] <- 3
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="10000"] <- 4
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="100000"] <- 5
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="1000000"] <- 6
# ps.dsmz.t@sam_data$Dilution[ps.dsmz.t@sam_data$Dilution=="10000000"] <- 7
# 
# 
# ps.dsmz.t.tv <- prune_taxa(taxa_names(ps.dsmz.t) %in% true_asvs, ps.dsmz.t)
# 
# psm <- psmelt(rename_taxa(ps.dsmz.t.tv, add_ASV_2tax = T))
# psm$ASV2 <- psm$ASV
# levels(psm$ASV2) <- make.unique(substr(levels(psm$ASV2),1,10))
# psm$ASV3 <- paste0(psm$Genus, psm$Species, "_", psm$OTU)
# 
# psm$TP <- psm$ASV
# levels(psm$TP) <- levels(psm$TP) %in% true_asvs
# 
# psm <- psm[!is.na(psm$Dilution), ]
# psm <- psm[grep("full",psm$Sample_Name),]
# psm$Dilution <- as.numeric(psm$Dilution)
# 
# ggplot(psm[!psm$Dilution %in% c(4,5,6,7),], aes(x=Polymerase, y=Abundance, color=ASV3)) + geom_point() +geom_line(aes(group=ASV)) + 
#   facet_grid(~Dilution) +
#   scale_y_log10() + 
#   scale_color_manual(values=colours) + geom_label(aes(label=OTU))
# 
# psm$length <- nchar(psm$OTU)
# 
# ggplot(psm[!psm$Dilution %in% c(4,5,6,7),], aes(x=Polymerase, y=Abundance, color=length)) + geom_point() + 
#   geom_line(aes(group=ASV)) + facet_grid(~Dilution) +
#   scale_color_gradient(low="red", high="green") +
#   scale_y_log10()
# 
# 
# taxa_facet_barplot_asv(rename_taxa(ps.dsmz.t.tv), Group = "Sample_Name",facet2 = "Project", facet1 = "Library", rank = "ASV", lumpNA = T)
# 
# ggplot(psm, aes(x=Sample_Name, y=Abundance, fill=ASV3)) + 
#   geom_bar(stat="identity") + 
#   facet_wrap(Project~Library, scales="free_x") + 
#   theme(text = element_text(size = 10),axis.text.x = element_text(angle = 90, hjust = 1))
