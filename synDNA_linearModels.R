#Task: annotate the script for my own understanding. 
#Primary goal is to understand the two key equations for the absolute quantification.
#eq 1, for each bacterial genome: readWeight = intercept + slope * normalizedReads #intercept: unitless. slope: log(counts/million)
#eq 2, for each bacterial genome: numCells = (readWeight * avogadro)/(genomeLength * molWeightBP) #Andy's question was: what are the units?


################
rm(list = ls())
#Install/set up packages
# List of required packages
packages <- c(
  "tidyverse", "data.table", "reshape2", "ggplot2", "RColorBrewer", "vegan",
  "plyr", "dplyr", "ggpubr", "Hmisc", "corrplot", "ggpmisc", "scales",
  "broom", "FactoMineR", "factoextra"
)

# Install any missing packages
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))
rm(packages, new_packages)

#Loading files (count table, metadata, synDNAs dilutions)
setwd("/Users/aal88/Documents/GitHub/SynDNA")
synDNAcounts <- read.delim("data/synDNA_Fwd_Rev_sam.biom.tsv", h=T) #each column is the number of reads aligned to each synDNA (the rows) for each (metagenomic sample, pool, and read direction) (the columns)
metadata <- read.delim("data/synDNA_metadata_updated.tsv", h=T) #contains the total number of reads for each sequencing run (metagenomic sample + pool) and read direction(forward/reverse)
dilutions <- read.delim("data/dilutions_metadata_pool_saliva.tsv") #contains the dilution factor for each of the synDNAs

#Calculation total number of reads aligned to the synDNAs

TotalAlign <- data.frame(apply(synDNAcounts[,-1], 2, sum))
TotalAlign <- data.frame(Files = rownames(TotalAlign), TotalAlign = TotalAlign$apply.synDNAcounts....1...2..sum.)


#Matrix to data.frame and merging with total number of aligned reads

synDNApoolsM <- melt(synDNAcounts) #reframes the data so that all numbers of aligned reads are in the same column, with metadata of the aligned synDNA and sample in the other columns
colnames(synDNApoolsM) <- c("SynDNAName", "Files", "Counts")
synDNApoolsMelt <- merge(synDNApoolsM, TotalAlign, by="Files") #combines the number of aligned reads for each synDNA + sample combo with the total number of reads aligned to any synDNA in the sample


#Merge table with metadata

synDNAMeta <- merge(synDNApoolsMelt, metadata, by = "Files") 

#Aggregate counts from same library, but different lane (Novaseq lanes)

synDNAMetaAg <- aggregate(cbind(Counts,TotalReads,TotalAlign)~Strand+ID+SynDNAName+Pool, synDNAMeta, sum) #not clear this did anything, because there are 420 rows before and after

# Calculate CPM and percentage
synDNAMetaAg$CPM <- (synDNAMetaAg$Counts / synDNAMetaAg$TotalReads) * 1e6 #calculates the percentage of reads aligned to a given synDNA over all reads aligned to synDNAs and multipies by a million
synDNAMetaAg$Percentage <- (synDNAMetaAg$Counts / synDNAMetaAg$TotalReads) / 100 #calculates the percentage of reads aligned to a given synDNA over all reads aligned to synDNAs and multipies by a hundred

# Log Transformation
synDNAMetaAg$CPMlog <- log10(synDNAMetaAg$CPM)

# Merge with dilutions metadata
synDNAMetaAgM <- merge(synDNAMetaAg, dilutions, by = c("SynDNAName", "Pool"))

# Remove samples with low quality and counts
synDNAMetaAgM <- synDNAMetaAgM[-grep("J1_", synDNAMetaAgM$ID), ]

#Visualizing the serial dilution

options(repr.plot.width=20, repr.plot.height=20)

ggplot(synDNAMetaAgM, aes(CPMlog, Dilution2)) +
  geom_point(shape = 15, size = 3) +
  #geom_errorbar(aes(xmin=CPMMeanlog-CPMSdlog, xmax=CPMMeanlog+CPMSdlog), width=.3,position=position_dodge(0.05)) +
  facet_wrap(~ Pool+Strand+ID, ncol = 6) +
  geom_smooth(method = "lm", color="red", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, aes(label = ..rr.label..), 
               parse=TRUE, label.x.npc = "right",, size = 5) +
  stat_fit_glance(method = 'lm', geom = 'text', aes(label = paste0('p = ', format(..p.value.., 3))), 
                  label.x = 3.3, label.y = 5, size = 5) + 
  labs(y = expression(bold(paste("Dilutions ", 10^"-n"))), x = "log10(Count Per Million)") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_bw() + guides(fill=guide_legend(ncol=1))  +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 16, face = "bold", colour = "black"),
        axis.text.y = element_text(angle = 0, size = 16, face = "bold", colour = "black"),
        axis.title=element_text(size=20, face = "bold", colour = "black"),
        strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
        strip.background = element_rect(fill="white"),
        legend.title=element_blank())



#Adjusting the libraries dilution. SynDNA pool follow the ratio 1:20 (5% proportion). 

DilutionLog <- synDNAMetaAgM$Dilution2
DilutionLog[which(DilutionLog == 1)] <- (1/20)  #So because the undiluted pool is 5% of the total concentration of the metagenomic sample, this is 5%
DilutionLog[which(DilutionLog == 2)] <- (0.1/20) #Because the next dilution in the series is 1:10, this is dilutionLog 1 /10
DilutionLog[which(DilutionLog == 3)] <- (0.01/20)
DilutionLog[which(DilutionLog == 4)] <- (0.001/20)
DilutionLog[which(DilutionLog == 5)] <- (0.0001/20)

synDNAMetaAgM$DilutionLog <- DilutionLog

options(repr.plot.width=18, repr.plot.height=20)

ggplot(synDNAMetaAgM, aes(CPMlog, log10((DilutionLog)))) +
  geom_point(size = 3) +
  #geom_errorbar(aes(xmin=log10(CPMMean)-log10(CPMStd), xmax=log10(CPMMean)+log10(CPMStd)), 
  #                  width=.3,position=position_dodge(0.05)) +
  facet_wrap(~ Pool+Strand+ID, ncol = 6) +
  geom_smooth(method = "lm", formula = y ~ x, se=FALSE) +
  stat_poly_eq(formula = y ~ x, aes(label = ..rr.label..), 
               parse=TRUE, label.x.npc = "left", size = 4) +
  stat_fit_glance(method = 'lm', geom = 'text', aes(label = paste0('p = ', format(..p.value.., 3))), 
                  label.x = 1.5, label.y = -1.8, size = 4) + 
  labs(y = expression(bold(paste("Dilutions ", 10^"n"))), x = "log10(CPM)") +  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_bw() + guides(fill=guide_legend(ncol=1))  +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14, face = "bold", colour = "black"),
        axis.text.y = element_text(angle = 0, size = 14, face = "bold", colour = "black"),
        axis.title=element_text(size=18, face = "bold", colour = "black"),
        strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
        strip.background = element_rect(fill="white"),
        legend.title=element_blank())


v1 <- unique(synDNAMetaAgM$Pool)
v2 <- unique(synDNAMetaAgM$Strand)
v3 <- unique(synDNAMetaAgM$ID)
v1v2v3 <- expand.grid(v1,v2,v3) #identify the set of all combination of pools (1:3), read direction (forward/verse), and sample ID

combinations <- v1v2v3[!duplicated(v1v2v3),] #verifies there are no duplications, which there weren't
sel <- c(intersect(grep("pool1", combinations$Var1), grep("pool1", combinations$Var3)), #identifies the set where the pool variable matches the sample ID, which makes sense because the pool variable is included in sample ID. Will reduce the set of combinations from 216 to 72
         intersect(grep("pool2", combinations$Var1), grep("pool2", combinations$Var3)),
         intersect(grep("pool3", combinations$Var1), grep("pool3", combinations$Var3)))

combinations <- combinations[sel,] #sorts the set of combinations sp that all the pool1's come first, then pool2, then pool3. 
colnames(combinations) <- c("Pool", "Strand", "ID")


coefall <- NULL

for(i in 1:dim(combinations)[1]){ #for each row in the set of combinations
  criteria <- as.vector(unlist(combinations[i,]))
  temp <- subset(synDNAMetaAgM, Pool %in% criteria & ID %in% criteria) #filters the dataset down to the forward and reverse reads for a given sample + spike-in combo. Because the sampleID also includes the directonality, in fact the model is calculated independently for the direction of the reads.
  temp <- temp[!(temp$Counts == 0), ] #remove the rows where a given synDNA had no reads that aligned to it
  lmtemp <- lm(temp$Dilution2~temp$CPMlog) #generate a linear model where stage of the serial dilution is related to the log(counts per million) of the reads aligned to a given synDNA/reads aligned to any synDNA
  coef <- c(criteria, lmtemp$coefficients[1], lmtemp$coefficients[2]) #records the intercept and slope respectively for the linear model
  coefall <- rbind(coefall, coef) #saves the intercept and slope for each set of parameters (pool, strand, and ID)
}

#at this point I can figure out the units of the equation: intercept: unitless. slope: log(counts/million)

colnames(coefall) <- c("Pool", "Strand", "ID", "a_intercept", "b_slope")
rownames(coefall) <- NULL

# Convert coefall to a data frame
coefall <- as.data.frame(coefall)

# Extract strand direction (Fwd/Rev) from the ID column
coefall$Strd <- gsub(".*_", "", coefall$ID)

# Remove rows where Strand is "reverse" but Strd is "Fwd"
coefall <- coefall[-which(coefall$Strand == "reverse" & coefall$Strd == "Fwd"), ]

# Remove rows where Strand is "forward" but Strd is "Rev"
coefall <- coefall[-which(coefall$Strand == "forward" & coefall$Strd == "Rev"), ]

write.table(coefall, "data/saliva_linear_models.tsv", quote = F, row.names = F, sep="\t")

