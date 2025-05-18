#Task: annotate the script for my own understanding. 
#Primary goal is to understand the two key equations for the absolute quantification.
#eq 1, for each bacterial genome: log(readWeight) = intercept + slope * normalizedReads #intercept: unitless. slope: log(counts/million). normalizedReads == log(CPM) where the CPM is readcount aligned to all MAGs ID'd as a given species with a coverage of 1x or higher/ readcount aligned to any MAG with a coverage of 1x or higher
#eq 2, for each bacterial genome: numCells = (readWeight * avogadro)/(genomeLength * molWeightBP) #Andy's question was: what are the units?

#numCells= readweight (countspermillion) * avogadro's number (units/mole) / the genome length (bp/cell) *1e6* the molecular weight of a based pair (g/bp * mol) 

#/mol / (g/mol*cell) = cell/g

log(readWeight) = intercept -log((ng/ul)) + slope (log(ng/ul) / log(cpm)) * normalized reads (log cpm)





rm(list = ls())

#load libraries:

library(tidyverse)
library(data.table)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(plyr)
library(dplyr) 
library(ggpubr)
library(Hmisc)
library(corrplot)
library(ggpmisc)
library(scales)
library(broom)
library(FactoMineR)
library(factoextra)

#Loading metadata, models, flow data and frequency table
setwd("/Users/aal88/Documents/GitHub/SynDNA")

linearmodels <- read.delim("data/saliva_linear_models.tsv", h=T) #this is the output of the linear models R script, recording the slope and interecept of the relationship between the aligned reads to a given synDNA/any synDNA and the dilution factor of the given synDNA for each pool, read direction, and sample
flowcounts <- read.delim("data/Saliva_Flow_data.tsv", head=T)
metadata <- read.delim("data/synDNA_metadata_updated.tsv", h=T)
freqtable <- fread("data/HMP_frequency_table.tsv", head=T)
freqtable <- as.data.frame(freqtable)

metadatagenomes <- fread("data/taxaID_length_fulllineage.tsv", sep = "\t", h=T)
metadatagenomes <- as.data.frame(metadatagenomes)


# Initialize output data frame
dfRawSpecies <- NULL

# Sum counts for each OTUID for each samples (pool * sample * read direction) in freqtable
dfsum <- ddply(freqtable, .(OTUID), colwise(sum))

# Loop over each sample column (starting from the 2nd column)
for(i in 2:ncol(dfsum)) {     #iterates over each sequencing sample (pool * sample * read direction)
  # Create a data frame with OTUID and raw counts for the current sample
  df <- data.frame(OTUID = dfsum$OTUID, RawCounts = dfsum[, i])   #create a dataframe that summarizes the number of reads aligned to each MAG in that sample
  SampleIDs <- colnames(dfsum) #currently unclear what this is for
  
  # Prepare metadata for merging
  metaspecies <- data.frame(     #summarizes for each MAG both its length and the species ID.
    OTUID = metadatagenomes$GenomeID,  
    Species = metadatagenomes$species,
    Length = metadatagenomes$Length
  )
  
  # Merge counts with metadata
  metaspecies <- merge(df, metaspecies, by = "OTUID") #combines the aligned reads for each MAG with the MAG length and species ID
  
  # Keep only rows with nonzero counts
  metaspecies <- metaspecies[which(metaspecies$RawCounts != 0), ] #remove MAGs that have no aligned reads (either not present or at very low abundance)
  
  # Aggregate raw counts by species
  dfaggreg <- aggregate(metaspecies$RawCounts, by = list(metaspecies$Species), FUN = sum) #seems to add together aligned reads from different OTUs that are ID'd as the same species
  colnames(dfaggreg) <- c("Species", "RawCounts")
  
  # Prepare output row for this sample 
  dfout <- data.frame(
    SampleID = SampleIDs[i],
    AlignedReads = sum(dfaggreg$RawCounts),  #adds up all of the reads that were aligned to any of the MAGs. 
    dfaggreg
  )
  
  # Append to the output data frame
  dfRawSpecies <- rbind(dfRawSpecies, dfout)
}

## Genome length (average) and coverage

# Initialize vectors to store genome lengths and total reads
len <- NULL
Tol <- NULL

# Loop through each row of dfRawSpecies
for(i in 1:nrow(dfRawSpecies)) { #for every sample/species pair
  # Get the mean genome length for the species
  len <- c(len, mean(metadatagenomes$Length[metadatagenomes$species == dfRawSpecies$Species[i]])) #the mean genome length is calculated for all MAGs belonging to the same species 
  # Get the total reads for the sample
  Tol <- c(Tol, metadata$TotalReads[metadata$ID == dfRawSpecies$SampleID[i]]) #pull the total number of reads (both aligned and unaligned) for each sample
}

# Add the calculated genome length and total reads to dfRawSpecies
dfRawSpecies$GenomeLength <- len
dfRawSpecies$TotalReads <- Tol

# Calculate coverage
dfRawSpecies$Cov <- (dfRawSpecies$RawCounts * 150) / dfRawSpecies$GenomeLength  #to calculate coverage, read counts are multiplied by read length (which is 150, presumably)

## Removing rare/low abundant species
dfCovSel <- dfRawSpecies[dfRawSpecies$Cov >= 1, ] #removes species with coverage of less than 1x, which cuts a ~23000 row dataframe down to ~500 rows

## Applying the linear models
dfall <- NULL

i <- 1
for(i in 1:nrow(linearmodels)) { #for each sample/read direction/synDNA combination:
  # Select data for the current sample
  df <- dfCovSel[dfCovSel$SampleID %in% linearmodels$ID[i], ]
  
  # Calculate relative abundance (%)
  df$Relab <- (df$RawCounts / sum(df$RawCounts)) * 100 #this is relative abundance calculated relative to all the other high abundance species? It seems it would be better calculated as relative to all other species-level MAGs
  
  # Normalize counts (log10 CPM)
  df$Norm <- log10((df$RawCounts / sum(df$RawCounts)) * 1e6)  # this is parallel to the synDNA log(CPM): log (aligned to given synDNA / aligned to any synDNA)
  
  # Apply linear model to estimate weight
  reads2weigth <- linearmodels$a_intercept[i] + (linearmodels$b_slope[i] * df$Norm)
  reads2weigth <- 10^(-reads2weigth)  # Convert from log scale
  
  # Calculate cell number
  GenomeSize <- df$GenomeLength
  #df$CellNumber <- (reads2weigth * 6.022e23) / (GenomeSize * 1e9* 650) #estimated cell number is the readweight (unitless) * avogadro's number (units/mole) / the genome length (bp) * the molecular weight of a based pair (g/bp * mol) * a billion for some unknown reason
  #my update:
  df$CellNumber <- (reads2weigth * 6.022e23) / (GenomeSize *1e6* 650) #estimated cell number is the readweight * avogadro's number/ the genome length * the molecular weight of a based pair
  
  # Predicted relative abundance (%)
  df$Predicted <- (df$CellNumber / sum(df$CellNumber)) * 100
  
  # Order by predicted abundance
  df <- df[order(df$Predicted, decreasing = TRUE), ]
  
  # Add flow cytometry cell counts
  df$FC_cells_per_ul_r1 <- flowcounts$FC_cells_per_ul_r1[which(flowcounts$SampleID == gsub("_.*", "", df$SampleID)[1])]
  
  # Calculate cell flow and percentage flow
  #df$CellFlow <- (df$RawCounts / df$GenomeLength) / df$FC_cells_per_ul_r1 #this is supposed to be the absolute number of cells I think. But maybe it should be by coverage instead of just raw reads/length?
  df$CellFlow <- df$Cov * df$FC_cells_per_ul_r1 * 1000  #need to multiply by 1000 to get cells/gram 
  
  df$PercFlow <- (df$CellFlow / sum(df$CellFlow)) * 100
  
  # Combine results
  dfall <- rbind(dfall, df)
}


## Select strand and pool

# Merge metadata with results
dfallmeta <- merge(metadata, dfall, by.x = "ID", by.y = "SampleID")

# Keep only rows with 'forward' strand
dfallmeta <- dfallmeta[grep("forward", dfallmeta$Strand), ] #why only forward?

# Keep only rows with 'pool2'
#dfallmeta <- dfallmeta[grep("pool2", dfallmeta$Pool), ] #why only pool2?

# Remove rows for Prevotella melaninogenica
dfallmeta <- dfallmeta[-grep("Prevotella melaninogenica", dfallmeta$Species), ] #why removed?


## Correlation analysis: SynDNA Predicted vs. Flow Cytometer

# Using TidyVerse to calculate Pearson Correlation between methods
# Predicted vs. Flow cytometer

library(tidyverse)
library(broom)

analysis <- dfallmeta %>%
  group_by(Sample, Pool) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(CellNumber ~ CellFlow, data = .)),
    cor = map(data, ~tidy(cor.test(.x$CellNumber, .x$CellFlow, method = "pearson"), 3))
  )

stats <- analysis %>%
  unnest(cor)

## Plotting correlation graphs

# Set plot dimensions (if using RStudio or Jupyter)
options(repr.plot.width = 15, repr.plot.height = 3)

# Define sample order for plotting
conditiongroup <- c("A1", "C1", "D1", "E1", "F1", "H1")
dfallmeta$Sample <- factor(dfallmeta$Sample, levels = conditiongroup, ordered = TRUE)
dfallmeta$Pool <- as.factor(dfallmeta$Pool)


#my version of the figure:
# Find the min and max across both axes
min_val <- min(c(dfallmeta$CellNumber, dfallmeta$CellFlow), na.rm = TRUE)
max_val <- max(c(dfallmeta$CellNumber, dfallmeta$CellFlow), na.rm = TRUE)

dfallmeta$synDNAToFlow <- dfallmeta$Predicted * dfallmeta$FC_cells_per_ul_r1 * 1000

analysis <- dfallmeta %>%
  group_by(Sample, Pool) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(synDNAToFlow ~ CellFlow, data = .)),
    cor = map(data, ~tidy(cor.test(.x$synDNAToFlow, .x$CellFlow, method = "pearson"), 3))
  )

stats <- analysis %>%
  unnest(cor)

p1 <- ggplot(dfallmeta, aes(synDNAToFlow/1000, CellFlow/1000)) +
  geom_point(shape = 20, size = 4) +
  facet_wrap(~ Sample, ncol = 6) +
  geom_smooth(method = "lm", color = "black", se = TRUE, formula = y ~ x) +
  #coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Predicted cell number according to synDNA", y = "Flow Cytometer-estimated cell numbers")+
  theme_bw()


p1 <- ggplot(dfallmeta, aes(synDNAToFlow/1000, CellFlow/1000)) +
  geom_point(shape = 20, size = 4) +
  facet_grid(Pool ~ Sample) +  # Rows = Pool, Columns = Sample
  geom_smooth(method = "lm", color = "black", se = TRUE, formula = y ~ x) +
  #coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Predicted cell number according to synDNA",
    y = "Flow Cytometer-estimated cell numbers"
  ) +
  theme_bw()

  
#check which species are a closer match:
dfallmeta$countMatch <- dfallmeta$synDNAToFlow/dfallmeta$CellFlow




# Create the plot
p1 <- ggplot(dfallmeta, aes(CellNumber, CellFlow)) +
  geom_point(shape = 20, size = 4) +
  facet_wrap(~ Sample, ncol = 6) +
  geom_smooth(method = "lm", color = "black", se = TRUE, formula = y ~ x) + 
  labs(x = "Predicted cell number according to synDNA", y = "Flow Cytometer-estimated cell numbers") #+ 
  geom_text(data = stats, aes(label = sprintf("r = %s", round(estimate, 3)), x = 9, y = 43), size = 5) +
  geom_text(data = stats, aes(label = sprintf("p = %s", round(p.value, 3)), x = 9, y = 37), size = 5) +
  +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_bw() + guides(fill = guide_legend(ncol = 1)) #+
  theme(
    axis.text.x = element_text(angle = 0, size = 16, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 16, colour = "black"),
    axis.text = element_text(size = 16, colour = "black"),
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
    strip.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "bottom", legend.box = "vertical", legend.margin = margin(),
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )

# To display the plot (uncomment if running interactively)
# print(p1)

  p1 <- ggplot(dfallmeta, aes(Predicted, PercFlow)) +
    facet_grid(Pool ~ Sample) +  # Rows = Pool, Columns = Sample
    geom_point(shape = 20, size = 4) +
    geom_smooth(method = "lm", color = "black", se = TRUE, formula = y ~ x)
