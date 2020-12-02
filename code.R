# Script to analyze distances between two dose groups of a chemical exposure using trinucleotide mutational patterns
library(ggplot2)
library(tidyverse)
library(spgs)
library(MutationalPatterns)
#library(SomaticSignatures)
library(deconstructSigs)
library(phyloseq)
library(vegan)
set.seed(3135)

# Import trinucleotide frequencies from spreadsheet saved as TSV
mutation.df <- read.table("./BaPDSTrinucleotideFrequenciesBaP.txt", header=T, sep="\t")

# Clean column names
colnames(mutation.df) <- c("mutation_type", "0", "12.5", "25", "50")

# Remove final row showing totals
mutation.df <- mutation.df[1:96,]

# Make data into long format
mutation.df.long <- tidyr::pivot_longer(mutation.df, cols=2:5, names_to="dose", values_to="frequency")

# Convert things to matrix format
mutation.matrix <- as.matrix(mutation.df[,2:5])
row.names(mutation.matrix) <- mutation.df[,1]
mutation.matrix.transposed <- t(mutation.matrix)

# Vegan package
dist <- vegdist(mutation.matrix, method="euclidean") # Could change distance here to others...
mds <- metaMDS(mutation.matrix, k=2)
stressplot(mds)
ordiplot(mds, type="n")
orditorp(mds, display="species",col="red",air=0.01)

dist_transposed <- vegdist(mutation.matrix.transposed, method="euclidean") # Could change distance here to others...
mds_transposed <- metaMDS(mutation.matrix.transposed, k=2) # Stress is nearly zero
stressplot(mds_transposed)
ordiplot(mds_transposed, type="n")
orditorp(mds_transposed, display="sites",cex=1.25,air=0.01)
orditorp(mds_transposed, display="species",col="red",air=0.01)

# Do more distance metrics in phyloseq

# Create phyloseq object
otu = otu_table(mutation.matrix.transposed, taxa_are_rows = F)
tax = as.matrix(mutation.df[,1])
rownames(tax) <- tax
tax = tax_table(tax)
sampledata <- data.frame(dose=factor(c("0","12.5","25","50")), chemical="BaP") # Need a dummy column of chemical...
row.names(sampledata) <- c("0","12.5","25","50")
sampledata_ps <- sample_data(sampledata)
# Put it together
ps <- phyloseq(otu, tax, sampledata_ps)

# Get distances metrics available
dist_methods <- unlist(distanceMethodList)
dist_methods <- dist_methods[-(1:3)]
dist_methods <- dist_methods[!dist_methods %in% c("cao","ANY")]
print(dist_methods)

# Loop over various distance methods, ordinate, plot
plist <- vector("list", length(dist_methods))
for( i in dist_methods ){
  message(i)
  # Calculate distance matrix
  iDist <- phyloseq::distance(ps, method=i)
  # Calculate ordination
  iMDS  <- ordinate(ps, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(ps, iMDS, type="samples", color="dose")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the plot
  plist[[i]] = p
}

# Create dataframe from ordination plot data
df = plyr::ldply(plist, function(x) x$data)
names(df)[1] <- "distance"

# Plot as faceted ggplot
ggplot(df, aes(Axis.1, Axis.2, color=dose)) +
  geom_point(size=3, alpha=0.5) +
  facet_wrap(~distance, scales="free") +
  ggtitle("MDS on various distance metrics for Duplex Sequencing BaP Bone Marrow dataset")

# Manhattan distance split CCA plot
plot_ordination(ps,
                ordinate(ps, method="CCA", distance="manhattan"),
                type="split",
                color="dose",
                title="Manhattan distance, split CCA plot") 

# Manhattan distance, but with different ordination methods
dist = "manhattan"
ord_meths = c("DCA", "CCA", "RDA",  "NMDS", "MDS", "PCoA")
plist = plyr::llply(as.list(ord_meths), function(i, ps, dist){
  ordi = ordinate(ps, method=i, distance=dist)
  plot_ordination(ps, ordi, "samples", color="dose")
}, ps, dist)

names(plist) <- ord_meths

pdataframe = plyr::ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})

names(pdataframe)[1] = "method"

# Plot for manhattant distance but various ordination methods
ggplot(pdataframe, aes(Axis_1, Axis_2, color=dose)) +
  geom_point(size=4) +
  facet_wrap(~method, scales="free") +
  scale_colour_brewer(palette="Set1")

# Cosine similarity
pairwise_cosine_similarity <- MutationalPatterns::cos_sim_matrix(mutation.matrix, mutation.matrix)
MutationalPatterns::plot_cosine_heatmap(pairwise_cosine_similarity)

#################################################
# Per mouse data ################################
#################################################

sampledata <- read.table("./sampledata.txt", sep="\t", header=T, row.names=1)
sampledata$Dose <- factor(sampledata$Dose)
sampledata$sample <- row.names(sampledata)
#per.nucleotide.data.only_mutations <- read.table("./VCFInformationMutationData_with_lacZ.txt", header=T, sep="\t") # Doesn't include all sites!
per.nucleotide.data.all <- read.table("./prj00068-20201202/mouse-processed-genome.mut", header=T, sep="\t")
per.nucleotide.data <- dplyr::left_join(per.nucleotide.data.all, sampledata)

# trinucleotide_frequencies_global <- per.nucleotide.data %>%
#   group_by(context, subtype) %>%
#   summarise(adjusted_depth_per_trinucleotide = sum(adjusted_depth))
# 
# trinucleotide_frequencies_dose_group <- per.nucleotide.data %>%
#   group_by(context, subtype, Dose) %>%
#   summarise(adjusted_depth_per_trinucleotide = sum(adjusted_depth))
# 
# trinucleotide_frequencies_per_mouse <- per.nucleotide.data %>%
#   group_by(context, subtype, sample) %>%
#   summarise(adjusted_depth_per_trinucleotide = sum(adjusted_depth))
# 
# trinucleotide_counts_per_dose <- per.nucleotide.data %>%
#   group_by(context, subtype, Dose) %>%
#   summarise(count = sum(final_somatic_alt_depth))

################
# Clean up data:
# Select only SNVs
# Remove sites where final somatic alt depth is zero
# Get reverse complement of sequence context where mutation is listed on purine context
# Change all purine substitutions to pyrimidine substitutions
# Make new column with COSMIC-style 96 base context
# Calculate depth for each sequence context and dose group
# Calculate frequency for each mouse within each 96 trinucleotide mutation
trinucleotide_frequencies_depth <- per.nucleotide.data %>%
  mutate(context = ifelse(subtype %in% c("G>T","G>A","G>C","A>T","A>C","A>G"),
                          mapply(function(x) spgs::reverseComplement(x, case="upper"), context),
                          context)) %>%
  mutate(subtype = str_replace(subtype, "G>T", "C>A")) %>%
  mutate(subtype = str_replace(subtype, "G>A", "C>T")) %>%
  mutate(subtype = str_replace(subtype, "G>C", "C>G")) %>%
  mutate(subtype = str_replace(subtype, "A>T", "T>A")) %>%
  mutate(subtype = str_replace(subtype, "A>C", "T>G")) %>%
  mutate(subtype = str_replace(subtype, "A>G", "T>C")) %>%
  mutate(context_with_mutation = paste0(str_sub(context, 1, 1),"[",subtype,"]",str_sub(context, 3, 3)) ) %>%
  group_by(context, Dose) %>%
  mutate(group_depth = sum(informative_total_depth)) %>%
  ungroup() %>%
  group_by(context_with_mutation, sample) %>%
  mutate(frequency = (sum(final_somatic_alt_depth)/group_depth) ) %>%
  ungroup() %>%
  filter(variation_type=="snv") %>%
  filter(!final_somatic_alt_depth==0)

# Convert table above to wide format
trinucleotide_frequencies_wide <- trinucleotide_frequencies_depth %>%
  select(context_with_mutation, frequency, sample) %>%
  distinct() %>%
  pivot_wider(names_from = context_with_mutation, values_from = frequency)

# Set NA to 0
trinucleotide_frequencies_wide[is.na(trinucleotide_frequencies_wide)] <- 0
# Revert to dataframe
trinucleotide_frequencies_wide <- as.data.frame(trinucleotide_frequencies_wide)
# Convert sample names to row names
row.names(trinucleotide_frequencies_wide) <- trinucleotide_frequencies_wide$sample
# Remove sample columns
trinucleotide_frequencies_wide <- trinucleotide_frequencies_wide %>% select(-sample)
# Reorder columns according to COSMIC
trinucleotide_frequencies_wide <- trinucleotide_frequencies_wide[colnames(signatures.cosmic)]
# Remove extra sample
trinucleotide_frequencies_wide <- trinucleotide_frequencies_wide[1:24,]

# Convert to proportion of row sums
trinucleotide_frequencies_proportions <- trinucleotide_frequencies_wide/rowSums(trinucleotide_frequencies_wide)

# No idea where their read depth is coming from.

#################################
# Per mouse re-analysis as above#
#################################

# Vegan package
dist <- vegdist(trinucleotide_frequencies_proportions, method="euclidean") # Could change distance here to others...
mds <- metaMDS(trinucleotide_frequencies_proportions, k=2)
stressplot(mds)
ordiplot(mds, type="n")
orditorp(mds, display="sites",col="red",air=0.01)
orditorp(mds, display="species",col="black",air=0.01)

trinucleotide_frequencies_proportions.transposed <- t(trinucleotide_frequencies_proportions)
dist_transposed <- vegdist(trinucleotide_frequencies_proportions.transposed, method="euclidean") # Could change distance here to others...
mds_transposed <- metaMDS(trinucleotide_frequencies_proportions.transposed, k=2) # Stress is nearly zero
stressplot(mds_transposed)
ordiplot(mds_transposed, type="n")
orditorp(mds_transposed, display="sites",cex=1.25,air=0.01)
orditorp(mds_transposed, display="species",col="red",air=0.01)

# Do more distance metrics in phyloseq

# Create phyloseq object
otu = otu_table(trinucleotide_frequencies_proportions, taxa_are_rows = F)
tax = as.matrix(colnames(trinucleotide_frequencies_proportions))
rownames(tax) <- tax
tax = tax_table(tax)
sampledata_ps <- sample_data(sampledata)
# Put it together
ps <- phyloseq(otu, tax, sampledata_ps)

# Get distances metrics available
dist_methods <- unlist(distanceMethodList)
dist_methods <- dist_methods[-(1:3)]
dist_methods <- dist_methods[!dist_methods %in% c("ANY")]
print(dist_methods)

# Loop over various distance methods, ordinate, plot
plist <- vector("list", length(dist_methods))
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(ps, method=i)
  # Calculate ordination
  iMDS  <- ordinate(ps, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(ps, iMDS, type="samples", color="Dose")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the plot
  plist[[i]] = p
}

# Create dataframe from ordination plot data
df = plyr::ldply(plist, function(x) x$data)
names(df)[1] <- "distance"

# Plot as faceted ggplot
ggplot(df, aes(Axis.1, Axis.2, color=Dose)) +
  geom_point(size=3, alpha=0.5) +
  facet_wrap(~distance, scales="free") +
  ggtitle("MDS on various distance metrics for Duplex Sequencing BaP Bone Marrow dataset")

# Binomial distance split NMDS plot
plot_ordination(ps,
                ordinate(ps, method="NMDS", distance="binomial"),
                type="biplot",
                color="Dose",
                title="Binomial distance, NMDS biplot") 

# Manhattan distance, but with different ordination methods
dist = "binomial"
ord_meths = c("DCA", "CCA", "RDA",  "NMDS", "MDS", "PCoA")
plist = plyr::llply(as.list(ord_meths), function(i, ps, dist){
  ordi = ordinate(ps, method=i, distance=dist)
  plot_ordination(ps, ordi, "samples", color="Dose")
}, ps, dist)

names(plist) <- ord_meths

pdataframe = plyr::ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})

names(pdataframe)[1] = "method"

# Plot for binomial distance but various ordination methods
ggplot(pdataframe, aes(Axis_1, Axis_2, color=Dose)) +
  geom_point(size=4) +
  facet_wrap(~method, scales="free") +
  scale_colour_brewer(palette="Set1")

# Cosine similarity
pairwise_cosine_similarity <- MutationalPatterns::cos_sim_matrix(t(trinucleotide_frequencies_proportions),
                                                                 t(trinucleotide_frequencies_proportions))

pheatmap::pheatmap(t(pairwise_cosine_similarity),
                   cluster_rows=F,
                   annotation=sampledata %>% select(Dose))

pheatmap::pheatmap(t(pairwise_cosine_similarity),
                   cluster_rows=F,
                   cluster_cols=F,
                   annotation=sampledata %>% select(Dose))

# Cosine similarity between known signatures
signatures <- get_known_signatures()
pairwise_cosine_similarity_sigs <- MutationalPatterns::cos_sim_matrix(t(trinucleotide_frequencies_proportions),
                                                                 signatures)

pheatmap::pheatmap(t(pairwise_cosine_similarity_sigs),
                   annotation=sampledata %>% select(Dose))

pheatmap::pheatmap(t(pairwise_cosine_similarity_sigs),
                   cluster_cols=F,
                   annotation=sampledata %>% select(Dose))

