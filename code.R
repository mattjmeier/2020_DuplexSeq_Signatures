# Script to analyze distances between two dose groups of a chemical exposure using trinucleotide mutational patterns
library(ggplot2)
library(tidyverse)
library(MutationalPatterns)
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
otu = otu_table(mutation.matrix.transposed, taxa_are_rows = TRUE)
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
print(dist_methods)

# Loop over various distance methods, ordinate, plot
plist <- vector("list", length(dist_methods))
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(ps, method=i)
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
