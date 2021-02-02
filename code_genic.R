# Script to analyze distances between two dose groups of a chemical exposure using trinucleotide mutational patterns
library(ggplot2)
library(tidyverse)
library(spgs)
library(MutationalPatterns)
#library(SomaticSignatures)
library(deconstructSigs)
library(phyloseq)
library(vegan)
library(IRanges)
library(GenomicRanges)
library(fuzzyjoin)
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

genic_regions <- read.table("genic_regions.txt", header=T, sep="\t")
sampledata <- read.table("./sampledata.txt", sep="\t", header=T, row.names=1)
sampledata$Dose <- factor(sampledata$Dose)
sampledata$sample <- row.names(sampledata)
#per.nucleotide.data.only_mutations <- read.table("./VCFInformationMutationData_with_lacZ.txt", header=T, sep="\t") # Doesn't include all sites!
bgzip -r 
per.nucleotide.data.all <- read.table(gzfile("./prj00068-20201202/mouse-processed-genome.mut.gz"), header=T, sep="\t")
per.nucleotide.data <- dplyr::left_join(per.nucleotide.data.all, sampledata)

per.nucleotide.data <- fuzzyjoin::genome_join(per.nucleotide.data, genic_regions, by=c("contig", "start", "end"), mode="inner")
colnames(per.nucleotide.data)[1:3] <- c("contig","start","end")

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

#################################################################
# From this point on, we split the data into genic vs. intergenic
#################################################################

trinucleotide_frequencies_intergenic <- trinucleotide_frequencies_depth %>%
  filter(location_relative_to_genes=="intergenic")

trinucleotide_frequencies_genic <- trinucleotide_frequencies_depth %>%
  filter(location_relative_to_genes=="genic")

# Convert table above to wide format
trinucleotide_frequencies_wide <- trinucleotide_frequencies_genic %>%
  select(context_with_mutation, frequency, sample) %>%
  distinct() %>%
  pivot_wider(names_from = context_with_mutation, values_from = frequency)

# Revert to dataframe
trinucleotide_frequencies_wide <- as.data.frame(trinucleotide_frequencies_wide)
# Convert sample names to row names
row.names(trinucleotide_frequencies_wide) <- trinucleotide_frequencies_wide$sample
# Remove sample columns
trinucleotide_frequencies_wide <- trinucleotide_frequencies_wide %>% select(-sample)
# Reorder columns according to COSMIC

actual_colnames <- colnames(signatures.cosmic)[colnames(signatures.cosmic) %in% colnames(trinucleotide_frequencies_wide)]
empty_matrix <- matrix(data = 0, nrow=nrow(trinucleotide_frequencies_wide), ncol=96)
colnames(empty_matrix) <- colnames(signatures.cosmic)
rownames(empty_matrix) <- rownames(trinucleotide_frequencies_wide)
for(i in actual_colnames) {
  print(i)
  empty_matrix[,i] <- trinucleotide_frequencies_wide[,i]
  }
trinucleotide_frequencies_mut_matrix <- empty_matrix
# Remove extra sample
trinucleotide_frequencies_mut_matrix <- trinucleotide_frequencies_mut_matrix[1:24,]

# Set NA to 0
trinucleotide_frequencies_mut_matrix[is.na(trinucleotide_frequencies_mut_matrix)] <- 0

# Convert to proportion of row sums
trinucleotide_frequencies_proportions <- trinucleotide_frequencies_mut_matrix/rowSums(trinucleotide_frequencies_mut_matrix)

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

#sorted <- dendextend::sort_dist_mat(pairwise_cosine_similarity_sigs) # Doesn't look better
sorted <- dendsort::dendsort(hclust(as.dist(pairwise_cosine_similarity_sigs), method = "average"))
pheatmap::pheatmap(t(pairwise_cosine_similarity_sigs),
                   cluster_cols = sorted,
                   cluster_rows = sorted,
                   annotation=sampledata %>% select(Dose))

pheatmap::pheatmap(t(pairwise_cosine_similarity_sigs),
                   cluster_cols=F,
                   annotation=sampledata %>% select(Dose))

################################################################
# Shrunken Centroid Analysis
################################################################
library(pamr)
library(heatmap3)

################################################################
# Training Set
################################################################
allDat <- as.matrix(t(trinucleotide_frequencies_proportions))
trainSet <- allDat[,sampledata$Group %in% c("C", "H")]
labels <- sampledata$Group[sampledata$Group %in% c("C", "H")]
p <- rownames(allDat)

################################################################
# Shrinkage Estimation
################################################################
train.data <- list(x = as.matrix(trainSet), y = as.factor(labels),
                   geneid=as.character(p), genenames=as.character(p))
my.train <- pamr.train(train.data)

################################################################
# 6-fold Cross Validation (CV) (6 fold because min class size is 6)
################################################################
results <- pamr.cv(my.train, train.data)
results

#Call:

#pamr.cv(fit = my.train, data = train.data)
#   threshold nonzero errors
#1  0.000     96      1    
#2  0.146     75      1    
#3  0.292     61      1    
#4  0.438     50      1    
#5  0.585     38      1    
#6  0.731     30      1    
#7  0.877     26      1    
#8  1.023     23      1     
#9  1.169     23      1    
#10 1.315     21      1    
#11 1.461     15      1    
#12 1.607     13      1    
#13 1.754     10      1    
#14 1.900      9      1    
#15 2.046      8      1    
#16 2.192      6      1    
#17 2.338      5      1    
#18 2.484      5      1    
#19 2.630      5      1    
#20 2.777      5      1    
#21 2.923      5      1    
#22 3.069      3      1    
#23 3.215      3      1    
#24 3.361      3      1    
#25 3.507      3      1    
#26 3.653      3      1    
#27 3.799      3      2    
#28 3.946      2      4    
#29 4.092      1      4    
#30 4.238      0      4

################################################################
# Choosing the largest delta with the lowest number of errors based on the CV
################################################################
# 3 mutations
cut <- 3.653
pamr.confusion(my.train, threshold = cut)

############################################################
# Writing out the centroids
################################################################
gs <- pamr.listgenes(my.train, train.data,  threshold=cut)[,1]
flag <- p %in% gs
x <- trainSet[flag,]
pp <- p[flag]
td <- list(x = as.matrix(x), y = as.factor(labels), geneid=as.character(pp), genenames=as.character(pp))
out <- pamr.train(td)
classifier <- data.frame(ID = pp)
classifier <- cbind.data.frame(classifier, out$centroids)
names(classifier) <- c("ID", "Control.score", "BaP.score")
classifier$std.dev <- out$sd

write.table(classifier,
            file="BaP Shunken Centroid Signature.txt", row.name = FALSE, col.name = TRUE, quote = FALSE, sep="\t")

################################################################
# Let's have a look at where the L and M groups get classified
################################################################
testDat <- cbind.data.frame(rownames(allDat), allDat)
names(testDat)[1] <- "ID"
testDat <- merge(classifier, testDat)

################################################################
#Predictions
################################################################
probBaP <- rep(0, ncol(testDat)-4)

for(k in 1:length(probBaP)){
  #Gausian Linear Discriminant Analysis
  BaP <- sum(((testDat[,k+4] - testDat$BaP.score)/testDat$std.dev)^2) - 2*log(0.5)
  Control <- sum(((testDat[,k+4] - testDat$Control.score)/testDat$std.dev)^2) - 2*log(0.5)
  probBaP[k] <- 1/(1 + exp(BaP/2 - Control/2))
}             

################################################################
#Heatmap
################################################################
heatDat <- as.matrix(testDat[,-c(1:4)])
rownames(heatDat) <- as.character(testDat$ID)
colnames(heatDat) <- as.character(sampledata$Group)
reOrder <- c(grep("C", sampledata$Group), grep("H", sampledata$Group),
             grep("L", sampledata$Group), grep("M", sampledata$Group))
heatDat <- heatDat[,reOrder]
probBaP <- probBaP[reOrder]

##############################################################
#Colours
my.col <- matrix(rep("white", 2*ncol(heatDat)), ncol = 2)
my.col[grep("C", colnames(heatDat)),1] <- "blue"
my.col[grep("H", colnames(heatDat)),1] <- "red"

flag <- probBaP > 0.9
my.col[flag,2] <- "red"
flag <- probBaP < 0.1
my.col[flag,2] <- "blue"

##############################################################
#Row center based on the mean of the controls
mn <- apply(heatDat[,1:6], 1, mean)
heatDat <- heatDat - mn

end.breaks <- c(0, 12, ncol(heatDat))
y <- NULL
my.col2 <- NULL
g <- NULL
col.blank <- rep(NA, nrow(heatDat))
pp <- NULL

for(k in 2:length(end.breaks)){
  y <- cbind(y, heatDat[, (end.breaks[k-1]+1):end.breaks[k]], col.blank)
  my.col2 <- rbind(my.col2, my.col[(end.breaks[k-1]+1):end.breaks[k],], c("white", "white"))
  g <- c(g, colnames(heatDat)[(end.breaks[k-1]+1):end.breaks[k]], " ")
}

y <- y[,-ncol(y)]
my.col2 <- my.col2[-nrow(my.col2),]
g <- g[-length(g)]
colnames(y) <- g

colnames(my.col) <- c("Class", "Prediction")
colnames(my.col2) <- c("Class", "Prediction")

heatmap3(y, Rowv = NA, Colv = NA, scale = "none", margins = c(10, 10), main = "",
         cexRow=1, cexCol=1, ColSideColors = my.col2)
