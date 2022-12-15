library(DESeq2)
library(tidyverse)
library(ggplot2)
library(here)
library(ggfortify)
library(ggbiplot)
library(lubridate)
library(utils)


# Data and metadata

# import count matrix
counts <- read.delim("data/raw_counts_matrix.tsv", row.names = 1)
head(counts, 10)

# import metadata
sampleinfo <- read_csv("data/sampleInfo.csv")
sampleinfo



# Quality control of the imported counts

# create DESeq object
dds <- DESeqDataSetFromMatrix(counts, 
                              colData = sampleinfo, 
                              design = ~Treated)
dds

# look at counts
counts(dds) %>% head(10)

# look at metadata
colData(dds)

# Visualize library size

# find library size
colSums(counts(dds))
apply(counts(dds), 2, sum)

# add library size to sampleinfo and create bar graph
sampleinfo$LibSize <- apply(counts(dds), 2, sum)
ggplot(data = sampleinfo,
       mapping = aes(x = Name, y = LibSize/1000000)) +
  geom_col(width = 0.85) +
  theme_bw() +
  labs(x = "Sample name", y = "Library size (million reads)")


# Visualize count distribution

boxplot(counts(dds))       # using boxplot() bc data not in long format required by ggplot

# perform variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = TRUE)     # blind the transformation to the experimental design

# check distribution using boxplot
boxplot(assay(vsd), 
        xlab = "", ylab = "Log2 counts per million", las = 2,
        main = "Normalised Distributions")
abline(h = median(assay(vsd)), col="blue")




# Principal Component Analysis (PCA) 

# Using plotPCA() 
plotPCA(vsd, intgroup = "condition")

# correct metadata
sampleinfo <- sampleinfo %>%
  mutate(condition = str_to_upper(condition))

# create data frame for PCA plot
PCA.plot <- plotPCA(vsd, intgroup = "condition", returnData = TRUE) %>%
  select(-group, -condition, Run = name) %>%
  mutate(Run = str_sub(Run, start = 2, end = -1)) %>%
  left_join(sampleinfo, by = "Run") %>%
  mutate(Replicate = as.factor(Replicate))

# recreate PCA plot
ggplot(data = PCA.plot, 
       mapping = aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "PC1: 61% variance", y = "PC2: 17% variance")

# examine batch effect
ggplot(data = PCA.plot, 
       mapping = aes(x = PC1, y = PC2, color = condition, shape = Replicate)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "PC1: 61% variance", y = "PC2: 17% variance")



# Using prcomp()
transformed.count.matrix <- assay(vsd)

# find top 500 genes with highest variance
row.sd <- apply(transformed.count.matrix, 1, sd)
row.sd.rank <- sort(row.sd, decreasing = TRUE)[1:500]
diff.gene <- names(row.sd.rank)

# create matrix for prcomp()
matrix.PCA <- t(transformed.count.matrix[diff.gene, ])

vsd.prcp <- prcomp(matrix.PCA)
summary(vsd.prcp)

autoplot(vsd.prcp)
ggbiplot::ggbiplot(vsd.prcp, scale = 0, var.axes = FALSE)