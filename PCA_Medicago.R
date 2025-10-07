# Load libraries
library(SNPRelate)
library(ggplot2)
library(RColorBrewer)

# Input / Output
vcf.fn <- "Joint_Medicago.vcf.gz"
gds.fn <- "Medicago.gds"

# Convert VCF to GDS
snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only")

# Load GDS file
genofile <- snpgdsOpen(gds.fn)

# Run PCA
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=FALSE)

# % variance explained
pc.percent <- round(pca$varprop*100, 2)

# Make dataframe for plotting
pca.df <- data.frame(
  sample.id = pca$sample.id,
  EV1 = pca$eigenvect[,1],
  EV2 = pca$eigenvect[,2],
  EV3 = pca$eigenvect[,3]
)

# 2D PCA plot
ggplot(pca.df, aes(x=EV1, y=EV2, label=sample.id)) +
  geom_point(size=5, shape=21, fill="steelblue", colour="black") +
  geom_text(vjust=1.5, hjust=1.2, size=3) +
  xlab(paste0("PC1 (", pc.percent[1], "%)")) +
  ylab(paste0("PC2 (", pc.percent[2], "%)")) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

ggsave("Medicago_PCA_2D.png", width=6, height=5)

# Optional: 3D scatterplot (if you want)
# Uncomment if you install scatterplot3d
# library(scatterplot3d)
# scatterplot3d(x=pca.df$EV1, y=pca.df$EV3, z=pca.df$EV2,
#               pch=21, bg="steelblue", cex.symbols=1.5,
#               xlab=paste0("PC1 (", pc.percent[1], "%)"),
#               ylab=paste0("PC3 (", pc.percent[3], "%)"),
#               zlab=paste0("PC2 (", pc.percent[2], "%)"))

