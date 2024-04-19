setwd("/work/larylab/NAYEMA/heatmap/")
df<-read.csv("snp-gene-heatmap.csv", header=TRUE)
head(df)


library(ggplot2)
library(dplyr)



# Convert P-Value to numeric
df$P.Value <- as.numeric(df$P.Value)


library(viridis)

heatmap<-ggplot(df, aes(x = Gene.Symbol, y = Tissue, fill = -log10(P.Value))) +
  geom_tile() +
  scale_fill_viridis() +  # Using Viridis color palette
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid
        panel.background = element_rect(fill = "white"),
        legend.title = element_text(face = "bold", size = 14),  # Bold legend title
        legend.text = element_text(face = "bold", size = 12),  # Increase legend text size
        axis.title = element_text(face = "bold", size = 14),  # Bold axis titles
        axis.text = element_text(face = "bold", size = 12)) +  # Increase axis label text size
  labs(title = "Heatmap of SNP-gene associations based on -log10(P-Value)",
       x = "Gene Symbol", y = "Tissue", fill = "-log10(P-Value)")  # Specify titles and labels
ggsave("result/heatmap-rs4690325.png", plot = heatmap, width = 10, height = 8, dpi = 300)

#####

df2<-read.csv("snp-gene-enrich.csv", header=TRUE)
head(df2)
# Convert P-Value to numeric
df2$P.Value <- as.numeric(df2$P.Value)

heatmap<-ggplot(df2, aes(x = Gene.Symbol, y = Tissue, fill = -log10(P.Value))) +
  geom_tile() +
  scale_fill_viridis() +  # Using Viridis color palette
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid
        panel.background = element_rect(fill = "white"),
        legend.title = element_text(face = "bold", size = 14),  # Bold legend title
        legend.text = element_text(face = "bold", size = 12),  # Increase legend text size
        axis.title = element_text(face = "bold", size = 14),  # Bold axis titles
        axis.text = element_text(face = "bold", size = 12)) +  # Increase axis label text size
  labs(title = "Heatmap of SNP-gene associations based on -log10(P-Value)",
       x = "Gene Symbol", y = "Tissue", fill = "-log10(P-Value)")  # Specify titles and labels
ggsave("result/heatmap-rs35220088.png", plot = heatmap, width = 10, height = 8, dpi = 300)


#########
#####

df2<-read.csv("DATA3.csv", header=TRUE)
head(df2)
# Convert P-Value to numeric
df2$P.Value <- as.numeric(df2$P.Value)

heatmap<-ggplot(df2, aes(x = Gene.Symbol, y = Tissue, fill = -log10(P.Value))) +
  geom_tile() +
  scale_fill_viridis() +  # Using Viridis color palette
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid
        panel.background = element_rect(fill = "white"),
        legend.title = element_text(face = "bold", size = 14),  # Bold legend title
        legend.text = element_text(face = "bold", size = 12),  # Increase legend text size
        axis.title = element_text(face = "bold", size = 14),  # Bold axis titles
        axis.text = element_text(face = "bold", size = 12)) +  # Increase axis label text size
  labs(title = "Heatmap of SNP-gene associations based on -log10(P-Value)",
       x = "Gene Symbol", y = "Tissue", fill = "-log10(P-Value)")  # Specify titles and labels
ggsave("result/heatmap-rs2292626.png", plot = heatmap, width = 10, height = 8, dpi = 300)



