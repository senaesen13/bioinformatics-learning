# Mini Analysis - Week 1
setwd("~/Desktop/bioinfo-learning/week1-day2-r-basics")
library(ggplot2)
library(dplyr)

# Load data
data <- read.csv("gene_expression.csv")

# Add color column
data$color <- "not significant"
data$color[data$significant == TRUE & data$log2FC > 1] <- "upregulated"
data$color[data$significant == TRUE & data$log2FC < -1] <- "downregulated"

# Volcano plot
ggplot(data, aes(x = log2FC, y = expression, color = color, label = gene)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  scale_color_manual(values = c("upregulated" = "red",
                                "downregulated" = "blue",
                                "not significant" = "gray")) +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "Expression") +
  theme_minimal()

# Top 3 genes
top3 <- data[data$significant == TRUE,]
top3 <- top3[order(-top3$expression),]
top3 <- head(top3, 3)
print(top3)
# Save the volcano plot as a PNG file
ggsave("volcano_plot.png", width = 8, height = 6)
log2(40/10)
