# ============================================
# Week 1 - Day 2
# Topic: Data Frames, CSV files, and Plots
# ============================================

# Set working directory - tell R where to find files
setwd("~/Desktop/bioinfo-learning/week1-day2-r-basics")

# Load ggplot2 library
library(ggplot2)

# ============================================
# PART 1 - Create a data frame manually
# ============================================

# data.frame() creates a table with rows and columns
# like an Excel spreadsheet in R
gene_table <- data.frame(
  gene = c("BRCA1", "TP53", "EGFR", "MYC", "PTEN"),  # text column
  expression = c(10, 25, 3, 47, 8),                    # number column
  significant = c(TRUE, TRUE, FALSE, TRUE, FALSE)       # TRUE/FALSE column
)

# Look at the table
gene_table

# How many rows and columns?
nrow(gene_table)   # 5 rows = 5 genes
ncol(gene_table)   # 3 columns

# Structure - shows data types of each column
# chr = text, num = numbers, logi = TRUE/FALSE
str(gene_table)

# Access one column with $
gene_table$expression

# Calculate mean of expression column
mean(gene_table$expression)

# Filter - show only significant genes
gene_table[gene_table$significant == TRUE, ]

# ============================================
# PART 2 - Read a CSV file into R
# ============================================

# read.csv() loads a CSV file as a data frame
# This is how you load real RNA-seq data
data <- read.csv("gene_expression.csv")

# Explore the data
data          # print the whole table
str(data)     # show structure
head(data)    # show first 6 rows

# ============================================
# PART 3 - Visualisation
# ============================================

# Basic barplot
barplot(data$expression,
        names.arg = data$gene,
        las = 2,
        main = "Gene Expression",
        ylab = "Expression")

# ggplot2 barplot - more beautiful
# aes() = aesthetics = what goes on x axis, y axis, and fill color
ggplot(data, aes(x = gene, y = expression, fill = significant)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression",
       x = "Gene",
       y = "Expression")
