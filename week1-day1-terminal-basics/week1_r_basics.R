# Week 1 - R Basics
# Day 1 - my first R script

# Variables
# <- means "assign" - store a value in a variable
x <- 10
x
x + 5
x * 2

# Vectors - a list of numbers
# c() = combine values together into one list
genes <- c(10, 25, 3, 47, 8)
genes

# Basic statistics on a vector
mean(genes)   # average = 18.6
sum(genes)    # total = 93
max(genes)    # highest value = 47
min(genes)    # lowest value = 3
# Strings - text vectors
# text in R always goes in quotes " "
gene_names <- c("BRCA1", "TP53", "EGFR", "MYC", "PTEN")
gene_names
# Data frame - a table connecting gene names to expression values
gene_table <- data.frame(
  gene = c("BRCA1", "TP53", "EGFR"),
  expression = c(10, 25, 3)
)
gene_table