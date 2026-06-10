# Week 1 - Day 2 Notes
## Data Frames, CSV Files, and Plots

---

### What is a data frame?
A data frame is a table with rows and columns — like Excel in R.
Each column can have a different data type:
- text (chr) = gene names
- numbers (num) = expression values
- TRUE/FALSE (logi) = significant or not

---

### Data types in R
| Type | Example | Meaning |
|------|---------|---------|
| chr | "BRCA1" | text — always in quotes |
| num | 150 | numbers |
| logi | TRUE / FALSE | yes or no |

---

### Key functions

| Function | Meaning | Example |
|----------|---------|---------|
| data.frame() | create a table | data.frame(gene=c("BRCA1")) |
| nrow() | number of rows | nrow(gene_table) |
| ncol() | number of columns | ncol(gene_table) |
| str() | show structure and data types | str(gene_table) |
| head() | show first 6 rows | head(data) |
| mean() | calculate average | mean(data$expression) |
| read.csv() | load a CSV file into R | read.csv("file.csv") |
| setwd() | set working directory | setwd("~/Desktop/folder") |
| library() | load a package | library(ggplot2) |

---

### What is a CSV file?
CSV = Comma Separated Values
A table saved as plain text — commas separate columns
R loads it with: data <- read.csv("file.csv")

---

### What is ggplot2?
The most popular R plotting package
Used in every RNA-seq paper
ggplot(data, aes(x=gene, y=expression, fill=significant)) +
  geom_bar(stat="identity") +
  labs(title="Gene Expression")

---

### Vocabulary
- package = a collection of R functions
- library() = loads a package so you can use it
- install.packages() = downloads a package for the first time
- setwd() = tells R which folder to look in
- working directory = the folder R is currently looking in
- $ = accesses one column from a data frame
---

## Statistics — log2FC and p-value

### log2FC — how much a gene changed
- Positive log2FC = gene went UP in disease
- Negative log2FC = gene went DOWN in disease
- log2FC = 1 means gene doubled (2x increase)
- log2FC = 2 means gene quadrupled (4x increase)
- log2FC = -1 means gene halved (2x decrease)
- We use |log2FC| > 1 as cutoff = only care about genes that changed at least 2x

### Example
MYC in healthy = 10, MYC in cancer = 40
log2(40/10) = 2 → MYC went up 4x in cancer

### p-value — is the change real or just noise?
- p-value < 0.05 = the change is statistically significant = real
- p-value > 0.05 = could be random chance = not significant
- In DESeq2 we use padj (adjusted p-value) because we test 20,000 genes at once

### The DEG cutoff in DESeq2
A gene is called differentially expressed when:
- padj < 0.05 (statistically significant)
- |log2FC| > 1 (biologically meaningful change)
Both conditions must be TRUE