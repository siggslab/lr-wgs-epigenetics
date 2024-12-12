library(rrcov)
library(ggplot2)
library(ggplotify)
library(dplyr)
library(data.table)


setwd("<INSERT WORKING DIRECTORY>")

##### Do stats
meth_freqs = fread("./methfreqs.tsv", header = TRUE)
rownames(meth_freqs) = meth_freqs$SITE
meth_freqs$SITE = NULL
meth_freqs = t(meth_freqs)

##### Remove any columns with NA entries
cf_DFinf2NA <- function(x)
{
  for (i in 1:ncol(x)){
    x[,i][is.infinite(x[,i])] = NA
  }
  return(x)
}
cleaned = cf_DFinf2NA(meth_freqs)
cleaned_2 = cleaned[ , colSums(is.na(cleaned))==0]


# Run PCA
rpca = PcaGrid(cleaned_2, k=2)

# Save outlier-y results
pdf("./outlier_map.pdf", width=12, height=12)
plot(rpca)
dev.off()
apply()
