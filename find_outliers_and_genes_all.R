### Finding regions under selection
# a quick R script to find genes near outliers

rm(list = ls())

library(tidyverse)

# read in divergence stats
div <- read_csv("Pundamilia_kivu_div_stats.csv")

# do some renaming
div <- rename(div, chr = scaffold)

# get the Fst data
div <- div %>% select(chr, mid, fst = Fst_NyerMak_PundMak)

# remove chr10 - sex chromosome
div <- div %>% filter(chr != "chr10")

# find outliers
a <- ggplot(div, aes(x = fst)) + geom_histogram(fill = "white", colour = "black")
a

# what about the quantiles?
thresholds <- quantile(div$fst, c(0.975, 0.995))

# plot the quantiles on the distribution
a + geom_vline(xintercept = thresholds, lty = 2, colour = "red")

# set a new outlier variable
div <- div %>% mutate(status = ifelse(fst >= thresholds[1], "outlier", "background"))

# get data exceeding thresholds
outliers <- filter(div, fst >= thresholds[1])

# read in the gene data
genes <- read_csv2("NCBI_tilapia_nyererei_mRNA.csv")
genes <- rename(genes, chr = CHR)

# calculate midpoints 
genes <- genes %>% mutate(mid = start + round(end/2))

gene_hits <- sapply(1:nrow(outliers), function(z){
  target <- outliers[z, ]
  # define variables
  my_chr <- target$chr
  my_mid <- target$mid
  # find gene hits
  hits <- genes %>% filter(chr == my_chr, abs(my_mid - mid) <= 100000)
  # get window id
  windows <- target %>% slice(rep(1, each = nrow(hits)))
  bind_cols(hits, windows)
}, simplify = FALSE)

gene_hits <- bind_rows(gene_hits)

gene_hits %>% filter(chr == "chr17") %>% pull(gene)



