require(metap)

rank_p <- function(null, switched)
{
  rank <- which(sort(c(null, switched)) == switched)/length(null)
  p = 1
  if (rank > 0.5)
  {
    p = 2*(1-rank)
  }
  else
  {
    p = 2*rank
  }
  
  return(p)
}

setwd("~/Documents/PhD/Maela carriage length/")
metadata <- read.table("tree_metadata.csv", header=T, sep=",", stringsAsFactors = F)

switch1 <- read.table("capsule_switch/serotype23F_NT_6841_7_5.labels.txt", stringsAsFactors = F)

null_pheno <- metadata[metadata$name %in% switch1$V1 & metadata$name != '6841_7_5', 3]
switch_pheno <- metadata[metadata$name == '6841_7_5', 3]

p1 <- rank_p(null_pheno, switch_pheno)

switch2 <- read.table("capsule_switch/serotype34_NT_6649_7_23.labels.txt", stringsAsFactors = F)

null_pheno <- metadata[metadata$name %in% switch2$V1 & metadata$name != '6649_7_23', 3]
switch_pheno <- metadata[metadata$name == '6649_7_23', 3]

p2 <- rank_p(null_pheno, switch_pheno)

switch3 <- read.table("capsule_switch/serotype23F_NT_6736_7_18.labels.txt", stringsAsFactors = F)

null_pheno <- metadata[metadata$name %in% switch3$V1 & metadata$name != '6841_7_5' & metadata$name != '6736_7_18', 3]
switch_pheno <- metadata[metadata$name == '6736_7_18', 3]

p3 <- rank_p(null_pheno, switch_pheno)

sumlog(c(p1,p2,p3))
