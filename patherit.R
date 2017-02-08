require(patherit)

fasttree <- read.tree("all_w_pheno.fasttree.tr")
common_fasttree <- read.tree("all_w_pheno.common_only.fasttree.tr")
maela_tree <- read.tree("~/Documents/PhD/claire_data/no_mitis_nj_tree.tre")

reg_pheno <- metadata_out$phenotype
names(reg_pheno) <- sub("#", "_", metadata_out$lane)

# H^2
pruned_tree <- drop.tip(fasttree, fasttree$tip.label[!(fasttree$tip.label %in% maela_tree$tip.label)])
maela_dists <- cophenetic.phylo(pruned_tree)
hist(maela_dists[upper.tri(maela_dists)], main="Maela children genetic distances", xlab="Patristic distance")
tree_phenotypes <- reg_pheno[match(fasttree$tip.label, names(reg_pheno))]

cpps <- analyseCPPs(tree_phenotypes, fasttree, CPPthr=0.04)
cpps$analysis.CPP[c('rA', 'bCI95lower', 'bCI95upper')]

#$rA
#[1] 0.6392637
#
#$bCI95lower
#[1] 0.5921743
#
#$bCI95upper
#[1] 0.6857109

