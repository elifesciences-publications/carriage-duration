require(phytools)
fasttree <- read.tree("all_w_pheno.fasttree.tr")
rooted_tree <- midpoint.root(fasttree)
long_branches <- c("6680_4_4", "6731_2_11", "6775_1_8", "6823_1_6", "6775_2_16", "7038_8_40")
new_tree <- drop.tip(rooted_tree, long_branches)

warpedlmm_newnt <- read.delim("~/Documents/PhD/Maela carriage length/warpedlmm_newnt.pheno", header=FALSE, stringsAsFactors=FALSE)
tree_phenotypes <- warpedlmm_newnt$V3
names(tree_phenotypes) <- warpedlmm_newnt$V1
new_tree <- drop.tip(new_tree,
                     new_tree$tip.label[!new_tree$tip.label %in% names(tree_phenotypes)])
tree_phenotypes <- tree_phenotypes[names(tree_phenotypes) %in% new_tree$tip.label]

new_mapped <- contMap(new_tree, tree_phenotypes)
new_colour <- setMap(new_mapped, colors=c("blue","yellow"))
pdf("trait_mapped_new.pdf", width = 12, height = 12)
plot(new_colour, fsize = c(0.1, 1), type="fan", lwd=1)
dev.off()
