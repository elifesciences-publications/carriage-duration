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

#mcmcfit <- mcmc.poumm(tree_phenotypes, pruned_tree, n.mcmc=2e5, n.adapt=20000, thin=100, acc.rate=0.1,scale=matrix(c(400, 0.00, 0.00, 0.00, 0.00,  2.00,  0.00,  0.00, 0.00,  0.00,  0.02,  0.00, 0.00,  0.00,  0.00,  0.02), nrow=4, ncol=4, byrow=T),distgr = 'maxlik')
#mcmcAnalysis <- analyseMCMCs(chains = mcmcfit$chains, stat = function(par) {H2e.poumm(z=tree_phenotypes, sigmae=sqrt(par[4]))}, statName='H2.OUe',start=1e5, end=2e5, thin=100)

mlfit <- ml.poumm(tree_phenotypes, pruned_tree, distgr = 'maxlik')
print(mlfit$par)
#alpha      theta      sigma     sigmae 
#22.2674543  4.3270402  3.3823372  0.7624445 
H2.poumm(mlfit$par[1], mlfit$par[3], mlfit$par[4], t=Inf, tm=0)
#[1] 0.3064671
H2e.poumm(z=tree_phenotypes, sigmae=mlfit$par[4])
#sigmae 
#0.3336053 

#> mcmcAnalysis$Mode
#[1] 0.2834562
#> mcmcAnalysis$HPD
#[[1]]
#lower     upper
#var1 0.2368475 0.3327728

saveRDS(tree_phenotypes, "tree_phenotypes.Rdata")
saveRDS(pruned_tree, "pruned_tree.Rdata")

#h^2
pruned_tree <- drop.tip(common_fasttree, common_fasttree$tip.label[!(common_fasttree$tip.label %in% maela_tree$tip.label)])
tree_phenotypes <- reg_pheno[match(common_fasttree$tip.label, names(reg_pheno))]

cpps <- analyseCPPs(tree_phenotypes, fasttree, CPPthr=0.04)
cpps$analysis.CPP[c('rA', 'bCI95lower', 'bCI95upper')]

#$rA
#[1] 0.438032
#
#$bCI95lower
#[1] 0.391244
#
#$bCI95upper
#[1] 0.4944363

#> mcmcAnalysis$Mode
#[1] 0.336443
#> mcmcAnalysis$HPD
#[[1]]
#lower    upper
#var1 0.2764059 0.403109
#attr(,"Probability")
#[1] 0.95

saveRDS(tree_phenotypes, "common_tree_phenotypes.Rdata")
saveRDS(pruned_tree, "common_pruned_tree.Rdata")
