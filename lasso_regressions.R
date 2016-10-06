require(glmnet)
require(parallel)
require(doMC)
require(lars)
require(covTest)

registerDoMC(cores=2)

setwd("~/Documents/PhD/Maela carriage length/post_test/")

# Sites
sites <- read.table("significant_sites.gen", sep = "\t", header = F, stringsAsFactors = F)
sites <- t(sites)
colnames(sites) <- sites[1,]
sites <- sites[-1,]

# Age and carried
covars <- read.table("baps.covars", header=T, stringsAsFactors = F)
covars <- covars[,c(2:4)]

# baps clusters
baps <- read.table("maela.pop.cluster1.phe", sep = "\t", header=F, stringsAsFactors = F)
colnames(baps) <- c("FID", "IID", "BAPS")

# phenotype (y)
phenotype <- read.table("length.pheno", sep = "\t", header=F, stringsAsFactors = F)
colnames(phenotype) <- c("IID", "pheno")

# Collect data
all_data <- merge(covars, baps)
all_data <- merge(all_data, phenotype)
all_data <- merge(all_data, sites, by.x="IID", by.y="ps")

post_pvals <- rep(1, ncol(all_data) - 6)

y <- all_data$pheno
for (i in 7:ncol(all_data))
{
  X <- model.matrix(y ~ factor(all_data$BAPS) + all_data$age +
                      all_data$carried + 
                      all_data[,i])[,-1]
  #glmmod <- glmnet(as.matrix(X), y, alpha = 1)
  #cv.glmmod <- cv.glmnet(as.matrix(X), y, alpha = 1, nfolds = nrow(X) - 1)

  # graphics
  #plot(glmmod,xvar="lambda",label=T)
  
  # to find out which order predictors enter the model
  #predict(glmmod,type="coef")
  #coef(glmmod, s=cv.glmmod$lambda.1se)
  
  # using significance test on results
  # https://projecteuclid.org/euclid.aos/1400592161
  elastic_new <- lars(X, y)
  lasso_sig <- covTest(elastic_new, X, y)
  
  post_pvals[i-6] <- lasso_sig$results[which(lasso_sig$results[,1] == ncol(X)),3]
}

post_test <- data.frame(sites=colnames(all_data[-1:-6]), pvals=post_pvals)
write.table(post_test[p.adjust(post_test$pvals) < 0.05,], 
            file = "post_test_significant.txt", quote = F, sep="\t", row.names = F,
            col.names = F)
