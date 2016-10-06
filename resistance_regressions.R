require(glmnet)
require(dplyr)

# Set up data
setwd("~/Documents/PhD/Maela carriage length/")

mds_components <- read.csv("~/Documents/PhD/Maela carriage length/struct.covars", 
                           sep="", stringsAsFactors=FALSE)
resistances <- read.delim("~/Documents/PhD/Maela carriage length/Drug_resistant.txt", stringsAsFactors=FALSE)
sequence_metadata <- readRDS("sequence_metadata.Rdata")

covariates <- read.delim("~/Documents/PhD/Maela carriage length/covariates.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(covariates) <- c("IID", "copy2", "age", "carried")

metadata_out <- readRDS("metadata_out.Rdata")

merge1<-merge(mds_components, resistances, by.x = "IID", by.y = "lane_id")
merge2<-merge(merge1, metadata_out,by.x = "IID", by.y = "lane")
merge3<-merge(merge2, sequence_metadata, by.x="IID", by.y = "lane")
merge4<-merge(merge3, covariates, by="IID")

merge4<-merge4[-3,] # Remove outlier

# From naive/hypoth-free GWAS, looks like Penicillin and Trimethoprim resistance elements have an effect

# Penicillin effect w/ pop structure
plot(as.factor(merge4$Penicillin), merge4$phenotype)
lm.pen <- lm(data = merge4, phenotype ~ Penicillin)
summary(lm.pen)

glm.pen <- glm(data = merge4, phenotype ~ as.factor(Penicillin) + MDS1 + MDS2 + MDS3)
summary(glm.pen)

penicillin <- model.matrix(merge4$phenotype ~ as.factor(merge4$Penicillin))[,-1]
glmmod.pen<-glmnet(as.matrix(data.frame(penicillin, merge4[,3:5])), merge4$phenotype, alpha = 1)
plot(glmmod.pen,xvar="lambda",label = T)
cv.glmmod.pen <- cv.glmnet(as.matrix(data.frame(penicillin, merge4[,3:5])), merge4$phenotype, alpha = 1)
plot(cv.glmmod.pen)

# All drug effects
glm.all <- glm(data = merge4, phenotype ~ as.factor(Penicillin) + 
                 as.factor(Chloramphenicol) + as.factor(Clindamycin) + 
                 as.factor(Erythromycin) + as.factor(Sulpha.trimethoprim) + 
                 as.factor(Tetracycline))
summary(glm.all)
plot(glm.all)





order <- c("SENSITIVE", "INTERMEDIATE", "RESISTANT")
drugs <- model.matrix(merge4$phenotype ~ 
                        ordered(merge4$Chloramphenicol, levels=c("SENSITIVE", "RESISTANT")) + 
                        ordered(merge4$Clindamycin, levels=order) + 
                        ordered(merge4$Erythromycin, levels=order) + 
                        ordered(merge4$Sulpha.trimethoprim, levels=order) + 
                        ordered(merge4$Penicillin, levels=c("SENSITIVE", "MIC")) + 
                        ordered(merge4$Tetracycline, levels=order))[,-1]
#drugs <- model.matrix(merge4$phenotype ~ 
#                        factor(merge4$Chloramphenicol, levels=c("SENSITIVE", "RESISTANT")) + 
#                        factor(merge4$Clindamycin, levels=order) + 
#                        factor(merge4$Erythromycin, levels=order) + 
#                        factor(merge4$Sulpha.trimethoprim, levels=order) + 
#                        factor(merge4$Penicillin, levels=c("SENSITIVE", "MIC")) + 
#                        factor(merge4$Tetracycline, levels=order))[,-1]
glmmod.drugs <- glmnet(as.matrix(drugs), merge4$phenotype, alpha = 1)
plot(glmmod.drugs,xvar="lambda",label=T)
cv.glmmod.drugs <- cv.glmnet(as.matrix(drugs), merge4$phenotype, alpha = 1)
plot(cv.glmmod.drugs)
coef(cv.glmmod.drugs, s="lambda.1se")

# If age is added in, it needs to be scaled appropriately (does seem to have an important effect after CV)
# Increased age -> increased carriage length
# Previous carriage smaller effect
# Carried before -> 

# Serotypes
serotype <- model.matrix(merge4$phenotype ~ as.factor(merge4$serotype))[,-1]
glmmod.sero <- glmnet(as.matrix(serotype), merge4$phenotype, alpha = 1)
plot(glmmod.sero,xvar="lambda",label=T)
cv.glmmod.sero <- cv.glmnet(as.matrix(serotype), merge4$phenotype, alpha = 1)
plot(cv.glmmod.sero)
coef(cv.glmmod.sero, s="lambda.1se")

# Resistances and serotypes
all <- model.matrix(merge4$phenotype ~ 
                      ordered(merge4$Chloramphenicol, levels=c("SENSITIVE", "RESISTANT")) + 
                      ordered(merge4$Clindamycin, levels=order) + 
                      ordered(merge4$Erythromycin, levels=order) + 
                      ordered(merge4$Sulpha.trimethoprim, levels=order) + 
                      ordered(merge4$Penicillin, levels=c("SENSITIVE", "MIC")) + 
                      ordered(merge4$Tetracycline, levels=order) + as.factor(merge4$serotype))[,-1]
glmmod.all <- glmnet(as.matrix(all), merge4$phenotype, alpha = 1)
plot(glmmod.all,xvar="lambda",label=T)
cv.glmmod.all <- cv.glmnet(as.matrix(all), merge4$phenotype, alpha = 1, nfolds = 2195)
plot(cv.glmmod.all)
coef(cv.glmmod.all, s="lambda.1se")
predict(glmmod.all,type="coef",s=0.1) # ab res enter model first

# to find out which order predictors enter the model
# predict(glmmod.all,type="coef",s=glmmod.all$lambda)
predict(glmmod.all,type="coef",s=glmmod.all)

# use these in regression
selected <- all[,which(coef(cv.glmmod.all, s="lambda.1se")[-1] != 0)]
lm.all <- lm(merge4$phenotype ~ selected)
summary(lm.all)

# print coefficients
effects = exp(coef(cv.glmmod.all, s="lambda.1se")[1,1]+coef(cv.glmmod.all, s="lambda.1se")) - 
  exp(coef(cv.glmmod.all, s="lambda.1se")[1,1])
effect_table = data.frame(factor=rownames(effects)[-1],effect_size=effects[-1,1])
write.table(effect_table, file="lineage_effect_sizes.txt", 
            quote = F, row.names = F, col.names = F, sep=",")

# multiple R-squared = 0.2279. Direction of effects in table

# (Could also do this with BAPS clusters rather than serotype)
maela.pop.cluster1 <- read.delim("~/Documents/PhD/Maela carriage length/maela.pop.cluster1.phe", header=FALSE, stringsAsFactors=FALSE)
merge5<-merge(merge4, maela.pop.cluster1, by.x="IID",by.y="V2")

baps <- model.matrix(merge5$phenotype ~ merge5$V3)[,-1]
glmmod.baps <- glmnet(as.matrix(baps), merge5$phenotype, alpha = 1)
plot(glmmod.baps,xvar="lambda",label=T)
cv.glmmod.baps <- cv.glmnet(as.matrix(baps), merge5$phenotype, alpha = 1, nfolds = 2157)
plot(cv.glmmod.baps)
coef(cv.glmmod.baps, s="lambda.1se")

selected.baps <- baps[,which(coef(cv.glmmod.baps, s="lambda.1se")[-1] != 0)]
lm.baps <- lm(merge5$phenotype ~ selected.baps)
summary(lm.baps)

# environmental effects; Multiple R-squared:  0.04139
summary(lm(phenotype ~ age + carried, data=merge4))

# with direct swab data
all_X <- readRDS("all_X.Rdata")
all_y <- readRDS("all_y.Rdata")

filtered_X <- all_X[apply(is.na(all_X), 1, sum) < 2,]
filtered_y <- all_y[apply(is.na(all_X), 1, sum) < 2,]

design_X <- model.matrix(filtered_y ~ 
                      ordered(filtered_X$Chloramphenicol, levels=order) + 
                      ordered(filtered_X$Clindamycin, levels=order) + 
                      ordered(filtered_X$Erythromycin, levels=order) + 
                      ordered(filtered_X$Sulpha.trimethoprim, levels=order) + 
                      ordered(filtered_X$Penicillin, levels=order) + 
                      ordered(filtered_X$Tetracycline, levels=order) + 
                      relevel(as.factor(filtered_X$serotype), ref = "6A/C"))[,-1]
all_lm <- lm(log(filtered_y) ~ as.matrix(design_X))

pdf("lm_diagnostics.pdf", width = 20, height = 14)
par(mfrow=c(2,2))
plot(all_lm)

dev.off()

outliers <- c(which(filtered_y < 3), 1653, 2013, 2066, 2129, 2147, 2148, 2149)
lasso_X <- as.matrix(design_X[-outliers,])
lasso_y <- log(filtered_y[-outliers])

glmmod.all <- glmnet(lasso_X, lasso_y, alpha = 1)
plot(glmmod.all,xvar="lambda",label=T)
cv.glmmod.all <- cv.glmnet(lasso_X, lasso_y, alpha = 1, nfolds = length(lasso_y))
plot(cv.glmmod.all)
coef(cv.glmmod.all, s="lambda.1se")

# use selected predictors in regression
selected <- lasso_X[,which(coef(cv.glmmod.all, s="lambda.1se")[-1] != 0)]
lm.all <- lm(lasso_y ~ selected)
summary(lm.all)

# R^2 = 0.1895

# print coefficients
effects = exp(coef(cv.glmmod.all, s="lambda.1se")[1,1]+coef(cv.glmmod.all, s="lambda.1se")) - 
  exp(coef(cv.glmmod.all, s="lambda.1se")[1,1])
effect_table = data.frame(factor=rownames(effects)[-1],effect_size=effects[-1,1])
write.table(effect_table, file="lineage_effect_sizes.txt", 
            quote = F, row.names = F, col.names = F, sep=",")

# environment
summary(lm(lasso_y ~ filtered_X$age_d[-outliers] + filtered_X$carried[-outliers]))

# R^2 = 0.04594
