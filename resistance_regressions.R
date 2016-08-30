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
# (Could also do this with BAPS clusters rather than serotype)
all <- model.matrix(merge4$phenotype ~ 
                      ordered(merge4$Chloramphenicol, levels=c("SENSITIVE", "RESISTANT")) + 
                      ordered(merge4$Clindamycin, levels=order) + 
                      ordered(merge4$Erythromycin, levels=order) + 
                      ordered(merge4$Sulpha.trimethoprim, levels=order) + 
                      ordered(merge4$Penicillin, levels=c("SENSITIVE", "MIC")) + 
                      ordered(merge4$Tetracycline, levels=order) + as.factor(merge4$serotype))[,-1]
glmmod.all <- glmnet(as.matrix(all), merge4$phenotype, alpha = 1)
plot(glmmod.all,xvar="lambda",label=T)
cv.glmmod.all <- cv.glmnet(as.matrix(all), merge4$phenotype, alpha = 1)
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

# adjusted R-squared = 0.2212. Direction of effects in table?
