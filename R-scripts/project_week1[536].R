# install.packages("ISLR")
# install.packages("tree")
# install.packages("randomForest")
# install.packages("e1071")
library(class)
library(MASS)
library(ISLR)
library(tree)
library(randomForest)
library(e1071)

set.seed(2)
# Read in the data
datafile <- "M:\\Period_4\\project\\Data\\expr4T.dat"
dataUnfiltered <- read.table(datafile)



# Last column contains tissue label
nn=ncol(dataUnfiltered)-1
m=sapply(dataUnfiltered[,1:nn],mean)
s=sapply(dataUnfiltered[,1:nn],sd)
# The standard deviation plotted with a max of 1.5, it shows a
plot(m[order(m, decreasing = T)], ylim=c(3,500))
sm=s/m
plot(sm[order(sm, decreasing = T)], ylim = c(1.5, 2))
# Filter genes based on minimum mean and stdev/mean
minsm = 1
minm = 15

# In the following, we make sure that tissue label is kept;
# by giving it an artificial sufficiently high value for m and sm
m=c(m,minm+1)
sm=c(sm,minsm+1)
dim(dataUnfiltered)
length(which(sm>minsm & m>minm))
dataToUse <- dataUnfiltered[,which(sm>minsm & m>minm)]
dim(dataToUse)

# The check if it contains the to genes used in the next week
containsGenes <- "ENSG00000271043.1_MTRNR2L2" %in% colnames(dataToUse) & "ENSG00000229344.1_RP5.857K21.7" %in% colnames(dataToUse)

# 2: models
# extract 2 different tissues
tissue.1 = "brain_amygdala"
tissue.2 = "brain_cerebellum"
amygdala=dataToUse[which(dataToUse$tissue==tissue.1),]
amygdala=droplevels(amygdala)
cerebellum=dataToUse[which(dataToUse$tissue==tissue.2),]
cerebellum=droplevels(cerebellum)

columnsToUse <- ncol(dataToUse) - 1
# Create training and test data
# Full
training <- sample(nrow(dataToUse), size = 600, replace = FALSE)
training_full <- dataToUse[training,]
test_full <- dataToUse[-training,1:columnsToUse]
# Amygdala
training_amyg <- sample(nrow(amygdala), size = 36, replace = FALSE)
training_amygdala <- amygdala[training_amyg,]
test_amygdala <- amygdala[-training_amyg,1:columnsToUse]
# Cerebellum
training_cere <- sample(nrow(cerebellum), size = 62, replace = FALSE)
training_cerebellum <- cerebellum[training_cere,]
test_cerebellum <- cerebellum[-training_cere,1:columnsToUse]

# ENSG00000271043.1_MTRNR2L2
# 12 linear regression models
##########
## 2ABC ##
##########

model.1.lm <- lm(ENSG00000271043.1_MTRNR2L2 ~ ., data=training_full[,1:columnsToUse])
summary(model.1.lm)
predictions.1.lm <- predict(model.1.lm, newdata = test_full[,1:columnsToUse])
#MSE
mean((test_full$ENSG00000271043.1_MTRNR2L2 - predictions.1.lm)[-training]^2)
plot(studres(model.1.lm) ~ model.1.lm$fitted, ylab= "Studentized residuals", xlab="Fitted values", main="Normal")

###############
model.1.lm.amygdala <- lm(ENSG00000271043.1_MTRNR2L2 ~ ., data=training_amygdala[,1:columnsToUse])
summary(model.1.lm.amygdala)
predictions.1.lm.amygdala <- predict(model.1.lm.amygdala, newdata=test_amygdala[,1:columnsToUse])
mean((test_amygdala$ENSG00000271043.1_MTRNR2L2 - predictions.1.lm.amygdala)[-training_amyg]^2)
###############
model.1.lm.cerebellum <- lm(ENSG00000271043.1_MTRNR2L2 ~ ., data=training_cerebellum[,1:columnsToUse])
summary(model.1.lm.cerebellum)
predictions.1.lm.cerebellum <- predict(model.1.lm.cerebellum, newdata = test_cerebellum[, 1:columnsToUse])
mean((test_cerebellum$ENSG00000271043.1_MTRNR2L2 - predictions.1.lm.cerebellum)[-training_cere]^2)
###############
model.2.lm <- lm(ENSG00000229344.1_RP5.857K21.7 ~ ., data=training_full[,1:columnsToUse])
summary(model.2.lm)
predictions.2.lm <- predict(model.2.lm, newdata = test_full[,1:columnsToUse])
mean((test_full$ENSG00000229344.1_RP5.857K21.7 - predictions.2.lm)[-training]^2)
###############
model.2.lm.amygdala <- lm(ENSG00000229344.1_RP5.857K21.7 ~ ., data=training_amygdala[,1:columnsToUse])
predictions.2.lm.amygdala <- predict(model.2.lm.amygdala, newdata = test_amygdala[,1:columnsToUse])
mean((test_amygdala$ENSG00000229344.1_RP5.857K21.7 - predictions.2.lm.amygdala)[-training_amyg]^2)
###############
model.2.lm.cerebellum <- lm(ENSG00000229344.1_RP5.857K21.7 ~ ., data=training_cerebellum[,1:columnsToUse])
predictions.2.lm.cerebellum <- predict(model.2.lm.cerebellum, newdata = test_cerebellum[,1:columnsToUse])
mean((test_cerebellum$ENSG00000229344.1_RP5.857K21.7 - predictions.2.lm.cerebellum)[-training_cere]^2)
###############
model.3.lm <- lm(ENSG00000225972.1_MTND1P23 ~ ., data=training_full[,1:columnsToUse])
predictions.3.lm <- predict(model.3.lm, newdata = test_full[,1:columnsToUse])
summary(model.3.lm)
mean((test_full$ENSG00000225972.1_MTND1P23 - predictions.3.lm)[-training]^2)
###############
model.3.lm.amygdala <- lm(ENSG00000225972.1_MTND1P23 ~ ., data=training_amygdala[,1:columnsToUse])
summary(model.3.lm.amygdala)
predictions.3.lm.amygdala <- predict(model.3.lm.amygdala, newdata=test_amygdala[,1:columnsToUse])
mean((test_amygdala$ENSG00000225972.1_MTND1P23 - predictions.3.lm.amygdala)[-training_amyg]^2)
###############
model.3.lm.cerebellum <- lm(ENSG00000225972.1_MTND1P23 ~ ., data=training_cerebellum[,1:columnsToUse])
summary(model.3.lm.cerebellum)
predictions.3.lm.cerebellum <- predict(model.3.lm.cerebellum, newdata = test_cerebellum[, 1:columnsToUse])
mean((test_cerebellum$ENSG00000225972.1_MTND1P23 - predictions.3.lm.cerebellum)[-training_cere]^2)
###############
model.4.lm <- lm(ENSG00000173110.6_HSPA6 ~ ., data=training_full[,1:columnsToUse])
summary(model.4.lm)
predictions.4.lm <- predict(model.4.lm, newdata = test_full[,1:columnsToUse])
mean((test_full$ENSG00000173110.6_HSPA6 - predictions.4.lm)[-training]^2)
###############
model.4.lm.amygdala <- lm(ENSG00000173110.6_HSPA6 ~ ., data=training_amygdala[,1:columnsToUse])
summary(model.4.lm.amygdala)
predictions.4.lm.amygdala <- predict(model.4.lm.amygdala, newdata=test_amygdala[,1:columnsToUse])
mean((test_amygdala$ENSG00000173110.6_HSPA6 - predictions.4.lm.amygdala)[-training_amyg]^2)
###############
model.4.lm.cerebellum <- lm(ENSG00000173110.6_HSPA6 ~ ., data=training_cerebellum[,1:columnsToUse])
summary(model.4.lm.cerebellum)
predictions.4.lm.cerebellum <- predict(model.4.lm.cerebellum, newdata = test_cerebellum[, 1:columnsToUse])
mean((test_cerebellum$ENSG00000173110.6_HSPA6 - predictions.4.lm.cerebellum)[-training_cere]^2)

########
## 2D ##
########
# Add +1 to all numbers to get rid of the zeros before log transforming
# Suggestion by D. De ridder
log.training <- log(training_full[,1:columnsToUse]+1)
log.testing <- log(test_full[,1:columnsToUse]+1)
log.model <- lm(ENSG00000271043.1_MTRNR2L2 ~ ., data=log.training[,1:columnsToUse])
log.test <- predict(log.model, newdata = log.testing)
summary(log.model)
plot(studres(log.model) ~ log.model$fitted.values ,ylab= "Studentized residuals", xlab="Fitted values", main="Log transformed")
##########

##############
##### 3 ######
##############
tissue1="brain_amygdala"
tissue2="brain_hippocampus"
dataNew=data.frame(dataToUse)
dataNew=dataNew[which(dataNew$tissue==tissue1|dataNew$tissue==tissue2),]
dataNew=droplevels(dataNew)

training_1 <- sample(nrow(dataNew), size = nrow(dataNew)/2, replace = FALSE)
training_set  <- dataNew[training_1,]
test_set <- dataNew[-training_1,]
### make the models
glm.fit <- glm(tissue ~ ., data=training_set, family = binomial)
summary(glm.fit)
coef(glm.fit)
glm.probs <- predict(glm.fit, newData=test_set, type="response")
glm.probs
glm.pred <- rep("brain_amygdala", nrow(test_set))
glm.pred[glm.probs>.5] = "brain_hippocampus"
table(glm.pred, test_set$tissue)
###
(lda.fit <- lda(tissue ~ ., training_set))
lda.pred <- predict(lda.fit, test_set)
table(lda.pred$class, test_set$tissue)

knn.pred <- knn(training_set[,1:ncol(training_set)-1], test_set[,1:ncol(test_set)-1], training_set$tissue, k=3)
table(knn.pred, test_set$tissue)


#3.2
tissue1="brain_hippocampus"
tissue2="brain_nucleusaccumbens"
dataNew=data.frame(dataToUse)
dataNew=dataNew[which(dataNew$tissue==tissue1|dataNew$tissue==tissue2),]
dataNew=droplevels(dataNew)

training_1 <- sample(nrow(dataNew), size = nrow(dataNew)/2, replace = FALSE)
training_set  <- dataNew[training_1,]
test_set <- dataNew[-training_1,]
### make the models
glm.fit <- glm(tissue ~ ., data=training_set, family = binomial)
summary(glm.fit)
coef(glm.fit)
glm.probs <- predict(glm.fit, newData=test_set, type="response")
glm.probs
glm.pred <- rep("brain_hippocampus", nrow(test_set))
glm.pred[glm.probs>.5] = "brain_nucleusaccumbens"
table(glm.pred, test_set$tissue)
###
(lda.fit <- lda(tissue ~ ., training_set))
lda.pred <- predict(lda.fit, test_set)
table(lda.pred$class, test_set$tissue)

knn.pred <- knn(training_set[,1:ncol(training_set)-1], test_set[,1:ncol(test_set)-1], training_set$tissue, k=3)
table(knn.pred, test_set$tissue)


#3.3
tissue1="brain_spinalcord"
tissue2="brain_substantianigra"
dataNew=data.frame(dataToUse)
dataNew=dataNew[which(dataNew$tissue==tissue1|dataNew$tissue==tissue2),]
dataNew=droplevels(dataNew)

training_1 <- sample(nrow(dataNew), size = nrow(dataNew)/2, replace = FALSE)
training_set  <- dataNew[training_1,]
test_set <- dataNew[-training_1,]
### make the models
glm.fit <- glm(tissue ~ ., data=training_set, family = binomial)
summary(glm.fit)
coef(glm.fit)
glm.probs <- predict(glm.fit, newData=test_set, type="response")
glm.probs
glm.pred <- rep("brain_spinalcord", nrow(test_set))
glm.pred[glm.probs>.5] = " brain_substantianigra"
table(glm.pred, test_set$tissue)
###
(lda.fit <- lda(tissue ~ ., training_set))
lda.pred <- predict(lda.fit, test_set)
table(lda.pred$class, test_set$tissue)

knn.pred <- knn(training_set[,1:ncol(training_set)-1], test_set[,1:ncol(test_set)-1], training_set$tissue, k=3)
table(knn.pred, test_set$tissue)


#3.4
tissue1="brain_cerebellarhemisphere"
tissue2="brain_cerebellum"
dataNew=data.frame(dataToUse)
dataNew=dataNew[which(dataNew$tissue==tissue1|dataNew$tissue==tissue2),]
dataNew=droplevels(dataNew)

training_1 <- sample(nrow(dataNew), size = nrow(dataNew)/2, replace = FALSE)
training_set  <- dataNew[training_1,]
test_set <- dataNew[-training_1,]
### make the models
glm.fit <- glm(tissue ~ ., data=training_set, family = binomial)
summary(glm.fit)
coef(glm.fit)
glm.probs <- predict(glm.fit, newData=test_set, type="response")
glm.probs
glm.pred <- rep("brain_cerebellarhemisphere", nrow(test_set))
glm.pred[glm.probs>.5] = "brain_cerebellum"
table(glm.pred, test_set$tissue)
###
(lda.fit <- lda(tissue ~ ., training_set))
lda.pred <- predict(lda.fit, test_set)
table(lda.pred$class, test_set$tissue)

knn.pred <- knn(training_set[,1:ncol(training_set)-1], test_set[,1:ncol(test_set)-1], training_set$tissue, k=3)
table(knn.pred, test_set$tissue)


#3.5
tissue1="brain_cerebellum"
tissue2="brain_amygdala"

dataNew=data.frame(dataToUse)
dataNew=dataNew[which(dataNew$tissue==tissue1|dataNew$tissue==tissue2),]
dataNew=droplevels(dataNew)

training_1 <- sample(nrow(dataNew), size = nrow(dataNew)/2, replace = FALSE)
training_set  <- dataNew[training_1,]
test_set <- dataNew[-training_1,]
### make the models
glm.fit <- glm(tissue ~ ., data=training_set, family = binomial)
summary(glm.fit)
coef(glm.fit)
glm.probs <- predict(glm.fit, newData=test_set, type="response")
glm.probs
glm.pred <- rep("brain_cerebellum", nrow(test_set))
glm.pred[glm.probs>.5] = "brain_amygdala"

table(glm.pred, test_set$tissue)
###
(lda.fit <- lda(tissue ~ ., training_set))
lda.pred <- predict(lda.fit, test_set)
table(lda.pred$class, test_set$tissue)

knn.pred <- knn(training_set[,1:ncol(training_set)-1], test_set[,1:ncol(test_set)-1], training_set$tissue, k=3)
table(knn.pred, test_set$tissue)

#3.6 1st of our choice
tissue1="brain_putamen"
tissue2="brain_spinalcord"

dataNew=data.frame(dataToUse)
dataNew=dataNew[which(dataNew$tissue==tissue1|dataNew$tissue==tissue2),]
dataNew=droplevels(dataNew)

training_1 <- sample(nrow(dataNew), size = nrow(dataNew)/2, replace = FALSE)
training_set  <- dataNew[training_1,]
test_set <- dataNew[-training_1,]
### make the models
glm.fit <- glm(tissue ~ ., data=training_set, family = binomial)
summary(glm.fit)
coef(glm.fit)
glm.probs <- predict(glm.fit, newData=test_set, type="response")
glm.probs
glm.pred <- rep("brain_putamen", nrow(test_set))
glm.pred[glm.probs>.5] = "brain_spinalcord"

table(glm.pred, test_set$tissue)
###
(lda.fit <- lda(tissue ~ ., training_set))
lda.pred <- predict(lda.fit, test_set)
table(lda.pred$class, test_set$tissue)

knn.pred <- knn(training_set[,1:ncol(training_set)-1], test_set[,1:ncol(test_set)-1], training_set$tissue, k=3)
table(knn.pred, test_set$tissue)

#3.7 2nd of our choice
tissue1="brain_anteriorcortex"
tissue2="brain_hypothalamus"

dataNew=data.frame(dataToUse)
dataNew=dataNew[which(dataNew$tissue==tissue1|dataNew$tissue==tissue2),]
dataNew=droplevels(dataNew)

training_1 <- sample(nrow(dataNew), size = nrow(dataNew)/2, replace = FALSE)
training_set  <- dataNew[training_1,]
test_set <- dataNew[-training_1,]
### make the models
glm.fit <- glm(tissue ~ ., data=training_set, family = binomial)
summary(glm.fit)
coef(glm.fit)
glm.probs <- predict(glm.fit, newData=test_set, type="response")
glm.probs
glm.pred <- rep("brain_anteriorcortex", nrow(test_set))
glm.pred[glm.probs>.5] = "brain_hypothalamus"
table(glm.pred, test_set$tissue)
###
(lda.fit <- lda(tissue ~ ., training_set))
lda.pred <- predict(lda.fit, test_set)
table(lda.pred$class, test_set$tissue)
knn.pred <- knn(training_set[,1:ncol(training_set)-1], test_set[,1:ncol(test_set)-1], training_set$tissue, k=3)
table(knn.pred, test_set$tissue)

#########
### 4 ###
#########
## non linear fit to the data
glm.fit.nl <- glm(tissue ~ . + ENSG00000004799.7_PDK4^2 + ENSG00000225972.1_MTND1P23^2, data=training_set, family = binomial)
summary(glm.fit.nl)
coef(glm.fit.nl)
glm.probs.nl <- predict(glm.fit.nl, newData=test_set, type="response")

glm.pred.nl <- rep("brain_anteriorcortex", nrow(test_set))
glm.pred.nl[glm.probs.nl>.5] = "brain_hypothalamus"
table(glm.pred.nl, test_set$tissue)

(lda.fit.nl <- lda(tissue ~ . + ENSG00000004799.7_PDK4^2 + ENSG00000225972.1_MTND1P23^2, training_set))
lda.pred.nl <- predict(lda.fit, test_set)
table(lda.pred.nl$class, test_set$tissue)

####
## It does not change anything because the amount of 
# predictors is probably to large for the non linear terms to have a real impact

#############################
########## WEEK 2 ###########
#############################
### Decision tree
### retrain data
#Week 2 question 1

# regression tree
tree.dataToUse1 <- tree(ENSG00000271043.1_MTRNR2L2 ~ ., data=dataToUse[,1:columnsToUse], subset = training)
#summary indicates the 4 variables that are used in constructing the tree
summary(tree.dataToUse1)
plot(tree.dataToUse1)
text(tree.dataToUse1, pretty = 0, cex=0.7)
tree.dataToUse1
#cv to see whether pruning the tree will improve performance
cv.dataToUse1 <- cv.tree(tree.dataToUse1)
plot(cv.dataToUse1$size, cv.dataToUse1$dev, type = "b")
#for pruning the tree we use the following, nothing changes though
pruned.dataToUse1 <- prune.tree(tree.dataToUse1, best = 5)
plot(pruned.dataToUse1)
text(pruned.dataToUse1, pretty = 0, cex=1)
#pruned tree to make predictions on the test set 
dataToUse1.pred <- predict(pruned.dataToUse1, newdata=dataToUse[-training,])
dataToUse1.test <- dataToUse[-training, "ENSG00000271043.1_MTRNR2L2"]
plot(dataToUse1.pred,  dataToUse1.test)
abline(0,1)
mean((dataToUse1.pred)^2)

#Random Forest -> classification, mtry is the number of predictors, correct ~ 0.80
rf.cl.train <- randomForest(as.factor(tissue) ~ . ,data = training_full, importance=T, ntree=1000)
rf.cl.test <- predict(rf.cl.train, newdata = test_full[,1:columnsToUse])
varImpPlot(rf.cl.train)
table(rf.cl.test, dataToUse[-training,]$tissue)

### SVM , linear
tune.svmfit <- tune(svm, tissue ~ ., data=training_full, kernel="linear", ranges=c(0.001, 0.01, 0.1, 1, 10, 1000))
summary(tune.svmfit$best.model)
svm.pred.cl <- predict(tune.svmfit$best.model, test_full)
table(svm.pred.cl, dataToUse[-training,]$tissue)

#polynomial
tune.svmfit.rd <- tune(svm, tissue ~ ., data=training_full, kernel="polynomial", ranges=c(0.001, 0.01, 0.1, 1, 10, 1000), gamma=c(0.5, 1, 2, 3, 4))
summary(tune.svmfit.rd)
svm.pred.rd <- predict(tune.svmfit.rd$best.model, test_full)
table(svm.pred.rd, dataToUse[-training,]$tissue)
