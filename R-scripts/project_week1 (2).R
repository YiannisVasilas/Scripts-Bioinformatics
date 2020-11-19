library(class)
library(MASS)
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
########## REMOVE THIS LATER
testing <- log(training_full[,1:columnsToUse])
testing_test <- log(test_full[,1:columnsToUse])
testing.model <- lm(ENSG00000271043.1_MTRNR2L2 ~ ., data=training_full[,1:columnsToUse])
predictions.test <- predict(model.1.lm, newdata = testing_test)
summary(testing.model)
mean((test_full$ENSG00000271043.1_MTRNR2L2 - predictions.test)[-training]^2)
#################################################################################################################
model.1.lm <- lm(ENSG00000271043.1_MTRNR2L2 ~ ., data=training_full[,1:columnsToUse])
summary(model.1.lm)
predictions.1.lm <- predict(model.1.lm, newdata = test_full[,1:columnsToUse])
#MSE
mean((test_full$ENSG00000271043.1_MTRNR2L2 - predictions.1.lm)[-training]^2)
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
# to save time predictions were instantly made instead of first the model.
#splitting up later for readability
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

##############
##### 3 ######
##############
set.seed(5)
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
