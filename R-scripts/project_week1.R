library(class)

# Read in the data
datafile <- "D:\\expr4T\\expr4T.dat"
dataUnfiltered <- read.table(datafile)

# Filter genes based on minimum mean and stdev/mean
minsm = 1.5
minm = 3

# Last column contains tissue label
nn=ncol(dataUnfiltered)-1
m=sapply(dataUnfiltered[,1:nn],mean)
s=sapply(dataUnfiltered[,1:nn],sd)
sm=s/m

# In the following, we make sure that tissue label is kept;
# by giving it an artificial sufficiently high value for m and sm
m=c(m,minm+1)
sm=c(sm,minsm+1)
dim(dataUnfiltered)
length(which(sm>minsm & m>minm))
dataToUse=dataUnfiltered[,which(sm>minsm & m>minm)]
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
model.1.lm <- lm(ENSG00000271043.1_MTRNR2L2 ~ . - ENSG00000271043.1_MTRNR2L2, data=training_full[,1:columnsToUse])
summary(model.1.lm)
predictions.1.lm <- predict(model.1.lm, newdata = test_full[,1:columnsToUse])
#MSE
mean((test_full$ENSG00000271043.1_MTRNR2L2 - predictions.1.lm)[-training]^2)
###############
model.1.lm.amygdala <- lm(ENSG00000271043.1_MTRNR2L2 ~ . - ENSG00000271043.1_MTRNR2L2, data=training_amygdala[,1:columnsToUse])
summary(model.1.lm.amygdala)
predictions.1.lm.amygdala <- predict(model.1.lm.amygdala, newdata=test_amygdala[,1:columnsToUse])
mean((test_amygdala$ENSG00000271043.1_MTRNR2L2 - predictions.1.lm.amygdala)[-training_amyg]^2)
###############
model.1.lm.cerebellum <- lm(ENSG00000271043.1_MTRNR2L2 ~ . - ENSG00000271043.1_MTRNR2L2, data=training_cerebellum[,1:columnsToUse])
summary(model.1.lm.cerebellum)
predictions.1.lm.cerebellum <- predict(model.1.lm.cerebellum, newdata = test_cerebellum[, 1:columnsToUse])
mean((test_cerebellum$ENSG00000271043.1_MTRNR2L2 - predictions.1.lm.cerebellum)[-training_cere]^2)
###############
model.2.lm <- lm(ENSG00000229344.1_RP5.857K21.7 ~ . - ENSG00000229344.1_RP5.857K21.7, data=training_full[,1:columnsToUse],)
predictions.2.lm <- predict(model.2.lm, newdata = test_full[,1:columnsToUse])
mean((test_full$ENSG00000229344.1_RP5.857K21.7 - predictions.2.lm)[-training]^2)
###############
model.2.lm.amygdala <- lm(ENSG00000229344.1_RP5.857K21.7 ~ . - ENSG00000229344.1_RP5.857K21.7, data=training_amygdala[,1:columnsToUse])
predictions.2.lm.amygdala <- predict(model.2.lm.amygdala, newdata = test_amygdala[,1:columnsToUse])
mean((test_amygdala$ENSG00000229344.1_RP5.857K21.7 - predictions.2.lm.amygdala)[-training_amyg]^2)
###############
model.2.lm.cerebellum <- lm(ENSG00000229344.1_RP5.857K21.7 ~ . - ENSG00000229344.1_RP5.857K21.7, data=training_cerebellum[,1:columnsToUse])
predictions.2.lm.cerebellum <- predict(model.2.lm.cerebellum, newdata = test_cerebellum[,1:columnsToUse])
mean((test_cerebellum$ENSG00000229344.1_RP5.857K21.7 - predictions.2.lm.amygdala)[-training_cere]^2)
###############
# to save time predictions were instantly made instead of first the model.
#splitting up later for readability
predictions.3.lm <- predict(lm(ENSG00000225972.1_MTND1P23 ~ . - ENSG00000225972.1_MTND1P23, data=training_full[,1:columnsToUse],), newdata = test_full[,1:columnsToUse])
mean((test_full$ENSG00000225972.1_MTND1P23 - predictions.3.lm)[-training]^2)
###############
predictions.3.lm.amygdala <- predict(lm(ENSG00000225972.1_MTND1P23 ~ . - ENSG00000225972.1_MTND1P23, data=training_amygdala[,1:columnsToUse]), newdata=test_amygdala[,1:columnsToUse])
mean((test_amygdala$ENSG00000225972.1_MTND1P23 - predictions.3.lm.amygdala)[-training_amyg]^2)
###############
predictions.3.lm.cerebellum <- predict(lm(ENSG00000225972.1_MTND1P23 ~ . - ENSG00000225972.1_MTND1P23, data=training_cerebellum[,1:columnsToUse]), newdata = test_cerebellum[, 1:columnsToUse])
mean((test_cerebellum$ENSG00000225972.1_MTND1P23 - predictions.3.lm.cerebellum)[-training_cere]^2)
###############
predictions.4.lm <- predict(lm(ENSG00000173110.6_HSPA6 ~ . - ENSG00000173110.6_HSPA6, data=training_full[,1:columnsToUse],), newdata = test_full[,1:columnsToUse])
mean((test_full$ENSG00000173110.6_HSPA6 - predictions.4.lm)[-training]^2)
###############
predictions.4.lm.amygdala <- predict(lm(ENSG00000173110.6_HSPA6 ~ . - ENSG00000173110.6_HSPA6, data=training_amygdala[,1:columnsToUse]), newdata=test_amygdala[,1:columnsToUse])
mean((test_amygdala$ENSG00000173110.6_HSPA6 - predictions.4.lm.amygdala)[-training_amyg]^2)
###############
predictions.4.lm.cerebellum <- predict(lm(ENSG00000173110.6_HSPA6 ~ . - ENSG00000173110.6_HSPA6, data=training_cerebellum[,1:columnsToUse]), newdata = test_cerebellum[, 1:columnsToUse])
mean((test_cerebellum$ENSG00000173110.6_HSPA6 - predictions.4.lm.cerebellum)[-training_cere]^2)

max(test_full$ENSG00000173110.6_HSPA6)

# model.1.glm <- glm(ENSG00000271043.1_MTRNR2L2 ~ . - ENSG00000271043.1_MTRNR2L2, data=dataToUse[,1:columnsToUse])
# model.1.knn <- knn(train = training_full, test = test_full, cl = as.factor(dataToUse$tissue[1:600]), k=13)
# brain_amygdala

# model.1.glm.amygdala <- glm(ENSG00000271043.1_MTRNR2L2 ~ . - ENSG00000271043.1_MTRNR2L2, data=amygdala[,1:columnsToUse])
# model.1.knn.amygdala <- knn(train = training_full, test = test_full, cl = as.factor(dataToUse$tissue[1:600]), k=1)
# brain_cerebellum

# model.1.glm.cerebellum <- glm(ENSG00000271043.1_MTRNR2L2 ~ . - ENSG00000271043.1_MTRNR2L2, data=cerebellum[,1:columnsToUse])
# model.1.knn.cerebellum <- knn(train = training_full, test = test_full, cl = as.factor(dataToUse$tissue[1:600]), k=1)

# ENSG00000229344.1_RP5.857K21.7
# non tissue specific

# model.2.glm <- glm(ENSG00000229344.1_RP5.857K21.7 ~ . - ENSG00000229344.1_RP5.857K21.7, data=dataToUse[,1:columnsToUse])
# model.2.knn <- knn(ENSG00000229344.1_RP5.857K21.7 ~ . - ENSG00000229344.1_RP5.857K21.7, data=dataToUse[,1:columnsToUse])
# brain_amygdala

# model.2.glm.amygdala <- glm(ENSG00000229344.1_RP5.857K21.7 ~ . - ENSG00000229344.1_RP5.857K21.7, data=amygdala[,1:columnsToUse])
# model.2.knn.amygdala <- knn(ENSG00000229344.1_RP5.857K21.7 ~ . - ENSG00000229344.1_RP5.857K21.7, data=amygdala[,1:columnsToUse])
# brain_cerebellum

# model.2.glm.cerebellum <- glm(ENSG00000229344.1_RP5.857K21.7 ~ . - ENSG00000229344.1_RP5.857K21.7, data=cerebellum[,1:columnsToUse])
# model.2.knn.cerebellum <- knn(ENSG00000229344.1_RP5.857K21.7 ~ . - ENSG00000229344.1_RP5.857K21.7, data=cerebellum[,1:columnsToUse])

# ENSG00000225972.1_MTND1P23
# non tissue specific

# model.3.glm <- glm(ENSG00000225972.1_MTND1P23 ~ . - ENSG00000225972.1_MTND1P23, data=dataToUse[,1:columnsToUse])
# model.3.knn <- knn(ENSG00000225972.1_MTND1P23 ~ . - ENSG00000225972.1_MTND1P23, data=dataToUse[,1:columnsToUse])
# brain_amygdala

# model.3.glm.amygdala <- glm(ENSG00000225972.1_MTND1P23 ~ . - ENSG00000225972.1_MTND1P23, data=amygdala[,1:columnsToUse])
# model.3.knn.amygdala <- knn(ENSG00000225972.1_MTND1P23 ~ . - ENSG00000225972.1_MTND1P23, data=amygdala[,1:columnsToUse])
# brain_cerebellum

# model.3.glm.cerebellum <- glm(ENSG00000225972.1_MTND1P23 ~ . - ENSG00000225972.1_MTND1P23, data=cerebellum[,1:columnsToUse])
# model.3.knn.cerebellum <- knn(ENSG00000225972.1_MTND1P23 ~ . - ENSG00000225972.1_MTND1P23, data=cerebellum[,1:columnsToUse])

# ENSG00000173110.6_HSPA6
# non tissue specific

# model.4.glm <- glm(ENSG00000173110.6_HSPA6 ~ . - ENSG00000173110.6_HSPA6, data=dataToUse[,1:columnsToUse])
# model.4.knn <- knn(ENSG00000173110.6_HSPA6 ~ . - ENSG00000173110.6_HSPA6, data=dataToUse[,1:columnsToUse])
# brain_amygdala

# model.4.glm.amygdala <- glm(ENSG00000173110.6_HSPA6 ~ . - ENSG00000173110.6_HSPA6, data=amygdala[,1:columnsToUse])
# model.4.knn.amygdala <- knn(ENSG00000173110.6_HSPA6 ~ . - ENSG00000173110.6_HSPA6, data=amygdala[,1:columnsToUse])
# brain_cerebellum

# model.4.glm.cerebellum <- glm(ENSG00000173110.6_HSPA6 ~ . - ENSG00000173110.6_HSPA6, data=cerebellum[,1:columnsToUse])
# model.4.knn.cerebellum <- knn(ENSG00000173110.6_HSPA6 ~ . - ENSG00000173110.6_HSPA6, data=cerebellum[,1:columnsToUse])

