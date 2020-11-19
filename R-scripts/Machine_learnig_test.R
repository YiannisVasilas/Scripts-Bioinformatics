#Check for updates 
install.packages("caret")
install.packages("caret", dependencies=c("Depends", "Suggests"))
#Load the package
library(caret)
#Set your data here
Filename='data.csv'
colnames=c('name_of_your_colums')
# create a list of 80% of the rows in the original dataset we can use for training
validation_index = createDataPartition(dataset$Species, p=0.80, list=FALSE)
# select 20% of the data for validation
validation = dataset[-validation_index,]
# use the remaining 80% of data to training and testing the models
dataset = dataset[validation_index,]
# create a list of 80% of the rows in the original dataset we can use for training
validation_index = createDataPartition(dataset$Species, p=0.80, list=FALSE)
# dimensions of dataset
dim(dataset)
# Check if the attributes are correct
sapply(dataset, class)
# summarize the class distribution
percentage = prop.table(table(dataset$Species)) * 100
cbind(freq=table(dataset$Species), percentage=percentage)
# summarize attribute distributions
summary(dataset)
#Visualize
x = dataset[,1:4]
y = dataset[,5]
#Box plot without color set col!
par(mfrow=c(1,4))
for(i in 1:4) {
  boxplot(x[,i], main=names(iris)[i])
}
#Barplot set col here also
plot(y)
# scatterplot matrix
featurePlot(x=x, y=y, plot="ellipse")
# Machine learnig algorithms
# a) linear algorithms
set.seed(7)
fit.lda = train(Species~., data=dataset, method="lda", metric=metric, trControl=control)
# b) nonlinear algorithms
# CART
set.seed(7)
fit.cart = train(Species~., data=dataset, method="rpart", metric=metric, trControl=control)
# kNN
set.seed(7)
fit.knn = train(Species~., data=dataset, method="knn", metric=metric, trControl=control)
# c) advanced algorithms
# SVM
set.seed(7)
fit.svm = train(Species~., data=dataset, method="svmRadial", metric=metric, trControl=control)
# Random Forest
set.seed(7)
fit.rf = train(Species~., data=dataset, method="rf", metric=metric, trControl=control)
#Accuracy of the models
results = resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf))
summary(results)
#Visualize
dotplot(results)
print(fit.lda)
#predictions here I check the LDA model 
# check any model you waqnt
predictions <- predict(fit.lda, validation)
confusionMatrix(predictions, validation$Species)
