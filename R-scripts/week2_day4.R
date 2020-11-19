# 9.6.1 - SVC

# Generate Gaussian data, slightly shift one class and plot it

set.seed(1)
x=matrix(rnorm(20*2), ncol=2)
y=c(rep(-1,10),rep(1,10))
x[y==1,]=x[y==1,]+1
plot(x,col=(3-y))

# Train a linear SVC and plot it
dat=data.frame(x=x,y=as.factor(y))
library(e1071)
svmfit=svm(y~.,data=dat,kernel="linear",cost=10,scale=FALSE)
plot(svmfit,dat)

# What are the support vectors?
svmfit$index
summary(svmfit)

# What if we lower the cost? Note: cost here is the opposite of C,
# the lower the cost (of errors), the wider the margin
svmfit=svm(y~.,data=dat,kernel="linear",cost=0.1,scale=FALSE)
plot(svmfit,dat)
svmfit$index

# Tune cost using 10-fold CV
set.seed(1)
tune.out=tune(svm,y~.,data=dat,kernel="linear",ranges=list(cost=c(0.001, 0.01, 0.1, 1,5,10,100)))
summary(tune.out)

# Best model is stored as well
bestmod=tune.out$best.model
summary(bestmod)

# Create some test data and test it
xtest=matrix(rnorm(20*2), ncol=2)
ytest=sample(c(-1,1), 20, rep=TRUE)
xtest[ytest==1,]=xtest[ytest==1,] + 1
testdat=data.frame(x=xtest, y=as.factor(ytest))
ypred=predict(bestmod ,testdat)
table(predict=ypred, truth=testdat$y)

# Test for a lower cost
svmfit=svm(y~.,data=dat,kernel="linear",cost=.01,scale=FALSE)
ypred=predict(svmfit,testdat)
table(predict=ypred,truth=testdat$y)

# Make classes linearly separable and try again (extremely high cost!)
x[y==1,]=x[y==1,]+0.5
plot(x, col=(y+5)/2, pch=19)
dat=data.frame(x=x,y=as.factor(y))
svmfit=svm(y~.,data=dat,kernel="linear",cost=1e5)
summary(svmfit)

# And with a more realistic cost
svmfit=svm(y~.,data=dat,kernel="linear",cost=1)
summary(svmfit)
plot(svmfit,dat)

# 9.6.2 - SVM

# Generate non-linearly separable data

set.seed (1)
x=matrix(rnorm(200*2),ncol=2)
x[1:100,]=x[1:100,]+2
x[101:150,]=x[101:150,]-2
y=c(rep(1,150),rep(2,50))
dat=data.frame(x=x,y=as.factor(y))
plot(x,col=y)

# Split into training and test set, train an RBF SVM

train=sample(200,100)
svmfit=svm(y~.,data=dat[train,],kernel="radial",gamma=1,cost =1)
plot(svmfit,dat[train,])
summary(svmfit)

# Overtrain
svmfit=svm(y~.,data=dat[train,],kernel="radial",gamma=1,cost=1e5)
plot(svmfit,dat[train,])

# We can also optimize gamma using CV
set.seed(1)
tune.out=tune(svm,y~.,data=dat[train,],kernel="radial",ranges=list(cost=c(0.1,1,10,100,1000),gamma=c(0.5,1,2,3,4)))
summary(tune.out)
plot(tune.out$best.model,dat[train,])

# And test the best model
table(true=dat[-train,"y"], pred=predict(tune.out$best.model, newx=dat[-train ,]))
