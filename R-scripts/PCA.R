#Install and load data
library(ggplot2)
Data = read.csv(file.choose())

myPr <- prcomp(iris[, -5])
myPr <- prcomp(iris[, -5], scale = TRUE)
plot(iris$Sepal.Length, iris$Sepal.Width)
plot(scale(iris$Sepal.Length), scale(iris$Sepal.Width))
plot((iris$Sepal.Length - mean(iris$Sepal.Length)) / sd(iris$Sepal.Length))
myPr
summary(myPr)
plot(myPr, type = "l")
biplot(myPr)
biplot(myPr, scale = 0)
str(myPr)
myPr$x
iris2 <- cbind(iris, myPr$x)


ggplot(iris2, aes(PC1, PC2, col = Species, fill = Species)) +
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")

cor(iris[, -5], iris2[, 6:9])