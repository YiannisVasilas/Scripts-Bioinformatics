# Intall Packages
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
install.packages("ggpubr")
# Load the packages
library("ggpubr")
library("ggpubr")
#Pearson and Spearman correlation 
# x and y have to be the same leght fix it!
cor(x, y, method = c("pearson","spearman"))
cor.test(x, y, method=c("pearson","spearman"))
# Load csv file
Data = read.csv(file.choose())
#Make a scatter plot
library("ggpubr")
ggscatter(Data, x = "", y = "", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#Pearson
res <- cor.test(my_data$wt, my_data$mpg, 
                method = "pearson")
res
#Spearman
res2 <-cor.test(my_data$wt, my_data$mpg,  method = "spearman")
res2