library(ggplot2)
dat2<-read.csv("for Correlationship.csv",header = T)
dat2_ICIM<- dat2[(1:9),]
row.names(dat2_ICIM)<-dat2_ICIM$ï..
dat2_ICIM<-dat2_ICIM[,-1]
class(dat2_ICIM)

dat2_ICIM<-as.data.frame(dat2_ICIM)

#plot with ggplot dataframe
dat_cor<-read.csv("tpm_filtered_ICM.csv",row.names = 1)
dat_cor<-as.data.frame(t(dat_cor))
qplot(CD8A, CXCL9, data=dat_cor, geom=c("point","smooth"), 
      method="lm", formula=y~x, 
      main="Regression of CXCL9 on CD8A", 
      xlab="CD8A", ylab="CXCL9")
ggsave("CD8A_CXCL9.pdf",width = 8,height = 5,path = "HighvsLow/")

# Data generation data as matrix
set.seed(1)
x <- as.numeric(dat_cor$CD8A)
          
y <- as.numeric(dat_cor$IFNG)

# Creating the plot
plot(x, y, pch = 19, col = "lightblue",xlab="CD8A", ylab="IFNG")
# Regression line
abline(lm(y ~ x), col = "red", lwd = 1)
# Pearson correlation
cor.test(x,y,method="pearson",use="complete.obs")

ggsave("correlationship_STAT1&CD8A.pdf",width = 8,height=5,path="exported pics/")


# Data generation data as matrix
dat2_ICIM<-as.matrix(dat2_ICIM)
set.seed(1)
x <- dat2_ICIM[,5]
y <- dat2_ICIM[,2]

# Creating the plot
plot(x, y, pch = 19, col = "lightblue",xlab="IFI30", ylab="CXCL9")

# Regression line
abline(lm(y ~ x), col = "red", lwd = 1)
# Pearson correlation
cor.test(x,y,method="pearson",use="complete.obs")

ggsave("correlationship_STAT1&CD8A.pdf",width = 8,height=5,path="exported pics/")
