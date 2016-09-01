setwd("")
tmpdata1 = read.csv("lambda_mean_206vars_2016-07-15.csv", sep=",")
tmpdata2 = read.csv("200_firms_returns_and_scaled_macro_2016-08-18.csv", sep=",")

lambda       = tmpdata1[, 2]
dates        = tmpdata1[, 1]
vix          = tmpdata2[128:2400, 202]
vix.dates    = tmpdata2[128:2400, 1]

par(mfrow = c(1,1))
vix.norm     = (vix - min(vix))/(max(vix) - min(vix))
lambda.norm  = (lambda - min(lambda))/(max(lambda) - min(lambda))
par(mar = c(5, 6, 1, 1))
plot(vix.norm, type = "l",  col =  "darkblue", axes = FALSE, 
     xlab = "Year", frame = TRUE, cex.main = 1.5, ylab = expression(paste("Average ", lambda)),
     cex.lab = 2)
at.tmp = c(grep("2008", dates)[1], grep("2009", dates)[1], grep("2010", dates)[1], grep("2011", dates)[1], 
           grep("2012", dates)[1], grep("2013", dates)[1], grep("2014", dates)[1], grep("2015", dates)[1],
           grep("2016", dates)[1])
axis(1, cex.axis = 1.5, labels = c(2008:2016), at = at.tmp)
axis(2, cex.axis = 1.5)
lines(lambda.norm, col = "red3")


