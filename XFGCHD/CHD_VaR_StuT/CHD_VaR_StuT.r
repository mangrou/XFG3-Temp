############################################### 1注意数据集，是五个列向量数据，五只股票，2需要增加exceeding ratio的计算
################# VaR_studentT     ########## 3可以画一个三维surface图，x:time，y:PnL，z:0.001-0.05 的1000个alpha分位数点，两个surface一个是基于alpha的一个是基于exceeding ratio的-有时间的话可以画出多个var面于一图
################# 20160312         ########## 20160312完成gaussian copula的作图和exceeding ratio的代码
############################################
library("foreach")
library("doParallel")

ptm <- proc.time() # 计算消耗时间
### 重新对126只股票做GARCH建模，求出1）eps，2）GARCH的4个估计参数，3）sigma.t
eps1=read.csv("C:/+++haindorf/eps1.csv") # eps1是时间升序第一个点是20070104
tikName=names(eps1)[-1]

combData=data.frame()
ticker100=tikName
tempo1 = read.csv(paste("C:/+++haindorf/sp500data/", ticker100[1],".csv", sep=""), header = T, stringsAsFactors = F)
combData1 = data.frame(tempo1[, 5]) # combData1 is the first ticker close price
names(combData1)[1]<- paste(ticker100[1], sep="")
for(i in 2:126){
#i=2
tempo = read.csv(paste("C:/+++haindorf/sp500data/", ticker100[i],".csv", sep=""), header = T, stringsAsFactors = F)
tempo = data.frame(tempo[, 5])
names(tempo)[1]<- paste(ticker100[i], sep="")
combData1 = data.frame(combData1, tempo)
} # combData1, construct data set with 126 ticker, matrix 2025*126 

combData1 = data.frame(read.csv(paste("C:/+++haindorf/sp500data/", "MMM.csv", sep=""), header = T, stringsAsFactors = F)[, 1], combData1) # add date uptrend in first column
names(combData1)[1]<- paste("date", sep="")
write.csv(combData1, file = paste("sp126.csv", sep=""))         ######################################################### sp126.csv是price时间序列第一点是20070103
sp.126=read.csv("sp126.csv", header=T)[,-c(1,2)] # sp.126 and sp126.csv, ticker matrix, 2025*126, without date
#sp.126=sp.126[,1:5]


### 求出log-return升序时间序列
p1=sp.126[-1,]
p2=sp.126[-length(sp.126[,1]),]
l.r=100*log(p1/p2) # l.r and logReturn126.csv, log-return, 2024*126, uptrend series
write.csv(l.r, file = paste("logReturn126.csv", sep=""))        ############################################### logReturn126.csv是log-return时间序列第一点是20070104

###The fit of a GARCH(1,1) model to the sample of log returns 估计GARCH模型
r=read.csv("logReturn126.csv", header=T)
r1=r[, -c(1)]
r2=r1 # r2, log return of 126 tickers, 2024*126
###########################################################################################################################################################################################################
### 合成一个5dim的数据集  r2 = data.frame(cbind(AVB, EQR, TXN, ADI, LLY))    ##############################################################################################################################
###########################################################################################################################################################################################################
#################### !!!!!!!!!!!!!!!!!!!!!!only first 2 tickers are employed for code test, i.e 2 dim, 300 window-width, 500 steps forecasting, 1000 simulation.
dims=5 #*************************************************************** 股票数目
attach(r2)
r2 = data.frame(cbind(AVB, EQR, TXN, ADI, LLY)) 
cor(r2, method="kendall", use="pairwise") 

############################# 窗宽300，循环计算fit GARCH
############################# 窗宽300，循环计算
############################# 窗宽300，循环计算


#install.packages("copula")
library(copula)
M=1000 #*************************************************************** 蒙特卡洛模拟次数
backtestNr=1000 #******************************************************** backtestNr
VaR=matrix(NA, backtestNr, 4)
slidingWindowLength=300 #********************************************** window width
install.packages("fGarch")
library(fGarch)
eps=matrix(NA, backtestNr, dim(r2)[2])
colnames(eps)=colnames(r2)
eps1=data.frame( ) # matrix 500*2, standard eps, eps_t-1
dat=eps1
paraMat=list() # matrix 500*2*4
para=matrix(NA, dims, 4) # matrix 2*4 parameter matrix of GARCH, para_t-1
sigma=matrix(NA, backtestNr, dims) # matrix 500*2, sigma from GARCH, sigma_t-1
# for(i in 1:dims){ ##################################################################################################### fit GARCH模型
    # fit = garchFit(~garch(1, 1), data = r2[,i], trace = F)
    # eps1[,i] = fit@residuals / fit@sigma.t
	# para[i,]= fit@fit$coef # 估计出的garch的参数                         
    # sigma[, i] = fit@sigma.t                                
	# write.csv(eps1, file = paste("eps1.csv", sep=""))

# } 
# write.csv(eps1, file = paste("eps1.csv", sep=""))              ######################################################### eps1.csv是标准化余差standardized residuals/para是估计的GARCH模型的4个参数/sigma是标准差

#############################做并行计算来估计GARCH参数-完成
#############################做并行计算来估计GARCH参数
#############################做并行计算来估计GARCH参数
#install.packages("foreach")
#install.packages("doParallel")

########### foreach+doParallel 并行计算

library(foreach)
library(doParallel)

############## 并行计算开始

cl <- makeCluster(37) # nr of kernels
registerDoParallel(cl) # 设定4核运算
getDoParWorkers()

#rhoSpace=as.list(rhoSpace)
### data matrix
	dat = r2 # dat, log return, 2024*2
	datMat=list()

for(i in 1:backtestNr){
	datMat[[i]]=dat[c(i:(i+(slidingWindowLength-1))), 1:dims] # plug blocked matrix of eps1 into list structure, datMat has 500 length from first 300 pts to slide 500 steps
} # datMat, log return matrix, 500*300*2

#dMat=datMat[[1]]
### objective function
lengthPara=4 # 估计的GARCH 参数个数
lengthEps1= slidingWindowLength# EPS和sigma个数
lengthSigma= slidingWindowLength
objFun=function(dMat){
# dMat=datMat[[1]]
	paraComb=list() # paraComb matrix dims*Nr of Para, 604*dims
    library(fGarch)
	for(i in 1:dims){
	#i=1
		fit = garchFit(~garch(1, 1), data = dMat[[i]], trace = F) # dMat[[i]] input first ticker's log-return
		eps1.loop= fit@residuals / fit@sigma.t
		para.loop= fit@fit$coef # 估计出的garch的参数                         
		sigma.loop = fit@sigma.t      
		lengthPara=length(para.loop)
		paraComb[[i]]=c(para.loop, sigma.loop, eps1.loop) # paraComb[[i]] matrix 604*1, first 4 are GARCH para, then first 300 sigma, then 300 eps1
	}
	return(paraComb)
}
   

resultD <- foreach(dMat=datMat) %dopar% objFun(dMat) # resultD; garch para list, 500*604*2, 500 lists for everyone contain 2 lists 
stopCluster(cl)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VERY IMPORTANT	
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VERY IMPORTANT
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VERY IMPORTANT
resultD  # resultD, para of GARCH为一个list, 500*dim*Nr of para, parameters of GARCH
#refer resultD[[5]][[1]], 第5期第1个ticker的参数列表， GARCH参数+300个sigma+300个eps1，其中resultD是一个list,长度500，每个element是一个604*2的矩阵
############## 并行计算结束65/30



write.csv(cbind(resultD), file = "garchPara500.604.dims.csv") # 把所有obj的结果写到表上



###
#2015-1-29下午
#############################做并行计算来估计copula参数-完成
#############################做并行计算来估计copula参数
#############################做并行计算来估计copula参数
#install.packages("foreach")
#install.packages("doParallel")

########### foreach+doParallel 并行计算

library(foreach)
library(doParallel)

############## 并行计算开始

cl <- makeCluster(37)
registerDoParallel(cl) # 设定4核运算
getDoParWorkers()

#rhoSpace=as.list(rhoSpace)
lengthPara=4 # 估计的GARCH 参数个数
lengthEps1= slidingWindowLength# EPS和sigma个数
lengthSigma= slidingWindowLength
totalLength=lengthPara+lengthEps1+lengthSigma
### data matrix
epsComb=list()
for(i in 1:backtestNr){
	#i=1
    k=pobs(as.data.frame(resultD[[i]]))
	epsComb[[i]]=k[-(1:(lengthPara+lengthEps1)), ] # epsComb, standard esp, 500*300*dims

}

datMatEps=epsComb # datMatEps, standard residuals, list, matrix 500*300*2
### objective function
objFun=function(dMat){
#dMat=datMatEps[[1]]
    library(copula)
    t.cop <- tCopula(0.3, dim = dims, dispstr = "ex", df = 2)

    
    
	#normal.cop = normalCopula(c(.91), dim = dims, dispstr = "ex")
	fit.mln=fitCopula(t.cop, dMat, method="mpl")
    para.normal=summary(fit.mln)$coefficient[1:2] # use eps1 data to estimate copula para for the 301 th eps' forecasting. i.e. Cop_300, Cop_301,..., Cop_799.
	return(para.normal)
}
   

resultDcopPara <- foreach(dMat=datMatEps,.combine='rbind') %dopar% objFun(dMat) # resultDcopPara， 500*1，估计出500个copula参数
stopCluster(cl)
resultDcopPara  # copula's 500 para

############## 并行计算结束65/30



write.csv(cbind(resultDcopPara), file = "C:/+++haindorf/resultDcopPara.csv") # 把所有obj的结果写到表上 


#??????????????????????????????????????????????????????????????????????????????????????????????????????????2015-1-29下午
#??????????????????????????????????????????????????????????????????????????????????????????????????????????2015-1-29下午
#??????????????????????????????????????????????????????????????????????????????????????????????????????????2015-1-29下午

#############################带入所有参数GARCH参数，COPULA参数计算VaR-未完成-仍然有问题
#############################带入所有参数GARCH参数，COPULA参数计算VaR
#############################带入所有参数GARCH参数，COPULA参数计算VaR
#r2 = data.frame(cbind(AVB, EQR, TXN, ADI, LLY)) 
###
spread.real.0=read.csv("sp126.csv", header=T)[,-c(1:2)]
attach(spread.real.0)
### real loss of portfolio
spread.real =data.frame(cbind(AVB, EQR, TXN, ADI, LLY))# 读取sp126数据2025*126,只有2个tickers, s, close price of 2 tickers, 2025*dims
#### portfolio的实际损失-real P&L
S=rowSums(spread.real)
S1=S[-1]
S2=S[-length(S)]
S3=S2-S1 # S3, portfolio value, 2024*1
L.real=S3[301:(300+backtestNr)]
    
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 129工作把此处重写
### 计算循环sliding window的VaR, 基于copula模拟的损失-simulated P&L
#install.packages("foreach")
#install.packages("doParallel")

########### foreach+doParallel 并行计算

library(foreach)
library(doParallel)

############## 并行计算开始

cl <- makeCluster(37)
registerDoParallel(cl) # 设定4核运算
getDoParWorkers()

#rhoSpace=as.list(rhoSpace)
lengthPara=4 # 估计的GARCH 参数个数
lengthEps1= slidingWindowLength# EPS和sigma个数
lengthSigma= slidingWindowLength
totalLength=lengthPara+lengthEps1+lengthSigma

### objective function
para.vec = read.csv("C:/+++haindorf/resultDcopPara.csv", header=T)[, -1] # 读取并行计算得到的估计参数数据，该数据由并行计算求得，下一步带入copula, para.vec, parameter vector, 500*1
VaR= matrix(NA, backtestNr, 4) # VaR, 500*4
#LRetSim=matrix(NA, 500, 2)

datMatIndex=c(1:backtestNr)
VaRt_store=matrix(777, backtestNr,4)
objFunVaR=function(dMat){
i=dMat
#i=1
library(copula)
	#normal.cop = normalCopula(c(.91), dim = dims, dispstr = "ex")         ############################################################################## construct normal cop in 126dim
	#dat = eps1
	#dat <- pobs(dat)
    #fit.mln <- fitCopula(normal.cop, dat[c(i:(i+499)), 1:dims], method="mpl")############################################################################### 估计copula参数，GAUSSIAN COPULA使用单个参数
    #fit.mln
    #para.normal = summary(fit.mln)$coefficients[1]
	para.normal=as.numeric(para.vec[i,])
    t.cop.2 <- tCopula(para.normal[1], dim = dims, dispstr = "ex", df = para.normal[2])
    
	#normal.cop.2 = normalCopula(para.normal, dim = dims, dispstr = "ex")         ####################################################################### construct normal cop in 126dim

	u = rCopula(M, t.cop.2) # u抽取的copula的随机数，matrix 1000*2, 1000*dims
	u = qnorm(u)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VERY IMPORTANT	
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VERY IMPORTANT
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VERY IMPORTANT resultD[[5]][[1]]
	k=as.data.frame(resultD[[i]]) # k, parameter of GARCH(para+sigma+epsilon), matrix 604*2
	k1=k[c(1:lengthPara,lengthPara+lengthEps1,totalLength),]# k1, 4 GARCH para+sig.301+eps.301 matrix 6*2
	#h = sqrt(params[sel.pairs,2] + params[sel.pairs,3] * sigma.t[dim(sigma.t)[1], sel.pairs]^2 * eps[dim(eps)[1], sel.pairs]^2 + params[sel.pairs,4] * sigma.t[dim(sigma.t)[1], sel.pairs]^2)    #############这个地方极其重要是使用GARCH对logreturn建模。                                                                                                                                                                       
   
    
	h= sqrt(k1[2, ] + k1[3, ] * k1[lengthPara+1, ]^2 * k1[lengthPara+2, ]^2 + k1[4, ] * k1[lengthPara+1, ]^2) # sigma_t+1的估计,matrix 2*1 
	k2=matrix(NA, M, dims)
	k3=matrix(NA, M, dims)
	k2[1,]=as.matrix(k1[1, ]) # mu
	mu.t=k2[rep(1, M), ] # mu, matrix 1000*2
	k3[1,]=as.matrix(k1[lengthPara+1, ]) # sigma 
	sig.t=k3[rep(1, M), ] # sigma, 1000*2
	
	R= mu.t+sig.t*u # forecasting of log return at t=302
	#LRetSim[i, ]=colMeans(R)
	#R = qnorm(u) # 此处得到的是simulation出来的log-return
    spread.real.loop=read.csv("sp126.csv")[,-c(1:2)]
attach(spread.real.loop)
### real loss of portfolio
spread =data.frame(cbind(AVB, EQR, TXN, ADI, LLY))
    
    
	#spread = read.csv("sp126.csv", header=T)[,-c(1:2)][,1:dims] # 读取spread数据,只有2个tickers
	s1=spread[i+(slidingWindowLength), ]##################################################### s1, spread of 2 ticker at t301, 1*nr of ticker
	
	st=data.frame(matrix(NA, M, dims)) # st, repeat of spread at 301 for 1000 times, 1000*dims
	#ma <- matrix(1:6, nrow  = 2)
#rma <- ma[rep(1:2, c(2,3)),]
    
	st[1,]=s1
	st=st[rep(1, M),]
	#for(j in 1:1000) st[j, ]=s1
	
	L.sim=matrix(NA, M, 1) # simulated value of portfolio at t302, portfolioValue_t302=P_t301*(exp(0.01*logreturn_t302)-1), matrix 1000*1, 
    #for(m in 1:M){
	#	ExpMinus1=exp(R[m, ])-1
	#	L.sim[m]=sum(s1*ExpMinus1)
	#}
	#m=1
    L.sim=rowSums(st*(exp(0.01*R)-1)) ################ 要明白copula模拟的是standard residual不是log return,而log-return需要用GARCH建模构建
		VaRt1005=quantile(L.sim, .05)
		VaRt1001=quantile(L.sim, .01)
		VaRt10005=quantile(L.sim, .005)
		VaRt10001=quantile(L.sim, .001)
		VaRt=c(VaRt1005, VaRt1001, VaRt10005, VaRt10001)


	return(VaRt) ##########################################################

	##
}

# for(dMat in 1:backtestNr){
# VaRt_store[dMat,]=objFunVaR(dMat)
# }


resultVaR <- foreach(dMat=datMatIndex,.combine='rbind') %dopar% objFunVaR(dMat) # resultDcopPara， 500*1，估计出500个copula参数
stopCluster(cl)
head(resultVaR)  # resultVaR, 500*4
VaR=resultVaR
#??????????????????????????????????????????????????????????????????????????????????????????????????????????2015-1-29下午
#??????????????????????????????????????????????????????????????????????????????????????????????????????????2015-1-29下午
#??????????????????????????????????????????????????????????????????????????????????????????????????????????2015-1-29下午
# ### 
# ### real loss of portfolio
# s=read.csv("sp100-128.csv", header=T)[,-1][,1:2]# 读取spread数据,只有2个tickers
# #### portfolio的实际损失-real P&L
# spread.real = read.csv("sp100-128.csv", header=T)[,-1][,1:2]# 读取spread数据,只有2个tickers
# S=rowSums(spread.real)
# S1=S[-1]
# S2=S[-length(S)]
# S3=S2-S1
# L.real=S3[300:(300+backtestNr-1)]

### combine L.real and VaR
resultComb=data.frame(L.real, VaR)
head(resultComb)
###


### exceeding ratio 
##  resultComb[,1]是real portfolio loss = l_t, resultComb[,2:5]是VaR_t值

Exceeding_Ratio=numeric(4)
for(alpha in 2:5){
nullVector=rep(0,backtestNr)
for(i in 1:length(resultComb[,alpha])){
  if(resultComb[,1][i]<resultComb[,alpha][i]){
    nullVector[i]=1
    }else{
    nullVector[i]=0
    }
}
Exceeding_Ratio[alpha]=sum(nullVector)/backtestNr
}
Exceeding_Ratio


# > Exceeding_Ratio
# [1] 0.000 0.034 0.013 0.011 0.006


###########################################################################################################################################################################################################
### 画图 ##################################################################################################################################################################################################
###########################################################################################################################################################################################################
#############################画图-未完成-仍然有问题
#############################画图-未完成
#############################画图-未完成
nameOfFigure="studentTtry"

#############################画图-未完成
alpha=2
filenamei <- paste("C:/xfg_figure/",nameOfFigure,alpha,".png",sep="")
png(filename=filenamei)

ptForPlot=c(min(resultComb[,1]), min(resultComb[,2]), min(resultComb[,3]), min(resultComb[,4]), min(resultComb[,5]))
lowPt=min(ptForPlot)
upPt=max(resultComb[,1])
Portfolio_Value=seq(lowPt,upPt, length.out=length(resultComb[,1]))
Time_Index=1:length(seq(lowPt,upPt, length.out=length(resultComb[,1])))
plot(Time_Index, Portfolio_Value, col = "white", pch = 19, cex = 2.5, xlab="Time Index", ylab="Profit and Loss of Portfolio")
lines(resultComb[,alpha], col = "gray", lwd=6)
#points(resultComb[,1], col = "black", pch = 19, cex = 1.5)
for(i in 1:length(resultComb[,alpha])){
  if(resultComb[,1][i]<resultComb[,alpha][i]){
    points(i,lowPt+1, col = "black", pch = 17, cex = 1.5) # solid triangle
    points(i,resultComb[,1][i], col = "black", pch = 3, cex = 2.5, lwd=1) # crossade
    points(i,resultComb[,1][i], col = "black", pch = 5, cex = 1.5, lwd=1) # crossade
    }else{
    points(i,resultComb[,1][i], col = "black", pch = 19, cex = 1) # solid points
    
    }
}


dev.off()

#############################画图-未完成
alpha=3
filenamei <- paste("C:/xfg_figure/",nameOfFigure,alpha,".png",sep="")
png(filename=filenamei)

ptForPlot=c(min(resultComb[,1]), min(resultComb[,2]), min(resultComb[,3]), min(resultComb[,4]), min(resultComb[,5]))
lowPt=min(ptForPlot)
upPt=max(resultComb[,1])
Portfolio_Value=seq(lowPt,upPt, length.out=length(resultComb[,1]))
Time_Index=1:length(seq(lowPt,upPt, length.out=length(resultComb[,1])))
plot(Time_Index, Portfolio_Value, col = "white", pch = 19, cex = 0.5, xlab="Time Index", ylab="Profit and Loss of Portfolio")
lines(resultComb[,alpha], col = "gray", lwd=6)
#points(resultComb[,1], col = "black", pch = 19, cex = 1.5)
for(i in 1:length(resultComb[,alpha])){
  if(resultComb[,1][i]<resultComb[,alpha][i]){
    points(i,lowPt+1, col = "black", pch = 17, cex = 1.5) # solid triangle
    points(i,resultComb[,1][i], col = "black", pch = 3, cex = 2.5, lwd=1) # crossade
    points(i,resultComb[,1][i], col = "black", pch = 5, cex = 1.5, lwd=1) # crossade
    }else{
    points(i,resultComb[,1][i], col = "black", pch = 19, cex = 1) # solid points
    
    }
}


dev.off()


#############################画图-未完成
alpha=4
filenamei <- paste("C:/xfg_figure/",nameOfFigure,alpha,".png",sep="")
png(filename=filenamei)

ptForPlot=c(min(resultComb[,1]), min(resultComb[,2]), min(resultComb[,3]), min(resultComb[,4]), min(resultComb[,5]))
lowPt=min(ptForPlot)
upPt=max(resultComb[,1])
Portfolio_Value=seq(lowPt,upPt, length.out=length(resultComb[,1]))
Time_Index=1:length(seq(lowPt,upPt, length.out=length(resultComb[,1])))
plot(Time_Index, Portfolio_Value, col = "white", pch = 19, cex = 0.5, xlab="Time Index", ylab="Profit and Loss of Portfolio")
lines(resultComb[,alpha], col = "gray", lwd=6)
#points(resultComb[,1], col = "black", pch = 19, cex = 1.5)
for(i in 1:length(resultComb[,alpha])){
  if(resultComb[,1][i]<resultComb[,alpha][i]){
    points(i,lowPt+1, col = "black", pch = 17, cex = 1.5) # solid triangle
    points(i,resultComb[,1][i], col = "black", pch = 3, cex = 2.5, lwd=1) # crossade
    points(i,resultComb[,1][i], col = "black", pch = 5, cex = 1.5, lwd=1) # crossade
    }else{
    points(i,resultComb[,1][i], col = "black", pch = 19, cex = 1) # solid points
    
    }
}


dev.off()

#############################画图-未完成
alpha=5
filenamei <- paste("C:/xfg_figure/",nameOfFigure,alpha,".png",sep="")
png(filename=filenamei)

ptForPlot=c(min(resultComb[,1]), min(resultComb[,2]), min(resultComb[,3]), min(resultComb[,4]), min(resultComb[,5]))
lowPt=min(ptForPlot)
upPt=max(resultComb[,1])
Portfolio_Value=seq(lowPt,upPt, length.out=length(resultComb[,1]))
Time_Index=1:length(seq(lowPt,upPt, length.out=length(resultComb[,1])))
plot(Time_Index, Portfolio_Value, col = "white", pch = 19, cex = 0.5, xlab="Time Index", ylab="Profit and Loss of Portfolio")
lines(resultComb[,alpha], col = "gray", lwd=6)
#points(resultComb[,1], col = "black", pch = 19, cex = 1.5)
for(i in 1:length(resultComb[,alpha])){
  if(resultComb[,1][i]<resultComb[,alpha][i]){
    points(i,lowPt+1, col = "black", pch = 17, cex = 1.5) # solid triangle
    points(i,resultComb[,1][i], col = "black", pch = 3, cex = 2.5, lwd=1) # crossade
    points(i,resultComb[,1][i], col = "black", pch = 5, cex = 1.5, lwd=1) # crossade
    }else{
    points(i,resultComb[,1][i], col = "black", pch = 19, cex = 1) # solid points
    
    }
}


dev.off()





proc.time() - ptm # 7m21s
##






