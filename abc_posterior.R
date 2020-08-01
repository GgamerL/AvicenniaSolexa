library(abc)
setwd("~/SOLEXA/A_marina/simulation")
dfrm=rep(c("mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"),c(100000,100000,100000,100000,100000,100000,100000,100000))
dfrm1 <- read.table("sumout_8mod4rplot.txt",header=T)
dfrm2 <- read.table("S2R9.obs",header=T)

################distinguish each model

pdf(paste("modsel_8modConfusionMatrix",".pdf",sep = ""),width=15,height=8)
par(mfrow=c(2,3))
cv.modsel1 <- cv4postpr(dfrm, dfrm1, nval=10, tol=.01, method="rejection")
mod_confuse1 <- summary(cv.modsel1) 
plot(cv.modsel1, names.arg=c("mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"),col=c("azure4","steelblue","coral","brown","green4","red","black","blue"))

cv.modsel2 <- cv4postpr(dfrm, dfrm1, nval=10, tol=.1, method="rejection")
mod_confuse2 <- summary(cv.modsel2)
plot(cv.modsel2, names.arg=c("mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"),col=c("azure4","steelblue","coral","brown","green4","red","black","blue"))

cv.modsel3 <- cv4postpr(dfrm, dfrm1, nval=10, tol=0.5, method="rejection")
mod_confuse3 <- summary(cv.modsel3)
plot(cv.modsel3, names.arg=c("mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"),col=c("azure4","steelblue","coral","brown","green4","red","black","blue"))

cv.modsel4 <- cv4postpr(dfrm,dfrm1,nval =10,tol=0.01,method = "neuralnet")
mod_confuse4 <- summary(cv.modsel4)
plot(cv.modsel4, names.arg=c("mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"),col=c("azure4","steelblue","coral","brown","green4","red","black","blue"))

cv.modsel5 <- cv4postpr(dfrm,dfrm1,nval =10,tol=0.1,method ="neuralnet")
mod_confuse5 <- summary(cv.modsel5)
plot(cv.modsel5, names.arg=c("mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"),col=c("azure4","steelblue","coral","brown","green4","red","black","blue"))

dev.off()
library(abc)

##################################

ptr.rej1 <- postpr(dfrm2,dfrm,dfrm1,tol = 0.01, method ="rejection")
summary(ptr.rej1)
ptr.rej2 <- postpr(dfrm2,dfrm,dfrm1,tol = 0.1, method ="rejection")
summary(ptr.rej2)
ptr.rej3 <- postpr(dfrm2,dfrm,dfrm1,tol = 0.8, method ="rejection")
summary(ptr.rej3)
ptr.mnlog1 <- postpr(dfrm2,dfrm,dfrm1,tol = 0.01,method="mnlogistic")
summary(ptr.mnlog1)
ptr.mnlog2 <- postpr(dfrm2,dfrm,dfrm1,tol = 0.1,method="mnlogistic")
summary(ptr.mnlog2)

ptr.nnet_tol0.8 <- postpr(dfrm2,dfrm,dfrm1,tol = 0.8,method = "neuralnet",maxit = 5000)
summary(ptr.nnet_tol0.8)
#######save environmen
save.image("sim_8mods.RData")
##################################
