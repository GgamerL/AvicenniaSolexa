## Code that has been used to simulate the human data set
## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html
## set to the path to the location of the programs "ms" and "sample_stats" on your computer
## if they are in the current directory you don't need to do anything
#written by WangZhengzhen
#This is a demo
#on 2020.1.1
## parameter declarations
numsim <- 20000#number of simulations
numloc <- 80#number of loci
mutrate <- 3.26*10^(-8)

##########################################################################
#########S2
#rep9
##model A no divergence from beginning,calculate order 1ma,2eu,3au 
head <- rep(numsim*numloc,each=numloc)
Ne <- rep(myN<-runif(numsim, 1000, 30000), each=numloc)
par.modS2A9 <- data.frame(formatC(head,digit=7),
			theta=Ne*4*mutrate*1000
			)
write.table(par.modS2A9, file="modS2A9", quote=F, row.names=F, col.names=F)
### run ms  
system(paste("ms 188 tbs -t tbs < modS2A9|", " perl calculate.pl 70 62 56 > afr_modS2A9.txt", sep=""))
#########S2
#rep9
### model B (ma+eu,au) calculate order 1ma,2eu,3au
Ne <- rep(myN <- runif(numsim,1000,30000), each = numloc)
t0 <- rep(myt0<-runif(numsim, 90001,110000),each = numloc)
head <- rep(numsim*numloc,each=numloc)
par.modS2B9 <- data.frame(formatC(head,digit=7),
			theta=Ne*4*mutrate*1000,
			t0=t0/(Ne*4))
write.table(par.modS2B9, file="modS2B9", quote=F, row.names=F, col.names=F)
### run ms  
system(paste("ms 188 tbs -t tbs -I 2 132 56 -n 1 1 -n 2 1 -ej tbs 1 2 < modS2B9 |", "perl calculate.pl 70 62 56 > afr_modS2B9.txt", sep=""))
#########S2
#rep9
##model C (ma+au,eu),calculate order 1ma, 2au,3eu
Ne <- rep(myN <- runif(numsim, 1000, 30000), each = numloc)
t0 <- rep(myt0 <-runif(numsim, 90001, 110000),each = numloc)
head <- rep(numsim*numloc,each=numloc)
par.modS2C9 <- data.frame(formatC(head,digit=7),
			theta =Ne*4*mutrate*1000,
			t0=t0/Ne/4)
write.table(par.modS2C9, file="modS2C9",quote=F,row.names=F,col.names=F)
system(paste("ms 188 tbs -t tbs -I 2 126 62 -n 1 0.9 -n 2 0.6 -ej tbs 1 2 < modS2C9 |", "perl calculate.pl 70 56 62 >afr_modS2C9.txt", sep=""))
#########S2
#rep9
#############################################3
##model D (au+eu,ma) calculate order 1au,2eu,3ma
Ne <- rep(myN <- runif(numsim, 1000, 30000), each = numloc)
t0 <- rep(myt0 <-runif(numsim, 90001, 110000),each = numloc)
head <- rep(numsim*numloc,each=numloc)
par.modS2D9 <- data.frame(formatC(head,digit=7),
			theta =Ne*4*mutrate*1000,
			t0=t0/Ne/4)
write.table(par.modS2D9, file="modS2D9",quote=F,row.names=F,col.names=F)
system(paste("ms 188 tbs -t tbs -I 2 118 70 -n 1 1.1 -n 2 0.7 -ej tbs 1 2 < modS2D9 |", "perl calculate.pl 56 62 70 >afr_modS2D9.txt", sep=""))
#########S2
#rep9
##model E ((ma,eu),au) calculate order 1ma,2eu,3au
Ne <- rep(myN <- runif(numsim, 1000, 30000), each = numloc)
t0 <- rep(myt0 <-runif(numsim, 90001, 110000),each = numloc)
t1 <- rep(myt1<-runif(numsim, 63001, 90000),each = numloc)
head <- rep(numsim*numloc,each=numloc)
par.modS2E9 <- data.frame(formatC(head,digit=7),
			theta =Ne*4*mutrate*1000,
			t1=t1/Ne/4,
			t0=t0/Ne/4)
write.table(par.modS2E9, file="modS2E9",quote=F,row.names=F,col.names=F)
system(paste("ms 188 tbs -t tbs -I 3 70 62 56 -n 1 0.7 -n 2 0.6 -n 3 1 -ej tbs 1 2 -ej tbs 2 3 < modS2E9 |", "perl calculate.pl 70 62 56 >afr_modS2E9.txt", sep=""))
#########S2
#rep9
##model F ((ma,au),eu) calculate order 1ma,2au,3eu
Ne <- rep(myN <- runif(numsim, 1000, 30000), each = numloc)
t0 <- rep(myt0 <-runif(numsim, 90001, 110000),each = numloc)
t1 <- rep(myt1<-runif(numsim, 63001, 90000),each = numloc)
head <- rep(numsim*numloc,each=numloc)
par.modS2F9 <- data.frame(formatC(head,digit=7),
			theta =Ne*4*mutrate*1000,
			t1=t1/Ne/4,
			t0=t0/Ne/4)
write.table(par.modS2F9, file="modS2F9",quote=F,row.names=F,col.names=F)
system(paste("ms 188 tbs -t tbs -I 3 70 56 62 -n 1 0.7 -n 2 1 -n 3 0.6 -ej tbs 1 2 -ej tbs 2 3 < modS2F9 |", "perl calculate.pl 70 56 62 >afr_modS2F9.txt", sep=""))
#########S2
#rep9
######################################################3
##model G((au,eu),ma) calculate order 1au,2eu,3ma
Ne <- rep(myN <- runif(numsim, 1000, 30000), each = numloc)
t0 <- rep(myt0 <-runif(numsim, 90001, 110000),each = numloc)
t1 <- rep(myt1<-runif(numsim, 63001, 90000),each = numloc)
head <- rep(numsim*numloc,each=numloc)
par.modS2G9 <- data.frame(formatC(head,digit=7),
			theta =Ne*4*mutrate*1000,
			t1=t1/Ne/4,
			t0=t0/Ne/4)
write.table(par.modS2G9, file="modS2G9",quote=F,row.names=F,col.names=F)
system(paste("ms 188 tbs -t tbs -I 3 56 62 70 -n 1 1 -n 2 0.6 -n 3 0.7 -ej tbs 1 2 -ej tbs 2 3 < modS2G9 |", "perl calculate.pl 56 62 70 >afr_modS2G9.txt", sep=""))
#########S2
#rep9
#model H (ma,eu,au) calculate order 1ma, 2eu 3au
Ne <- rep(myN <- runif(numsim, 1000, 30000), each = numloc)
t0 <- rep(myt0 <-runif(numsim, 63001, 110000),each = numloc)
head <- rep(numsim*numloc,each=numloc)
par.modS2H9 <- data.frame(formatC(head,digit=7),#
			theta = Ne*4*mutrate*1000,
			ta = t0/Ne/4,
			tb = t0/Ne/4)
write.table(par.modS2H9, file="modS2H9",quote=F,row.names=F,col.names=F)
system(paste("ms 188 tbs -t tbs -I 3 70 62 56 -n 1 0.7 -n 2 0.6 -n 3 1 -ej tbs 1 2 -ej tbs 3 2 < modS2H9 >flow.modS2H9 ", "&& perl calculate.pl 70 62 56 <flow.modS2H9 >afr_modS2H9.txt ", "&& rm flow.modS2H9", sep=""))
