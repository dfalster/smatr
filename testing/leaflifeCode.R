#first load smatr from CRAN into R. For R version 2 or above, this can be done using the menus: "Packages...Install Packages" then "Packages...Load Package"
require(smatr)

setwd('F:\\_resources\\C\\SMATR\\test')
#source('F:\\_resources\\C\\_projects\\SMATR\\R\\shift.com.R')
dat=read.table("leaflife.txt",header=T)

#------------------------------------------------
#Without Measurement Error
#checks of slope.test:
ft=slope.test(log10(longev), log10(lma), test.value=1.3, data=dat, method="SMA", intercept=TRUE)
c(ft$b,ft$ci,ft$test.value, ft$r, ft$p)
ft=slope.test(log10(longev), log10(lma), test.value=1.3, data=dat, method="MA", intercept=TRUE )
c(ft$b,ft$ci,ft$test.value, ft$r, ft$p)
ft=slope.test(log10(longev), log10(lma), test.value=0.1, data=dat, method="SMA", intercept=FALSE)
c(ft$b,ft$ci,ft$test.value, ft$r, ft$p)
ft=slope.test(log10(longev), log10(lma), test.value=0.1, data=dat, method="MA", intercept=FALSE)
c(ft$b,ft$ci,ft$test.value, ft$r, ft$p)

#checks of elev.test:
ft=elev.test(log10(longev), log10(lma), -2.5,data = dat,method="MA")
c(ft$a,ft$a.ci,ft$test.value, ft$t, ft$p)
#note mistake in elev.test

#checks of line.cis:
ft=line.cis(log10(longev), log10(lma), data = dat,method="MA")

#checks of slope.com:
ft=slope.com(log10(longev), log10(lma), site, data = dat,method="SMA")
c(ft$b,ft$ci,ft$LR,ft$p)
ft=slope.com(log10(longev), log10(lma), site, data = dat,method="MA")
c(ft$b,ft$ci,ft$LR,ft$p)
ft=slope.com(log10(longev), log10(lma), site, data = dat, method="SMA", intercept=F)
c(ft$b,ft$ci,ft$LR,ft$p)
ft=slope.com(log10(longev), log10(lma), site, data = dat, method="MA", intercept=F)
c(ft$b,ft$ci,ft$LR,ft$p)

#checks of elev.com:
ft=elev.com(log10(longev), log10(lma), site, data = dat,method="SMA")
c(ft$a,ft$ci,ft$stat,ft$p,ft$as)
ft=elev.com(log10(longev), log10(lma), site, data = dat, method="MA")
c(ft$a,ft$ci,ft$stat,ft$p,ft$as)

#checks of shift.com:
source('F:\\_resources\\C\\_projects\\SMATR\\R\\shift.com.R')
source('F:\\_resources\\C\\_projects\\SMATR\\R\\slope.com.R')

ft=shift.com(log10(longev), log10(lma), site, data = dat,method="SMA")
c(ft$stat,ft$p,ft$f.mean)
ft=shift.com(log10(longev), log10(lma), site, data = dat,method="MA")
c(ft$stat,ft$p,ft$f.mean)
ft=shift.com(log10(longev), log10(lma), site, data = dat, method="SMA", intercept=F)
c(ft$stat,ft$p,ft$f.mean)
ft=shift.com(log10(longev), log10(lma), site, data = dat, method="MA", intercept=F)
c(ft$stat,ft$p,ft$f.mean)

#------------------------------------------------
#With Measurement Error
#measurement error
VV=c(0.002,0,0,0.001); 
dim(VV)=c(2,2)

#checks of slope.test:
ft=slope.test(log10(longev), log10(lma), test.value=1.3, data=dat, method="SMA", intercept=TRUE, V=VV)
c(ft$b,ft$ci,ft$test.value, ft$r, ft$p)
ft=slope.test(log10(longev), log10(lma), test.value=1.3, data=dat, method="MA", intercept=TRUE, V=VV)
c(ft$b,ft$ci,ft$test.value, ft$r, ft$p)
ft=slope.test(log10(longev), log10(lma), test.value=0.1, data=dat, method="SMA", intercept=FALSE, V=VV)
c(ft$b,ft$ci,ft$test.value, ft$r, ft$p)
ft=slope.test(log10(longev), log10(lma), test.value=0.1, data=dat, method="MA", intercept=FALSE, V=VV)
c(ft$b,ft$ci,ft$test.value, ft$r, ft$p)

#checks of elev.test:
ft=elev.test(log10(longev), log10(lma), -2.5,data = dat,method="MA", V=VV)
c(ft$a,ft$a.ci,ft$test.value, ft$t, ft$p)
#note mistake in elev.test

#checks of line.cis:
line.cis(log10(longev), log10(lma), data = dat,method="SMA", V=VV)
line.cis(log10(longev), log10(lma), data = dat,method="MA", V=VV)
line.cis(log10(longev), log10(lma), data = dat,method="SMA", V=VV,intercept=FALSE)
line.cis(log10(longev), log10(lma), data = dat,method="MA", V=VV,intercept=FALSE)

#measurement error
V=c(0.002,0,0,0.001); VV=c(V,V,V,V)
dim(VV)=c(2,2,4)

#checks of slope.com:
ft=slope.com(log10(longev), log10(lma), site, data = dat,method="SMA", V=VV)
c(ft$b,ft$ci,ft$LR,ft$p)
ft=slope.com(log10(longev), log10(lma), site, data = dat,method="MA", V=VV)
c(ft$b,ft$ci,ft$LR,ft$p)
ft=slope.com(log10(longev), log10(lma), site, data = dat, method="SMA", intercept=F, V=VV)
c(ft$b,ft$ci,ft$LR,ft$p)
ft=slope.com(log10(longev), log10(lma), site, data = dat, method="MA", intercept=F, V=VV)
c(ft$b,ft$ci,ft$LR,ft$p)

#checks of elev.com:
ft=elev.com(log10(longev), log10(lma), site, data = dat,method="SMA", V=VV)
c(ft$a,ft$ci,ft$stat,ft$p,ft$as)
ft=elev.com(log10(longev), log10(lma), site, data = dat, method="MA", V=VV)
c(ft$a,ft$ci,ft$stat,ft$p,ft$as)

#checks of shift.com:

ft=shift.com(log10(longev), log10(lma), site, data = dat,method="SMA", V=VV)
c(ft$stat,ft$p,ft$f.mean)
ft=shift.com(log10(longev), log10(lma), site, data = dat,method="MA", V=VV)
c(ft$stat,ft$p,ft$f.mean)
ft=shift.com(log10(longev), log10(lma), site, data = dat, method="SMA", intercept=F, V=VV)
c(ft$stat,ft$p,ft$f.mean)
ft=shift.com(log10(longev), log10(lma), site, data = dat, method="MA", intercept=F, V=VV)
c(ft$stat,ft$p,ft$f.mean)





# to do measurement error model checks of slope.com:
V=c(0.002,0,0,0.001)
VV=c(V,V,V,V)
dim(VV)=c(2,2,4)
ft=slope.com(log10(longev), log10(lma), site, data = dat,method="MA",intercept=T,V=VV)
c(ft$b,ft$ci,ft$LR,ft$p)
ft=slope.com(log10(longev), log10(lma), site, data = dat,method="SMA",intercept=T,V=VV)
c(ft$b,ft$ci,ft$LR,ft$p)
ft=slope.com(log10(longev), log10(lma), site, data = dat,method="SMA",intercept=F,V=VV)
c(ft$b,ft$ci,ft$LR,ft$p)
ft=slope.com(log10(longev), log10(lma), site, data = dat,method="MA",intercept=F,V=VV)
c(ft$b,ft$ci,ft$LR,ft$p)
