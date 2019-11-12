require(mgcv)
#require(lme4)
#require(xtable)
require(itsadug)
#require(plotfunctions)
#require(devtools)
#require(RCurl)
#require(akima)
#require(RColorBrewer)
#require(colorRamps)
require(dichromat)

#source(textConnection(getURL("https://gist.github.com/mages/5339689/raw/576263b8f0550125b61f4ddba127f5aa00fa2014/add.alpha.R")))
#source(textConnection(getURL("https://raw.githubusercontent.com/soskuthy/gamm_intro/master/gamm_hacks.r")))

load('GAMM_data.Rda')


# Monophthongs: Differences between /i:/ and /a:/
vowels <- c('i:','a:')
subdata <- vt.data[vt.data$word %in% c("bat","biete","Dieter","Rate","Rita","Tat"),]
subdata <- subdata[subdata$stress=="N",]
subdata$vowel <- droplevels(subdata$vowel)
subdata$vowel <- factor(subdata$vowel,levels=vowels)
contrasts(subdata$vowel) <- "contr.treatment"

subdata$AR.start  <- subdata$gridline==1





# NB: optimal knot (k) parameters were determined by model diagnosis using gamm.check()
m1 <- mgcv::bam(aperture ~ vowel
                + te(time.norm, gridline.norm, k=20)
                + te(time.norm, gridline.norm, by=vowel, k=15)
                + s(time.norm, speaker, bs='fs', m=1, xt="tp", k=4)
                + s(gridline.norm, speaker, bs='fs', m=1, xt="tp", k=4)
                + s(word, bs='re', m=1),
                data=subdata, method="fREML")

# optimize rho parameter
valRho <- itsadug::acf_resid(m1, split_pred = c("gridline.norm"),plot = F)[2]

m1.rho <- mgcv::bam(aperture ~ vowel
                    + te(time.norm, gridline.norm, k=20)
                    + te(time.norm, gridline.norm, by=vowel, k=15)
                    + s(time.norm, speaker, bs='fs', m=1, xt="tp", k=4)
                    + s(gridline.norm, speaker, bs='fs', m=1, xt="tp", k=4)
                    + s(word, bs='re', m=1),
                    AR.start=AR.start, rho=valRho,
                    data=subdata, method="fREML")




par(mfrow=c(1,2), cex=0.9)
p1  <- acf_resid(m1, split_pred = c("time.norm"),plot=F)
p2  <- acf_resid(m1.rho, split_pred = c("time.norm"),plot=F)
p   <- plot(p1,col='blue',type='h',lwd=5,main="GAMM autocorrelation: normalized time",xlab="Lag",ylab="ACF")
lines(p2, col='red', type='h', lwd=5)
legend("topright", legend=c("Before correction", "After correction"),
       col=c("blue", "red"), lty=1, cex=1.2, lwd=5)

p1  <- acf_resid(m1, split_pred = c("gridline.norm"),plot=F)
p2  <- acf_resid(m1.rho, split_pred = c("gridline.norm"),plot=F)
p   <- plot(p1,col='blue',type='h',lwd=5,main="GAMM autocorrelation: vocal tract grid lines",xlab="Lag",ylab="ACF")
lines(p2, col='red',type='h',lwd=5)
legend("topright", legend=c("Before correction", "After correction"),
       col=c("blue", "red"), lty=1, cex=1.2, lwd=5)




mylims    <- c(min(m1.rho$fitted.values), max(m1.rho$fitted.values))
mapcols   <- rev(colorRamps::matlab.like2(100))
diffcols <- rev(colorRampPalette(colorschemes$BluetoOrange.12, space = "Lab")(12))
#labels for y-axis plotting
ticknames<-c('glottis','hypo-\n pharynx','hyper-\n pharynx','velum','palate','alveolar\n ridge')
ticks <- c(0,0.2,0.4,0.6,0.8,1)
ticks2 <- c(1,5.6,11.2,16.8,22.4,28)


tiff(filename=paste0("GAMM_heatmaps_",vowels[1],"-",vowels[2],".tiff"),h=6,w=14,units="in",res=280,pointsize=24)

par(mfrow=c(1,2), cex=1, mar=c(3.5,4,1.5,2), mgp=c(1.5,0.75,0))

fvisgam(m1.rho, view=c("time.norm", "gridline.norm"),
        cond=list(vowel=vowels[1]),add.color.legend=F,
        main=paste0("VT aperture of /",vowels[1],"/"),ylab='',yaxt='n',xlab='Time (normalized)',
        rm.ranef=T,zlim = c(-1,17),color=mapcols,hide.label=T,
        n.grid=30,nCol=100,cex.lab=0.8,cex.axis=0.7,cex.main=0.8)
gradientLegend(c(0,16), pos=.875, side=4, color=mapcols, inside = F)
axis(2, at=ticks2, labels=ticknames, las=2, cex.axis=0.7)

par(cex=1, mar=c(3.5,3.5,1.5,2.5), mgp=c(1.5,0.75,0))
fvisgam(m1.rho, view=c("time.norm", "gridline.norm"),
        cond=list(vowel=vowels[2]),add.color.legend=F,
        main=paste0("VT aperture of /",vowels[2],"/"),ylab='',yaxt='n',xlab='Time (normalized)',
        rm.ranef=T,zlim = c(-1,17),color=mapcols,hide.label=T,
        n.grid=30,nCol=100,cex.lab=0.8,cex.axis=0.7,cex.main=0.8)
gradientLegend(c(0,16), pos=.875, side=4, color=mapcols, inside = F)
axis(2, at=ticks2, labels=ticknames, las=2, cex.axis=0.7)
dev.off()



# highlight significant differences
tiff(filename=paste0("GAMM_difference_",vowels[1],"-",vowels[2],".tiff"),h=7,w=9,units="in",res=300,pointsize=24)

par(mfrow=c(1,1), cex=1, mar=c(4,5,2,3), mgp=c(1.5,0.75,0))
plot_diff2(m1.rho, view=c("time.norm", "gridline.norm"), 
           comp=list(vowel=vowels),
           main=paste0("VT aperture difference of /",vowels[1],"/ − /",vowels[2],"/"),ylab='',yaxt='n',xlab='Time (normalized)',
           zlim=c(-14,14),show.diff = T,color=diffcols,add.color.legend = F,alpha.diff = 0.4,rm.ranef = F,hide.label = T,col="black",
           n.grid=30,cex.lab=0.8,cex.axis=0.7,cex.main=0.8)
gradientLegend(c(-14,14), pos=.875, side=4, color=diffcols, inside=F)
abline(v=0.2,lwd=1.5,lty=2)
abline(v=0.8,lwd=1.5,lty=2)
axis(2, at=ticks2, labels=ticknames, las = 1, cex.axis=0.8)
dev.off()


# Diphthongs: differences between /aI/ and /a:/
vowels <- c('aI','a:')

subdata <- vt.data[vt.data$word %in% c("bahne","bat","wate","weine","weinte","weihte"),]
subdata <- subdata[subdata$stress=="N",]
subdata$vowel <- droplevels(subdata$vowel)
subdata$vowel <- factor(subdata$vowel,levels=vowels)
contrasts(subdata$vowel) <- "contr.treatment"

subdata$AR.start  <- subdata$gridline==1


# NB: optimal knot (k) parameters were determined by model diagnosis using gamm.check()
m2 <- mgcv::bam(aperture ~ vowel
                + te(time.norm, gridline.norm, k=15)
                + te(time.norm, gridline.norm, by=vowel, k=15)
                + s(time.norm, speaker, bs='fs', m=1, xt="tp", k=3)
                + s(gridline.norm, speaker, bs='fs', m=1, xt="tp", k=3)
                + s(word, bs='re', m=1),
                data=subdata, method="fREML")

# optimize rho parameter
valRho <- itsadug::acf_resid(m2, split_pred = c("gridline.norm"),plot = F)[2]

m2.rho <- mgcv::bam(aperture ~ vowel
                    + te(time.norm, gridline.norm, k=15)
                    + te(time.norm, gridline.norm, by=vowel, k=15)
                    + s(time.norm, speaker, bs='fs', m=1, xt="tp", k=3)
                    + s(gridline.norm, speaker, bs='fs', m=1, xt="tp", k=3)
                    + s(word, bs='re', m=1),
                    AR.start=AR.start, rho=valRho,
                    data=subdata, method="fREML")




par(mfrow=c(1,2), cex=0.9)
p1  <- acf_resid(m2, split_pred = c("time.norm"),plot=F)
p2  <- acf_resid(m2.rho, split_pred = c("time.norm"),plot=F)
p   <- plot(p1,col='blue',type='h',lwd=5,main="GAMM autocorrelation: normalized time",xlab="Lag",ylab="ACF")
lines(p2, col='red', type='h', lwd=5)
legend("topright", legend=c("Before correction", "After correction"),
       col=c("blue", "red"), lty=1, cex=1.2, lwd=5)

p1  <- acf_resid(m2, split_pred = c("gridline.norm"),plot=F)
p2  <- acf_resid(m2.rho, split_pred = c("gridline.norm"),plot=F)
p   <- plot(p1,col='blue',type='h',lwd=5,main="GAMM autocorrelation: vocal tract grid lines",xlab="Lag",ylab="ACF")
lines(p2, col='red',type='h',lwd=5)
legend("topright", legend=c("Before correction", "After correction"),
       col=c("blue", "red"), lty=1, cex=1.2, lwd=5)




# highlight significant differences
tiff(filename=paste0("GAMM_difference_",vowels[1],"-",vowels[2],".tiff"),h=7,w=9,units="in",res=300,pointsize=24)

par(mfrow=c(1,1), cex=1, mar=c(4,5,2,3), mgp=c(1.5,0.75,0))
plot_diff2(m2.rho, view=c("time.norm", "gridline.norm"), 
           comp=list(vowel=vowels),
           main=paste0("VT aperture difference of /",vowels[1],"/ − /",vowels[2],"/"),ylab='',yaxt='n',xlab='Time (normalized)',
           zlim=c(-8,8),show.diff = T,color=diffcols,add.color.legend = F,alpha.diff = 0.4,rm.ranef = F,hide.label = T,col="black",
           n.grid=30,cex.lab=0.8,cex.axis=0.7,cex.main=0.8)
gradientLegend(c(-8,8), pos=.875, side=4, color=diffcols, inside=F)
abline(v=0.2,lwd=1.5,lty=2)
abline(v=0.8,lwd=1.5,lty=2)
axis(2, at=ticks2, labels=ticknames, las = 1, cex.axis=0.8)
dev.off()



# Stress: Differences between accentuaed and neutral vowels
stresses <- c('A','N')

subdata <- vt.data[vt.data$word %in% c('ahnde','ahnte','sahnst','sahnt','sahst'),]
subdata <- subdata[subdata$stress %in% stresses,]
subdata$stress <- droplevels(subdata$stress)
subdata$stress <- factor(subdata$stress,ordered=F)
contrasts(subdata$stress) <- "contr.treatment"

subdata$AR.start  <- subdata$gridline==1


# NB: optimal knot (k) parameters were determined by model diagnosis using gamm.check()
m3 <- mgcv::bam(aperture ~ stress
                + te(time.norm, gridline.norm, k=12)
                + te(time.norm, gridline.norm, by=stress, k=12)
                + s(time.norm, speaker, bs='fs', m=1, xt="tp", k=3)
                + s(gridline.norm, speaker, bs='fs', m=1, xt="tp", k=3)
                + s(word, bs='re', m=1),
                data=subdata, method="fREML")

# optimize rho parameter
valRho <- itsadug::acf_resid(m3, split_pred = c("gridline.norm"),plot = F)[2]

m3.rho <- mgcv::bam(aperture ~ stress
                    + te(time.norm, gridline.norm, k=12)
                    + te(time.norm, gridline.norm, by=stress, k=12)
                    + s(time.norm, speaker, bs='fs', m=1, xt="tp", k=3)
                    + s(gridline.norm, speaker, bs='fs', m=1, xt="tp", k=3)
                    + s(word, bs='re', m=1),
                    AR.start=AR.start, rho=valRho,
                    data=subdata, method="fREML")


# highlight significant differences
tiff(filename=paste0("GAMM_difference_stress_",stresses[1],"-",stresses[2],".tiff"),h=7,w=9,units="in",res=300,pointsize=24)

par(mfrow=c(1,1), cex=1, mar=c(4,5,2,3), mgp=c(1.5,0.75,0))
plot_diff2(m3.rho, view=c("time.norm", "gridline.norm"), 
           comp=list(stress=stresses),
           main=paste0(stresses[1]," − ",stresses[2]," stress difference in /a:/"),ylab='',yaxt='n',xlab='Time (normalized)',
           zlim=c(-8,8),show.diff = T,color=diffcols,add.color.legend = F,alpha.diff = 0.4,rm.ranef = F,hide.label = T,col="black",
           n.grid=30,cex.lab=0.8,cex.axis=0.7,cex.main=0.8)
gradientLegend(c(-8,8), pos=.875, side=4, color=diffcols, inside=F)
abline(v=0.2,lwd=1.5,lty=2)
abline(v=0.8,lwd=1.5,lty=2)
axis(2, at=ticks2, labels=ticknames, las = 1, cex.axis=0.8)
dev.off()