# Filename: LabPhon_rtMRI_GAMMs.R
# Date: 2019-11-14
# Author: Christopher Carignan
# Email: c.carignan@ucl.ac.uk
# Institution: Speech, Hearing & Phonetic Sciences, University College London
# Description:
#   Code to recreate the generalized additive mixed models (GAMMs) and associated figures appearing in:
#   Carignan, C., Hoole, P., Kunay, E., Pouplier, M., Joseph, A., Voit, D., Frahm, J., & Harrington, J. (under revision), 
#     "Analyzing speech in both time and space: Generalized additive mixed models can uncover systematic patterns of variation 
#     in vocal tract shape in real-time MRI", Laboratory Phonology.


# Load required packages
require(mgcv)
require(itsadug)
require(dichromat)

# Load MRI vocal tract aperture data
load("GAMM_data.Rda")


## plotting parameters to use throughout the script
# colors
mapcols  <- rev(colorRamps::matlab.like2(100))
diffcols <- rev(colorRampPalette(colorschemes$BluetoOrange.12, space="Lab")(12))
# labels for axes
ticknames <- c("glottis","hypo-\n pharynx","hyper-\n pharynx","velum","palate","alveolar\n ridge")
tickvals  <- c(1,5.6,11.2,16.8,22.4,28)



# # # # # # # # # # # # # # # # # # # # # # # # # #
# Monophthongs: Differences between /a:/ and /i:/ #
# # # # # # # # # # # # # # # # # # # # # # # # # #

# Subset the data 
vowels  <- c("a:","i:")
subdata <- vt.data[vt.data$word %in% c("bat","biete","Dieter","Rate","Rita","Tat"),]
subdata <- subdata[subdata$stress=="N",]
subdata$vowel <- droplevels(subdata$vowel)
subdata$vowel <- factor(subdata$vowel, levels=vowels, ordered=T)
contrasts(subdata$vowel) <- "contr.treatment"


## Create the GAM models
# NB: optimal knot (k) parameters were determined by model diagnosis using gamm.check()

# Base model
m1 <- mgcv::bam(aperture ~ vowel
                + te(time.norm, gridline.norm, k=12)
                + te(time.norm, gridline.norm, by=vowel, k=12)
                + s(time.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                + s(gridline.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                + s(word, bs="re", m=1),
                data=subdata, method="fREML")

# Optimize the rho parameter to correct for correlation throughout the vocal tract
valRho <- itsadug::acf_resid(m1, split_pred=c("gridline.norm"), plot=F)[2]

# New model, corrected for autocorrelation
m1.rho <- mgcv::bam(aperture ~ vowel
                    + te(time.norm, gridline.norm, k=12)
                    + te(time.norm, gridline.norm, by=vowel, k=12)
                    + s(time.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                    + s(gridline.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                    + s(word, bs="re", m=1),
                    AR.start=grid.start, rho=valRho,
                    data=subdata, method="fREML")

# Summary of final model
summary(m1.rho)

# Inspect random effects (Figure 5 in paper)
par(mfrow=c(1,2),mar=c(4,5,1,1))
plot(m1.rho, select=3)
plot(m1.rho, select=4)


## Plot separate articulatory heatmaps for /a:/ and /i:/, save as TIFF (Figure 3 in paper)
tiff(filename="GAMM_heatmaps_a:-i:.tiff", h=6, w=14, units="in", res=280, pointsize=24)

# Articulatory heatmap figure for /a:/
# NB: change the n.grid value from 30 (default) to 500 to create the high-resolution figures shown in the paper
par(mfrow=c(1,2), cex=1, mar=c(3,3.5,1.5,2), mgp=c(1.5,0.75,0))
itsadug::fvisgam(m1.rho, view=c("time.norm", "gridline.norm"),
                 cond=list(vowel=vowels[1]), add.color.legend=F,
                 main="VT aperture of /a:/", xlab="Time (normalized)", ylab="", yaxt="n",
                 rm.ranef=F, zlim=c(-1.5,17), color=mapcols, hide.label=T,
                 n.grid=30, nCol=100, cex.lab=0.8, cex.axis=0.7, cex.main=0.8)
gradientLegend(c(0,16), pos=.875, side=4, color=mapcols, inside = F)
axis(2, at=tickvals, labels=ticknames, las=2, cex.axis=0.7)

# Articulatory heatmap figure for /i:/
# NB: change the n.grid value from 30 (default) to 500 to create the high-resolution figures shown in the paper
par(cex=1, mar=c(3,3.5,1.5,2.25), mgp=c(1.5,0.75,0))
itsadug::fvisgam(m1.rho, view=c("time.norm", "gridline.norm"),
                 cond=list(vowel=vowels[2]), add.color.legend=F,
                 main="VT aperture of /i:/", xlab="Time (normalized)", ylab="", yaxt="n",
                 rm.ranef=F, zlim = c(-1.5,17), color=mapcols, hide.label=T,
                 n.grid=30, nCol=100, cex.lab=0.8, cex.axis=0.7, cex.main=0.8)
gradientLegend(c(0,16), pos=.875, side=4, color=mapcols, inside = F)
axis(2, at=tickvals, labels=ticknames, las=2, cex.axis=0.7)
dev.off()


## Plot articulatory differences between /a:/ and /i:/, highlighting areas of significant difference, save as TIFF (Figure 4 in paper)
tiff(filename="GAMM_difference_i:-a:.tiff", h=7, w=9, units="in", res=300, pointsize=24)

# NB: change the n.grid value from 30 (default) to 500 to create the high-resolution figures shown in the paper
par(mfrow=c(1,1), cex=1, mar=c(3,3.5,1.5,2.5), mgp=c(1.5,0.75,0))
itsadug::plot_diff2(m1.rho, view=c("time.norm", "gridline.norm"), comp=list(vowel=vowels),
                    main="VT aperture difference of /a:/ − /i:/", xlab="Time (normalized)", ylab="", yaxt="n",
                    rm.ranef=F, show.diff=T, zlim=c(-14.2,14.2), alpha.diff=0.4, hide.label=T,
                    color=diffcols, add.color.legend=F,  col="black",
                    n.grid=30, cex.lab=0.8, cex.axis=0.7, cex.main=0.8)
gradientLegend(c(-14,14), pos=.875, side=4, color=diffcols, inside=F)
abline(v=0.2, lwd=1.5, lty=2)
abline(v=0.8, lwd=1.5, lty=2)
axis(2, at=tickvals, labels=ticknames, las=1, cex.axis=0.8)
dev.off()



# # # # # # # # # # # # # # # # # # # # # # # # #
# Diphthongs: Differences between /a:/ and /aI/ #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Subset the data
vowels  <- c("a:","aI")
subdata <- vt.data[vt.data$word %in% c("bahne","bat","wate","weine","weinte","weihte"),]
subdata <- subdata[subdata$stress=="N",]
subdata$vowel <- droplevels(subdata$vowel)
subdata$vowel <- factor(subdata$vowel, levels=vowels, ordered=T)
contrasts(subdata$vowel) <- "contr.treatment"


## Create the GAM models
# NB: optimal knot (k) parameters were determined by model diagnosis using gamm.check()

# Base model
m2 <- mgcv::bam(aperture ~ vowel
                + te(time.norm, gridline.norm, k=15)
                + te(time.norm, gridline.norm, by=vowel, k=15)
                + s(time.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                + s(gridline.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                + s(word, bs="re", m=1),
                data=subdata, method="fREML")

# Optimize the rho parameter to correct for correlation throughout the vocal tract
valRho <- itsadug::acf_resid(m2, split_pred=c("gridline.norm"), plot=F)[2]

# New model, corrected for autocorrelation
m2.rho <- mgcv::bam(aperture ~ vowel
                    + te(time.norm, gridline.norm, k=15)
                    + te(time.norm, gridline.norm, by=vowel, k=15)
                    + s(time.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                    + s(gridline.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                    + s(word, bs="re", m=1),
                    AR.start=grid.start, rho=valRho,
                    data=subdata, method="fREML")

# Summary of final model
summary(m2.rho)

# Inspect random effects (Figure 8 in paper)
par(mfrow=c(1,2),mar=c(4,5,1,1))
plot(m2.rho, select=3)
plot(m2.rho, select=4)


## Plot articulatory differences between /a:/ and /aI/, highlighting areas of significant difference, save as TIFF (Figure 7 in paper)
tiff(filename="GAMM_difference_a:-aI.tiff", h=7, w=9, units="in", res=300, pointsize=24)

# NB: change the n.grid value from 30 (default) to 500 to create the high-resolution figures shown in the paper
par(mfrow=c(1,1), cex=1, mar=c(3,3.5,1.5,2.5), mgp=c(1.5,0.75,0))
itsadug::plot_diff2(m2.rho, view=c("time.norm", "gridline.norm"), comp=list(vowel=vowels),
                    main="VT aperture difference of /a:/ - /aɪ/", xlab="Time (normalized)", ylab="", yaxt="n",
                    rm.ranef=F, show.diff=T, zlim=c(-8,8), alpha.diff=0.4, hide.label=T, 
                    color=diffcols, add.color.legend=F,  col="black",
                    n.grid=30, cex.lab=0.8, cex.axis=0.7, cex.main=0.8)
gradientLegend(c(-8,8), pos=.875, side=4, color=diffcols, inside=F)
abline(v=0.2, lwd=1.5, lty=2)
abline(v=0.8, lwd=1.5, lty=2)
axis(2, at=tickvals, labels=ticknames, las=1, cex.axis=0.8)
dev.off()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Stress: Differences between accentuaed and neutral vowels #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Subset the data
stresses <- c("A","N")
subdata  <- vt.data[vt.data$word %in% c('ahnde','ahnte','sahnst','sahnt','sahst'),]
subdata  <- subdata[subdata$stress %in% stresses,]
subdata$stress <- droplevels(subdata$stress)
subdata$stress <- factor(subdata$stress, levels=stresses, ordered=T)
contrasts(subdata$stress) <- "contr.treatment"


## Create the GAM models
# NB: optimal knot (k) parameters were determined by model diagnosis using gamm.check()

# Base model
m3 <- mgcv::bam(aperture ~ stress
                + te(time.norm, gridline.norm, k=15)
                + te(time.norm, gridline.norm, by=stress, k=15)
                + s(time.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                + s(gridline.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                + s(word, bs="re", m=1),
                data=subdata, method="fREML")

# Optimize the rho parameter to correct for correlation throughout the vocal tract
valRho <- itsadug::acf_resid(m3, split_pred=c("gridline.norm"), plot=F)[2]

# New model, corrected for autocorrelation
m3.rho <- mgcv::bam(aperture ~ stress
                    + te(time.norm, gridline.norm, k=15)
                    + te(time.norm, gridline.norm, by=stress, k=15)
                    + s(time.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                    + s(gridline.norm, speaker, bs="fs", m=1, xt="tp", k=4)
                    + s(word, bs="re", m=1),
                    AR.start=grid.start, rho=valRho,
                    data=subdata, method="fREML")

# Summary of final model
summary(m3.rho)

# Inspect random effects (Figure 11 in paper)
par(mfrow=c(1,2),mar=c(4,5,1,1))
plot(m3.rho, select=3)
plot(m3.rho, select=4)


## Plot articulatory differences between /i:/ and /a:/, highlighting areas of significant difference, save as TIFF (Figure 10 in paper)
tiff(filename="GAMM_difference_stress_A-N.tiff", h=7, w=9, units="in", res=300, pointsize=24)

# NB: change the n.grid value from 30 (default) to 500 to create the high-resolution figures shown in the paper
par(mfrow=c(1,1), cex=1, mar=c(3,3.5,1.5,2.5), mgp=c(1.5,0.75,0))
itsadug::plot_diff2(m3.rho, view=c("time.norm", "gridline.norm"), comp=list(stress=stresses),
                    main="A − N stress difference in /a:/", xlab="Time (normalized)", ylab="", yaxt="n",
                    rm.ranef=F, show.diff=T, zlim=c(-8,8), alpha.diff=0.4, hide.label=T, 
                    color=diffcols, add.color.legend=F,  col="black",
                    n.grid=30, cex.lab=0.8, cex.axis=0.7, cex.main=0.8)
gradientLegend(c(-8,8), pos=.875, side=4, color=diffcols, inside=F)
abline(v=0.2, lwd=1.5, lty=2)
abline(v=0.8, lwd=1.5, lty=2)
axis(2, at=tickvals, labels=ticknames, las=1, cex.axis=0.8)
dev.off()