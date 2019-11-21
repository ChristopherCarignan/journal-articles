# Filename: LabPhon_rtMRI_FLMMs.R
# Date: 2019-11-14
# Author: Christopher Carignan
# Email: c.carignan@ucl.ac.uk
# Institution: Speech, Hearing & Phonetic Sciences, University College London
# Description:
#   Code to recreate the functional linear mixed models (FLMMs) and associated figures appearing in:
#   Carignan, C., Hoole, P., Kunay, E., Pouplier, M., Joseph, A., Voit, D., Frahm, J., & Harrington, J. (under revision), 
#     "Analyzing speech in both time and space: Generalized additive mixed models can uncover systematic patterns of variation 
#     in vocal tract shape in real-time MRI", Laboratory Phonology.


# Load required packages
library(sparseFLMM)
library(data.table)

# Load MRI vocal tract aperture data for 20% and 80% of the vowel interval
load("FLMM_data.Rda")


## plotting parameters to use throughout the script
cols    <- rainbow(12)
xticks  <- seq(0,1, by=0.2)
xlabs   <- c("glottis","hypo-pharynx","hyper-pharynx","velum","palate","alveolar ridge")
yticks  <- seq(0,20, by=2)
miny    <- 0
maxy    <- 16



# # # # # # # # # # # # # # # # # # # # # # # # # #
# Monophthongs: Differences between /i:/ and /a:/ #
# # # # # # # # # # # # # # # # # # # # # # # # # #

# Subset the data 
vowels  <- c("a:","i:")
subdata <- vt.data[vt.data$word %in% c("bat","biete","Dieter","Rate","Rita","Tat"),]
subdata <- subdata[subdata$stress=="N",]
subdata$vowel <- droplevels(subdata$vowel)
subdata$vowel <- factor(subdata$vowel, levels=vowels, ordered=T)
contrasts(subdata$vowel) <- "contr.treatment"


# Prepare FLMM table for 20% of the vowel interval
FLMM.table <- subdata[,c("aperture.20","gridline.scale","ident","speaker","word","vowel")]
names(FLMM.table) <- c("y_vec","t","n_long","subject_long","word_long","covariate.1")

# Identity coding for covariates
x <- 1
for (i in unique(FLMM.table$n_long)){
  FLMM.table$curve1[FLMM.table$n_long==i] <- x
  x <- x+1
}

x <- 1
for (i in unique(FLMM.table$subject_long)){
  FLMM.table$subject1[FLMM.table$subject_long==i] <- x
  x <- x+1
}

x <- 1
for (i in unique(FLMM.table$word_long)){
  FLMM.table$word1[FLMM.table$word_long==i] <- x
  x <- x+1
}

# Dummy code factor levels
FLMM.table$covariate.new[FLMM.table$covariate.1==vowels[1]] <- 0
FLMM.table$covariate.new[FLMM.table$covariate.1==vowels[2]] <- 1

# Create the FLMM table
FLMM.table <- cbind(FLMM.table[,c("y_vec","t","curve1","subject1","word1","covariate.new")])
names(FLMM.table) <- c("y_vec","t","n_long","subject_long","word_long","covariate.1")
FLMM.table$n_long <- as.integer(FLMM.table$n_long)
FLMM.table$subject_long <- as.integer(FLMM.table$subject_long)
FLMM.table <- data.table::as.data.table(FLMM.table)


# Build the FLMM for 20% of the vowel interval
flmm.20 <- sparseFLMM::sparseFLMM(curve_info=FLMM.table,
                                  use_RI=T, 
                                  use_simple=F,
                                  covariate=T,
                                  num_covariates=1,
                                  covariate_form="by",
                                  interaction=T,
                                  bf_covs=c(10,10),
                                  m_covs=list(c(2,3), c(2,3)),
                                  use_famm=T
)


# Build the FLMM for 80% of the vowel interval
FLMM.table$y_vec <- subdata$aperture.80

flmm.80 <- sparseFLMM::sparseFLMM(curve_info=FLMM.table,
                                  use_RI=T, 
                                  use_simple=F,
                                  covariate=T,
                                  num_covariates=1,
                                  covariate_form="by",
                                  interaction=T,
                                  bf_covs=c(10,10),
                                  m_covs=list(c(2,3), c(2,3)),
                                  use_famm=T
)


# Extract curve info for plotting (20% of vowel interval)
intercept.20  <- flmm.20$fpc_famm_hat_tri_constr$intercept
y_mean.20     <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_mean$value
se_mean.20    <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_mean$se
my_grid.20    <- flmm.20$my_grid

# Effect of covariate (20% of vowel interval)
y_cov1.20     <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_covariate.1$value
se_cov1.20    <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_covariate.1$se

# Sum the intercept, mean and intercept, mean, covariate effect (20% of vowel interval)
ref.fit.20    <- intercept.20 + y_mean.20
cov.fit.20    <- intercept.20 + y_mean.20 + y_cov1.20


# Extract curve info for plotting (80% of vowel interval)
intercept.80  <- flmm.80$fpc_famm_hat_tri_constr$intercept
y_mean.80     <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_mean$value
se_mean.80    <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_mean$se
my_grid.80    <- flmm.80$my_grid

# Effect of covariate (80% of vowel interval)
y_cov1.80     <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_covariate.1$value
se_cov1.80    <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_covariate.1$se

# Sum the intercept, mean and intercept, mean, covariate effect (80% of vowel interval)
ref.fit.80    <- intercept.80 + y_mean.80
cov.fit.80    <- intercept.80 + y_mean.80 + y_cov1.80


## PLOTTING ##

# Figure 6 in paper
cairo_pdf("FLMM_a:-i:.pdf", h=7, w=26, onefile=T)

# Prepare the 2-panel plot layout
par(mfrow=c(1,2))


# Plot the reference mean and then the covariate effect (20% of vowel interval)
par(mar=c(5,6,4,2)+0.1, mgp=c(4,1,0))
plot(my_grid.20, ref.fit.20, col=cols[1], lty=1, lwd=4, t="l",
     ylim=c(miny, maxy),
     ylab="", xlab="", xaxt="n", yaxt="n")
abline(h=yticks, lwd=1.5, lty=2)
par(new=T)
plot(my_grid.20, cov.fit.20, col=cols[10], lty=4, lwd=4,
     ylim=c(miny, maxy), xaxt="n", yaxt="n",
     ylab="Average VT aperture (mm)", xlab="",
     main="Vocal tract aperture at 20% of vowel interval", t="l",
     cex.main=2.5, font.main=1, cex.lab=1.8)
axis(1, at=xticks, labels=xlabs, las=1, cex.axis=1.7)
axis(2, at=yticks, labels=yticks, las=1, cex.axis=1.7)

# Add confidence interval bands to the plot
lines(x=my_grid.20, y=(y_mean.20 + se_mean.20 + intercept.20), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.20, y=(y_mean.20 - se_mean.20 + intercept.20), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.20, y=(cov.fit.20 + se_cov1.20), lty=1, lwd=2,
      col=cols[10])
lines(x=my_grid.20, y=(cov.fit.20 - se_cov1.20), lty=1, lwd=2,
      col=cols[10])

# Create legend
legend("topleft", legend=c(vowels[1], vowels[2]), bty="o", col=c(cols[1], cols[10]),
       lty=c(1,4), lwd=4, cex=1.5, box.col="black")


# Plot the reference mean and then the covariate effect (80% of vowel interval)
par(mar=c(5,6,4,2)+0.1, mgp=c(4,1,0))
plot(my_grid.80, ref.fit.80, col=cols[1], lty=1, lwd=4, t="l",
     ylim=c(miny, maxy),
     ylab="", xlab="", xaxt="n", yaxt="n")
abline(h=yticks, lwd=1.5, lty=2)
par(new=T)
plot(my_grid.80, cov.fit.80, col=cols[10], lty=4, lwd=4,
     ylim=c(miny, maxy), xaxt="n", yaxt="n",
     ylab="Average VT aperture (mm)", xlab="",
     main="Vocal tract aperture at 80% of vowel interval", t="l",
     cex.main=2.5, font.main=1, cex.lab=1.8)
axis(1, at=xticks, labels=xlabs, las=1, cex.axis=1.7)
axis(2, at=yticks, labels=yticks, las=1, cex.axis=1.7)

# Add confidence interval bands to the plot
lines(x=my_grid.80, y=(y_mean.80 + se_mean.80 + intercept.80), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.80, y=(y_mean.80 - se_mean.80 + intercept.80), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.80, y=(cov.fit.80 + se_cov1.80), lty=1, lwd=2,
      col=cols[10])
lines(x=my_grid.80, y=(cov.fit.80 - se_cov1.80), lty=1, lwd=2,
      col=cols[10])

# Create legend
legend("topleft", legend=c(vowels[1], vowels[2]), bty="o", col=c(cols[1], cols[10]),
       lty=c(1,4), lwd=4, cex=1.5, box.col="black")

dev.off()



# # # # # # # # # # # # # # # # # # # # # # # # #
# Diphthongs: Differences between /aI/ and /a:/ #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Subset the data
vowels  <- c("a:","aI")
subdata <- vt.data[vt.data$word %in% c("bahne","bat","wate","weine","weinte","weihte"),]
subdata <- subdata[subdata$stress=="N",]
subdata$vowel <- droplevels(subdata$vowel)
subdata$vowel <- factor(subdata$vowel, levels=vowels, ordered=T)
contrasts(subdata$vowel) <- "contr.treatment"


# Prepare FLMM table for 20% of the vowel interval
FLMM.table <- subdata[,c("aperture.20","gridline.scale","ident","speaker","word","vowel")]
names(FLMM.table) <- c("y_vec","t","n_long","subject_long","word_long","covariate.1")

# Identity coding for covariates
x <- 1
for (i in unique(FLMM.table$n_long)){
  FLMM.table$curve1[FLMM.table$n_long==i] <- x
  x <- x+1
}

x <- 1
for (i in unique(FLMM.table$subject_long)){
  FLMM.table$subject1[FLMM.table$subject_long==i] <- x
  x <- x+1
}

x <- 1
for (i in unique(FLMM.table$word_long)){
  FLMM.table$word1[FLMM.table$word_long==i] <- x
  x <- x+1
}

# Dummy code factor levels
FLMM.table$covariate.new[FLMM.table$covariate.1==vowels[1]] <- 0
FLMM.table$covariate.new[FLMM.table$covariate.1==vowels[2]] <- 1

# Create the FLMM table
FLMM.table <- cbind(FLMM.table[,c("y_vec","t","curve1","subject1","word1","covariate.new")])
names(FLMM.table) <- c("y_vec","t","n_long","subject_long","word_long","covariate.1")
FLMM.table$n_long <- as.integer(FLMM.table$n_long)
FLMM.table$subject_long <- as.integer(FLMM.table$subject_long)
FLMM.table <- data.table::as.data.table(FLMM.table)


# Build the FLMM for 20% of the vowel interval
flmm.20 <- sparseFLMM::sparseFLMM(curve_info=FLMM.table,
                                  use_RI=T, 
                                  use_simple=F,
                                  covariate=T,
                                  num_covariates=1,
                                  covariate_form="by",
                                  interaction=T,
                                  bf_covs=c(10,10),
                                  m_covs=list(c(2,3), c(2,3)),
                                  use_famm=T
)


# Build the FLMM for 80% of the vowel interval
FLMM.table$y_vec <- subdata$aperture.80

flmm.80 <- sparseFLMM::sparseFLMM(curve_info=FLMM.table,
                                  use_RI=T, 
                                  use_simple=F,
                                  covariate=T,
                                  num_covariates=1,
                                  covariate_form="by",
                                  interaction=T,
                                  bf_covs=c(10,10),
                                  m_covs=list(c(2,3), c(2,3)),
                                  use_famm=T
)


# Extract curve info for plotting (20% of vowel interval)
intercept.20  <- flmm.20$fpc_famm_hat_tri_constr$intercept
y_mean.20     <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_mean$value
se_mean.20    <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_mean$se
my_grid.20    <- flmm.20$my_grid

# Effect of covariate (20% of vowel interval)
y_cov1.20     <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_covariate.1$value
se_cov1.20    <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_covariate.1$se

# Sum the intercept, mean and intercept, mean, covariate effect (20% of vowel interval)
ref.fit.20    <- intercept.20 + y_mean.20
cov.fit.20    <- intercept.20 + y_mean.20 + y_cov1.20


# Extract curve info for plotting (80% of vowel interval)
intercept.80  <- flmm.80$fpc_famm_hat_tri_constr$intercept
y_mean.80     <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_mean$value
se_mean.80    <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_mean$se
my_grid.80    <- flmm.80$my_grid

# Effect of covariate (80% of vowel interval)
y_cov1.80     <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_covariate.1$value
se_cov1.80    <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_covariate.1$se

# Sum the intercept, mean and intercept, mean, covariate effect (80% of vowel interval)
ref.fit.80    <- intercept.80 + y_mean.80
cov.fit.80    <- intercept.80 + y_mean.80 + y_cov1.80


## PLOTTING ##

# Figure 9 in paper
cairo_pdf("FLMM_a:-aI.pdf", h=7, w=26, onefile=T)

# Prepare the 2-panel plot layout
par(mfrow=c(1,2))


# Plot the reference mean and then the covariate effect (20% of vowel interval)
par(mar=c(5,6,4,2)+0.1, mgp=c(4,1,0))
plot(my_grid.20, ref.fit.20, col=cols[1], lty=1, lwd=4, t="l",
     ylim=c(miny, maxy),
     ylab="", xlab="", xaxt="n", yaxt="n")
abline(h=yticks, lwd=1.5, lty=2)
par(new=T)
plot(my_grid.20, cov.fit.20, col=cols[10], lty=4, lwd=4,
     ylim=c(miny, maxy), xaxt="n", yaxt="n",
     ylab="Average VT aperture (mm)", xlab="",
     main="Vocal tract aperture at 20% of vowel interval", t="l",
     cex.main=2.5, font.main=1, cex.lab=1.8)
axis(1, at=xticks, labels=xlabs, las=1, cex.axis=1.7)
axis(2, at=yticks, labels=yticks, las=1, cex.axis=1.7)

# Add confidence interval bands to the plot
lines(x=my_grid.20, y=(y_mean.20 + se_mean.20 + intercept.20), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.20, y=(y_mean.20 - se_mean.20 + intercept.20), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.20, y=(cov.fit.20 + se_cov1.20), lty=1, lwd=2,
      col=cols[10])
lines(x=my_grid.20, y=(cov.fit.20 - se_cov1.20), lty=1, lwd=2,
      col=cols[10])

# Create legend
legend("topleft", legend=c(vowels[1], vowels[2]), bty="o", col=c(cols[1], cols[10]),
       lty=c(1,4), lwd=4, cex=1.5, box.col="black")


# Plot the reference mean and then the covariate effect (80% of vowel interval)
par(mar=c(5,6,4,2)+0.1, mgp=c(4,1,0))
plot(my_grid.80, ref.fit.80, col=cols[1], lty=1, lwd=4, t="l",
     ylim=c(miny, maxy),
     ylab="", xlab="", xaxt="n", yaxt="n")
abline(h=yticks, lwd=1.5, lty=2)
par(new=T)
plot(my_grid.80, cov.fit.80, col=cols[10], lty=4, lwd=4,
     ylim=c(miny, maxy), xaxt="n", yaxt="n",
     ylab="Average VT aperture (mm)", xlab="",
     main="Vocal tract aperture at 80% of vowel interval", t="l",
     cex.main=2.5, font.main=1, cex.lab=1.8)
axis(1, at=xticks, labels=xlabs, las=1, cex.axis=1.7)
axis(2, at=yticks, labels=yticks, las=1, cex.axis=1.7)

# Add confidence interval bands to the plot
lines(x=my_grid.80, y=(y_mean.80 + se_mean.80 + intercept.80), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.80, y=(y_mean.80 - se_mean.80 + intercept.80), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.80, y=(cov.fit.80 + se_cov1.80), lty=1, lwd=2,
      col=cols[10])
lines(x=my_grid.80, y=(cov.fit.80 - se_cov1.80), lty=1, lwd=2,
      col=cols[10])

# Create legend
legend("topleft", legend=c(vowels[1], vowels[2]), bty="o", col=c(cols[1], cols[10]),
       lty=c(1,4), lwd=4, cex=1.5, box.col="black")

dev.off()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Stress: Differences between accentuaed and neutral vowels #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Subset the data
stresses <- c("A","N")
subdata  <- vt.data[vt.data$word %in% c("ahnde","ahnte","sahnst","sahnt","sahst"),]
subdata  <- subdata[subdata$stress %in% stresses,]
subdata$stress <- droplevels(subdata$stress)
subdata$stress <- factor(subdata$stress, levels=stresses, ordered=T)
contrasts(subdata$stress) <- "contr.treatment"


# Prepare FLMM table for 20% of the vowel interval
FLMM.table <- subdata[,c("aperture.20","gridline.scale","ident","speaker","word","stress")]
names(FLMM.table) <- c("y_vec","t","n_long","subject_long","word_long","covariate.1")

# Identity coding for covariates
x <- 1
for (i in unique(FLMM.table$n_long)){
  FLMM.table$curve1[FLMM.table$n_long==i] <- x
  x <- x+1
}

x <- 1
for (i in unique(FLMM.table$subject_long)){
  FLMM.table$subject1[FLMM.table$subject_long==i] <- x
  x <- x+1
}

x <- 1
for (i in unique(FLMM.table$word_long)){
  FLMM.table$word1[FLMM.table$word_long==i] <- x
  x <- x+1
}

# Dummy code factor levels
FLMM.table$covariate.new[FLMM.table$covariate.1==stresses[1]] <- 0
FLMM.table$covariate.new[FLMM.table$covariate.1==stresses[2]] <- 1

# Create the FLMM table
FLMM.table <- cbind(FLMM.table[,c("y_vec","t","curve1","subject1","word1","covariate.new")])
names(FLMM.table) <- c("y_vec","t","n_long","subject_long","word_long","covariate.1")
FLMM.table$n_long <- as.integer(FLMM.table$n_long)
FLMM.table$subject_long <- as.integer(FLMM.table$subject_long)
FLMM.table <- data.table::as.data.table(FLMM.table)


# Build the FLMM for 20% of the vowel interval
flmm.20 <- sparseFLMM::sparseFLMM(curve_info=FLMM.table,
                                  use_RI=T, 
                                  use_simple=F,
                                  covariate=T,
                                  num_covariates=1,
                                  covariate_form="by",
                                  interaction=T,
                                  bf_covs=c(10,10),
                                  m_covs=list(c(2,3), c(2,3)),
                                  use_famm=T
)


# Build the FLMM for 80% of the vowel interval
FLMM.table$y_vec <- subdata$aperture.80

flmm.80 <- sparseFLMM::sparseFLMM(curve_info=FLMM.table,
                                  use_RI=T, 
                                  use_simple=F,
                                  covariate=T,
                                  num_covariates=1,
                                  covariate_form="by",
                                  interaction=T,
                                  bf_covs=c(10,10),
                                  m_covs=list(c(2,3), c(2,3)),
                                  use_famm=T
)


# Extract curve info for plotting (20% of vowel interval)
intercept.20  <- flmm.20$fpc_famm_hat_tri_constr$intercept
y_mean.20     <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_mean$value
se_mean.20    <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_mean$se
my_grid.20    <- flmm.20$my_grid

# Effect of covariate (20% of vowel interval)
y_cov1.20     <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_covariate.1$value
se_cov1.20    <- flmm.20$fpc_famm_hat_tri_constr$famm_cb_covariate.1$se

# Sum the intercept, mean and intercept, mean, covariate effect (20% of vowel interval)
ref.fit.20    <- intercept.20 + y_mean.20
cov.fit.20    <- intercept.20 + y_mean.20 + y_cov1.20


# Extract curve info for plotting (80% of vowel interval)
intercept.80  <- flmm.80$fpc_famm_hat_tri_constr$intercept
y_mean.80     <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_mean$value
se_mean.80    <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_mean$se
my_grid.80    <- flmm.80$my_grid

# Effect of covariate (80% of vowel interval)
y_cov1.80     <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_covariate.1$value
se_cov1.80    <- flmm.80$fpc_famm_hat_tri_constr$famm_cb_covariate.1$se

# Sum the intercept, mean and intercept, mean, covariate effect (80% of vowel interval)
ref.fit.80    <- intercept.80 + y_mean.80
cov.fit.80    <- intercept.80 + y_mean.80 + y_cov1.80


## PLOTTING ##

# Figure 12 in paper
cairo_pdf("FLMM_stress_A-N.pdf", h=7, w=26, onefile=T)

# Prepare the 2-panel plot layout
par(mfrow=c(1,2))


# Plot the reference mean and then the covariate effect (20% of vowel interval)
par(mar=c(5,6,4,2)+0.1, mgp=c(4,1,0))
plot(my_grid.20, ref.fit.20, col=cols[1], lty=1, lwd=4, t="l",
     ylim=c(miny, maxy),
     ylab="", xlab="", xaxt="n", yaxt="n")
abline(h=yticks, lwd=1.5, lty=2)
par(new=T)
plot(my_grid.20, cov.fit.20, col=cols[10], lty=4, lwd=4,
     ylim=c(miny, maxy), xaxt="n", yaxt="n",
     ylab="Average VT aperture (mm)", xlab="",
     main="Vocal tract aperture at 20% of vowel interval", t="l",
     cex.main=2.5, font.main=1, cex.lab=1.8)
axis(1, at=xticks, labels=xlabs, las=1, cex.axis=1.7)
axis(2, at=yticks, labels=yticks, las=1, cex.axis=1.7)

# Add confidence interval bands to the plot
lines(x=my_grid.20, y=(y_mean.20 + se_mean.20 + intercept.20), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.20, y=(y_mean.20 - se_mean.20 + intercept.20), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.20, y=(cov.fit.20 + se_cov1.20), lty=1, lwd=2,
      col=cols[10])
lines(x=my_grid.20, y=(cov.fit.20 - se_cov1.20), lty=1, lwd=2,
      col=cols[10])

# Create legend
legend("topleft", legend=c(stresses[1], stresses[2]), bty="o", col=c(cols[1], cols[10]),
       lty=c(1,4), lwd=4, cex=1.5, box.col="black")


# Plot the reference mean and then the covariate effect (80% of vowel interval)
par(mar=c(5,6,4,2)+0.1, mgp=c(4,1,0))
plot(my_grid.80, ref.fit.80, col=cols[1], lty=1, lwd=4, t="l",
     ylim=c(miny, maxy),
     ylab="", xlab="", xaxt="n", yaxt="n")
abline(h=yticks, lwd=1.5, lty=2)
par(new=T)
plot(my_grid.80, cov.fit.80, col=cols[10], lty=4, lwd=4,
     ylim=c(miny, maxy), xaxt="n", yaxt="n",
     ylab="Average VT aperture (mm)", xlab="",
     main="Vocal tract aperture at 80% of vowel interval", t="l",
     cex.main=2.5, font.main=1, cex.lab=1.8)
axis(1, at=xticks, labels=xlabs, las=1, cex.axis=1.7)
axis(2, at=yticks, labels=yticks, las=1, cex.axis=1.7)

# Add confidence interval bands to the plot
lines(x=my_grid.80, y=(y_mean.80 + se_mean.80 + intercept.80), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.80, y=(y_mean.80 - se_mean.80 + intercept.80), lty=1, lwd=2,
      col=cols[1])
lines(x=my_grid.80, y=(cov.fit.80 + se_cov1.80), lty=1, lwd=2,
      col=cols[10])
lines(x=my_grid.80, y=(cov.fit.80 - se_cov1.80), lty=1, lwd=2,
      col=cols[10])

# Create legend
legend("topleft", legend=c(stresses[1], stresses[2]), bty="o", col=c(cols[1], cols[10]),
       lty=c(1,4), lwd=4, cex=1.5, box.col="black")

dev.off()