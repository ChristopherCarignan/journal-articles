library(sparseFLMM)

load('FLMM_data.Rda')

# Monophthongs: Differences between /i:/ and /a:/
vowels <- c('a:','i:')

subdata <- subdata[subdata$word %in% c("bat","biete","Dieter","Rate","Rita","Tat"),]
subdata$vowel <- droplevels(subdata$vowel)
subdata$vowel <- factor(subdata$vowel,levels=vowels)
contrasts(subdata$vowel) <- "contr.treatment"
subdata <- subdata[complete.cases(subdata),]


FLMM.table <- subdata[,c('aperture','vt.norm','ident','speaker','word','vowel')]
names(FLMM.table) <- c('y_vec','t','n_long','subject_long','word_long','covariate.1')


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


cov0 <- vowels[1]
cov1 <- vowels[2]

FLMM.table$covariate.new[FLMM.table$covariate.1==cov0] <- 0
FLMM.table$covariate.new[FLMM.table$covariate.1==cov1] <- 1

FLMM.table <- FLMM.table[complete.cases(FLMM.table),]

FLMM.table <- cbind(FLMM.table[,c('y_vec','t','curve1','subject1','word1','covariate.new')])
names(FLMM.table) <- c('y_vec','t','n_long','subject_long','word_long','covariate.1')
FLMM.table$n_long <- as.integer(FLMM.table$n_long)
FLMM.table$subject_long <- as.integer(FLMM.table$subject_long)
FLMM.table <- as.data.table(FLMM.table)

results <- sparseFLMM(curve_info=FLMM.table,
                      use_RI=T,use_simple=F,
                      covariate=T,
                      num_covariates=1,
                      covariate_form='by',
                      interaction=F,
                      bf_covs=c(10,10),
                      m_covs=list(c(2,3), c(2,3)),
                      use_famm=T
)

intercept <- results$fpc_famm_hat_tri_constr$intercept
y_mean <- results$fpc_famm_hat_tri_constr$famm_cb_mean$value
se_mean <- results$fpc_famm_hat_tri_constr$famm_cb_mean$se
my_grid <- results$my_grid

#effect of covariate
y_cov1<-results$fpc_famm_hat_tri_constr$famm_cb_covariate.1$value
se_cov1<-results$fpc_famm_hat_tri_constr$famm_cb_covariate.1$se



cols <- rainbow(12)
ticks <- c(0,0.2,0.4,0.6,0.8,1)

#labels for x-axis plotting
xtick<-c('glottis','hypo-pharynx','hyper-pharynx','velum','palate','alveolar ridge')
#xtick <- c('glottis','hypo-\n pharynx','hyper-\n pharynx','velum','palate','alveolar\n ridge')



cairo_pdf(paste0("FLMM_VT_",timepoint,"_",vowels[1],"-",vowels[2],".pdf"),h=6,w=18,onefile=T)


# summed effect plots
# define colors and scaling of the y-axis
cairo_pdf(paste0("FLMM_VT_",timepoint,"_",vowels[1],"-",vowels[2],"_summary.pdf"),h=7,w=14,onefile=T)

miny <-0
maxy <-16


# sum the intercept, mean and intercept, mean, covariate effect
none <- intercept+y_mean
Status <- intercept+y_mean+y_cov1
#plot first the reference mean and then the covariate effect
par(mfrow=c(1,1))
par(mar=c(5,6,4,2)+0.1,mgp=c(4,1,0))
p <- plot(my_grid, none, col = cols[1], lty = 1, lwd = 4, t = "l",
          ylim = c(miny, maxy),
          ylab = "", xlab = "", xaxt = "n", yaxt="n")
abline(h=0,lwd=1.5,lty=2)
#abline(h=5/resolution,lwd=1.5,lty=2)
abline(h=5,lwd=1.5,lty=2)
#abline(h=10/resolution,lwd=1.5,lty=2)
abline(h=10,lwd=1.5,lty=2)
#abline(h=15/resolution,lwd=1.5,lty=2)
abline(h=15,lwd=1.5,lty=2)
#abline(h=20/resolution,lwd=1.5,lty=2)
abline(h=20,lwd=1.5,lty=2)
abline(h=25,lwd=1.5,lty=2)
abline(h=30,lwd=1.5,lty=2)
par(new = TRUE)
plot(my_grid, Status, col = cols[10], lty = 4, lwd = 4,
     ylim = c(miny, maxy),xaxt = "n", yaxt = "n",
     #ylab = expression("Average VT aperture " (mm^{2})), xlab = "",
     ylab = "Average VT aperture (mm)", xlab = "",
     main = paste0("Summed Effect Curves"), t = "l",
     cex.main = 1.8, font.main = 1, cex.lab = 1.8)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1),labels=xtick, las = 1, cex.axis = 1.7)
#axis(2, at=ticks,labels=round(ticks*resolution,1), las = 1, cex.axis = 1.5)
#axis(2, at=c(0,5,10,15,20)/resolution,labels=c(0,5,10,15,20), las = 1, cex.axis = 1.7)
axis(2, at=c(0,5,10,15,20,25,30),labels=c(0,5,10,15,20,25,30), las = 1, cex.axis = 1.7)

#add confidence bands to the plot
lines(x = my_grid, y = (y_mean+se_mean+intercept), lty = 1, lwd = 2,
      col = cols[1])
lines(x = my_grid, y = (y_mean-se_mean+intercept), lty = 1, lwd = 2,
      col = cols[1])
lines(x = my_grid, y = (Status+se_cov1), lty = 1, lwd = 2,
      col = cols[10])
lines(x = my_grid, y = (Status-se_cov1), lty = 1, lwd = 2,
      col = cols[10])
# legend text
#leg.text <- c(expression(paste(f[0](t))),
#              expression(paste(f[0](t)+f[1](t))))
leg.text <- c(cov0,cov1)
legend("topleft", legend = leg.text, bty = "o", col = c(cols[1], cols[10]),
       lty = c(1,4), lwd = 4,cex = 1.5, box.col='black')

print(p)
dev.off()