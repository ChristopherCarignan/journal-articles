require(brms)
require(tidyverse)
require(lme4)
require(MuMIn)
require(parallel)

# setwd("/home/chris/Documents/research/journal-articles/rtMRI-velum/")

#### Prepare data ####
options(mc.cores=parallel::detectCores())

my.seed <- 123
set.seed(my.seed)

matdat <- read.csv('./rtMRI-velum/velum_data.csv', header=T)
matdat$word <- paste0(matdat$prev,matdat$vowel)

# separate coda contexts by alveolar voiced stop vs. alveolar voiceless stop
voiceless <- c('nt__','nt_@','nt_6','nt_a')
voiced    <- c('nd_@','nd_6','nd_a')

# only include alveolar nasal items preceding a voiced or voiceless stop consonant
coda <- c(voiceless,voiced)

# only include neutrally stressed utterances
stresses <- c("N")

# subset the data
subdat <- matdat[matdat$post %in% coda & matdat$stress %in% stresses, ]
subdat$voicing <- c()
subdat$voicing[subdat$post %in% voiceless] <- "voiceless"
subdat$voicing[subdat$post %in% voiced] <- "voiced"



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 3                                                        #
# Relation between N duration and duration of vowel nasalization  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# duration of nasal consonant: interval between vowel offset and velum gesture offset
subdat$dur.NN <- 1000*(subdat$velumopening_maxvel_off-subdat$Vokal_off)
# duration of vowel nasalization: interval between velum gesture onset and vowel offset
subdat$dur.VN <- 1000*(subdat$Vokal_off-subdat$velumopening_maxvel_on)

# subset the data by voicing context
nt.dat <- subdat[subdat$voicing=="voiceless",]
nd.dat <- subdat[subdat$voicing=="voiced",]


## Estimating R^2 of the fixed effect (all data combined)
# first create the full model and reduced model (without fixed effect)
mod.all     <- lme4::lmer(dur.VN ~ dur.NN + (1|speaker), data=subdat) # full model
mod.all.red <- lme4::lmer(dur.VN ~ (1|speaker), data=subdat) # reduced model

# compare the two models
anova(mod.all,mod.all.red)

# partition the R^2 of the fixed effect
r2.all  <- MuMIn::r.squaredGLMM(mod.all)[1] # R^2
int.all <- as.numeric(coef(lm(dur.VN ~ dur.NN, data=subdat))[1]) # intercept

# get ranges of data (for plotting)
all.xlims <- range(subdat$dur.NN)
all.ylims <- all.xlims*-sqrt(r2.all) + int.all


## Estimating R^2 of the fixed effect (only /nt/ data)
# first create the full model and reduced model (without fixed effect)
mod.nt      <- lmer(dur.VN ~ dur.NN + (1|speaker), data=nt.dat) # full model
mod.nt.red  <- lmer(dur.VN ~ (1|speaker), data=nt.dat) # reduced model

# compare the two models
anova(mod.nt,mod.nt.red)

# partition the R^2 of the fixed effect
r2.nt  <- r.squaredGLMM(mod.nt)[1] # R^2
int.nt <- as.numeric(coef(lm(dur.VN ~ dur.NN, data=nt.dat))[1]) # intercept

# get ranges of data (for plotting)
nt.xlims <- range(nt.dat$dur.N)
nt.ylims <- nt.xlims*-sqrt(r2.nt) + int.nt

## Estimating R^2 of the fixed effect (only /nt/ data)
# first create the full model and reduced model (without fixed effect)
mod.nd      <- lmer(dur.VN ~ dur.NN + (1|speaker), data=nd.dat) # full model
mod.nd.red  <- lmer(dur.VN ~ (1|speaker), data=nd.dat) # reduced model

# compare the two models
anova(mod.nd,mod.nd.red)

# partition the R^2 of the fixed effect
r2.nd  <- r.squaredGLMM(mod.nd)[1] # R^2
int.nd <- as.numeric(coef(lm(dur.VN ~ dur.NN, data=nd.dat))[1]) # intercept

# get ranges of data (for plotting)
nd.xlims <- range(nd.dat$dur.NN)
nd.ylims <- nd.xlims*-sqrt(r2.nd) + int.nd


## create Figure 3
pdf(file="rtMRI-velum/plots/trading_relation.pdf",width=5,height=3,onefile=T,pointsize=16)
ggplot(subdat,aes(x=dur.NN, y=dur.VN, col=voicing, shape=voicing)) + 
  geom_point(alpha=0.7,cex=2) + scale_shape_manual(values=c(2,4)) +
  geom_segment(x=nt.xlims[1],xend=nt.xlims[2],y=nt.ylims[1],yend=nt.ylims[2],lty=2,lwd=0.6,col='black') +
  geom_segment(x=nd.xlims[1],xend=nd.xlims[2],y=nd.ylims[1],yend=nd.ylims[2],lty=3,lwd=0.6,col='black') +
  scale_x_continuous(breaks=seq(-100,300,50)) +
  scale_y_continuous(breaks=seq(-200,300,50)) +
  xlab("Duration of consonant nasalization (ms)") + ylab("Duration of vowel nasalization (ms)") + theme_bw()
dev.off()



# # # # # # # # # # # # # # # # # # # # # # #
# Bayesian regression models (Figures 5-11) #
# # # # # # # # # # # # # # # # # # # # # # #

# we will use these same priors for all of the BRMs in Figures 5-9

# full model priors:
priors <- c(
  prior(normal(0, 100), class = Intercept),
  prior(normal(0, 100), class = b, coef = voicingvoiceless),
  prior(cauchy(0, 10), class = sd),
  prior(cauchy(0, 10), class = sigma),
  prior(lkj(2), class = cor)
)

# null model priors:
priors_null <- c(
  prior(normal(0, 100), class = Intercept),
  prior(cauchy(0, 10), class = sd),
  prior(cauchy(0, 10), class = sigma)
)


# # # # # # # # # # # # # #
# Figure 5                #
# Velum gesture duration  #
# # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- subdat$velumopening_maxvel_dur*1000

# full model
m1 <- brms::brm(
  formula = DV ~ voicing + (1+voicing|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors,
  file = "models/dur",
  save_all_pars = TRUE
)

# model summary
summary(m1)

# reduced/null model
m1_null <- brms::brm(
  formula = DV ~ (1|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors_null,
  file = "./rtMRI-velum/models/dur_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(m1, m1_null)

# calculate the marginal posteriors of the full model
m1_post <- brms::posterior_samples(m1, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 90% CI of /nd/ items
nd.ci <- quantile(m1_post$DV[m1_post$context=="nd"], probs=c(0.05,0.95))
# 90% CI of /nt/ items
nt.ci <- quantile(m1_post$DV[m1_post$context=="nt"], probs=c(0.05,0.95))

## create Figure 5
pdf(file="rtMRI-velum/plots/duration.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m1_plot <- m1_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture duration") + 
  ylab("Posterior probability (density)") + xlab("Duration (ms)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m1_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m1_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(100,200,5)) + theme_bw()
dev.off()



# # # # # # # # # # # # #
# Figure 6              #
# Velum gesture offset  #
# # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- (subdat$velumopening_maxvel_off - subdat$Vokal_off)*1000

# full model
m2 <- brms::brm(
  formula = DV ~ voicing + (1+voicing|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors,
  file = "./rtMRI-velum/models/offset",
  save_all_pars = TRUE
)

# model summary
summary(m2)

# reduced/null model
m2_null <- brms::brm(
  formula = DV ~ (1|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors_null,
  file = "./rtMRI-velum/models/offset_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(m2, m2_null)

# calculate the marginal posteriors of the full model
m2_post <- brms::posterior_samples(m2, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 90% CI of /nd/ items
nd.ci <- quantile(m2_post$DV[m2_post$context=="nd"], probs=c(0.05,0.95))
# 90% CI of /nt/ items
nt.ci <- quantile(m2_post$DV[m2_post$context=="nt"], probs=c(0.05,0.95))

## create Figure 6
pdf(file="rtMRI-velum/plots/offset.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m2_plot <- m2_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture offset") + 
  ylab("Posterior probability (density)") + xlab("Time (ms) relative to acoustic vowel offset")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m2_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m2_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(50,150,5)) + theme_bw()
dev.off()



# # # # # # # # # # # # #
# Figure 7              #
# Velum gesture offset  #
# # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- (subdat$velumopening_maxvel_on - subdat$Vokal_off)*1000

# full model
m3 <- brms::brm(
  formula = DV ~ voicing + (1+voicing|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors,
  file = "./rtMRI-velum/models/onset",
  save_all_pars = TRUE
)

# model summary
summary(m3)

# reduced/null model
m3_null <- brms::brm(
  formula = DV ~ (1|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors_null,
  file = "./rtMRI-velum/models/onset_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(m3, m3_null)

# calculate the marginal posteriors of the full model
m3_post <- brms::posterior_samples(m3, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 90% CI of /nd/ items
nd.ci <- quantile(m3_post$DV[m3_post$context=="nd"], probs=c(0.05,0.95))
# 90% CI of /nt/ items
nt.ci <- quantile(m3_post$DV[m3_post$context=="nt"], probs=c(0.05,0.95))

## create Figure 7
pdf(file="rtMRI-velum/plots/onset.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m3_plot <- m3_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture onset") + 
  ylab("Posterior probability (density)") + xlab("Time (ms) relative to acoustic vowel offset")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m3_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m3_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(-100,0,5)) + theme_bw()
dev.off()



# # # # # # # # # # # # # # # #
# Figure 8                    #
# Velum gesture peak (timing) #
# # # # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- (subdat$velumopening_maxcon_on - subdat$Vokal_off)*1000

# full model
m4 <- brms::brm(
  formula = DV ~ voicing + (1+voicing|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors,
  file = "./rtMRI-velum/models/gest_max",
  save_all_pars = TRUE
)

# model summary
summary(m4)

# reduced/null model
m4_null <- brms::brm(
  formula = DV ~ (1|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors_null,
  file = "./rtMRI-velum/models/gest_max_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(m4, m4_null)

# calculate the marginal posteriors of the full model
m4_post <- brms::posterior_samples(m4, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 90% CI of /nd/ items
nd.ci <- quantile(m4_post$DV[m4_post$context=="nd"], probs=c(0.05,0.95))
# 90% CI of /nt/ items
nt.ci <- quantile(m4_post$DV[m4_post$context=="nt"], probs=c(0.05,0.95))

## create Figure 8
pdf(file="rtMRI-velum/plots/gest_max.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m4_plot <- m4_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture peak (timing)") + 
  ylab("Posterior probability (density)") + xlab("Time (ms) relative to acoustic vowel offset")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m4_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m4_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(0,100,5)) + theme_bw()
dev.off()



# # # # # # # # # # # # # # # # #
# Figure 9                      #
# Velum gesture peak (magnitude)#
# # # # # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- subdat$velum2US_velumopening_maxcon_onset

# full model
m5 <- brms::brm(
  formula = DV ~ voicing + (1+voicing|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors,
  file = "./rtMRI-velum/models/gest_max_mag",
  save_all_pars = TRUE
)

# model summary
summary(m5)

# reduced/null model
m5_null <- brms::brm(
  formula = DV ~ (1|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors_null,
  file = "./rtMRI-velum/models/gest_max_mag_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(m5, m5_null)

# calculate the marginal posteriors of the full model
m5_post <- brms::posterior_samples(m5, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 90% CI of /nd/ items
nd.ci <- quantile(m5_post$DV[m5_post$context=="nd"], probs=c(0.05,0.95))
# 90% CI of /nt/ items
nt.ci <- quantile(m5_post$DV[m5_post$context=="nt"], probs=c(0.05,0.95))

## create Figure 9
pdf(file="rtMRI-velum/plots/gest_max_mag.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m5_plot <- m5_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture peak (magnitude)") + 
  ylab("Posterior probability (density)") + xlab("Velum opening magnitude (speaker-scaled)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m5_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m5_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(0,1,0.05)) + theme_bw()
dev.off()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Figure 10                                               #
# Effect of speech rate on duration of vowel nasalization #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# we will rename the priors for Figures 10 and 11, but we're using the same prior values as from Figures 5-9

# full model priors:
priors <- c(
  prior(normal(0, 100), class = Intercept),
  prior(normal(0, 100), class = b, coef = stressF),
  prior(cauchy(0, 10), class = sd),
  prior(cauchy(0, 10), class = sigma),
  prior(lkj(2), class = cor)
)


# include both neutral/normal and fast speech utterances
stresses <- c("N","F")

# subset the data
subdat <- matdat[matdat$post %in% coda & matdat$stress %in% stresses, ]
subdat$voicing <- c()
subdat$voicing[subdat$post %in% voiceless] <- "voiceless"
subdat$voicing[subdat$post %in% voiced] <- "voiced"

# only include /nt/ items
subdat <- subdat[subdat$voicing=="voiceless",]

# reorder the speech rate factor order, so that the normal speech rate is the baseline
subdat$stress <- factor(subdat$stress, levels=c("N","F"))

# duration of nasal consonant: interval between vowel offset and velum gesture offset
subdat$dur.NN <- 1000*(subdat$velumopening_maxvel_off - subdat$Vokal_off)
# duration of vowel nasalization: interval between velum gesture onset and vowel offset
subdat$dur.VN <- 1000*(subdat$Vokal_off - subdat$velumopening_maxvel_on)


# first, we look at the duration of the vowel nasalization
subdat$DV <- subdat$dur.VN

# full model
m6 <- brms::brm(
  formula = DV ~ stress + (1+stress|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors,
  file = "./rtMRI-velum/models/rate_v_nasal",
  save_all_pars = TRUE
)

# model summary
summary(m6)

# reduced/null model
m6_null <- brms::brm(
  formula = DV ~ (1|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors_null,
  file = "./rtMRI-velum/models/rate_v_nasal_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(m6, m6_null)

# calculate the marginal posteriors of the full model
m6_post <- brms::posterior_samples(m6, pars="b_") %>%
  dplyr::mutate(
    normal = b_Intercept,
    fast = b_Intercept + b_stressF
  ) %>%
  dplyr::select(normal, fast) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 90% CI of normal rate items
norm.ci <- quantile(m6_post$DV[m6_post$context=="normal"], probs=c(0.05,0.95))
# 90% CI of fast rate items
fast.ci <- quantile(m6_post$DV[m6_post$context=="fast"], probs=c(0.05,0.95))

## create Figure 10
pdf(file="rtMRI-velum/plots/speech_rate_v_nasal.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m6_plot <- m6_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Effect of speech rate on duration of vowel nasalization") + 
  ylab("Posterior probability (density)") + xlab("Duration (ms)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m6_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m6_plot + 
  # whisker for context 1
  geom_segment(aes(x=norm.ci[1], xend=norm.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=norm.ci[1], xend=norm.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=norm.ci[2], xend=norm.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=fast.ci[1], xend=fast.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=fast.ci[1], xend=fast.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=fast.ci[2], xend=fast.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(0,100,5)) + theme_bw()
dev.off()


# now, we look at the duration of the consonant nasalization
subdat$DV <- subdat$dur.NN

# full model
m7 <- brms::brm(
  formula = DV ~ stress + (1+stress|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors,
  file = "./rtMRI-velum/models/rate_n_nasal",
  save_all_pars = TRUE
)

# model summary
summary(m7)

# reduced/null model
m7_null <- brms::brm(
  formula = DV ~ (1|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors_null,
  file = "./rtMRI-velum/models/rate_n_nasal_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(m7, m7_null)

# calculate the marginal posteriors of the full model
m7_post <- brms::posterior_samples(m7, pars="b_") %>%
  dplyr::mutate(
    normal = b_Intercept,
    fast = b_Intercept + b_stressF
  ) %>%
  dplyr::select(normal, fast) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 90% CI of normal rate items
norm.ci <- quantile(m7_post$DV[m7_post$context=="normal"], probs=c(0.05,0.95))
# 90% CI of fast rate items
fast.ci <- quantile(m7_post$DV[m7_post$context=="fast"], probs=c(0.05,0.95))


## create Figure 11
pdf(file="rtMRI-velum/plots/speech_rate_n_nasal.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m7_plot <- m7_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Effect of speech rate on duration of consonant nasalization") + 
  ylab("Posterior probability (density)") + xlab("Duration (ms)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m7_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m7_plot + 
  # whisker for context 1
  geom_segment(aes(x=norm.ci[1], xend=norm.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=norm.ci[1], xend=norm.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=norm.ci[2], xend=norm.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=fast.ci[1], xend=fast.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=fast.ci[1], xend=fast.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=fast.ci[2], xend=fast.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(0,100,5)) + theme_bw()
dev.off()
