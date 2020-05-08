require(brms)
require(bayesplot)
require(tidyverse)
require(tidybayes)
require(lme4)
require(parallel)
require(patchwork)

#### Prepare data ####
core.num <- parallel::detectCores()
options(mc.cores=core.num)

my.seed <- 123
set.seed(my.seed)

matdat <- read.csv('./rtMRI-velum/velum_data.csv', header=T)

# get rid of items that gesture onsets that begin *after* the vowel offset (only 5/7152 total items)
matdat <- matdat[(matdat$velumopening_gesture_on - matdat$Vokal_off)<0,]

# separate coda contexts by alveolar voiced stop vs. alveolar voiceless stop
voiceless <- c('nt__','nt_@','nt_6','nt_a')
voiced    <- c('nd_@','nd_6','nd_a')

matdat$voicing <- c()
matdat$voicing[matdat$post %in% voiceless]  <- "voiceless"
matdat$voicing[matdat$post %in% voiced]     <- "voiced"

# only include alveolar nasal items preceding a voiced or voiceless stop consonant
coda <- c(voiceless,voiced)

# only include neutrally stressed utterances
stresses <- c("N")

# subset the data
subdat <- matdat[matdat$post %in% coda & matdat$stress %in% stresses, ]



# # # # # # # # # # # # # #
# Figure X                #
# Velum gesture duration  #
# # # # # # # # # # # # # #

The BRM for velum gesture duration was built using a log-normal distribution \citep{Rosen2005, Ratnikova2017, Gahl2019}. The following distributions were used as weakly informative priors (on the log-odds scale): for the intercept of duration (corresponding to \ips{nd} and overall mean speech rate), a normal distribution with mean 0 and standard deviation 3 (Normal(0, 3)); for the effect of voicing (when \ips{nt}) and centred speech rate, Normal(0, 1). These roughly correspond to a belief that the intercept is between 0 and 400 ms, and the duration changes (increase or decreases) by a factor of 0 to 7 in the \ips{nd} context, at 95\% confidence. For the model standard deviation and the random intercept standard deviation we used a half Cauchy distribution with location 0 and scale 0.1. For the correlation between random effects, an LKJ(2) distribution. The same prior specification was used for the other models in this section (velum gesture onset time and offset time).

# create the dependent variable
subdat$DV <- subdat$velumopening_gesture_dur*1000

# full model
dur <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(normal(0, 1), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 0.1), class = sd),
            prior(cauchy(0, 0.1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/dur",
  save_all_pars = TRUE
)

# model summary
summary(dur)

# posterior predictive check
pp_check(dur, nsamples = 100)

# reduced/null model
dur_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(cauchy(0, 0.1), class = sd),
            prior(cauchy(0, 0.1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/dur_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(dur, dur_null)

# calculate the marginal posteriors of the full model
dur_post <- brms::posterior_samples(dur, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 95% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(dur_post$DV[dur_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(dur_post$DV[dur_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 5
pdf(file="./rtMRI-velum/plots/duration.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
post_plot <- dur_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture duration") + 
  ylab("Posterior probability (density)") + xlab("Duration (ms)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(post_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
post_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(0,1000,5)) + theme_bw()
dev.off()



# # # # # # # # # # # #
# Figure X            #
# Velum gesture onset #
# # # # # # # # # # # #

# create the dependent variable
subdat$DV <- -(subdat$velumopening_gesture_on - subdat$Vokal_off)*1000

# full model
onset <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(normal(0, 1), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 0.1), class = sd),
            prior(cauchy(0, 0.1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/onset",
  save_all_pars = TRUE
)

# model summary
summary(onset)

# reduced/null model
onset_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(cauchy(0, 0.1), class = sd),
            prior(cauchy(0, 0.1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/onset_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(onset, onset_null)

# calculate the marginal posteriors of the full model
onset_post <- brms::posterior_samples(onset, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(-onset_post$DV[onset_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(-onset_post$DV[onset_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 7
pdf(file="./rtMRI-velum/plots/onset.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
post_plot <- onset_post %>%
  ggplot(aes(-DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture onset") + 
  ylab("Posterior probability (density)") + xlab("Time (ms) relative to acoustic vowel offset")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(post_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
post_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(-1000,1000,5)) + theme_bw()
dev.off()



# # # # # # # # # # # # #
# Figure X              #
# Velum gesture offset  #
# # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- (subdat$velumopening_gesture_off - subdat$Vokal_off)*1000

# full model
offset <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(normal(0, 1), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 0.1), class = sd),
            prior(cauchy(0, 0.1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/offset",
  save_all_pars = TRUE
)

# model summary
summary(offset)

# reduced/null model
offset_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(cauchy(0, 0.1), class = sd),
            prior(cauchy(0, 0.1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/offset_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(offset, offset_null)

# calculate the marginal posteriors of the full model
offset_post <- brms::posterior_samples(offset, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(offset_post$DV[offset_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(offset_post$DV[offset_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 6
pdf(file="./rtMRI-velum/plots/offset.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
post_plot <- offset_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture offset") + 
  ylab("Posterior probability (density)") + xlab("Time (ms) relative to acoustic vowel offset")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(post_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
post_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(0,1000,5)) + theme_bw()
dev.off()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Figure X                                                    #
# Integral of velum movement in vowel (area under the curve)  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- subdat$vowel.integ

# full model
integ <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = gaussian(),
  prior = c(prior(normal(0, 50), class = Intercept),
            prior(normal(0, 25), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 1), class = sd),
            prior(cauchy(0, 1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4, 
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/vowel_integ",
  save_all_pars = TRUE
)

# model summary
summary(integ)

# reduced/null model
integ_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = gaussian(),
  prior = c(prior(normal(0, 50), class = Intercept),
            prior(cauchy(0, 1), class = sd),
            prior(cauchy(0, 1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/vowel_integ_null",
  save_all_pars = TRUE
)


# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(integ, integ_null)

# calculate the marginal posteriors of the full model
integ_post <- brms::posterior_samples(integ, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(integ_post$DV[integ_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(integ_post$DV[integ_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 9
pdf(file="./rtMRI-velum/plots/vowel_integ.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
post_plot <- integ_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Integral of velum movement in vowel interval") + 
  ylab("Posterior probability (density)") + xlab("Velum movement integral (time X magnitude)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(post_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
post_plot + 
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



# # # # # # # # # # # # # # # #
# Figure 8                    #
# Velum gesture peak (timing) #
# # # # # # # # # # # # # # # #

The BRM for the time point of the velum gesture peak was built using a Gaussian distribution; unlike the other measures of timing, the time point of the velum gesture peak is not expected to follow a one-sided distribution, since the peak can potentially occur before or after the vowel offset. The following distributions were used as weakly informative priors (on the milliseconds scale): for the intercept (corresponding to \ips{nd} and overall mean speech rate), Normal(0, 200); for the effect of voicing (when \ips{nt}) and centred speech rate, Normal(0, 100). These correspond to a belief that the intercept is between -400 and 400 ms, and that the time changes by -200 to 200 ms in the \ips{nt} context, at 95\% confidence. For the model standard deviation and the random intercept standard deviation we used HalfCauchy(0, 5). For the correlation between random effects, an LKJ(2) distribution.

# create the dependent variable
subdat$DV <- (subdat$velumopening_maxcon_on - subdat$Vokal_off)*1000

# full model
gest.max <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = gaussian(),
  prior = c(prior(normal(0, 200), class = Intercept),
            prior(normal(0, 100), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 5), class = sd),
            prior(cauchy(0, 5), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4, 
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/gest_max",
  save_all_pars = TRUE
)

# model summary
summary(gest.max)

# reduced/null model
gest.max_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = gaussian(),
  prior = c(prior(normal(0, 200), class = Intercept),
            prior(cauchy(0, 5), class = sd),
            prior(cauchy(0, 5), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/gest_max_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(gest.max, gest.max_null)

# calculate the marginal posteriors of the full model
gest.max_post <- brms::posterior_samples(gest.max, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%    
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(gest.max_post$DV[gest.max_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(gest.max_post$DV[gest.max_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 8
pdf(file="./rtMRI-velum/plots/gest_max.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
post_plot <- gest.max_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture peak (timing)") + 
  ylab("Posterior probability (density)") + xlab("Time (ms) relative to acoustic vowel offset")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(post_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
post_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(-1000,1000,5)) + theme_bw()
dev.off()



# # # # # # # # # # # # # # # # #
# Figure 9                      #
# Velum gesture peak (magnitude)#
# # # # # # # # # # # # # # # # #

The BRM for the velum gesture magnitude was built using a Beta distribution, since the magnitude values are on a 0-1 scale. The following weakly informative priors were used: Normal(0, 10) for the intercept; Normal(0, 5) for voicing, speech rate, and the random effects standard deviations; the \textit{brms} default prior for the $\phi$ parameter of the beta distribution (gamma(0.01, 0.01)); and LKJ(2) prior for the random effects correlation.

# create the dependent variable
subdat$DV <- subdat$velum2US_velumopening_maxcon_onset

# full model
gest.max.mag <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = Beta(),
  prior = c(prior(normal(0, 10), class = Intercept),
            prior(normal(0, 5), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 5), class = sd),
            prior(gamma(0.01, 0.01), class = phi),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4, 
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/gest_max_mag",
  save_all_pars = TRUE
)

# model summary
summary(gest.max.mag)

# reduced/null model
gest.max.mag_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = Beta(),
  prior = c(prior(normal(0, 10), class = Intercept),
            prior(cauchy(0, 5), class = sd),
            prior(gamma(0.01, 0.01), class = phi),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/gest_max_mag_null",
  save_all_pars = TRUE
)


# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(gest.max.mag, gest.max.mag_null)

# calculate the marginal posteriors of the full model
gest.max.mag_post <- brms::posterior_samples(gest.max.mag, pars="b_") %>%
  dplyr::mutate(
    nd = plogis(b_Intercept),
    nt = plogis(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(gest.max.mag_post$DV[gest.max.mag_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(gest.max.mag_post$DV[gest.max.mag_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 9
pdf(file="./rtMRI-velum/plots/gest_max_mag.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
post_plot <- gest.max.mag_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture peak (magnitude)") + 
  ylab("Posterior probability (density)") + xlab("Velum opening magnitude (speaker-scaled)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(post_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
post_plot + 
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



# # # # # # # # # # # # # # # # # # # # # #
# Figure 11                               #
# Velum gesture plateau acceleration peak #
# # # # # # # # # # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- scale(subdat$accel.peak)

# full model
accel <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = gaussian(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(normal(0, 3), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 1), class = sd),
            prior(cauchy(0, 1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num, 
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/accel_peak",
  save_all_pars = TRUE
)

# model summary
summary(accel)

# reduced/null model
accel_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = gaussian(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(cauchy(0, 1), class = sd),
            prior(cauchy(0, 1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/accel_peak_null",
  save_all_pars = TRUE
)


# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(accel, accel_null)

# calculate the marginal posteriors of the full model
accel_post <- brms::posterior_samples(accel, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(accel_post$DV[accel_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(accel_post$DV[accel_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 9
pdf(file="./rtMRI-velum/plots/plateau_accel_peak.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
post_plot <- accel_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Maximum acceleration within gesture plateau") + 
  ylab("Posterior probability (density)") + xlab("Velum movement acceleration (normalized mm/s^2)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(post_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
post_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(-1,1,0.1)) + theme_bw()
dev.off()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Figure XX                                                   #
# Integral of velum movement in vowel (area under the curve)  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- scale(abs(subdat$quad.coef))

# full model
quad <- brms::brm(
  DV ~ voicing +
    (1 + voicing| speaker) +
    (1 + voicing | word),
  data = subdat,
  family = gaussian(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(normal(0, 3), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 1), class = sd),
            prior(cauchy(0, 1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4, 
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/quad_coef",
  save_all_pars = TRUE
)

# model summary
summary(quad)

# reduced/null model
quad_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = gaussian,
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(cauchy(0, 1), class = sd),
            prior(cauchy(0, 1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/quad_coef_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(quad, quad_null)

# calculate the marginal posteriors of the full model
quad_post <- brms::posterior_samples(quad, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(quad_post$DV[quad_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(quad_post$DV[quad_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 9
pdf(file="./rtMRI-velum/plots/quad_coef.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
post_plot <- quad_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Quadratic coefficient of second-order polynomial") + 
  ylab("Posterior probability (density)") + xlab("Quadratic coefficient (normalized)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(post_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
post_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(-1,1,0.1)) + theme_bw()
dev.off()



# # # # # # # # # # # # # # #
# Figure XX                 #
# Gesture onset stiffness   #
# # # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- subdat$stiff.ons

# full model
stiff.ons <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(normal(0, 1), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 0.1), class = sd),
            prior(cauchy(0, 0.1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4, 
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/stiff_ons",
  save_all_pars = TRUE
)

# model summary
summary(stiff.ons)

# reduced/null model
stiff.ons_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(cauchy(0, 0.1), class = sd),
            prior(cauchy(0, 0.1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/stiff_ons_null",
  save_all_pars = TRUE
)


# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(stiff.ons, stiff.ons_null)

# calculate the marginal posteriors of the full model
stiff.ons_post <- brms::posterior_samples(stiff.ons, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(stiff.ons_post$DV[stiff.ons_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(stiff.ons_post$DV[stiff.ons_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 9
pdf(file="./rtMRI-velum/plots/stiff_ons.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
post_plot <- stiff.ons_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Stiffness of the velum gesture onset") + 
  ylab("Posterior probability (density)") + xlab("Stiffness")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(post_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
post_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(-1,1,0.00125)) + theme_bw()
dev.off()



# # # # # # # # # # # # # # #
# Figure XX                 #
# Gesture offset stiffness  #
# # # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- subdat$stiff.off

# full model
stiff.off <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(normal(0, 1), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 0.1), class = sd),
            prior(cauchy(0, 0.1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4, 
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/stiff_off",
  save_all_pars = TRUE
)

# model summary
summary(stiff.off)

# reduced/null model
stiff.off_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(cauchy(0, 0.1), class = sd),
            prior(cauchy(0, 0.1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/stiff_off_null",
  save_all_pars = TRUE
)


# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(stiff.off, stiff.off_null)

# calculate the marginal posteriors of the full model
stiff.off_post <- brms::posterior_samples(stiff.off, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(stiff.off_post$DV[stiff.off_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(stiff.off_post$DV[stiff.off_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 9
pdf(file="./rtMRI-velum/plots/stiff_off.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
post_plot <- stiff.off_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Stiffness of the velum gesture offset") + 
  ylab("Posterior probability (density)") + xlab("Stiffness")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(post_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
post_plot + 
  # whisker for context 1
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1], xend=nd.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2], xend=nd.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[1]) +
  # whisker for context 2
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1], xend=nt.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2], xend=nt.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[2]) +
  scale_x_continuous(breaks=seq(-1,1,0.00125)) + theme_bw()
dev.off()




# Composite marginal post plot ####
my.cols <- c("#2c7fb8","#7fcdbb")

time_post <- bind_rows(
  dur_post %>% rename(`duration` = DV) %>% pivot_longer(`duration`, names_to = "outcome", values_to = "estimate"),
  offset_post %>% rename(`offset` = DV) %>% pivot_longer(`offset`, names_to = "outcome", values_to = "estimate"),
  onset_post %>% rename(`onset` = DV) %>% pivot_longer(`onset`, names_to = "outcome", values_to = "estimate"),
  gest.max_post %>% rename(`peak` = DV) %>% pivot_longer(`peak`, names_to = "outcome", values_to = "estimate")
) %>%
  mutate(outcome = factor(outcome, levels = c("duration", "onset", "offset", "peak")))

time_post$estimate[time_post$outcome=="onset"] <- -time_post$estimate[time_post$outcome=="onset"]

time_post$outcome <- factor(time_post$outcome, levels = c("offset","peak","onset","duration"))


# separate plots

p1 <- time_post %>%
  filter(outcome == "duration") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,500,5)) +
  coord_cartesian(ylim = c(1.3, 1.5), xlim = c(250,327)) +
  scale_fill_manual(values = my.cols) +
  labs(x = "Duration (ms)", y = element_blank()) + theme_bw() + theme(legend.position = "none",axis.text=element_text(size=12),axis.title=element_text(size=13))

p2 <- time_post %>%
  filter(outcome %in% c("onset","peak","offset")) %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  geom_vline(xintercept = 0, lty=2) +
  scale_x_continuous(breaks=seq(-200,300,20)) +
  coord_cartesian(xlim = c(-125,205)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Time (ms) relative to acoustic vowel offset", y = element_blank(), fill = "Context") + theme_bw() + theme(axis.text=element_text(size=12),axis.title=element_text(size=13),panel.grid.major.y=element_blank())

pdf(file="./rtMRI-velum/plots/time_plots.pdf",width=9.5,height=5,onefile=T,pointsize=14)
p1 + p2 + patchwork::plot_layout(ncol = 1, guides = "collect") + theme(legend.position = "right")
dev.off()



# Composite marginal post plot ####

time_post <- bind_rows(
  gest.max.mag_post %>% rename(`peak    \nmagnitude` = DV) %>% pivot_longer(`peak    \nmagnitude`, names_to = "outcome", values_to = "estimate"),
  integ_post %>% rename(`gesture\nintegral` = DV) %>% pivot_longer(`gesture\nintegral`, names_to = "outcome", values_to = "estimate"),
  accel_post %>% rename(`maximum  \nacceleration` = DV) %>% pivot_longer(`maximum  \nacceleration`, names_to = "outcome", values_to = "estimate"),
  quad_post %>% rename(`quadratic \ncoefficient` = DV) %>% pivot_longer(`quadratic \ncoefficient`, names_to = "outcome", values_to = "estimate"),
  stiff.ons_post %>% rename(`opening\nstiffness` = DV) %>% pivot_longer(`opening\nstiffness`, names_to = "outcome", values_to = "estimate"),
  stiff.off_post %>% rename(`closing \nstiffness` = DV) %>% pivot_longer(`closing \nstiffness`, names_to = "outcome", values_to = "estimate")
) %>%
  mutate(outcome = factor(outcome, levels = c("peak    \nmagnitude", "gesture\nintegral", "maximum  \nacceleration", "quadratic \ncoefficient", "closing \nstiffness", "opening\nstiffness")))



# separate plots

p1 <- time_post %>%
  filter(outcome == "peak    \nmagnitude") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,1,0.025)) +
  coord_cartesian(xlim = c(0.58,0.775), ylim = c(1.4,1.4)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Velum opening magnitude (speaker-scaled)", y = element_blank(), fill = "Context") + theme_bw() + theme(axis.text=element_text(size=12),axis.title=element_text(size=13))

p2 <- time_post %>%
  filter(outcome == "gesture\nintegral") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(10,50,1)) +
  coord_cartesian(xlim = c(18,32), ylim = c(1.4,1.4)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Velum displacement integral (time X magnitude)", y = element_blank(), fill = "Context") + theme_bw() + theme(axis.text=element_text(size=12),axis.title=element_text(size=13))

p3 <- time_post %>%
  filter(outcome %in% c("maximum  \nacceleration","quadratic \ncoefficient")) %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(-1,1,0.1)) +
  coord_cartesian(xlim = c(-0.5,0.45), ylim = c(1.5,1.95)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "z-scale normalized values", y = element_blank(), fill = "Context") + theme_bw() + theme(legend.position = "none",axis.text=element_text(size=12),axis.title=element_text(size=13),panel.grid.major.y=element_blank())

p4 <- time_post %>%
  filter(outcome %in% c("opening\nstiffness","closing \nstiffness")) %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,30,0.5)) +
  coord_cartesian(xlim = c(12,17.5), ylim = c(1.4,2.3)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Stiffness", y = element_blank(), fill = "Context") + theme_bw() + theme(axis.text=element_text(size=12),axis.title=element_text(size=13),panel.grid.major.y=element_blank())

pdf(file="./rtMRI-velum/plots/shape_plots.pdf",width=9,height=9,onefile=T,pointsize=14)
p1 + p2 + p3 + p4 + patchwork::plot_layout(ncol = 1, guides = "collect") + theme(legend.position = "right")
dev.off()


pdf(file="./rtMRI-velum/plots/magnitude_plot.pdf",width=9,height=3,onefile=T,pointsize=14)
p1 + patchwork::plot_layout(ncol = 1, guides = "collect") + theme(legend.position = "right")
dev.off()

pdf(file="./rtMRI-velum/plots/integral_plot.pdf",width=9,height=3,onefile=T,pointsize=14)
p2 + patchwork::plot_layout(ncol = 1, guides = "collect") + theme(legend.position = "right")
dev.off()

pdf(file="./rtMRI-velum/plots/peakiness_plots.pdf",width=9,height=4,onefile=T,pointsize=14)
p3 + patchwork::plot_layout(ncol = 1, guides = "collect") + theme(legend.position = "right")
dev.off()


pdf(file="./rtMRI-velum/plots/stiffness_plots.pdf",width=9,height=4,onefile=T,pointsize=14)
p4 + patchwork::plot_layout(ncol = 1, guides = "collect") + theme(legend.position = "right")
dev.off()
