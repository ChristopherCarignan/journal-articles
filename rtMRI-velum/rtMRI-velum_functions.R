require(brms)
require(bayesplot)
require(tidyverse)
require(lme4)
require(parallel)

#### Prepare data ####
core.num <- parallel::detectCores()
options(mc.cores=core.num)

my.seed <- 123
set.seed(my.seed)

matdat <- read.csv('./rtMRI-velum/velum_data.csv', header=T)

# get rid of items that gesture onsets that begin *after* the vowel offset (only 5/7152 total items)
matdat <- matdat[(matdat$velumopening_gesture_on - matdat$Vokal_off)<0,]

# Center the speech rate
matdat <- matdat %>% mutate(speech_rate_c = sp_rate - mean(sp_rate))

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



# # # # # # # # # # # # # # # # # # # # # # #
# Bayesian regression models (Figures 5-11) #
# # # # # # # # # # # # # # # # # # # # # # #

# we will use these same priors for all of the BRMs in Figures 5-9

# full model priors:
priors <- c(
  prior(normal(0, 3), class = Intercept),
  prior(normal(0, 1), class = b, coef = voicingvoiceless),
  prior(normal(0, 1), class = b, coef = speech_rate_c),
  prior(cauchy(0, 0.1), class = sd),
  prior(cauchy(0, 0.1), class = sigma),
  prior(lkj(2), class = cor)
)

# null model priors:
priors_null <- c(
  prior(normal(0, 3), class = Intercept),
  prior(normal(0, 1), class = b, coef = speech_rate_c),
  prior(cauchy(0, 0.1), class = sd),
  prior(cauchy(0, 0.1), class = sigma),
  prior(lkj(2), class = cor)
)


# # # # # # # # # # # # # #
# Figure 5                #
# Velum gesture duration  #
# # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- subdat$velumopening_gesture_dur*1000

# full model
m1 <- brms::brm(
  DV ~ voicing +
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = lognormal(),
  prior = priors,
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
summary(m1)

# posterior predictive check
pp_check(m1, nsamples = 100)

# reduced/null model
m1_null <- brms::brm(
  DV ~
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = lognormal(),
  prior = priors_null,
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
#brms::bayes_factor(m1, m1_null)

# calculate the marginal posteriors of the full model
m1_post <- brms::posterior_samples(m1, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 95% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(m1_post$DV[m1_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(m1_post$DV[m1_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 5
pdf(file="./rtMRI-velum/plots/duration.pdf",width=7,height=3,onefile=T,pointsize=16)
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
  scale_x_continuous(breaks=seq(0,1000,5)) + theme_bw()
dev.off()

# sensitivity analysis (Betancourt 2018)

m1_fixed <- fixef(m1) %>% as_tibble(rownames = "term")

m1_fixed %>%
  mutate(
    theta = c(0, 0, 0),
    sigma_prior = c(3, 1, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250)


# # # # # # # # # # # # #
# Figure 6              #
# Velum gesture offset  #
# # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- (subdat$velumopening_gesture_off - subdat$Vokal_off)*1000

# full model
m2 <- brms::brm(
  DV ~ voicing +
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = lognormal(),
  prior = priors,
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
summary(m2)

# posterior predictive check
pp_check(m2, nsamples = 100)

# reduced/null model
m2_null <- brms::brm(
  DV ~
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = lognormal(),
  prior = priors_null,
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
#brms::bayes_factor(m2, m2_null)

# calculate the marginal posteriors of the full model
m2_post <- brms::posterior_samples(m2, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(m2_post$DV[m2_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(m2_post$DV[m2_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 6
pdf(file="./rtMRI-velum/plots/offset.pdf",width=7,height=3,onefile=T,pointsize=16)
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
  scale_x_continuous(breaks=seq(0,1000,5)) + theme_bw()
dev.off()

# sensitivity analysis (Betancourt 2018)

m2_fixed <- fixef(m2) %>% as_tibble(rownames = "term")

m2_fixed %>%
  mutate(
    theta = c(0, 0, 0),
    sigma_prior = c(3, 1, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 200)

# # # # # # # # # # # #
# Figure 7            #
# Velum gesture onset #
# # # # # # # # # # # #

# create the dependent variable
subdat$DV <- -(subdat$velumopening_gesture_on - subdat$Vokal_off)*1000

# full model
m3 <- brms::brm(
  DV ~ voicing +
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = lognormal(),
  prior = priors,
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
summary(m3)

# posterior predictive check
pp_check(m3, nsamples = 100)

# reduced/null model
m3_null <- brms::brm(
  DV ~
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = lognormal(),
  prior = priors_null,
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
#brms::bayes_factor(m3, m3_null)

# calculate the marginal posteriors of the full model
m3_post <- brms::posterior_samples(m3, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(-m3_post$DV[m3_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(-m3_post$DV[m3_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 7
pdf(file="./rtMRI-velum/plots/onset.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m3_plot <- m3_post %>%
  ggplot(aes(-DV, fill = context)) +
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
  scale_x_continuous(breaks=seq(-1000,1000,5)) + theme_bw()
dev.off()

# sensitivity analysis (Betancourt 2018)

m3_fixed <- fixef(m3) %>% as_tibble(rownames = "term")

m3_fixed %>%
  mutate(
    theta = c(0, 0, 0),
    sigma_prior = c(3, 1, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 100)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Figure X                                                    #
# Integral of velum movement in vowel (area under the curve)  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# full model priors:
priors <- c(
  prior(normal(0, 50), class = Intercept),
  prior(normal(0, 25), class = b, coef = voicingvoiceless),
  prior(normal(0, 25), class = b, coef = speech_rate_c),
  prior(cauchy(0, 1), class = sd),
  prior(cauchy(0, 1), class = sigma),
  prior(lkj(2), class = cor)
)

# null model priors:
priors_null <- c(
  prior(normal(0, 50), class = Intercept),
  prior(normal(0, 25), class = b, coef = speech_rate_c),
  prior(cauchy(0, 1), class = sd),
  prior(cauchy(0, 1), class = sigma),
  prior(lkj(2), class = cor)
)

# create the dependent variable
subdat$DV <- subdat$vowel.integ

# full model
m4 <- brms::brm(
  DV ~ voicing +
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = gaussian,
  prior = priors,
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
summary(m4)

# posterior predictive check
pp_check(m4, nsamples = 100)

# reduced/null model
m4_null <- brms::brm(
  DV ~
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = gaussian,
  prior = priors_null,
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
#brms::bayes_factor(m9, m9_null)

# calculate the marginal posteriors of the full model
m4_post <- brms::posterior_samples(m4, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(m4_post$DV[m4_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(m4_post$DV[m4_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 9
pdf(file="./rtMRI-velum/plots/vowel_integ.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m4_plot <- m4_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Integral of velum movement in vowel interval") + 
  ylab("Posterior probability (density)") + xlab("Velum movement integral (time X magnitude)")
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

# sensitivity analysis (Betancourt 2018)

m4_fixed <- fixef(m4) %>% as_tibble(rownames = "term")

m4_fixed %>%
  mutate(
    theta = c(0, 0, 0),
    sigma_prior = c(50, 25, 25),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 15)

# # # # # # # # # # # # # # # #
# Figure 8                    #
# Velum gesture peak (timing) #
# # # # # # # # # # # # # # # #

# full model priors:
priors <- c(
  prior(normal(0, 200), class = Intercept),
  prior(normal(0, 100), class = b, coef = voicingvoiceless),
  prior(normal(0, 100), class = b, coef = speech_rate_c),
  prior(cauchy(0, 5), class = sd),
  prior(cauchy(0, 5), class = sigma),
  prior(lkj(2), class = cor)
)

# null model priors:
priors_null <- c(
  prior(normal(0, 200), class = Intercept),
  prior(normal(0, 100), class = b, coef = speech_rate_c),
  prior(cauchy(0, 5), class = sd),
  prior(cauchy(0, 5), class = sigma),
  prior(lkj(2), class = cor)
)

# create the dependent variable
subdat$DV <- (subdat$velumopening_maxcon_on - subdat$Vokal_off)*1000

# full model
m5 <- brms::brm(
  DV ~ voicing +
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = gaussian(),
  prior = priors,
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
summary(m5)

# posterior predictive check
pp_check(m5, nsamples = 100)

# reduced/null model
m5_null <- brms::brm(
  DV ~
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = gaussian(),
  prior = priors_null,
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
#brms::bayes_factor(m5, m5_null)

# calculate the marginal posteriors of the full model
m5_post <- brms::posterior_samples(m5, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%    
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(m5_post$DV[m5_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(m5_post$DV[m5_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 8
pdf(file="./rtMRI-velum/plots/gest_max.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m5_plot <- m5_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture peak (timing)") + 
  ylab("Posterior probability (density)") + xlab("Time (ms) relative to acoustic vowel offset")
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
  scale_x_continuous(breaks=seq(-1000,1000,5)) + theme_bw()
dev.off()

# sensitivity analysis (Betancourt 2018)

m5_fixed <- fixef(m5) %>% as_tibble(rownames = "term")

m5_fixed %>%
  mutate(
    theta = c(0, 0, 0),
    sigma_prior = c(200, 100, 100),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 15)


# # # # # # # # # # # # # # # # #
# Figure 9                      #
# Velum gesture peak (magnitude)#
# # # # # # # # # # # # # # # # #

# full model priors:
priors <- c(
  prior(normal(0, 10), class = Intercept),
  prior(normal(0, 5), class = b, coef = voicingvoiceless),
  prior(normal(0, 5), class = b, coef = speech_rate_c),
  prior(cauchy(0, 5), class = sd),
  prior(gamma(0.01, 0.01), class = phi),
  prior(lkj(2), class = cor)
)

# null model priors:
priors_null <- c(
  prior(normal(0, 10), class = Intercept),
  prior(normal(0, 5), class = b, coef = speech_rate_c),
  prior(cauchy(0, 5), class = sd),
  prior(gamma(0.01, 0.01), class = phi),
  prior(lkj(2), class = cor)
)

# create the dependent variable
subdat$DV <- subdat$velum2US_velumopening_maxcon_onset

# full model
m6 <- brms::brm(
  DV ~ voicing +
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = Beta(),
  prior = priors,
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
summary(m6)

# posterior predictive check
pp_check(m6, nsamples = 100)

# reduced/null model
m6_null <- brms::brm(
  DV ~
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = Beta(),
  prior = priors_null,
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
#brms::bayes_factor(m6, m6_null)

# calculate the marginal posteriors of the full model
m6_post <- brms::posterior_samples(m6, pars="b_") %>%
  dplyr::mutate(
    nd = plogis(b_Intercept),
    nt = plogis(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(m6_post$DV[m6_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(m6_post$DV[m6_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 9
pdf(file="./rtMRI-velum/plots/gest_max_mag.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m6_plot <- m6_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture peak (magnitude)") + 
  ylab("Posterior probability (density)") + xlab("Velum opening magnitude (speaker-scaled)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m6_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m6_plot + 
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

# sensitivity analysis (Betancourt 2018)

m6_fixed <- fixef(m6) %>% as_tibble(rownames = "term")

m6_fixed %>%
  mutate(
    theta = c(0, 0, 0),
    sigma_prior = c(10, 5, 5),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 10)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Figure 10                                               #
# Effect of speech rate on duration of vowel nasalization #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# we will rename the priors for Figures 10 and 11, but we're using the same prior values as from Figures 5-9

# full model priors:
priors <- c(
  prior(normal(0, 3), class = Intercept),
  prior(normal(0, 1), class = b, coef = stressN),
  prior(normal(0, 1), class = b, coef = speech_rate_c),
  prior(cauchy(0, 0.1), class = sd),
  prior(cauchy(0, 0.1), class = sigma),
  prior(lkj(2), class = cor)
)

# null model priors:
priors_null <- c(
  prior(normal(0, 3), class = Intercept),
  prior(normal(0, 1), class = b, coef = speech_rate_c),
  prior(cauchy(0, 0.1), class = sd),
  prior(cauchy(0, 0.1), class = sigma),
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

# duration of nasal consonant: interval between vowel offset and velum gesture offset
subdat$dur.NN <- 1000*(subdat$velumopening_gesture_off - subdat$Vokal_off)
# duration of vowel nasalization: interval between velum gesture onset and vowel offset
subdat$dur.VN <- 1000*(subdat$Vokal_off - subdat$velumopening_gesture_on)


# first, we look at the duration of the vowel nasalization
subdat$DV <- subdat$dur.VN

# full model
m7 <- brms::brm(
  DV ~ stress +
    speech_rate_c +
    (1 + stress + speech_rate_c | speaker) +
    (1 + stress + speech_rate_c | word),
  data = subdat,
  family = lognormal(),
  prior = priors,
  seed = my.seed,
  iter = 10000,
  warmup = 5000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/rate_v_nasal",
  save_all_pars = TRUE
)


# model summary
summary(m7)

# posterior predictive check
pp_check(m7, nsamples = 100)

# reduced/null model
m7_null <- brms::brm(
  DV ~
    speech_rate_c +
    (1 + stress + speech_rate_c | speaker) +
    (1 + stress + speech_rate_c | word),
  data = subdat,
  family = lognormal(),
  prior = priors_null,
  seed = my.seed,
  iter = 10000,
  warmup = 5000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/rate_v_nasal_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(m7, m7_null)

# calculate the marginal posteriors of the full model
m7_post <- brms::posterior_samples(m7, pars="b_") %>%
  dplyr::mutate(
    fast = exp(b_Intercept),
    normal = exp(b_Intercept) * exp(b_stressN)
  ) %>%
  dplyr::select(normal, fast) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of normal rate items
norm.ci <- quantile(m7_post$DV[m7_post$context=="normal"], probs=c(0.025,0.975))
# 95% CI of fast rate items
fast.ci <- quantile(m7_post$DV[m7_post$context=="fast"], probs=c(0.025,0.975))

## create Figure 10
pdf(file="./rtMRI-velum/plots/speech_rate_v_nasal.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m7_plot <- m7_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Effect of speech rate on duration of vowel nasalization") + 
  ylab("Posterior probability (density)") + xlab("Duration (ms)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m7_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m7_plot + 
  # whisker for context 1
  geom_segment(aes(x=norm.ci[1], xend=norm.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[2]) +
  geom_segment(aes(x=norm.ci[1], xend=norm.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=norm.ci[2], xend=norm.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[2]) +
  # whisker for context 2
  geom_segment(aes(x=fast.ci[1], xend=fast.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[1]) +
  geom_segment(aes(x=fast.ci[1], xend=fast.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=fast.ci[2], xend=fast.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[1]) +
  scale_x_continuous(breaks=seq(-1000,1000,5)) + theme_bw()
dev.off()

# sensitivity analysis (Betancourt 2018)

m7_fixed <- fixef(m7) %>% as_tibble(rownames = "term")

m7_fixed %>%
  mutate(
    theta = c(0, 0, 0),
    sigma_prior = c(3, 1, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 100)


# now, we look at the duration of the consonant nasalization
subdat$DV <- subdat$dur.NN

# full model
m8 <- brms::brm(
  DV ~ stress +
    speech_rate_c +
    (1 + stress + speech_rate_c | speaker) +
    (1 + stress + speech_rate_c | word),
  data = subdat,
  family = lognormal(),
  prior = priors,
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/rate_n_nasal",
  save_all_pars = TRUE
)

# model summary
summary(m8)

# posterior predictive check
pp_check(m8, nsamples = 100)

# reduced/null model
m8_null <- brms::brm(
  DV ~
    speech_rate_c +
    (1 + stress + speech_rate_c | speaker) +
    (1 + stress + speech_rate_c | word),
  data = subdat,
  family = lognormal(),
  prior = priors_null,
  seed = my.seed,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/rate_n_nasal_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
#brms::bayes_factor(m8, m8_null)

# calculate the marginal posteriors of the full model
m8_post <- brms::posterior_samples(m8, pars="b_") %>%
  dplyr::mutate(
    fast = exp(b_Intercept),
    normal = exp(b_Intercept) * exp(b_stressN)
  ) %>%
  dplyr::select(normal, fast) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of normal rate items
norm.ci <- quantile(m8_post$DV[m8_post$context=="normal"], probs=c(0.025,0.975))
# 95% CI of fast rate items
fast.ci <- quantile(m8_post$DV[m8_post$context=="fast"], probs=c(0.025,0.975))


## create Figure 11
pdf(file="./rtMRI-velum/plots/speech_rate_n_nasal.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m8_plot <- m8_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Effect of speech rate on duration of consonant nasalization") + 
  ylab("Posterior probability (density)") + xlab("Duration (ms)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m8_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m8_plot + 
  # whisker for context 1
  geom_segment(aes(x=norm.ci[1], xend=norm.ci[2], y=-yrange/10, yend=-yrange/10), col=my.cols[2]) +
  geom_segment(aes(x=norm.ci[1], xend=norm.ci[1], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[2]) +
  geom_segment(aes(x=norm.ci[2], xend=norm.ci[2], y=(-yrange/10+yrange/40), yend=(-yrange/10-yrange/40)), col=my.cols[2]) +
  # whisker for context 2
  geom_segment(aes(x=fast.ci[1], xend=fast.ci[2], y=-yrange/20, yend=-yrange/20),col=my.cols[1]) +
  geom_segment(aes(x=fast.ci[1], xend=fast.ci[1], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[1]) +
  geom_segment(aes(x=fast.ci[2], xend=fast.ci[2], y=(-yrange/20+yrange/40), yend=(-yrange/20-yrange/40)), col=my.cols[1]) +
  scale_x_continuous(breaks=seq(-1000,1000,5)) + theme_bw()
dev.off()

# sensitivity analysis (Betancourt 2018)

m8_fixed <- fixef(m8) %>% as_tibble(rownames = "term")

m8_fixed %>%
  mutate(
    theta = c(0, 0, 0),
    sigma_prior = c(3, 1, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 200)


# # # # # # # # # # # # # # # # # # # # # #
# Figure 11                               #
# Velum gesture plateau acceleration peak #
# # # # # # # # # # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- scale(subdat$accel.peak)


# full model priors:
priors <- c(
  prior(normal(0, 3), class = Intercept),
  prior(normal(0, 3), class = b, coef = voicingvoiceless),
  prior(normal(0, 3), class = b, coef = speech_rate_c),
  prior(cauchy(0, 1), class = sd),
  prior(cauchy(0, 1), class = sigma),
  prior(lkj(2), class = cor)
)

# null model priors:
priors_null <- c(
  prior(normal(0, 3), class = Intercept),
  prior(normal(0, 3), class = b, coef = speech_rate_c),
  prior(cauchy(0, 1), class = sd),
  prior(cauchy(0, 1), class = sigma),
  prior(lkj(2), class = cor)
)


# full model
m10 <- brms::brm(
  DV ~ voicing +
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = gaussian,
  prior = priors,
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
summary(m10)

# posterior predictive check
pp_check(m10, nsamples = 100)

# reduced/null model
m10_null <- brms::brm(
  DV ~
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = gaussian,
  prior = priors_null,
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
#brms::bayes_factor(m10, m10_null)

# calculate the marginal posteriors of the full model
m10_post <- brms::posterior_samples(m10, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(m10_post$DV[m10_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(m10_post$DV[m10_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 9
pdf(file="./rtMRI-velum/plots/plateau_accel_peak.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m10_plot <- m10_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Maximum acceleration within gesture plateau") + 
  ylab("Posterior probability (density)") + xlab("Velum movement acceleration (normalized mm/s^2)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m10_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m10_plot + 
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

# sensitivity analysis (Betancourt 2018)

m10_fixed <- fixef(m10) %>% as_tibble(rownames = "term")

m10_fixed %>%
  mutate(
    theta = c(0, 0, 0),
    sigma_prior = c(3, 3, 3),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 5)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Figure XX                                                   #
# Integral of velum movement in vowel (area under the curve)  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# create the dependent variable
subdat$DV <- scale(subdat$quad.coef)

# full model
m11 <- brms::brm(
  DV ~ voicing +
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = gaussian,
  prior = priors,
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
summary(m11)

# posterior predictive check
pp_check(m11, nsamples = 100)

# reduced/null model
m11_null <- brms::brm(
  DV ~
    speech_rate_c +
    (1 + voicing + speech_rate_c | speaker) +
    (1 + voicing + speech_rate_c | word),
  data = subdat,
  family = gaussian,
  prior = priors_null,
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
#brms::bayes_factor(m11, m11_null)

# calculate the marginal posteriors of the full model
m11_post <- brms::posterior_samples(m11, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)

# calculate the 90% credible intervals
# 95% CI of /nd/ items
nd.ci <- quantile(m11_post$DV[m11_post$context=="nd"], probs=c(0.025,0.975))
# 95% CI of /nt/ items
nt.ci <- quantile(m11_post$DV[m11_post$context=="nt"], probs=c(0.025,0.975))

## create Figure 9
pdf(file="./rtMRI-velum/plots/quad_coef.pdf",width=7,height=3,onefile=T,pointsize=16)
# perceptually distinct colors that are also safe for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")
# base plot
m11_plot <- m11_post %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Quadratic coefficient of second-order polynomial") + 
  ylab("Posterior probability (density)") + xlab("Quadratic coefficient (normalized)")
# get y-range values of the base plot (for determining height of CI whisker bars)
yrange <- ggplot_build(m11_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
# add CI whiskers to base plot
m11_plot + 
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


# sensitivity analysis (Betancourt 2018)

m11_fixed <- fixef(m11) %>% as_tibble(rownames = "term")

m11_fixed %>%
  mutate(
    theta = c(0, 0, 0),
    sigma_prior = c(3, 3, 3),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 5)
