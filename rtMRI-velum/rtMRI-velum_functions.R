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
brms::bayes_factor(dur, dur_null)

# calculate the marginal posteriors of the full model
dur_post <- brms::posterior_samples(dur, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)



# # # # # # # # # # # #
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
brms::bayes_factor(onset, onset_null)

# calculate the marginal posteriors of the full model
onset_post <- brms::posterior_samples(onset, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)



# # # # # # # # # # # # # # # #
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
brms::bayes_factor(gest.max, gest.max_null)

# calculate the marginal posteriors of the full model
gest.max_post <- brms::posterior_samples(gest.max, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)



# # # # # # # # # # # # #
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
brms::bayes_factor(offset, offset_null)

# calculate the marginal posteriors of the full model
offset_post <- brms::posterior_samples(offset, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)



# # # # # # # # # # # # # # # # #
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
brms::bayes_factor(gest.max.mag, gest.max.mag_null)

# calculate the marginal posteriors of the full model
gest.max.mag_post <- brms::posterior_samples(gest.max.mag, pars="b_") %>%
  dplyr::mutate(
    nd = plogis(b_Intercept),
    nt = plogis(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)



# # # # # # # # # # # # # #
# Gesture onset stiffness #
# # # # # # # # # # # # # #

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
brms::bayes_factor(stiff.ons, stiff.ons_null)

# calculate the marginal posteriors of the full model
stiff.ons_post <- brms::posterior_samples(stiff.ons, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)



# # # # # # # # # # # # # # #
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
brms::bayes_factor(stiff.off, stiff.off_null)

# calculate the marginal posteriors of the full model
stiff.off_post <- brms::posterior_samples(stiff.off, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)



# # # # # # #
# Kurtosis  #
# # # # # # #

# create the dependent variable
subdat$DV <- subdat$kurtosis

# full model
kurt <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(normal(0, 3), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 1), class = sd),
            prior(cauchy(0, 1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/kurtosis",
  save_all_pars = TRUE
)

# model summary
summary(kurt)

# reduced/null model
kurt_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(cauchy(0, 1), class = sd),
            prior(cauchy(0, 1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/kurtosis_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(kurt, kurt_null)

# calculate the marginal posteriors of the full model
kurt_post <- brms::posterior_samples(kurt, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)



# # # # # # # # #
# Crest factor  #
# # # # # # # # #

# create the dependent variable
subdat$DV <- subdat$crest.fact

# full model
crest.fact <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(normal(0, 3), class = b, coef = voicingvoiceless),
            prior(cauchy(0, 1), class = sd),
            prior(cauchy(0, 1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/crest_factor",
  save_all_pars = TRUE
)

# model summary
summary(crest.fact)

# reduced/null model
crest.fact_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 3), class = Intercept),
            prior(cauchy(0, 1), class = sd),
            prior(cauchy(0, 1), class = sigma),
            prior(lkj(2), class = cor)),
  seed = my.seed,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = core.num,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  file = "./rtMRI-velum/models/crest_factor_null",
  save_all_pars = TRUE
)

# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(crest.fact, crest.fact_null)

# calculate the marginal posteriors of the full model
crest.fact_post <- brms::posterior_samples(crest.fact, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
brms::bayes_factor(integ, integ_null)

# calculate the marginal posteriors of the full model
integ_post <- brms::posterior_samples(integ, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)



# # # # # # # # # # # #
# Plotting posteriors #
# # # # # # # # # # # #


# colors to be used for plotting (suitable for B/W printing
my.cols <- c("#2c7fb8","#7fcdbb")


## Figure 5: duration, onset, peak (timing), offset ##

# gather posteriors
fig5_post <- bind_rows(
  dur_post %>% rename(`duration` = DV) %>% pivot_longer(`duration`, names_to = "outcome", values_to = "estimate"),
  onset_post %>% rename(`onset` = DV) %>% pivot_longer(`onset`, names_to = "outcome", values_to = "estimate"),
  gest.max_post %>% rename(`peak` = DV) %>% pivot_longer(`peak`, names_to = "outcome", values_to = "estimate"),
  offset_post %>% rename(`offset` = DV) %>% pivot_longer(`offset`, names_to = "outcome", values_to = "estimate")
) %>%
  mutate(outcome = factor(outcome, levels = c("offset","peak","onset","duration")))

# flip sign of velum gesture onset posteriors for interpretability
fig5_post$estimate[fig5_post$outcome=="onset"] <- -fig5_post$estimate[fig5_post$outcome=="onset"]

# make separate plots
p1 <- fig5_post %>%
  filter(outcome == "duration") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,500,5)) +
  coord_cartesian(ylim = c(1.3, 1.5), xlim = c(250,327)) +
  scale_fill_manual(values = my.cols) +
  labs(x = "Duration (ms)", y = element_blank()) + theme_bw() + theme(legend.position = "none",axis.text=element_text(size=12),axis.title=element_text(size=13))

p2 <- fig5_post %>%
  filter(outcome %in% c("onset","peak","offset")) %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  geom_vline(xintercept = 0, lty=2) +
  scale_x_continuous(breaks=seq(-200,300,20)) +
  coord_cartesian(xlim = c(-125,205)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Time (ms) relative to acoustic vowel offset", y = element_blank(), fill = "Context") + theme_bw() + theme(axis.text=element_text(size=12),axis.title=element_text(size=13),panel.grid.major.y=element_blank())

# save coposite plot
pdf(file="./rtMRI-velum/plots/time_plots.pdf",width=9.5,height=5,onefile=T,pointsize=14)
(p1 + p2) + patchwork::plot_layout(ncol = 1, guides = "collect") + theme(legend.position = "right")
dev.off()


## Figure 6: peak (magnitude) ##

# gather posteriors
fig6_post <- bind_rows(
  gest.max.mag_post %>% rename(`peak    \nmagnitude` = DV) %>% pivot_longer(`peak    \nmagnitude`, names_to = "outcome", values_to = "estimate")
)

# make plot
p1 <- fig6_post %>%
  filter(outcome == "peak    \nmagnitude") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,1,0.025)) +
  coord_cartesian(xlim = c(0.58,0.775), ylim = c(1.4,1.4)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Velum opening magnitude (speaker-scaled)", y = element_blank(), fill = "Context") + theme_bw() + theme(axis.text=element_text(size=12),axis.title=element_text(size=13))

# save plot
pdf(file="./rtMRI-velum/plots/magnitude_plot.pdf",width=9,height=3,onefile=T,pointsize=14)
p1
dev.off()


## Figure 7: opening stiffness, closing stiffness ##

# gather posteriors
fig7_post <- bind_rows(
  stiff.ons_post %>% rename(`opening\nstiffness` = DV) %>% pivot_longer(`opening\nstiffness`, names_to = "outcome", values_to = "estimate"),
  stiff.off_post %>% rename(`closing \nstiffness` = DV) %>% pivot_longer(`closing \nstiffness`, names_to = "outcome", values_to = "estimate")
)

# make composite plot on same scale
p1 <- fig7_post %>%
  filter(outcome %in% c("opening\nstiffness","closing \nstiffness")) %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,30,0.5)) +
  coord_cartesian(xlim = c(12,17.5), ylim = c(1.4,2.3)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Stiffness", y = element_blank(), fill = "Context") + theme_bw() + theme(axis.text=element_text(size=12),axis.title=element_text(size=13),panel.grid.major.y=element_blank())

# save plot
pdf(file="./rtMRI-velum/plots/stiffness_plots.pdf",width=9,height=4,onefile=T,pointsize=14)
p1
dev.off()


## Figure 9: kurtosis, crest factor ##

# gather posteriors
fig9_post <- bind_rows(
  kurt_post %>% rename(`kurtosis` = DV) %>% pivot_longer(`kurtosis`, names_to = "outcome", values_to = "estimate"),
  crest.fact_post %>% rename(`crest \nfactor` = DV) %>% pivot_longer(`crest \nfactor`, names_to = "outcome", values_to = "estimate")
)

# make separate plots
p1 <- fig9_post %>%
  filter(outcome == "kurtosis") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,5,0.05)) +
  coord_cartesian(xlim = c(2.81,3.18), ylim = c(1.4,1.4)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Kurtosis", y = element_blank(), fill = "Context") + theme_bw() + theme(legend.position = "none",axis.text=element_text(size=12),axis.title=element_text(size=13),panel.grid.major.y=element_blank())

p2 <- fig9_post %>%
  filter(outcome == "crest \nfactor") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,3,0.01)) +
  coord_cartesian(xlim = c(1.805,1.885), ylim = c(1.4,1.4)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Crest factor (ratio of peak to average displacement)", y = element_blank(), fill = "Context") + theme_bw() + theme(legend.position = "none",axis.text=element_text(size=12),axis.title=element_text(size=13),panel.grid.major.y=element_blank())

# save composite plot
pdf(file="./rtMRI-velum/plots/peakedness_plots.pdf",width=9,height=5,onefile=T,pointsize=14)
(p1 + p2) + patchwork::plot_layout(ncol = 1, guides = "collect") + theme(legend.position = "right")
dev.off()


## Figure 10: velum displacement integral

# gather posteriors
fig10_post <- bind_rows(
  integ_post %>% rename(`gesture\nintegral` = DV) %>% pivot_longer(`gesture\nintegral`, names_to = "outcome", values_to = "estimate")
)

# make plot
p1 <- fig10_post %>%
  filter(outcome == "gesture\nintegral") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(10,50,1)) +
  coord_cartesian(xlim = c(18,32), ylim = c(1.4,1.4)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Velum displacement integral (time X magnitude)", y = element_blank(), fill = "Context") + theme_bw() + theme(axis.text=element_text(size=12),axis.title=element_text(size=13))

# save plot
pdf(file="./rtMRI-velum/plots/integral_plot.pdf",width=9,height=3,onefile=T,pointsize=14)
p1
dev.off()
