## ----setup, include=TRUE, message=FALSE--------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = here::here())
library(tidyverse)
library(patchwork)
library(brms)
library(tidybayes)
library(extraDistr)
library(HDInterval)

core.num <- parallel::detectCores()
options(mc.cores = core.num)

my.seed <- 123
set.seed(my.seed)


## ----read-data---------------------------------------------------------------------------------------
matdat <- read.csv("./rtMRI-velum/velum_data.csv", header = T)

# get rid of items that gesture onsets that begin *after* the vowel offset (only 5/7152 total items)
matdat <- matdat[(matdat$velumopening_gesture_on - matdat$Vokal_off) < 0,]

# separate coda contexts by alveolar voiced stop vs. alveolar voiceless stop
voiceless <- c("nt__", "nt_@", "nt_6", "nt_a")
voiced    <- c("nd_@", "nd_6", "nd_a")

matdat$voicing <- c()
matdat$voicing[matdat$post %in% voiceless]  <- "voiceless"
matdat$voicing[matdat$post %in% voiced]     <- "voiced"

# only include alveolar nasal items preceding a voiced or voiceless stop consonant
coda <- c(voiceless, voiced)

# only include neutrally stressed utterances
stresses <- c("N")

# subset the data
subdat <- matdat[matdat$post %in% coda & matdat$stress %in% stresses, ]


## ----half-cauchy-------------------------------------------------------------------------------------
inverseCDF(c(0.025, 0.975), phcauchy, sigma = 0.1)


## ----dur, cache=TRUE---------------------------------------------------------------------------------
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


## ----dur-summary-------------------------------------------------------------------------------------
dur


## ----dur-null, cache=TRUE----------------------------------------------------------------------------
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


## ----dur-bf, cache=TRUE, dependson=c("dur", "dur-null")----------------------------------------------
# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(dur, dur_null)


## ----dur-post----------------------------------------------------------------------------------------
# calculate the marginal posteriors of the full model
dur_post <- brms::posterior_samples(dur, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)


## ----dur-pp------------------------------------------------------------------------------------------
pp_check(dur, nsamples = 50) + theme_minimal()


## ----dur-sens----------------------------------------------------------------------------------------
dur_fixed <- fixef(dur) %>% as_tibble(rownames = "term")

dur_fixed %>%
  mutate(
    theta = c(0, 0),
    sigma_prior = c(3, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250) +
  theme_minimal()


## ----onset, cache=TRUE-------------------------------------------------------------------------------
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


## ----onset-summary-----------------------------------------------------------------------------------
onset


## ----onset-null, cache=TRUE--------------------------------------------------------------------------
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


## ----onset-bf, cache=TRUE, dependson=c("onset", "onset-null")----------------------------------------
# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(onset, onset_null)


## ----onset-post--------------------------------------------------------------------------------------
# calculate the marginal posteriors of the full model
onset_post <- brms::posterior_samples(onset, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)


## ----onset-pp----------------------------------------------------------------------------------------
pp_check(onset, nsamples = 50) + theme_minimal()


## ----onset-sens--------------------------------------------------------------------------------------
onset_fixed <- fixef(onset) %>% as_tibble(rownames = "term")

onset_fixed %>%
  mutate(
    theta = c(0, 0),
    sigma_prior = c(3, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250) +
  theme_minimal()


## ----half-cauchy-get-max-----------------------------------------------------------------------------
inverseCDF(c(0.025, 0.975), phcauchy, sigma = 5)


## ----gest-max, cache=TRUE----------------------------------------------------------------------------
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


## ----gest-max-summary--------------------------------------------------------------------------------
gest.max


## ----gest-max-null, cache=TRUE-----------------------------------------------------------------------
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


## ----gest-max-bf, cache=TRUE, dependson=c("gest-max", "gest-max-null")-------------------------------
# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(gest.max, gest.max_null)


## ----gest-max-post-----------------------------------------------------------------------------------
# calculate the marginal posteriors of the full model
gest.max_post <- brms::posterior_samples(gest.max, pars="b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%    
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)


## ----gest-max-pp-------------------------------------------------------------------------------------
pp_check(gest.max, nsamples = 50) + theme_minimal()


## ----gest-max-sens-----------------------------------------------------------------------------------
gest.max_fixed <- fixef(gest.max) %>% as_tibble(rownames = "term")

gest.max_fixed %>%
  mutate(
    theta = c(0, 0),
    sigma_prior = c(200, 100),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250) +
  theme_minimal()


## ----offset, cache=TRUE------------------------------------------------------------------------------
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


## ----offset-summary----------------------------------------------------------------------------------
offset


## ----offset-null, cache=TRUE-------------------------------------------------------------------------
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


## ----offset-bf, cache=TRUE, dependson=c("offset", "offset-null")-------------------------------------
# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(offset, offset_null)


## ----offset-post-------------------------------------------------------------------------------------
# calculate the marginal posteriors of the full model
offset_post <- brms::posterior_samples(offset, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)


## ----offset-pp---------------------------------------------------------------------------------------
pp_check(offset, nsamples = 50) + theme_minimal()


## ----offset-sens-------------------------------------------------------------------------------------
offset_fixed <- fixef(offset) %>% as_tibble(rownames = "term")

offset_fixed %>%
  mutate(
    theta = c(0, 0),
    sigma_prior = c(3, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250) +
  theme_minimal()


## ----gest-max-mag, cache=TRUE------------------------------------------------------------------------
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


## ----gest-max-mag-summary----------------------------------------------------------------------------
gest.max.mag


## ----gest-max-mag-null, cache=TRUE-------------------------------------------------------------------
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


## ----gest-max-mag-bf, cache=TRUE, dependson=c("gest-max-mag", "gest-max-mag-null")-------------------
# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(gest.max.mag, gest.max.mag_null)


## ----gest-max-mag-post-------------------------------------------------------------------------------
# calculate the marginal posteriors of the full model
gest.max.mag_post <- brms::posterior_samples(gest.max.mag, pars="b_") %>%
  dplyr::mutate(
    nd = plogis(b_Intercept),
    nt = plogis(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)


## ----gest-max-mag-pp---------------------------------------------------------------------------------
pp_check(gest.max.mag, nsamples = 50) + theme_minimal()


## ----gest-max-mag-sens-------------------------------------------------------------------------------
gest.max.mag_fixed <- fixef(gest.max.mag) %>% as_tibble(rownames = "term")

gest.max.mag_fixed %>%
  mutate(
    theta = c(0, 0),
    sigma_prior = c(10, 5),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250) +
  theme_minimal()


## ----stiff-ons, cache=TRUE---------------------------------------------------------------------------
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


## ----stiff-ons-summary-------------------------------------------------------------------------------
stiff.ons


## ----stiff-ons-null, cache=TRUE----------------------------------------------------------------------
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


## ----stiff-ons-bf, cache=TRUE, dependson=c("stiff-ons", "stiff-ons-null")----------------------------
# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(stiff.ons, stiff.ons_null)


## ----stiff-ons-post----------------------------------------------------------------------------------
# calculate the marginal posteriors of the full model
stiff.ons_post <- brms::posterior_samples(stiff.ons, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)


## ----stiff-ons-pp------------------------------------------------------------------------------------
pp_check(stiff.ons, nsamples = 50) + theme_minimal()


## ----stiff-ons-sens----------------------------------------------------------------------------------
stiff.ons_fixed <- fixef(stiff.ons) %>% as_tibble(rownames = "term")

stiff.ons_fixed %>%
  mutate(
    theta = c(0, 0),
    sigma_prior = c(3, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250) +
  theme_minimal()


## ----stiff-off, cache=TRUE---------------------------------------------------------------------------
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


## ----stiff-off-summary-------------------------------------------------------------------------------
stiff.off


## ----stiff-off-null, cache=TRUE----------------------------------------------------------------------
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


## ----stiff-off-bf, cache=TRUE, dependson=c("stiff-off", "stiff-off-null")----------------------------
# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(stiff.off, stiff.off_null)


## ----stiff-off-post----------------------------------------------------------------------------------
# calculate the marginal posteriors of the full model
stiff.off_post <- brms::posterior_samples(stiff.off, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept) * exp(b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)


## ----stiff-off-pp------------------------------------------------------------------------------------
pp_check(stiff.off, nsamples = 50) + theme_minimal()


## ----stiff-off-sens----------------------------------------------------------------------------------
stiff.off_fixed <- fixef(stiff.off) %>% as_tibble(rownames = "term")

stiff.off_fixed %>%
  mutate(
    theta = c(0, 0),
    sigma_prior = c(3, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250) +
  theme_minimal()


## ----kurt, cache=TRUE--------------------------------------------------------------------------------
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


## ----kurt-summary------------------------------------------------------------------------------------
kurt


## ----kurt-null, cache=TRUE---------------------------------------------------------------------------
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


## ----kurt-bf, cache=TRUE, dependson=c("kurt", "kurt-null")-------------------------------------------
# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(kurt, kurt_null)


## ----kurt-post---------------------------------------------------------------------------------------
# calculate the marginal posteriors of the full model
kurt_post <- brms::posterior_samples(kurt, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)


## ----kurt-pp-----------------------------------------------------------------------------------------
pp_check(kurt, nsamples = 50) + theme_minimal()


## ----kurt-sens---------------------------------------------------------------------------------------
kurt_fixed <- fixef(kurt) %>% as_tibble(rownames = "term")

kurt_fixed %>%
  mutate(
    theta = c(0, 0),
    sigma_prior = c(3, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250) +
  theme_minimal()


## ----crest-fact, cache=TRUE--------------------------------------------------------------------------
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


## ----crest-fact-summary------------------------------------------------------------------------------
crest.fact


## ----crest-fact-null, cache=TRUE---------------------------------------------------------------------
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


## ----creast-fact-bf, cache=TRUE, dependson=c("crest-fact", "crest-fact-null")------------------------
# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(crest.fact, crest.fact_null)


## ----crest-fact-post---------------------------------------------------------------------------------
# calculate the marginal posteriors of the full model
crest.fact_post <- brms::posterior_samples(crest.fact, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)


## ----crest-fact-pp-----------------------------------------------------------------------------------
pp_check(crest.fact, nsamples = 50) + theme_minimal()


## ----crest-fact-sens---------------------------------------------------------------------------------
crest.fact_fixed <- fixef(crest.fact) %>% as_tibble(rownames = "term")

crest.fact_fixed %>%
  mutate(
    theta = c(0, 0),
    sigma_prior = c(3, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250) +
  theme_minimal()


## ----integ, cache=TRUE-------------------------------------------------------------------------------
# create the dependent variable
subdat$DV <- subdat$vowel.integ

# full model
integ <- brms::brm(
  DV ~ voicing +
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 5), class = Intercept),
            prior(normal(0, 5), class = b, coef = voicingvoiceless),
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


## ----integ-summary-----------------------------------------------------------------------------------
integ


## ----integ-null, cache=TRUE--------------------------------------------------------------------------
# reduced/null model
integ_null <- brms::brm(
  DV ~
    (1 + voicing | speaker) +
    (1 + voicing | word),
  data = subdat,
  family = lognormal(),
  prior = c(prior(normal(0, 5), class = Intercept),
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


## ----integ-bf, cache=TRUE, dependson=c("integ", "integ-null")----------------------------------------
# calculate the Bayes factor of the difference between the full and null models
brms::bayes_factor(integ, integ_null)


## ----integ-post--------------------------------------------------------------------------------------
# calculate the marginal posteriors of the full model
integ_post <- brms::posterior_samples(integ, pars="b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept + b_voicingvoiceless)
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV)


## ----integ-pp----------------------------------------------------------------------------------------
pp_check(integ, nsamples = 50) + theme_minimal()


## ----integ-sens--------------------------------------------------------------------------------------
integ_fixed <- fixef(integ) %>% as_tibble(rownames = "term")

integ_fixed %>%
  mutate(
    theta = c(0, 0),
    sigma_prior = c(3, 1),
    # it's called here std.error but is the standard deviation
    z = abs((Estimate - theta) / Est.Error),
    s = 1 - (Est.Error^2 / sigma_prior^2)
  ) %>%
  ggplot(aes(s, z, label = term)) +
  geom_point() +
  geom_label(nudge_x = -0.1, size = 2, alpha = 0.5) +
  xlim(0, 1) + ylim(0, 250) +
  theme_minimal()


## ----cols--------------------------------------------------------------------------------------------
# colors to be used for plotting (suitable for B/W printing)
my.cols <- c("#2c7fb8","#7fcdbb")


## ----figure-5----------------------------------------------------------------------------------------
# gather posteriors
fig5_post <- bind_rows(
  dur_post %>%
    rename(`duration` = DV)
  %>% pivot_longer(`duration`, names_to = "outcome", values_to = "estimate"),
  onset_post %>%
    rename(`onset` = DV) %>%
    pivot_longer(`onset`, names_to = "outcome", values_to = "estimate"),
  gest.max_post %>%
    rename(`peak` = DV) %>%
    pivot_longer(`peak`, names_to = "outcome", values_to = "estimate"),
  offset_post %>%
    rename(`offset` = DV) %>%
    pivot_longer(`offset`, names_to = "outcome", values_to = "estimate")
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
  labs(x = "Duration (ms)", y = element_blank()) + theme_bw() +
  theme(legend.position = "none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=13))

p2 <- fig5_post %>%
  filter(outcome %in% c("onset","peak","offset")) %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  geom_vline(xintercept = 0, lty=2) +
  scale_x_continuous(breaks=seq(-200,300,20)) +
  coord_cartesian(xlim = c(-125,205)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Time (ms) relative to acoustic vowel offset",
       y = element_blank(), fill = "Context") + theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=13),
        panel.grid.major.y=element_blank())


## ----figure-5-save-----------------------------------------------------------------------------------
# save coposite plot
pdf(file="./rtMRI-velum/plots/time_plots.pdf",width=9.5,height=5,
    onefile=T,pointsize=14)
(p1 + p2) + patchwork::plot_layout(ncol = 1, guides = "collect") +
  theme(legend.position = "right")
dev.off()


## ----figure-6----------------------------------------------------------------------------------------
# gather posteriors
fig6_post <- bind_rows(
  gest.max.mag_post %>%
    rename(`peak    \nmagnitude` = DV) %>%
    pivot_longer(`peak    \nmagnitude`, names_to = "outcome",
                 values_to = "estimate")
) 

# make plot
p1 <- fig6_post %>%
  filter(outcome == "peak    \nmagnitude") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,1,0.025)) +
  coord_cartesian(xlim = c(0.58,0.775), ylim = c(1.4,1.4)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Velum opening magnitude (speaker-scaled)", y = element_blank(),
       fill = "Context") + theme_bw() +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=13))


## ----figure-6-save-----------------------------------------------------------------------------------
# save plot
pdf(file="./rtMRI-velum/plots/magnitude_plot.pdf",width=9,height=3,onefile=T,
    pointsize=14)
p1
dev.off()


## ----figure-7----------------------------------------------------------------------------------------
# gather posteriors
fig7_post <- bind_rows(
  stiff.ons_post %>% rename(`opening\nstiffness` = DV) %>%
    pivot_longer(`opening\nstiffness`, names_to = "outcome", values_to = "estimate"),
  stiff.off_post %>% rename(`closing \nstiffness` = DV) %>%
    pivot_longer(`closing \nstiffness`, names_to = "outcome", values_to = "estimate")
) 

# make composite plot on same scale
p1 <- fig7_post %>%
  filter(outcome %in% c("opening\nstiffness","closing \nstiffness")) %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,30,0.5)) +
  coord_cartesian(xlim = c(12,17.5), ylim = c(1.4,2.3)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Stiffness", y = element_blank(), fill = "Context") + theme_bw() +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=13),
        panel.grid.major.y=element_blank())


## ----figure-7-save-----------------------------------------------------------------------------------
# save plot
pdf(file="./rtMRI-velum/plots/stiffness_plots.pdf",width=9,height=4,onefile=T,
    pointsize=14)
p1
dev.off()


## ----figure-9----------------------------------------------------------------------------------------
# gather posteriors
fig9_post <- bind_rows(
  kurt_post %>% rename(`kurtosis` = DV) %>%
    pivot_longer(`kurtosis`, names_to = "outcome", values_to = "estimate"),
  crest.fact_post %>% rename(`crest \nfactor` = DV) %>%
    pivot_longer(`crest \nfactor`, names_to = "outcome", values_to = "estimate")
) 

# make separate plots
p1 <- fig9_post %>%
  filter(outcome == "kurtosis") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,5,0.05)) +
  coord_cartesian(xlim = c(2.81,3.18), ylim = c(1.4,1.4)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Kurtosis", y = element_blank(), fill = "Context") +
  theme_bw() +
  theme(legend.position = "none",axis.text=element_text(size=12),
        axis.title=element_text(size=13),panel.grid.major.y=element_blank())

p2 <- fig9_post %>%
  filter(outcome == "crest \nfactor") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,3,0.01)) +
  coord_cartesian(xlim = c(1.805,1.885), ylim = c(1.4,1.4)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Crest factor (ratio of peak to average displacement)",
       y = element_blank(), fill = "Context") + theme_bw() +
  theme(legend.position = "none",axis.text=element_text(size=12),
        axis.title=element_text(size=13),panel.grid.major.y=element_blank())


## ----figure-9-save-----------------------------------------------------------------------------------
# save composite plot
pdf(file="./rtMRI-velum/plots/peakedness_plots.pdf",width=9,height=5,
    onefile=T,pointsize=14)
(p1 + p2) + patchwork::plot_layout(ncol = 1, guides = "collect") +
  theme(legend.position = "right")
dev.off()


## ----figure-10---------------------------------------------------------------------------------------
# gather posteriors
fig10_post <- bind_rows(
  integ_post %>% rename(`gesture\nintegral` = DV) %>%
    pivot_longer(`gesture\nintegral`, names_to = "outcome",
                 values_to = "estimate")
) 

# make plot
p1 <- fig10_post %>%
  filter(outcome == "gesture\nintegral") %>%
  ggplot(aes(estimate, outcome, group = context, fill = context)) +
  geom_halfeyeh(slab_color="black", slab_alpha=0.5) +
  scale_x_continuous(breaks=seq(0,100,2)) +
  coord_cartesian(xlim = c(31,64.5), ylim = c(1.4,1.4)) +
  scale_fill_manual(values=my.cols) +
  labs(x = "Velum displacement integral: time (ms) X magnitude (normalized)",
       y = element_blank(), fill = "Context") + theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=13))


## ----figure-10-save----------------------------------------------------------------------------------
# save plot
pdf(file="./rtMRI-velum/plots/integral_plot.pdf",width=9,height=3,onefile=T,
    pointsize=14)
p1
dev.off()

