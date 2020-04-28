library(tidyverse)
library(brms)

# Calculations ----

m1 <- readRDS("./rtMRI-velum/models/dur.rds")

m1

m1_fixed <- fixef(m1)

# Duration of velum gesture when C is voiced (Intercept)
# Lower
exp(m1_fixed[1, 3])
# Upper
exp(m1_fixed[1, 4])

# Effect of voicing in odds
# 
# Lower
exp(m1_fixed[2, 3])
# Upper
exp(m1_fixed[2, 4])

# Effect of voicing in percentage decrease
(1 - exp(m1_fixed[2, 4])) * 100
(1 - exp(m1_fixed[2, 3])) * 100

# The effect is a decrease in duration between 6.5% and 14%

m1_post <- brms::posterior_samples(m1, pars = "b_") %>%
  dplyr::mutate(
    nd = exp(b_Intercept),
    nt = exp(b_Intercept + b_voicingvoiceless),
    diff = nt - nd
  )

# Posterior distribution of the effect of voicing in ms
m1_post %>%
  ggplot(aes(diff)) + geom_density()

# 95% CI of effect of voicing in ms. Use this for reporting.
quantile(m1_post$diff, probs = c(0.025, 0.975))

# The effect of voicing is a decrease of duration between 19 and 42 ms.

# Predicted probabilities of mean duration in voiced and voiceless
m1_post %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, DV) %>%
  ggplot(aes(DV, fill = context)) +
  geom_density(alpha = 0.5)

# 95% CI of the probability of mean duration in voiceless
quantile(m1_post$nt, probs = c(0.025, 0.975))

# Reporting ----

# According to the model intercept, the velum gesture duration when C is voiced is between 286 and 318 ms at 95% confidene ($\hat{theta}$ = ..., SE = ... [in log-odds]). When C is voiceless, the duration of the gesture decreases by 6.5-15% at 95% confidence ($\hat{theta}$ = ..., SE = ... [in log-odds]). This change corresponds to a decrese of 19-42 ms. Speech rate has ...

# Somewhere in the paper we will have to be explicit about how we made the calculations above, to aid the reader.

# What really matters is the effect of voicing, not the predicted probability of the mean duration of the gesture when C is voiceless.
