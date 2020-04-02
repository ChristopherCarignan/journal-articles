require(brms)
require(tidyverse)
require(HDInterval)
require(faintr) # devtools::install_github('michael-franke/bayes_mixed_regression_tutorial/faintr', build_vignettes = TRUE)
require(tidybayes)
require(parallel)
require(viridisLite)
require(ggpmisc)
require(lme4)
require(lmerTest)
require(MuMIn)

setwd("/home/chris/Documents/research/Language_Nasality/")

#### Prepare data ####

options(mc.cores = parallel::detectCores())
my.seed <- 123

set.seed(my.seed)

matdat <- read.csv('velum_data.csv',header = T)
matdat$word <- paste0(matdat$prev,matdat$vowel)


vowels  <- c('E_', 'i:', 'o:', '2:', 'aU', 'U_', 'a:', 'a_', 'y:', 'y_', 'I_', 'E:', 'O_', 'Y_', 'u:', 'aI', 'e:', '9_')
mono    <- c('E_', 'i:', 'o:', '2:', 'U_', 'a:', 'a_', 'y:', 'y_', 'I_', 'E:', 'O_', 'Y_', 'u:', 'e:', '9_')

high  <- c('i:', 'o:', 'y:', 'y_', 'I_', 'Y_', 'u:', 'e:')
low   <- c('a:', 'a_')

voiceless <- c('nt__','nt_@','nt_6','nt_a')
voiced    <- c('nd_@','nd_6','nd_a')

coda <- c(voiceless,voiced)
stresses <- c("N")

#subdat <- matdat[matdat$post %in% coda & matdat$vowel %in% low,]
subdat <- matdat[matdat$post %in% coda & matdat$stress %in% stresses, ]
subdat$voicing <- c()
subdat$voicing[subdat$post %in% voiceless] <- "voiceless"
subdat$voicing[subdat$post %in% voiced] <- "voiced"
subdat$height <- c()
subdat$height[subdat$vowel %in% high] <- "high"
subdat$height[subdat$vowel %in% low] <- "low"

table(subdat$voicing)



# Variable names:
# Vokal_off: time point (s) of the acoustic vowel offset
# velumopening_gesture_dur: duration (s) of the velum gesture
# velumopening_gesture_on: time point (s) of the onset of the velum gesture
# velumopening_gesture_off: time point (s) of the offset of the velum gesture
# velumopening_maxcon_on: time point (s) of the peak of the velum gesture
# velum2US_velumopening_maxcon_onset: magnitude of the peak of the velum gesture


#subdat$DV <- subdat$velumopening_gesture_dur*1000
#subdat$DV <- subdat$velumopening_maxvel_dur*1000
subdat$DV <- (subdat$velumopening_maxvel_on - subdat$Vokal_off)*1000
#subdat$DV <- subdat$velum2US_velumopening_maxcon_onset




#### Model m1 ####

# Use get_prior() to print the default priors of the model
# get_prior(DV ~ voicing + (1+voicing|speaker) + (1|word), data = subdat)

priors <- c(
  prior(normal(0, 100), class = Intercept),
  prior(normal(0, 100), class = b, coef = voicingvoiceless),
  prior(cauchy(0, 10), class = sd),
  prior(cauchy(0, 10), class = sigma),
  prior(lkj(2), class = cor)
)

m1 <- brm(
  formula = DV ~ voicing + (1+voicing|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors,
  file = "./cache/max_mag_N",
  save_all_pars = TRUE
)

m1

# Plot posteriors
plot(m1)

# Posterior predictive checks
pp_check(m1, nsamples = 100)

# Residuals
subdat %>%
  add_predicted_draws(m1) %>%
  summarise(
    p_residual = mean(.prediction < DV),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()


# Null model (i.e., main effect not included)
priors_null <- c(
  prior(normal(0, 100), class = Intercept),
  prior(cauchy(0, 10), class = sd),
  prior(cauchy(0, 10), class = sigma)
)

m1_null <- brm(
  formula = DV ~ (1|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors_null,
  file = "./cache/onset_vel_null_N",
  save_all_pars = TRUE
)

bayes_factor(m1, m1_null)



#### Plot marginal posteriors ####

m1_post <- brms::posterior_samples(m1, pars = "b_") %>%
  dplyr::mutate(
    nd = b_Intercept,
    nt = b_Intercept + b_voicingvoiceless
  ) %>%
  dplyr::select(nd, nt) %>%
  tidyr::gather(context, velum_on)

# # If you add vowel height in the model, assuming high is the reference level
# m1_post <- posterior_samples(m1, pars = "b_") %>%
#   mutate(
#     nd_high = b_Intercept,
#     nd_low = b_Intercept + b_heightlow,
#     nt_high = b_Intercept + b_voicingvoiceless,
#     nt_low = b_Intercept + b_heightlow + b_voicingvoiceless,
#   ) %>%
#   select(nd_high:nt_low) %>%
#   gather(context, velum_on) %>%
#   separate(context, c("voicing", "height"))

# You can get the 90% CI with quantile()
# 90% CI of "nd"
nd.ci <- quantile(m1_post[m1_post$context == "nd",]$velum_on, probs = c(0.05, 0.95))
# 90% CI of "nt"
nt.ci <- quantile(m1_post[m1_post$context == "nt",]$velum_on, probs = c(0.05, 0.95))

# Plot marginal posterior distributions
pdf(file="plots/gesture_max_mag_N.pdf",width=7,height=3,onefile=T,pointsize=16)
#my.cols <- c("red","blue")
#my.cols <- viridisLite::magma(n = 3)
my.cols <- c("#2c7fb8","#7fcdbb")
m1_plot <- m1_post %>%
  ggplot(aes(velum_on, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Velum gesture peak (magnitude)") + 
  ylab("Posterior probability (density)") + xlab("Velum-opening magnitude (speaker-scaled)")
yrange <- ggplot_build(m1_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
m1_plot + 
  geom_segment(aes(x=nd.ci[1],xend=nd.ci[2],y=-yrange/10,yend=-yrange/10),col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1],xend=nd.ci[1],y=(-yrange/10 + yrange/40),yend=(-yrange/10 - yrange/40)),col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2],xend=nd.ci[2],y=(-yrange/10 + yrange/40),yend=(-yrange/10 - yrange/40)),col=my.cols[1]) +
  geom_segment(aes(x=nt.ci[1],xend=nt.ci[2],y=-yrange/20,yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1],xend=nt.ci[1],y=(-yrange/20 + yrange/40),yend=(-yrange/20 - yrange/40)),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2],xend=nt.ci[2],y=(-yrange/20 + yrange/40),yend=(-yrange/20 - yrange/40)),col=my.cols[2]) +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,1,0.05))
dev.off()




## Relation between N duration and duration of vowel nasalization
subdat$dur.N  <- 1000*(subdat$velumopening_maxvel_off-subdat$Vokal_off)
subdat$dur.VN <- 1000*(subdat$Vokal_off-subdat$velumopening_maxvel_on)

cor.test(subdat$dur.VN, subdat$dur.N) 

nt.dat <- subdat[subdat$voicing=="voiceless",]
cor.test(nt.dat$dur.VN, nt.dat$dur.N) 

nd.dat <- subdat[subdat$voicing=="voiced",]
cor.test(nd.dat$dur.VN, nd.dat$dur.N) 


# Estimating R^2 (Beddor, 2009, p.791)
mod.all <- lmer(dur.VN ~ dur.N + (1|speaker), data=subdat)
mod.all.red <- lmer(dur.VN ~ (1|speaker), data=subdat)
anova(mod.all,mod.all.red)
r2.all  <- r.squaredGLMM(mod.all)[1]
int.all <- as.numeric(coef(lm(dur.VN ~ dur.N, data=subdat))[1])
all.xlims <- range(subdat$dur.N)
all.ylims <- all.xlims*-sqrt(r2.all) + int.all

mod.nt <- lmer(dur.VN ~ dur.N + (1|speaker), data=nt.dat)
mod.nt.red <- lmer(dur.VN ~ (1|speaker), data=nt.dat)
anova(mod.nt,mod.nt.red)
r2.nt  <- r.squaredGLMM(mod.nt)[1]
int.nt <- as.numeric(coef(lm(dur.VN ~ dur.N, data=nt.dat))[1])
nt.xlims <- range(nt.dat$dur.N)
nt.ylims <- nt.xlims*-sqrt(r2.nt) + int.nt

mod.nd <- lmer(dur.VN ~ dur.N + (1|speaker), data=nd.dat)
mod.nd.red <- lmer(dur.VN ~ (1|speaker), data=nd.dat)
anova(mod.nd,mod.nd.red)
r2.nd  <- r.squaredGLMM(mod.nd)[1]
int.nd <- as.numeric(coef(lm(dur.VN ~ dur.N, data=nd.dat))[1])
nd.xlims <- range(nd.dat$dur.N)
nd.ylims <- nd.xlims*-sqrt(r2.nd) + int.nd


pdf(file="plots/trading_relation_N.pdf",width=5,height=3,onefile=T,pointsize=16)
ggplot(subdat,aes(x=dur.N, y=dur.VN, col=voicing, shape=voicing)) + 
  geom_point(alpha=0.7,cex=2) + scale_shape_manual(values=c(2,4)) +
  #geom_abline(aes(intercept=int.nt,slope=-sqrt(r2.nt)),lty=2,lwd=0.8) +
  #geom_abline(aes(intercept=int.nd,slope=-sqrt(r2.nd)),lty=3,lwd=0.8) +
  #geom_segment(x=all.xlims[1],xend=all.xlims[2],y=all.ylims[1],yend=all.ylims[2],lty=1,lwd=0.6,col='black') +
  geom_segment(x=nt.xlims[1],xend=nt.xlims[2],y=nt.ylims[1],yend=nt.ylims[2],lty=2,lwd=0.6,col='black') +
  geom_segment(x=nd.xlims[1],xend=nd.xlims[2],y=nd.ylims[1],yend=nd.ylims[2],lty=3,lwd=0.6,col='black') +
  #geom_smooth(method='lm',se=F,aes(lty=voicing),col='black') + 
  #scale_linetype_manual(values=c(2,3)) +
  #geom_smooth(method='lm',se=F,col='black',lty=1) +
  scale_x_continuous(breaks=seq(-100,300,50)) +
  scale_y_continuous(breaks=seq(-200,300,50)) +
  xlab("Duration of nasal consonant (ms)") + ylab("Duration of vowel nasalization (ms)") +
  theme_bw()
dev.off()




## Does speech rate affect the coarticulatory source and effect?
coda <- c(voiceless,voiced)
stresses <- c("N","F")

subdat <- matdat[matdat$post %in% coda & matdat$stress %in% stresses, ]
subdat$voicing <- c()
subdat$voicing[subdat$post %in% voiceless] <- "voiceless"
subdat$voicing[subdat$post %in% voiced] <- "voiced"
subdat[subdat$voicing=="voiceless",]
subdat$height <- c()
subdat$height[subdat$vowel %in% high] <- "high"
subdat$height[subdat$vowel %in% low] <- "low"

subdat$dur.N  <- 1000*(subdat$velumopening_maxvel_off-subdat$Vokal_off)
subdat$dur.VN <- 1000*(subdat$Vokal_off-subdat$velumopening_maxvel_on)


subdat$DV <- subdat$dur.N

priors <- c(
  prior(normal(0, 100), class = Intercept),
  prior(normal(0, 100), class = b, coef = stressN),
  prior(cauchy(0, 10), class = sd),
  prior(cauchy(0, 10), class = sigma),
  prior(lkj(2), class = cor)
)

m1 <- brm(
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
  file = "./cache/speech_rate_vnasal",
  save_all_pars = TRUE
)

m1


# Null model (i.e., main effect not included)
priors_null <- c(
  prior(normal(0, 100), class = Intercept),
  prior(cauchy(0, 10), class = sd),
  prior(cauchy(0, 10), class = sigma)
)

m1_null <- brm(
  formula = DV ~ (1|speaker) + (1|word),
  data = subdat,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  seed = my.seed,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 20),
  prior = priors_null,
  file = "./cache/speech_rate_vnasal_null",
  save_all_pars = TRUE
)

bayes_factor(m1, m1_null)




#### Plot marginal posteriors ####

m1_post <- brms::posterior_samples(m1, pars = "b_") %>%
  dplyr::mutate(
    fast = b_Intercept,
    neutral = b_Intercept + b_stressN
  ) %>%
  dplyr::select(fast, neutral) %>%
  tidyr::gather(context, velum_on)

# You can get the 90% CI with quantile()
# 90% CI of "nd"
nd.ci <- quantile(m1_post[m1_post$context == "fast",]$velum_on, probs = c(0.05, 0.95))
# 90% CI of "nt"
nt.ci <- quantile(m1_post[m1_post$context == "neutral",]$velum_on, probs = c(0.05, 0.95))

# Plot marginal posterior distributions
pdf(file="plots/speech_rate_vnasal.pdf",width=7,height=3,onefile=T,pointsize=16)
#my.cols <- c("red","blue")
#my.cols <- viridisLite::magma(n = 3)
my.cols <- c("#2c7fb8","#7fcdbb")
m1_plot <- m1_post %>%
  ggplot(aes(velum_on, fill = context)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values=my.cols) +
  ggtitle("Effect of speech rate on duration of vowel nasalization") + 
  ylab("Posterior probability (density)") + xlab("Duration (ms)")
yrange <- ggplot_build(m1_plot)$layout$panel_scales_y[[1]]$range$range
yrange <- yrange[2] - yrange[1]
m1_plot + 
  geom_segment(aes(x=nd.ci[1],xend=nd.ci[2],y=-yrange/10,yend=-yrange/10),col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[1],xend=nd.ci[1],y=(-yrange/10 + yrange/40),yend=(-yrange/10 - yrange/40)),col=my.cols[1]) +
  geom_segment(aes(x=nd.ci[2],xend=nd.ci[2],y=(-yrange/10 + yrange/40),yend=(-yrange/10 - yrange/40)),col=my.cols[1]) +
  geom_segment(aes(x=nt.ci[1],xend=nt.ci[2],y=-yrange/20,yend=-yrange/20),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[1],xend=nt.ci[1],y=(-yrange/20 + yrange/40),yend=(-yrange/20 - yrange/40)),col=my.cols[2]) +
  geom_segment(aes(x=nt.ci[2],xend=nt.ci[2],y=(-yrange/20 + yrange/40),yend=(-yrange/20 - yrange/40)),col=my.cols[2]) +
  theme_bw() +
  scale_x_continuous(breaks=seq(-100,300,5))
dev.off()