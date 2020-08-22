# %% LOAD ----------------------------------------------------------------------

source("01-Import-Private-Data.R")

pacman::p_load(
  data.table,
  summarytools,
  ciTools,
  dplyr,
  stringr,
  lubridate,
  ggplot2,
  ggthemes,
  viridis,
  MASS,
  rms,
  Hmisc,
  mgcv,
  survey
)

# function to calculate number of unique levels in a variable
unql <- function(data) length(unique(data))

# import cleaned wide dataset
pd_path <- Sys.getenv("ARTNET_PATH")

an <- fread(file.path(pd_path, "artnet-wide-cleaned.csv"))
anl <- fread(file.path(pd_path, "artnet-long-cleaned.csv"))
imp <- readRDS(file.path(pd_path, "artnet-long-imputed.Rds"))

sort(names(an))
sort(names(anl))

# set default ggplot themes
theme_set(theme_base())


# %% STI TESTING FREQUENCY - MSM NEVER ON PREP  --------------------------------

ra_grid <-
  as.data.table(expand.grid(
    race.cat = unique(an$race.cat),
    age.grp = unique(an$age.grp),
    stringsAsFactors = FALSE
  ))[order(age.grp, race.cat)]

print(ra_grid)


## NEVER ON PREP

noprep <- an[prep_revised != 1, .(
  id, age, race.cat, age.grp, prep_revised, stitest_2yr
  )][order(id)]

print(noprep)

prep <- an[prep_revised == 1, .(
  id, age, race.cat, age.grp, prep_revised, stitest_2yr_prep
  )][order(id)]

print(prep)

# 2 missing stitest_2yr info
mice::md.pattern(noprep)
mice::md.pattern(prep)

# distribution of number of STI tests (any), past 2 years
hist(noprep$stitest_2yr)
summary(noprep$stitest_2yr)

noprep[, .(
  mean = mean(stitest_2yr, na.rm = T),
  var = var(stitest_2yr, na.rm = T)
)]

fit_stitest_2yr_noprep <- MASS::glm.nb(
  stitest_2yr ~ race.cat + factor(age.grp),
  data = noprep
)

summary(fit_stitest_2yr_noprep)

print(broom::tidy(fit_stitest_2yr_noprep, conf.int = TRUE))

stitest_noprep_pred <- cbind(ra_grid,
  pred_2yr = predict(
    fit_stitest_2yr_noprep,
    newdata = ra_grid,
    type = "response"
    )
  ) %>%
  ciTools::add_ci(.,
    fit_stitest_2yr_noprep,
    names = c("pred_2yr_ll95", "pred_2yr_ul95")
  ) %>%
  setDT %>%
  .[, ":="(
      wkrate = pred_2yr / (52 * 2),
      wkrate_ll95 = pred_2yr_ll95 / (52 * 2),
      wkrate_ul95 = pred_2yr_ul95 / (52 * 2)
    )]

print(stitest_noprep_pred)

ggplot(
  stitest_noprep_pred,
  aes(x = factor(age.grp),
      y = pred_2yr,
      group = race.cat,
      color = race.cat,
      fill = race.cat)
  ) +
  geom_pointrange(
    aes(ymin = pred_2yr_ll95,
        ymax = pred_2yr_ul95),
    size = 1) +
  scale_color_viridis_d() +
  geom_line(size = 1) +
  labs(title = "STI tests among never-PrEP MSM")


# %% STI TESTING FREQUENCY - MSM EVER ON PREP  ---------------------------------

# STI tests outside of PrEP follow-up
fit_stitest_2yr_prep <- glm.nb(
  stitest_2yr_prep ~ race.cat + factor(age.grp),
  data = prep
)

summary(fit_stitest_2yr_prep)

stitest_prep_pred <- cbind(ra_grid,
  pred_2yr = predict(
    fit_stitest_2yr_prep,
    newdata = ra_grid,
    type = "response"
    )
  ) %>%
  ciTools::add_ci(.,
    fit_stitest_2yr_prep,
    names = c("pred_2yr_ll95", "pred_2yr_ul95")
  ) %>%
  setDT %>%
  .[, ":="(
      wkrate = pred_2yr / (52 * 2),
      wkrate_ll95 = pred_2yr_ll95 / (52 * 2),
      wkrate_ul95 = pred_2yr_ul95 / (52 * 2)
    )]

print(stitest_prep_pred)

ggplot(
  stitest_prep_pred,
  aes(x = factor(age.grp),
      y = pred_2yr,
      group = race.cat,
      color = race.cat,
      fill = race.cat)
  ) +
  geom_pointrange(
    aes(ymin = pred_2yr_ll95,
        ymax = pred_2yr_ul95),
    size = 1,
    alpha = 0.4) +
  scale_color_viridis_d() +
  geom_line(size = 1) +
  labs(title = "Non-PrEP STI tests among ever-PrEP MSM")


# %% SEX ACT RATES ------------------------------------------

plot_actrates_main <- function(data, var, filterby) {

  selcols <- c("id", "pid", "ptype", filterby, var)

  data[get(filterby) == 1, ..selcols] %>%
    ggplot(aes_string(x = var)) +
    geom_bar(stat = "bin", color = "white", fill = "#990000") +
    geom_rug(color = "gray", alpha = 0.5) +
    facet_wrap(~ ptype)

}

p.recai.rate <- plot_actrates_main(anl[ptype %in% 1:2], "recai.rate", "p_rai")
p.insai.rate <- plot_actrates_main(anl[ptype %in% 1:2], "insai.rate", "p_iai")
p.recoi.rate <- plot_actrates_main(anl[ptype %in% 1:2], "recoi.rate", "p_roi")
p.insoi.rate <- plot_actrates_main(anl[ptype %in% 1:2], "insoi.rate", "p_ioi")

# 14.4% of eligible responses missing
anl[ptype %in% 1:2 & p_rai == 1, .(
  missing = sum(is.na(recai.rate)),
  n_elig = .N
  )][, pct := round(missing / n_elig * 100, digits = 2)][]

# 14.8% of eligible responses missing
anl[ptype %in% 1:2 & p_iai == 1, .(
  missing = sum(is.na(insai.rate)),
  n_elig = .N
)][, pct := round(missing / n_elig * 100, digits = 2)][]

# 15.3% of eligible responses missing
anl[ptype %in% 1:2 & p_roi == 1, .(
  missing = sum(is.na(recoi.rate)),
  n_elig = .N
)][, pct := round(missing / n_elig * 100, digits = 2)][]

# 35.5% of eligible responses missing
anl[ptype %in% 1:2 & p_ioi == 1, .(
  missing = sum(is.na(insoi.rate)),
  n_elig = .N
)][, pct := round(missing / n_elig * 100, digits = 2)][]


# %% SEX ACT PROBABILITY, ONE-TIME CONTACTS ------------------------------------

sort(names(anl))
nrow(anl[ptype == 3])

nmiss_op_rai <- anl[ptype == 3, sum(is.na(p_rai_once))]
nmiss_op_iai <- anl[ptype == 3, sum(is.na(p_iai_once))]
nmiss_op_roi <- anl[ptype == 3, sum(is.na(p_roi_once))]
nmiss_op_ioi <- anl[ptype == 3, sum(is.na(p_ioi_once))]

data.table(
  nmiss_op_rai,
  nmiss_op_iai,
  nmiss_op_roi,
  nmiss_op_ioi
)[]

# all missing _once variables accompanied by numeric entry to
#  p_actsprefernot or p_actsdk (meaning non-structurally missing)
anl[ptype == 3, .N, keyby = .(p_rai_once, p_actsprefernot, p_actsdk)]
anl[ptype == 3, .N, keyby = .(p_iai_once, p_actsprefernot, p_actsdk)]
anl[ptype == 3, .N, keyby = .(p_roi_once, p_actsprefernot, p_actsdk)]
anl[ptype == 3, .N, keyby = .(p_ioi_once, p_actsprefernot, p_actsdk)]

anl[ptype == 3 & !is.na(p_rai_once), .N, keyby = p_rai_once
  ][, P := round(N / sum(N), 3)][]

anl[ptype == 3 & !is.na(p_iai_once), .N, keyby = p_iai_once
  ][, P := round(N / sum(N), 3)][]

anl[ptype == 3 & !is.na(p_roi_once), .N, keyby = p_roi_once
  ][, P := round(N / sum(N), 3)][]

anl[ptype == 3 & !is.na(p_ioi_once), .N, keyby = p_ioi_once
  ][, P := round(N / sum(N), 3)][]


# %% CONDOM USE ---------------------------------------------

ggplot(anl[ptype %in% 1:2]) +
  geom_density(aes(x = insuai.prob, color = "insuai.prob")) +
  geom_density(aes(x = recuai.prob, color = "recuai.prob")) +
  geom_density(aes(x = cond.prob, color = "cond.prob"), size = 1) +
  facet_wrap(~ ptype) +
  scale_color_colorblind() +
  labs(
    x = "probability of condom use",
    title = "Main and casual partnerships"
  )


# %% MODELS --------------------------------------------------------------------

## convert character vars to match epidemic model
anhiv <- copy(an)

anhiv[, ":="(
  race.cat = as.numeric(as.factor(race.cat)),
  age.grp = as.numeric(as.factor(age.grp))
)]

## HIV diagnosis status
anhiv[, .N, .(hiv2, hiv3, hiv.ego)]
anhiv[, .N, .(artnetstatus, artnetrcntrslt, artnetevrpos, hiv.ego)]

### set unknown HIV status to missing for the HIV model only
anhiv[hiv.ego == 2, hiv.ego := NA]
anhiv[, .N, hiv.ego]

### missing data weight numerator
ipw.num.s1 <- anhiv[, sum(!is.na(hiv.ego)) / nrow(an)]
ipw.num.s0 <- 1 - ipw.num.s1

# missing data weight denominator model
hivfit.cube.age <- glm(
  !is.na(hiv.ego) ~ race.cat * age + I(age^2) + I(age^3),
  family = binomial(),
  data = anhiv
)

hivfit.quad.age <- glm(
  !is.na(hiv.ego) ~ race.cat * age + I(age^2),
  family = binomial(),
  data = anhiv
)

# cubic term for age gives lower AIC
hivfit.quad.age$aic
hivfit.cube.age$aic

hiv.s1.pred <- predict(hivfit.cube.age, newdata = anhiv, type = "response")

anhiv[, ":="(
  ipw.num = ifelse(is.na(hiv.ego), ipw.num.s0, ipw.num.s1),
  ipw.den = ifelse(is.na(hiv.ego), 1 - hiv.s1.pred, hiv.s1.pred)
)][, ipsw := ipw.num / ipw.den]

### check weights
anhiv[!is.na(hiv.ego), summary(ipsw)]
anhiv[!is.na(hiv.ego), sum(ipsw)]

### run outcome model
hivdes <- svydesign(
  data = anhiv[!is.na(hiv.ego), .(id, hiv.ego, age, race.cat, ipsw)],
  ids = ~id,
  weights = ~ipsw
)

hiv.mod.cube.age <- svyglm(
  hiv.ego ~ factor(race.cat) * age + I(age^2) + I(age^3),
  family = binomial(),
  design = hivdes
)

hiv.mod.quad.age <- svyglm(
  hiv.ego ~ factor(race.cat) + age + I(age^2),
  family = binomial(),
  design = hivdes
)

# quadratic term for age gives lower AIC
hiv.mod.quad.age$aic
hiv.mod.cube.age$aic

summary(hiv.mod.quad.age)

# compare results of complete case analysis with IPW
cc.hiv.quad <- glm(
  hiv.ego ~ factor(race.cat) + age + I(age^2),
  family = binomial(),
  data = anhiv
)

summary(cc.hiv.quad)

data.table(
  complete_cases = coef(cc.hiv.quad),
  ipw = coef(hiv.mod.quad.age)
)

pred.hiv.df <- as.data.table(
  expand.grid(race.cat = 1:4, age = 18:65)
)

pred.hiv.df2 <- data.table(
  race.cat = pred.hiv.df$race.cat,
  age = pred.hiv.df$age,
  ccpreds.hiv = predict(
    cc.hiv.quad,
    newdata = pred.hiv.df,
    type = "response"
  ),
  ipwpreds.hiv = predict(
    hiv.mod.quad.age,
    newdata = pred.hiv.df,
    type = "response"
  )
)

ggplot(pred.hiv.df2, aes(x = age, color = factor(race.cat))) +
  geom_line(aes(y = ccpreds.hiv, linetype = "CC")) +
  geom_line(aes(y = ipwpreds.hiv, linetype = "IPW")) +
  labs(y = "prediction") +
  ggtitle("HIV status") +
  scale_color_colorblind() +
  theme_tufte()

hiv.mod <- hiv.mod.quad.age


## %% SEX ACT MODEL REUSABLES --------------------------------------------------

anl[ptype == 1, summary(ai.rate)]
anl[ptype == 2, summary(ai.rate)]

anl.newvars <- copy(anl)

anl.newvars[p_race.cat == "", p_race.cat := NA]

anl.newvars <- anl.newvars[, ":="(
  race.i = as.numeric(as.factor(ego.race.cat)),
  race.j = as.numeric(as.factor(p_race.cat)),
  age.i = ego.age,
  age.j = p_age_imputed,
  diag.status.i = ego.hiv,
  diag.status.j = p_hiv2
)]

# Check race/eth factor coding
anl.newvars[, .N, keyby = .(ego.race.cat, race.i)]
anl.newvars[, .N, keyby = .(p_race.cat, race.j)]


# Subset to ongoing main and casual partnerships (18 missing).
# Name variables according to use in epidemic simulation modules.
maincas <- anl.newvars[ptype %in% 1:2 & p_ongoing_ind == 1, .(
  ai.rate,
  oi.rate,
  cond.prob,
  race.i,
  age.i,
  diag.status.i,
  race.j,
  age.j,
  diag.status.j,
  abs_sqrt_agediff,
  ptype,
  durat_wks
)]


plot(density(maincas$ai.rate * 52, na.rm = TRUE))
summary(maincas$ai.rate * 52)

plot(density(log(maincas$ai.rate * 52), na.rm = TRUE))
boxplot(maincas[(ai.rate * 52) < 3000, ai.rate * 52])

plot(density(maincas$oi.rate * 52, na.rm = TRUE))
summary(maincas$oi.rate * 52)


# Model formulas

ai.outcome <- "floor(ai.rate * 52)"
ai.log.outcome <- "log(ai.rate * 52)"

oi.outcome <- "floor(oi.rate * 52)"
oi.log.outcome <- "log(oi.rate * 52)"

base.fml <-
  "factor(race.i) +
   factor(race.j) +
   age.i +
   age.j +
   abs_sqrt_agediff +
   factor(ptype) +
   factor(diag.status.i) +
   factor(diag.status.j) +
   durat_wks"

poly.fml <- . ~ . -
   I(age.i^2) + I(age.i^3) +
   I(age.j^2) + I(age.j^3) +
   I(abs_sqrt_agediff^2) + I(abs_sqrt_agediff^3) +
   I(durat_wks^2) + I(durat_wks^3)

spline.fml <- . ~ . -
  rcs(age.i) +
  rcs(age.j) +
  rcs(abs_sqrt_agediff) +
  rcs(durat_wks)


## %% ANAL SEX ACT MODEL -------------------------------------------------------

# NOTE: See Harrell /Regression Modeling Strategies/ for recommendation re:
#       backward selection.


#### Negative binomial, pick between base model and polynomials
base.fit.ai.nb <- glm.nb(
  as.formula(paste(ai.outcome, "~", base.fml)),
  data = maincas
)

poly.fit.ai.nb <- update(base.fit.ai.nb, poly.fml)

nb.ai.stepaic <- stepAIC(
  poly.fit.ai.nb,
  direction = "backward",
  scope = list(
    upper = poly.fit.ai.nb,
    lower = base.fit.ai.nb
  )
)

#### Negative binomial, pick between base model and spline
spline.fit.ai.nb <- update(base.fit.ai.nb, spline.fml)

nb.ai.stepaic2 <- stepAIC(
  spline.fit.ai.nb,
  direction = "backward",
  scope = list(
    upper = spline.fit.ai.nb,
    lower = base.fit.ai.nb
  )
)

#### Poisson, pick between base and polynomial model
base.fit.ai.pois <- glm(
  as.formula(paste(ai.outcome, "~", base.fml)),
  data = maincas,
  family = poisson()
)

poly.fit.ai.pois <- update(base.fit.ai.pois, poly.fml)

pois.ai.stepaic <- stepAIC(
  poly.fit.ai.pois,
  direction = "backward",
  scope = list(
    upper = poly.fit.ai.pois,
    lower = base.fit.ai.pois
  )
)

#### Poisson, pick between base model and spline
spline.fit.ai.pois <- update(base.fit.ai.pois, spline.fml)

pois.ai.stepaic2 <- stepAIC(
  spline.fit.ai.pois,
  direction = "backward",
  scope = list(
    upper = spline.fit.ai.pois,
    lower = base.fit.ai.pois
  )
)

#### Gaussian with transformed outcome, base vs. polynomial model
base.fit.ai.loglm <- glm(
  as.formula(paste(ai.log.outcome, "~", base.fml)),
  data = maincas,
  family = gaussian
)

poly.fit.ai.loglm <- update(base.fit.ai.loglm, poly.fml)

loglm.ai.stepaic <- stepAIC(
  poly.fit.ai.loglm,
  direction = "backward",
  scope = list(
    upper = poly.fit.ai.loglm,
    lower = base.fit.ai.loglm
  )
)

spline.fit.ai.loglm <- update(base.fit.ai.loglm, spline.fml)

loglm.ai.stepaic2 <- stepAIC(
  spline.fit.ai.loglm,
  direction = "backward",
  scope = list(
    upper = spline.fit.ai.loglm,
    lower = base.fit.ai.loglm
  )
)

ai.modlist <- list(
  nb_polyAIC_selected      = nb.ai.stepaic,
  nb_splineAIC_selected    = nb.ai.stepaic2,
  pois_polyAIC_selected    = pois.ai.stepaic,
  pois_splineAIC_selected  = pois.ai.stepaic2,
  loglm_polyAIC_selected   = loglm.ai.stepaic,
  loglm_splineAIC_selected = loglm.ai.stepaic2
)

## Plot residuals
par(mfrow = c(3, 2))
for (i in seq_along(ai.modlist)) {
  hist(resid(ai.modlist[[i]]), main = names(ai.modlist[i]), xlab = "residuals")
}

## Plot simulated responses
s.ai <- list()

par(mfrow = c(3, 2))
for (i in seq_along(ai.modlist)) {
  s.ai[[i]] <- unlist(simulate(ai.modlist[[i]], seed = 2000, type = "response"))

  if (grepl("^loglm", names(ai.modlist[i]))) {
    s.ai[[i]] <- exp(s.ai[[i]])
  }

  hist(s.ai[[i]] / 52, main = names(ai.modlist[i]), xlab = "simulated (n = 1)")
}

names(s.ai) <- names(ai.modlist)

sapply(s.ai, function(x) {
  transf <- round(x / 52, 3)
  summary(transf)
})

# Compare models on AIC
ai.aic <- sapply(ai.modlist, AIC)

# View summary for selected model
selected.ai.mod <- names(ai.aic)[ai.aic == min(ai.aic)]
summary(ai.modlist[[selected.ai.mod]])


## %% ORAL ACTS MODEL ----------------------------------------------------------

#### Pick between base (main effects) and polynomial model
base.fit.oi.nb <- glm.nb(
  as.formula(paste(oi.outcome, "~", base.fml)),
  data = maincas
)

poly.fit.oi.nb <- update(base.fit.oi.nb, poly.fml)

nb.oi.stepaic <- stepAIC(
  poly.fit.oi.nb,
  direction = "backward",
  scope = list(
    upper = poly.fit.oi.nb,
    lower = base.fit.oi.nb
  )
)

#### Pick between base (main effects) and spline model
spline.fit.oi.nb <- update(base.fit.oi.nb, spline.fml)

nb.oi.stepaic2 <- stepAIC(
  spline.fit.oi.nb,
  direction = "backward",
  scope = list(
    upper = spline.fit.oi.nb,
    lower = base.fit.oi.nb
  )
)

#### Poisson, pick between base and polynomial model
base.fit.oi.pois <- glm(
  as.formula(paste(oi.outcome, "~", base.fml)),
  data = maincas,
  family = poisson()
)

poly.fit.oi.pois <- update(base.fit.oi.pois, poly.fml)

pois.oi.stepaic <- stepAIC(
  poly.fit.oi.pois,
  direction = "backward",
  scope = list(
    upper = poly.fit.oi.pois,
    lower = base.fit.oi.pois
  )
)

#### Poisson, pick between base model and spline
spline.fit.oi.pois <- update(base.fit.oi.pois, spline.fml)

pois.oi.stepaic2 <- stepAIC(
  spline.fit.oi.pois,
  direction = "backward",
  scope = list(
    upper = spline.fit.oi.pois,
    lower = base.fit.oi.pois
  )
)

#### Gaussian model with log-transformed outcome

# Two records have weekly oral sex rate estimate of 0. Set to 0.001 for
# log-transformation.

maincas[oi.rate == 0, oi.rate := 0.001]

base.fit.oi.loglm <- glm(
  as.formula(paste(oi.log.outcome, "~", base.fml)),
  data = maincas,
  family = gaussian
)

poly.fit.oi.loglm <- update(base.fit.oi.loglm, poly.fml)

loglm.oi.stepaic <- stepAIC(
  poly.fit.oi.loglm,
  direction = "backward",
  scope = list(
    upper = poly.fit.oi.loglm,
    lower = base.fit.oi.loglm
  )
)

spline.fit.oi.loglm <- update(base.fit.oi.loglm, spline.fml)

loglm.oi.stepaic2 <- stepAIC(
  spline.fit.oi.loglm,
  direction = "backward",
  scope = list(
    upper = spline.fit.oi.loglm,
    lower = base.fit.oi.loglm
  )
)

## Compare and select models

oi.modlist <- list(
  nb_polyAIC_selected      = nb.oi.stepaic,
  nb_splineAIC_selected    = nb.oi.stepaic2,
  pois_polyAIC_selected    = pois.oi.stepaic,
  pois_splineAIC_selected  = pois.oi.stepaic2,
  loglm_polyAIC_selected   = loglm.oi.stepaic,
  loglm_splineAIC_selected = loglm.oi.stepaic2
)

## Plot residuals
par(mfrow = c(3, 2))
for (i in seq_along(oi.modlist)) {
  hist(resid(oi.modlist[[i]]), main = names(oi.modlist[i]), xlab = "residuals")
}

## Plot simulated responses
s.oi <- list()

par(mfrow = c(3, 2))
for (i in seq_along(oi.modlist)) {
  s.oi[[i]] <- unlist(simulate(oi.modlist[[i]], seed = 2000, type = "response"))
  if (grepl("^loglm", names(oi.modlist[i]))) s.oi[[i]] <- exp(s.oi[[i]])

  hist(s.oi[[i]] / 52, main = names(oi.modlist[i]), xlab = "simulated (n = 1)")
}

## Compare oral act models
names(s.oi) <- names(oi.modlist)

sapply(s.oi, function(x) {
  transf <- x / 52
  summary(transf)
})

## Summary of selected model
oi.aic <- sapply(oi.modlist, AIC)
selected.oi.mod <- names(oi.aic)[oi.aic == min(oi.aic)]
summary(oi.modlist[[selected.oi.mod]])


## %% CONDOM MODELS ------------------------------------------------------------

#### Pick between base (main effects) and polynomial model
base.fit.cp <- glm(
  as.formula(paste("cond.prob", "~", base.fml)),
  data = maincas,
  family = binomial
)

poly.fit.cp <- update(base.fit.cp, poly.fml)

spline.fit.cp <- update(base.fit.cp, spline.fml)

gamk40.fit.cp <- gam(
  cond.prob ~
    factor(race.i) +
    factor(race.j) +
    s(age.i, k = 40) +
    s(age.j, k = 40) +
    s(abs_sqrt_agediff, k = 40) +
    s(durat_wks, k = 40) +
    factor(diag.status.i) +
    factor(diag.status.j) +
    factor(ptype),
  data = maincas,
  family = binomial
)

gamk50.fit.cp <- gam(
  cond.prob ~
    factor(race.i) +
    factor(race.j) +
    s(age.i, k = 50) +
    s(age.j, k = 50) +
    s(abs_sqrt_agediff, k = 50) +
    s(durat_wks, k = 50) +
    factor(diag.status.i) +
    factor(diag.status.j) +
    factor(ptype),
  data = maincas,
  family = binomial
)

cp.modlist <- list(
  cp.mod.quad = base.fit.cp,
  cp.mod.cubic = poly.fit.cp,
  cp.mod.spline = spline.fit.cp,
  cp.mod.gamk40 = gamk40.fit.cp,
  cp.mod.gamk50 = gamk50.fit.cp
)

# Restricted cubic spline model has lowest AIC
cp.aic <- sapply(setNames(cp.modlist, names(cp.modlist)), AIC)
selected.cp.mod <- names(cp.aic)[cp.aic == min(cp.aic)]
summary(cp.modlist[[selected.cp.mod]])


################################################################################
## ONE-TIME PARTNERSHIPS ##
################################################################################

otp <- anl.newvars[ptype == 3][, .(
  id, pid, pid_unique,
  race.i, race.j, age.i, age.j,
  p_iai_once, p_rai_once,
  p_ioi_once, p_roi_once,
  p_insuai_once, p_recuai_once,
  abs_sqrt_agediff
)]

head(otp)
nrow(otp)

otp[, .N, keyby = p_iai_once]
otp[, .N, keyby = p_rai_once]
otp[, .N, keyby = .(p_iai_once, p_rai_once)]
otp[, .N, keyby = .(p_ioi_once, p_roi_once)]
otp[, .N, keyby = .(p_insuai_once, p_recuai_once)]

sapply(otp, function(x) sum(is.na(x)) / nrow(otp))

## Create multinomial outcome
## (1 = ego insertive only, 2 = ego receptive only, 3 = both acts)
otp[, ":="(
  ai.acts = p_iai_once + p_rai_once,
  oi.acts = p_ioi_once + p_roi_once
)]

otp[, .N, keyby = ai.acts]
otp[, .N, keyby = oi.acts]

otp[, .N, keyby = .(p_iai_once, p_rai_once, ai.acts)]
otp[, .N, keyby = .(p_ioi_once, p_roi_once, oi.acts)]

otp[, .N, keyby = .(p_iai_once, p_rai_once, p_insuai_once, p_recuai_once)]



################################################################################
## WRITE TO FILE
################################################################################

epistats <- list()

epistats$ai.acts.mod <- ai.modlist[[selected.ai.mod]]
epistats$oi.acts.mod <- oi.modlist[[selected.oi.mod]]
epistats$cp.mod <- cp.modlist[[selected.cp.mod]]
epistats$hiv.mod <- hiv.mod

saveRDS(epistats, file.path("netstats", "epistats.Rds"))


# debug versions of models without spline terms
epistats_debug <- list()

epistats_debug$ai.acts.mod <- ai.modlist[["loglm_polyAIC_selected"]]
epistats_debug$oi.acts.mod <- oi.modlist[["loglm_polyAIC_selected"]]
epistats_debug$cp.mod <- cp.modlist[["cp.mod.cubic"]]
epistats_debug$hiv.mod <- hiv.mod

saveRDS(epistats_debug, file.path("netstats", "epistats_debug.Rds"))
