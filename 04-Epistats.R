#[allow(non_snake_case)]

################################################################################
                       ## HELPER OBJECTS AND FUNCTIONS ##
################################################################################

## h/t: https://stackoverflow.com/questions/36344415/is-there-some-way-to-not-
##      show-a-warning-for-non-snake-case-identifiers

source("01-Import-Private-Data.R")

pacman::p_load(
  data.table,
  ciTools,
  dplyr,
  stringr,
  lubridate,
  ggplot2,
  ggthemes,
  ggridges,
  viridis,
  MASS,
  rms,
  Hmisc,
  mgcv,
  survey,
  mice,
  pscl
)

# function to calculate number of unique levels in a variable
unql <- function(data) length(unique(data))

# import cleaned wide dataset
pd_path <- Sys.getenv("ARTNET_PATH")

an <- fread(file.path(pd_path, "artnet-wide-cleaned.csv"))
anl <- fread(file.path(pd_path, "artnet-long-cleaned.csv"))

imp_mc <- readRDS(file.path(pd_path, "artnet-imputed-mc-augmented.Rds"))
imp_otp <- readRDS(file.path(pd_path, "artnet-imputed-otp-augmented.Rds"))
mcdt <- as.data.table(complete(imp_mc, action = "long"))
otpdt <- as.data.table(complete(imp_otp, action = "long"))

sort(names(an))
sort(names(anl))

# set default ggplot themes
theme_set(theme_base())


################################################################################
                       ## HELPER OBJECTS AND FUNCTIONS ##
################################################################################

## NOTE
## Majority vote method described in:
## Van Buuren, S. (2018). Flexible imputation of missing data, second edition.
## http://dx.doi.org/10.1201/9780429492259
##
## Limitation: method used here is subject to overfitting.

## Identify the average values at these quantiles across imputations and use
## as common knots across imputation data sets. Very little difference in
## quantile values across the data sets, but need to use static values to
## predict on new data without error.
knots3 <- c(0.1, 0.5, 0.9)
knots5 <- c(0.05, 0.275, 0.5, 0.725, 0.95)
contvars <- c("age.i", "age.j", "abs_sqrt_agediff", "durat_wks")

get_qts <- function(data, contvars, k) {
  lapply(setNames(contvars, contvars), function(x) {
    a <- split(data[.imp > 0], by = ".imp")
    sapply(a, function(y) {
      if (x == "age.i") {
        quantile(y[, first(get(x)), id][, V1], k)
      } else {
        quantile(y[, get(x)], k)
      }
    })
  })
}

mc_qts3 <- get_qts(mcdt, contvars = contvars, knots3)
mc_qts3 <- as.data.frame(sapply(mc_qts3, rowMeans))

mc_qts5 <- get_qts(mcdt, contvars = contvars, knots5)
mc_qts5 <- as.data.frame(sapply(mc_qts5, rowMeans))

otp_qts3 <- get_qts(otpdt, contvars = contvars[-4], knots3)
otp_qts3 <- as.data.frame(sapply(otp_qts3, rowMeans))

otp_qts5 <- get_qts(otpdt, contvars = contvars[-4], knots5)
otp_qts5 <- as.data.frame(sapply(otp_qts5, rowMeans))

## Set up the combinations of spline and polynomial terms

fml_pieces_mc <- list(
  contvars = contvars,
  poly3 = paste0("poly(", contvars, ", 3)"),
  poly2 = paste0("poly(", contvars, ", 2)"),
  spline5k = paste0("rms::rcs(", contvars, ", mc_qts5$", contvars, ")"),
  spline3k = paste0("rms::rcs(", contvars, ", mc_qts3$", contvars, ")")
)

fml_pieces_otp <- list(
  contvars = contvars[-4],
  poly3 = paste0("poly(", contvars[-4], ", 3)"),
  poly2 = paste0("poly(", contvars[-4], ", 2)"),
  spline5k = paste0(
    "rms::rcs(", contvars[-4], ", otp_qts5$", contvars[-4], ")"
  ),
  spline3k = paste0("rms::rcs(", contvars[-4], ", otp_qts3$", contvars[-4], ")")
)

base_fml <- "~ factor(race.combo) + factor(hiv.concord) + any.prep"
base_fml_mc <- paste(base_fml, "+ ptype")

build_fml <- function(base, addons) {
  paste0(base, " + ", paste(addons, collapse = " + "))
}

fml_mc_upper <- sapply(
  fml_pieces_mc, . %>% build_fml(base = base_fml_mc, addons = .)
)

fml_otp_upper <- sapply(
  fml_pieces_otp, . %>% build_fml(base = base_fml, addons = .)
)

fml_mc_upper
fml_otp_upper

make_lower_fmls <- function(fml_list) {
  sp <- unlist(str_split(fml_list, " \\+ "))
  lowfml <- paste(sp[grepl("race|age\\.", sp)], collapse = " + ")
  lowfml
}

fml_mc_lower <- sapply(fml_mc_upper, make_lower_fmls)
fml_otp_lower <- sapply(fml_otp_upper, make_lower_fmls)

fml_mc_lower
fml_otp_lower

fml_mc <- mapply(
  FUN = function(x, y) list(x, y),
  x = fml_mc_upper,
  y = fml_mc_lower,
  SIMPLIFY = FALSE
)

fml_otp <- mapply(
  FUN = function(x, y) list(x, y),
  x = fml_otp_upper,
  y = fml_otp_lower,
  SIMPLIFY = FALSE
)

fml_mc
fml_otp

## This function calculates the mean AIC value for each best-fitting model
## selected by stepwise selection. The model with the lowest AIC will be chosen.
get_AIC <- function(fitlist) {
  sapply(fitlist, function(x) {
    mean(sapply(x$analyses, AIC))
  })
}

## This function tabulates the number of times predictors are chosen in the
## backwards selection procedure across imputed datasets.
count_votes <- function(models) {
  # models = the analyses item from a mira object (MICE package)
  lapply(models$analyses, formula) %>%
  lapply(., terms) %>%
  lapply(., labels) %>%
  unlist %>%
  table
}

## This function modifies an existing model object to create a new object that
## stores pooled coefficient estimates from multiply imputed datasets. These
## exported ojbects should ONLY be used with the predict() function, as other
## model components (e.g., AIC, theta, R) ARE NOT updated.
export_pooled <- function(fits, ests) {
  fit <- fits$analyses[[1]]
  coef <- ests$pooled$estimate

  fit$coefficients <- coef
  names(fit$coefficients) <- ests$pooled$term

  fit
}


################################################################################
                 ## STI TESTING FREQUENCY - NON-IMPUTED DATA ##
################################################################################

## NOTE: Only a few respondents were missing STI testing rate after data
##       cleaning. Don't need to use the imputed datasets.

## Subset.
dsti <- an[, .(id, race.string = race.cat, age, stitest_perweek_all)]
setkey(dsti, "id")

setnames(dsti, c("stitest_perweek_all"), c("stitest_perweek"))

dsti[, race.cat := match(race.string, c("black", "hispanic", "other", "white"))]
dsti[, .N, keyby = .(race.string, race.cat)]
dsti[, stitest.52 := round(stitest_perweek * 52)]

## EDA.
dsti[!is.na(stitest.52), .(mean = mean(stitest.52), var = var(stitest.52))]

dsti[
  !is.na(stitest.52) & stitest.52 > 0,
  .(mean = mean(stitest.52), var = var(stitest.52))
]

ggplot(dsti, aes(x = stitest.52)) +
  geom_histogram(color = "white")

ggplot(dsti, aes(x = factor(race.cat), y = stitest.52)) +
  geom_boxplot()

set.seed(1971)
ggplot(dsti, aes(x = age, y = stitest.52)) +
  geom_point(position = "jitter", alpha = 0.3) +
  stat_smooth()

## Models. Fit a zero-inflated Poisson model for yearly STI rate.
fit_stitest_5k <- zeroinfl(
  stitest.52 ~ factor(race.cat) + rms::rcs(age, 5),
  data = dsti,
  x = TRUE
)

fit_stitest_3k <- zeroinfl(
  stitest.52 ~ factor(race.cat) + rms::rcs(age, 3),
  data = dsti,
  x = TRUE
)

fit_stitest_3p <- zeroinfl(
  stitest.52 ~ factor(race.cat) + poly(age, 2),
  data = dsti,
  x = TRUE
)

fit_stitest_2p <- zeroinfl(
  stitest.52 ~ factor(race.cat) + poly(age, 3),
  data = dsti,
  x = TRUE
)

sapply(
  list(fit_stitest_5k, fit_stitest_3k, fit_stitest_3p, fit_stitest_2p),
  AIC
)

summary(fit_stitest_5k)

pred_stitest <- as.data.table(expand.grid(
  race.cat = unique(dsti$race.cat),
  age = unique(dsti$age)
))

setkeyv(pred_stitest, c("race.cat", "age"))

pred_stitest[, ":="(
  pred_rate = predict(
    fit_stitest_5k, newdata = pred_stitest, type = "response"
  ),
  pred_zero = predict(fit_stitest_5k, newdata = pred_stitest, type = "zero"),
  pred_pois = predict(fit_stitest_5k, newdata = pred_stitest, type = "count")
)]

psti <- melt(
  pred_stitest,
  measure.vars = c("pred_rate", "pred_zero", "pred_pois"),
  variable.name = "pred_type",
  value.name = "pred_val"
  )

ggplot(psti, aes(x = age, y = pred_val)) +
  geom_line(aes(color = factor(race.cat)), size = 1.5) +
  scale_color_viridis_d(option = "magma", end = 0.8) +
  facet_wrap(~ pred_type, scales = "free")


################################################################################
        ## ANAL SEX ACT RATES (VIZ ONLY), MAIN/CASUAL - IMPUTED DATA ##
################################################################################

## EDA. (Uses the first imputed dataset)
dmc <- mcdt[.imp == 1]

hist(dmc$ai.rate.52)
hist(sqrt(dmc$ai.rate.52))
hist(dmc$ai.rate.52^(1 / 3))

dmc[, air52sqrt := sqrt(ai.rate.52)]
dmc[, air52cbrt := ai.rate.52^(1 / 3)]

## NOTE Cube root makes the yearly anal sex act rate variable approximately
##      Poisson-distributed (but semi-continuous, so use Gamma family).
dmc[, .(
  mean.sqrt = mean(air52sqrt),
  var.sqrt = var(air52sqrt),
  mean.cbrt = mean(air52cbrt),
  var.cbrt = var(air52cbrt)
)]

ggplot(dmc, aes(x = age.i, y = age.j)) +
  geom_point(
    aes(color = air52cbrt),
    position = "jitter"
  ) +
  stat_smooth()

ggplot(dmc, aes(x = abs_sqrt_agediff, y = air52cbrt)) +
  geom_density_2d_filled() +
  geom_point(
    alpha = 0.4,
    size = 0.7,
    position = "jitter",
    color = "white"
  ) +
  stat_smooth() +
  scale_fill_viridis_d(option = "magma")

ggplot(dmc, aes(y = factor(race.combo), x = air52cbrt)) +
  geom_density_ridges()

ggplot(dmc, aes(y = factor(hiv.concord), x = air52cbrt)) +
  geom_density_ridges()

ggplot(dmc, aes(x = air52cbrt)) +
  geom_histogram(color = "white") +
  facet_wrap(~any.prep, ncol = 1)


################################################################################
                        ## MODEL SELECTION FUNCTIONS ##
################################################################################

select_model <- function(impobj, fml, outcome,
                         glmselect = c("negbin", "logit"), threshold = 10) {

  if (glmselect == "negbin") {
    fits <- lapply(fml, function(x) {
      expr <- expression(
        f1 <- glm.nb(
          as.formula(paste(outcome, x[[1]])),
          model = FALSE,
          y = FALSE
        ),
        f2 <- stepAIC(
          f1,
          scope = list(lower = x[[2]], upper = formula(f1)),
          direction = "backward"
        )
      )

      fit <- with(impobj, expr)
      fit
    })
  }

  if (glmselect == "logit") {
    fits <- lapply(fml, function(x) {
      expr <- expression(
        f1 <- glm(
          as.formula(paste(outcome, x[[1]])),
          family = binomial,
          model = FALSE,
          y = FALSE,
          control = glm.control(maxit = 100)
        ),
        f2 <- stepAIC(
          f1,
          scope = list(lower = x[[2]], upper = formula(f1)),
          direction = "backward")
      )

      fit <- with(impobj, expr)
      fit
    })
  }

  majority_votes <- lapply(fits, function(x) {
    v <- count_votes(x)
    names(v)[v > threshold]
  })

  bestfits <- lapply(majority_votes, function(x, gs = glmselect) {

    if (gs == "negbin") {
      imp <- with(
        impobj,
        glm.nb(
          as.formula(paste(outcome, "~", paste(x, collapse = " + "))),
          model = FALSE,
          y = FALSE
        )
      )
    }

    if (gs == "logit") {
      imp <- with(
        impobj,
        glm(
          as.formula(paste(outcome, "~", paste(x, collapse = " + "))),
          model = FALSE,
          y = FALSE,
          family = binomial
        ))
    }

    imp
  })

  aics <- sort(get_AIC(bestfits))
  winner <- names(aics)[aics == min(aics)]
  pooled <- pool(bestfits[[winner]])

  fit_out <- export_pooled(bestfits[[winner]], pooled)

  ## Before we save epistats, need workaround to prevent saveRDS from inflating
  ## the file size (makes a huge .Rds object).
  ##
  ## Solutions at:
  ## stackoverflow.com/questions/42230920/saverds-inflating-size-of-object
  ## blogs.oracle.com/r/is-the-size-of-your-lm-model-causing-you-headaches
  fit_out$fitted.values <- NULL
  fit_out$data <- NULL
  attr(fit_out$terms, ".Environment") <- NULL
  attr(fit_out$formula, ".Environment") <- NULL

  if (glmselect == "negbin") {
    thetas <- unlist(lapply(bestfits[[winner]]$analyses, . %>% .$theta))
    theta <- mean(thetas)
  }

  out <- list(
    stepfits = fits,
    majority_votes = majority_votes,
    bestfits = bestfits,
    AIC = aics,
    winner = winner,
    pooled = pooled,
    fit_out = fit_out
  )

  if (exists("theta")) out <- c(out, list(theta = theta))
  out
}

## This function checks that the pooled fit object generates a different
## prediction than one of the individual fit objects in an imputed data set.
check_pooled_fit <- function(selmod, pred_df) {

  fit1 <- with(selmod, bestfits[[winner]]$analyses[[1]])
  fit2 <- selmod$fit_out
  attr(fit2$terms, ".Environment") <- globalenv()

  pred1 <- predict(fit1, newdata = pred_df, type = "response")
  pred2 <- predict(fit2, newdata = pred_df, type = "response")

  list(pred_single_fit = unname(pred1), pred_pooled_fit = unname(pred2))
}


################################################################################
                 ## ANAL ACT RATES, MAIN/CAUSAL PARTNERSHIPS ##
################################################################################

ai52 <- select_model(imp_mc, fml_mc, "ai.rate.52", "negbin")

## Confirm that replacement coefficients are used correctly with predict().
pred_ai52 <- data.table(
  ptype = 2,
  any.prep = 1,
  race.combo = 11,
  age.i = 45,
  age.j = 30,
  abs_sqrt_agediff = abs(sqrt(45) - sqrt(30)),
  hiv.concord = 3,
  durat_wks = 30
)

## Prediction from model fit to first imputation
check_pooled_fit(ai52, pred_ai52)

################################################################################
                   ## CONDOM USE, MAIN/CASUAL PARTNERSHIPS ##
################################################################################

cond_mc <- select_model(imp_mc, fml_mc, outcome = "cond.prob", "logit")

pred_condmc <- data.table(
  ptype = 1,
  race.combo = 23,
  hiv.concord = 2,
  age.i = 23,
  age.j = 21,
  any.prep = 0,
  durat_wks = 4
)

## Create a model object to use with predict() in the epidemic model.
## USE ***ONLY*** FOR THAT PURPOSE!
check_pooled_fit(cond_mc, pred_condmc)


################################################################################
                 ## SEX ACT PROBABILITIES, ONE-TIME CONTACTS ##
################################################################################

ai_once <- select_model(imp_otp, fml_otp, outcome = "ai_once", "logit")

pred_aionce <- data.table(
  race.combo = 11,
  age.i = 45,
  age.j = 60,
  any.prep = 0,
  hiv.concord = 2
)

## Create a model object to use with predict() in the epidemic model.
## USE ***ONLY*** FOR THAT PURPOSE!
check_pooled_fit(ai_once, pred_aionce)


################################################################################
                 ## CONDOM USE, ONE-TIME CONTACTS ##
################################################################################

cond_otp <- select_model(imp_otp, fml_otp, outcome = "cond.prob", "logit")

pred_condotp <- data.table(
  race.combo = 14,
  age.i = 30,
  age.j = 34,
  abs_sqrt_agediff = abs(sqrt(30) - sqrt(34)),
  any.prep = 1,
  hiv.concord = 3
)

## Create a model object to use with predict() in the epidemic model.
## USE ***ONLY*** FOR THAT PURPOSE!
check_pooled_fit(cond_otp, pred_condotp)


################################################################################
             ## ORAL ACT PROBABILITIES, MAIN/CASUAL PARTNERSHIPS ##
################################################################################

oi52 <- select_model(imp_mc, fml_mc, "oi.rate.52", "negbin")

pred_oi52 <- data.table(
  any.prep = 0,
  race.combo = 14,
  ptype = 2,
  age.i = 18,
  age.j = 20,
  abs_sqrt_agediff = abs(sqrt(18) - sqrt(20)),
  durat_wks = 50
)

## Create a model object to use with predict() in the epidemic model.
## USE ***ONLY*** FOR THAT PURPOSE!
check_pooled_fit(oi52, pred_oi52)


################################################################################
                ## ORAL ACT PROBABILITIES, ONE-TIME CONTACTS ##
################################################################################

oi_once <- select_model(imp_otp, fml_otp, "oi_once", "logit")

pred_oionce <- data.table(
  hiv.concord = 1,
  race.combo = 11,
  age.i = 45,
  age.j = 37,
  abs_sqrt_agediff = abs(sqrt(45) - sqrt(37)),
  any.prep = 1
)

## Create a model object to use with predict() in the epidemic model.
## USE ***ONLY*** FOR THAT PURPOSE!
check_pooled_fit(oi_once, pred_oionce)


################################################################################
                    ## HIV DIAGNOSIS STATUS MODEL (SEED) ##
################################################################################

## Select a model with which to seed the population. Because HIV diagnosis
## status in the seeded population is used to calculate network target
## statistics, we should ## aim for an accurate enough seed diagnosis prevalence
## so that those network statistics are well-estimated.

hiv_matchrole <- anl[, .(id, role.class = ego.anal.role)]
hiv_matchrole[role.class == "", role.class := NA]

hiv_mr <- dcast(
  hiv_matchrole,
  id ~ role.class,
  function(x) as.numeric(length(x) > 0)
)

hiv_mr[, role.class := fcase(
           Versatile == 1 | (Insertive + Receptive > 1), "Versatile",
           Insertive == 1 & Receptive == 0, "Insertive",
           Receptive == 1 & Insertive == 0, "Receptive",
           default = NA
         )]

hiv_mr[, .N, role.class]
hiv_mr[, .N, keyby = .(role.class, Insertive, Receptive, Versatile)]

hivdt <- an[
  hiv_mr[, .(id, role.class)],
  on = "id"
][, .(id, hiv2, race.cat, age, age.grp, role.class)]

hivdt[, .(N = .N, pr_hiv = mean(hiv2)), .(role.class, race.cat, age)] %>%
  ggplot(aes(x = age, y = pr_hiv, color = race.cat, fill = race.cat)) +
  stat_smooth(se = FALSE) +
  geom_point(aes(size = N), shape = 21, color = "white") +
  facet_wrap(~role.class) +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

hiv_mod_5k <- glm(
  hiv2 ~ race.cat + rms::rcs(age, 5) + factor(role.class),
  data = hivdt,
  family = binomial
)

hiv_mod_3k <- glm(
  hiv2 ~ race.cat + rms::rcs(age, 3) + factor(role.class),
  data = hivdt,
  family = binomial
)

hiv_mod_poly3 <- glm(
  hiv2 ~ race.cat + poly(age, 3) + factor(role.class),
  data = hivdt,
  family = binomial
)

hiv_mod_poly2 <- glm(
  hiv2 ~ race.cat + poly(age, 2) + factor(role.class),
  data = hivdt,
  family = binomial
)

hiv_mod_linear <- glm(
  hiv2 ~ race.cat + age + factor(role.class),
  data = hivdt,
  family = binomial
)

hiv_mod_cat <- glm(
  hiv2 ~ race.cat + factor(age.grp) + factor(role.class),
  data = hivdt,
  family = binomial
)

hivmods <- list(
  hiv_mod_5k = hiv_mod_5k,
  hiv_mod_3k = hiv_mod_3k,
  hiv_mod_poly3 = hiv_mod_poly3,
  hiv_mod_poly2 = hiv_mod_poly2,
  hiv_mod_linear = hiv_mod_linear,
  hiv_mod_cat = hiv_mod_cat
)

## Running into probable sparsity issues with 5-knot model. Compare against
## a 3-knot model.

## Based on likelihood ratio tests, poly2 seems to be the most parsimonious.
lrtest(hiv_mod_5k, hiv_mod_3k)
lrtest(hiv_mod_3k, hiv_mod_poly3)
lrtest(hiv_mod_poly3, hiv_mod_poly2)
lrtest(hiv_mod_poly2, hiv_mod_linear)
lrtest(hiv_mod_linear, hiv_mod_cat)

## AIC selects the 3-knot model (poly2 a close second)
sapply(hivmods, AIC, simplify = FALSE)

## To avoid having to export the knot quantiles, and because the models
## performed similarly, choose the poly2 model.

pred_hiv <- data.table(
  expand.grid(
    race.cat = unique(hivdt$race.cat),
    age = unique(hivdt$age),
    role.class = c("Versatile", "Receptive", "Insertive")
  )
)

setorder(pred_hiv, "age")

pred_hiv[, predhiv := predict(hiv_mod_poly2, newdata = pred_hiv, type = "response")]
pred_hiv

ggplot(pred_hiv, aes(x = age, y = predhiv, color = race.cat)) +
  geom_line() +
  facet_wrap(~ role.class) +
  scale_color_viridis_d(option = "magma", end = 0.9)


################################################################################
## WRITE TO FILE
################################################################################

epistats <- list()

epistats$ai.acts.mc <- ai52$fit_out
epistats$ai.acts.mc.theta <- ai52$theta
epistats$ai.acts.oo <- ai_once$fit_out
epistats$mc_qts <- list(q5 = mc_qts5, q3 = mc_qts3)
epistats$oi.acts.mc <- oi52$fit_out
epistats$oi.acts.mc.theta <- oi52$theta
epistats$oi.acts.oo <- oi_once$fit_out
epistats$otp_qts <- list(q5 = otp_qts5, q3 = otp_qts3)
epistats$cond.mc <- cond_mc$fit_out
epistats$cond.oo <- cond_otp$fit_out
epistats$stitest <- fit_stitest_5k
epistats$hiv.mod <- hiv_mod_poly2


saveRDS(epistats, file.path("netstats", "epistats.Rds"))
