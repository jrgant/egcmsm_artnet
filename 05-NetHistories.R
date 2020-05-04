# %% Setup --------------------------------------------------------------------

source("01-Import-Private-Data.R")

pacman::p_load(
  ggplot2,
  ggthemes,
  magrittr,
  data.table,
  MASS,
  tidyverse,
  summarytools
)

an <- fread(paste0(pd_path, "/artnet-wide-cleaned.csv"))
names(an)

# Set ggplot theme globally
theme_set(theme_classic())

# %% Partnerships -------------------------------------------------------------

# Rationale for negative binomial and BIC:

# Hamilton DT, Handcock MS, Morris M. Degree Distributions in Sexual Networks: # A Framework for Evaluating Evidence. Sex Transm Dis. 01/2008
# [cited 2016 Jun 5];35(1):30–40.

# Check Data

varsel <- c(
  "id", "race.cat", "age", "age.grp", "hiv.ego", names(an)[grep("pn", names(an))]
)

st_options(plain.ascii = T,
           dfSummary.graph.col = F,
           dfSummary.varnumbers = F)

dfSummary(an[, ..varsel])
dfSummary(an[hiv.ego == 1, ..varsel])
dfSummary(an[hiv.ego == 0, ..varsel])

an[!is.na(hiv.ego) & !is.na(pn_ongoing), summarytools::descr(pn_ongoing)]

# Function to test best-fitting distribution
# Types: Poisson, quasipoisson, Negative binomial
degree_est <- function(data, model_fml) {

  require(MASS)

  fml <- as.formula(model_fml)

  # "null" model
  pois_fit <- glm(fml, data = data, family = "poisson")

  # quasipoisson model
  qpois_fit <- glm(fml, data = data, family = "quasipoisson")

  # negative binomial model
  nbin_fit <- MASS::glm.nb(fml, data = data)

  # store
  fits <- list(pois_fit  = pois_fit,
               qpois_fit = qpois_fit,
               nbin_fit  = nbin_fit)

  bics <- c(pois_bic = BIC(pois_fit),
            qpois_bic = BIC(qpois_fit), # TODO qpois doesn't admit BIC
            nbin_bic = BIC(nbin_fit))

  bestfit <- bics[bics == min(bics)]

  results <- list(fits = fits, bic = bics, best = bestfit)
  results

}


# Generate model formulas
indvars  <- c("race.cat", "age.grp")
outcomes <- c("pnoa_12m", "pna_12m", "pnua_12m", "pn_ongoing")

rhs <- glue::glue("{ indvars[1] } * factor({ indvars[2] })")

# Subset data
an_sub <- an[, .SD, .SDcols = c(indvars, outcomes)]
names(an_sub)

# ... Estimate
deg_ests <- lapply(setNames(outcomes, outcomes), function(x) {
  fit <- degree_est(an_sub, paste(x, rhs, sep = " ~ "))
  fit
  })

names(deg_ests)


# View model summaries across partner number outcomes

summary(deg_ests$pnoa_12m$fits$qpois_fit)        # cumulative
summary(deg_ests$pna_12m$fits$qpois_fit)         # cumulative
summary(deg_ests$pnua_12m$fits$qpois_fit)        # cumulative
summary(deg_ests$pn_ongoing$fits$qpois_fit)      # ongoing

dfSummary(an[, .(race.cat, age.grp, pn_ongoing)])

# Function to predict degree for specified combinations of traits
predict_pnum <- function(yvec = outcomes, model_list, pf) {

  require(ciTools)

  pred_df <- expand.grid(
    race.cat = unique(an$race.cat),
    age.grp = unique(an$age.grp),
    hiv.ego = 0:1) %>%
    setDT %>%
    .[order(race.cat, age.grp)]

  plot_facet <- "~ pn_measure"

  # calculate confidence intervals
  predictions <- list()

  for (i in 1:length(yvec)) {

    # generate prediction dataset
    predictions[[yvec[i]]] <- as.data.table(
      add_ci(tb = pred_df,
             fit = model_list[[i]],
             names = c("ll95", "ul95"))
    )

    predictions[[yvec[i]]]$pn_measure <- yvec[i]
  }

  preds_out <- list()
  preds_out[["predictions"]] <- rbindlist(predictions)

  # plot predicted means
  preds_out[["plot"]] <- preds_out[["predictions"]] %>%
    ggplot(aes(x = age.grp,
               y = pred,
               col = race.cat)) +
      geom_line(size = 2) +
      scale_color_viridis_d() +
      facet_wrap(as.formula(pf), ncol = 2) +
      theme_clean()

  preds_out
}

# %% Partner Number Predictions ----------------------------------------

# @NOTE
# - Switched to Jones et al's age groups
# - Now including age * race.cat interaction term
# - Replaces prior TODO
# - Revisit CI

pred_ongoing <- predict_pnum(
  yvec = "pn_ongoing",
  model_list = list(deg_ests$pn_ongoing$age_categorical$fits$nbin_fit),
  pf = quo(hiv.ego)
 )

pred_ongoing$plot

#
# # %% Simplify Models -----------------------------------------------------------
#
# # ONGOING
#
# # @NOTE:
# # - No evidence of race * age interaction for ongoing partnership
# #   number (pn_ongoing)
#
# fit.ong <- MASS::glm.nb(pn_ongoing ~ race.cat + age + hiv.ego, data = an_sub)
# summary(fit.ong)
# nrow(fit.ong$model)
#
# pl.ong <- predict_pnum("pn_ongoing", list(fit.ong), quo(hiv.ego))
#
# pl.ong$plot +
#   labs(title = "Ongoing partnerships by Ego HIV status")
#
#
# temp <- "C:/Users/jason/Desktop/"
#
# # ORAL OR ANAL, 12 MONTHS
#
# fit.pnoa12m <- MASS::glm.nb(pnoa_12m ~ race.cat * age + hiv.ego, data = an_sub)
#
# summary(fit.pnoa12m)
# ccn <- nrow(fit.pnoa12m$model)
#
# pl.pnoa12m <- predict_pnum("pnoa_12m", list(fit.pnoa12m), quo(hiv.ego))
#
# pl.pnoa12m$plot +
#   labs(title = "Oral or Anal Partnerships (past 12 months), by Ego HIV Status",
#        caption = glue::glue("Based on { ccn } complete cases"))
#
# ggsave(paste0(temp, "pnoa12m.png"), width = 7, height = 4, units = "in")
#
#
# # ANAL PARTNERS, 12 MONTHS
#
# fit.pna12m <- MASS::glm.nb(pna_12m ~ race.cat * age + hiv.ego, data = an_sub)
# summary(fit.pna12m)
# ccn <- nrow(fit.pna12m$model)
#
# pl.pna12m <- predict_pnum("pna_12m", list(fit.pna12m), quo(hiv.ego))
#
# pl.pna12m$plot +
#   labs(title = "Anal Partnerships (past 12 months), by Ego HIV Status",
#        caption = glue::glue("Based on { ccn } complete cases"))
#
# ggsave(paste0(temp, "pna12m.png"), width = 7, height = 4, units = "in")
#
# # UNPROTECTED ANAL PARTNERS, 12 MONTHS
#
# fit.pnua12m <- MASS::glm.nb(pnua_12m ~ race.cat + age + hiv.ego, data = an_sub)
# summary(fit.pnua12m)
# ccn <- nrow(fit.pnua12m$model)
#
# pl.pnua12m <- predict_pnum("pnua_12m", list(fit.pno12m), quo(hiv.ego))
#
# pl.pnua12m$plot +
#   labs(
#     title = "Unprotected Anal Partnerships (past 12 months),\nby Ego HIV Status",
#     caption = glue::glue("Based on { ccn } complete cases"))
# ggsave(paste0(temp, "pnua12m.png"), width = 7, height = 4, units = "in")
