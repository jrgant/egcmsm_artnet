# %% Setup -------------------------------------------------------------------

pacman::p_load(
  data.table,
  ggplot2,
  ggthemes,
  ggridges,
  magrittr,
  summarytools,
  dplyr,
  haven
)

theme_set(theme_classic())

an <- fread("artnet-cleaned.csv")
head(an)

anyDuplicated(names(an))

# %% PrEP ----------------------------------------------------------------------

an[prep_revised == 1, freq(prep_hivtestfreq)]   %>% print
an[prep_revised == 1, freq(prep_stithroatfreq)] %>% print
an[prep_revised == 1, freq(prep_stirectfreq)]   %>% print
an[prep_revised == 1, freq(prep_stiurethfreq)]  %>% print

# %% Descriptive Statistics and Derived Variables ------------------------------

# Age
an[, descr(age)] %>% print
plot(density(an$age))

an[, freq(artnetevertest)] %>% print


# Create and view race/ethnicity variable
an[, .N, by = .(race, race.cat)][order(race.cat)] %>% print

freq(an$race.cat) %>% print



# %% HIV by Demo --------------------------------------------------------------

# .. Race
freq(an$hiv) %>% print
an[, stby(hiv, race.cat, freq)]


# .. Age
hiv_by_age <- an %>%
  .[!is.na(hiv)] %>%
  .[, .(.N, hiv_n = sum(hiv)), keyby = age]

print(hiv_by_age)

hiv_by_age_sum <- cbind(hiv_by_age,
      hiv_by_age[, Hmisc::binconf(x = hiv_n, n = N, return.df = T)])

print(hiv_by_age_sum)

hiv_by_age_sum %>%
  ggplot(aes(x = age, y = PointEst)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3)

fit <- glm(hiv ~ age * race, data = an, family = "binomial")
summary(fit)


# %% STI Testing by Demo - Never PrEP Takers ----------------------------------

an[prep_revised == 0, stby(stitest_2yr_cat, race.cat, freq)]

an[prep_revised == 0, stby(stitest_2yr_sympt_pct, race.cat, freq)]

an[prep_revised == 0, stby(stitestfreq_cat, race.cat, freq)]


# Overall number of STI tests sought due to presence of symptoms

stisympt <- an[
  prep_revised == 0 & !is.na(stitest_2yr) & !is.na(stitest_2yr_sympt),
  .(id, stitest_2yr, stitest_2yr_sympt)
  ]

bootsti <- function(data, n) {
  sampdat <- data[sample(1:n, replace = T)]
  anlysdat <- sampdat[, .(.N,
                        pct_sympt = sum(stitest_2yr_sympt) / sum(stitest_2yr))]
  anlysdat
}


set.seed(1971)
bootdist_stisympt <- replicate(
  n = 10000,
  expr = do.call("bootsti",
                 args = list(data = stisympt,
                             n = nrow(stisympt))),
  simplify = F) %>%
  rbindlist

estlabs <- c("Prob", "LL95", "UL95")

sti_sympt_prob_param <-
    c(stisympt[, sum(stitest_2yr_sympt) / sum(stitest_2yr)],
    quantile(bootdist_stisympt$pct_sympt, c(0.025, 0.975)))

names(sti_sympt_prob_param) <- estlabs

sti_sympt_prob_param

# %% STI Testing by Demo - Ever PrEP Takers ----------------------------------

# STI tests sought:
#  - while on PrEP
#  - outside of a regular PrEP follow-up visit

stisympt_p <- an[
  prep_revised == 1 &
  !is.na(stitest_2yr_prep) &
  !is.na(stitest_2yr_sympt_prep),
  .(id,
    stitest_2yr = stitest_2yr_prep,
    stitest_2yr_sympt = stitest_2yr_sympt_prep
  )]

print(stisympt_p)

set.seed(1971)
bootdist_stisympt_p <- replicate(
  n = 10000,
  expr = do.call("bootsti",
                 args = list(data = stisympt_p,
                             n = nrow(stisympt_p))),
  simplify = F) %>%
  rbindlist

sti_sympt_prob_prep_param <-
  c(stisympt_p[, sum(stitest_2yr_sympt) / sum(stitest_2yr)],
    quantile(bootdist_stisympt_p$pct_sympt, c(0.025, 0.975)))

names(sti_sympt_prob_prepusers_param) <- estlabs
sti_sympt_prob_prepusers_param


# STI test frequency at PrEP follow-up visits

# .. Throat
an[prep_revised == 1, freq(prep_stithroatfreq)] %>% print

# .. Anus
an[prep_revised == 1, freq(prep_stirectfreq)] %>% print

# .. Urine or genital swab
an[prep_revised == 1, freq(prep_stiurethfreq)] %>% print


# %% Mean (12-month) Degree by Race/Ethnicity ----------------------------------

rc <- unique(an$race.cat)
p <- ggplot(an, aes(group = race.cat))


sapply(rc, function(x) an[race.cat == x, summary(pnoa_12m)])
p + geom_boxplot(aes(y = pnoa_12m))

sapply(rc, function(x) an[race.cat == x, summary(pnua_12m)])
p + geom_boxplot(aes(y = pnua_12m))

sapply(rc, function(x) an[race.cat == x, summary(pna_12m)])
p + geom_boxplot(aes(y = pna_12m))

# %% Mean Degree by Age ---------------------------------------------

plotdeg_by_age <- function(yvar) {

  ggplot(an, aes(x = age, y = !!yvar)) +
    geom_point(size = 0.7) +
    geom_smooth()

}


plotdeg_by_age(quo(pnoa_12m)) +
  labs(title = "Oral or Anal Partners (12 months)")

plotdeg_by_age(quo(pna_12m)) +
  labs(title = "Anal-only Partners (12 months)")

plotdeg_by_age(quo(pnua_12m)) +
  labs(title = "Unprotected Anal Partners (12 months)")

plotdeg_by_age(quo(pn_ongoing)) +
  labs(title = "Ongoing Partners (Oral or Anal)")
