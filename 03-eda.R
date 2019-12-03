# %% PrEP Variables ----------------------------------------------------------

prepvars <- varselect[grepl("prep", varselect)]

prepl <- lapply(set_names(prepvars, prepvars), function(x) {
  avs[artnetprep_current == 1 | prep_revised == 1, .N,
        keyby = x][, p := N / sum(N)]
  })

prepl

# %% Descriptive Statistics and Derived Variables ----------------------------

# Age
avs[, descr(age)] %>% print


avs[, freq(artnetevertest)] %>% print


# Create and view race/ethnicity variable
avs[, .N, by = .(race, race.cat)] %>% print
raceth_lookup <- c("white",
                   "other",
                   "black",
                   "other",
                   "other",
                   "other",
                   "other")

names(raceth_lookup) <- unique(avs$race)

avs[, raceth := ifelse(race.cat == "hispanic",
                         "hispanic",
                         raceth_lookup[race]
                         )]

# .. check race/ethnicity coding
avs[, .N, .(race, race.cat, raceth)] %>%
  .[, pct := N / sum(N) * 100] %>%
  .[order(race, race.cat, raceth)] %>%
  print

freq(avs$raceth) %>% print


# %% HIV by Demo --------------------------------------------------------------

freq(avs$hiv)
avs[, stby(hiv, raceth, freq)]

# HIV status by race, age
library(rms)

avhiv <- avs[!is.na(hiv), .(id, age, raceth, hiv)]
print(avhiv)
nrow(avhiv)

bootit <- function(n) avhiv[sample(1:n, replace = T)]

hiv_bsamps <- replicate(n = 10000, simplify = F,
                    expr = do.call("bootit",
                                   args = list(n = nrow(avhiv))))

print(hiv_bsamps[1:5])

testix <- function(bsamps, dep, ind, dat, fam = "binomial",
                   keepreg = FALSE) {

  fml1 <- as.formula(paste(dep, "~", paste(ind, collapse = "+")))
  fml2 <- as.formula(paste(dep, "~", paste(ind, collapse = "*")))

  fit1 <- glm(fml1, data = dat, family = fam)
  fit2 <- glm(fml2, data = dat, family = fam)

  lrt <- lrtest(fit1, fit2)

  if (keepreg) {
    list(fit1, fit2, lrt)
  } else {
    list(lrt)
  }
}

hiv_lrt <- lapply(hiv_bsamps, function(x) {
  testix(dat = x, dep = "hiv", ind = c("raceth", "age"))
})

hiv_lrt_p <- lapply(hiv_lrt, function(x) x[[1]]$stats[["P"]]) %>% unlist
print(descr(hiv_lrt_p))
summary(hiv_lrt_p)
hiv_aic <- lapply(hiv_bsamps, function(x) {
  fit <- glm("hiv ~ raceth * age", data = x, family = "binomial")
  MASS::stepAIC(fit)
  })

hiv_aic_varsets <- lapply(hiv_aic, function(x) coef(x))


predhra <- expand.grid(raceth = unique(avs$raceth),
                       age = unique(avs$age)) %>%
                       setDT %>%
                       .[order(raceth, age)]


hiv_preds <- lapply(hiv_bsamps)
predhra$hiv_hat <- predict(object = hra,
                           newdata = predhra,
                           type = "response")

print(predhra)

ggplot(predhra) +
  geom_line(
    aes(x = age, y = hiv_hat, color = raceth),
    size = 1
    ) +
  scale_color_viridis_d() +
  theme_clean()

# %% PrEP by Demo -------------------------------------------------------------

freq(avs$prep_revised) %>% print


avs[prep_revised == 1, freq(artnetprep_current)] %>%
  print

avs[, an_prep_current := ifelse(
  !is.na(prep_revised) & is.na(artnetprep_current),
  0, artnetprep_current
  )]

freq(avs$an_prep_current)

# %% Simple Summaries by Race/ethnicity ---------------------------------------

print(sort(names(avs)))

# .. continuous variables
avs[, stby(age, raceth, descr)] %>%
  print

avs[, stby(zap_labels(stitestfreq), raceth, descr)] %>%
  print

# Source for raincloud plots:

# Allen M, Poggiali D, Whitaker K, Marshall TR, Kievit RA. Raincloud plots: a
# multi-platform tool for robust data visualization. Wellcome Open Res
# 2019;4:63. doi:10.12688/wellcomeopenres.15191.1.

# @DEV 2019-10-15: Make a raincloud plot function

make_it_rain <- function(data, xvar, yvar, fillvar, colvar,
                         groupvar = NULL, legend_title) {
   avs %>%
     ggplot(aes(x = !!xvar,
                y = !!yvar,
                fill = !!fillvar,
                color = !!colvar,
                group = !!groupvar)) +
     geom_density_ridges(scale = 0.7, alpha = 0.2) +
     geom_point(position = position_jitter(width = 0.4, height = 0.1),
                size = 0.25) +
     scale_color_viridis_d() +
     scale_fill_viridis_d(name = legend_title) +
     guides(color = F) +
     theme_ridges() +
     theme(panel.grid.major.y = element_blank(),
           axis.ticks.y = element_blank())
}

make_it_rain(
  data = avs,
  xvar = quo(age),
  yvar = quo(raceth),
  fillvar = quo(raceth),
  colvar = quo(raceth),
  groupvar = quo(raceth),
  legend_title = "Race/ethnicity"
)

make_it_rain(
  data = avs,
  xvar = quo(zap_labels(stitestfreq)),
  yvar = quo(raceth),
  fillvar = quo(raceth),
  colvar = quo(raceth),
  groupvar = quo(raceth),
  legend_title = "Race/ethnicity"
)


# categorical variables
avs[, stby(stireg, raceth, freq)] %>% print
avs[, stby(artnetevertest, raceth, freq)]
avs[, stby(stitestfreq, raceth, freq)]
avs[, stby(stitest_2yr, raceth, freq)]
avs[, stby(, raceth, freq)]
avs[, stby(, raceth, freq)]


by_race <- c(varselect, "raceth")
avs[, ..by_race]

stby(artnetevertest)
