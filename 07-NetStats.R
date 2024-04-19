################################################################################
                                  ## SETUP ##
################################################################################

pacman::p_load(
  data.table,
  EpiModelHIV,
  readxl,
  magrittr,
  ggplot2,
  ggthemes,
  rms,
  nnet
)


################################################################################
                                  ## INPUTS ##
################################################################################

network_size <- 20000

# Degree distributions
## degdist <- readRDS("netstats/aggregate_degree_summaries.Rds")
## str(degdist)

# Predictions based on ART-Net Partnership Data
pdat <- readRDS("netstats/predictions.Rds")
str(pdat)

epistats <- readRDS("netstats/epistats.Rds")


################################################################################
                             ## POPULATION SIZE ##
################################################################################

num <- network_size


################################################################################
                            ## RACE/ETHNICITY ##
################################################################################

# Source: Grey JA, Bernstein KT, Sullivan PS, Kidd SE, Gift TL, Hall EW, et al.
# Rates of Primary and Secondary Syphilis Among White and Black Non-Hispanic
# Men Who Have Sex With Men, United States, 2014. J Acquir Immune Defic Syndr.
# 2017;76: e65-e73. doi:10.1097/QAI.0000000000001508

race.lvls <- c("W", "B", "H", "O")
race.prob <- c(0.586, 0.135, 0.191, 0.088)
race.num <- round(num * race.prob)
sum(race.num) == num

race.dist <- data.table(race.lvls, race.prob, race.num)[order(race.lvls)]

print(race.dist)


################################################################################
                                ## AGE GROUP ##
################################################################################

# Source: Jones J, Grey JA, Purcell DW, Bernstein KT, Sullivan PS, Rosenberg
# ES. Estimating Prevalent Diagnoses and Rates of New Diagnoses of HIV at the
# State-Level by Age Group among Men Who Have Sex with Men in the United
# States. Open Forum Infect Dis. 2018 May 29

t.unit <- 52

age.grp.lvls <- list(
  "18-24" = list(min = 18 * t.unit, max = 25 * t.unit - 1),
  "25-34" = list(min = 25 * t.unit, max = 34 * t.unit - 1),
  "35-44" = list(min = 35 * t.unit, max = 44 * t.unit - 1),
  "45-54" = list(min = 45 * t.unit, max = 54 * t.unit - 1),
  "55+"   = list(min = 55 * t.unit, max = 64 * t.unit - 1)
)

age.grp.prob <- c(0.129, 0.242, 0.239, 0.216, 0.174)
age.grp.num <- round(num * age.grp.prob)
sum(age.grp.num) == num

age.grp.dist <- data.table(
  age.grp.lvls = names(age.grp.lvls),
  age.grp.prob,
  age.grp.num
)


################################################################################
                             ## ANAL ROLE CLASS ##
################################################################################

# Source: ArtNet
role.class.lvls <- names(pdat$demo$ai.role.pr)[3:5]
role.class.prob <- pdat$demo$ai.role.pr


################################################################################
                        ## INITIALIZE NODE ATTRIBUTES ##
################################################################################

# NOTE:
# These are used to seed the model.
# Some attributes are seeded independently of the others, but the burn-in
# period should ensure the population distributions are independent of these
# initial seeded values.

set.seed(19410524)

## Assign age
attr_age.grp <- sample(
  seq_len(length(names(age.grp.lvls))),
  size = num,
  prob = age.grp.prob,
  replace = TRUE
)

attr_age.wk <- vapply(attr_age.grp, FUN = function(x) {
  round(runif(n = 1, min = age.grp.lvls[[x]]$min, age.grp.lvls[[x]]$max))
  },
  FUN.VALUE = 23.56) %>%
  unname

attr_age.yr <- attr_age.wk / 52

## Assign race/ethnicity
attr_race <- sample(
  1:4,
  size = num,
  prob = race.dist[, race.prob],
  replace = TRUE
)

## Assign anal sex role
# Source: ARTNet
rolepreds <- lapply(seq_len(num), function(x) {
  r <- as.character(sort(unique(role.class.prob$race.i))[attr_race[x]])
  a <- round(attr_age.yr[x])
  rd <- role.class.prob[age.i == a & race.i == r]
  rd[, .(Insertive, Receptive, Versatile)]
}) %>% rbindlist

attr_role.class_c <- sapply(seq_len(num), function(x) {
  sample(role.class.lvls, 1, replace = TRUE, prob = unlist(rolepreds[x, ]))
})

# reformat for use in netsim
attr_role.class <- NULL
attr_role.class[attr_role.class_c == "Insertive"] <- 0
attr_role.class[attr_role.class_c == "Receptive"] <- 1
attr_role.class[attr_role.class_c == "Versatile"] <- 2
attr_role.class <- as.numeric(attr_role.class)

prop.table(table(attr_role.class))

## Assign anal insertativity quotients
vers <- which(attr_role.class == 2)
attr_ins.quot <- rep(NA, num)
attr_ins.quot[vers] <- runif(length(vers))

## Assign oral insertativity quotient
attr_ins.quot.oral <- runif(num)

## Assign diagnosis status
race_char <- c("black", "hispanic", "other", "white")
role.class.char <- c("Insertive", "Receptive", "Versatile")

attr_diag.status.prob <- predict(
  epistats$hiv.mod,
  newdata = data.table(
    race.cat = race_char[attr_race],
    age = attr_age.yr,
    role.class = role.class.char[attr_role.class + 1]
    ),
  type = "response"
)

attr_diag.status.prob[1:50]

set.seed(94084354)
attr_diag.status <- rbinom(num, 1, attr_diag.status.prob)
mean(attr_diag.status)

## Main degree
attr_deg.main.prob <- predict(
  pdat$main$degprob,
  newdata = data.table(
    race.cat = race_char[attr_race],
    age.grp = attr_age.grp,
    hiv.ego = attr_diag.status
  ),
  type = "probs"
)

set.seed(1800)
attr_deg.main <- vapply(1:num, FUN.VALUE = 1, FUN = function(x) {
  degprob <- as.vector(attr_deg.main.prob[x, ])
  draw <- sample(0:2, 1, prob = degprob)
  draw
})

prop.table(table(attr_deg.main))

## Casual degree
attr_deg.casl.prob <- predict(
  pdat$casl$degprob,
  newdata = data.table(
    race.cat = race_char[attr_race],
    age.grp = attr_age.grp,
    hiv.ego = attr_diag.status,
    degmain_trunc2 = attr_deg.main
  ),
  type = "probs"
)

set.seed(5882300)
attr_deg.casl <- vapply(1:num, FUN.VALUE = 1, FUN = function(x) {
  degprob <- as.vector(attr_deg.casl.prob[x, ])
  draw <- sample(0:5, 1, prob = degprob)
  draw
})

prop.table(table(attr_deg.casl))


################################################################################
                  ## AGE- AND RACE-SPECIFIC MORTALITY RATES ##
################################################################################

# Source: Arias E, Heron M, Xu J. United states life tables, 2014. Natl Vital
# Stat Rep. 2017;66:1–64. Tables 11, 14, and 17.

lt <- here::here("data", list.files("data", pattern = "lifetables"))

pull_lt <- function(data,
                    agerange = "A22:B68",
                    colselect = c("age", "qx")) {

  read_excel(
    data,
    skip = 2,
    range = agerange,
    col_names = colselect
  ) %>%
  setDT %>%
  .[, race := stringr::str_extract(data, "(?<=table[0-9]{2}-).*(?=-)")]

}

lt_long <- lapply(lt, function(x) pull_lt(data = x)) %>% rbindlist

lt_wide <- dcast(lt_long, age ~ race, value.var = "qx") %>%
  .[, agestart := stringr::str_extract(age, "^[0-9]{2}")] %>%
  .[, age.grp := dplyr::case_when(
    agestart %in% c(18:24) ~ "[18, 25)",
    agestart %in% c(25:34) ~ "[25, 35)",
    agestart %in% c(35:44) ~ "[35, 45)",
    agestart %in% c(45:54) ~ "[45, 55)",
    agestart >= 55 ~ "55+"
    )]

lt_wide[, ":=" (
  vec.asmr.B = 1 - (1 - nhblack)^(1 / 52),
  vec.asmr.H = 1 - (1 - hisp)^(1 / 52),
  vec.asmr.O = 1 - (1 - nhwhite)^(1 / 52),
  vec.asmr.W = 1 - (1 - nhwhite)^(1 / 52)
)]

lt_wide

# export weekly mortality rates
asmr <- lt_wide[, .(
  age = as.numeric(agestart),
  vec.asmr.B,
  vec.asmr.H,
  vec.asmr.O,
  vec.asmr.W
)]

asmr <- rbind(
  asmr,
  data.table(
    age = 65,
    vec.asmr.B = 1,
    vec.asmr.H = 1,
    vec.asmr.O = 1,
    vec.asmr.W = 1
  ))

## # calculate population margin mortality rate (accounting for deterministic "death" at 65)
mort_rates <- rbind(
  lt_wide[, .(age, nhblack, hisp, nhwhite, age.grp)],
  data.table(age = 65, nhblack = 1, hisp = 1, nhwhite = 1, age.grp = "55+")
)[, .(B = mean(nhblack),
      H = mean(hisp),
      W = mean(nhwhite)), by = age.grp]

# Assign Other the same mortality rate as Whites
mort_rates[, O := W]
mort_rates_annual_raspec <- as.data.frame(mort_rates)

print(mort_rates_annual_raspec)

asmr[, age.grp := cut(
         age,
         breaks = c(18, 25, 35, 45, 55, 66),
         right = FALSE,
         labels = FALSE
       )]

asmr[, paste(range(age), collapse = ","), age.grp]

mort_rates_long <- melt(
  as.data.table(mort_rates_annual_raspec),
  id.vars = "age.grp",
  measure.vars = c("B", "H", "O", "W")
)[, .(age.grp, race = variable, rate_percap = value)
  ][, ":="(
   age.grp_numeric = rep(1:5, 4),
   race.cat_numeric = fcase(
     race == "B", 1,
     race == "H", 2,
     race == "O", 3,
     race == "W", 4
   )
 )]

mort_rates_long[, ":=" (
  rprob = race.dist$race.prob[match(race, race.dist$race.lvls)],
  aprob = age.grp.dist$age.grp.prob[age.grp_numeric]
)][, jt_prob := rprob * aprob]

pop <- data.table(
  age.grp_numeric = attr_age.grp,
  race.cat_numeric = attr_race
)

pop <- pop[
  mort_rates_long[, .(age.grp_numeric, race.cat_numeric, rate_percap)],
  on = c("age.grp_numeric", "race.cat_numeric")
]

print(pop)

## Calculate population margin departure rate, weighted by distribution of
## race/ethnicity and age group in the virtual population.
mort_rate_annual_popmargin <-
  mort_rates_long[, weighted.mean(rate_percap, jt_prob)]

print(mort_rate_annual_popmargin)


################################################################################
                              ## CIRCUMCISION ##
################################################################################

library(RNHANES)
library(survey)

circ_in <- as.data.table(
  nhanes_load_data("SXQ_I", "2015-2016", demographics = TRUE)
)

circ <- circ_in[, .(
  SEQN, DMDHRGND, DMDHRAGE, RIDRETH3,
  WTINT2YR, SDMVPSU, SDMVSTRA, SXQ280
)]

names(circ) <- tolower(names(circ))

circ <- circ[dmdhrgnd == 1]

# missing values
circ[sxq280 %in% c(7, 9), sxq280 := NA]
circ[sxq280 %in% c(NA, 7, 9), missing := 1]
circ[sxq280 %in% c(1, 2), missing := 0]

# race/ethnicity recode
circ[, race4 := fcase(
         ridreth3 == 3, "White",
         ridreth3 == 4, "Black",
         ridreth3 %in% c(1, 2), "Hispanic",
         ridreth3 %in% c(6, 7), "Other"
       )]

# selection weight numerator
design <- svydesign(
  ids = ~ seqn,
  data = circ,
  weights = ~ wtint2yr
)

selwt_mod <- svyglm(
  missing ~ race4 * dmdhrage,
  design,
  family = "binomial"
)

selwt_mod_p2 <- svyglm(
  missing ~ race4 * poly(dmdhrage, 2),
  design,
  family = "binomial"
)

selwt_mod_p3 <- svyglm(
  missing ~ race4 * poly(dmdhrage, 3),
  design,
  family = "binomial"
)

selwt_mod_s3 <- svyglm(
  missing ~ race4 * rcs(dmdhrage, 3),
  design,
  family = "binomial"
)

selwt_mod_s5 <- svyglm(
  missing ~ race4 * rcs(dmdhrage, 5),
  design,
  family = "binomial"
)

AIC(selwt_mod, selwt_mod_p2, selwt_mod_p3, selwt_mod_s3, selwt_mod_s5)

selwt_pred <- predict(selwt_mod_s5, newdata = circ, type = "response")
selwt_pred <- selwt_pred[seq_len(length(selwt_pred))]

summary(selwt_pred)
boxplot(selwt_pred)

circ[, selwt := fcase(
         missing == 1, mean(missing) / selwt_pred,
         missing == 0, (1 - mean(missing)) / (1 - selwt_pred)
       )][, comb_wt := wtint2yr * selwt]

# combined weight
comb_design <- svydesign(
  ids = ~ seqn,
  data = circ[missing == 0],
  weights = ~ comb_wt
)

comb_design_num <- svydesign(
  ids = ~ seqn,
  data = circ[missing == 0],
  weights = ~ comb_wt_numcd
)

circ_mod_pl <- svyglm(
  sxq280 == 1 ~ race4 * dmdhrage,
  comb_design,
  family = "binomial"
)

circ_mod_p2 <- svyglm(
  sxq280 == 1 ~ race4 * poly(dmdhrage, 2),
  comb_design,
  family = "binomial"
)

circ_mod_p3 <- svyglm(
  sxq280 == 1 ~ race4 * poly(dmdhrage, 3),
  comb_design,
  family = "binomial"
)

circ_mod_s3 <- svyglm(
  sxq280 == 1 ~ race4 * rcs(dmdhrage, 3),
  comb_design,
  family = "binomial"
)

circ_mod_s5 <- svyglm(
  sxq280 == 1 ~ race4 * rcs(dmdhrage, 5),
  comb_design,
  family = "binomial"
)

### pick poly3 model
AIC(circ_mod_pl, circ_mod_p2, circ_mod_p3, circ_mod_s3, circ_mod_s5)

rlabs <- c("Black", "Hispanic", "Other", "White")
nd2 <- data.table(race4 = rlabs[attr_race], dmdhrage = attr_age.yr)

circmods <- ls(pattern = "circ_mod")

circpreds <- lapply(
  setNames(circmods, circmods),
  function(.x) {
    cbind(
      nd2,
      as.data.table(predict(get(.x), newdata = nd2, type = "response"))
    )
  }) %>% rbindlist(., idcol = "model")

ggplot(
  circpreds,
  aes(x = dmdhrage, y = response, color = model, fill = model)
  ) +
  geom_line() +
  ## geom_ribbon(
  ##   aes(ymin = response - SE, ymax = response + SE),
  ##   alpha = 0.2, color = "white") +
  facet_wrap(~ race4) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_base(base_size = 24)

circ_pred_compare <- dcast(
  circpreds[, .(circ_prob = mean(response)), .(model, race4)],
  race4 ~ model
)

circ.probs <- circ_pred_compare[, circ_mod_p3]
names(circ.probs) <- rlabs


################################################################################
                        ## NODEFACTOR HELPER OBJECTS ##
################################################################################

## This function calculates nodefactor terms for categorical agent attributes.
calc_nodefactor <- function(outcome, groupvar, data = attr_dt) {
  dt <- data[, .(.N, mean = mean(get(outcome))), keyby = groupvar]
  dt[, target := mean * N]
  dt[, target]
}

## This function checks that the nodefactor terms generated by calc_nodefactor()
## imply the correct number of total edges in the network.
check_nodefactor <- function(nf, edges) {
  ifelse(
    all.equal(sum(nf) / 2, edges),
    "Good to go. Target stats imply the same number of total edges.",
    "Squash a bug somewhere."
  )
}

## Create three-level deg.casl attribute for use in concurrency model.
attr_deg.casl3 <- ifelse(attr_deg.casl >= 3, 3, attr_deg.casl)
table(attr_deg.casl, attr_deg.casl3)

## Data frame to use to calculate nodefactor terms.
attr_dt <- data.table(
  attr_race,
  attr_age.grp,
  attr_age.yr,
  attr_deg.casl,
  attr_deg.casl3,
  attr_deg.main,
  attr_role.class
)


################################################################################
               ## MAIN PARTNERSHIPS, NETWORK TARGET STATISTICS ##
################################################################################

## Data frame to use to predict from nodematch models. Includes only agents
## initialized to have at least 1 main partnership.
pred_main_nodematch <- data.table(
  race.i = race_char[attr_race],
  age.grp.i = attr_age.grp,
  diag.status.i = attr_diag.status,
  ptype = 1,
  dm = attr_deg.main
)[dm > 0][, dm := NULL]

## Data frame to use to predict from nodemix models. Includes only agents
## initialized to have at least 1 main partnership.
pred_main_nodemix <- data.table(
  race.i = race_char[attr_race],
  age.grp.i = attr_age.grp,
  diag.status.i = attr_diag.status,
  ptype = 1,
  dm = attr_deg.main
)[dm > 0][, dm := NULL]

## Calculate total number of edges in main network.
edges_main <- mean(attr_deg.main) * num / 2
edges_main

## nodefactor_race_main
nodefactor_race_main <- calc_nodefactor("attr_deg.main", "attr_race")
check_nodefactor(nodefactor_race_main, edges_main)

## nodefactor_age.grp_main
nodefactor_age.grp_main <- calc_nodefactor("attr_deg.main", "attr_age.grp")
check_nodefactor(nodefactor_age.grp_main, edges_main)

## nodefactor_diag.status_main
nodefactor_diag.status_main <- calc_nodefactor(
  "attr_deg.main", "attr_diag.status"
)
check_nodefactor(nodefactor_diag.status_main, edges_main)

## nodefactor_deg.casl_main
nodefactor_deg.casl_main <- calc_nodefactor("attr_deg.main", "attr_deg.casl")
check_nodefactor(nodefactor_deg.casl_main, edges_main)

## nodematch_race.eth_main
samerace_prob_main <- predict(
  pdat$mc$racematch,
  newdata = pred_main_nodematch,
  type = "response"
)

nodematch_race.eth_main <- mean(samerace_prob_main) * edges_main
nodematch_race.eth_main

## nodematch_age.grp_main
sameage_prob_main <- predict(
  pdat$mc$age.grpmatch,
  newdata = pred_main_nodematch,
  type = "response"
)

nodematch_age.grp_main <- mean(sameage_prob_main) * edges_main
nodematch_age.grp_main

## nodematch_diag.status_main
samediag_prob_main <- 1 - predict(
  pdat$mc$hiv.discord,
  newdata = pred_main_nodematch,
  type = "response"
)

nodematch_diag.status_main <- mean(samediag_prob_main) * edges_main
nodematch_diag.status_main

## nodemix_race.eth_main
racemix_probs_main <- rowMeans(sapply(
  pdat$mc$racemix,
  function(.x) {
    colMeans(predict(.x, newdata = pred_main_nodemix, type = "probs"))
  }
))

sum(racemix_probs_main) == 1
nodemix_race.eth_main <- racemix_probs_main * edges_main
sum(nodemix_race.eth_main) == edges_main

## nodemix_age.grp_main
agemix_probs_main <- rowMeans(sapply(
  pdat$mc$agemix,
  function(.x) {
    colMeans(predict(.x, newdata = pred_main_nodemix, type = "probs"))
  }
))

sum(agemix_probs_main) == 1
nodemix_age.grp_main <- agemix_probs_main * edges_main
sum(nodemix_age.grp_main) == edges_main

## Set up prediction data for concurrency term.
pred_main_conc <- data.table(
  race.cat = race_char,
  age.grp = attr_age.grp,
  degcasl3 = attr_deg.casl3,
  hiv.ego = attr_diag.status
)

conc_prob_main <- predict(
  pdat$main$concurrent,
  newdata = pred_main_conc,
  type = "response"
)

concurrent_main <- mean(conc_prob_main) * num

## Average partnership duration among main partnerships.
pred_main_durat <- data.table(
  ego.race.cat = race_char[attr_race],
  ego.age.grp = attr_age.grp,
  hiv2 = attr_diag.status,
  dm = attr_deg.main
)[dm > 0][, dm := NULL]

durat_wks_main_preds <- predict(
  pdat$main$durat_wks,
  newdata = pred_main_durat,
  type = "response"
)

durat_wks_main <- mean(durat_wks_main_preds)
durat_wks_main

pred_main_durat[, durat_preds := durat_wks_main_preds]

## average duration by age group (ego population)
durat_wks_main_byage <-
  pred_main_durat[, .(
    mean_duration = mean(durat_preds)
  ), keyby = ego.age.grp][, mean_duration]

durat_wks_main_byage

## average duratoin by age group combo (within partnership)
pred_main_durat_byagec <- data.table(
  age_combo = c(
    "11", "12", "13", "14", "15",
    "22", "23", "24", "25",
    "33", "34", "35",
    "44", "45",
    "55"
  )
)

durat_wks_main_byagec <- unname(rowMeans(sapply(
  pdat$main$durat_wks_byagec$analyses,
  function(.x) predict(.x, newdata = pred_main_durat_byagec, type = "response")
)))

durat_wks_main_byagec


################################################################################
              ## CASUAL PARTNERSHIPS, NETWORK TARGET STATISTICS ##
################################################################################

## Data frame to use to predict from nodematch models. Includes only agents
## initialized to have at least 1 casual partnership.
pred_casl_nodematch <- data.table(
  race.i = race_char[attr_race],
  age.grp.i = attr_age.grp,
  diag.status.i = attr_diag.status,
  ptype = 2,
  dc = attr_deg.casl
)[dc > 0][, dc := NULL]

## Data frame to use to predict from nodemix models. Includes only agents
## initialized to have at least 1 main partnership.
pred_casl_nodemix <- data.table(
  race.i = race_char[attr_race],
  age.grp.i = attr_age.grp,
  diag.status.i = attr_diag.status,
  ptype = 2,
  dc = attr_deg.casl
)[dc > 0][, dc := NULL]


## Calculate total number of edges in casual network.
edges_casl <- mean(attr_deg.casl) * num / 2
edges_casl

## nodefactor_race_casl
nodefactor_race_casl <- calc_nodefactor("attr_deg.casl", "attr_race")
check_nodefactor(nodefactor_race_casl, edges_casl)

## nodefactor_age.grp_casl
nodefactor_age.grp_casl <- calc_nodefactor("attr_deg.casl", "attr_age.grp")
check_nodefactor(nodefactor_age.grp_casl, edges_casl)

## nodefactor_diag.status_casl
nodefactor_diag.status_casl <- calc_nodefactor(
  "attr_deg.casl", "attr_diag.status"
)
check_nodefactor(nodefactor_diag.status_casl, edges_casl)

## nodefactor_deg.main_casl
nodefactor_deg.main_casl <- calc_nodefactor("attr_deg.casl", "attr_deg.main")
check_nodefactor(nodefactor_deg.main_casl, edges_casl)

## nodematch_race.eth_casl
samerace_prob_casl <- predict(
  pdat$mc$racematch,
  newdata = pred_casl_nodematch,
  type = "response"
)

nodematch_race.eth_casl <- mean(samerace_prob_casl) * edges_casl
nodematch_race.eth_casl

## nodematch_age.grp_casl
sameage_prob_casl <- predict(
  pdat$mc$age.grpmatch,
  newdata = pred_casl_nodematch,
  type = "response"
)

nodematch_age.grp_casl <- mean(sameage_prob_casl) * edges_casl
nodematch_age.grp_casl

## nodematch_diag.status_casl
samediag_prob_casl <- 1 - predict(
  pdat$mc$hiv.discord,
  newdata = pred_casl_nodematch,
  type = "response"
)

nodematch_diag.status_casl <- mean(samediag_prob_casl) * edges_casl
nodematch_diag.status_casl

## nodemix_race.eth_casl
racemix_probs_casl <- rowMeans(sapply(
  pdat$mc$racemix,
  function(.x) {
    colMeans(predict(.x, newdata = pred_casl_nodemix, type = "probs"))
  }
))

sum(racemix_probs_casl) == 1
nodemix_race.eth_casl <- racemix_probs_casl * edges_casl

## nodemix_age.grp_casl
agemix_probs_casl <- rowMeans(sapply(
  pdat$mc$agemix,
  function(.x) {
    colMeans(predict(.x, newdata = pred_casl_nodemix, type = "probs"))
  }
))

sum(agemix_probs_casl) == 1
nodemix_age.grp_casl <- agemix_probs_casl * edges_casl


## Set up prediction data for concurrency term.
pred_casl_conc <- data.table(
  race.cat = race_char,
  age.grp = attr_age.grp,
  degmain_trunc2 = attr_deg.main,
  hiv.ego = attr_diag.status
)

conc_prob_casl <- predict(
  pdat$casl$concurrent,
  newdata = pred_casl_conc,
  type = "response"
)

concurrent_casl <- mean(conc_prob_casl) * num

## Average partnership duration among casual partnerships.
pred_casl_durat <- data.table(
  ego.race.cat = race_char[attr_race],
  ego.age.grp = attr_age.grp,
  hiv2 = attr_diag.status,
  dc = attr_deg.casl
)[dc > 0][, dc := NULL]

durat_wks_casl_preds <- predict(
  pdat$casl$durat_wks,
  newdata = pred_casl_durat,
  type = "response"
)

durat_wks_casl <- mean(durat_wks_casl_preds)
durat_wks_casl

pred_casl_durat[, durat_preds := durat_wks_casl_preds]

## casual duration by age group (ego population)
durat_wks_casl_byage <-
  pred_casl_durat[, .(
    mean_duration = mean(durat_preds)
  ), keyby = ego.age.grp][, mean_duration]

durat_wks_casl_byage

## casual duration by age group combo (per partnership)
pred_casl_durat_byagec <- data.table(
  age_combo = c(
    "11", "12", "13", "14", "15",
    "22", "23", "24", "25", "33",
    "34", "35",
    "44", "45",
    "55"
  )
)

durat_wks_casl_byagec <- unname(rowMeans(sapply(
  pdat$casl$durat_wks_byagec$analyses,
  function(.x) predict(.x, newdata = pred_casl_durat_byagec, type = "response")
)))

durat_wks_casl_byagec


################################################################################
## ONE-TIME CONTACTS, NETWORK TARGET STATISTICS ##
################################################################################

## Set up prediction data set for instantaneous partnerships.
pred_inst_nodefact <- data.table(
  race.cat = race_char[attr_race],
  age.grp = attr_age.grp,
  degcasl = attr_deg.casl,
  degmain_trunc2 = attr_deg.main,
  hiv.ego = attr_diag.status
)

## Calculate total number of instantaneous partnerships. We divide by
## 52 here because the instantaneous partnership rates were originally
## estimated on the yearly time scale. So here, we calculate the
## expected number of instantaneous partnerships per week.
attr_deg.inst <- predict(
  pdat$inst$instrate,
  pred_inst_nodefact,
  type = "response"
) / 52

edges_inst <- mean(attr_deg.inst) * num / 2
edges_inst

## nodefactor_race
nodefactor_race_i <- calc_nodefactor("attr_deg.inst", "attr_race")
check_nodefactor(nodefactor_race_i, edges_inst)

## nodefactor_age.grp
nodefactor_age.grp_i <- calc_nodefactor("attr_deg.inst", "attr_age.grp")
check_nodefactor(nodefactor_age.grp_i, edges_inst)

## nodefactor_degmain
nodefactor_deg.main_i <- calc_nodefactor("attr_deg.inst", "attr_deg.main")
check_nodefactor(nodefactor_deg.main_i, edges_inst)

## nodefactor_degcasl
nodefactor_deg.casl_i <- calc_nodefactor("attr_deg.inst", "attr_deg.casl")
check_nodefactor(nodefactor_deg.casl_i, edges_inst)

## nodefactor_diagstatus
nodefactor_diag.status_i <- calc_nodefactor("attr_deg.inst", "attr_diag.status")
check_nodefactor(nodefactor_diag.status_i, edges_inst)


################################################################################
## SAVE PARAMETERS TO FILE ##
################################################################################

out <- list()

# ... STORE INPUT PARAMETERS

out$inputs <- list()
out$inputs$race.dist <- race.dist[, -c("race.num")]
out$inputs$age.grp.dist <- age.grp.dist[, -c("age.grp.num")]
out$inputs$role.class.dist <- prop.table(table(attr_role.class_c))
out$inputs$circ.probs <- circ.probs

out$inputs

# ... STORE DEMOGRAPHICS

out$demog <- list()
out$demog$num <- num
out$demog$asmr <- asmr
out$demog$age.breaks <- c(18, 25, 35, 45, 55, 66)

# Race/ethnicity
for (i in seq_len(length(race.lvls))) {
  out$demog[paste0(
        "num.",
        race.dist[i, sort(race.lvls)]
      )] <- race.dist[i, race.num]
}

# Age.Grp
for (i in seq_len(length(age.grp.lvls))) {
  out$demog[paste0("num.age.grp_", i)] <- age.grp.dist[i, age.grp.num]
}

# Age limits
out$demog$ages <- c(18, 65)

# Mortality Rates (Weekly, Age-and-Race-specific)
for (i in seq_len(length(mort_rates_annual_raspec) - 1)) {
  for (j in 1:nrow(mort_rates_annual_raspec)) {
    out$demog[paste0("mortrate.",
    "race.", names(mort_rates_annual_raspec)[i + 1],
    ".age.grp_", j)] <-
    round(1 - (1 - mort_rates_annual_raspec[j, i + 1]) ^ (1 / 52), 6)
  }
}

# Mortality/Departure Rate (weekly, population margin)
out$demog$mortrate.marginal <- round(
  1 - (1 - mort_rate_annual_popmargin) ^ (1 / 52), digits = 6
)

# Anal Role Class
for (i in seq_len(length(role.class.lvls))) {
  out$demog[paste0(
        "role.class.", role.class.lvls[i]
      )] <- sum(attr_role.class_c == role.class.lvls[i])
}

print(out$demog)

# ... NETWORK MODEL TARGETS (MAIN PARTNERSHIPS)
out$netmain <- list()
out$netmain$edges <- edges_main
out$netmain$nodefactor_race <- nodefactor_race_main
out$netmain$nodefactor_age.grp <- nodefactor_age.grp_main
out$netmain$nodefactor_degcasl <- nodefactor_deg.casl_main
out$netmain$nodefactor_diagstatus <- nodefactor_diag.status_main
out$netmain$nodematch_race <- nodematch_race.eth_main
out$netmain$nodematch_age.grp <- nodematch_age.grp_main
out$netmain$nodematch_diagstatus <- nodematch_diag.status_main
out$netmain$nodemix_race <- nodemix_race.eth_main
out$netmain$nodemix_age.grp <- nodemix_age.grp_main
out$netmain$concurrent <- concurrent_main
out$netmain$durat_wks <- durat_wks_main
out$netmain$durat_wks_byage <- durat_wks_main_byage
out$netmain$durat_wks_byagec <- durat_wks_main_byagec

# ... NETWORK MODEL TARGETS (CASUAL PARTNERSHIPS)
out$netcasl <- list()
out$netcasl$edges <- edges_casl
out$netcasl$nodefactor_race <- nodefactor_race_casl
out$netcasl$nodefactor_age.grp <- nodefactor_age.grp_casl
out$netcasl$nodefactor_degmain <- nodefactor_deg.main_casl
out$netcasl$nodefactor_diagstatus <- nodefactor_diag.status_casl
out$netcasl$nodematch_race <- nodematch_race.eth_casl
out$netcasl$nodematch_age.grp <- nodematch_age.grp_casl
out$netcasl$nodematch_diagstatus <- nodematch_diag.status_casl
out$netcasl$nodemix_race <- nodemix_race.eth_casl
out$netcasl$nodemix_age.grp <- nodemix_age.grp_casl
out$netcasl$concurrent <- concurrent_casl
out$netcasl$durat_wks <- durat_wks_casl
out$netcasl$durat_wks_byage <- durat_wks_casl_byage
out$netcasl$durat_wks_byagec <- durat_wks_casl_byagec

# ... NETWORK MODEL TARGETS (INSTANTANEOUS PARTNERSHIPS)
out$netinst <- list()
out$netinst$edges <- edges_inst
out$netinst$nodefactor_race <- nodefactor_race_i
out$netinst$nodefactor_age.grp <- nodefactor_age.grp_i
out$netinst$nodefactor_degmain <- nodefactor_deg.main_i
out$netinst$nodefactor_degcasl <- nodefactor_deg.casl_i
out$netinst$nodefactor_diagstatus <- nodefactor_diag.status_i

# ... STORE ATTRIBUTES

out$attr <- list()
out$attr$age <- attr_age.yr
out$attr$sqrt.age <- sqrt(attr_age.yr)
out$attr$age.wk <- attr_age.wk
out$attr$age.grp <- attr_age.grp
out$attr$race <- attr_race
out$attr$deg.main <- attr_deg.main
out$attr$deg.casl <- attr_deg.casl
out$attr$role.class <- attr_role.class
out$attr$ins.quot <- attr_ins.quot
out$attr$ins.quot.oral <- attr_ins.quot.oral
out$attr$diag.status <- attr_diag.status


# ... WRITE
saveRDS(out, here::here("netstats", "netstats.Rds"))
