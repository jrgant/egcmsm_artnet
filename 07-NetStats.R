# %% SETUP ---------------------------------------------------------------------

pacman::p_load(
  data.table,
  EpiModelHIV,
  readxl,
  magrittr,
  ggplot2,
  ggthemes
)


# %% INPUTS --------------------------------------------------------------------

network_size <- 20000

# Degree distributions
degdist <- readRDS("netstats/aggregate_degree_summaries.Rds")
str(degdist)

# Predictions based on ART-Net Partnership Data
pdat <- readRDS("netstats/predictions.Rds")
str(pdat)

epistats <- readRDS("netstats/epistats.Rds")


# %% INITIALIZE DEMOGRAPHICS ---------------------------------------------------

## population size
num <- network_size

# ... RACE and ETHNICITY

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

# ... AGE GROUP

# Source: Jones J, Grey JA, Purcell DW, Bernstein KT, Sullivan PS, Rosenberg
# ES. Estimating Prevalent Diagnoses and Rates of New Diagnoses of HIV at the
# State-Level by Age Group among Men Who Have Sex with Men in the United
# States. Open Forum Infect Dis. 2018 May 29

age.grp.lvls <- list(
  "18-24" = list(min = 18, max = 24),
  "25-34" = list(min = 25, max = 34),
  "35-44" = list(min = 35, max = 44),
  "45-54" = list(min = 45, max = 54),
  "55+"   = list(min = 55, max = 64)
)

age.grp.prob <- c(0.129, 0.242, 0.239, 0.216, 0.174)
age.grp.num <- round(num * age.grp.prob)
sum(age.grp.num) == num

age.grp.dist <- data.table(
  age.grp.lvls = names(age.grp.lvls),
  age.grp.prob,
  age.grp.num
)

# ... JOINT RACExAGE

raceage.dist <- cbind(
  expand.grid(
    race = race.dist$race.lvls,
    age.grp = age.grp.dist$age.grp.lvls
  ),
  expand.grid(
    race.prob = race.dist$race.prob,
    age.grp.prob = age.grp.dist$age.grp.prob
  )) %>%
  setDT %>%
  .[, raceage.prob := race.prob * age.grp.prob]

print(raceage.dist)
sum(raceage.dist$raceage.prob)


# ... MAIN DEGREE DISTRIBUTION

maindeg.dist <- pdat$main$degprob[
  order(race.cat, age.grp),
  .(maindeg.prob = weighted.mean(
      preds,
      raceage.dist[order(race, age.grp), raceage.prob])
    ),
  outcome]

cbind(degdist$main_artnet_sum$sum_main_trunc, maindeg.dist)
print(maindeg.dist)


# ... CASUAL DEGREE DISTRIBUTION

casldeg.dist <- pdat$casl$degprob[
  order(race.cat, age.grp),
  .(casldeg.prob = weighted.mean(
      preds,
      raceage.dist[order(race, age.grp), raceage.prob])
    ),
  outcome]

cbind(degdist$casl_artnet_sum$sum_casl_total, casldeg.dist)
print(casldeg.dist)


# ... JOINT DEMOGRAPHIC PROB (INCORPORATING MAIN/CASUAL DEGREE)

ram_grid <- cbind(
  expand.grid(
    race.cat = race.dist$race.lvls,
    age.grp.cat = age.grp.dist$age.grp.lvls,
    degmain.cat = maindeg.dist$outcome
  ),
  expand.grid(
    race.prob = race.dist$race.prob,
    age.grp.prob = age.grp.dist$age.grp.prob,
    degmain.prob = maindeg.dist$maindeg.prob
  )) %>%
  setDT  %>%
  .[, jt_prob := race.prob * age.grp.prob * degmain.prob]

print(ram_grid)

rac_grid <- cbind(
  expand.grid(
    race.cat = race.dist$race.lvls,
    age.grp.cat = age.grp.dist$age.grp.lvls,
    degcasl.cat = casldeg.dist$outcome
  ),
  expand.grid(
    race.prob = race.dist$race.prob,
    age.grp.prob = age.grp.dist$age.grp.prob,
    degcasl.prob = casldeg.dist$casldeg.prob
  )
  ) %>%
  setDT  %>%
  .[, jt_prob := race.prob * age.grp.prob * degcasl.prob]

print(rac_grid)

racm_grid <- as.data.table(cbind(
  expand.grid(
    race.cat = race.dist$race.lvls,
    age.grp.cat = age.grp.dist$age.grp.lvls,
    degcasl.cat = casldeg.dist$outcome,
    degmain.cat = maindeg.dist$outcome
    ),
  expand.grid(
    race.prob = race.dist$race.prob,
    age.grp.prob = age.grp.dist$age.grp.prob,
    degcasl.prob = casldeg.dist$casldeg.prob,
    degmain.prob = maindeg.dist$maindeg.prob
    )
  ))[, jt_prob := race.prob * age.grp.prob * degcasl.prob * degmain.prob]

print(racm_grid)


# ... ANAL SEX ROLE

# Source: ArtNet
role.class.lvls <- pdat$demo$ai.role.pr$anal.sex.role
role.class.prob <- pdat$demo$ai.role.pr$P
role.class.num <- round(num * role.class.prob)
sum(role.class.num) == num

role.class.dist <- data.table(
  role.class.lvls,
  role.class.prob,
  role.class.num
)[order(role.class.lvls)]


# ... AGE-AND-RACE/ETHNICITY-SPECIFIC MORTALITY

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
  vec.asmr.B = 1 - (1 - nhblack)^(1/52),
  vec.asmr.H = 1 - (1 - hisp)^(1/52),
  vec.asmr.O = 1 - (1 - nhwhite)^(1/52),
  vec.asmr.W = 1 - (1 - nhwhite)^(1/52)
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

# calculate population margin mortality rate
mort_rates <- lt_wide[, .(B = mean(nhblack),
                          H = mean(hisp),
                          W = mean(nhwhite)),
                          by = age.grp]

# Assign Other the same mortality rate as Whites
mort_rates[, O := W]
mort_rates_annual_raspec <- as.data.frame(mort_rates)

print(mort_rates_annual_raspec)

mort_rates_long <- melt(
  as.data.table(mort_rates_annual_raspec),
  id.vars = "age.grp",
  measure.vars = c("B", "H", "O", "W")
  ) %>%
  .[, .(age.grp, race = variable, rate_percap = value)] %>%
  .[, age.grp_num := rep(1:5, 4)]

mort_rates_long <- cbind(
  mort_rates_long,
  raceage.dist[order(race, age.grp), .(age.grp, race, raceage.prob)]
)

print(mort_rates_long)
sum(mort_rates_long$raceage.prob) == 1

mort_rates

mort_rate_annual_popmargin <-
  mort_rates_long[, weighted.mean(rate_percap, raceage.prob)]

print(mort_rate_annual_popmargin)


# %% INITIALIZE NODE ATTRIBUTES ------------------------------------------------

# NOTE:
# These are used to seed the model.
# Some attributes are seeded independently of the others, but the burn-in
# period should ensure the population distributions are independent of these
# initial seeded values.

set.seed(19410524)

## Assign age
attr_age.grp <- sample(
  1:length(names(age.grp.lvls)),
  size = num,
  prob = age.grp.prob,
  replace = T
)

attr_age.yr <- vapply(attr_age.grp, FUN = function(x) {
  runif(n = 1, min = age.grp.lvls[[x]]$min, age.grp.lvls[[x]]$max + 0.99)
  },
  FUN.VALUE = 23.56) %>%
  unname

## Assign race/ethnicity
attr_race <- sample(
  1:4,
  size = num,
  prob = race.dist[, race.prob],
  replace = TRUE
)

## Main degree
attr_deg.main <- sample(
  0:2,
  size = num,
  replace = T,
  prob = maindeg.dist[, maindeg.prob]
) %>% as.integer

## Casual degree
attr_deg.casl <- sample(
  0:5,
  size = num,
  replace = T,
  prob = casldeg.dist[, casldeg.prob]
) %>% as.integer

## Anal Sex Role
# Source: ARTNet
attr_role.class <- sample(
  role.class.dist[, role.class.lvls],
  size = num,
  prob = role.class.dist[, role.class.prob],
  replace = TRUE
)

# reformat for use in netsim
attr_role.class[attr_role.class == "I"] <- 0
attr_role.class[attr_role.class == "R"] <- 1
attr_role.class[attr_role.class == "V"] <- 2
attr_role.class <- as.numeric(attr_role.class)

prop.table(table(attr_role.class))

## HIV status
hiv.df <- data.frame(race.cat = attr_race, age = attr_age.yr)

attr_diag.status.prob <- predict(
  epistats$hiv.mod,
  newdata = hiv.df,
  type = "response"
)

attr_diag.status.prob[1:50]

attr_diag.status <- rbinom(num, 1, attr_diag.status.prob)


# %% PARTNERSHIP REUSABLES -----------------------------------------------------

race_char <- c("black", "hispanic", "other", "white")

race.wt <- raceage.dist[,
  weight := raceage.prob / sum(raceage.prob),
  keyby = race
][, .(race, age.grp, weight)]

print(race.wt)

age.grp.wt <- raceage.dist[,
  weight := raceage.prob / sum(raceage.prob),
  keyby = age.grp
][, .(race, age.grp, weight)]

print(age.grp.wt)


# %% MAIN PARTNERSHIPS ---------------------------------------------------------

## mdeg_preds <- pdat$main$degpred_joint[hiv.ego == 0]
## setkeyv(mdeg_preds, c("race.cat", "age.grp", "degcasl"))
## mdeg_preds[, degcasl := as.factor(degcasl)]

## setkeyv(rac_grid, c("race.cat", "age.grp.cat", "degcasl.cat"))
## setnames(
##   rac_grid[, age.grp.cat := as.numeric(age.grp.cat)], "age.grp.cat", "age.grp"
## )
## setnames(rac_grid, "degcasl.cat", "degcasl")
## rac_grid[, degcasl := rep(1:5, 24)]

## mlook <- mdeg_preds[rac_grid, on = c("race.cat", "age.grp", "degcasl")]

## main_lookup <- as.data.table(cbind(
##   pdat$main$degpred_joint[hiv.ego == 0][order(race.cat, age.grp, degcasl)],
##   rac_grid[order(race.cat, age.grp.cat, degcasl.cat)]
## ))

## main_lookup[, .N, .(race.cat, race.cat)]
## main_lookup[, .N, .(age.grp, age.grp.cat)]

## sum(main_lookup$jt_prob)
## print(main_lookup$jt_prob)

## calc mean degree by race
pred_race_newdata <- data.table(
  race.cat = race_char[attr_race],
  age.grp = attr_age.grp,
  degcasl = attr_deg.casl,
  hiv.ego = attr_diag.status
)


predict(pdat$main$degpred_joint, newdata = pred_race_newdata)

race_preds_m <-  mlook %>%
  .[, .(pred = weighted.mean(pred_mdeg_joint, jt_prob)), race.cat]

print(race_preds_m)

## calc mean degree by age grp (marginalized over race/ethnicity & casl.deg)
age.grp_preds_m <- main_lookup %>%
  .[, .(pred = weighted.mean(pred_mdeg_joint, jt_prob)), age.grp]

print(age.grp_preds_m)

## calculate mean degree by casual degree
casl_preds_m <- main_lookup %>%
  .[, .(pred = weighted.mean(pred_mdeg_joint, jt_prob)), degcasl.cat]

print(casl_preds_m)

## edges_main
edges_main <- main_lookup %>%
  .[, round(weighted.mean(pred_mdeg_joint, jt_prob) * num / 2)]

print(edges_main)

## nodefactor_race_main
nodefactor_race_main <- round(race_preds_m$pred * race.dist$race.prob * num)

print(sum(nodefactor_race_main) / 2)

## nodefactor_age.grp
nodefactor_age.grp_main <-
  round(age.grp_preds_m$pred * age.grp.dist$age.grp.prob * num)

print(sum(nodefactor_age.grp_main) / 2)

## nodefactor_degcasl
nodefacter_degcasl_main <-
  round(casl_preds_m$pred * casldeg.dist$casldeg.prob * num)

print(sum(nodefacter_degcasl_main) / 2)

## nodematch_race.eth
nodematch_race.eth_main <-
  round(pdat$main$racematch[samerace == 1, P] * edges_main)

print(nodematch_race.eth_main)

## nodematch_age.grp
nodematch_age.grp_main <-
  round(pdat$main$age.grpmatch[sameage == 1, P] * edges_main)

print(nodematch_age.grp_main)


## concurrent
main_concpr <- pdat$main$concurrent[order(race.cat, age.grp, degcasl)] %>%
  cbind(., main_lookup[, .(race.cat, age.grp.cat, degcasl.cat, jt_prob)])

print(main_concpr)

concurrent_main <-
  round(main_concpr[, weighted.mean(pred_main_concurrent, jt_prob)] * num)

print(concurrent_main)

attr_race.cat <- rep(NA, length(attr_race))
attr_race.cat[which(attr_race == 1)] <- "black"
attr_race.cat[which(attr_race == 2)] <- "hispanic"
attr_race.cat[which(attr_race == 3)] <- "other"
attr_race.cat[which(attr_race == 4)] <- "white"

conc_newdat <- data.table(
  race.cat = attr_race.cat,
  age.grp = as.factor(attr_age.grp),
  hiv.ego = attr_diag.status,
  degcasl3 = ifelse(attr_deg.casl >= 3, 3, attr_deg.casl)
)

pred_mconcur_prob <- predict(
  fit_main_concurrent,
  newdata = conc_newdat,
  type = "response"
)

concurrent_main <- round(mean(pred_mconcur_prob) * num)


## duration
durat_wks_main <- as.data.table(pdat$main$durat_wks)

durat_wks_main[, ego.race.cat2 := ego.race.cat]
durat_wks_main[, ego.race.cat := as.numeric(as.factor(ego.race.cat2))]
durat_wks_main[, .N, .(ego.race.cat, ego.race.cat2)]

durat_pred_df <- data.table(
  ego.race.cat = attr_race,
  ego.age.grp = attr_age.grp,
  ego.hiv = attr_diag.status
)

durat_main_join <- durat_pred_df[
  durat_wks_main,
  on = c("ego.race.cat", "ego.age.grp", "ego.hiv")
  ][, -c("ego.race.cat2")]

durat_wks_main <-
  durat_main_join[, .(durat_wks_main = mean(pred_main_durat_wks))]


# %% CASUAL PARTNERSHIPS -------------------------------------------------------

casl_lookup <- as.data.table(cbind(
  pdat$casl$degpred_joint[order(race.cat, age.grp, degmain_trunc2)],
  ram_grid[order(race.cat, age.grp.cat, degmain.cat)]
))

# calculate mean degree by race (marginalized over age.grp and main degree)
race_preds_c <- casl_lookup %>%
  .[, .(pred = weighted.mean(pred_cdeg_joint, jt_prob)), race.cat]

print(race_preds_c)

# calculate mean degree by age group (marginalized over race and main degree)
age.grp_preds_c <- casl_lookup %>%
  .[, .(pred = weighted.mean(pred_cdeg_joint, jt_prob)), age.grp]

print(age.grp_preds_c)

# calculate mean degree by casual degree
main_preds_c <- casl_lookup %>%
  .[, .(pred = weighted.mean(pred_cdeg_joint, jt_prob)), degmain_trunc2]

print(main_preds_c)

## edges_casl
edges_casl <- casl_lookup %>%
  .[, round(weighted.mean(pred_cdeg_joint, jt_prob) * num / 2)]

print(edges_casl)

## nodefactor_race
nodefactor_race_casl <-
  round(race_preds_c$pred * race.dist$race.prob * num)

print(sum(nodefactor_race_casl) / 2)

## nodefactor_age.grp
nodefactor_age.grp_casl <-
  round(age.grp_preds_c$pred * age.grp.dist$age.grp.prob * num)

print(sum(nodefactor_age.grp_casl) / 2)

## nodefactor_degmain
nodefacter_degmain_casl <-
  round(main_preds_c$pred * maindeg.dist$maindeg.prob * num)

print(sum(nodefacter_degmain_casl) / 2)

## nodematch_race.eth
nodematch_race.eth_casl <-
  round(pdat$casl$racematch[samerace == 1, P] * edges_casl)

print(nodematch_race.eth_casl)

## nodematch_age.grp
nodematch_age.grp_casl <-
  round(pdat$casl$age.grpmatch[sameage == 1, P] * edges_casl)

print(nodematch_age.grp_casl)

## concurrent
casl_concpr <-
  pdat$casl$concurrent[order(race.cat, age.grp, degmain_trunc2)] %>%
  cbind(., casl_lookup[, .(race.cat, age.grp.cat, degmain.cat, jt_prob)])

print(casl_concpr)

concurrent_casl <-
  round(casl_concpr[, weighted.mean(pred_casl_concurrent, jt_prob)] * num)

print(concurrent_casl)


## duration
durat_wks_casl <- as.data.table(pdat$casl$durat_wks)

durat_wks_casl[, ego.race.cat2 := ego.race.cat]
durat_wks_casl[, ego.race.cat := as.numeric(as.factor(ego.race.cat2))]
durat_wks_casl[, .N, .(ego.race.cat, ego.race.cat2)]

durat_casl_join <- durat_pred_df[
  durat_wks_casl,
  on = c("ego.race.cat", "ego.age.grp", "ego.hiv")
][, -c("ego.race.cat2")]

durat_wks_casl <-
  durat_casl_join[, .(durat_wks_casl = mean(pred_casl_durat_wks))]


# %% ONE-TIME PARTNERSHIPS -----------------------------------------------------

inst_lookup <- cbind(
  pdat$inst$inst_joint[order(race.cat, age.grp, degcasl, degmain_trunc2)],
  racm_grid[order(race.cat, age.grp.cat, degcasl.cat, degmain.cat)]
) %>% setDT

inst_lookup[, .N, .(race.cat, race.cat)]
inst_lookup[, .N, .(age.grp, age.grp.cat)]
inst_lookup[, .N, .(degcasl, degcasl.cat)]
inst_lookup[, .N, .(degmain_trunc2, degmain.cat)]

## instantaneous partnerships by race
race_preds_i <- inst_lookup %>%
  .[, .(pred = weighted.mean(pred_inst_joint, jt_prob)), race.cat]

print(race_preds_i)

## instantaneous partnerships by age.grp
age.grp_preds_i <- inst_lookup %>%
  .[, .(pred = weighted.mean(pred_inst_joint, jt_prob)), age.grp]

print(age.grp_preds_i)

## instantaneous partnerships by main degree
main_preds_i <- inst_lookup %>%
  .[, .(pred = weighted.mean(pred_inst_joint, jt_prob)), degmain_trunc2]

print(main_preds_i)

## instantaneous partnerships by casual degree
casl_preds_i <- inst_lookup %>%
  .[, .(pred = weighted.mean(pred_inst_joint, jt_prob)), degcasl]

print(casl_preds_i)

## edges_casl
edges_inst <- round(
  inst_lookup[, weighted.mean(pred_inst_joint, jt_prob)] * num / 2
)

print(edges_inst)

## nodefactor_race
nodefactor_race_i <-
  round(race_preds_i$pred * race.dist$race.prob * num)

print(sum(nodefactor_race_i) / 2)

## nodefactor_age.grp
nodefactor_age.grp_i <-
  round(age.grp_preds_i$pred * age.grp.dist$age.grp.prob * num)

print(sum(nodefactor_age.grp_i) / 2)

## nodefactor_degmain
nodefactor_degmain_i <-
  round(main_preds_i$pred * maindeg.dist$maindeg.prob * num)

print(sum(nodefactor_degmain_i) / 2)

## nodefactor_degcasl
nodefactor_degcasl_i <-
  round(casl_preds_i$pred * casldeg.dist$casldeg.prob * num)

print(sum(nodefactor_degcasl_i) / 2)


# %% SAVE PARAMETERS TO FILE ---------------------------------------------------

out <- list()

# ... STORE INPUT PARAMETERS

out$inputs <- list()
out$inputs$race.dist <- race.dist[, -c("race.num")]
out$inputs$age.grp.dist <- age.grp.dist[, -c("age.grp.num")]
out$inputs$main.dist <- maindeg.dist[, .(degmain_trunc2 = outcome,
                                         prob = maindeg.prob)]
out$inputs$casl.dist <- casldeg.dist[, .(degcasl = outcome,
                                         prob = casldeg.prob)]
out$inputs$role.class.dist <- role.class.dist[, -c("role.class.num")]

out$inputs

# ... STORE DEMOGRAPHICS

out$demog <- list()
out$demog$num <- num
out$demog$asmr <- asmr

# Race/ethnicity
for (i in 1:length(race.lvls)) {
  out$demog[paste0("num.", race.dist[i, sort(race.lvls)])] <- race.dist[i, race.num]
}

# Age.Grp
for (i in 1:length(age.grp.lvls)) {
  out$demog[paste0("num.age.grp_", i)] <- age.grp.dist[i, age.grp.num]
}

# Age limits
out$demog$ages <- c(18, 65)

# Mortality Rates (Weekly, Age-and-Race-specific)
for (i in 1:(length(mort_rates_annual_raspec) - 1)) {
  for (j in 1:nrow(mort_rates_annual_raspec)) {
    out$demog[paste0("mortrate.",
    "race.", names(mort_rates_annual_raspec)[i + 1],
    ".age.grp_", j)] <-
    round(1 - (1 - mort_rates_annual_raspec[j, i + 1]) ^ (1 / 52), 6)
  }
}

# Mortality Rate (Weekly, population margin)
out$demog$mortrate.marginal <- round(
  1 - (1 - mort_rate_annual_popmargin) ^ (1 / 52), digits = 6
)

# Anal Role Class
for (i in 1:length(role.class.lvls)) {
  out$demog[paste0("role.class.", role.class.dist[i, role.class.lvls])] <-
    role.class.dist[i, role.class.num]
}

print(out$demog)

# ... NETWORK MODEL TARGETS (MAIN PARTNERSHIPS)

out$netmain <- list()
out$netmain$edges <- edges_main
out$netmain$nodefactor_race <- nodefactor_race_main
out$netmain$nodefactor_age.grp <- nodefactor_age.grp_main
out$netmain$nodefactor_degcasl <- nodefacter_degcasl_main
out$netmain$nodematch_race <- nodematch_race.eth_main
out$netmain$nodematch_age.grp <- nodematch_age.grp_main
out$netmain$concurrent <- concurrent_main
out$netmain$durat_wks <- durat_wks_main

# ... NETWORK MODEL TARGETS (CASUAL PARTNERSHIPS)

out$netcasl <- list()
out$netcasl$edges <- edges_casl
out$netcasl$nodefactor_race <- nodefactor_race_casl
out$netcasl$nodefactor_age.grp <- nodefactor_age.grp_casl
out$netcasl$nodefactor_degmain <- nodefacter_degmain_casl
out$netcasl$nodematch_race <- nodematch_race.eth_casl
out$netcasl$nodematch_age.grp <- nodematch_age.grp_casl
out$netcasl$concurrent <- concurrent_casl
out$netcasl$durat_wks <- durat_wks_casl

# ... NETWORK MODEL TARGETS (INSTANTANEOUS PARTNERSHIPS)
out$netinst <- list()
out$netinst$edges <- edges_inst
out$netinst$nodefactor_race <- nodefactor_race_i
out$netinst$nodefactor_age.grp <- nodefactor_age.grp_i
out$netinst$nodefactor_degmain <- nodefactor_degmain_i
out$netinst$nodefactor_degcasl <- nodefactor_degcasl_i

# ... STORE ATTRIBUTES

out$attr <- list()
out$attr$age <- attr_age.yr
out$attr$sqrt.age <- sqrt(attr_age.yr)
out$attr$age.wk <- attr_age.yr * 52
out$attr$age.grp.char <- attr_age.grp
out$attr$age.grp <- as.numeric(
  cut(attr_age.yr, c(18, 25, 35, 45, 55, 65), right = FALSE)
)
out$attr$race <- attr_race
out$attr$deg.main <- attr_deg.main
out$attr$deg.casl <- attr_deg.casl
out$attr$role.class <- attr_role.class
out$attr$diag.status <- attr_diag.status


# ... WRITE

saveRDS(out, here::here("netstats", "netstats.Rds"))
