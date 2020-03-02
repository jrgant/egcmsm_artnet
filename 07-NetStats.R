# %% SETUP ---------------------------------------------------------------------

pacman::p_load(data.table,
               EpiModelHIV,
               readxl,
               magrittr)


# %% INPUTS --------------------------------------------------------------------

network_size <- 20000

# Degree distributions
degdist <- readRDS("netstats/aggregate_degree_summaries.Rds")
str(degdist)

# Predictions based on ART-Net Partnership Data
pdat <- readRDS("netstats/predictions.Rds")
str(pdat)


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

age5.lvls <- list("18-24" = list(min = 18, max = 24),
                  "25-34" = list(min = 25, max = 34),
                  "35-44" = list(min = 35, max = 44),
                  "45-54" = list(min = 45, max = 54),
                  "55+"   = list(min = 55, max = 64))

age5.prob <- c(0.129, 0.242, 0.239, 0.216, 0.174)
age5.num <- round(num * age5.prob)
sum(age5.num) == num

age5.dist <- data.table(age5.lvls = names(age5.lvls), age5.prob, age5.num)

# ... JOINT RACExAGE

raceage.dist <- cbind(
  expand.grid(race = race.dist$race.lvls,
              age5 = age5.dist$age5.lvls),
  expand.grid(race.prob = race.dist$race.prob,
              age5.prob = age5.dist$age5.prob)) %>%
  setDT %>%
  .[, raceage.prob := race.prob * age5.prob]

print(raceage.dist)
sum(raceage.dist$raceage.prob)

# ... ANAL SEX ROLE

# Source: ArtNet

role.class.lvls <- pdat$demo$ai.role.pr$anal.sex.role
role.class.prob <- pdat$demo$ai.role.pr$P
role.class.num <- round(num * role.class.prob)
sum(role.class.num) == num

role.class.dist <- data.table(role.class.lvls,
                              role.class.prob,
                              role.class.num)[order(role.class.lvls)]


# ... MAIN DEGREE DISTRIBUTION
maindeg.dist <- pdat$main$degprob[
  order(race.cat, age5),
  .(maindeg.prob = weighted.mean(
      preds,
      raceage.dist[order(race, age5), raceage.prob])
    ),
  outcome]

cbind(degdist$main_artnet_sum$sum_main_trunc, maindeg.dist)
print(maindeg.dist)

# ... CASUAL DEGREE DISTRIBUTION

# casldeg.dist <- degdist$casl_summaries$sum_casl_total[, -c("N")] %>%
#                 .[, num := round(P * num)]

casldeg.dist <- pdat$casl$degprob[
  order(race.cat, age5),
  .(casldeg.prob = weighted.mean(
      preds,
      raceage.dist[order(race, age5), raceage.prob])
    ),
  outcome]

cbind(degdist$casl_artnet_sum$sum_casl_total, casldeg.dist)
print(casldeg.dist)


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
  .[, age5 := dplyr::case_when(
    agestart %in% c(18:24) ~ "[18, 25)",
    agestart %in% c(25:34) ~ "[25, 35)",
    agestart %in% c(35:44) ~ "[35, 45)",
    agestart %in% c(45:54) ~ "[45, 55)",
    agestart >= 55 ~ "55+"
    )]

mort_rates <- lt_wide[, .(B = mean(nhblack),
                          H = mean(hisp),
                          W = mean(nhwhite)),
                          by = age5]

# Assign Other the same mortality rate as Whites
mort_rates[, O := W]
mort_rates_annual_raspec <- as.data.frame(mort_rates)

print(mort_rates_annual_raspec)

mort_rates_long <- melt(
  as.data.table(mort_rates_annual_raspec),
  id.vars = "age5",
  measure.vars = c("B", "H", "O", "W")
  ) %>%
  .[, .(age5, race = variable, rate_percap = value)] %>%
  .[, age5_num := rep(1:5, 4)]

mort_rates_long <- cbind(
  mort_rates_long,
  raceage.dist[order(race, age5), .(age5, race, raceage.prob)]
)

print(mort_rates_long)
sum(mort_rates_long$raceage.prob) == 1

mort_rate_annual_popmargin <-
  mort_rates_long[, weighted.mean(rate_percap, raceage.prob)]

print(mort_rate_annual_popmargin)


# %% INITIALIZE NODE ATTRIBUTES ------------------------------------------------

set.seed(19410524)

## Assign age
attr_age5 <- sample(names(age5.lvls),
                    size = num,
                    prob = age5.prob,
                    replace = T)

attr_age.yr <- vapply(attr_age5, FUN = function(x) {
  runif(n = 1, min = age5.lvls[[x]]$min, age5.lvls[[x]]$max + 0.99)
  },
  FUN.VALUE = 23.56) %>%
  unname

## Assign race/ethnicity
attr_race <- sample(1:4,
                    size = num,
                    prob = race.dist[, race.prob],
                    replace = T)

## Main degree
attr_deg.main <- sample(0:2,
                        size = num,
                        replace = T,
                        prob = maindeg.dist[, maindeg.prob]) %>%
                 as.integer

## Casual degree
attr_deg.casl <- sample(0:5,
                        size = num,
                        replace = T,
                        prob = casldeg.dist[, casldeg.prob]) %>%
                  as.integer

## Anal Sex Role
# Source: ARTNet
attr_role.class <- sample(role.class.dist[, role.class.lvls],
                          size = num,
                          prob = role.class.dist[, role.class.prob],
                          replace = T)

# %% PARTNERSHIP REUSABLES -----------------------------------------------------

race.wt <- raceage.dist[, weight := raceage.prob / sum(raceage.prob),
keyby = race][, .(race, age5, weight)]
print(race.wt)

age5.wt <- raceage.dist[, weight := raceage.prob / sum(raceage.prob),
keyby = age5][, .(race, age5, weight)]
print(age5.wt)


# %% MAIN PARTNERSHIPS ---------------------------------------------------------

# calculate mean degree by race (marginalized over age5)
race_preds_m <- cbind(pdat$main$degpred_byra[order(race.cat, age5)],
                    race.wt[order(race, age5)]) %>%
      .[, .(pred_margage = weighted.mean(pred_mdeg_byra, weight)), .(race)]

print(race_preds_m)

# calculate mean degree by age group (marginalized over race/ethnicity)
age5_preds_m <- cbind(pdat$main$degpred_byra[order(race.cat, age5)],
                    age5.wt[order(race, age5)]) %>%
      .[, .(pred_margrace = weighted.mean(pred_mdeg_byra, weight)), .(age5)]

print(age5_preds_m)

# calculate mean degree by casual degree
main_bycasl_tab <- pdat$main$degpred_bycasl
print(main_bycasl_tab)

## edges_main

edge_pred_m <- cbind(
  pred = pdat$main$degpred_byra[order(race.cat, age5)] %>%
    .[, pred_mdeg_byra],
  raceage.prob = raceage.dist[order(race, age5)][, raceage.prob]) %>%
      as.data.table

print(edge_pred_m)
sum(edge_pred_m$raceage.prob) == 1

edges_main <- round(
  weighted.mean(
    edge_pred_m$pred,
    edge_pred_m$raceage.prob
    ) * num / 2)

print(edges_main)

## nodefactor_race

nodefactor_race_main <-
  round(race_preds_m$pred_margage * race.dist$race.prob * num)

print(sum(nodefactor_race_main) / 2)

## nodefactor_age5
nodefactor_age5_main <-
  round(age5_preds_m$pred_margrace * age5.dist$age5.prob * num)

print(sum(nodefactor_age5_main) / 2)

## nodefactor_degcasl
nodefacter_degcasl_main <-
  round(pdat$main$degpred_bycasl$pred_mdeg_bycasl *
            casldeg.dist$casldeg.prob * num)

print(sum(nodefacter_degcasl_main) / 2)

## nodematch_race.eth
nodematch_race.eth_main <-
  round(pdat$main$racematch[samerace == 1, P] * edges_main)

print(nodematch_race.eth_main)

## nodematch_age5
nodematch_age5_main <-
  round(pdat$main$age5match[sameage == 1, P] * edges_main)

print(nodematch_age5_main)

## concurrent
concurrent_main <- round(pdat$main$concurrent * num)
print(concurrent_main)

## duration
durat_wks_main <- pdat$main$durat_wks
print(durat_wks_main)


# %% CASUAL PARTNERSHIPS -------------------------------------------------------

# calculate mean degree by race (marginalized over age5)
race_preds_c <- cbind(pdat$casl$degpred_byra[order(race.cat, age5)],
                    race.wt[order(race, age5)]) %>%
      .[, .(pred_margage = weighted.mean(pred_cdeg_byra, weight)), .(race)]

print(race_preds_c)

# calculate mean degree by age group (marginalized over race/ethnicity)
age5_preds_c <- cbind(pdat$casl$degpred_byra[order(race.cat, age5)],
                      age5.wt[order(race, age5)]) %>%
      .[, .(pred_margrace = weighted.mean(pred_cdeg_byra, weight)), .(age5)]

print(age5_preds_c)

# calculate mean degree by casual degree
casl_bymain_tab <- pdat$casl$degpred_bymain
print(casl_bymain_tab)

## edges_casl

edge_pred_c <- cbind(
  pred = pdat$casl$degpred_byra[order(race.cat, age5)] %>%
    .[, pred_cdeg_byra],
  raceage.prob = raceage.dist[order(race, age5)][, raceage.prob]) %>%
      as.data.table

print(edge_pred_c)
sum(edge_pred_c$raceage.prob) == 1

edges_casl <- round(
  weighted.mean(
    edge_pred_c$pred,
    edge_pred_c$raceage.prob
    ) * num / 2)

print(edges_casl)

## nodefactor_race

nodefactor_race_casl <-
  round(race_preds_c$pred_margage * race.dist$race.prob * num)

print(sum(nodefactor_race_casl) / 2)

## nodefactor_age5
nodefactor_age5_casl <-
  round(age5_preds_c$pred_margrace * age5.dist$age5.prob * num)

print(sum(nodefactor_age5_casl) / 2)

## nodefactor_degcasl
nodefacter_degmain_casl <-
  round(pdat$casl$degpred_bymain$pred_cdeg_bymain *
          maindeg.dist$maindeg.prob * num)

print(sum(nodefacter_degmain_casl) / 2)

## nodematch_race.eth
nodematch_race.eth_casl <-
  round(pdat$casl$racematch[samerace == 1, P] * edges_casl)

print(nodematch_race.eth_casl)

## nodematch_age5
nodematch_age5_casl <-
  round(pdat$casl$age5match[sameage == 1, P] * edges_casl)

print(nodematch_age5_casl)

## concurrent
concurrent_casl <- round(pdat$casl$concurrent * num)
print(concurrent_casl)

## duration
durat_wks_casl <- pdat$casl$durat_wks
print(durat_wks_casl)


# %% ONE-TIME PARTNERSHIPS -----------------------------------------------------



# %% SAVE PARAMETERS TO FILE ---------------------------------------------------

out <- list()

# ... STORE INPUT PARAMETERS

out$inputs <- list()
out$inputs$race.dist <- race.dist[, -c("race.num")]
out$inputs$age5.dist <- age5.dist[, -c("age5.num")]
out$inputs$main.dist <- maindeg.dist[, .(degmain_trunc2 = outcome,
                                         prob = maindeg.prob)]
out$inputs$casl.dist <- casldeg.dist[, .(degcasl = outcome,
                                         prob = casldeg.prob)]
out$inputs$role.class.dist <- role.class.dist[, -c("role.class.num")]

out$inputs

# ... STORE DEMOGRAPHICS

out$demog <- list()
out$demog$num <- num

# Race/ethnicity
for (i in 1:length(race.lvls)) {
  out$demog[paste0("num.", race.dist[i, race.lvls])] <- race.dist[i, race.num]
}

# Age5
for (i in 1:length(age5.lvls)) {
  out$demog[paste0("num.age5_", i)] <- age5.dist[i, age5.num]
}

# Mortality Rates (Weekly, Age-and-Race-specific)
for (i in 1:(length(mort_rates_annual_raspec) - 1)) {
  for (j in 1:nrow(mort_rates_annual_raspec)) {
    out$demog[paste0("mortrate.",
    "race.", names(mort_rates_annual_raspec)[i + 1],
    ".age5_", j)] <-
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
out$netmain$nodefactor_age5 <- nodefactor_age5_main
out$netmain$nodefactor_degcasl <- nodefacter_degcasl_main
out$netmain$nodematch_race.eth <- nodematch_race.eth_main
out$netmain$nodematch_age5 <- nodematch_age5_main
out$netmain$concurrent <- concurrent_main
out$netmain$durat_wks <- durat_wks_main

# ... NETWORK MODEL TARGETS (CASUAL PARTNERSHIPS)

out$netcasl <- list()
out$netcasl$edges <- edges_casl
out$netcasl$nodefactor_race <- nodefactor_race_casl
out$netcasl$nodefactor_age5 <- nodefactor_age5_casl
out$netcasl$nodefactor_degmain <- nodefacter_degmain_casl
out$netcasl$nodematch_race.eth <- nodematch_race.eth_casl
out$netcasl$nodematch_age5 <- nodematch_age5_casl
out$netcasl$concurrent <- concurrent_casl
out$netcasl$durat_wks <- durat_wks_casl


# ... STORE ATTRIBUTES

out$attr <- list()
out$attr$age <- attr_age.yr
out$attr$age.wk <- attr_age.yr * 52
out$attr$age5 <- attr_age5
out$attr$race <- attr_race
out$attr$deg.main <- attr_deg.main
out$attr$deg.casl <- attr_deg.casl
out$attr$role.class <- attr_role.class

# print(out$attr)

# ... WRITE

saveRDS(out, here::here("netstats", "netstats.Rds"))
