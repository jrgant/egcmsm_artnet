# %% PACKAGES ------------------------------------------------------------------

library(magrittr)
library(data.table)

degdist <- readRDS("netstats/aggregate_degree_summaries.Rds")

# %% INITIALIZE POPULATION -----------------------------------------------------

# @TODO 2020-01-29:
# - Move population initilization into another file that conforms to
#   the out$demo lists from ARTNet

pop.N <- 20000
set.seed(1971)

# ... RACE and ETHNICITY

# Source: Grey JA, Bernstein KT, Sullivan PS, Kidd SE, Gift TL, Hall EW, et al.
# Rates of Primary and Secondary Syphilis Among White and Black Non-Hispanic
# Men Who Have Sex With Men, United States, 2014. J Acquir Immune Defic Syndr.
# 2017;76: e65–e73. doi:10.1097/QAI.0000000000001508

race.lvls <- c("W", "B", "H", "O")
race.prob <- c(0.586, 0.135, 0.191, 0.088)

race.nums <- cbind(race.lvls, race.prob,
                   num = pop.N * race.prob) %>%
  as.data.frame %>%
  setDT %>%
  .[order(race.lvls)] %>%
  .[, num := as.numeric(as.character(num))]

race.nums
sum(race.prob) == 1

race <- sample(race.lvls, size = pop.N, replace = T, prob = race.prob)
prop.table(table(race))

# ... AGE GROUP

# Source: Jones J, Grey JA, Purcell DW, Bernstein KT, Sullivan PS, Rosenberg
# ES. Estimating Prevalent Diagnoses and Rates of New Diagnoses of HIV at the
# State-Level by Age Group among Men Who Have Sex with Men in the United
# States. Open Forum Infect Dis. 2018 May 29

age.lvls <- list("18-24" = list(min = 18, max = 24),
                 "25-34" = list(min = 25, max = 34),
                 "35-44" = list(min = 35, max = 44),
                 "45-54" = list(min = 45, max = 54),
                 "55+"   = list(min = 55, max = 64))

age.prob <- c(0.129, 0.242, 0.239, 0.216, 0.174)
sum(age.prob) == 1    # check probabilities sum to 1

age5 <- sample(x = names(age.lvls),
                  size = pop.N,
                  replace = T,
                  prob = age.prob)

prop.table(table(age5))

age.yr <- vapply(age5, FUN = function(x) {
  runif(n = 1, min = age.lvls[[x]]$min, age.lvls[[x]]$max + 0.99)
  },
  FUN.VALUE = 23.56)

age.wk <- floor(age.yr * 52)
summary(age.wk)

library(summarytools)
stby(data = age.yr, INDICES = age5, FUN = descr)


# ... INITIAL DEGREE ASSIGNMENTS

maindat <- degdist$main_summaries$sum_main_trunc

deg.main <- sample(maindat[, degmain_trunc2],
                   size = pop.N,
                   replace = T,
                   prob = maindat[, P])

summary(factor(deg.main))

casldat <- degdist$casl_summaries$sum_casl_total

deg.casl <- sample(casldat[, degcasl],
                   size = pop.N,
                   replace = T,
                   prob = casldat[, P])

summary(factor(deg.casl))



# %% TECHNICAL PARAMETERS ------------------------------------------------------

mcmc.maxiterations <- 500
use_ncores <- 4


# %% NETWORK ESTIMATION --------------------------------------------------------

suppressMessages(library("EpiModelHIV"))

# Load ART-Net Partnership Data
pdat <- readRDS("netstats/predictions.Rds")

# Initialize network
nw <- network.initialize(pop.N, directed = F)
nw <- set.vertex.attribute(nw, "race.eth", race)
nw <- set.vertex.attribute(nw, "age5", age5)
nw <- set.vertex.attribute(nw, "age.wk", age.wk)
nw <- set.vertex.attribute(nw, "deg.main", deg.main)
nw <- set.vertex.attribute(nw, "deg.casl", deg.casl)

# %% MAIN PARTNERSHIPS ---------------------------------------------------------

library(stringr)
library(data.table)
library(magrittr)

wts <- prop.table(table(race))

edges_main_calc <-
  as.data.table(pdat$main$deg_byrace) %>%
  .[, match := toupper(str_extract(race.cat, "^[a-z]{1}"))] %>%
  .[order(match)]

edges_main_calc %>%
  .[, weight := wts[match(match, names(wts))]] %>%
  .[, num := weight * pop.N]

print(edges_main_calc)

edges_main <- round(weighted.mean(edges_main_calc$pred,
                                  edges_main_calc$weight) * pop.N / 2)
print(edges_main)

main_concurrent <- round(pdat$main$concurrent_prob * pop.N)

formation <- ~edges +
             nodefactor("race.eth", levels = -1) +
             nodefactor("age5", levels = -1) +
             nodefactor("deg.casl", levels = -1) +
             concurrent +
             degrange(from = 3)


netstats_main <- c(
  edges = edges_main,
  nodefactor_race = round((pdat$main$deg_byrace$pred[-1] *
                           table(race)[-1])),
  nodefactor_age5 = round((pdat$main$deg_byage5$pred[-1] *
                           table(age5)[-1])),
  nodefactor_degcasl = round((pdat$main$deg_bycasltot$pred[-1] *
                               table(deg.casl)[-1])),
  concurrent = main_concurrent,
  degrange = 0)

netstats_main <- unname(netstats_main)
netstats_main

# @TODO 2020-01-28
# add mortality rates
coef_diss_main <- dissolution_coefs(~offset(edges),
                                    pdat$main$mean_durat_wks)

netest_main <- netest(
       nw = nw,
       formation = formation,
       target.stats = netstats_main,
       coef.diss = coef_diss_main,
       set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
       )

summary(netest_main)
saveRDS(netest_main, "netest/netest_main.Rds")

mcmc.diagnostics(netest_main$fit)

dx_main <- netdx(netest_main,
                 nsims = 10,
                 nsteps = 52 * 3,
                 ncores = use_ncores)

saveRDS(dx_main, "netest/dx_main.Rds")
plot(dx_main)

# %% CASUAL PARTNERSHIPS -------------------------------------------------------

edge_casl_calc <-
  as.data.table(pdat$casl$deg_byrace) %>%
  .[, match := toupper(str_extract(race.cat, "^[a-z]{1}"))] %>%
  .[order(match)]

edge_casl_calc %>%
  .[, weight := wts[match(match, names(wts))]] %>%
  .[, num := weight * pop.N]

print(edge_casl_calc)

edges_casl <- round(weighted.mean(edge_casl_calc$pred,
                                  edge_casl_calc$weight) * pop.N / 2)
print(edges_casl)

casl_concurrent <- round(pdat$casl$concurrent_prob * pop.N)
print(casl_concurrent)

formation <- ~edges +
             nodefactor("race.eth", levels = -1) +
             nodefactor("age5", levels = -1) +
             nodefactor("deg.main", levels = -1) +
             concurrent +
             degrange(from = 5)


netstats_casl <- c(
  edges = edges_casl,
  nodefactor_race = round((pdat$casl$deg_byrace$pred[-1] *
                          table(race)[-1])),
  nodefactor_age5 = round((pdat$casl$deg_byage5$pred[-1] *
                          table(age5)[-1])),
  nodefactor_degmain = round((pdat$casl$deg_bymaindegt2$pred[-1] *
                               table(deg.main)[-1])),
  concurrent = casl_concurrent,
  degrange = 0)

netstats_casl <- unname(netstats_casl)
print(netstats_casl)

# @TODO 2020-01-28
# add mortality rates
coef_diss_casl <- dissolution_coefs(~offset(edges),
                                    pdat$casl$mean_durat_wks)

netest_casl <- netest(
       nw = nw,
       formation = formation,
       target.stats = netstats_casl,
       coef.diss = coef_diss_casl,
       set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
       )

summary(netest_casl)
saveRDS(netest_casl, "netest/netest_casl.Rds")

mcmc.diagnostics(netest_casl$fit)

dx_casl <- netdx(netest_casl,
                 nsims = 10,
                 nsteps = 52 * 3,
                 ncores = use_ncores)

plot(dx_casl)
saveRDS(dx_casl, "netest/dx_casl.Rds")

# %% ONE-TIME PARTNERSHIPS -----------------------------------------------------
