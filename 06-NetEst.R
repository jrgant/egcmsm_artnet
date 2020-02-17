# %% PACKAGES ------------------------------------------------------------------

pacman::p_load(magrittr,
               data.table,
               summarytools,
               stringr)

# Load degree distributions
degdist <- readRDS("netstats/aggregate_degree_summaries.Rds")

# Load ART-Net Partnership Data
pdat <- readRDS("netstats/predictions.Rds")

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
# 2017;76: e65-e73. doi:10.1097/QAI.0000000000001508

race.lvls <- c("W", "B", "H", "O")
race.prob <- c(0.586, 0.135, 0.191, 0.088)

race.nums <- as.data.frame(
  cbind(race.lvls,
        race.prob,
        num = pop.N * race.prob)
        ) %>%
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

sqrt.age.yr <- sqrt(age.yr)

age.wk <- floor(age.yr * 52)
summary(age.wk)

stby(data = age.yr, INDICES = age5, FUN = descr)


# ... INITIAL ANAL SEX ROLE PROPORTIONS
# Source: ARTNet
ai.role.lvls <- pdat$demo$ai.role$anal.sex.role
ai.role.pr <- pdat$demo$ai.role.pr$P

ai.role <- sample(ai.role.lvls, size = pop.N, replace = T, prob = ai.role.pr)

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

# Initialize network
nw <- network.initialize(pop.N, directed = F)
nw <- set.vertex.attribute(nw, "race.eth", race)
nw <- set.vertex.attribute(nw, "age5", age5)
nw <- set.vertex.attribute(nw, "age.yr", age.yr)
nw <- set.vertex.attribute(nw, "sqrt.age.yr", sqrt.age.yr)
nw <- set.vertex.attribute(nw, "age.wk", age.wk)
nw <- set.vertex.attribute(nw, "deg.main", deg.main)
nw <- set.vertex.attribute(nw, "deg.casl", deg.casl)
nw <- set.vertex.attribute(nw, "ai.role", ai.role)

# %% MAIN PARTNERSHIPS ---------------------------------------------------------

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

sort(unique(nw %v% "age5"))

agecomb_levels <- c("11",
                    "12", "22",
                    "13", "23", "33",
                    "14", "24", "34", "44",
                    "15", "25", "35", "45", "55")

main_formation <- ~edges +
                  nodefactor("race.eth", levels = -1) +
                  nodefactor("age5", levels = -1) +
                  nodefactor("deg.casl", levels = -1) +
                  concurrent +
                  nodematch("race.eth") +
                  nodematch("age5") +
                  nodematch("ai.role", levels = c("R", "I"))
                  

netstats_main <- c(
  edges = edges_main,
  nodefactor_race =
    round(pdat$main$deg_byrace$pred[-1] * table(race)[-1]),
  nodefactor_age5 =
    round(pdat$main$deg_byage5$pred[-1] * table(age5)[-1]),
  nodefactor_degcasl =
    round(pdat$main$deg_bycasltot$pred[-1] * table(deg.casl)[-1]),
  concurrent = main_concurrent,
  nodematch_race.eth = round(pdat$main$racematch[samerace == 1, P] * edges_main),
  nodematch_age5 = round(pdat$main$age5match[sameage == 1, P] * edges_main),
  nodematch_ai.role = 0
)

netstats_main <- unname(netstats_main)
netstats_main

# @TODO 2020-01-28
# add mortality rates
coef_diss_main <- dissolution_coefs(~offset(edges),
                                    pdat$main$mean_durat_wks)

netest_main <- netest(
       nw = nw,
       formation = main_formation,
       target.stats = netstats_main,
       coef.diss = coef_diss_main,
       # edapprox = T,
       set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
       )

summary(netest_main)

mcmc.diagnostics(netest_main$fit)

dx_main <- netdx(netest_main,
                 nsims = 10,
                 nsteps = 52 * 3,
                 ncores = use_ncores)

plot(dx_main, qnts = 0.95)

saveRDS(netest_main, "netest/netest_main.Rds")
saveRDS(dx_main, "netest/dx_main.Rds")


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

# cp_agemix <- pdat$casl$agemix[, agecomb := factor(agecomb, agecomb_levels)] %>%
#   .[order(agecomb)]
# 
# print(cp_agemix)
# sum(cp_agemix$P)


casl_formation <- ~edges +
                  nodefactor("race.eth", levels = -1) +
                  nodefactor("age5", levels = -1) +
                  nodefactor("deg.main", levels = -1) +
                  concurrent +
                  degrange(from = 4) +
                  nodematch("race.eth")

netstats_casl <- c(
  edges = edges_casl,
  nodefactor_race =
    round((pdat$casl$deg_byrace$pred[-1] * table(race)[-1])),
  nodefactor_age5 =
    round((pdat$casl$deg_byage5$pred[-1] * table(age5)[-1])),
  nodefactor_degmain =
    round((pdat$casl$deg_bymaindegt2$pred[-1] * table(deg.main)[-1])),
  concurrent = casl_concurrent,
  degrange = round(
    sum(degdist$casl_summaries$sum_casl_total$P[4:5]) * edges_casl
    ),
  nodematch_race.eth = round(pdat$casl$racematch[samerace == 1, P] * edges_casl)
  )
#nodematch_age5 = round(pdat$casl$age5match[sameage == 1, P] * edges_casl)

netstats_casl <- unname(netstats_casl)
print(netstats_casl)

# @TODO 2020-01-28
# add mortality rates
coef_diss_casl <- dissolution_coefs(~offset(edges),
                                    pdat$casl$mean_durat_wks)

netest_casl <- netest(
       nw = nw,
       formation = casl_formation,
       target.stats = netstats_casl,
       coef.diss = coef_diss_casl,
       edapprox = T,
       set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
       )

summary(netest_casl)

mcmc.diagnostics(netest_casl$fit)

dx_casl <- netdx(netest_casl,
                 nsims = 10,
                 nsteps = 52 * 3,
                 ncores = use_ncores)

plot(dx_casl, qnts = 0.95)

saveRDS(netest_casl, "netest/netest_casl.Rds")
saveRDS(dx_casl, "netest/dx_casl.Rds")


# %% ONE-TIME PARTNERSHIPS -----------------------------------------------------
