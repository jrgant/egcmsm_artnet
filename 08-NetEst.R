# %% PACKAGES ------------------------------------------------------------------

pacman::p_load(magrittr,
               data.table,
               summarytools,
               stringr)

netstats <- readRDS(here::here("netstats", "netstats.Rds"))


# %% TECHNICAL PARAMETERS ------------------------------------------------------

mcmc.maxiterations <- 500
use_ncores <- 6


# %% NETWORK ESTIMATION --------------------------------------------------------

suppressMessages(library("EpiModelHIV"))

# Initialize network
nw <- network.initialize(netstats$demog$num, directed = F)
nw <- set.vertex.attribute(nw, "race.eth", netstats$attr$race)
nw <- set.vertex.attribute(nw, "age", netstats$attr$age)
nw <- set.vertex.attribute(nw, "age.wk", netstats$attr$age.wk)
nw <- set.vertex.attribute(nw, "age5", netstats$attr$age5)
nw <- set.vertex.attribute(nw, "deg.main", netstats$attr$deg.main)
nw <- set.vertex.attribute(nw, "deg.casl", netstats$attr$deg.casl)
nw <- set.vertex.attribute(nw, "role.class", netstats$attr$role.class)


# %% MAIN PARTNERSHIPS ---------------------------------------------------------

main_formation_full <-
  ~ edges +
  nodefactor("race.eth", levels = NULL) +
  nodefactor("age5", levels = NULL) +
  nodefactor("deg.casl", levels = NULL) +
  concurrent +
  nodematch("race.eth") +
  nodematch("age5") +
  nodematch("role.class", diff = TRUE, levels = c("R", "I"))

main_formation <-
  ~ edges +
  nodefactor("race.eth", levels = -1) +
  nodefactor("age5", levels = -1) +
  nodefactor("deg.casl", levels = -1) +
  concurrent +
  nodematch("race.eth") +
  nodematch("age5") +
  nodematch("role.class", diff = TRUE, levels = c("R", "I"))

netstats_main <-
  with(netstats$netmain,
    c(edges = edges,
      nodefactor_race = nodefactor_race[-1],
      nodefactor_age5 = nodefactor_age5[-1],
      nodefactor_degcasl = nodefactor_degcasl[-1],
      concurrent = concurrent,
      nodematch_race.eth = nodematch_race.eth,
      nodematch_age5 = nodematch_age5,
      nodematch_ai.role = c(0, 0)
))

netstats_main <- unname(netstats_main)
netstats_main

coef_diss_main <- dissolution_coefs(
  dissolution = ~offset(edges),
  duration = netstats$netmain$durat_wks,
  d.rate = netstats$demog$mortrate.marginal
)

netest_main <- netest(
  nw = nw,
  formation = main_formation,
  target.stats = netstats_main,
  coef.diss = coef_diss_main,
  set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
)

summary(netest_main)

mcmc.diagnostics(netest_main$fit)

dx_main <- netdx(
  netest_main,
  nsims = 10,
  nsteps = 500,
  ncores = use_ncores,
  skip.dissolution = FALSE,
  nwstats.formula = main_formation_full
)

plot(dx_main, qnts = 0.95)

# need to burn in for a long time to assess duration match
# will be biased downward for early time steps
plot(dx_main, type = "duration")
print(dx_main)

saveRDS(netest_main, "netest/netest_main.Rds")
saveRDS(dx_main, "netest/dx_main.Rds")


# %% CASUAL PARTNERSHIPS -------------------------------------------------------

casl_formation_full <-
  ~ edges +
  nodefactor("race.eth", levels = NULL) +
  nodefactor("age5", levels = NULL) +
  nodefactor("deg.main", levels = NULL) +
  concurrent +
  nodematch("race.eth", levels = NULL) +
  nodematch("age5") +
  nodematch("role.class", diff = TRUE, levels = c("R", "I"))

casl_formation <-
  ~ edges +
  nodefactor("race.eth", levels = -1) +
  nodefactor("age5", levels = -1) +
  nodefactor("deg.main", levels = -1) +
  concurrent +
  nodematch("race.eth") +
  nodematch("age5") +
  nodematch("role.class", diff = TRUE, levels = c("R", "I"))

netstats_casl <-
  with(netstats$netcasl,
    c(edges = edges,
      nodefactor_race = nodefactor_race[-1],
      nodefactor_age5 = nodefactor_age5[-1],
      nodefactor_degmain = nodefactor_degmain[-1],
      concurrent = concurrent,
      nodematch_race.eth = nodematch_race.eth,
      nodematch_age5 = nodematch_age5,
      nodematch_ai.role = c(0, 0)
))


netstats_casl <- unname(netstats_casl)
print(netstats_casl)

coef_diss_casl <- dissolution_coefs(
  dissolution = ~offset(edges),
  duration = netstats$netcasl$durat_wks,
  d.rate = netstats$demog$mortrate.marginal
)

netest_casl <- netest(
  nw = nw,
  formation = casl_formation,
  target.stats = netstats_casl,
  coef.diss = coef_diss_casl,
  set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
)

summary(netest_casl)

mcmc.diagnostics(netest_casl$fit)

dx_casl <- netdx(
  netest_casl,
  nsims = 10,
  nsteps = 500,
  ncores = use_ncores,
  skip.dissolution = FALSE,
  nwstats.formula = casl_formation_full
)

plot(dx_casl, qnts = 0.95)
print(dx_casl)

saveRDS(netest_casl, "netest/netest_casl.Rds")
saveRDS(dx_casl, "netest/dx_casl.Rds")


# %% ONE-TIME PARTNERSHIPS -----------------------------------------------------
