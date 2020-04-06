# %% PACKAGES ------------------------------------------------------------------

pacman::p_load(
  magrittr,
  data.table,
  summarytools,
  stringr
)

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
  nodematch("role.class", levels = c("R", "I"))

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
  nsteps = 2000,
  ncores = use_ncores,
  skip.dissolution = FALSE,
  nwstats.formula = main_formation_full
)

plot(dx_main, qnts = 0.95)

# need to burn in for a long time to assess duration match
# will be biased downward for early time steps
plot(dx_main, type = "duration")
print(dx_main)


# %% CASUAL PARTNERSHIPS -------------------------------------------------------

casl_formation_full <-
  ~ edges +
  nodefactor("race.eth", levels = NULL) +
  nodefactor("age5", levels = NULL) +
  nodefactor("deg.main", levels = NULL) +
  concurrent +
  nodematch("race.eth", levels = NULL) +
  nodematch("age5") +
  nodematch("role.class", levels = NULL)

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
  nsteps = 2000,
  ncores = use_ncores,
  skip.dissolution = FALSE,
  nwstats.formula = casl_formation_full
)

plot(dx_casl, qnts = 0.95)
print(dx_casl)


# %% ONE-TIME PARTNERSHIPS -----------------------------------------------------

inst_formation_full <-
  ~ edges +
  nodefactor("race.eth", levels = NULL) +
  nodefactor("age5", levels = NULL) +
  nodefactor("deg.main", levels = NULL) +
  nodefactor("deg.casl", levels = NULL) +
  nodematch("role.class", levels = NULL)

inst_formation <-
  ~ edges +
  nodefactor("race.eth", levels = -1) +
  nodefactor("age5", levels = -1) +
  nodefactor("deg.main", levels = -1) +
  nodefactor("deg.casl", levels = -1) +
  nodematch("role.class", diff = TRUE, levels = c("R", "I"))

netstats_inst <-
  with(netstats$netinst,
    c(edges = edges,
      nodefactor_race = nodefactor_race[-1],
      nodefactor_age5 = nodefactor_age5[-1],
      nodefactor_degmain = nodefactor_degmain[-1],
      nodefactor_degcasl = nodefactor_degcasl[-1],
      nodefactor_ai.role = c(0, 0)
))

netstats_inst <- unname(netstats_inst)
print(netstats_inst)

coef_diss_inst <- dissolution_coefs(
  dissolution = ~offset(edges),
  duration = 1,
  d.rate = netstats$demog$mortrate.marginal
)

netest_inst <- netest(
  nw = nw,
  formation = inst_formation,
  target.stats = netstats_inst,
  coef.diss = coef_diss_inst,
  set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
)

summary(netest_inst)
mcmc.diagnostics(netest_inst$fit)

dx_inst <- netdx(
  netest_inst,
  nsims = 10000,
  ncores = use_ncores,
  skip.dissolution = FALSE,
  nwstats.formula = inst_formation_full,
  dynamic = FALSE
)

plot(dx_inst, qnts = 0.95)
print(dx_inst)


# %% Write ------------------------------------------------------------------

netest_out <- list(
  fit_main = netest_main,
  fit_casl = netest_casl,
  fit_inst = netest_inst
)

dx_out <- list(dx_main, dx_casl, dx_inst)

saveRDS(netest_out, "netest/netest.Rds")
saveRDS(dx_out, "netest/netdx.Rds")
