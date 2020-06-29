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
nw <- network.initialize(netstats$demog$num, directed = FALSE)
nw <- set.vertex.attribute(nw, "race", netstats$attr$race)
nw <- set.vertex.attribute(nw, "age", netstats$attr$age)
nw <- set.vertex.attribute(nw, "age.wk", netstats$attr$age.wk)
nw <- set.vertex.attribute(nw, "age.grp", netstats$attr$age.grp)
nw <- set.vertex.attribute(nw, "diag.status", netstats$attr$diag.status)
nw <- set.vertex.attribute(nw, "deg.main", netstats$attr$deg.main)
nw <- set.vertex.attribute(nw, "deg.casl", netstats$attr$deg.casl)
nw <- set.vertex.attribute(nw, "role.class", netstats$attr$role.class)


# %% NOTE ------------------------------------------------------------------

# Role class coding in model formulas:
#   The levels argument in nodematch("role.class") indicates
#   1 for insertive-only agents (coded as 0 in role.class) and
#   2 for receptive-only agents (coded as 1 in role.class)


# %% MAIN PARTNERSHIPS ---------------------------------------------------------

main_formation_full <-
  ~ edges +
  nodefactor("race", levels = NULL) +
  nodefactor("age.grp", levels = NULL) +
  nodefactor("deg.casl", levels = NULL) +
  degrange(from = 3) +
  concurrent +
  nodematch("race") +
  nodematch("age.grp") +
  nodematch("role.class", levels = c(1, 2))

main_formation <-
  ~ edges +
  nodefactor("race", levels = I(2:4)) +
  nodefactor("age.grp", levels = I(2:5)) +
  nodefactor("deg.casl", levels = I(1:5)) +
  degrange(from = 3) +
  concurrent +
  nodematch("race") +
  nodematch("age.grp") +
  nodematch("role.class", diff = TRUE, levels = c(1, 2))

netstats_main <-
  with(netstats$netmain,
    c(edges = edges,
      nodefactor_race = nodefactor_race[-1],
      nodefactor_age.grp = nodefactor_age.grp[-1],
      nodefactor_degcasl = nodefactor_degcasl[-1],
      degrange = 0,
      concurrent = concurrent,
      nodematch_race = nodematch_race,
      nodematch_age.grp = nodematch_age.grp,
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
  nsteps = 100,
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
  nodefactor("race", levels = NULL) +
  nodefactor("age.grp", levels = NULL) +
  nodefactor("deg.main", levels = NULL) +
  degrange(from = 6) +
  concurrent +
  nodematch("race", levels = NULL) +
  nodematch("age.grp") +
  nodematch("role.class", levels = c(1, 2))

casl_formation <-
  ~ edges +
  nodefactor("race", levels = I(2:4)) +
  nodefactor("age.grp", levels = I(2:5)) +
  nodefactor("deg.main", levels = I(1:2)) +
  degrange(from = 6) +
  concurrent +
  nodematch("race") +
  nodematch("age.grp") +
  nodematch("role.class", diff = TRUE, levels = c(1, 2))

netstats_casl <-
  with(netstats$netcasl,
    c(edges = edges,
      nodefactor_race = nodefactor_race[-1],
      nodefactor_age.grp = nodefactor_age.grp[-1],
      nodefactor_degmain = nodefactor_degmain[-1],
      degrange = 0,
      concurrent = concurrent,
      nodematch_race = nodematch_race,
      nodematch_age.grp = nodematch_age.grp,
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
  nsteps = 100,
  ncores = use_ncores,
  skip.dissolution = FALSE,
  nwstats.formula = casl_formation_full
)

plot(dx_casl, qnts = 0.95)
print(dx_casl)


# %% ONE-TIME PARTNERSHIPS -----------------------------------------------------

inst_formation_full <-
  ~ edges +
  nodefactor("race", levels = NULL) +
  nodefactor("age.grp", levels = NULL) +
  nodefactor("deg.main", levels = NULL) +
  nodefactor("deg.casl", levels = NULL) +
  nodematch("role.class", levels = c(1, 2))

inst_formation <-
  ~ edges +
  nodefactor("race", levels = I(2:4)) +
  nodefactor("age.grp", levels = I(2:5)) +
  nodefactor("deg.main", levels = I(1:2)) +
  nodefactor("deg.casl", levels = I(1:5)) +
  nodematch("role.class", diff = TRUE, levels = c(1, 2))

netstats_inst <-
  with(netstats$netinst,
    c(edges = edges,
      nodefactor_race = nodefactor_race[-1],
      nodefactor_age.grp = nodefactor_age.grp[-1],
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
  nsims = 1000,
  ncores = use_ncores,
  skip.dissolution = FALSE,
  nwstats.formula = inst_formation_full,
  dynamic = FALSE
)

plot(dx_inst, sim.lines = TRUE, alpha = 0.4)
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
