################################################################################
                                 ## PACKAGES ##
################################################################################

pacman::p_load(
  magrittr,
  data.table,
  summarytools,
  stringr
)

netstats <- readRDS(here::here("netstats", "netstats.Rds"))


################################################################################
                           ## TECHNICAL PARAMETERS ##
################################################################################

mcmc.maxiterations <- 500


################################################################################
                            ## NETWORK ESTIMATION ##
################################################################################

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


################################################################################
                            ## MAIN PARTNERSHIPS ##
################################################################################

main_formation <-
  ~ edges +
    nodefactor("age.grp", levels = I(2:5)) +
    nodefactor("race", levels = I(2:4)) +
    nodefactor("deg.casl", levels = I(1:5)) +
    nodefactor("diag.status", levels = I(1)) +
    degrange(from = 3) +
    concurrent +
    nodematch("race") +
    nodematch("age.grp") +
    nodematch("diag.status") +
    nodematch("role.class", diff = TRUE, levels = c(1, 2))

netstats_main <-
  with(netstats$netmain,
    c(edges = edges,
      nodefactor_age.grp = nodefactor_age.grp[-1],
      nodefactor_race = nodefactor_race[-1],
      nodefactor_degcasl = nodefactor_degcasl[-1],
      nodefactor_diag.status = nodefactor_diagstatus[-1],
      degrange = 0,
      concurrent = concurrent,
      nodematch_race = nodematch_race,
      nodematch_age.grp = nodematch_age.grp,
      nodematch_diag.status = nodematch_diagstatus,
      nodematch_ai.role = c(0, 0)
))

netstats_main <- unname(netstats_main)
netstats_main

coef_diss_main <- dissolution_coefs(
  dissolution = ~offset(edges) + offset(nodefactor("age.grp", levels = I(2:5))),
  duration = c(
    netstats$netmain$durat_wks,
    netstats$netmain$durat_wks_byage[-1]
  ),
  d.rate = netstats$demog$mortrate.marginal
)

netest_main <- netest(
  nw = nw,
  formation = main_formation,
  target.stats = netstats_main,
  coef.diss = coef_diss_main,
  set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
)


################################################################################
                           ## CASUAL PARTNERSHIPS ##
################################################################################

casl_formation <-
  ~ edges +
    nodefactor("age.grp", levels = I(2:5)) +
    nodefactor("race", levels = I(2:4)) +
    nodefactor("deg.main", levels = I(1:2)) +
    nodefactor("diag.status", level = I(1)) +
    degrange(from = 6) +
    concurrent +
    nodematch("race") +
    nodematch("age.grp") +
    nodematch("diag.status") +
    nodematch("role.class", diff = TRUE, levels = c(1, 2))

netstats_casl <-
  with(netstats$netcasl,
    c(edges = edges,
      nodefactor_age.grp = nodefactor_age.grp[-1],
      nodefactor_race = nodefactor_race[-1],
      nodefactor_degmain = nodefactor_degmain[-1],
      nodefactor_diag.status = nodefactor_diagstatus[-1],
      degrange = 0,
      concurrent = concurrent,
      nodematch_race = nodematch_race,
      nodematch_age.grp = nodematch_age.grp,
      nodematch_diag.status = nodematch_diagstatus,
      nodematch_ai.role = c(0, 0)
))

netstats_casl <- unname(netstats_casl)
print(netstats_casl)

coef_diss_casl <- dissolution_coefs(
  dissolution = ~offset(edges) + offset(nodefactor("age.grp", levels = I(2:5))),
  duration = c(
    netstats$netcasl$durat_wks,
    netstats$netcasl$durat_wks_byage[-1]
  ),
  d.rate = netstats$demog$mortrate.marginal
)

netest_casl <- netest(
  nw = nw,
  formation = casl_formation,
  target.stats = netstats_casl,
  coef.diss = coef_diss_casl,
  set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
)


################################################################################
                          ## ONE-TIME PARTNERSHIPS ##
################################################################################

inst_formation <-
  ~ edges +
    nodefactor("race", levels = I(2:4)) +
    nodefactor("age.grp", levels = I(2:5)) +
    nodefactor("deg.main", levels = I(1:2)) +
    nodefactor("deg.casl", levels = I(1:5)) +
    nodefactor("diag.status", levels = I(1)) +
    nodematch("role.class", diff = TRUE, levels = c(1, 2))

netstats_inst <-
  with(netstats$netinst,
    c(edges = edges,
      nodefactor_race = nodefactor_race[-1],
      nodefactor_age.grp = nodefactor_age.grp[-1],
      nodefactor_degmain = nodefactor_degmain[-1],
      nodefactor_degcasl = nodefactor_degcasl[-1],
      nodefactor_diagstatus = nodefactor_diagstatus[-1],
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


# %% Write ------------------------------------------------------------------

netest_out <- list(
  fit_main = netest_main,
  fit_casl = netest_casl,
  fit_inst = netest_inst
)

saveRDS(netest_out, "netest/netest.Rds")
