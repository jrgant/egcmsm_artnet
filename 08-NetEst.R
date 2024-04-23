################################################################################
                                 ## PACKAGES ##
################################################################################

pacman::p_load(
  magrittr,
  data.table,
  summarytools,
  stringr,
  foreach,
  doParallel
)

netstats <- readRDS(here::here("netstats", "netstats.Rds"))


################################################################################
                           ## MORTALITY CORRECTION ##
################################################################################

## Tweak the exit rate based on epidemic model simulations that resulted
## in population sizes close to N = 20,000. We want dissolution rates to
## reflect that exit rate.
mort_correct <- 1.285 / 20000


################################################################################
                            ## NETWORK ESTIMATION ##
################################################################################

suppressMessages(library("EpiModelHIVxgc"))

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

# Specification of node mixing terms for age group and race/ethnicity
#   - These specs drop age groups involving age group 5 (55+) and race/eth
#     group 4 (white).
#   - Additionally, nodefactor and nodematch terms are dropped because they
#     are redundant (and lead to model fit issues) in the presence of the
#     nodemix terms. The original model specifications are left commented out
#     for reference.


################################################################################
               ## MAIN PARTNERSHIPS: FORMULAE AND TARGET STATS ##
################################################################################

main_formation <- ~ edges +
  nodemix("age.grp") +
  nodemix("race") +
  nodefactor("deg.casl", levels = -1) +
  nodefactor("diag.status", levels = -1) +
  degrange(from = 3) +
  concurrent +
  nodematch("diag.status") +
  nodematch("role.class", diff = TRUE, levels = I(0:1))


netstats_main <- with(netstats$netmain, c(
  edges = edges,
  mix.age.grp = nodemix_age.grp[-1],
  mix.race = nodemix_race[-1],
  nodefactor.deg.casl = nodefactor_degcasl[-1],
  nodefactor.diag.status = nodefactor_diagstatus[-1],
  `deg3+` = 0,
  concurrent = concurrent,
  nodematch.diag.status = nodematch_diagstatus,
  nodematch.role.class = c("0" = 0, "1" = 0)
))

#netstats_main <- unname(netstats_main)
netstats_main

coef_diss_main <- dissolution_coefs(
  dissolution = ~offset(edges) + offset(nodemix("age.grp")),
  duration = c(netstats$netmain$durat_wks,
               netstats$netmain$durat_wks_byagec[-1]),
  d.rate = netstats$demog$mortrate.marginal + mort_correct
)


################################################################################
              ## CASUAL PARTNERSHIPS: FORMULAE AND TARGET STATS ##
################################################################################

casl_formation <- ~ edges +
  nodemix("age.grp") +
  nodemix("race") +
  nodefactor("deg.main", levels = -1) +
  nodefactor("diag.status", levels = -1) +
  degrange(from = 6) +
  concurrent +
  nodematch("diag.status") +
  nodematch("role.class", diff = TRUE, levels = I(0:1))

netstats_casl <- with(netstats$netcasl, c(
  edges = edges,
  mix.age.grp = nodemix_age.grp[-1],
  mix.race = nodemix_race[-1],
  nodefactor.deg.main = nodefactor_degmain[-1],
  nodefactor.diag.status = nodefactor_diagstatus[-1],
  `deg6+` = 0,
  concurrent = concurrent,
  nodematch.diag.status = nodematch_diagstatus,
  nodematch.role.class = c("0" = 0, "1" = 0)
))

netstats_casl

coef_diss_casl <- dissolution_coefs(
  dissolution =
    ~offset(edges) + offset(nodemix("age.grp")),
  duration = c(
    netstats$netcasl$durat_wks,
    netstats$netcasl$durat_wks_byagec[-1]
  ),
  d.rate = netstats$demog$mortrate.marginal + mort_correct
)


################################################################################
               ## ONE-TIME CONTACTS: FORMULA AND TARGET STATS ##
################################################################################

inst_formation <- ~ edges +
  nodefactor("race", levels = -1) +
  nodefactor("age.grp", levels = -1) +
  nodefactor("deg.main", levels = -1) +
  nodefactor("deg.casl", levels = -1) +
  nodefactor("diag.status", levels = -1) +
  nodematch("role.class", diff = TRUE, levels = I(0:1))

netstats_inst <- with(netstats$netinst, c(
  edges = edges,
  nodefactor.race = nodefactor_race[-1],
  nodefactor.age.grp = nodefactor_age.grp[-1],
  nodefactor.deg.main = nodefactor_degmain[-1],
  nodefactor.deg.casl = nodefactor_degcasl[-1],
  nodefactor.diag.status = nodefactor_diagstatus[-1],
  nodematch.role.class = c("0" = 0, "1" = 0)
))

netstats_inst

coef_diss_inst <- dissolution_coefs(
  dissolution = ~offset(edges),
  duration = 1,
  d.rate = netstats$demog$mortrate.marginal
)


################################################################################
## FIT NETWORK MODELS ##
################################################################################

# List model components in order of: main, casl, inst
formation_formulae <- c("main_formation", "casl_formation", "inst_formation")
target_stats <- c("netstats_main", "netstats_casl", "netstats_inst")
dissolution_coef_fits <- c("coef_diss_main", "coef_diss_casl", "coef_diss_inst")
model_names <- c("fit_main", "fit_casl", "fit_inst")

# Set maximumn number of MCMC iterations
mcmc.maxiterations <- 500

# Fit models in parallel
registerDoParallel(cores = 3)
netest_out <- foreach(i = 1:3, .inorder = FALSE) %dopar% {
  netest(
    nw = nw,
    formation = get(formation_formulae[i]),
    target.stats = get(target_stats[i]),
    coef.diss = get(dissolution_coef_fits[i]),
    set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
  )
}
stopImplicitCluster()

names(netest_out) <- model_names


# %% Write ------------------------------------------------------------------
saveRDS(netest_out, "netest/netest.Rds")
