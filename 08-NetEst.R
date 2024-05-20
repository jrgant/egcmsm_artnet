################################################################################
                                 ## PACKAGES ##
################################################################################

pacman::p_load(
  magrittr,
  data.table,
  summarytools,
  stringr,
  foreach,
  doParallel,
  EpiModel
)

netstats <- readRDS(here::here("netstats", "netstats.Rds"))


################################################################################
                           ## MORTALITY CORRECTION ##
################################################################################

## Tweak the exit rate based on epidemic model simulations that resulted
## in population sizes close to N = 20,000. We want dissolution rates to
## reflect that exit rate.
mort_correct <- 1.285 / netstats$demog$num


################################################################################
                            ## NETWORK ESTIMATION ##
################################################################################

# Initialize network
nw <- EpiModel::network_initialize(netstats$demog$num)
nw <- EpiModel::set_vertex_attribute(nw, "race", netstats$attr$race)
nw <- EpiModel::set_vertex_attribute(nw, "age", netstats$attr$age)
nw <- EpiModel::set_vertex_attribute(nw, "age.wk", netstats$attr$age.wk)
nw <- EpiModel::set_vertex_attribute(nw, "age.grp", netstats$attr$age.grp)
nw <- EpiModel::set_vertex_attribute(nw, "diag.status", netstats$attr$diag.status)
nw <- EpiModel::set_vertex_attribute(nw, "deg.main", netstats$attr$deg.main)
nw <- EpiModel::set_vertex_attribute(nw, "deg.casl", netstats$attr$deg.casl)
nw <- EpiModel::set_vertex_attribute(nw, "role.class", netstats$attr$role.class)


################################################################################
## MIXING TERM ORDERING ##
################################################################################

## Specify the order of the mixing target statistics (lexicographic order)
## Remember that age and race group labels, while numeric, indicate nominal
## categories.
age_mix_order <- c("1.2", "2.2",
                   "1.3", "2.3", "3.3",
                   "1.4", "2.4", "3.4", "4.4",
                   "1.5", "2.5", "3.5", "4.5", "5.5")

race_mix_order <- c("1.2", "2.2",
                    "1.3", "2.3", "3.3",
                    "1.4", "2.4", "3.4", "4.4")


################################################################################
               ## MAIN PARTNERSHIPS: FORMULAS AND TARGET STATS ##
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
  mix.age.grp = nodemix_age.grp[-1][match(age_mix_order, names(nodemix_age.grp[-1]))],
  mix.race = nodemix_race[-1][match(race_mix_order, names(nodemix_race[-1]))],
  nodefactor.deg.casl = nodefactor_degcasl[-1],
  nodefactor.diag.status = nodefactor_diagstatus[-1],
  `deg3+` = 0,
  concurrent = concurrent,
  nodematch.diag.status = nodematch_diagstatus,
  nodematch.role.class = c("0" = 0, "1" = 0)
))

netstats_main <- unname(netstats_main)
netstats_main

coef_diss_main <- dissolution_coefs(
  dissolution = ~offset(edges) + offset(nodemix("age.grp")),
  duration = c(
    netstats$netmain$durat_wks,
    with(netstats$netmain,
         durat_wks_byagec[-1][match(age_mix_order, names(durat_wks_byagec[-1]))])
  ),
  d.rate = netstats$demog$mortrate.marginal + mort_correct
)


################################################################################
              ## CASUAL PARTNERSHIPS: FORMULAS AND TARGET STATS ##
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
  mix.age.grp = nodemix_age.grp[-1][match(age_mix_order, names(nodemix_age.grp[-1]))],
  mix.race = nodemix_race[-1][match(race_mix_order, names(nodemix_race[-1]))],
  nodefactor.deg.main = nodefactor_degmain[-1],
  nodefactor.diag.status = nodefactor_diagstatus[-1],
  `deg6+` = 0,
  concurrent = concurrent,
  nodematch.diag.status = nodematch_diagstatus,
  nodematch.role.class = c("0" = 0, "1" = 0)
))

netstats_casl <- unname(netstats_casl)
netstats_casl

coef_diss_casl <- dissolution_coefs(
  dissolution =
    ~offset(edges) + offset(nodemix("age.grp")),
  duration = c(
    netstats$netcasl$durat_wks,
    with(netstats$netcasl,
         durat_wks_byagec[-1][match(age_mix_order, names(durat_wks_byagec[-1]))])
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

netstats_inst <- unname(netstats_inst)
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
nwlabs                  <- c("main", "casl", "inst")
formation_formulas      <- paste0(nwlabs, "_formation")
target_stats            <- paste0("netstats_", nwlabs)
dissolution_coef_fits   <- paste0("coef_diss_", nwlabs)
model_names             <- paste0("fit_", nwlabs)

# Set maximumn number of MCMC iterations
mcmc.maxiterations <- 500

# Fit models in parallel
registerDoParallel(cores = 3)
netest_out <- foreach(i = 1:3, .inorder = FALSE) %dopar% {
  EpiModel::netest(
    nw = nw,
    formation = get(formation_formulas[i]),
    target.stats = get(target_stats[i]),
    coef.diss = get(dissolution_coef_fits[i]),
    set.control.ergm = control.ergm(MCMLE.maxit = mcmc.maxiterations)
  )
}
stopImplicitCluster()

names(netest_out) <- model_names


# %% Write ------------------------------------------------------------------
saveRDS(netest_out, "netest/netest.Rds")
