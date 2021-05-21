# %% SETUP ---------------------------------------------------------------------

# NOTE: If not resimulating networks, go to last section to view plots of
#       previously simulations.

suppressMessages(library(EpiModel))
est <- readRDS("netest/netest.Rds")

use_ncores <- as.numeric(Sys.getenv("SLURM_NPROCS"))
n_networks <- 50
n_timesteps <- 52 * 10


# %% MAIN PARTNERSHIPS ---------------------------------------------------------

main_formation_full <-
  ~ edges +
    nodemix("age.grp", levels = NULL) +
    nodemix("race", levels = NULL) +
    ## nodefactor("age.grp", levels = NULL) +
    ## nodefactor("race", levels = NULL) +
    nodefactor("deg.casl", levels = NULL) +
    nodefactor("diag.status", levels = NULL) +
    degrange(from = 3) +
    concurrent +
    ## nodematch("race") +
    ## nodematch("age.grp") +
    nodematch("diag.status") +
    nodematch("role.class", levels = c(1, 2))

dx_main <- netdx(
  est[[1]],
  nsims = n_networks,
  nsteps = n_timesteps,
  ncores = use_ncores,
  skip.dissolution = FALSE,
  nwstats.formula = main_formation_full,
)

saveRDS(dx_main, "netest/netdx_main.Rds")
dx_main <- NULL
gc()

# %% CASUAL PARTNERSHIPS -------------------------------------------------------

casl_formation_full <-
  ~ edges +
    nodemix("age.grp", levels = NULL) +
    nodemix("race", levels = NULL) +
    ## nodefactor("age.grp", levels = NULL) +
    ## nodefactor("race", levels = NULL) +
    nodefactor("deg.main", levels = NULL) +
    nodefactor("diag.status", levels = NULL) +
    degrange(from = 6) +
    concurrent +
    ## nodematch("race", levels = NULL) +
    ## nodematch("age.grp") +
    nodematch("diag.status") +
    nodematch("role.class", levels = c(1, 2))

dx_casl <- netdx(
  est[[2]],
  nsims = n_networks,
  nsteps = n_timesteps,
  ncores = use_ncores,
  skip.dissolution = FALSE,
  nwstats.formula = casl_formation_full
)

saveRDS(dx_casl, "netest/netdx_casl.Rds")
dx_casl <- NULL
gc()


# %% ONE-TIME CONTACTS ---------------------------------------------------------

inst_formation_full <-
  ~ edges +
    nodefactor("race", levels = NULL) +
    nodefactor("age.grp", levels = NULL) +
    nodefactor("deg.main", levels = NULL) +
    nodefactor("deg.casl", levels = NULL) +
    nodefactor("diag.status", levels = NULL) +
    nodematch("role.class", levels = c(1, 2))

dx_inst <- netdx(
  est[[3]],
  nsims = 10000,
  nsteps = n_timesteps,
  ncores = use_ncores,
  skip.dissolution = FALSE,
  nwstats.formula = inst_formation_full,
  dynamic = FALSE
)

saveRDS(dx_inst, "netest/netdx_inst.Rds")
dx_inst <- NULL
gc()
