# %% SETUP ---------------------------------------------------------------------

suppressMessages(library(EpiModel))
est <- readRDS("netest/netest.Rds")

if (Sys.getenv("SLURM_NPROCS") == "") {
  ncores_avail <- 12
} else {
  ncores_avail <- as.numeric(Sys.getenv("SLURM_NPROCS"))
}

## Wrapper around netdx.
run_netdx <- function(nw_est, n_networks = 50, n_timesteps = 52 * 60,
                      use_ncores = ncores_avail, formation_fml,
                      dynset = TRUE) {
  netdx(
    nw_est,
    nsims = n_networks,
    nsteps = n_timesteps,
    ncores = use_ncores,
    dynamic = dynset,
    sequential = FALSE,
    verbose = FALSE,
    skip.dissolution = FALSE,
    nwstats.formula = formation_fml,
  )
} 


# %% MAIN PARTNERSHIPS ---------------------------------------------------------

main_formation_full <- ~ edges +
  nodemix("age.grp", levels = NULL) +
  nodemix("race", levels = NULL) +
  nodefactor("deg.casl", levels = NULL) +
  nodefactor("diag.status", levels = NULL) +
  degrange(from = 3) +
  concurrent +
  nodematch("diag.status") +
  nodematch("role.class", diff = TRUE, levels = I(0:1))

dx_main <- run_netdx(
  nw_est = est[[1]],
  formation_fml = main_formation_full
)

saveRDS(dx_main, "netest/netdx_main.Rds")
dx_main <- NULL
gc()


# %% CASUAL PARTNERSHIPS -------------------------------------------------------

casl_formation_full <- ~ edges +
  nodemix("age.grp", levels = NULL) +
  nodemix("race", levels = NULL) +
  nodefactor("deg.main", levels = NULL) +
  nodefactor("diag.status", levels = NULL) +
  degrange(from = 6) +
  concurrent +
  nodematch("diag.status") +
  nodematch("role.class", levels = I(0:1))

dx_casl <- run_netdx(
  nw_est = est[[2]],
  formation_fml = casl_formation_full
)

saveRDS(dx_casl, "netest/netdx_casl.Rds")
dx_casl <- NULL
gc()


# %% ONE-TIME CONTACTS ---------------------------------------------------------

inst_formation_full <- ~ edges +
  nodefactor("race", levels = NULL) +
  nodefactor("age.grp", levels = NULL) +
  nodefactor("deg.main", levels = NULL) +
  nodefactor("deg.casl", levels = NULL) +
  nodefactor("diag.status", levels = NULL) +
  nodematch("role.class", levels = I(0:1))

dx_inst <- run_netdx(
  nw_est = est[[3]],
  formation_fml = inst_formation_full,
  dynset = FALSE
)

saveRDS(dx_inst, "netest/netdx_inst.Rds")
dx_inst <- NULL
gc()
