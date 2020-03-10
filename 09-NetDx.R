# %% SETUP ---------------------------------------------------------------------

suppressMessages(library(EpiModel))

mdx <- readRDS("netest/dx_main.Rds")
cdx <- readRDS("netest/dx_casl.Rds")
idx <- readRDS("netest/dx_inst.Rds")



# %% PLOT SIMULATED NETWORKS ---------------------------------------------------

# Quick Diagnostics
plot(mdx, qnts = 0.95, sim.col = "cornflowerblue")
plot(cdx, qnts = 0.95, sim.col = "cornflowerblue")
plot(idx, qnts = 0.95, sim.col = "cornflowerblue")


print(mdx)
print(cdx)
print(idx)
