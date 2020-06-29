# %% SETUP ---------------------------------------------------------------------

suppressMessages(library(EpiModel))

dx <- readRDS("netest/netdx.Rds")

# %% PLOT SIMULATED NETWORKS ---------------------------------------------------

# Quick Diagnostics
plot(dx[[1]], qnts = 0.95, sim.col = "cornflowerblue")
plot(cdx, qnts = 0.95, sim.col = "cornflowerblue")
plot(idx, qnts = 0.95, sim.col = "cornflowerblue")


print(mdx)
print(cdx)
print(idx)
