# %% SETUP ---------------------------------------------------------------------

suppressMessages(library(EpiModel))

dx <- readRDS("netest/netdx.Rds")

# %% PLOT SIMULATED NETWORKS ---------------------------------------------------

plotdx <- function(network, color = "skyblue") {
  nt <- c("main", "casual", "inst")
  plot(dx[[grep(network, nt)]], qnts = 0.95, sim.col = color)
}

plotdx("main")
plotdx("casual")
plotdx("inst")


# Get network summaries

print(mdx)
print(cdx)
print(idx)
quickdag::qd_swig()
