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

plot(dx[[1]], "duration")

# Get network summaries
dx[[1]]
dx[[2]]
dx[[3]]
