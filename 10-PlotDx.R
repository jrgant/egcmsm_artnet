# %% PLOT SIMULATED NETWORKS ---------------------------------------------------

suppressMessages(library(EpiModelHIV))
ns <- readRDS("netstats/netstats.Rds")

dx <- list(
  readRDS("netest/netdx_main.Rds"),
  readRDS("netest/netdx_casl.Rds"),
  readRDS("netest/netdx_inst.Rds")
)

plotdx <- function(network, color = "skyblue") {
  nt <- c("main", "casual", "inst")
  plot(dx[[grep(network, nt)]], qnts = 0.95, sim.col = color)
}

plotdx("main")
plotdx("casual")
plotdx("inst")

plot(dx[[1]], sim.lines = TRUE)
plot(dx[[2]], sim.lines = TRUE)
plot(dx[[3]], sim.lines = TRUE)
