# %% SETUP ---------------------------------------------------------------------

suppressMessages(library(EpiModel))

mdx <- readRDS("netest/dx_main.Rds")
cdx <- readRDS("netest/dx_casl.Rds")


# %% PLOT SIMULATED NETWORKS ---------------------------------------------------

# Quick Diagnostics

plot(mdx, qnts = 0.95)
plot(cdx, qnts = 0.95)

print(mdx)
print(cdx)
