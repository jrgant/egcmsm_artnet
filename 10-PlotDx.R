# %% PLOT SIMULATED NETWORKS ---------------------------------------------------

suppressMessages(library(EpiModel))
library(data.table)
library(ggplot2)
library(ggthemes)
library(stringr)

# Fitted network models.
ne <- readRDS("netest/netest.Rds")
ns <- readRDS("netstats/netstats.Rds")

# Simulated diagnostics.
dx <- list(
  main = readRDS("netest/netdx_main.Rds"),
  casl = readRDS("netest/netdx_casl.Rds"),
  inst = readRDS("netest/netdx_inst.Rds")
)

# Plot simulated networks (EpiModel default visualization)
plot(dx$main, qnts = 0.95)
plot(dx$casl, qnts = 0.95)
plot(dx$inst, sim.lines = TRUE)

# Function: Bring in original network statistics as a double-check.
organize_targets <- function(netstats) {
  lapply(names(netstats), \(x) {
    dt <- as.data.table(netstats[[x]], keep.rownames = TRUE)
    dt[, term := x]
    if (x %in% c("edges", "concurrent", "durat_wks_byage", "durat_wks") |
          x %like% "nodematch") {
      setnames(dt, "V1", "specified_target")
      dt[, names := x]
    } else {
      setnames(dt, c("V1", "V2"), c("group", "specified_target"))
      dt[, names := paste0(str_extract(x, "nodefactor|mix"), ".",
                           str_extract(x, "race|age.grp|degmain|degcasl|diagstatus"), ".",
                           group)]
      if (x %like% "byagec") {
        dt[, names := paste0("mix.age.grp", group)]
      }
    }
    dt[, .(names, specified_target)]
  }) |> rbindlist()
}

# Function: Retrieve the simulated means, as well as the target statistics
#           used to fit the network models.
extract_sim_stats <- function(dxobj) {
  dt <- as.data.table(dxobj$stats.table.formation, keep.rownames = TRUE)
  dt[, .(names = rn, fitted_target = Target, simulated_target = `Sim Mean`,
         pctdiff = `Pct Diff`)]
}

# Format the original target objects to match the diagnosis output.
nsmain_targets <- organize_targets(ns$netmain)
nscasl_targets <- organize_targets(ns$netcasl)

# Extract the simulated statistics.
main_simstats <- extract_sim_stats(dx$main)
casl_smistats <- extract_sim_stats(dx$casl)

# Merge the original specified targets with the model output
main_dx_merged <- merge(main_simstats, nsmain_targets, by = "names")
casl_dx_merged <- merge(casl_smistats, nscasl_targets, by = "names")

plot_merged <- function(data, title) {
  data |>
    ggplot(aes(x = names, y = pctdiff)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_segment(aes(xend = names, yend = 0)) +
    geom_point(shape = 21, fill = "white", size = 3) +
    labs(x = "Names", y = "Percent Difference") +
    ggtitle(title) +
    coord_flip() +
    theme_pander(base_size = 16) +
    theme()
}

plot_merged(main_dx_merged, "Main Partnership Network")
plot_merged(casl_dx_merged, "Casual Partnership Network")
