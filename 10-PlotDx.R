# %% PLOT SIMULATED NETWORKS ---------------------------------------------------

suppressMessages(library(EpiModel))
library(data.table)
library(ggplot2)
library(ggthemes)
library(mgcv)

# Fitted network models.
ne <- readRDS("netest/netest.Rds")

# Simulated diagnostics.
dx <- list(
  main = readRDS("netest/netdx_main.Rds"),
  casl = readRDS("netest/netdx_casl.Rds"),
  inst = readRDS("netest/netdx_inst.Rds")
)

main_dx <- melt(as.data.table(get_nwstats(dx$main)), id.vars = c("time", "sim"))
casl_dx <- melt(as.data.table(get_nwstats(dx$casl)), id.vars = c("time", "sim"))
## inst_dx <- melt(as.data.table(get_nwstats(dx[[3]])), id.vars = c("time", "sim"))

# Simulated means.
## main_gam <- mgcv::bam(value ~ s(time, bs = "cs", by = variable),
##                       family = gaussian(),
##                       data = main_dx)

main_dx_means <- unique(main_dx[, .(variable, time)])
main_dx_means[, pred := predict(main_gam, newdata = main_dx_means)]

# Read in target statistics.
main_targs <- as.data.table(ne$fit_main$target.stats, keep.rownames = TRUE)
casl_targs <- as.data.table(ne$fit_casl$target.stats, keep.rownames = TRUE)
## inst_targs <- as.data.table(ne$fit_inst$target.stats, keep.rownames = TRUE)

targ_varnames <- c("variable", "value")
setnames(main_targs, new = targ_varnames)
setnames(casl_targs, new = targ_varnames)
## setnames(inst_targs, new = targ_varnames)

# Function: Plot estimated network statistics against the corresponding target
#           statistics
plotnet <- function(data_dx, data_targ, time_start = 1, sim_lines = TRUE) {

  plot <- data_dx[time >= time_start] |>
    ggplot(aes(x = time, y = value))

  if (sim_lines == TRUE) {
    plot <- plot +
      geom_line(
        aes(group = sim),
        alpha = 0.2,
        color = "firebrick"
      )
  } else {
    plot <- plot +
      geom_smooth(aes(x = time, y = value), method = "loess", se = FALSE)
  }
  plot +
    geom_hline(
      data = data_targ,
      aes(yintercept = value),
      linetype = "dashed",
      size = 1
    ) +
    facet_wrap(
      ~ variable,
      scales = "free_y"
    ) +
    ggtitle("Diagnostic Statistics",
            subtitle = paste("Time steps", time_start, "--3120")) +
    theme_few(base_size = 9)
}

mplot <- plotnet(main_dx, main_targs, time_start = 2860)
cplot <- plotnet(casl_dx, casl_targs, time_start = 2860)
# iplot <- plotnet(inst_dx, inst_targs)

ggsave(
  "plotdx_main.png",
  mplot,
  width = 10,
  height = 7
)

ggsave(
  "plotdx_casl.png",
  cplot,
  width = 10,
  height = 7
)
