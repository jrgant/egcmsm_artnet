# %% PLOT SIMULATED NETWORKS ---------------------------------------------------

suppressMessages(library(EpiModelHIV))
ns <- readRDS("netstats/netstats.Rds")

dx <- list(
  readRDS("netest/netdx_main.Rds"),
  readRDS("netest/netdx_casl.Rds"),
  readRDS("netest/netdx_inst.Rds")
)

main <- melt(as.data.table(get_nwstats(dx[[1]])), id.vars = c("time", "sim"))
casl <- melt(as.data.table(get_nwstats(dx[[2]])), id.vars = c("time", "sim"))
# inst <- melt(as.data.table(get_nwstats(dx[[3]])), id.vars = c("time", "sim"))

main_targs <- as.data.table(dx[[1]]$target.stats)
casl_targs <- as.data.table(dx[[2]]$target.stats)
# inst_targs <- as.data.table(dx[[3]]$target.stats)

setnames(main_targs, "names", "variable")
setnames(casl_targs, "names", "variable")
# setnames(inst_targs, "names", "variable")

plotnet <- function(data, targ) {
  data %>%
    ggplot(aes(x = time, y = value)) +
    geom_line(
      aes(group = sim),
      alpha = 0.2,
      color = "firebrick"
    ) +
    geom_hline(
      data = targ,
      aes(yintercept = targets),
      linetype = "dashed",
      size = 1
    ) +
    facet_wrap(
      ~ variable,
      scales = "free_y"
    ) +
    theme_few(base_size = 9)
}

mplot <- plotnet(main, main_targs)
cplot <- plotnet(casl, casl_targs)
# iplot <- plotnet(inst, inst_targs)

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
