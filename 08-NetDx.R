suppressMessages(library(EpiModel))


mdx <- readRDS("netest/dx_main.Rds")
cpx <- readRDS("netest/dx_casl.Rds")

plot(mdx, qnts = 0.95)
plot(cpx, qnts = 0.95)
