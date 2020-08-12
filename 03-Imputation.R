################################################################################
                                  ## SETUP ##
################################################################################

pacman::p_load(
  data.table,
  magrittr,
  ggplot2,
  ggthemes,
  viridis,
  MASS,
  rms,
  Hmisc,
  mice,
  miceadds,
  parallel,
  tictoc
)

an <- fread(paste0(Sys.getenv("ARTNET_PATH"), "/artnet-wide-cleaned.csv"))
anl <- fread(paste0(Sys.getenv("ARTNET_PATH"), "/artnet-long-cleaned.csv"))

# default ggplot theme
theme_set(theme_base())

## NOTE:
## I am using a slightly modified version of the mice() package that I altered
## in order to be able to use a function from the miceadds package in
## parlmice(). Without this change, the clusters could not see miceadds at all.

################################################################################
                           ## SET UP EGO DATA SET ##
################################################################################

ansub <- copy(an)
ansub <- ansub[, .(id, cuml.pnum, prep_revised)]

## Fix NA coding
anl[ego.anal.role == "", ego.anal.role := NA]
anl[p_race.cat == "", p_race.cat := NA]


################################################################################
                             ## UNIVERSAL FIXES ##
################################################################################

## Convert character variables to factor
anl[, ego.anal.role := as.factor(ego.anal.role)]
anl[, ego.race.cat := as.factor(ego.race.cat)]
anl[, p_race.cat := as.factor(p_race.cat)]


################################################################################
                       ## SET UP PARTNERSHIP DATA SETS ##
################################################################################

anlsub <- copy(anl)

## NOTE:
##   PrEP and ART use are coded a little differently than most partnership
##   variables, where the "p_" prefix indicates the alter's characteristic.
##   For PrEP and ART, p_artuse, for instance, indicates the ego's ART use with
##   the given alter, and the "_part" suffix indicates the alter's status.

anlsub <- anlsub[, .(
  id, pid, pid_unique,
  ego.race.cat, ego.age, ego.hiv, ego.anal.role,
  p_race.cat, p_age_i1 = p_age_imputed, abs_sqrt_agediff,
  p_hiv2, p_ongoing_ind, durat_wks, ptype,
  cond.prob, p_rai, p_iai, p_roi, p_ioi,
  p_prepuse = p_prepuse2, p_prepuse_part = p_prepuse_part2,
  oi.rate.52, ai.rate.52,
  p_artuse_bin, p_artuse_part_bin
)]

## Merge selected ego-level variables with the subset of the partnership
## dataset.
mcong_otp <- ansub[
  anlsub[
  (ptype %in% 1:2 & p_ongoing_ind == 1) | ptype == 3],
  on = "id"
]

mcong_otp[, p_ongoing_ind := NULL]

str(mcong_otp)


################################################################################
                          ## MISSING DATA SUMMARIES ##
################################################################################

dropfromviz <- c("id", "pid", "pid_unique")
dropformc <- c(dropfromviz, paste0("p_", c("rai", "iai", "roi", "ioi")))
dropforotp <- c(dropfromviz, "ai.rate.52", "oi.rate.52")

## Main/casual pratnerships
md.pattern(mcong_otp[ptype %in% 1:2, -..dropformc], rotate.names = TRUE)
misscount <- sapply(
  mcong_otp[ptype %in% 1:2, -..dropformc],
  function(x) sum(is.na(x))
)
misscount
round(misscount / mcong_otp[ptype %in% 1:2, .N], 3)

## One-time partnerships
md.pattern(mcong_otp[ptype == 3, -..dropforotp], rotate.names = TRUE)
misscount_otp <- sapply(
  mcong_otp[ptype == 3, -..dropforotp],
  function(x) sum(is.na(x))
)
misscount_otp
round(misscount_otp / mcong_otp[ptype == 3, .N], 3)

## Few complete cases exist, partly due to complex skip patterns.
ncc(mcong_otp)
nic(mcong_otp)

## Calculate proportion of missing data for each variable.
## Unit of analysis: sexual partnerships
missing_summary <- function(data, pship_type) {
  s <- match(
    x = pship_type,

    table = c("main", "casual", "onetime"))

  lapply(data[ptype == s], function(x) {
    data.table(
      misscount = sum(is.na(x)),
      veclength = length(x)
    )[, misspct := misscount / veclength * 100][, misspct := round(misspct, 2)]
  }) %>% rbindlist(idcol = "var")
}

missing_summary(mcong_otp, "main")
missing_summary(mcong_otp, "casual")
missing_summary(mcong_otp, "onetime")

mdp.main <- md.pattern(mcong_otp[ptype == 1], plot = FALSE)
mdp.casl <- md.pattern(mcong_otp[ptype == 2], plot = FALSE)

mdp.main.cnt <- as.numeric(rownames(mdp.main))
mdp.casl.cnt <- as.numeric(rownames(mdp.casl))

# Check if sum of missingness pattern counts equals number of partnerships.
# Should return TRUE.
sum(mdp.main.cnt, na.rm = TRUE) == nrow(mcong_otp[ptype == 1])
sum(mdp.casl.cnt, na.rm = TRUE) == nrow(mcong_otp[ptype == 2])

anymiss_main <- 1 - (mdp.main.cnt[1] / sum(mdp.main.cnt, na.rm = TRUE))
anymiss_casl <- 1 - (mdp.casl.cnt[1] / sum(mdp.casl.cnt, na.rm = TRUE))

round(anymiss_main * 100, 2)
round(anymiss_casl * 100, 2)


################################################################################
                          ## BASE IMPUTATION SETUP ##
################################################################################

## NOTE:
## Imputation approaches and decisions rely heavily on advice provided by Van
## Buuren in the following sources:
##
## 1. van Buuren, S., & Groothuis-Oudshoorn, K. (2011). Mice: multivariate
##    imputation by chained equations in R. J Stat Software, 45(3),
##    1–67. http://dx.doi.org/10.18637/jss.v045.i03
##
## 2. van Buuren, S. (2018). Flexible imputation of missing data, 2nd edition,
##     http://dx.doi.org/10.1201/9780429492259


## FLUX CHECK ------------------------------------------------------------------

anflux <- copy(an)
anflux <- anflux[, grep("part|m_|age|race", names(an)) := NULL]
round(sapply(anflux, function(x) sum(is.na(x))) / nrow(anflux), 3)

anflux_ong <- anflux[mcong_otp, on = "id"]
anflux_ong[, grep("^i\\.", names(anflux_ong)) := NULL]  # drop duplicated vars
anflux_ong
names(anflux_ong)

as.data.table(
  flux(anflux_ong), keep.rownames = TRUE
)[order(outflux), .(rn, pobs, outflux)]

## Select ego-level (i.e., level-2) variables with high outflux to use in
## imputation models. However, STITEST_PERWEEK is also one of the outcomes
## of interest and so should be included a priori in the imputation models.
check_flux <- anflux[, .(
  id,
  an_prep_current, prep_revised,
  pnua_12m, cuml.pnum, stitest_perweek, mmconc
), keyby = id]

check_flux

sapply(check_flux, function(x) sum(is.na(x))) / nrow(check_flux)

## Reassign mcong_otp with additional variables
droplevel1 <- setdiff(names(anflux), names(check_flux))
droplevel1

mcong_otp <- copy(anflux_ong)
mcong_otp[, c(droplevel1) := NULL]

names(mcong_otp)


## SPECIFY PREDICTOR VARIABLES -------------------------------------------------

pred <- make.predictorMatrix(mcong_otp)

completevars <- c(
  "id", "pid", "pid_unique",
  "ego.race.cat", "ego.age", "ego.hiv",
  "p_hiv2", "ptype", "cuml.pnum"
)

## Don't impute complete variables.
pred[completevars, ] <- 0

## Don't impute using ID variables.
pred[, c("id", "pid", "pid_unique")] <- 0

## Identify respondent ID (mice will treat as class/cluster ID)
pred[, "id"] <- -2

## Use passive imputation only for abs_sqrt_agediff
pred["abs_sqrt_agediff", ] <- 0

## Don't use variables that are related determinstically due to skip patterns.
pred["prep_revised", c("ego.hiv", "prep_revised", "p_artuse_bin")] <- 0
pred["an_prep_current", c("prep_revised", "ego.hiv", "p_artuse_bin")] <- 0
pred["p_prepuse", c("prep_revised", "an_prep_current", "ego.hiv")] <- 0
pred["p_prepuse_part", "p_hiv2"] <- 0
pred["p_artuse_bin", c("prep_revised", "an_prep_current", "ego.hiv")] <- 0
pred["p_artuse_part_bin", "p_hiv2"] <- 0


## Prevent ptype-specific variables from imputing outside of relevant ptype.
## For instance, in main and casual partnerships, we want to impute ai.rate.52
## and oi.rate.52 only. For one-time partnerships, we don't want or need these
## rates and are interested in p_rai, p_iai, p_roi, and p_ioi only (indicators)
## of which sex acts took place during the encounter).
sexact.ind <- c("p_rai", "p_iai", "p_roi", "p_ioi")
rate.vars <- c("ai.rate.52", "oi.rate.52")

pred

## OUTPUT PREDICTION MATRIX VIZ ------------------------------------------------
pdt <- as.data.table(pred, keep.rownames = TRUE)
pdt <- melt(pdt, id.vars = "rn")
idvars <- c("id", "pid", "pid_unique")

pdt <- pdt[!rn %in% idvars & !variable %in% idvars]
pdt[, ":=" (
  rn = factor(rn, levels = unique(rn)),
  variable = factor(variable, levels = unique(rn))
)]

dmc <- c("p_rai", "p_iai", "p_roi", "p_ioi")
mc <- pdt[!rn %in% dmc & !variable %in% dmc]

dmo <- c("ego.anal.role", "ptype", "durat_wks", "ai.rate.52", "oi.rate.52")
otp <- pdt[!rn %in% dmo & !variable %in% dmo]

fullyobs <- c(
  "ptype", "p_hiv2", "ego.hiv",
  "ego.age", "ego.race.cat", "cuml.pnum"
)

pmplot <- function(data) {

  ggplot(data, aes(x = variable, y = rn)) +
    geom_point(
      aes(shape = factor(value), fill = factor(value)),
      size = 2
    ) +
    scale_fill_viridis_d(
      name = "Predictor Used to\nImpute Outcome",
      option = "magma",
      begin = 0.5,
      direction = -1
    ) +
    scale_shape_manual(
      name = "Predictor Used to\nImpute Outcome",
      values = c(4, 21)
    ) +
    labs(
      x = "Predictor",
      y = "Outcome"
    ) +
    theme_minimal() +
    theme(
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust = 1
      ))
}


pmp_mc <- pmplot(mc[!rn %in% fullyobs])
pmp_otp <- pmplot(otp[!rn %in% fullyobs])

ggsave(
  "imputation_diagnostics/predmat_plot_maincas.pdf",
  plot = pmp_mc,
  device = "pdf",
  height = 6,
  width = 7.5
)

ggsave(
  "imputation_diagnostics/predmat_plot_onetime.pdf",
  plot = pmp_otp,
  device = "pdf",
  height = 6,
  width = 7.5
)



## SPECIFY IMPUTATION METHODS --------------------------------------------------

meth <- make.method(mcong_otp)
length(meth)

## Level 1 variables to impute
l1_cat <- c(
  "ego.anal.role",
  "p_race.cat"
)

l1_bin <- c(
  "prep_revised",
  "p_prepuse",
  "p_prepuse_part",
  "p_artuse_bin",
  "p_artuse_part_bin",
  sexact.ind
)

l1_count <- c(
  "cond.prob",
  "p_age_i1",
  "durat_wks",
  "ai.rate.52",
  "oi.rate.52"
)

## Level 2 variables to impute
l2 <- c(
  "an_prep_current",
  "prep_revised",
  "stitest_perweek",
  "pnua_12m",
  "mmconc"
)

## Methods
meth["abs_sqrt_agediff"] <- "~ I(abs(sqrt(ego.age) - sqrt(p_age_i1)))"
meth[c(l1_cat, l1_count)] <- "2l.pmm"
meth[l1_bin] <- "2l.bin"
meth[l2] <- "2lonly.pmm"

meth

sapply(mcong_otp, class)

## Fix classes for a couple of variables.
fixclass <- c("p_hiv2", "ego.hiv")
mcong_otp[, (fixclass) := lapply(.SD, as.factor), .SDcols = fixclass]
mcong_otp[, oi.rate.52 := as.integer(oi.rate.52)]

## NOTE: The conversion of oi.rate.52 above throws a warning about introducing
##       missing values. Values appear to match their originals,
##       so considering this warning to be benign.


################################################################################
                ## SPLIT IMPUTATION TASKS BY PARTNERSHIP TYPE ##
################################################################################

## Make two subsets of the data, predictors, and methods for imputing values by
## partnerships type. Main/casual in one group and one-time contacts in the
## other.


# MAIN/CASUAL PARTNERSHIPS -----------------------------------------------------

mcong <- mcong_otp[ptype %in% 1:2][, (sexact.ind) := NULL]

predmc <- pred[
  -grep(paste0(sexact.ind, collapse = "|"), rownames(pred)),
  -grep(paste0(sexact.ind, collapse = "|"), colnames(pred))
  ]

methmc <- meth[!names(meth) %in% sexact.ind]

## Initialize the main/casual MICE object.
inimc <- mice(
  mcong,
  maxit = 0,
  predictorMatrix = predmc,
  method = methmc
)

## Check visit sequence to make sure that abs_sqrt_agediff is passively imputed
## right after p_age_i1 is imputed. The variable subsetting conducted at the
## beginning of this script ordered these variables intentionally.
vismc <- inimc$visitSequence
vismc


# ONE-TIME CONTACTS -----------------------------------------------------------

dropfromotp <- c(rate.vars, "ego.anal.role", "durat_wks", "ptype")

otp <- mcong_otp[ptype == 3][, (dropfromotp) := NULL]

predotp <- pred[
  -grep(paste0(dropfromotp, collapse = "|"), rownames(pred)),
  -grep(paste0(dropfromotp, collapse = "|"), colnames(pred))
]

methotp <- meth[!names(meth) %in% dropfromotp]

iniotp <- mice(
  otp,
  maxit = 0,
  predictorMatrix = predotp,
  method = methotp
)

## Check visit sequence to make sure that abs_sqrt_agediff is passively imputed
## right after p_age_i1 is imputed. The variable subsetting conducted at the
## beginning of this script ordered these variables intentionally.
visotp <- iniotp$visitSequence
visotp

################################################################################
                 ## POST-PROCESSING (ALL PARTNERSHIP TYPES) ##
################################################################################

## Post-process PrEP and ART use variables to avoid violating skip patterns
## enforced in the ARTnet questionnaire.

## Skip pattern summary:
##
## PrEP use within partnership:
##
##   + Only HIV-NEGATIVE EGOS reported on their own PrEP use during a
##     partnership (P_PREPUSE) or ever (PREP_REVISE).
##
##   + Only HIV-NEGATIVE or -UNKNOWN PARTNERS had information regarding their
##     PrEP use during a partnership.
##
## ART use within partnership:
##
##   + Only HIV-POSITIVE EGOS reported their own ART use during a partnership.
##
##   + Only HIV-POSITIVE PARTNERS had information regarding their ART use during
##     a partnership.
##
## Ever PrEP use:
##
##   + Respondents were asked about ever taking PrEP only if they reported
##     being HIV-negative.
##
## Current PrEP use:
##
##   + Respondents were asked about current PrEP use only if they reported ever
##     having taken PrEP.

postmc <- inimc$post
postotp <- iniotp$post

## Formulas.
post_p_prepuse <- "imp[[j]][data$prep_revised[!r[, j]] != 1, i] <- 0"
post_p_artuse_bin <- "imp[[j]][data$ego.hiv[!r[, j]] %in% c(0, 2), i] <- 0"

post_p_prepuse_part <- "imp[[j]][data$p_hiv2[!r[, j]] == 1, i] <- 0"
post_p_artuse_part_bin <- "imp[[j]][data$p_hiv2[!r[, j]] %in% c(0, 2), i] <- 0"

post_prep_revised <- "imp[[j]][data$ego.hiv[!r[, j]] == 1, i] <- 0"
post_an_prep_current <- "imp[[j]][data$prep_revised[!r[, j]] == 0, i] <- 0"

## Assignment.
postmc["p_prepuse"] <- post_p_prepuse
postotp["p_prepuse"] <- post_p_prepuse

postmc["p_artuse_bin"] <- post_p_artuse_bin
postotp["p_artuse_bin"] <- post_p_artuse_bin

postmc["p_prepuse_part"] <- post_p_prepuse_part
postotp["p_prepuse_part"] <- post_p_prepuse_part

postmc["p_artuse_part_bin"] <- post_p_artuse_part_bin
postotp["p_artuse_part_bin"] <- post_p_artuse_part_bin

postmc["prep_revised"] <- post_prep_revised
postotp["prep_revised"] <- post_prep_revised

postmc["an_prep_current"] <- post_an_prep_current
postotp["an_prep_current"] <- post_an_prep_current

postmc
postotp


################################################################################
                     ## IMPUTE MAIN/CASUAL PARTNERSHIPS ##
################################################################################

## NOTE: Donor pool size selected based on:

##       Morris TP, White IR, Royston P. Tuning multiple imputation by
##       predictive mean matching and local residual draws. BMC Med. Res.
##       Methodol. 2014;14:75. http://dx.doi.org/10.1186/1471-2288-14-75

## Make sure lengths are all equal. Should return TRUE.
all.equal(
  length(names(mcong)),
  length(colnames(predmc)),
  length(methmc),
  length(vismc),
  length(postmc),
  tolerance = 0
)

check.varlist <- as.data.table(cbind(
  sort(names(mcong)),
  sort(colnames(predmc)),
  sort(names(methmc)),
  sort(vismc),
  sort(names(postmc))
))

check.varlist[, .N, names(check.varlist)][, check := N == 1][]

# FIXME
# Getting the "glmer does not run" error messages.
# See org TODO about this issue.

imp_mc <- parlmice(
  mcong,
  maxit = 100,
  n.core = 5,
  m = 4,
  predictorMatrix = predmc,
  method = methmc,
  visitSequence = vismc,
  post = postmc,
  cluster.seed = 45345,
  donors = 10L,
  # printFlag = TRUE,
  cl.loads = TRUE
)

imp_mc$loggedEvents


################################################################################
                      ## WRITE MAIN/CASUAL IMPUTATIONS ##
################################################################################
saveRDS(imp_mc, file.path(Sys.getenv("ARTNET_PATH"), "artnet-imputed-mc.Rds"))


################################################################################
                     ## DIAGNOSE MAIN/CASUAL IMPUTATIONS ##
################################################################################

imp_mc <- readRDS(
  file.path(Sys.getenv("ARTNET_PATH"), "artnet-imputed-mc.Rds")
)

## Check for impossible combinations (violations of skip patterns from ARTnet)
cimp_mc <- as.data.table(complete(imp_mc, "long"))

### PrEP, ego
mcong[, .N, keyby = .(ego.hiv, prep_revised)]   # original data
cimp_mc[, .N, keyby = .(ego.hiv, prep_revised)]     # imp_mcuted data

mcong[, .N, keyby = .(ego.hiv, p_prepuse)]      # original data
cimp_mc[, .N, keyby = .(ego.hiv, p_prepuse)]        # imp_mcuted data

mcong[, .N, keyby = .(an_prep_current, prep_revised)]  # original data
cimp_mc[, .N, keyby = .(an_prep_current, prep_revised)]    # imp_mcuted data


### PrEP, alter
mcong[, .N, keyby = .(p_hiv2, p_prepuse_part)]  # original data
cimp_mc[, .N, keyby = .(p_hiv2, p_prepuse_part)]    # imp_mcuted data

### ART, ego
mcong[, .N, keyby = .(ego.hiv, p_artuse_bin)]      # original data
cimp_mc[, .N, keyby = .(ego.hiv, p_artuse_bin)]        # imp_mcuted data

### ART, alter
mcong[, .N, keyby = .(p_hiv2, p_artuse_part_bin)]  # original data
cimp_mc[, .N, keyby = .(p_hiv2, p_artuse_part_bin)]    # imp_mcuted data


imp_dx_dir <- "imputation_diagnostics"

### Trace plots
pdf(file.path(imp_dx_dir, paste0("imp_mc_trace_", Sys.Date(), ".pdf")), onefile = TRUE)
plot(imp_mc)
dev.off()

### Density plots
pdf(file.path(imp_dx_dir, paste0("imp_mc_density_", Sys.Date(), ".pdf")), onefile = TRUE)
densityplot(
  imp_mc,
  thicker = 3,
  scales = list(
    x = list(relation = "free"),
    y = list(relation = "free")
  )
)
dev.off()

### closer look at oral and anal act rates
## pdf(file.path(imp_mcdx_dir, "imp_mc_density_actrates-only.pdf"))
## densityplot(imp_mc, ~ oi.rate.52 + ai.rate.52)
## dev.off()


################################################################################
                       ## IMPUTE ONE-TIME PARTNERSHIPS ##
################################################################################


## Make sure lengths are all equal. Should return TRUE.
all.equal(
  length(names(otp)),
  length(colnames(predotp)),
  length(methotp),
  length(visotp),
  length(postotp),
  tolerance = 0
)

check.varlist.otp <- as.data.table(cbind(
  sort(names(otp)),
  sort(colnames(predotp)),
  sort(names(methotp)),
  sort(visotp),
  sort(names(postotp))
))

check.varlist.otp[, .N, names(check.varlist.otp)][, check := N == 1][]

# FIXME
# Getting the "glmer does not run" error messages.
# See org TODO about this issue.

imp_otp <- parlmice(
  otp,
  maxit = 50,
  n.core = 5,
  m = 4,
  predictorMatrix = predotp,
  method = methotp,
  visitSequence = visotp,
  post = postotp,
  cluster.seed = 45345,
  donors = 10L,
  # printFlag = TRUE,
  cl.loads = TRUE
)

imp_otp$loggedEvents


################################################################################
                   ## WRITE ONE-TIME CONTACTS IMPUTATIONS ##
################################################################################
saveRDS(imp_otp, file.path(Sys.getenv("ARTNET_PATH"), "artnet-imputed-otp.Rds"))


################################################################################
                     ## DIAGNOSE MAIN/CASUAL IMPUTATIONS ##
################################################################################

imp_otp <- readRDS(
  file.path(Sys.getenv("ARTNET_PATH"), "artnet-imputed-otp.Rds")
)

## Check for impossible combinations (violations of skip patterns from ARTnet)
cimp_otp <- as.data.table(complete(imp_otp, "long"))

### PrEP, ego
otp[, .N, keyby = .(ego.hiv, prep_revised)]       # original data
cimp_otp[, .N, keyby = .(ego.hiv, prep_revised)]  # imp_otputed data

otp[, .N, keyby = .(ego.hiv, p_prepuse)]          # original data
cimp_otp[, .N, keyby = .(ego.hiv, p_prepuse)]     # imp_otputed data

otp[, .N, keyby = .(an_prep_current, prep_revised)]       # original data
cimp_otp[, .N, keyby = .(an_prep_current, prep_revised)]  # imp_otputed data

### PrEP, alter
otp[, .N, keyby = .(p_hiv2, p_prepuse_part)]       # original data
cimp_otp[, .N, keyby = .(p_hiv2, p_prepuse_part)]  # imp_otputed data

### ART, ego
otp[, .N, keyby = .(ego.hiv, p_artuse_bin)]        # original data
cimp_otp[, .N, keyby = .(ego.hiv, p_artuse_bin)]   # imp_otputed data

### ART, alter
otp[, .N, keyby = .(p_hiv2, p_artuse_part_bin)]       # original data
cimp_otp[, .N, keyby = .(p_hiv2, p_artuse_part_bin)]  # imp_otputed data

imp_dx_dir <- "imputation_diagnostics"

### Trace plots
pdf(file.path(imp_dx_dir, paste0("imp_otp_trace_", Sys.Date(), ".pdf")), onefile = TRUE)
plot(imp_otp)
dev.off()

### Density plots
pdf(file.path(imp_dx_dir, paste0("imp_otp_density_", Sys.Date(), ".pdf")), onefile = TRUE)
densityplot(
  imp_otp,
  thicker = 3,
  scales = list(
    x = list(relation = "free"),
    y = list(relation = "free")
  )
)
dev.off()

### closer look at oral and anal act rates
## pdf(file.path(imp_otpdx_dir, "imp_otp_density_actrates-only.pdf"))
## densityplot(imp_otp, ~ oi.rate.52 + ai.rate.52)
## dev.off()


################################################################################
                              ## TEST ANALYSIS ##
################################################################################

## NOTE: Remove this section when done.

library(data.table)
library(MASS)
library(mice)
library(ggplot2)
library(magrittr)

## NOTE:
## Reference for variance estimation of model predictions using the "combine
## then predict" approach:
##
## Miles A. Obtaining Predictions from Models Fit to Multiply Imputed Data.
## Sociol. Methods Res. 2016;45(1):175–185.
## https://doi.org/10.1177/0049124115610345

imp2 <- readRDS(
  file.path(Sys.getenv("ARTNET_PATH"), "artnet-maincas-imputed.rds")
)

test.ai <- with(imp2, glm.nb(
  ai.rate.52 ~
    ptype +
    ego.race.cat + p_race.cat +
    ego.age + p_age_i1 +
    abs_sqrt_agediff
))

test.ai.pooled <- pool(test.ai)

sapply(1:20, function(x) test.ai$analyses[[x]]$theta)

test.ai.cc <- glm.nb(
  ai.rate.52 ~
    ptype + ego.race.cat +
    p_race.cat + ego.age +
    p_age_i1 + abs_sqrt_agediff,
  data = mcong_otp
)

pool(with(imp2, lm(p_age_i1 ~ 1))) # quick test to get pooled means with standard errors

cbind(test.ai.pooled$pooled$estimate, coef(test.ai.cc))

lapply(test.ai$analyses, predict, type = "response")

class(test.ai.pooled)
summary(test.ai.pooled)
predict(test.ai.pooled)
