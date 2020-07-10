################################################################################
                                  ## SETUP ##
################################################################################

pacman::p_load(
  data.table,
  summarytools,
  ciTools,
  dplyr,
  stringr,
  lubridate,
  ggplot2,
  ggthemes,
  viridis,
  MASS,
  rms,
  Hmisc,
  mgcv,
  survey,
  mice,
  miceadds
)

an <- fread(paste0(Sys.getenv("ARTNET_PATH"), "/artnet-wide-cleaned.csv"))
anl <- fread(paste0(Sys.getenv("ARTNET_PATH"), "/artnet-long-cleaned.csv"))

# default ggplot theme
theme_set(theme_base())


################################################################################
                           ## SET UP EGO DATA SET ##
################################################################################

ansub <- copy(an)
ansub <- ansub[, .(id, cuml.pnum, prep_revised)]


################################################################################
                             ## UNIVERSAL FIXES ##
################################################################################

## Convert character variables to factor
anl[, ego.anal.role := as.factor(ego.anal.role)]
anl[, ego.race.cat := as.factor(ego.race.cat)]
anl[, p_race.cat := as.factor(p_race.cat)]

## Fix NA coding
anl[ego.anal.role == "", ego.anal.role := NA]
anl[mc_ong[p_race.cat == "", p_race.cat := NA]


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
  ego.race.cat, ego.age, ego.hiv, ego.ongoing, ego.anal.role,
  p_race.cat, p_age_i1 = p_age_imputed, p_hiv2, p_ongoing_ind,
  durat_wks, ptype, ai.rate, oi.rate, cond.prob,
  p_prepuse, p_prepuse_part, p_artuse, p_artuse_part
)]

# create main/casual partnership data set and merge in level 2 (ego) variables
mc_ong <- ansub[anlsub[ptype != 3 & p_ongoing_ind == 1], on = "id"]

dropcols <- c("ego.ongoing", "p_ongoing_ind")
mc_ong <- mc_ong[, (dropcols) := NULL]
mc_ong

md.pattern(mc_ong, rotate.names = TRUE)
ncc(mc_ong)
nic(mc_ong)

mc_ong[, ai.rate.52 := floor(ai.rate * 52)]
mc_ong[, oi.rate.52 := floor(oi.rate * 52)]

## View and update PrEP use within partnership
mc_ong[, .N, keyby = p_prepuse]
mc_ong[, .N, keyby = p_prepuse_part]

## RECODE EGO PREP USE

## If respondent reported never being on PrEP, set ego PrEP use within
## partnership to "Never".
mc_ong[prep_revised == 0 | p_prepuse == 3, p_prepuse := 0]
mc_ong[p_prepuse %in% c(1, 2), p_prepuse := 1]  # any PrEP use in relationship

## RECODE ALTER PREP USE
mc_ong[p_prepuse_part %in% c(1, 2), p_prepuse_part := 1]
mc_ong[p_prepuse_part == 3, p_prepuse_part := 0]

mc_ong[, .N, keyby = p_prepuse]
mc_ong[, .N, keyby = p_prepuse_part]

## RECODE EGO ART USE
mc_ong[, .N, .(ego.hiv, p_artuse)]

## any ego ART use with given partner
mc_ong[ego.hiv == 1, p_artuse_bin := ifelse(
  p_artuse %in% 1:2, 1, p_artuse
)]

mc_ong[p_artuse_bin == 3, p_artuse_bin := 0]

mc_ong[, .N, keyby = .(p_artuse, p_artuse_bin)]

## RECODE ALTER ART USE
mc_ong[, .N, keyby = .(p_artuse_part, p_hiv2)]

mc_ong[p_hiv2 == 1, p_artuse_part_bin := ifelse(
  p_artuse_part %in% 1:2, 1, p_artuse_part
)]

mc_ong[p_artuse_part_bin == 3, p_artuse_part_bin := 0]

mc_ong[, .N, keyby = .(p_artuse_part, p_artuse_part_bin)]

mc_ong[, c("p_artuse", "p_artuse_part") := NULL]

## CALCULATE DIFFERENCE IN AGE SQUARE ROOT DIFFERENCES
mc_ong[, abs_sqrt_agediff := abs(sqrt(ego.age) - sqrt(p_age_i1))]

weekly_rates <- c("ai.rate", "oi.rate")
mc_ong[, (weekly_rates) := NULL]

mc_ong[, .N, keyby = ego.hiv]
mc_ong[, .N, keyby = p_hiv2]

sapply(mc_ong, class)


################################################################################
                          ## MISSING DATA SUMMARIES ##
################################################################################

## Calculate proportion of missing data for each variable
## Unit of analysis, sexual partnerships
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

missing_summary(mc_ong, "main")
missing_summary(mc_ong, "casual")

mdp.main <- md.pattern(mc_ong[ptype == 1], plot = FALSE)
mdp.casl <- md.pattern(mc_ong[ptype == 2], plot = FALSE)

mdp.main.cnt <- as.numeric(rownames(mdp.main))
mdp.casl.cnt <- as.numeric(rownames(mdp.casl))

# Check if sum of missingness pattern counts equals number of partnerships.
# Should return TRUE.
sum(mdp.main.cnt, na.rm = TRUE) == nrow(mc_ong[ptype == 1])
sum(mdp.casl.cnt, na.rm = TRUE) == nrow(mc_ong[ptype == 2])

anymiss_main <- 1 - (mdp.main.cnt[1] / sum(mdp.main.cnt, na.rm = TRUE))
anymiss_casl <- 1 - (mdp.casl.cnt[1] / sum(mdp.casl.cnt, na.rm = TRUE))

round(anymiss_main * 100, 2)
round(anymiss_casl * 100, 2)


################################################################################
                           ## IMPUTE MISSING DATA ##
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

ncc(mc_ong) # no complete cases register, due to complex skip patterns

anflux <- copy(an)
anflux <- anflux[, grep("part|m_|age|race", names(an)) := NULL]
round(sapply(anflux, function(x) sum(is.na(x))) / nrow(anflux), 3)

anflux_ong <- anflux[mc_ong, on = "id"]
anflux_ong[, grep("^i\\.", names(anflux_ong)) := NULL]  # drop duplicated vars
anflux_ong
names(anflux_ong)

as.data.table(flux(anflux_ong), keep.rownames = TRUE)[order(outflux)]

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

## Reassign mc_ong with additional variables
droplevel1 <- setdiff(names(anflux), names(check_flux))
droplevel1

mc_ong <- copy(anflux_ong)
mc_ong[, c(droplevel1) := NULL]

names(mc_ong)


## SPECIFY PREDICTOR VARIABLES -------------------------------------------------

pred <- make.predictorMatrix(mc_ong)

completevars <- c(
  "id", "pid", "pid_unique",
  "ego.race.cat", "ego.age", "ego.hiv",
  "p_hiv2", "ptype"
)

## Don't impute complete variables
pred[completevars, ] <- 0

## Don't impute using ID variables
pred[, c("id", "pid", "pid_unique")] <- 0

## Identify respondent ID (mice will treat as class/cluster ID)
pred[, "id"] <- -2

## Use passive imputation only for abs_sqrt_agediff
pred["abs_sqrt_agediff", ] <- 0

pred


## SPECIFY IMPUTATION METHODS --------------------------------------------------

meth <- make.method(mc_ong)
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
  "p_artuse_part_bin"
)

l1_count <- c(
  "p_age_i1",
  "durat_wks",
  "cond.prob",
  "ai.rate.52",
  "oi.rate.52"
)

## Level 2 variables to impute
l2 <- c(
  "an_prep_current",
  "prep_revised",
  "stitest_perweek",
  "pnua_12m",
  "cuml.pnum",
  "mmconc"
)

## Methods
meth["abs_sqrt_agediff"] <- "~I(abs(sqrt(ego.age) - sqrt(p_age_i1)))"
meth[c(l1_cat, l1_count)] <- "2l.pmm"
meth[l1_bin] <- "2l.bin"
meth[l2] <- "2lonly.pmm"

meth

sapply(mc_ong, class)

# initialize a mice object
ini <- mice(
  mc_ong,
  maxit = 0,
  predictorMatrix = pred,
  methods = meth
)

## Alter visit sequence so that abs_sqrt_agediff is passively
## imputed right after p_age_i1 is imputed.
vis <- ini$visitSequence
vis <- c(vis[1:15], vis[26], vis[16:25])
vis

## POST-PROCESSING -------------------------------------------------------------

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
## Current PrEP use:
##
##   + Respondents were asked about current PrEP use only if they reported ever
##     having taken PrEP.

post <- ini$post

post["p_prepuse"] <-
  "imp[[j]][data$ego.hiv[!r[, j]] %in% 1:2, i] <- 0"

post["p_artuse_bin"] <-
  "imp[[j]][data$ego.hiv[!r[, j]] %in% c(0, 2), i] <- 0"

post["p_prepuse_part"] <-
  "imp[[j]][data$p_hiv2[!r[, j]] == 1, i] <- 0"

post["p_artuse_part_bin"] <-
  "imp[[j]][data$p_hiv2[!r[, j]] %in% c(0, 2), i] <- 0"

post["an_prep_current"] <-
  "imp[[j]][data$prep_revised[!r[, j]] == 0, i] <- 0"

post["prep_revised"] <-
  "imp[[j]][data$ego.hiv[!r[, j]] == 1, i] <- 0"


## IMPUTE ----------------------------------------------------------------------

## NOTE: Donor pool size selected based on:

##       Morris TP, White IR, Royston P. Tuning multiple imputation by
##       predictive mean matching and local residual draws. BMC Med. Res.
##       Methodol. 2014;14:75. http://dx.doi.org/10.1186/1471-2288-14-75

imp <- parlmice(
  mc_ong,
  maxit = 50,
  n.core = 5,
  n.imp.core = 4,
  predictorMatrix = pred,
  method = meth,
  visitSequence = vis,
  post = post,
  cluster.seed = 45345,
  donors = 10L,
  cl.loads = TRUE
)

## Imputation models kick out some variables with perfect prediction due to
## skip patterns. Dealt with via post-processing above.
imp$loggedEvents


################################################################################
                           ## DIAGNOSE IMPUTATIONS ##
################################################################################


## Check for impossible combinations (violations of skip patterns from ARTnet)
cimp <- as.data.table(complete(imp, "long"))

### PrEP, ego
mc_ong[, .N, keyby = .(ego.hiv, prep_revised)]   # original data
cimp[, .N, keyby = .(ego.hiv, prep_revised)]     # imputed data

mc_ong[, .N, keyby = .(ego.hiv, p_prepuse)]      # original data
cimp[, .N, keyby = .(ego.hiv, p_prepuse)]        # imputed data

mc_ong[, .N, keyby = .(an_prep_current, prep_revised)]  # original data
cimp[, .N, keyby = .(an_prep_current, prep_revised)]    # imputed data


### PrEP, alter
mc_ong[, .N, keyby = .(p_hiv2, p_prepuse_part)]  # original data
cimp[, .N, keyby = .(p_hiv2, p_prepuse_part)]    # imputed data

### ART, ego
mc_ong[, .N, keyby = .(ego.hiv, p_artuse_bin)]      # original data
cimp[, .N, keyby = .(ego.hiv, p_artuse_bin)]        # imputed data

### ART, alter
mc_ong[, .N, keyby = .(p_hiv2, p_artuse_part_bin)]  # original data
cimp[, .N, keyby = .(p_hiv2, p_artuse_part_bin)]    # imputed data


impdx_dir <- "imputation_diagnostics"

### Trace plots
pdf(file.path(impdx_dir, "imp_trace.pdf"), onefile = TRUE)
plot(imp)
dev.off()

### Density plots
pdf(file.path(impdx_dir, "imp_density.pdf"), onefile = TRUE)
densityplot(
  imp,
  thicker = 3,
  scales = list(
    x = list(relation = "free"),
    y = list(relation = "free")
  )
)
dev.off()

### closer look at oral and anal act rates
pdf(file.path(impdx_dir, "imp_density_actrates-only.pdf"))
densityplot(imp, ~ oi.rate.52 + ai.rate.52)
dev.off()


################################################################################
                                  ## WRITE ##
################################################################################

saveRDS(imp, file.path(Sys.getenv("ARTNET_PATH"), "artnet-maincas-imputed.Rds"))


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
  data = mc_ong
)

pool(with(imp2, lm(p_age_i1 ~ 1))) # quick test to get pooled means with standard errors

cbind(test.ai.pooled$pooled$estimate, coef(test.ai.cc))

lapply(test.ai$analyses, predict, type = "response")

class(test.ai.pooled)
summary(test.ai.pooled)
predict(test.ai.pooled)
