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
  mice
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
## SET UP PARTNERSHIP DATA SET ##
################################################################################

anlsub <- copy(anl)

## NOTE: PrEP and ART use are coded a little differently than most partnership
## variables, where the "p_" prefix indicates the alter's characteristic.
## For PrEP and ART, p_artuse, for instance, indicates the ego's ART use with the given alter, and the "_part" suffix indicates the alter's status.

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
mc_ong[p_prepuse %in% c(1,2), p_prepuse := 1]  # any PrEP use in relationship

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

## Fix NA coding 
mc_ong[ego.anal.role == "", ego.anal.role := NA]
mc_ong[p_race.cat == "", p_race.cat := NA]

weekly_rates <- c("ai.rate", "oi.rate")
mc_ong[, (weekly_rates) := NULL]

mc_ong[, .N, keyby = ego.hiv]
mc_ong[, .N, keyby = p_hiv2]

## Convert character variables to factor
mc_ong[, ego.anal.role := as.factor(ego.anal.role)]
mc_ong[, ego.race.cat := as.factor(ego.race.cat)]
mc_ong[, p_race.cat := as.factor(p_race.cat)]

sapply(mc_ong, class)


################################################################################
## MISSING DATA SUMMARIES ##
################################################################################

## Calculate proportion of missing data for each variable
## N = partnerships
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

## SPECIFY PREDICTOR VARIABLES
pred <- make.predictorMatrix(mc_ong)

completevars <- c(
  "id", "pid", "pid_unique",
  "ego.race.cat", "ego.age", "ego.hiv",
  "p_hiv2", "ptype"
)

## Don't impute complete variables
pred[completevars, ] <- 0
pred

## Don't impute using ID variables
pred[, c("id", "pid", "pid_unique")] <- 0

## Identify respondent ID (mice will treat as class/cluster ID)
pred[, "id"] <- -2

## Use passive imputation only for abs_sqrt_agediff
pred["abs_sqrt_agediff", ] <- 0

pred

## SPECIFY IMPUTATION METHODS
meth <- make.method(mc_ong)

# TODO: Decide on final imputation model forms and additional covariates
#       to use in prediction

imp_cat <- c("ego.anal.role", "p_race.cat")

imp_bin <- c(
  "prep_revised", "p_prepuse", "p_prepuse_part",
  "p_artuse_bin", "p_artuse_part_bin"
)

meth["abs_sqrt_agediff"] <- "~I(abs(sqrt(ego.age) - sqrt(p_age_i1)))"
meth[imp_cat] <- "polyreg"
meth[imp_bin] <- "logreg"

meth

sapply(mc_ong, class)

## IMPUTE

## NOTE: Donor pool size selected based on:

##       Morris TP, White IR, Royston P. Tuning multiple imputation by
##       predictive mean matching and local residual draws. BMC Med. Res.
##       Methodol. 2014;14:75. http://dx.doi.org/10.1186/1471-2288-14-75

imp <- parlmice(
  mc_ong,
  maxit = 20,
  n.core = 5,
  n.imp.core = 4,
  predictorMatrix = pred,
  method = meth,
  cluster.seed = 45345,
  donors = 10L
)

imp$loggedEvents


################################################################################
## WRITE ##
################################################################################

saveRDS(imp, file.path(Sys.getenv("ARTNET_PATH"), "artnet-long-imputed.Rds"))


################################################################################
## TEST ANALYSIS ##
################################################################################

imp2 <- readRDS(file.path(Sys.getenv("ARTNET_PATH"), "artnet-long-imputed.Rds"))

imp_complete <- as.data.table(mice::complete(imp2, "long"))
sapply(imp_complete, function(x) sum(is.na(x)))

# make sure no incompatible PrEP/HIV and ART/HIV combos exist
imp_complete[, .N, keyby = .(ego.hiv, p_prepuse)]
imp_complete[, .N, keyby = .(p_hiv2, p_prepuse_part)]
imp_complete[, .N, keyby = .(ego.hiv, p_artuse_bin)]
imp_complete[, .N, keyby = .(p_hiv2, p_artuse_part_bin)]

plot(imp2, layout = c(4, 8))
densityplot(imp2, layout = c(3, 5))
stripplot(imp2, ai.rate.52 + oi.rate.52 + p_age_i1 + abs_sqrt_agediff ~ .imp)

## Ref. for variance estimation of model predictions
## Miles A. Obtaining Predictions from Models Fit to Multiply Imputed Data. Sociol. Methods Res. [electronic article]. 2016;45(1):175–185. (https://doi.org/10.1177/0049124115610345)

test.ai <- with(imp2, glm.nb(
  ai.rate.52 ~
    ego.race.cat + p_race.cat +
    ego.age + p_age_i1 +
    abs_sqrt_agediff
))

lapply(test.ai$analyses, predict, type = "response")

test.ai.pooled <- pool(test.ai)
summary(test.ai.pooled)
