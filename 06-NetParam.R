# %% SETUP ---------------------------------------------------------------------

source("01-Import-Private-Data.R")

pacman::p_load(data.table,
               summarytools,
               ciTools,
               dplyr,
               stringr,
               lubridate)

# function to calculate number of unique levels in a variable
unql <- function(data) length(unique(data))

# import cleaned wide dataset
an <- fread(paste0(Sys.getenv("ARTNET_PATH"), "/artnet-wide-cleaned.csv"))
anl <- fread(paste0(Sys.getenv("ARTNET_PATH"), "/artnet-long-cleaned.csv"))


# %% MAIN AND CAUSAL PARTNERSHIPS ----------------------------------------------

egos <- length(unique(anl$id))

total_main <- unname(table(anl$ptype == 1)[2])
total_casl <- unname(table(anl$ptype == 2)[2])
total_pships <- anl[!is.na(ptype), length(unique(paste0(id, pid)))]
sum(total_main, total_casl)

cat(
  sep = "",
  "AMONG MSM WITH PARTNERSHIPS\n", rep("=", 30), "\n\n",
  "Unique egos\n", format(egos, big.mark = ","), "\n\n",
  "Unique main partnerships\n", format(total_main, big.mark = ","), "\n\n",
  "Unique casual partnerships\n", format(total_casl, big.mark = ","), "\n\n",
  "Overall partnerships\n", format(total_pships, big.mark = ",")
)

# ongoing main or casual partnerships
ong_cols <- c("sub_date", "id", "pid",
              "ptype", "psubtype",
              "p_startyyyy", "p_startyyyydk", "p_startmm", "p_startdt",
              "durat_days", "durat_wks",
              "ego.race.cat", "ego.age5",
              "ego.anal.role")

maincas_ong <- anl[ptype %in% 1:2 & p_ongoing_ind == 1, .N, ong_cols]
setkeyv(maincas_ong, cols = c("id", "pid"))

maincas_ong[, table(ptype, psubtype, exclude = NULL)]

# convert startdt to date
maincas_ong[, p_startdt := ymd(p_startdt)]
str(maincas_ong)

# ongoing main and casual partnerships
sum(maincas_ong$ptype == 1)
sum(maincas_ong$ptype == 2)


# drops 1 partnership with no sex
# drops 12 partnerships with missing data on psubtype
ong_ptype_cts <- dcast(data = maincas_ong[!psubtype %in% c("no sex", NA)],
                       formula = id ~ psubtype + ptype,
                       fun.aggregate = sum,
                       value.var = "N")

print(ong_ptype_cts)


# %% PARTNERSHIP DURATIONS -----------------------------------------------------

maincas_ong[, .(
  id, pid, ego.race.cat, ego.age5,
  ptype, p_startyyyy, p_startmm, p_startdt, sub_date,
  durat_days, durat_wks
)]

maincas_ong[, table(durat_wks < 0)]

ggplot(maincas_ong,
       aes(x = durat_wks)) +
       geom_histogram(col = "black", fill = "aliceblue") +
       facet_wrap(~ ptype, nrow = 2) +
       theme_clean()

fit_main_dur <- maincas_ong[ptype == 1,
                            glm(durat_wks ~ 1, family = "quasipoisson")]

fit_casl_dur <- maincas_ong[ptype == 2,
                            glm(durat_wks ~ 1, family = "quasipoisson")]

durats <- rbind(
  add_ci(fit_main_dur, tb = data.frame(1)),
  add_ci(fit_casl_dur, tb = data.frame(1))
)

durats$X1 <- NULL
durats$ptype <- c("main", "casual")

setDT(durats)

print(durats)


# %% IMPUTE AGE --------------------------------------------------------------

# Source for age imputation:

# Weiss KM, Goodreau SM, Morris M, Prasad P, Ramaraju R, Sanchez T, et al.
# Egocentric Sexual Networks of Men Who Have Sex with Men in the United States:
# Results from the ARTnet Study. medRxiv. 2019. doi:10.1101/19010579

# @NOTE: My method varied slightly from theirs

plot(density(anl$p_age, na.rm = T))
anl[, .N, key = p_age]

# set age == 0 to missing
anl[, p_age := ifelse(p_age == 0, NA, p_age)]
print(anl[, .(id, p_age)])

sum(is.na(anl$p_age))

anl[is.na(p_age) & is.na(p_relage), sum(.N)]
anl[, .N, key = p_relage] %>% print
class(anl$p_age)

# Based on PARTXRELAGE (ART-Net survey)
set.seed(1998)

anl[, p_age_imputed := dplyr::case_when(
          p_relage == 1 ~ ego.age - sample(size = 1, x = 10:15),
          p_relage == 2 ~ ego.age - sample(size = 1, x = 2:10),
          p_relage == 3 ~ ego.age + sample(size = 1, x = -1:1),
          p_relage == 4 ~ ego.age + sample(size = 1, x = 2:10),
          p_relage == 5 ~ ego.age + sample(size = 1, x = 10:15),
          !is.na(p_age) ~ p_age,
          TRUE ~ NA_integer_),
      .(id, pid)] %>%
  .[, p_age5 := dplyr::case_when(
          p_age_imputed <= 24 ~ 1,
          p_age_imputed >= 25 & p_age_imputed <= 34 ~ 2,
          p_age_imputed >= 35 & p_age_imputed <= 44 ~ 3,
          p_age_imputed >= 45 & p_age_imputed <= 54 ~ 4,
          p_age_imputed >= 55 ~ 5,
          TRUE ~ NA_real_
    )]

print(anl[, .(p_age, p_age_imputed)])

anl[, .(agediff = p_age_imputed - ego.age),
      key = p_relage] %>%
  .[, .(min = min(agediff),
        max = max(agediff)),
      key = p_relage]


# %% CALCULATE AGE DIFFERENCES -------------------------------------------------

# calculate age differences
anl[, ":="(abs_agediff = abs(p_age_imputed - ego.age),
           abs_sqrt_agediff = abs(sqrt(p_age_imputed) - sqrt(ego.age))
           )]

anl[, .(abs_agediff, abs_sqrt_agediff)]

# @TODO: the paste0(id, pid) doesn't drop rows that should be dropped
# missing imputed partner age, among egos with partnerships
# 1.8% of eligible observations
anl[!is.na(nn), sum(is.na(p_age_imputed)) / unql(paste0(id, pid))]

print(names(anl))


# %% RACE/ETHNICITY ------------------------------------------------------------

# Create partner version of race.cat
racematch <- c("1" = "other",
               "2" = "black",
               "3" = "white",
               "4" = "other",
               "5" = "other",
               "6" = "other")

# @NOTE:
#  - Default to reported race/ethnicity if hispanic indicator missing
anl[, p_race.cat := ifelse(
        p_hisp == 0 | is.na(p_hisp),
        racematch[p_race],
        "hispanic"
        )]

anl[, .(id,
        pid,
        nn,
        p_hisp,
        p_race,
        p_race.cat)]

anl[, .N, .(id)][, summary(N)]

dfSummary(anl, plain.ascii = T, graph.col = F)


# %% IMPUTE RACE/ETHNICITY ---------------------------------------------------

# @NOTE 2020-01-28
# holding off for now on this; might not be necessary


# %% CALCULATE MAIN AND CASUAL DEGREE ------------------------------------------

# Main and causal partnerships (ongoing)

# join ongoing partner counts with ego traits
deg_data <- ong_ptype_cts[an[, .(id, pn_ongoing,
                                 race.cat, age5)], on = "id"]

setkey(deg_data, id)

# should be 0
sum(is.na(deg_data$race.cat))
sum(is.na(deg_data$age5))

# if partnership type count missing, set to 0
deg_data[, ":="(
      analonly_1 = ifelse(is.na(analonly_1), 0, analonly_1), # main, anal-only
      oralanal_1 = ifelse(is.na(oralanal_1), 0, oralanal_1), # main, oral-anal
      oralonly_1 = ifelse(is.na(oralonly_1), 0, oralonly_1), # main, oral-only
      analonly_2 = ifelse(is.na(analonly_2), 0, analonly_2), # casual, anal-only
      oralanal_2 = ifelse(is.na(oralanal_2), 0, oralanal_2), # casual, oral-anal
      oralonly_2 = ifelse(is.na(oralonly_2), 0, oralonly_2)  # casual, oral-only
    )] %>%
    # main degree
    .[, degmain := analonly_1 + oralanal_1 + oralonly_1] %>%
    .[, degmain_trunc2 := ifelse(degmain >= 2, 2, degmain)] %>%
    # casual degree
    .[, degcasl := analonly_2 + oralanal_2 + oralonly_2] %>%
    # total degree (main + casual)
    .[, degtotal := degmain + degcasl] %>%
    # concurrency indicators
    .[, main_conc_ind := ifelse(degmain >= 2, 1, 0)] %>%
    .[, casl_conc_ind := ifelse(degcasl >= 2, 1, 0)] %>%
    .[, maincasl_conc_ind := degtotal > 1]

str(deg_data)

deg_data[, .N, degmain][, P := N / sum(N)] %>% print
degmain_t2 <- deg_data[, .N, degmain_trunc2][, P := N / sum(N)]
print(degmain_t2)

deg_data[, .N, main_conc_ind][, P := N / sum(N)] %>% print
deg_data[, .N, keyby = degcasl][, P := N / sum(N)] %>% print
deg_data[, .N, casl_conc_ind][, P := N / sum(N)] %>% print
deg_data[, .N, maincasl_conc_ind][, P := N / sum(N)] %>% print

names(deg_data)

# create dummy variables for main and casl degrees
deg_data %>%
  .[, ":="(
      casl0 = dplyr::case_when(degcasl == 0 ~ 1,
                               degcasl %in% 1:5 ~ 0,
                               TRUE ~ NA_real_),
      casl1 = dplyr::case_when(degcasl == 1 ~ 1,
                               degcasl %in% c(0, 2:5) ~ 0,
                               TRUE ~ NA_real_),
      casl2 = dplyr::case_when(degcasl == 2 ~ 1,
                               degcasl %in% c(0:1, 3:5) ~ 0,
                               TRUE ~ NA_real_),
      casl3 = dplyr::case_when(degcasl == 3 ~ 1,
                               degcasl %in% c(0:2, 4:5) ~ 0,
                               TRUE ~ NA_real_),
      casl4 = dplyr::case_when(degcasl == 4 ~ 1,
                               degcasl %in% c(0:3, 5) ~ 0,
                               TRUE ~ NA_real_),
      casl5 = dplyr::case_when(degcasl == 5 ~ 1,
                               degcasl %in% 0:4 ~ 0,
                               TRUE ~ NA_real_),
      main0 = dplyr::case_when(degmain_trunc2 == 0 ~ 1,
                               degmain_trunc2 %in% 1:2 ~ 0,
                               TRUE ~ NA_real_),
      main1 = dplyr::case_when(degmain_trunc2 == 1 ~ 1,
                               degmain_trunc2 %in% c(0, 2) ~ 0,
                               TRUE ~ NA_real_),
      main2 = dplyr::case_when(degmain_trunc2 == 2 ~ 1,
                               degmain_trunc2 %in% 0:1 ~ 0,
                               TRUE ~ NA_real_)
    )]

print(deg_data)


# %% INSPECT DEGREE/ONETIME DATASETS -------------------------------------------

# oa = oral and anal sex partnership
# ao = anal-only partnership
# oo = oral-only partnership

deg_data[, .(mn_degmain = mean(degmain),
             mn_oa_1 = mean(oralanal_1),
             mn_ao_1 = mean(analonly_1),
             mn_oo_1 = mean(oralonly_1))]

sum_main_total <- deg_data[, .N, key = degmain][, P := round(N / sum(N), 3)]
sum_main_oa    <- deg_data[, .N, key = oralanal_1][, P := round(N / sum(N), 3)]
sum_main_ao    <- deg_data[, .N, key = analonly_1][, P := round(N / sum(N), 3)]
sum_main_oo    <- deg_data[, .N, key = oralonly_1][, P := round(N / sum(N), 3)]

sum_casl_total <- deg_data[, .N, key = degcasl][, P := round(N / sum(N), 3)]
sum_casl_oa    <- deg_data[, .N, key = oralanal_2][, P := round(N / sum(N), 3)]
sum_casl_ao    <- deg_data[, .N, key = analonly_2][, P := round(N / sum(N), 3)]
sum_casl_oo    <- deg_data[, .N, key = oralonly_2][, P := round(N / sum(N), 3)]

pcols <- c("id",
           "pid",
           "ptype",
           "psubtype",
           "p_race.cat",
           "p_age5",
           "p_hiv",
           names(anl)[grep("^ego", names(anl))])


# %% RACE MIXING ---------------------------------------------------------------

mp <- anl[ptype == 1 & p_ongoing_ind == 1, ..pcols]
cp <- anl[ptype == 2 & p_ongoing_ind == 1, ..pcols]

# Mixing matrix, main partnership

mp_racemix <- mp[, ctable(ego.race.cat, p_race.cat, prop = "n", useNA = "no")]
mp_mixmat <- with(mp, round(prop.table(table(ego.race.cat, p_race.cat)), 3))
print(mp_racemix)

# Mixing matrix, casual partnership

cp_racemix <- cp[, ctable(ego.race.cat, p_race.cat, prop = "n", useNA = "no")]
cp_mixmat <- with(cp, round(prop.table(table(ego.race.cat, p_race.cat)), 3))
print(cp_racemix)

# Partner race probabilities, main partnerships

mp_racematch <- mp[!is.na(p_race.cat)] %>%
 .[, samerace := ifelse(ego.race.cat == p_race.cat, 1, 0)] %>%
 .[, .N, by = "samerace"] %>%
 .[, P := round(N / sum(N), 4)]

print(mp_racematch)

# Partner race probabilities, casual partnerships

cp_racematch <- cp[!is.na(p_race.cat)] %>%
 .[, samerace := ifelse(ego.race.cat == p_race.cat, 1, 0)] %>%
 .[, .N, by = "samerace"] %>%
 .[, P := round(N / sum(N), 4)]

print(cp_racematch)

# plot mixing proportions by ego raceeth
# plotmixing <- function(mixdat, title) {
#   ggplot(data = mixdat,
#          mapping = aes(x = ego.race.cat,
#                        y = P,
#                        group = p_race.cat,
#                        fill = p_race.cat)) +
#     geom_col(position = "dodge",
#              col = "white",
#              width = 0.3) +
#     scale_fill_viridis_d(name = "Partner\nRace/ethnicity") +
#     labs(x = "Ego Race/ethnicity",
#          y = "Probability",
#          title = title) +
#     theme_classic()
# }
#
# plotmixing(mp_racesummary, "Mixing in main partnerships")
# plotmixing(cp_racesummary, "Mixing in casual partnerships")


# %% AGE MIXING ----------------------------------------------------------------

# Mixing

mp_agemix <- mp[!is.na(p_age5)] %>%
  .[, .N, keyby = .(ego.age5, p_age5)] %>%
  .[, ":=" (minage = min(ego.age5, p_age5),
            maxage = max(ego.age5, p_age5)), by = 1:nrow(.)] %>%
  .[, agecomb := paste0(minage, maxage)] %>%
  .[, .(N = sum(N)), agecomb] %>%
  .[, P := round(N / sum(N), 4)]

print(mp_agemix)

cp_agemix <- cp[!is.na(p_age5)] %>%
  .[, .N, keyby = .(ego.age5, p_age5)] %>%
  .[, ":=" (minage = min(ego.age5, p_age5),
            maxage = max(ego.age5, p_age5)), by = 1:nrow(.)] %>%
  .[, agecomb := paste0(minage, maxage)] %>%
  .[, .(N = sum(N)), agecomb] %>%
  .[, P := round(N / sum(N), 4)]

print(cp_agemix)

# Matching

mp_age5match <- mp[!is.na(p_age5)] %>%
  .[, .N, keyby = .(ego.age5, p_age5)] %>%
  .[, sameage := ifelse(ego.age5 == p_age5, 1, 0)] %>%
  .[, .(N = sum(N)), sameage] %>%
  .[, P := round(N / sum(N), 4)]

print(mp_age5match)

cp_age5match <- cp[!is.na(p_age5)] %>%
  .[, .N, keyby = .(ego.age5, p_age5)] %>%
  .[, sameage := ifelse(ego.age5 == p_age5, 1, 0)] %>%
  .[, .(N = sum(N)), sameage] %>%
  .[, P := round(N / sum(N), 4)]

print(cp_age5match)


# %% EGO SEX ROLE --------------------------------------------------------------

# @TODO:
# - Impute missing anal.sex.role

# Overall (among ids with main or casual partnerships)

mc_analrole <- dcast(anl[ptype %in% 1:2],
                     id ~ ego.anal.role,
                     value.var = "ego.anal.role") %>%
  .[, anal.sex.role := case_when(
    Versatile > 0 | (Insertive == 1 & Receptive == 1) ~ "V",
    Receptive > 0 & (Insertive == 0 & Versatile == 0) ~ "R",
    Insertive > 0 & (Receptive == 0 & Versatile == 0) ~ "I",
    TRUE ~ NA_character_
  )]

# about 11% missing, ego-wise
mc_analrole[, .N, anal.sex.role][, P := N / sum(N)] %>% print

mc_analrole <- mc_analrole[!is.na(anal.sex.role),
                           .N, anal.sex.role][, P := N / sum(N)]

print(mc_analrole)

# Main Partnerships

# between 9-10% missing, partnership-wise
table(is.na(mp$ego.anal.role))


mp_analrole <- dcast(mp, id ~ ego.anal.role, value.var = "ego.anal.role") %>%
  .[, anal.sex.role := case_when(
    Versatile > 0 | (Insertive == 1 & Receptive == 1) ~ "V",
    Receptive > 0 & (Insertive == 0 & Versatile == 0) ~ "R",
    Insertive > 0 & (Receptive == 0 & Versatile == 0) ~ "I",
    TRUE ~ NA_character_
  )]

# about 9% missing, ego-wise
mp_analrole[, .N, anal.sex.role][, P := N / sum(N)] %>% print

mp_analrole <- mp_analrole[!is.na(anal.sex.role),
                           .N, anal.sex.role][, P := N / sum(N)]

print(mp_analrole)


# Casual Partnerships

# about 21% missing, partnership-wise
table(is.na(cp$ego.anal.role))

cp_analrole <- dcast(cp, id ~ ego.anal.role, value.var = "ego.anal.role") %>%
  .[, anal.sex.role := case_when(
    Versatile > 0 | (Insertive == 1 & Receptive == 1) ~ "V",
    Receptive > 0 & (Insertive == 0 & Versatile == 0) ~ "R",
    Insertive > 0 & (Receptive == 0 & Versatile == 0) ~ "I",
    TRUE ~ NA_character_
  )]

# about 16% missing, ego-wise
cp_analrole[, .N, anal.sex.role][, P := N / sum(N)] %>% print

cp_analrole <- cp_analrole[!is.na(anal.sex.role),
                           .N, anal.sex.role][, P := N / sum(N)]

print(cp_analrole)


# %% PARTNERSHIP REUSABLES -----------------------------------------------------

# store unique levels for categorical variables
race_unq <- data.table(race.cat = sort(unique(deg_data$race.cat)))
age5_unq <- data.table(age5 = sort(unique(deg_data$age5)))

main_degt2_unq <- data.table(
  degmain_trunc2 = sort(unique(deg_data$degmain_trunc2))
)

casual_deg_unq <- data.table(degcasl = sort(unique(deg_data$degcasl)))

# labels for ciTools::add_ci output
cilabs <- c("ll95", "ul95")

# set up tables for predictor level combinations
ra_grid <- expand.grid(
  race.cat = unlist(race_unq),
  age5 = unlist(age5_unq)
)

rac_grid <- expand.grid(
  race.cat = unlist(race_unq),
  age5 = unlist(age5_unq),
  degcasl = unlist(casual_deg_unq)
)

ram_grid <- expand.grid(
  race.cat = unlist(race_unq),
  age5 = unlist(age5_unq),
  degmain_trunc2 = unlist(main_degt2_unq)
)

racm_grid <- expand.grid(
  race.cat = unlist(race_unq),
  age5 = unlist(age5_unq),
  degcasl = unlist(casual_deg_unq),
  degmain_trunc2 = unlist(main_degt2_unq)
)


# %% MAIN PARTNERSHIPS ---------------------------------------------------------

## ... PROBS FOR SEEDING MAIN DEGREE

# main degree dummies
main_bin <- names(deg_data)[grepl("main[0-2]{1}", names(deg_data))]

# function to generate predictions
predict_degree_prob <- function(yvar, xvar, dat, newdat) {
  y <- yvar
  x <- paste(xvar, collapse = "+")

  fit <- glm(paste(y, "~", x), data = dat, family = "binomial")
  preds <- predict(fit, newdat, type = "response")

  as.data.table(cbind(newdat, preds, outcome = y))

}

# used to estimate the distribution of main degree
# marginalized over race and age (used for seeding population)
mainprob_by_ra <- lapply(main_bin, function(x) {
  predict_degree_prob(
    yvar = x,
    xvar = c("race.cat", "factor(age5)"),
    dat = deg_data,
    newdat = ra_grid
  )
})


## ... EXPECTED MAIN DEGREE (BY RACE, AGE, AND CASUAL DEGREE)

fit_mdeg_joint <- glm(
  degmain_trunc2 ~ race.cat + factor(age5) + factor(degcasl),
  data = deg_data,
  family = "quasipoisson"
)

summary(fit_mdeg_joint)

pred_mdeg_joint <- predict(
  fit_mdeg_joint,
  newdata = rac_grid,
  type = "response"
)

pred_mdeg_joint <- cbind(rac_grid, pred_mdeg_joint)
print(pred_mdeg_joint)

# ... MAIN CONCURRENCY

# proportion of individuals with concurrent main partnerships
# main_concurrent_prob <- deg_data[, mean(main_conc_ind)]
# round(main_concurrent_prob, 3)

head(deg_data[, .(id, race.cat, age5, degmain, degcasl, main_conc_ind)])

fit_main_concurrent <- glm(
  main_conc_ind ~ race.cat + factor(age5) + factor(degcasl),
  family = "binomial",
  data = deg_data
)

pred_main_concurrent <- predict(
  fit_main_concurrent,
  newdata = rac_grid,
  type = "response"
)

pred_main_concurrent <- cbind(rac_grid, pred_main_concurrent)

# ... MAIN RELATIONSHIP DURATION

# average duration (in weeks)
main_durat_wks <- durats[ptype == "main", pred]
print(main_durat_wks)


# %% CASUAL PARTNERSHIPS -------------------------------------------------------

## ... PROBS FOR SEEDING CASUAL DEGREE

# casual degree dummies
casl_bin <- names(deg_data)[grepl("casl[0-5]{1}", names(deg_data))]

caslprob_by_ra <- lapply(casl_bin, function(x) {
    predict_degree_prob(yvar = x,
                        xvar = c("race.cat",
                                 "factor(age5)"),
                        dat = deg_data,
                        newdat = ra_grid)
  })

caslprob_by_ra %>% print

# used to estimate the distribution of casual degree
# marginalized over race and age (used for seeding population)
caslprob_by_ra <- lapply(casl_bin, function(x) {
  predict_degree_prob(
    yvar = x,
    xvar = c("race.cat", "factor(age5)"),
    dat = deg_data,
    newdat = ra_grid
  )
})

print(caslprob_by_ra)


## ... EXPECTED MAIN DEGREE (BY RACE, AGE, AND MAIN DEGREE)

fit_cdeg_joint <- glm(
  degcasl ~ race.cat + factor(age5) + factor(degmain_trunc2),
  data = deg_data,
  family = "quasipoisson"
)

summary(fit_cdeg_joint)

pred_cdeg_joint <- predict(
  fit_cdeg_joint,
  newdata = ram_grid,
  type = "response"
)

pred_cdeg_joint <- cbind(ram_grid, pred_cdeg_joint)
print(pred_cdeg_joint)


# ... CASUAL CONCURRENCY

# proportion of individuals with concurrent casual partnerships
# casl_concurrent_prob <- deg_data[, mean(casl_conc_ind)]
# round(casl_concurrent_prob, 3)

fit_casl_concurrent <- glm(
  casl_conc_ind ~ race.cat + factor(age5) + factor(degmain_trunc2),
  family = "binomial",
  data = deg_data
)

pred_casl_concurrent <- predict(
  fit_casl_concurrent,
  newdata = ram_grid,
  type = "response"
)

pred_casl_concurrent <- cbind(ram_grid, pred_casl_concurrent)
print(pred_casl_concurrent)

# ... CASUAL RELATIONSHIP DURATION

# average duration (in weeks)
casl_durat_wks <- durats[ptype == "casual", pred]
print(casl_durat_wks)


# %% ONE-TIME PARTNERSHIPS -----------------------------------------------------

ponetime <- anl[
  ptype == 3,
  .(id, pid,
    ego.age, ego.race.cat,
    p_age_imputed, p_race.cat,
    p_rai_once, p_iai_once, p_roi_once, p_ioi_once,
    ptype, psubtype)
    ]

print(ponetime)
ponetime[, .N, psubtype]

ponetime <- anl[ptype == 3, .N, .(id, psubtype, ptype)]

pinst_rate <- an[order(id), .(id, race.cat, age, pnoa_12m, pna_12m, pno_12m)]

# 19 missing
mice::md.pattern(pinst_rate)

pinst_rate[, inst_wkrate := pnoa_12m / 52]

pinst_rate <- pinst_rate[!is.na(inst_wkrate)]

pinst_rate[, plot(density(pnoa_12m))]
pinst_rate[, .(mean = mean(pnoa_12m), var = var(pnoa_12m))]

pinst_rate[, plot(density(inst_wkrate))]
pinst_rate[, .(mean = mean(inst_wkrate), var = var(inst_wkrate))]

# join deg_data and pinst_rate
pinst_rate <- deg_data[
  pinst_rate,
  .(id, age, age5, race.cat,
    degmain_trunc2, degcasl,
    inst_wkrate)
    ]

print(pinst_rate)


## ... EXPECTED WEEKLY PARTNERSHIP RATE (BY RACE, AGE, MAIN, AND CASUAL DEGREE)

fit_inst_joint <- glm(
  inst_wkrate ~
    race.cat + factor(age5) + factor(degcasl) + factor(degmain_trunc2),
  data = pinst_rate,
  family = "quasipoisson"
)

summary(fit_inst_joint)

pred_inst_joint <- predict(
  fit_inst_joint,
  newdata = racm_grid,
  type = "response"
)

pred_inst_joint <- cbind(racm_grid, pred_inst_joint)
print(pred_inst_joint)


# %% WRITE SUMMARIES -----------------------------------------------------------

main_summaries <- list(sum_main_total = sum_main_total,
                       sum_main_trunc = degmain_t2,
                       sum_main_oa    = sum_main_oa,
                       sum_main_ao    = sum_main_ao,
                       sum_main_oo    = sum_main_oo)

casl_summaries <- list(sum_casl_total = sum_casl_total,
                       sum_casl_oa    = sum_casl_oa,
                       sum_casl_ao    = sum_casl_ao,
                       sum_casl_oo    = sum_casl_oo)

saveRDS(list(main_artnet_sum = main_summaries,
             casl_artnet_sum = casl_summaries),
        file = "netstats/aggregate_degree_summaries.Rds")

# Main predictions
main_pred_probs <- mainprob_by_ra %>% rbindlist
main_pred_joint <- as.data.table(pred_mdeg_joint)

# Casual predictions
casl_pred_probs <- caslprob_by_ra %>% rbindlist
casl_pred_joint <- as.data.table(pred_cdeg_joint)

# Instantaneous predictions
inst_pred_joint <- as.data.table(pred_inst_joint)


nparams <- list(

  demo = list(ai.role.pr = mc_analrole),

  main = list(degprob = main_pred_probs,
              degpred_joint = main_pred_joint,
              concurrent = as.data.table(pred_main_concurrent),
              racematch = mp_racematch,
              age5match = mp_age5match,
              durat_wks = main_durat_wks,
              role.class = mp_analrole
            ),

  casl = list(degprob = casl_pred_probs,
              degpred_joint = casl_pred_joint,
              concurrent = as.data.table(pred_casl_concurrent),
              racematch = cp_racematch,
              age5match = cp_age5match,
              durat_wks = casl_durat_wks,
              role.class = cp_analrole
            ),

  inst = list(inst_joint = inst_pred_joint)
)

print(nparams)

saveRDS(nparams, file = "netstats/predictions.Rds")
