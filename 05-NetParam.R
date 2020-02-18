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
              "p_startyyyy", "p_startyyyydk", "p_startmm",
              "ego.anal.role")

maincas_ong <- anl[ptype %in% 1:2 & p_ongoing_ind == 1, .N, ong_cols]
setkeyv(maincas_ong, cols = c("id", "pid"))

maincas_ong[, table(ptype, psubtype, exclude = NULL)]

# convert subdate to date
maincas_ong[, sub_date := ymd(sub_date)]
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

# impute start month
maincas_ong[is.na(p_startmm),
            p_startmm := sample(1:12, size = 1),
            by = .(id, pid)]

set.seed(283798)

maincas_ong[,
  p_startdt := ymd(
    paste(
      p_startyyyy, p_startmm,
      # sample day of the month (all rows)
      ((p_startmm %in% c(1, 3, 5, 7, 8, 10, 12)) * sample(1:31, size = 1)) +
        ((p_startmm %in% c(4, 6, 9, 11)) * sample(1:30, size = 1)) +
        ((p_startmm == 2) * sample(1:28, size = 1)),
      sep = "-")),
  by = .(id, pid)]

maincas_ong[is.na(p_startdt)] %>% print

datecols <- c("ptype", "p_startyyyy", "p_startyyyydk", "p_startmm", "p_startdt")

# if provided year range of partnership start date, impute
# set upper limit to 15
startdt_match <- list(ydk1 = 1:364,
                      ydk2 = 365 : (2 * 365 - 1),
                      ydk3 = (2 * 365) : (5 * 365 - 1),
                      ydk4 = (5 * 365) : (10 * 365 - 1),
                      ydk5 = (10 * 365) : (15 * 364 - 1))

sapply(startdt_match, range)

maincas_ong[!is.na(p_startyyyydk),
            p_startdt := sub_date - sample(startdt_match[[p_startyyyydk]],
                                           size = 1),
            by = .(id, pid)]

maincas_ong[is.na(p_startdt)]

# @TODO 2020-01-28: fix this issue more naturally
# set abs to account for imputed dates that were after the sub date
maincas_ong[, durat_days := abs(as.numeric(sub_date - p_startdt))]
maincas_ong[, durat_wks := round(durat_days / 7)]
class(maincas_ong$durat_days)
class(maincas_ong$durat_wks)

maincas_ong[, .(id, pid, p_startyyyy, p_startmm,
                p_startdt, sub_date, durat_days, durat_wks)]

maincas_ong[, table(durat_wks < 0)]

ggplot(maincas_ong,
       aes(x = durat_wks)) +
       geom_density() +
       facet_wrap(~ ptype) +
       theme_clean()

fit_all_dur <- maincas_ong[, glm(durat_wks ~ 1, family = "quasipoisson")]

fit_main_dur <- maincas_ong[ptype == 1,
                            glm(durat_wks ~ 1, family = "quasipoisson")]

fit_casl_dur <- maincas_ong[ptype == 2,
                            glm(durat_wks ~ 1, family = "quasipoisson")]

durats <- rbind(
  add_ci(fit_all_dur, tb = data.frame(1)),
  add_ci(fit_main_dur, tb = data.frame(1)),
  add_ci(fit_casl_dur, tb = data.frame(1))
)

durats$X1 <- NULL
durats$ptype <- c("overall", "main", "casual")

setDT(durats)

print(durats)

# %% ONE-TIME PARTNERSHIPS -----------------------------------------------------

# onetime partnerships
ponetime <- anl[ptype == 3, .N, .(id, psubtype, ptype)]


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

# # if initial casual estimate missing, assign casual_init degree as 0
# deg_data[is.na(casl_init), casl_init := 0]
#
# # use self-reported ongoing partnerships to calculate ongoing causal,
# # assuming any non-main, ongoing partnership is casual
# deg_data[!is.na(pn_ongoing),
#          casl_upper :=
#           (casl_init > (pn_ongoing - main)) * casl_init +
#           (casl_init < (pn_ongoing - main)) * (pn_ongoing - main)]
#
# deg_data[is.na(pn_ongoing), casl_upper := casl_init]
#
# deg_data[, tot_pships_init := main + casl_init]
# deg_data[, tot_pships_upper := main + casl_upper]
#
# print(deg_data)
#
# # One-time partnerships
# onetime_data <- merge(ponetime[, .(id, onetime = N)],
#                       an[, .(id, race.cat, age5)],
#                       by = "id")
#
# dfSummary(onetime_data, graph.col = F)


# %% INSPECT DEGREE/ONETIME DATASETS ---------------------------------------------

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


# %% RACE MIXING -----------------------------------------------------------------

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


# %% AGE MIXING ------------------------------------------------------------------

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


# %% EGO SEX ROLE ----------------------------------------------------------------

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


# %% PARTNERSHIP REUSABLES -------------------------------------------------------

race_unq <- data.table(race.cat = sort(unique(deg_data$race.cat)))
age5_unq <- data.table(age5 = sort(unique(deg_data$age5)))

main_degt2_unq <- data.table(
  degmain_trunc2 = sort(unique(deg_data$degmain_trunc2))
)

casual_deg_unq <- data.table(degcasl = sort(unique(deg_data$degcasl)))

cilabs <- c("ll95", "ul95")


# %% MAIN PARTNERSHIPS -----------------------------------------------------------

# Fits
fit_main_deg_byrace <- glm(degmain ~ race.cat,
                           data = deg_data,
                           family = "quasipoisson")

fit_main_deg_byage5 <- glm(degmain ~ factor(age5),
                           data = deg_data,
                           family = "quasipoisson")

fit_main_deg_bycasl <- glm(degmain ~ factor(degcasl),
                           data = deg_data,
                           family = "quasipoisson")

# Predictions

## nodefactor("race")
main_deg_byrace_preds <- add_ci(race_unq,
                                fit_main_deg_byrace,
                                names = cilabs)
print(main_deg_byrace_preds)

## nodefactor("age5")
main_deg_byage5_preds <- add_ci(age5_unq,
                                fit_main_deg_byage5,
                                names = cilabs)
print(main_deg_byage5_preds)

## nodefactor("degcasl")
main_deg_bycasldeg_preds <- add_ci(casual_deg_unq,
                                   fit_main_deg_bycasl,
                                   names = cilabs)
print(main_deg_bycasldeg_preds)


main_concurrent_prob <- deg_data[, mean(main_conc_ind)]
round(main_concurrent_prob, 3)


# %% CASUAL PARTNERSHIPS -------------------------------------------------------

# Fits
fit_casl_deg_byrace <- glm(degcasl ~ race.cat,
                           data = deg_data,
                           family = "quasipoisson")

fit_casl_deg_byage5 <- glm(degcasl ~ factor(age5),
                           data = deg_data,
                           family = "quasipoisson")

fit_casl_deg_bymaint2 <- glm(degcasl ~ factor(degmain_trunc2),
                             data = deg_data,
                             family = "quasipoisson")

# Predictions

## nodefactor("race")
casl_deg_byrace_preds <- add_ci(race_unq,
                                fit_casl_deg_byrace,
                                names = cilabs)
print(casl_deg_byrace_preds)

## nodefactor("age5")
casl_deg_byage5_preds <- add_ci(age5_unq,
                                fit_casl_deg_byage5,
                                names = cilabs)
print(casl_deg_byage5_preds)

## nodefactor("degmain")
casl_deg_bymaindegt2_preds <- add_ci(main_degt2_unq,
                                   fit_casl_deg_bymaint2,
                                   names = cilabs)
print(casl_deg_bymaindegt2_preds)

casl_concurrent_prob <- deg_data[, mean(casl_conc_ind)]
round(casl_concurrent_prob, 3)


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

saveRDS(list(main_summaries = main_summaries,
             casl_summaries = casl_summaries),
        file = "netstats/aggregate_degree_summaries.Rds")

# Main predictions

main_predictions <- list()

# degree
main_predictions$deg_byrace <- main_deg_byrace_preds
main_predictions$deg_byage5 <- main_deg_byage5_preds
main_predictions$deg_bycasltot <- main_deg_bycasldeg_preds
main_predictions$concurrent_prob <- main_concurrent_prob
main_predictions$mean_durat_wks <- durats[ptype == "main", pred]
main_predictions$absdiff_sqrtage <- anl[ptype == 1,
                                        mean(abs_sqrt_agediff, na.rm = T)]

print(main_predictions)

# Casual predictions

casl_predictions <- list()

# degree
casl_predictions$deg_byrace <- casl_deg_byrace_preds
casl_predictions$deg_byage5 <- casl_deg_byage5_preds
casl_predictions$deg_bymaindegt2 <- casl_deg_bymaindegt2_preds
casl_predictions$mean_durat_wks <- durats[ptype == "casual", pred]
casl_predictions$concurrent_prob <- casl_concurrent_prob
casl_predictions$absdiff_sqrtage <- anl[ptype == 2,
                                        mean(abs_sqrt_agediff, na.rm = T)]

print(casl_predictions)

netstats <- list(
  demo = list(ai.role.pr = mc_analrole),
  main = main_predictions,
  casl = casl_predictions)



# Other network stats

# race matches
netstats$main$racematch <- mp_racematch
netstats$casl$racematch <- cp_racematch

# age mixing
netstats$main$agemix <- mp_agemix
netstats$casl$agemix <- cp_agemix

# age5 matches
netstats$main$age5match <- mp_age5match
netstats$casl$age5match <- mp_age5match

# anal sex role
netstats$main$ai.role <- mp_analrole
netstats$casl$ai.role <- cp_analrole

print(netstats)

saveRDS(netstats, file = "netstats/predictions.Rds")
