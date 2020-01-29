# %% Setup ---------------------------------------------------------------------

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


# %% INDIVIDUAL PARTNER DATA ---------------------------------------------------

# Set up ID table - one row per partnership

pid_tab <- melt(
  an[, c("id", paste0("partnn", 1:5))],
  measure = patterns("partnn")
  ) %>%
  .[, pid := stringr::str_remove(variable, "nn")] %>%
  .[order(id, pid), .(id, pid, nn = value)]

setkey(pid_tab, id, pid)

print(pid_tab)


# %% CONVERT WIDE TO LONG ------------------------------------------------------

# Select partner variables

pc_select <- names(an)[grep("^part[1-5]{1}", names(an))] %>% .[. != "part2013"]
partner_cols <- c("id", pc_select)

print(partner_cols)

# Partner variable suffixes
slugs <- unique(str_extract(partner_cols[-1], "(?<=part[1-5]{1}).*"))
sum(is.na(slugs))
print(sort(slugs))

# Specify sets of columns [(ego id + 5 partners) * 1 variable]
colsets <- lapply(slugs, function(x) {

  pcols <- partner_cols[
    grep(paste0("^id|(?<=part[1-5]{1})", x, "(?!.)"),
         partner_cols,
         perl = T)]

  pcols
  })

print(colsets)

# Create long-format tables for each variable type

# @NOTE:
#  - data.table will kick a warning that some of the columns beings collapsed
#    are not of the same type (some double, some integer).
#  - Output column is converted to double.
#  - Interger-to-double conversion should be safe.

sapply(slugs[grepl("recai|recuai|insai|insuai|once|acts", slugs)], print)

integer_coerce <- slugs[grepl("recai|recuai|insai|insuai|recoi|insoi|once|acts",
                              slugs)]

melted_vars <- lapply(setNames(slugs, slugs), function(x) {

  selpartcols <- grep(paste0("(?<=part[1-5]{1})", x, "(?!.)"),
                      partner_cols,
                      perl = T)

  currcols <- c("id", partner_cols[selpartcols])

  dat <- an[, ..currcols][order(id)]

  # coerce certain variables to integer format (avoid melt conflicts)
  if (x %in% integer_coerce) {
    dat[, (currcols[-1]) := lapply(.SD, as.integer), .SDcols = currcols[-1]]
  }

  dmelt <- melt(dat, measure = patterns("part"))[order(id)]
  dmelt[, pid := str_extract(variable, "part[1-5]{1}")]

  names(dmelt)[names(dmelt) == "value"] <- paste0("p_", x)

  out <- dmelt[, c("id", "pid",
                   paste0("p_", x)), with = F][order(pid)]
  setkey(out, id, pid)

  return(out)
  })

lapply(melted_vars[integer_coerce], function(x) sapply(x[, 3], class))

length(melted_vars)

checknrow <- sapply(melted_vars, nrow) %>% unlist

table(checknrow == nrow(pid_tab))

# Bind partner data into long-format data.table
plong <- Reduce(
  function(x, y) merge(x, y, by = c("id", "pid"), all = T),
  melted_vars
)

# remove double underscores
names(plong) <- str_replace_all(names(plong), "__", "_")
plong <- pid_tab[plong, on = .(id, pid)]

str(plong)
unql(plong$id)


# %% SET UP LONG DATASET ------------------------------------------------------

# pull in ego variables
plong2 <- plong[an[, .(id,
                       sub_date,
                       ego.age = age,
                       ego.age5 = age5,
                       ego.race.cat = race.cat,
                       ego.hiv = hiv.ego,
                       ego.ongoing = pn_ongoing)],
             on = "id"] %>%
  .[order(id, pid)] %>%
  .[, nn := ifelse(nn == "", NA, nn)]

names(plong2)
str(plong2)

# variables coded with 88/99 as "don't know" or "prefer not to answer"
na_8899 <- paste0("p_", c("hisp",
                          "race",
                          "hiv",
                          "hivknow",
                          "startyyyydk",
                          "startmm",
                          "maincasendmm",
                          "ofendmm",
                          "ongoing",
                          "main_ong",
                          "main_end",
                          "iev",
                          "concurr",
                          "concurr2",
                          na.omit(
                            stringr::str_extract(names(anl), "(?<=^p_).*_once$")
                            ),
                          "stidiag",
                          "stitx",
                          "prepuse",
                          "prepstart",
                          "prepuse_part",
                          "prepimp_start",
                          "prepimp_start_part",
                          "prepimp_condoms",
                          "prepimp_condoms_part",
                          "relage",
                          "artuse",
                          "artuse_part",
                          "artstart_part",
                          "artimp_start",
                          "artimp_start_part",
                          "artimp_condoms",
                          "artimp_condoms_part",
                          "artimp_prep"))

# set 88/99 to NA for all variables in na_8899
plong2[, (na_8899) := lapply(.SD, function(x) ifelse(x == 88 | x == 99, NA, x)),
         .SDcols = na_8899]

print(unql(plong2$id))

plong2[is.na(nn) & (!is.na(p_age) | !is.na(p_race) | !is.na(p_hisp)),
    nn := "Nickname"]

# still contains filler rows from the dcast procedure
str(plong2)


# %% CLASSIFY PARTNERSHIP TYPES ------------------------------------------------

# label main partnerships
plong2[p_main_ong == 1 | p_main_end == 1, ptype := 1]

# label casual partnerships
plong2[p_once == 2 &
    (p_main_ong %in% c(NA, 0) &
     p_main_end %in% c(NA, 0)),
    ptype := 2]

# label one-time partnerships
plong2[p_once == 1, ptype := 3]

# remove filler rows and rename dataset
anl <- plong2[!is.na(ptype)]

str(anl)
anl[, length(unique(id))]

table(anl$ptype)

print(anl[, .(id,
              pid,
              ptype,
              p_ongoing,
              p_main_ong,
              p_once)])

anl[, .N, .(p_once,
            p_main_ong,
            p_main_end,
            ptype)] %>%
  .[order(ptype, p_once)] %>%
  .[, P := round(N / sum(N), 3)] %>%
  print

# create ongoing partnership indicator based on response to
# p_ongoing and p_main_ongoing
anl[ptype %in% 1:2,
    p_ongoing_ind := ifelse(p_ongoing == 1 | p_main_ong == 1, 1, NA)]

# regardless of reported intentions of ongoingness, partnerships with one
# sex act treated as one-time partnerships
anl[is.na(p_ongoing_ind), p_ongoing_ind := 0]

# check coding by ptype
anl[, .N, keyby = .(ptype, p_ongoing_ind)] %>%
  .[, P := N / sum(N), by = ptype] %>%
  print

# check main/casual ongoing indicator
anl[ptype %in% 1:2, .N, keyby = .(ptype,
                                  p_ongoing,
                                  p_main_ong,
                                  p_ongoing_ind)]


# create ongoing partnership indicator based on reported relationship
# start and end date
anl[ptype %in% 1:2, .N, .(p_startyyyy, p_startyyyydk)]
anl[, unique(p_startyyyy)]
anl[, unique(p_startyyyydk)]

# replace invalid years with NA
anl[ptype %in% 1:2, p_startyyyy :=
    ifelse(!grepl("^[19|20][0-9]{2}", p_startyyyy),
           NA, p_startyyyy)]

anl[ptype %in% 1:2, prop.table(table(is.na(p_startyyyy)))]
anl[ptype %in% 1:2, prop.table(table(is.na(p_startyyyydk)))]

# where startyyyy missing and startyyyydk estiamte present,
# impute partnership start date

set.seed(3000)

# capture year-month of submission date
anl[, sub_ym := floor_date(as.Date(sub_date), "month")]
class(anl$sub_ym)
ymd(as_date(anl$sub_ym[1:10]) - sample(1:365, size = 1))

anl[ptype %in% 1:2, .N,
    keyby = .(yyyy_miss = is.na(p_startyyyy),
              yyyydk_miss = is.na(p_startyyyydk),
              mm_miss = is.na(p_startmm))]

anl[ptype %in% 1:2, .N, keyby = .(ptype, p_ongoing_ind)]

anl[ptype %in% 1:2,
    date_impute_type := case_when(
      !is.na(p_startyyyy) & !is.na(p_startmm) ~ "none",
      !is.na(p_startyyyy) &  is.na(p_startmm) ~ "m",
       is.na(p_startyyyy) & !is.na(p_startyyyydk) & !is.na(p_startmm) ~ "y",
       is.na(p_startyyyy) & !is.na(p_startyyyydk) &  is.na(p_startmm) ~ "ym",
       TRUE ~ NA_character_
    ), .(id, pid)]

anl[ptype %in% 1:2, .N, keyby = .(p_ongoing_ind, date_impute_type)]

# FALSE-FALSE cell should equal 0
anl[, table(is.na(p_rai), is.na(p_rai_once))]
anl[, table(is.na(p_iai), is.na(p_iai_once))]
anl[, table(is.na(p_ioi), is.na(p_ioi_once))]
anl[, table(is.na(p_roi), is.na(p_roi_once))]
anl[, table(is.na(p_actsprefernot), is.na(p_actsprefernot_once))]
anl[, table(is.na(p_actsdk), is.na(p_actsdk_once))]


# %% SEXUAL ACTS PRESENT AND ANAL/ORAL SEX ROLE --------------------------------

# set partner subtype, main and casual partnerships
anl[ptype %in% 1:2,
    psubtype := dplyr::case_when(
      (p_rai == 1 | p_iai == 1) & (p_roi == 1 | p_ioi == 1) ~ "oa",
      (p_rai == 0 & p_iai == 0) & (p_roi == 1 | p_ioi == 1) ~ "oo",
      (p_rai == 1 | p_iai == 1) & (p_roi == 0 & p_ioi == 0) ~ "ao",
      p_rai + p_iai + p_roi + p_ioi +
        p_actsprefernot + p_actsdk == 0 ~ "nosex",
      (p_rai + p_iai + p_roi + p_ioi == 0) &
        (p_actsprefernot == 1 | p_actsdk == 1) ~ NA_character_,
      TRUE ~ "999999"
      )]

# set partner subtype, one-time contacts
anl[ptype == 3,
    psubtype := dplyr::case_when(
      (p_rai_once == 1 | p_iai_once == 1) &
        (p_roi_once == 1 | p_ioi_once == 1) ~ "oa",
      (p_rai_once == 0 & p_iai_once == 0) &
        (p_roi_once == 1 | p_ioi_once == 1) ~ "oo",
      (p_rai_once == 1 | p_iai_once == 1) &
        (p_roi_once == 0 & p_ioi_once == 0) ~ "ao",
      p_rai_once + p_iai_once + p_roi_once + p_ioi_once +
        p_actsprefernot_once + p_actsdk_once == 0 ~ "nosex",
      (p_rai_once + p_iai_once + p_roi_once + p_ioi_once == 0) &
        (p_actsprefernot_once == 1 | p_actsdk_once == 1) ~ NA_character_,
      TRUE ~ NA_character_
      )]

# set ego sex roles, main and casual partnerships
anl[ptype %in% 1:2, ":=" (
    ego.anal.role = dplyr::case_when(
      p_rai == 1 & p_iai == 0 ~ "Receptive",
      p_rai == 0 & p_iai == 1 ~ "Insertive",
      p_rai == 1 & p_iai == 1 ~ "Versatile",
      TRUE ~ NA_character_
    ),
    ego.oral.role = dplyr::case_when(
      p_roi == 1 & p_ioi == 0 ~ "Receptive",
      p_roi == 0 & p_ioi == 1 ~ "Insertive",
      p_roi == 1 & p_iai == 1 ~ "Versatile",
      TRUE ~ NA_character_
    )
)]

anl[ptype %in% 1:2, .(id, pid, ptype, psubtype,
                      p_rai, p_iai, p_roi, p_ioi,
                      ego.anal.role, ego.oral.role)]

sc <- names(anl)[grepl("once", names(anl))][2:7]
anl[, table(ptype, psubtype, exclude = NULL)]

egos <- length(unique(anl$id))

total_main <- unname(table(anl$ptype == 1)[2])
total_casl <- unname(table(anl$ptype == 2)[2])
total_pships <- anl[!is.na(ptype), length(unique(paste0(id, pid)))]
sum(total_main, total_casl)


cat(
  sep = "",
  "Unique egos\n", format(egos, big.mark = ","), "\n\n",
  "Unique main partnerships\n", format(total_main, big.mark = ","), "\n\n",
  "Unique casual partnerships\n", format(total_casl, big.mark = ","), "\n\n",
  "Overall partnerships\n", format(total_pships, big.mark = ",")
)

# ongoing main or casual partnerships
ong_cols <- c("id", "pid", "ptype", "psubtype",
              "startyyyy", "startyyyydk", "startmm",
              "ego.anal.role")

maincas_ong <- anl[ptype %in% 1:2 & p_ongoing_ind == 1, .N, ong_cols]

print(maincas_ong)
maincas_ong[, table(ptype, psubtype, exclude = NULL)]



# onetime partnerships
ponetime <- anl[ptype == 3, .N, .(id, psubtype, ptype)]

# drops 1 partnership with no sex
# drops 12 partnerships with missing data on psubtype
ong_counts <- dcast(maincas_ong[!psubtype %in% c("nosex", NA)],
                    id ~ psubtype + ptype,
                    value.var = "N")

print(ong_counts)


# %% IMPUTE AGE ----------------------------------------------------

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


# # %% IMPUTE RACE/ETHNICITY -----------------------------------------------------
#
# # check missing partner race/ethnicity, among egos with partnerships
# anl[!is.na(nn), sum(is.na(p_race.cat)) / unql(paste0(id, pid))]


# %% CALCULATE MAIN AND CASUAL DEGREE ------------------------------------------

# Main and causal partnerships (ongoing)

deg_data <- ong_counts[an[, .(id,
                              pn_ongoing,
                              race.cat,
                              age5)],
                       on = "id"][order(id)]

sum(is.na(deg_data$race.cat))
sum(is.na(deg_data$age5))

# if partnership type count missing, set to 0
deg_data[, ":=" (
  ao_1 = ifelse(is.na(ao_1), 0, ao_1),
  ao_2 = ifelse(is.na(ao_2), 0, ao_2),
  oa_1 = ifelse(is.na(oa_1), 0, oa_1),
  oa_2 = ifelse(is.na(oa_2), 0, oa_2),
  oo_1 = ifelse(is.na(oo_1), 0, oo_1),
  oo_2 = ifelse(is.na(oo_2), 0, oo_2)
  )] %>%
  .[, degmain := ao_1 + oa_1 + oo_1] %>%
  .[, degmain_trunc2 := ifelse(degmain >= 2, 2, degmain)] %>%
  .[, main_conc := ifelse(degmain >= 2, 1, 0)] %>%
  .[, degcasl := ao_2 + oa_2 + oo_2] %>%
  .[, casl_conc := ifelse(degcasl >= 2, 1, 0)]

print(deg_data)

deg_data[, .N, degmain][, P := N / sum(N)] %>% print
deg_data[, .N, main_conc][, P := N / sum(N)] %>% print
deg_data[, .N, keyby = degcasl][, P := N / sum(N)] %>% print
deg_data[, .N, casl_conc][, P := N / sum(N)] %>% print


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

# %% MAIN AND CASUAL DURATIONS -------------------------------------------------



# %% INSPECT DEGREE/ONETIME DATASETS -------------------------------------------

# oa = oral and anal sex partnership
# ao = anal-only partnership
# oo = oral-only partnership

deg_data[, .(mn_degmain = mean(degmain),
             mn_oa_1 = mean(oa_1),
             mn_ao_1 = mean(ao_1),
             mn_oo_1 = mean(oo_1))]

sum_main_total <- deg_data[, .N, key = degmain][, P := round(N / sum(N), 3)]
sum_main_oa    <- deg_data[, .N, key = oa_1][, P := round(N / sum(N), 3)]
sum_main_ao    <- deg_data[, .N, key = ao_1][, P := round(N / sum(N), 3)]
sum_main_oo    <- deg_data[, .N, key = oo_1][, P := round(N / sum(N), 3)]

sum_casl_total <- deg_data[, .N, key = degcasl][, P := round(N / sum(N), 3)]
sum_casl_oa    <- deg_data[, .N, key = oa_2][, P := round(N / sum(N), 3)]
sum_casl_ao    <- deg_data[, .N, key = ao_2][, P := round(N / sum(N), 3)]
sum_casl_oo    <- deg_data[, .N, key = oo_2][, P := round(N / sum(N), 3)]

pcols <- c("id",
           "pid",
           "ptype",
           "psubtype",
           "p_race.cat",
           "p_hiv",
           names(anl)[grep("^ego", names(anl))])

mp <- anl[ptype == 1 & p_ongoing_ind == 1, ..pcols]
cp <- anl[ptype == 2 & p_ongoing_ind == 1, ..pcols]

print(mp)
print(cp)

mp_racemix <- mp[, ctable(ego.race.cat, p_race.cat, prop = "n", useNA = "no")]
cp_racemix <- cp[, ctable(ego.race.cat, p_race.cat, prop = "n", useNA = "no")]

print(mp_racemix)
print(cp_racemix)

mp_racesummary <- mp[, .N, .(ego.race.cat, p_race.cat)] %>%
 .[, P := N / sum(N), .(ego.race.cat)] %>%
 na.omit

mp_racesummary[order(ego.race.cat, p_race.cat)] %>%
  .[, P := round(P, 3)] %>%
  print

ggplot(data = mp_racesummary,
       mapping = aes(x = ego.race.cat,
                     y = P,
                     group = p_race.cat,
                     fill = p_race.cat)) +
  geom_col(position = "dodge",
           col = "white",
           width = 0.3) +
  scale_fill_viridis_d(name = "Partner\nRace/ethnicity") +
  labs(x = "Ego Race/ethnicity",
       y = "Probability") +
  theme_classic()


# %% PARTNERSHIP REUSABLES -----------------------------------------------------

race_unq <- data.table(race.cat = sort(unique(deg_data$race.cat)))
age5_unq <- data.table(age5 = sort(unique(deg_data$age5)))

main_degt2_unq <- data.table(
  degmain_trunc2 = sort(unique(deg_data$degmain_trunc2))
)
casual_deg_unq <- data.table(degcasl = sort(unique(deg_data$degcasl)))

cilabs <- c("ll95", "ul95")

# %% MAIN PARTNERSHIPS ---------------------------------------------------------

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

main_concurrent_prob <- deg_data[, mean(main_conc)]
main_concurrent_prob


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
                                   fit_casl_deg_bymain,
                                   names = cilabs)
print(casl_deg_bymaindegt2_preds)

casl_concurrent_prob <- deg_data[, mean(casl_conc)]
casl_concurrent_prob


# %% WRITE DATA ----------------------------------------------------------------

# write summaries

main_summaries <- list(sum_main_total = sum_main_total,
                       sum_main_oa    = sum_main_oa,
                       sum_main_ao    = sum_main_ao,
                       sum_main_oo    = sum_main_oo)

casl_summaries <- list(sum_casl_total = sum_casl_total,
                       sum_casl_oa    = sum_casl_oa,
                       sum_casl_ao    = sum_casl_ao,
                       sum_casl_oo    = sum_casl_oo)

saveRDS(list(main_summaries = main_summaries,
             casl_summaries = casl_summaries),
        file = "netstats/degree_summaries.Rds")

# Main predictions

main_predictions <- list()

main_predictions$main_deg_byrace <- list(
  fits  = fit_main_deg_byrace,
  preds = main_deg_byrace_preds
)

main_predictions$main_deg_byage5 <- list(
  fits  = fit_main_deg_byage5,
  preds = main_deg_byage5_preds
)

main_predictions$main_deg_bycasltot <- list(
  fits  = fit_main_deg_bycasl,
  preds = main_deg_bycasldeg_preds
)

# Casual predictions

casl_predictions <- list()

casl_predictions$casl_deg_byrace <- list(
  fits = fit_casl_deg_byrace,
  preds = casl_deg_byrace_preds
)

casl_predictions$casl_deg_byage5 <- list(
  fits  = fit_casl_deg_byage5,
  preds = casl_deg_byage5_preds
)

casl_predictions$casl_deg_bymaindegt2 <- list(
  fits  = fit_casl_deg_bymaint2,
  preds = casl_deg_bymaindegt2_preds
)


pship_predictions <- list(main = main_predictions,
                          casl = casl_predictions)
pship_predictions

saveRDS(pship_predictions, file = "netstats/pship_predictions.Rds")
