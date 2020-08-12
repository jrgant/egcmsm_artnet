# %% Setup ---------------------------------------------------------------------

source("01-Import-Private-Data.R")

pacman::p_load(
  data.table,
  summarytools,
  ciTools,
  dplyr,
  stringr,
  lubridate,
  forcats
)

# function to calculate number of unique levels in a variable
unql <- function(data) length(unique(data))

# import cleaned wide dataset
an <- fread(paste0(Sys.getenv("ARTNET_PATH"), "/artnet-wide-cleaned.csv"))

# default ggplot2 theme
theme_set(theme_base())


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
#  - Data.table kicks a warning that NAs are introduced due to coercion.
#  - A result of converting a subset of variables to integer format.
#  - Some were saved as character with missing values recorded as "".

sapply(slugs[grepl("recai|recuai|insai|insuai|once|acts|unit", slugs)], print)

integer_coerce <-
  slugs[grepl("recai|recuai|insai|insuai|recoi|insoi|once|acts|unit",
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
lapply(melted_vars[integer_coerce], function(x) unique(unlist(x[, 3])))
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
plong2 <- plong[
  an[, .(
    id,
    sub_date,
    ego.age = age,
    ego.age.grp = age.grp,
    ego.race.cat = race.cat,
    ego.hiv = hiv.ego,
    ego.ongoing = pn_ongoing,
    prep_revised
  )],
  on = "id"][
    order(id, pid)][
    , nn := ifelse(nn == "", NA, nn)
    ]

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
                            stringr::str_extract(
                              names(plong2),
                              "(?<=^p_).*_once$")
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

plong2[p_main_ong == 1 | p_main_end == 1, ptype := 1] # label main partnerships
plong2[p_once == 2 & is.na(ptype), ptype := 2] # label casual partnerships
plong2[p_once == 1 & is.na(ptype), ptype := 3] # label one-time partnerships


# remove filler rows and rename dataset
anl <- plong2[!is.na(ptype)]

str(anl)
anl[, length(unique(id))]

table(anl$ptype)
table(is.na(anl$sub_date))

anl[, .(id, pid, ptype, p_ongoing, p_main_ong, p_once)][]

anl[, .N, .(p_once, p_main_ong, p_main_end, ptype)][
  order(ptype, p_once)][,
  P := round(N / sum(N), 3)][]

# create ongoing partnership indicator based on response to
# p_ongoing and p_main_ongoing
anl[ptype %in% 1:2,
    p_ongoing_ind := ifelse(p_ongoing == 1 | p_main_ong == 1, 1, NA)]

# regardless of reported intentions of ongoingness, partnerships with one
# sex act treated as one-time partnerships
anl[is.na(p_ongoing_ind), p_ongoing_ind := 0]

# check coding by ptype
anl[, .N, keyby = .(ptype, p_ongoing_ind)
   ][, P := N / sum(N), by = ptype][]

# check main/casual ongoing indicator
anl[ptype %in% 1:2, .N, keyby = .(ptype, p_ongoing, p_main_ong, p_ongoing_ind)]


# %% PARTNERSHIP DATES ---------------------------------------------------------

anl[ptype %in% 1:2, .N, .(p_startyyyy, p_startyyyydk)]
anl[, unique(p_startyyyy)]
anl[, unique(p_startyyyydk)]

# replace invalid years with NA
anl[ptype %in% 1:2,
    p_startyyyy := ifelse(
      !grepl("^(19|20)[0-9]{2}", p_startyyyy),
      NA, p_startyyyy
    )]

# check overall start dates
anl[ptype %in% 1:2, .N, .(startyyyy_miss = is.na(p_startyyyy))
   ][, P := N / sum(N)][]

anl[ptype %in% 1:2, .N, .(startyyyydk_miss = is.na(p_startyyyydk))
   ][, P := N / sum(N)][]

anl[ptype %in% 1:2, .N, .(startmm_miss = is.na(p_startmm))
   ][, P := N / sum(N)][]

anl[ptype %in% 1:2, .N, .(maincasendmm_miss = is.na(p_maincasendmm))
   ][, P := N / sum(N)][]

# where startyyyy missing and startyyyydk estimate present,
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
      !is.na(p_startyyyy) & is.na(p_startmm) ~ "month",
      is.na(p_startyyyy) & !is.na(p_startyyyydk) &
        !is.na(p_startmm) ~ "year",
      is.na(p_startyyyy) & !is.na(p_startyyyydk) &
        is.na(p_startmm) ~ "year and month",
      is.na(p_startyyyy) & is.na(p_startyyyydk) &
        !is.na(p_startmm) ~ "year",
      TRUE ~ NA_character_
    ), .(id, pid)]

anl[ptype %in% 1:2, .N, keyby = .(p_ongoing_ind,
                                  date_impute_type)]


# %% SEXUAL ACTS PRESENT AND ANAL/ORAL SEX ROLE --------------------------------

# FALSE-FALSE cell should equal 0
anl[, table(is.na(p_rai), is.na(p_rai_once))]
anl[, table(is.na(p_iai), is.na(p_iai_once))]
anl[, table(is.na(p_ioi), is.na(p_ioi_once))]
anl[, table(is.na(p_roi), is.na(p_roi_once))]
anl[, table(is.na(p_actsprefernot), is.na(p_actsprefernot_once))]
anl[, table(is.na(p_actsdk), is.na(p_actsdk_once))]

anl[ptype == 3, .N, keyby = .(p_rai, p_rai_once)]
anl[ptype == 3, .N, keyby = .(p_iai, p_iai_once)]
anl[ptype == 3, .N, keyby = .(p_roi, p_roi_once)]
anl[ptype == 3, .N, keyby = .(p_ioi, p_ioi_once)]

## Consolidate RAI, IAI, ROI, and IOI variables, following what Sam & co.
## did in the ARTnet data cleaning code.
anl[is.na(p_rai) & !is.na(p_rai_once), p_rai := p_rai_once]
anl[is.na(p_iai) & !is.na(p_iai_once), p_iai := p_iai_once]
anl[is.na(p_roi) & !is.na(p_roi_once), p_roi := p_roi_once]
anl[is.na(p_ioi) & !is.na(p_ioi_once), p_ioi := p_ioi_once]

## Recheck
anl[ptype == 3, .N, keyby = .(p_rai, p_rai_once)]
anl[ptype == 3, .N, keyby = .(p_iai, p_iai_once)]
anl[ptype == 3, .N, keyby = .(p_roi, p_roi_once)]
anl[ptype == 3, .N, keyby = .(p_ioi, p_ioi_once)]


## Create a partnership/contact subtype categorizing partners as:
##  + oral-only
##  + anal-only
##  + oral-anal
anl[,
    psubtype := dplyr::case_when(
    (p_rai == 1 | p_iai == 1) & (p_roi == 1 | p_ioi == 1) ~ "oralanal",
    (p_rai == 0 & p_iai == 0) & (p_roi == 1 | p_ioi == 1) ~ "oralonly",
    (p_rai == 1 | p_iai == 1) & (p_roi == 0 & p_ioi == 0) ~ "analonly",
    p_rai + p_iai + p_roi + p_ioi +
      p_actsprefernot + p_actsdk == 0 ~ "nosex",
    (p_rai + p_iai + p_roi + p_ioi == 0) &
      (p_actsprefernot == 1 | p_actsdk == 1 |
         p_actsprefernot_once == 1 | p_actsdk_once == 1) ~ NA_character_,
    TRUE ~ "999999")
    ]

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

anl[, .N, keyby = .(ptype, ego.anal.role, p_rai, p_iai)]

anl[, .(
  id, pid, ptype, psubtype,
  p_rai, p_iai, p_roi, p_ioi,
  ego.anal.role, ego.oral.role
)]

sc <- names(anl)[grepl("once", names(anl))][2:7]
anl[, table(ptype, psubtype, exclude = NULL)]


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

anl[, .N, .(p_race, p_hisp, p_race.cat)]
anl[, .(id, pid, nn, p_hisp, p_race, p_race.cat)]
anl[, .N, .(id)][, summary(N)]


# %% IMPUTE PARTNER AGE (FIRST PASS BASED ON RELAGE RESPONSE) ------------------

## NOTE: Source for age imputation:
## Weiss KM, Goodreau SM, Morris M, Prasad P, Ramaraju R, Sanchez T, et al.
## Egocentric Sexual Networks of Men Who Have Sex with Men in the United States:
## Results from the ARTnet Study. medRxiv. 2019. doi:10.1101/19010579

## NOTE:
## - I used a different method to impute these values than in the Weiss paper

plot(density(anl$p_age, na.rm = T))
anl[, .N, key = p_age]

# set age == 0 to missing
anl[, p_age := ifelse(p_age == 0, NA, p_age)]
anl[, .(id, p_age)][]

# 8.8% missing
total_propmissing_page <- anl[, sum(is.na(p_age)) / .N]
total_propmissing_page

# 7.0% missing
anl[is.na(p_age) & !is.na(p_relage), sum(.N)] / nrow(anl)
anl[, .N, key = p_relage]

class(anl$p_age)

# Based on PARTXRELAGE (ART-Net survey)

set.seed(1998)

# impute ages in partnerships where ego provided binned relative age of partner

anl[, p_age_imputed := dplyr::case_when(
          p_relage == 1 ~ ego.age - sample(size = 1, x = 10:15),
          p_relage == 2 ~ ego.age - sample(size = 1, x = 2:10),
          p_relage == 3 ~ ego.age + sample(size = 1, x = -1:1),
          p_relage == 4 ~ ego.age + sample(size = 1, x = 2:10),
          p_relage == 5 ~ ego.age + sample(size = 1, x = 10:15),
          !is.na(p_age) ~ p_age,
          TRUE ~ NA_integer_),
      .(id, pid)] %>%
  # set lower bound to minimum partner age reported in dataset
  .[, p_age_imputed := ifelse(p_age_imputed < 11, 11, round(p_age_imputed))] %>%
  .[, p_age.grp := dplyr::case_when(
          p_age_imputed <= 24 ~ 1,
          p_age_imputed >= 25 & p_age_imputed <= 34 ~ 2,
          p_age_imputed >= 35 & p_age_imputed <= 44 ~ 3,
          p_age_imputed >= 45 & p_age_imputed <= 54 ~ 4,
          p_age_imputed >= 55 ~ 5,
          TRUE ~ NA_real_
    )]

anl[, .(p_age, p_age_imputed)][]
anl[!is.na(p_age), summary(p_age)]
anl[, summary(p_age_imputed)]

missing_after_relage_impute <- anl[, sum(is.na(p_age_imputed))]
missing_after_relage_impute

anl[p_age_imputed < 15, .(
  id, pid, ego.age, p_age, p_relage, p_age_imputed
)]


# %% CALCULATE PARTNERSHIP DURATIONS -------------------------------------------

# convert subdate to date
anl[, sub_date := ymd(sub_date)]
str(anl)

set.seed(283798)

# impute start month
anl[ptype %in% 1:2 & p_ongoing_ind == 1 & is.na(p_startmm),
    p_startmm := sample(1:12, size = 1),
    by = .(id, pid)]

anl[ptype %in% 1:2 & p_ongoing_ind == 1, p_startdt := ymd(
  paste(
    p_startyyyy, p_startmm,
    # sample day of the month (all rows)
    ((p_startmm %in% c(1, 3, 5, 7, 8, 10, 12)) * sample(1:31, size = 1)) +
      ((p_startmm %in% c(4, 6, 9, 11)) * sample(1:30, size = 1)) +
      ((p_startmm == 2) * sample(1:28, size = 1)),
    sep = "-")),
  by = .(id, pid)]

anl[ptype %in% 1:2 & p_ongoing_ind == 1 & is.na(p_startdt)][]
anl[ptype %in% 1:2 & p_ongoing_ind == 1 & !is.na(p_startdt)][]

datecols <- c(
  "ptype",
  "p_startyyyy", "p_startyyyydk",
  "p_startmm", "p_startdt"
)

# if provided year range of partnership start date, impute
# set upper limit to 15
startdt_match <- list(ydk1 = 1:364,
                      ydk2 = 365 : (2 * 365 - 1),
                      ydk3 = (2 * 365) : (5 * 365 - 1),
                      ydk4 = (5 * 365) : (10 * 365 - 1),
                      ydk5 = (10 * 365) : (15 * 365 - 1))

sapply(startdt_match, range)

anl[!is.na(p_startyyyydk),
  p_startdt := sub_date - sample(startdt_match[[p_startyyyydk]], size = 1),
  by = .(id, pid)]

anl[ptype %in% 1:2 & p_ongoing_ind == 1 & !is.na(p_startdt), .(
    sub_date, p_startdt
    )]


## NOTE:
## Here we take the absolute value of the difference in dates to account
## for imputed dates that were after the sub date.

anl[, durat_days := abs(as.numeric(sub_date - p_startdt))]
anl[, durat_wks := round(durat_days / 7)]
anl[ptype == 3, durat_wks := 0]

class(anl$durat_days)
class(anl$durat_wks)

anl[, .(
  id, pid,
  ptype, ego.race.cat, ego.age.grp,
  p_startyyyy, p_startmm, p_startdt, sub_date,
  durat_days, durat_wks
)]


# %% REFORMAT HAVEN-LABELED UNIT VARIABLES -------------------------------------

## reformat sex act rate unit to calculate weekly sex act rates

# receptive anal sex time units
anl[, p_unitrai_fix := dplyr::case_when(
  p_unitrai %in% c(11647, 11685, 11688, 11691, 11694) ~ 4,  # monthly
  p_unitrai %in% c(11648, 11686, 11689, 11692, 11695) ~ 52, # yearly
  p_unitrai %in% c(11646, 11684, 11687, 11690, 11693) ~ 1,  # weekly
  p_unitrai %in% c(13061, 13062, 13063, 13064, 13065) ~ durat_wks,
  TRUE ~ NA_real_
  )]

# insertive anal sex time units
anl[, p_unitiai_fix := dplyr::case_when(
  p_unitiai %in% c(11655, 11700, 11703, 11706, 11801) ~ 4,  # monthly
  p_unitiai %in% c(11656, 11701, 11704, 11707, 11802) ~ 52, # yearly
  p_unitiai %in% c(11654, 11699, 11702, 11705, 11800) ~ 1,  # weekly
  p_unitiai %in% c(13066, 13069, 13072, 13074, 13077) ~ durat_wks, # overall
  TRUE ~ NA_real_
  )]

# receptive oral sex time units
anl[, p_unitroi_fix := dplyr::case_when(
  p_unitroi %in% c(11658, 11709, 11712, 11715, 11718) ~ 4,  # monthly
  p_unitroi %in% c(11659, 11710, 11713, 11716, 11719) ~ 52, # yearly
  p_unitroi %in% c(11657, 11708, 11711, 11714, 11717) ~ 1,  # weekly
  p_unitroi %in% c(13067, 13070, 13073, 13075, 13078) ~ durat_wks, # overall
  TRUE ~ NA_real_
  )]

# insertive oral sex time units
anl[, p_unitioi_fix := dplyr::case_when(
  p_unitioi %in% c(11661, 11733, 11736, 11739, 11742) ~ 4,  # monthly
  p_unitioi %in% c(11662, 11734, 11737, 11740, 11743) ~ 52, # yearly
  p_unitioi %in% c(11660, 11732, 11735, 11738, 11741) ~ 1,  # weekly
  p_unitroi %in% c(13068, 13071, 13076, 13079, 13080) ~ durat_wks, # overall
  TRUE ~ NA_real_
  )]

anl[ptype %in% 1:2, .N, keyby = p_unitrai_fix]
anl[ptype %in% 1:2, .N, keyby = p_unitiai_fix]
anl[ptype %in% 1:2, .N, keyby = p_unitroi_fix]
anl[ptype %in% 1:2, .N, keyby = p_unitioi_fix]

## Set units that equal 0 to 0.5 to avoid infinite values in rate calculations.
anl[ptype %in% 1:2 & p_rai == 1 & p_unitrai_fix == 0, p_unitrai_fix := 0.5]
anl[ptype %in% 1:2 & p_iai == 1 & p_unitiai_fix == 0, p_unitiai_fix := 0.5]
anl[ptype %in% 1:2 & p_roi == 1 & p_unitroi_fix == 0, p_unitroi_fix := 0.5]
# NOTE: p_unitioi_fix had no 0 time units

anl[ptype %in% 1:2 & p_rai == 1, summary(p_unitrai_fix)]
anl[ptype %in% 1:2 & p_iai == 1, summary(p_unitiai_fix)]
anl[ptype %in% 1:2 & p_roi == 1, summary(p_unitroi_fix)]
anl[ptype %in% 1:2 & p_ioi == 1, summary(p_unitioi_fix)]

## Calculate act rates within main and casual partnerships. After we translate
## the overall and unprotected act rates to the weekly scale, we use that
## rate ratio as the per-act probability of condomless sex.
anl[ptype %in% 1:2, ":=" (
  recai.rate = p_recai / p_unitrai_fix,
  insai.rate = p_insai / p_unitiai_fix,
  recoi.rate = p_recoi / p_unitroi_fix,
  insoi.rate = p_insoi / p_unitioi_fix,
  recuai.rate = p_recuai1 / p_unitrai_fix,
  insuai.rate = p_insuai1 / p_unitiai_fix
)][, incompatible.uai := ifelse(
       recuai.rate > recai.rate | insuai.rate > insai.rate,
       1, NA
)][recuai.rate > recai.rate, recuai.rate := NA
 ][insuai.rate > insai.rate, insuai.rate := NA
 ][, ":=" (
    recuai.prob = recuai.rate / recai.rate,
    insuai.prob = insuai.rate / insai.rate
 )][]

anl[ptype == 3, ":=" (
  recuai.prob = p_recuai_once,
  insuai.prob = p_insuai_once
)][]

# number of condom probabilities set to NA due to incompatible UAI response
# n = 7
anl[, sum(incompatible.uai, na.rm = TRUE)]

# calculate summary act rate (aggregate over insertive/receptive)
anl[, pid_unique := .I
   ][ptype %in% 1:2, ":=" (
       ai.rate = mean(c(recai.rate, insai.rate), na.rm = TRUE),
       oi.rate = mean(c(recoi.rate, insoi.rate), na.rm = TRUE)
    ), pid_unique]

# check missing
anl[ptype %in% 1:2, .N]
anl[ptype %in% 1:2, sum(is.na(ai.rate))]
anl[ptype %in% 1:2, sum(is.na(oi.rate))]

# check for infinite rate calculations
anl[
  ptype %in% 1:2,
  lapply(.SD, function(x) sum(x == Inf, na.rm = TRUE)),
  .SDcols = c(
    "recai.rate", "insai.rate",
    "recoi.rate", "insoi.rate",
    "recuai.rate", "insuai.rate"
  )]

# add variables for yearly act rates
anl[, ":="(
  ai.rate.52 = floor(ai.rate * 52),
  oi.rate.52 = floor(oi.rate * 52)
)]

# calculate proportion of sex acts condom-protected (used as per-act prob.)
anl[
  ptype %in% 1:2 & !is.na(recuai.prob) & !is.na(insuai.prob),
  cond.prob := 1 - ((recuai.prob + insuai.prob) / 2)
][]

anl[
  !is.na(recuai.prob) & is.na(insuai.prob),
  cond.prob := 1 - recuai.prob
][]

anl[
  is.na(recuai.prob) & !is.na(insuai.prob),
  cond.prob := 1 - insuai.prob
][]

anl[ptype %in% 1:2, sum(is.na(cond.prob)) / .N]


# %% EDGE CHARACTERISTICS ------------------------------------------------------

## make factors of race.cat and age.grp varables
convert2factor <- c(
  "ego.race.cat",
  "ego.age.grp",
  "p_race.cat",
  "p_age_imputed"
)

anl[, (convert2factor) := lapply(.SD, as.factor), .SDcols = convert2factor]


# %% PARTNER HIV STATUS --------------------------------------------------------

anl[p_hiv == 1, p_hiv2 := 1]  # HIV-positive
anl[p_hiv == 2, p_hiv2 := 0]  # HIV-negative
anl[p_hiv %in% 3:4 | is.na(p_hiv), p_hiv2 := 2]  # HIV status unknown

anl[, .N, .(p_hiv, p_hiv2)]
anl[ptype %in% 1:2, .N, keyby = p_hiv2]
anl[ptype == 3, .N, keyby = p_hiv2]

anl[, .N, ego.hiv]

anl[ego.hiv == 0 & p_hiv2 == 0, hiv.concord := 0] # both HIV-negative
anl[ego.hiv == 1 & p_hiv2 == 1, hiv.concord := 1] # both HIV-positive
anl[ego.hiv == 2 & p_hiv2 == 2, hiv.concord := 2] # both HIV unknown
anl[ego.hiv != p_hiv2, hiv.concord := 3] # HIV serodiscordant

anl[, .N, keyby = .(ego.hiv, p_hiv2, hiv.concord)]


# %% ART AND PREP USE RECODING ------------------------------------------------

## Respondent PrEP use
##
##   + Respondents were asked about current PrEP use only if they reported ever
##     having taken PrEP.
##
##   + Respondents were asked about ever taking PrEP only if the reported
##     being HIV-negative.

anl[ego.hiv %in% c(1, 2), prep_revised := 0]
anl[, .N, .(ego.hiv, prep_revised)]

anl[prep_revised == 0, an_prep_current := 0]
anl[, .N, .(prep_revised, an_prep_current)]

## PrEP use skip patterns (within partnerships):
##
##   + Only HIV-NEGATIVE EGOS reported on their own PrEP use during a
##     partnership (P_PREPUSE) or ever (PREP_REVISED). Set to 0.
##
##   + Only HIV-NEGATIVE or -UNKNOWN PARTNERS had information regarding their
##     PrEP use during a partnership. Set to 0.

## Recode ego PrEP use within partnership as ever/never:
##   - If respondent reported never being on PrEP, set ego PrEP use within
##     partnership to "Never".

anl[prep_revised == 0 | p_prepuse == 3 | ego.hiv %in% 1:2, p_prepuse2 := 0]
anl[p_prepuse %in% c(1, 2), p_prepuse2 := 1]  # any PrEP use in partnership
# anl[ego.hiv %in% 1:2, p_prepuse2 := 0] # see skip pattern summary

## Recode partner PrEP use within partnership as ever/never.
anl[p_prepuse_part %in% c(1, 2), p_prepuse_part2 := 1]
anl[p_prepuse_part == 3 | p_hiv2 == 1, p_prepuse_part2 := 0] # any PrEP use in partnership
# anl[p_hiv2 == 1, p_prepuse_part2 := 0]  # see skip pattern summary

anl[, .N, keyby = p_prepuse]
anl[, .N, keyby = p_prepuse_part]

## ART use skip patterns (within partnerships):
##
##   + Only HIV-POSITIVE EGOS reported their own ART use during a partnership.
##
##   + Only HIV-POSITIVE PARTNERS had information regarding their ART use during
##     a partnership.
##

## Recode ego ART use within partnership as ever/never.
anl[
  ego.hiv == 1,
  p_artuse_bin := ifelse(
    p_artuse %in% 1:2, 1, p_artuse
  )]

anl[p_artuse_bin == 3 | ego.hiv %in% c(0, 2), p_artuse_bin := 0]
# anl[ego.hiv %in% c(0, 2), p_artuse_bin := 0] # See skip pattern summary.

anl[, .N, keyby = .(p_artuse, p_artuse_bin)]

## Recode partner ART use within partnership as ever/never.
anl[
  p_hiv2 == 1,
  p_artuse_part_bin := ifelse(
    p_artuse_part %in% 1:2, 1, p_artuse_part
  )]

anl[p_artuse_part_bin == 3 | p_hiv2 %in% c(0, 2), p_artuse_part_bin := 0]
# anl[p_hiv2 %in% c(0, 2), p_artuse_part_bin := 0] # See skip pattern summary.

anl[, .N, keyby = .(p_artuse_part, p_artuse_part_bin)]


# %% CALCULATE AGE DIFFERENCES -------------------------------------------------

# calculate age differences
anl[, p_age_imputed := as.numeric(as.character(p_age_imputed))]

anl[, ":="(
  abs_agediff = abs(p_age_imputed - ego.age),
  abs_sqrt_agediff = abs(sqrt(p_age_imputed) - sqrt(ego.age))
)]

anl[, .(abs_agediff, abs_sqrt_agediff)]
sum(is.na(anl$abs_sqrt_agediff))


# %% VARIABLE CLASS FIXES ------------------------------------------------------

## Fix NA coding
anl[ego.anal.role == "", ego.anal.role := NA]
anl[p_race.cat == "", p_race.cat := NA]


# %% WRITE LONG DATASET --------------------------------------------------------

fwrite(anl, paste0(Sys.getenv("ARTNET_PATH"), "/artnet-long-cleaned.csv"))
