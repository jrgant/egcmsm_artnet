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

print(names(plong2))
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
table(is.na(anl$sub_date))

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


# %% PARTNERSHIP DATES ---------------------------------------------------------

anl[ptype %in% 1:2, .N, .(p_startyyyy, p_startyyyydk)]
anl[, unique(p_startyyyy)]
anl[, unique(p_startyyyydk)]

# replace invalid years with NA
anl[ptype %in% 1:2,
    p_startyyyy := ifelse(!grepl("^(19|20)[0-9]{2}", p_startyyyy),
                          NA, p_startyyyy)]

# check overall start dates
anl[ptype %in% 1:2, .N, .(startyyyy_miss = is.na(p_startyyyy))] %>%
  .[, P := N / sum(N)] %>%
  print

anl[ptype %in% 1:2, .N, .(startyyyydk_miss = is.na(p_startyyyydk))] %>%
  .[, P := N / sum(N)] %>%
  print

anl[ptype %in% 1:2, .N, .(startmm_miss = is.na(p_startmm))] %>%
  .[, P := N / sum(N)] %>%
  print

anl[ptype %in% 1:2, .N, .(maincasendmm_miss = is.na(p_maincasendmm))] %>%
  .[, P := N / sum(N)] %>%
  print

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

# set partner subtype, main and casual partnerships
anl[ptype %in% 1:2,
    psubtype := dplyr::case_when(
      (p_rai == 1 | p_iai == 1) & (p_roi == 1 | p_ioi == 1) ~ "oralanal",
      (p_rai == 0 & p_iai == 0) & (p_roi == 1 | p_ioi == 1) ~ "oralonly",
      (p_rai == 1 | p_iai == 1) & (p_roi == 0 & p_ioi == 0) ~ "analonly",
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
        (p_roi_once == 1 | p_ioi_once == 1) ~ "oralanal",
      (p_rai_once == 0 & p_iai_once == 0) &
        (p_roi_once == 1 | p_ioi_once == 1) ~ "oralonly",
      (p_rai_once == 1 | p_iai_once == 1) &
        (p_roi_once == 0 & p_ioi_once == 0) ~ "analonly",
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

# dfSummary(anl, plain.ascii = T, graph.col = F)


# %% IMPUTE RACE/ETHNICITY -----------------------------------------------------

# @NOTE 2020-01-28
# holding off for now on this; might not be necessary


# %% IMPUTE PARTNER AGE --------------------------------------------------------

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
  .[, p_age_imputed := ifelse(p_age_imputed < 11, 11, p_age_imputed)] %>%
  .[, p_age5 := dplyr::case_when(
          p_age_imputed <= 24 ~ 1,
          p_age_imputed >= 25 & p_age_imputed <= 34 ~ 2,
          p_age_imputed >= 35 & p_age_imputed <= 44 ~ 3,
          p_age_imputed >= 45 & p_age_imputed <= 54 ~ 4,
          p_age_imputed >= 55 ~ 5,
          TRUE ~ NA_real_
    )]

print(anl[, .(p_age, p_age_imputed)])

anl[!is.na(p_age), summary(p_age)]

anl[, summary(p_age_imputed)]

anl[p_age_imputed < 15, .(id, pid, ego.age, p_age, p_relage, p_age_imputed)]

anl[, .(agediff = p_age_imputed - ego.age),
      key = p_relage] %>%
  .[, .(min = min(agediff),
        max = max(agediff)),
      key = p_relage]

anl[, table(is.na(p_age_imputed))]



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

anl[ptype %in% 1:2 & p_ongoing_ind == 1 & is.na(p_startdt)] %>% print
anl[ptype %in% 1:2 & p_ongoing_ind == 1 & !is.na(p_startdt)] %>% print

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
                      ydk5 = (10 * 365) : (15 * 364 - 1))

sapply(startdt_match, range)

anl[!is.na(p_startyyyydk),
  p_startdt := sub_date - sample(startdt_match[[p_startyyyydk]], size = 1),
  by = .(id, pid)]

anl[ptype %in% 1:2 & p_ongoing_ind == 1 & !is.na(p_startdt), .(sub_date, p_startdt)]

# @TODO 2020-01-28: fix this issue more naturally
# set abs to account for imputed dates that were after the sub date
anl[, durat_days := abs(as.numeric(sub_date - p_startdt))]
anl[, durat_wks := round(durat_days / 7)]
class(anl$durat_days)
class(anl$durat_wks)

anl[ptype %in% 1:2, .(
  id, pid,
  ptype, ego.race.cat, ego.age5,
  p_startyyyy, p_startmm, p_startdt, sub_date,
  durat_days, durat_wks
)]


# %% REFORMAT HAVEN-LABELED UNIT VARIABLES -------------------------------------

## reformat sex act rate unit to calculate weekly rates in EpiStats.R

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

# set units that equal 0 to 0.5 to avoid infinite values in rate calculations

anl[ptype %in% 1:2 & p_rai == 1 & p_unitrai_fix == 0, p_unitrai_fix := 0.5]
anl[ptype %in% 1:2 & p_iai == 1 & p_unitiai_fix == 0, p_unitiai_fix := 0.5]
anl[ptype %in% 1:2 & p_roi == 1 & p_unitroi_fix == 0, p_unitroi_fix := 0.5]
# p_unitioi_fix had no 0 time units

anl[ptype %in% 1:2 & p_rai == 1, summary(p_unitrai_fix)]
anl[ptype %in% 1:2 & p_iai == 1, summary(p_unitiai_fix)]
anl[ptype %in% 1:2 & p_roi == 1, summary(p_unitroi_fix)]
anl[ptype %in% 1:2 & p_ioi == 1, summary(p_unitioi_fix)]

# %% WRITE LONG DATASET --------------------------------------------------------

fwrite(anl, paste0(Sys.getenv("ARTNET_PATH"), "/artnet-long-cleaned.csv"))
