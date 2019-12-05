# %% Setup --------------------------------------------------------------------

# Running this source file should result in a printout of the variables
# included in each of the datasets the import script returns.
source("01-import-private-data.R")

pacman::p_load(
  magrittr,
  data.table,
  ggplot2,
  ggthemes,
  ggridges,
  viridis,
  haven,
  rlang,
  purrr,
  dplyr,
  rms
)

# View PrEP-related variables
names(av)[grep("prep", names(av))] %>% sort %>% print

# View Ego ART Use Variables
names(av)[grep("(?<!p)art(?!net)", names(av), perl = T)] %>%
  .[-grep("\\_part", .)] %>%
  sort %>%
  print

# View Partner ART Use Variables
names(av)[grep("(?<!p)art(?!net)", names(av), perl = T)] %>%
  .[grep("\\_part", .)] %>%
  sort %>%
  print

# View all objects in global environment
ls.str()

av_data
print(av)
names(av) %>% sort %>% print

# How many variables are haven-labelled?
sapply(av, function(x) is.labelled(x)) %>% table


# %% Select Variables --------------------------------------------------------

varselect <- c(
  # Demographic
  "id",
  "age",
  "race",
  "race.cat",
  # HIV
  "artnetevertest",
  "artnetrcntrslt",
  "artnetevrpos",
  "artnetstatus",
  # PrEP
  "prep_revised",
  "artnetprep_current",
  "prep_hivtestfreq",
  # STI testing
  "prep_stithroatfreq",
  "prep_stirectfreq",
  "prep_stiurethfreq",
  "stireg",
  names(av)[grep("stitest", names(av))],
  # Partner variables
  names(av)[grep("part", names(av))],
  # Sexual behavior
  names(av)[grep("^m_|mmconc", names(av))]
  )

print(varselect)

# select variable subset
avs <- av[, ..varselect] %>%
  # fix duplicate name error
  as_tibble(.name_repair = "minimal") %>%
  setDT

glimpse(avs)

# view summary of data frame
dfSummary(avs, graph.col = F) %>% print


# %% HIV Status --------------------------------------------------------

freq(avs$artnetevertest) %>% print

# Logic: if artnetevertest == 1
freq(avs$artnetrcntrslt) %>% print

# Logic: if artnetrcntrslt == 1, 3, 4, 88, 99
freq(avs$artnetevrpos) %>% print

# Logic: if artnetevrpos == 1 & artnetrcntrslt != 2
freq(avs$artnetstatus) %>% print

avs[, .N, .(artnetrcntrslt, artnetevrpos, artnetstatus)] %>% print

# create HIV status variable
avs[, hiv := case_when(
    !is.na(artnetstatus) ~ artnetstatus - 1,
    artnetevrpos == 1 ~ 1,
    artnetrcntrslt == 2 ~ 1,
    artnetrcntrslt == 1 ~ 0,
    TRUE ~ NA_real_)]

# drop inherited label
label(avs$hiv) <- ""
freq(avs$hiv) %>% print

# check
avs[, .N, keyby = .(hiv, artnetrcntrslt, artnetevrpos, artnetstatus)] %>%
  print

# %% PrEP Status --------------------------------------------------------------

freq(avs$prep_revised) %>% print

avs[prep_revised == 1, freq(artnetprep_current)] %>%
  print

avs[, an_prep_current := ifelse(
  !is.na(prep_revised) & is.na(artnetprep_current),
  0, artnetprep_current
  )]

freq(avs$an_prep_current) %>% print

# %% STI Testing Variables (Ever PrEP Users) -----------------------------------

avs[, .N, .(prep_revised, artnetprep_current)]

ep <- avs[prep_revised == 1]

ep[, .N, keyby = .(prep_stithroatfreq, prep_stiurethfreq, prep_stirectfreq)]
ep[, prep_any_missing_stifreq := ifelse(
       is.na(prep_stithroatfreq + prep_stiurethfreq + prep_stirectfreq), 1, 0)]

ep[, .N, keyby = .(prep_any_missing_stifreq,
                   prep_stiurethfreq,
                   prep_stithroatfreq,
                   prep_stirectfreq)]

ep[, .N, keyby = prep_any_missing_stifreq][, p := N / sum(N)] %>% print

# Non-PrEP-related STI Tests among Ever PrEP Users

ep[, .N, keyby = stitest_2yr_prep]

# check extreme response
ep[stitest_2yr_prep == 2000, .(id, age, race.cat,
                               stitest_2yr_prep,
                               stitest_2yr_sympt_prep,
                               stitest_2yr_notif_prep)]

# set extreme value to missing (rationale: impossible value)
avs$stitest_2yr_prep[avs$stitest_2yr_prep == 2000] <- NA
avs[prep_revised == 1, .N, keyby = stitest_2yr_prep]

# create categorical verstion of stitest_2yr_prep
avs[, stitest_2yr_pcat := ifelse(stitest_2yr_prep >= 7, "7+", stitest_2yr_prep)]
avs[prep_revised == 1, .N, keyby = stitest_2yr_pcat]

# proportion of STI tests sought due to presence of symptoms
avs[!is.na(stitest_2yr_prep), .N, keyby = stitest_2yr_sympt_prep]
avs[, stitest_2yr_psympt_pct := stitest_2yr_sympt_prep / stitest_2yr_prep]


# %% STI Testing Variables (Never PrEP Users) ----------------------------------

avs[prep_revised != 1 & !is.na(prep_revised), freq(stitest_2yr)]

# check extreme response
avs[stitest_2yr == 2015, .(id, age, race.cat,
                           stitest_2yr,
                           stitest_2yr_sympt,
                           stitest_2yr_notif)]

# set extreme value to missing (rationale: impossible value)
avs$stitest_2yr[avs$stitest_2yr == 2015] <- NA
avs[prep_revised != 1, .N, keyby = stitest_2yr]

# create categorical version of stitest_2yr
avs[, stitest_2yr_cat := ifelse(stitest_2yr >= 7, "7+", stitest_2yr)]
avs[prep_revised != 1, .N, keyby = stitest_2yr_cat]

# proportion of STI tests sought due to presence of symptoms
avs[prep_revised != 1, freq(stitest_2yr_sympt)]

avs[, stitest_2yr_sympt_pct := stitest_2yr_sympt / stitest_2yr]

# regular STI testing
avs[prep_revised != 1 & stitest_2yr > 0, freq(stireg)] %>% print

avs[stireg == 1, freq(stitestfreq)] %>% print
avs[, stitestfreq_cat := ifelse(stitestfreq %in% c(NA, 9), NA, stitestfreq)]
avs[stireg == 1, freq(stitestfreq_cat)]


# %% Sexual Behavior Variables ------------------------------------------------

# Naming convention for derivative variables
# p - partner, n - number, o - oral, a -  anal, Xm - past X months

names(av)[grep("m_", names(av))] %>% sort %>% print

# .. Oral and Anal Partners
avs[, .N, .(m_mp12oanum)]
boxplot(avs$m_mp12oanum)

# drop one impossible partner number
avs[, pnoa_12m := ifelse(m_mp12oanum > 7000, NA, m_mp12oanum)]
boxplot(avs$pnoa_12m)

# compare original with derived
cbind(orig = summary(avs$m_mp12oanum), new = summary(avs$pnoa_12m)) %>% print

# .. Ongoing Partnerships
cbind(
  gt1part_ongoing = summary(avs$mmconc),
  eq1part_ongoing = summary(avs$mmconc_onepart)
  ) %>%
  print

boxplot(avs$mmconc)
boxplot(avs$mmconc_onepart)

avs[, pn_ongoing := ifelse(!is.na(mmconc_onepart), mmconc_onepart, mmconc)]

summary(avs$pn_ongoing)

freq(avs$pn_ongoing) %>% print

# .. Anal-only Partners

# all anal partners
freq(avs$m_mp12anum2) %>% print
freq(avs$m_mp12anum2_onepart) %>% print

avs[, pna_12m := ifelse(!is.na(m_mp12anum2_onepart),
                        m_mp12anum2_onepart, m_mp12anum2)]

freq(avs$pna_12m) %>% print

# unprotected anal partners
freq(avs$m_mp12uanum2) %>% print
freq(avs$m_mp12uanum2_onepart) %>% print

avs[, pnua_12m := ifelse(!is.na(m_mp12uanum2_onepart),
                         m_mp12uanum2_onepart, m_mp12uanum2)]

freq(avs$pnua_12m) %>% print

# single-contact anal partners
freq(avs$m_mp12instanum2) %>% print


# %% Write Cleaned Dataset ---------------------------------------------------

fwrite(avs, "artnet-cleaned.csv", row.names = F)
