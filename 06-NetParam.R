# %% SETUP ---------------------------------------------------------------------

source("01-Import-Private-Data.R")

pacman::p_load(
  data.table,
  summarytools,
  ciTools,
  dplyr,
  stringr,
  lubridate,
  mice,
  MASS,
  rms,
  nnet
)

# function to calculate number of unique levels in a variable
unql <- function(data) length(unique(data))

# path to private data
pd_path <- Sys.getenv("ARTNET_PATH")

# import cleaned wide dataset
an <- fread(file.path(pd_path, "/artnet-wide-cleaned.csv"))
anl <- fread(file.path(pd_path, "/artnet-long-cleaned.csv"))

# import imputed datasets
imp_mc <- readRDS(file.path(pd_path, "artnet-imputed-mc-augmented.Rds"))
imp_otp <- readRDS(file.path(pd_path, "artnet-imputed-otp-augmented.Rds"))


################################################################################
                         ## SET UP NON-IMPUTED DATA ##
################################################################################

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
ong_cols <- c(
  "sub_date", "id", "pid",
  "ptype", "psubtype",
  "p_startyyyy", "p_startyyyydk", "p_startmm", "p_startdt",
  "durat_days", "durat_wks",
  "ego.race.cat", "ego.age.grp", "ego.age", "ego.hiv", "hiv2",
  "ego.anal.role", "hiv.concord"
)

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
ong_ptype_cts <- dcast(
  data = maincas_ong[!psubtype %in% c("no sex", NA)],
  formula = id ~ psubtype + ptype,
  fun.aggregate = sum,
  value.var = "N"
)

print(ong_ptype_cts)


################################################################################
                          ## PARTNERSHIP DURATIONS ##
################################################################################

maincas_ong[, .(
  id, pid, ego.race.cat, ego.age, ego.age.grp, ego.hiv,
  ptype, p_startyyyy, p_startmm, p_startdt, sub_date,
  durat_days, durat_wks
)]

maincas_ong[, table(durat_wks < 0)]

ggplot(maincas_ong, aes(x = durat_wks)) +
  geom_histogram(col = "black", fill = "aliceblue") +
  facet_wrap(~ ptype, nrow = 2) +
  theme_base()

# 18 missing partnership durations (among ongoing)
maincas_ong[, .(
  mn = mean(durat_wks, na.rm = TRUE),
  var = var(durat_wks, na.rm = TRUE),
  missing_durat = sum(is.na(durat_wks)),
  pct_missing_durat = sum(is.na(durat_wks)) / .N
), ptype]

fit_main_dur <- maincas_ong[
  ptype == 1,
  glm.nb(durat_wks ~ ego.race.cat + factor(ego.age.grp) + factor(hiv2))
]

fit_main_dur_ix <- update(
  fit_main_dur, . ~ + ego.race.cat * factor(ego.age.grp) + factor(hiv2)
)

## No support for interaction between race and age.
summary(fit_main_dur_ix)

pmdur <- as.data.table(
  expand.grid(
    ego.race.cat = unique(maincas_ong$ego.race.cat),
    hiv2 = c(0, 1),
    ego.age.grp = 1:5
  )
)

pmdur[, pred := predict(fit_main_dur, newdata = pmdur, type = "response")]

pmdur %>%
  ggplot(aes(x = ego.age.grp, y = pred, color = ego.race.cat)) +
  geom_line(size = 1) +
  geom_point(shape = 21, fill = "white", size = 5) +
  facet_wrap(
    ~ hiv2,
    labeller = labeller(c("HIV undiagnosed" = 0, "HIV diagnosed" = 1))
  ) +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  ggtitle("Predicted ongoing partnership durations") +
  labs(subtitle = "Main partnerships") +
  theme_clean()


## Casual partnerships

fit_casl_dur <- maincas_ong[
  ptype == 2,
  glm.nb(durat_wks ~ ego.race.cat + factor(ego.age.grp) + factor(hiv2))
]

## No support for interaction between race and age.
fit_casl_dur_ix <- update(
  fit_casl_dur, . ~ ego.race.cat * factor(ego.age.grp) + factor(hiv2)
)

summary(fit_casl_dur_ix)

pcdur <- as.data.table(
  expand.grid(
    ego.race.cat = unique(maincas_ong$ego.race.cat),
    hiv2 = c(0, 1),
    ego.age.grp = 1:5
  )
)

pcdur[, pred := predict(fit_casl_dur, newdata = pcdur, type = "response")]

pcdur %>%
  ggplot(aes(x = ego.age.grp, y = pred, color = ego.race.cat)) +
  geom_line(size = 1) +
  geom_point(shape = 21, fill = "white", size = 5) +
  facet_wrap(~ hiv2) +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  ggtitle("Predicted ongoing partnership durations") +
  labs(subtitle = "Casual partnerships") +
  theme_clean()


################################################################################
                    ## MAKE DATASET TO CALCULATE DEGREES ##
################################################################################

# Main and causal partnerships (ongoing)

# join ongoing partner counts with ego traits
deg_data <- ong_ptype_cts[
  an[, .(id, pn_ongoing, hiv.ego, race.cat, age.grp)],
  on = "id"
]

setkey(deg_data, id)

# should be 0
sum(is.na(deg_data$race.cat))
sum(is.na(deg_data$age.grp))

# if partnership type count missing, set to 0
deg_data[, ":="(
  analonly_1 = ifelse(is.na(analonly_1), 0, analonly_1), # main, anal-only
  oralanal_1 = ifelse(is.na(oralanal_1), 0, oralanal_1), # main, oral-anal
  oralonly_1 = ifelse(is.na(oralonly_1), 0, oralonly_1), # main, oral-only
  analonly_2 = ifelse(is.na(analonly_2), 0, analonly_2), # casual, anal-only
  oralanal_2 = ifelse(is.na(oralanal_2), 0, oralanal_2), # casual, oral-anal
  oralonly_2 = ifelse(is.na(oralonly_2), 0, oralonly_2)  # casual, oral-only
 )][, degmain := analonly_1 + oralanal_1 + oralonly_1
  ][, degmain_trunc2 := ifelse(degmain >= 2, 2, degmain)
  ][, degcasl := analonly_2 + oralanal_2 + oralonly_2
  ][, degtotal := degmain + degcasl
  ][, main_conc_ind := ifelse(degmain >= 2, 1, 0)
  ][, casl_conc_ind := ifelse(degcasl >= 2, 1, 0)
  ][, maincasl_conc_ind := degtotal > 1]

str(deg_data)

deg_data[, .N, degmain][, P := N / sum(N)][]
degmain_t2 <- deg_data[, .N, degmain_trunc2][, P := N / sum(N)]
print(degmain_t2)

deg_data[, .N, main_conc_ind][, P := N / sum(N)][]
deg_data[, .N, keyby = degcasl][, P := N / sum(N)][]
deg_data[, .N, casl_conc_ind][, P := N / sum(N)][]
deg_data[, .N, maincasl_conc_ind][, P := N / sum(N)][]

names(deg_data)

# create dummy variables for main and casl degrees
deg_data[, ":="(
  casl0 = dplyr::case_when(
    degcasl == 0 ~ 1,
    degcasl %in% 1:5 ~ 0,
    TRUE ~ NA_real_
  ),
  casl1 = dplyr::case_when(
    degcasl == 1 ~ 1,
    degcasl %in% c(0, 2:5) ~ 0,
    TRUE ~ NA_real_
  ),
  casl2 = dplyr::case_when(
    degcasl == 2 ~ 1,
    degcasl %in% c(0:1, 3:5) ~ 0,
    TRUE ~ NA_real_
  ),
  casl3 = dplyr::case_when(
    degcasl == 3 ~ 1,
    degcasl %in% c(0:2, 4:5) ~ 0,
    TRUE ~ NA_real_
  ),
  casl4 = dplyr::case_when(
    degcasl == 4 ~ 1,
    degcasl %in% c(0:3, 5) ~ 0,
    TRUE ~ NA_real_
  ),
  casl5 = dplyr::case_when(
    degcasl == 5 ~ 1,
    degcasl %in% 0:4 ~ 0,
    TRUE ~ NA_real_
  ),
  main0 = dplyr::case_when(
    degmain_trunc2 == 0 ~ 1,
    degmain_trunc2 %in% 1:2 ~ 0,
    TRUE ~ NA_real_
  ),
  main1 = dplyr::case_when(
    degmain_trunc2 == 1 ~ 1,
    degmain_trunc2 %in% c(0, 2) ~ 0,
    TRUE ~ NA_real_
  ),
  main2 = dplyr::case_when(
    degmain_trunc2 == 2 ~ 1,
    degmain_trunc2 %in% 0:1 ~ 0,
    TRUE ~ NA_real_
  )
)]

print(deg_data)


# %% INSPECT DEGREE/ONETIME DATASETS -------------------------------------------

# oa = oral and anal sex partnership
# ao = anal-only partnership
# oo = oral-only partnership

deg_data[, .(
  mn_degmain = mean(degmain),
  mn_oa_1 = mean(oralanal_1),
  mn_ao_1 = mean(analonly_1),
  mn_oo_1 = mean(oralonly_1)
)]

sum_main_total <- deg_data[, .N, key = degmain][, P := round(N / sum(N), 3)]
sum_main_oa    <- deg_data[, .N, key = oralanal_1][, P := round(N / sum(N), 3)]
sum_main_ao    <- deg_data[, .N, key = analonly_1][, P := round(N / sum(N), 3)]
sum_main_oo    <- deg_data[, .N, key = oralonly_1][, P := round(N / sum(N), 3)]

sum_casl_total <- deg_data[, .N, key = degcasl][, P := round(N / sum(N), 3)]
sum_casl_oa    <- deg_data[, .N, key = oralanal_2][, P := round(N / sum(N), 3)]
sum_casl_ao    <- deg_data[, .N, key = analonly_2][, P := round(N / sum(N), 3)]
sum_casl_oo    <- deg_data[, .N, key = oralonly_2][, P := round(N / sum(N), 3)]

pcols <- c(
  "id",
  "pid",
  "ptype",
  "psubtype",
  "p_race.cat",
  "p_age.grp",
  "p_hiv",
  "hiv.concord",
  names(anl)[grep("^ego", names(anl))]
)


################################################################################
                                ## RACE MATCH ##
################################################################################

mcdat <- as.data.table(complete(imp_mc, action = "long", include = TRUE))
str(mcdat)

age_breaks.i <- c(15, 25, 35, 45, 55, 66)
age_breaks.j <- c(14, 25, 35, 45, 55, 88)

# drop partners not aged between [18, 65]
mcdat[, ":=" (
  samerace = ifelse(race.i == race.j, 1, 0),
  age.grp.i = as.numeric(cut(age.i, age_breaks.i, right = FALSE)),
  age.grp.j = as.numeric(cut(age.j, age_breaks.j, right = FALSE)),
  serodisc = ifelse(hiv.concord == 2, 1, 0)
)][, sameage := ifelse(age.grp.i == age.grp.j, 1, 0)]

mcdat[.imp > 0, .N, keyby = .(race.combo, samerace)]
mcdat[.imp > 0, .(min = min(age.i), max = max(age.i)), keyby = age.grp.i]
mcdat[.imp > 0, .(min = min(age.j), max = max(age.j)), keyby = age.grp.j]
mcdat[.imp > 0, .N, keyby = .(sameage, age.grp.i, age.grp.j)]
mcdat[.imp > 0, .N, keyby = .(hiv.concord, serodisc)]

imp_mc <- as.mids(mcdat)

fit_racematch <- with(
  imp_mc,
  glm(
    samerace ~ ptype + race.i + factor(age.grp.i) + factor(diag.status.i),
    family = binomial,
    subset = age.j >= 18 & age.j <= 65
  )
)

pool_racematch <- pool(fit_racematch)
summary(pool_racematch)

## Save pooled estimates to a model object. USE ONLY FOR PREDICTING
## PROBABILITY OF SAME-RACE PARTNERSHIPS!!
fit_racematch_out <- fit_racematch$analyses[[1]]
fit_racematch_out$coefficients <- pool_racematch$pooled$estimate
names(fit_racematch_out$coefficients) <- pool_racematch$pooled$term

prmatch <- data.table(
  expand.grid(
    ptype = 1:2,
    age.grp.i = 1:5,
    race.i = unique(mcdat$race.i),
    diag.status.i = 0:1
  )
)

prmatch[, pred := predict(
            fit_racematch_out,
            newdata = prmatch,
            type = "response"
          )]

prmatch %>%
  ggplot(
    aes(x = age.grp.i, y = pred,
        color = race.i,
        linetype = factor(diag.status.i))
  ) +
  geom_line(size = 1) +
  geom_point(shape = 21, size = 5, fill = "white") +
  facet_wrap(~ ptype) +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  theme_clean()


################################################################################
                            ## AGE GROUP MATCHING ##
################################################################################

fit_agematch <- with(
  imp_mc,
  glm(
    sameage ~ ptype + race.i + factor(age.grp.i) + factor(diag.status.i),
    family = binomial,
    subset = age.j >= 18 & age.j <= 65
  )
)

pool_agematch <- pool(fit_agematch)
summary(pool_agematch)

fit_agematch_out <- fit_agematch$analyses[[1]]
fit_agematch_out$coefficients <- pool_agematch$pooled$estimate
names(fit_agematch_out$coefficients) <- pool_agematch$pooled$term

pamatch <- as.data.table(
  expand.grid(
    ptype = 1:2,
    race.i = unique(mcdat$race.i),
    age.grp.i = 1:5,
    diag.status.i = 0:1
  )
)

pamatch[, pred := predict(
            fit_agematch_out,
            newdata = pamatch,
            type = "response"
          )]

pamatch %>%
  ggplot(aes(x = age.grp.i, y = pred, color = race.i)) +
  geom_line(aes(linetype = factor(diag.status.i))) +
  geom_point(shape = 21, size = 5, fill = "white") +
  facet_wrap(~ ptype) +
  ggtitle("Predicted probability of same-age partnerships") +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  theme_clean()


################################################################################
                             ## HIV CONCORDANCE  ##
################################################################################

fit_serodisc <- with(
  imp_mc,
  glm(
    serodisc ~ ptype + race.i + age.grp.i + factor(diag.status.i),
    family = binomial,
    subset = age.j >= 18 & age.j <= 65
  )
)

pool_serodisc <- pool(fit_serodisc)
summary(pool_serodisc)

fit_serodisc_out <- fit_serodisc$analyses[[1]]
fit_serodisc_out$coefficients <- pool_serodisc$pooled$estimate
names(fit_serodisc_out$coefficients) <- pool_serodisc$pooled$term

pserodisc <- as.data.table(
  expand.grid(
    ptype = 1:2,
    race.i = unique(mcdat$race.i),
    age.grp.i = 1:5,
    diag.status.i = 0:1
  )
)

pserodisc[, pred := predict(
              fit_serodisc_out,
              newdata = pserodisc,
              type = "response"
            )]

pserodisc %>%
  ggplot(aes(x = age.grp.i, y = pred, color = race.i)) +
  geom_line(aes(linetype = factor(diag.status.i))) +
  geom_point(shape = 21, size = 5, fill = "white") +
  facet_wrap(~ ptype) +
  ggtitle("Predicted probability of HIV-serodiscordant partnership") +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  theme_clean()



################################################################################
                               ## EGO SEX ROLE ##
################################################################################

# Overall (among ids with main or casual partnerships)
# Main Partnerships

mp <- maincas_ong[ptype == 1]
mp[ego.anal.role == "", ego.anal.role := "Missing"]
mp_analrole <- dcast(
  mp,
  id + ego.race.cat + ego.age.grp + ego.hiv ~ ego.anal.role
)

## Create a role class variable based on available responses.
mp_analrole[, role.class := fcase(
    Versatile > 0 | (Insertive == 1 & Receptive == 1), "Versatile",
    Receptive > 0 & (Insertive == 0 & Versatile == 0), "Receptive",
    Insertive > 0 & (Receptive == 0 & Versatile == 0), "Insertive",
    default = NA_character_
  )]

mp_analrole[, .N, role.class][, p := N / sum(N)][]
mp_analrole[, missing := ifelse(is.na(role.class), 1, 0)]

ipsmod <- glm(
  missing ~ ego.race.cat + factor(ego.age.grp) + ego.hiv,
  data = mp_analrole,
  family = binomial
)

mp_analrole[, ips_den := predict(ipsmod, mp_analrole, type = "response")]
mp_analrole[missing == 0, ips_den := 1 - ips_den]
mp_analrole[,
  ips_num := fcase(
    missing == 1, sum(missing == 1) / nrow(mp_analrole),
    missing == 0, 1 - (sum(missing == 1) / nrow(mp_analrole))
  )]
mp_analrole[, ipw := ips_num / ips_den]

mp_analrole[, summary(ipw)]
boxplot(mp_analrole$ipw)

fit_mrole <- multinom(
  role.class ~ ego.race.cat + factor(ego.age.grp),
  data = mp_analrole
)

fit_mrole_ipw <- multinom(
  role.class ~ ego.race.cat + factor(ego.age.grp),
  weights = ipw,
  data = mp_analrole
)

pmrole <- as.data.table(
  expand.grid(
    ego.race.cat = unique(mp_analrole$ego.race.cat),
    ego.age.grp = 1:5
  )
)

mrole_preds <- predict(fit_mrole, newdata = pmrole, type = "probs")
mrole_preds_ipw <- predict(fit_mrole_ipw, newdata = pmrole, type = "probs")

pmrole[, ":="(
  Insertive = mrole_preds[, 1],
  Receptive = mrole_preds[, 2],
  Versatile = mrole_preds[, 3],
  Insertive_ipw = mrole_preds_ipw[, 1],
  Receptive_ipw = mrole_preds_ipw[, 2],
  Versatile_ipw = mrole_preds_ipw[, 3]
)]

pmrole

pmrole_melt <- melt(
  pmrole,
  id.vars = c("ego.race.cat", "ego.age.grp"),
  measure.vars = c(
    "Insertive", "Receptive", "Versatile",
    "Insertive_ipw", "Receptive_ipw", "Versatile_ipw"
  )
)

pmrole_melt[, ":="(
  model = fcase(grepl("ipw", variable), "ipw", default = "cc"),
  role.class = fcase(
    grepl("Insertive", variable), "Insertive",
    grepl("Receptive", variable), "Receptive",
    grepl("Versatile", variable), "Versatile"
  )
)]

pmrole_melt %>%
  ggplot(aes(x = ego.age.grp, y = value, color = ego.race.cat)) +
  geom_line(aes(linetype = model)) +
  geom_point(shape = 21, size = 5, fill = "white") +
  facet_wrap(~ role.class) +
  ggtitle("Anal role class probabilities") +
  labs(subtitle = "Main partnerships", y = "Membership probability") +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  theme_clean()


# Casual Partnerships

cp <- maincas_ong[ptype == 2]
cp[ego.anal.role == "", ego.anal.role := "Missing"]
cp_analrole <- dcast(
  cp,
  id + ego.race.cat + ego.age.grp + ego.hiv ~ ego.anal.role
)

## Create a role class variable based on available responses.
cp_analrole[, role.class := fcase(
    Versatile > 0 | (Insertive == 1 & Receptive == 1), "Versatile",
    Receptive > 0 & (Insertive == 0 & Versatile == 0), "Receptive",
    Insertive > 0 & (Receptive == 0 & Versatile == 0), "Insertive",
    default = NA_character_
  )]

cp_analrole[, .N, role.class][, p := N / sum(N)][]
cp_analrole[, missing := ifelse(is.na(role.class), 1, 0)]

ipsmod <- glm(
  missing ~ ego.race.cat + factor(ego.age.grp) + ego.hiv,
  data = cp_analrole,
  family = binomial
)

cp_analrole[, ips_den := predict(ipsmod, cp_analrole, type = "response")]
cp_analrole[missing == 0, ips_den := 1 - ips_den]
cp_analrole[,
  ips_num := fcase(
    missing == 1, sum(missing == 1) / nrow(cp_analrole),
    missing == 0, 1 - (sum(missing == 1) / nrow(cp_analrole))
  )]
cp_analrole[, ipw := ips_num / ips_den]

cp_analrole[, summary(ipw)]
boxplot(cp_analrole$ipw)

fit_crole <- multinom(
  role.class ~ ego.race.cat + factor(ego.age.grp),
  data = cp_analrole
)

fit_crole_ipw <- multinom(
  role.class ~ ego.race.cat + factor(ego.age.grp),
  weights = ipw,
  data = cp_analrole
)

pcrole <- as.data.table(
  expand.grid(
    ego.race.cat = unique(cp_analrole$ego.race.cat),
    ego.age.grp = 1:5
  )
)

crole_preds <- predict(fit_crole, newdata = pcrole, type = "probs")
crole_preds_ipw <- predict(fit_crole_ipw, newdata = pcrole, type = "probs")

pcrole[, ":="(
  Insertive = crole_preds[, 1],
  Receptive = crole_preds[, 2],
  Versatile = crole_preds[, 3],
  Insertive_ipw = crole_preds_ipw[, 1],
  Receptive_ipw = crole_preds_ipw[, 2],
  Versatile_ipw = crole_preds_ipw[, 3]
)]

pcrole

pcrole_melt <- melt(
  pcrole,
  id.vars = c("ego.race.cat", "ego.age.grp"),
  measure.vars = c(
    "Insertive", "Receptive", "Versatile",
    "Insertive_ipw", "Receptive_ipw", "Versatile_ipw"
  )
)

pcrole_melt[, ":="(
  model = fcase(grepl("ipw", variable), "ipw", default = "cc"),
  role.class = fcase(
    grepl("Insertive", variable), "Insertive",
    grepl("Receptive", variable), "Receptive",
    grepl("Versatile", variable), "Versatile"
  )
)]

pcrole_melt %>%
  ggplot(aes(x = ego.age.grp, y = value, color = ego.race.cat)) +
  geom_line(aes(linetype = model)) +
  geom_point(shape = 21, size = 5, fill = "white") +
  facet_wrap(~ role.class) +
  ggtitle("Anal role class probabilities") +
  labs(subtitle = "Main partnerships", y = "Membership probability") +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  theme_clean()


################################################################################
                  ## SEED PROBABILITIES, MAIN PARTNERSHIPS ##
################################################################################

## ... PROBS FOR SEEDING MAIN DEGREE

# main degree dummies
main_bin <- names(deg_data)[grepl("main[0-2]{1}", names(deg_data))]

## Because main degree will be seeded first, leave degcasl out. Seed casual
## degree will be generated subsequently, conditional on demographics
## and main degree.
fit_mdeg <- multinom(
  degmain_trunc2 ~ race.cat + factor(age.grp) + hiv.ego,
  data = deg_data
)

summary(fit_mdeg)
coef(fit_mdeg)
confint(fit_mdeg)

pmdeg <- as.data.table(
  expand.grid(
    race.cat = unique(deg_data$race.cat),
    age.grp = 1:5,
    hiv.ego = 0:1,
    degcasl = unique(deg_data$degcasl)
  )
)

mdeg_probs <- as.data.table(predict(fit_mdeg, newdata = pmdeg, type = "probs"))
names(mdeg_probs) <- paste0("maindeg", 0:2)
pmdeg <- cbind(pmdeg, mdeg_probs)

melt(pmdeg, measure.vars = c("maindeg0", "maindeg1", "maindeg2")) %>%
  ggplot(aes(x = age.grp, y = value, color = race.cat)) +
  geom_line(aes(linetype = factor(hiv.ego))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  facet_wrap(~ factor(degcasl) + variable) +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  theme_clean()


################################################################################
                      ## CONCURRENCY, MAIN PARTNERSHIPS ##
################################################################################

# proportion of individuals with concurrent main partnerships

deg_data[, .(id, race.cat, age.grp, degmain, degcasl, main_conc_ind)][]
deg_data[, degcasl3 := ifelse(degcasl >= 3, 3, degcasl)]

# sparse at degcasl > 3 (group for the target stat estimation)
fit_main_concurrent <- glm(
  as.numeric(main_conc_ind) ~ race.cat + factor(age.grp) + hiv.ego + factor(degcasl3),
  family = "binomial",
  data = deg_data
)

summary(fit_main_concurrent)
coef(fit_main_concurrent)
confint(fit_main_concurrent)


################################################################################
                     ## SEED PROBABILITY, CASUAL DEGREE ##
################################################################################

## ... PROBS FOR SEEDING CASUAL DEGREE

# casual degree dummies
casl_bin <- names(deg_data)[grepl("^casl[0-5]{1}", names(deg_data))]

# used to estimate the distribution of casual degree
# marginalized over race and age (used for seeding population)

## ... EXPECTED CASL  DEGREE (BY RACE, AGE, HIV STATUS, AND MAIN DEGREE)

fit_cdeg <- multinom(
  degcasl ~ race.cat + factor(age.grp) + hiv.ego + factor(degmain_trunc2),
  data = deg_data
)

summary(fit_cdeg)
confint(fit_cdeg)

pcdeg <- as.data.table(
  expand.grid(
    race.cat = unique(deg_data$race.cat),
    age.grp = 1:5,
    hiv.ego = 0:1,
    degmain_trunc2 = 0:2
  )
)

cdeg_probs <- as.data.table(predict(fit_cdeg, newdata = pcdeg, type = "probs"))

names(cdeg_probs) <- paste0("casldeg", 0:5)

pcdeg <- cbind(pcdeg, cdeg_probs)

melt(pcdeg, measure.vars = names(cdeg_probs)) %>%
  ggplot(aes(x = age.grp, y = value, color = race.cat)) +
  geom_line(aes(linetype = factor(hiv.ego))) +
  geom_point(shape = 21, size = 5, fill = "white") +
  facet_wrap(~ degmain_trunc2 + variable, ncol = 6) +
  ggtitle("Casual degree probabilities") +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  theme_clean()


################################################################################
                     ## CONCURRENCY, CASUAL PARTNERSHIPS ##
################################################################################

# proportion of individuals with concurrent casual partnerships
# casl_concurrent_prob <- deg_data[, mean(casl_conc_ind)]
# round(casl_concurrent_prob, 3)

fit_casl_concurrent <- glm(
  as.numeric(casl_conc_ind) ~ race.cat + factor(age.grp) + hiv.ego + factor(degmain_trunc2),
  family = "binomial",
  data = deg_data
)

summary(fit_casl_concurrent)
coef(fit_casl_concurrent)
confint(fit_casl_concurrent)


################################################################################
                 ## INSTANTANEOUS PARTNERSHIP FORMATION RATE ##
################################################################################

instvars <- c(
  "id", "race.cat", "age", "hiv.ego",
  "pnoa_12m", "m_mp12anum2", "m_mp12instanum2",
  "m_mp12anum2_onepart", "part1once", "pn_ongoing"
)

pinst_rate <- an[, ..instvars]

setkey(pinst_rate, "id")
sapply(pinst_rate, class)

pinst_rate[,
           gtonetime := max(
             pn_ongoing, (m_mp12anum2 - m_mp12instanum2), na.rm = TRUE
           ), id]

pinst_rate[, ":="(
  # calculate number of one-time partnerships as total oral/anal partners
  # minus total anal partners with > 1 contact.
  # NOTE: ARTNet didn't ask about one- or multi- oral sexual contacts.
  pnoa_12m_onetime = fcase(
    pnoa_12m == 0, as.integer(0),
    pnoa_12m == 1 & part1once == 1, pnoa_12m,
    pnoa_12m == 1 & part1once > 1, as.integer(0),
    pnoa_12m > 1 & m_mp12anum2_onepart == 0, pnoa_12m,
    pnoa_12m > 1 & !is.na(gtonetime), pnoa_12m - gtonetime
  ),
  # calculate percentage of oral-or-anal partnerships with whom respondent
  # had anal sex
  pnoa_12m_pct_anal = m_mp12anum2 / pnoa_12m
  )]

pinst_rate[, summary(pnoa_12m_pct_anal)]
pinst_rate[, summary(pnoa_12m_onetime)]
pinst_rate[, .N, .(pnoa_12m_onetime > pnoa_12m)]

ggplot(pinst_rate[pnoa_12m > 0], aes(x = pnoa_12m - pnoa_12m_onetime)) +
  geom_histogram(color = "white", fill = "firebrick", binwidth = 1) +
  theme_clean()

# 34 missing
mice::md.pattern(pinst_rate)

pinst_rate[, inst_wkrate := pnoa_12m_onetime / 52]

pinst_rate <- pinst_rate[!is.na(inst_wkrate)]

pinst_rate[, plot(density(pnoa_12m_onetime))]
pinst_rate[, plot(density(inst_wkrate))]

pinst_rate[, .(mean = mean(pnoa_12m_onetime), var = var(pnoa_12m_onetime))]
pinst_rate[, .(mean = mean(inst_wkrate), var = var(inst_wkrate))]

# join deg_data and pinst_rate
pinst_rate <- deg_data[
  pinst_rate,
  .(id, age, age.grp, race.cat, hiv.ego,
    degmain_trunc2, degcasl,
    pnoa_12m_onetime, inst_wkrate
 ), on = "id"]

print(pinst_rate)


## ... EXPECTED WEEKLY ONE-TIME RATE (BY RACE, AGE, HIV, & MAIN/CASUAL DEGREE)

fit_instrate <- glm.nb(
  pnoa_12m_onetime ~
    race.cat +
    factor(age.grp) +
    factor(degcasl) +
    hiv.ego +
    factor(degmain_trunc2),
  data = pinst_rate
)

summary(fit_instrate)
coef(fit_instrate)
confint(fit_instrate)

pinstrate <- as.data.table(
  expand.grid(
    race.cat = unique(pinst_rate$race.cat),
    age.grp = 1:5,
    degcasl = unique(pinst_rate$degcasl),
    hiv.ego = 0:1,
    degmain_trunc2 = unique(pinst_rate$degmain_trunc2)
  )
)

pinstrate[,
  pred := predict(fit_instrate,
                  newdata = pinstrate,
                  type = "response"
                  )]

pinstrate %>%
  ggplot(aes(x = age.grp, y = pred, color = race.cat)) +
  geom_line(aes(linetype = factor(hiv.ego))) +
  geom_point(shape = 21, size = 4, fill = "white") +
  facet_grid(degmain_trunc2 ~ degcasl) +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  theme_clean()


################################################################################
                             ## ANAL ROLE CLASS ##
################################################################################

mc <- as.data.table(complete(imp_mc, "long", include = TRUE))
otp <- as.data.table(complete(imp_otp, "long", include = TRUE))

mc <- mc[, .(.imp, .id, id, race.i, age.i, role.class)]
otp <- otp[, .(.imp, .id, id, race.i, age.i, role.class)]

all <- rbind(mc, otp)
setkeyv(all, c(".imp", "id"))

all_oneid <- all[, head(.SD, 1), .(.imp, id)]
all_oneid[, .id := seq_len(nrow(all_oneid))]

imp_all <- as.mids(all_oneid)

# Look at relationship betwee age and role.class probability.
# Doesn't look like a spline is needed.
rolebyage <-
  all_oneid[.imp != 0, .N, .(.imp, age.i, role.class)
    ][, P := N/sum(N), keyby = .(.imp, age.i)
    ][, .(P = mean(P)), .(age.i, role.class)]

ggplot(rolebyage, aes(x = age.i, fill = role.class)) +
  geom_col(
    aes(y = P),
    color = "black"
  ) +
  scale_fill_viridis_d(option = "magma")

# Summarize role.class probability across imputed datasets,
# conditional on race/ethnicity and age.
fit_universal_role <-
  with(imp_all, multinom(role.class ~ race.i * age.i))

fit_universal_role_pooled <- summary(pool(fit_universal_role))

roledat <- expand.grid(
  age.i = 18:65,
  race.i = c("black", "hispanic", "other", "white")
)

role_class_probs <- lapply(
  fit_universal_role$analyses, function(x) {
    as.data.table(
      cbind(roledat, predict(x, roledat, type = "probs"))
    )
  }) %>% rbindlist(., idcol = ".imp")

role_class_probs <-
  role_class_probs[, .(Insertive = mean(Insertive),
                       Receptive = mean(Receptive),
                       Versatile = mean(Versatile)), .(age.i, race.i)]

rcp <- melt(role_class_probs, id.vars = c("age.i", "race.i"))

ggplot(rcp, aes(x = age.i, y = value)) +
  geom_line(aes(color = variable), size = 1) +
  facet_wrap(~race.i) +
  scale_color_viridis_d() +
  theme_base()


################################################################################
                             ## WRITE SUMMARIES ##
################################################################################

main_summaries <- list(
  sum_main_total = sum_main_total,
  sum_main_trunc = degmain_t2,
  sum_main_oa    = sum_main_oa,
  sum_main_ao    = sum_main_ao,
  sum_main_oo    = sum_main_oo
)

casl_summaries <- list(
  sum_casl_total = sum_casl_total,
  sum_casl_oa    = sum_casl_oa,
  sum_casl_ao    = sum_casl_ao,
  sum_casl_oo    = sum_casl_oo
)

saveRDS(
  list(
    main_artnet_sum = main_summaries,
    casl_artnet_sum = casl_summaries
  ),
  file = "netstats/aggregate_degree_summaries.Rds"
)

nparams <- list(
  demo = list(
    ai.role.fit = fit_universal_role_pooled,
    ai.role.pr = role_class_probs
  ),
  mc = list(
    racematch = fit_racematch_out,
    age.grpmatch = fit_agematch_out,
    hiv.discord = fit_serodisc_out
  ),
  main = list(
    degprob = fit_mdeg,
    concurrent = fit_main_concurrent,
    durat_wks = fit_main_dur,
    role.class = fit_mrole
  ),
  casl = list(
    degprob = fit_cdeg,
    concurrent = fit_casl_concurrent,
    durat_wks = fit_casl_dur,
    role.class = fit_crole
  ),
  inst = list(
    instrate = fit_instrate
  )
)

print(nparams)

saveRDS(nparams, file = "netstats/predictions.Rds")
