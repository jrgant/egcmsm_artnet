# %% Setup -------------------------------------------------------------------

source("01-import-private-data.R")

pacman::p_load(data.table)
an <- fread(paste0(pd_path, "/artnet-wide-cleaned.csv"))

# %% Individual-Level Partner Data ---------------------------------------------

# Set up ID table - one row per partnership

pid_tab <- melt(
  an[, c("id", paste0("partnn", 1:5))],
  measure = patterns("partnn")
  ) %>%
  .[, pid := stringr::str_remove(variable, "nn")] %>%
  .[order(id, pid), .(id, pid, nn = value)]

setkey(pid_tab, id, pid)

print(pid_tab)


# Convert wide to long

# Select partner variables

pc_select <- names(an)[grep("^part[1-5]{1}", names(an))] %>% .[. != "part2013"]
partner_cols <- c("id", pc_select)

partner_cols

# Partner variable suffixes
slugs <- unique(stringr::str_extract(partner_cols[-1], "(?<=part[1-5]{1}).*"))
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

colsets

# Create long-format tables for each variable type

# @NOTE:
#  - data.table will kick a warning that some of the columns beings collapsed
#    are not of the same type (some double, some integer).
#  - Should be fine. Output column is converted to double.

melted_vars <- lapply(slugs, function(x) {

  selpartcols <- grep(paste0("(?<=part[1-5]{1})", x, "(?!.)"),
                    partner_cols,
                    perl = T)

  currcols <- c("id", partner_cols[selpartcols])

  dat <- an[, ..currcols][order(id)]

  dmelt <- melt(dat, measure = patterns("part"))[order(id)]
  dmelt[, pid := stringr::str_extract(variable, "part[1-5]{1}")]

  names(dmelt)[names(dmelt) == "value"] <- paste0("p_", x)

  out <- dmelt[, c("id", "pid", paste0("p_", x)), with = F][order(pid)]
  setkey(out, id, pid)

  return(out)
  })

length(melted_vars)

checknrow <- sapply(melted_vars, nrow) %>% unlist

table(checknrow == nrow(pid_tab))

# Bind partner data into long-format data.table
plong <- Reduce(
  function(x, y) merge(x, y, by = c("id", "pid"), all = T),
  melted_vars
)

# remove double underscores
names(plong) <- stringr::str_replace_all(names(plong), "__", "_")
head(plong)

plong <- pid_tab[plong, on = .(id, pid)]
head(plong)

# %% Clean/Create variables ----------------------------------------------------

# pull in ego variables
anl <- plong[an[, .(
  id,
  ego.age = age,
  ego.race.cat = race.cat,
  ego.hiv = hiv.ego,
  ego.ongoing = pn_ongoing)], on = .(id)] %>%
  .[order(id, pid)] %>%
  .[nn != ""]  # drop extraneous partner slots

head(anl)
str(anl)

anl[, p_age := ifelse(p_age == 0, NA, p_age)]
head(anl[, .(id, p_age)])
names(anl)

# variables coded with 88/99 as "don't know" or "prefer not to answer"
na_8899 <- paste0("p_", c(
  "hisp", "race", "hiv", "hivknow", "startyyyydk", "startmm", "maincasendmm",
  "ofendmm", "ongoing", "main_ong", "main_end", "iev", "concurr", "concurr2",
  na.omit(stringr::str_extract(names(anl), "(?<=^p_).*_once$")), "stidiag",
  "stitx", "prepuse",
  "prepstart", "prepuse_part", "prepimp_start", "prepimp_start_part",
  "prepimp_condoms", "prepimp_condoms_part", "artuse", "artuse_part",
  "artstart_part", "artimp_start", "artimp_start_part", "artimp_condoms",
  "artimp_condoms_part", "artimp_prep"))

# set 88/99 to NA for all variables in na_8899
anl[, (na_8899) := lapply(.SD, function(x) ifelse(x == 88 | x == 99, NA, x)),
      .SDcols = na_8899]

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
  p_hisp == 0 | is.na(p_hisp), racematch[p_race], "hispanic")]

head(anl[, .(id, pid, p_hisp, p_race, p_race.cat)], n = 50)
anl[, .N, .(id)][, summary(N)]
dfSummary(anl, plain.ascii = T, graph.col = F)

# Recreate partnership type
anl[p_main_ong == 1 | p_main_end == 1, ptype := 1]
anl[p_once == 1, ptype := 3]

anl[(p_main_ong %in% c(0, NA) & p_main_end %in% c(NA, 0)) &
        p_once == 2,
        ptype := 2]


table(is.na(anl$ptype))

head(anl[, .(id, pid, ptype, p_main_ong, p_once, p_ongoing)])

anl[, .N, .(p_once, p_ongoing, p_main_ong, p_main_end, ptype)] %>%
  .[order(ptype, p_once)] %>%
  .[, p := round(N / sum(N), 3)] %>%
  print

anl[, p_ongoing := ifelse(is.na(p_ongoing), 0, p_ongoing)]

anl[, .N, .(p_ongoing)]

maincas_ong <- anl[ptype %in% 1:2 & p_ongoing == 1, .N, .(id, ptype)]
ponetime <- anl[ptype == 3, .N, .(id, ptype)]


pcounts <- dcast(maincas_ong, id ~ ptype, value.var = "N")

names(pcounts)[2:3] <- c("main", "casl_init")
head(pcounts)


deg_data <- pcounts[an[, .(id, pn_ongoing, race.cat, age)]][order(id)]

deg_data[is.na(main), main := 0]
deg_data[is.na(casl_init), casl_init := 0]

deg_data

deg_data[!is.na(pn_ongoing),
         casl_corr :=
          (casl_init > (pn_ongoing - main)) * casl_init +
          (casl_init < (pn_ongoing - main)) * (pn_ongoing - main)]

deg_data[is.na(pn_ongoing), casl_corr := casl_init]

deg_data

mean(deg_data$main)
mean(deg_data$casl_init)
mean(deg_data$casl_corr)

library(summarytools)
deg_data[, ctable(main, casl_init)]
deg_data[, ctable(main, casl_corr)]
hist(deg_data$casl_init)
hist(deg_data$casl_corr)
hist(deg_data$main)

pcols <- c("id", "pid", "ptype", "p_race.cat", "p_hiv",
           names(anl)[grep("^ego", names(anl))])

mp <- anl[ptype == 1 & p_ongoing == 1, ..pcols]
cp <- anl[ptype == 2 & p_ongoing == 1, ..pcols]

head(mp)

mp_racemix <- mp[, ctable(ego.race.cat, p_race.cat, prop = "n")]
cp_racemix <- cp[, ctable(ego.race.cat, p_race.cat, prop = "n")]

mp[, .N, .(ego.race.cat, p_race.cat)] %>%
 .[, P := N / sum(N), .(ego.race.cat)] %>%
 na.omit %>%
 ggplot(aes(x = ego.race.cat, y = P, group = p_race.cat, fill = p_race.cat)) +
  geom_col(position = "dodge", col = "white", width = 0.3) +
  scale_fill_viridis_d(name = "Partner\nRace/ethnicity") +
  labs(x = "Ego Race/ethnicity", y = "Probability") +
  theme_classic()

descr(deg_data$main)
fit_main <- glm(main ~ race.cat + age,
                data = deg_data,
                family = "quasipoisson")

summary(fit_main)

fit_casl_init <- glm(casl_init ~ race.cat + age,
                     data = deg_data,
                     family = "quasipoisson")

summary(fit_casl_init)

fit_casl_corr <- glm(casl_corr ~ race.cat + age,
                     data = deg_data,
                     family = "quasipoisson")

summary(fit_casl_corr)

degpred <- expand.grid(
  race.cat = unique(an$race.cat),
  age = unique(an$age)) %>%
  setDT %>%
  .[order(race.cat, age)]

degpred

main_deg <- cbind(
  degpred,
  pred = predict(fit_main,
                 newdata = degpred,
                 type = "response")
  )

main_deg

casl_init_deg <- cbind(
  degpred,
  pred = predict(fit_casl_init,
                 newdata = degpred,
                 type = "response")
  )

casl_init_deg

casl_corr_deg <- cbind(
  degpred,
  pred = predict(fit_casl_corr,
                 newdata = degpred,
                 type = "response")
  )

casl_corr_deg


broom::tidy(fit_main, conf.int = T, exponentiate = T)

plotdeg_predict <- function(data) {
  ggplot(data, aes(x = age, y = pred, col = race.cat)) +
  geom_line(size = 1) +
  scale_fill_viridis_d() +
  labs(caption = "SOURCE: ART-Net") +
  theme_classic()
}

temp <- "C:/Users/jason/Desktop/"

pm <- plotdeg_predict(main_deg) +
  labs(title = "Main partnerships (quasipoisson)")

pm

# ggsave(paste0(temp, "plotmain.png"),
#        plot = pm,
#        width = 6,
#        height = 4,
#        units = "in")

pci <- plotdeg_predict(casl_init_deg) +
  labs(title = "Casual ongoing, direct (quasipoisson)")

pci

# ggsave(paste0(temp, "plotci.png"),
#        plot = pci,
#        width = 6,
#        height = 4,
#        units = "in")

pcc <- plotdeg_predict(casl_corr_deg) +
  labs(title = "Casual ongoing, indirect (quasipoisson)")

pcc

# ggsave(paste0(temp, "plotcc.png"),
#        plot = pcc,
#        width = 6,
#        height = 4,
#        units = "in")
