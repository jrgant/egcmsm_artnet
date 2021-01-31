# %% Setup --------------------------------------------------------------------

pacman::p_load(data.table,
               ggplot2,
               ggthemes,
               magrittr,
               summarytools)

# %% Import Private Data ------------------------------------------------------

# PD_PATH is a system environment variable on my local machine
artdat <- file.path(
  Sys.getenv("ARTNET_PATH"),
  list.files(Sys.getenv("ARTNET_PATH"), pattern = "ARTNet-Merged")
)

print(artdat)

# av = ARTNet survey data
load(artdat)

# convert variable names to lowercase
names(av) <- tolower(names(av))
names(av) %>% sort %>% print
# convert to data.table
setDT(av)

# Drop unneeded columns
str_drop <- function(string, data = av) {
  cols <- names(data)[grepl(string, names(data))]
  return(cols)
}

drop_cols <- c(
  "surveyeligstarttime",
  "consent",
  "artnetconsent_time",
  "assent",
  "artnetassent_time",
  "surveystarttime",
  "homeless_fr",
  "homeless_p12m",
  str_drop("completion"),
  str_drop("niuse"),
  str_drop("inj"),
  str_drop("suicide"),
  str_drop("k6"),
  str_drop("stigma"),
  str_drop("typ_"),
  str_drop("dchoice"),
  str_drop("mess"),
  str_drop("truvada"),
  str_drop("fp_"),
  # all reported male
  str_drop("gender"),
  "refsites",
  "mj_med",
  "seehcp",
  str_drop("future")
)


hr <- paste0(rep("-", 40), collapse = "")

cat(hr, "COLUMNS DROPPED DURING IMPORT", hr,
  drop_cols, hr,
  paste("TOTAL:", length(drop_cols)), "\n\n",
  sep = "\n"
)

makenull <- function(x) return(NULL)
av[, (drop_cols) := lapply(.SD, makenull), .SDcols = drop_cols]

head(av)

# %% Compare Dataset Attributes -----------------------------------------------

# Function: summarise dataset characteristics
df_metadata <- function(data) {

  sdf <- data.table(
    nrows = nrow(data), unique.ids = length(unique(data$id)),
    numvars = length(data), anydups = anyDuplicated(data)
  )

  return(sdf)
}

# display dataset summary
av_data <- df_metadata(av)
print(av_data)
