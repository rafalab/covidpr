library(tidyverse)
library(lubridate)
library(splines)

if(grepl("fermat|leo", Sys.info()["nodename"])){
  rda_path <- "/homes10/rafa/dashboard/covidpr/rdas"
} else{
  rda_path <- "rdas"
}

ma7 <- function(d, y, k = 7) 
  tibble(date = d, moving_avg = as.numeric(stats::filter(y, rep(1/k, k), side = 1)))

sum7 <- function(d, y, k = 7) 
  tibble(date = d, moving_sum = as.numeric(stats::filter(y, rep(1, k), side = 1)))

# -- Fixed values
pr_pop <- 3285874 ## population of puerto rico

icu_beds <- 229 #if available beds is missing change to this

first_day <- make_date(2020, 3, 12)

last_complete_day <- today() - 1

the_years <- seq(2020, year(today()))

age_levels <-  c("0 to 9", "10 to 19", "20 to 29", "30 to 39", "40 to 49", "50 to 59", "60 to 69", 
                 "70 to 79", "80 to 89", "90 to 99", "100 to 109", "110 to 119", "120 to 129")

imputation_delay  <- 2

alpha <- 0.05

# Computing positivity rate by lab ----------------------------------------

url <- "https://bioportal.salud.gov.pr/api/administration/reports/tests-by-collected-date-and-entity"

message("Reading lab data.")

all_labs_data <- jsonlite::fromJSON(url)

labs <- all_labs_data %>%
  select(-molecular, -serological, -antigens) %>%
  rename(Laboratorio = entityName,
         date = collectedDate) %>%
  mutate(date = as_date(date),
         Laboratorio = str_trim(str_remove_all(tolower(Laboratorio), "\t|inc|\\.|\\,")))


##check the most common labs
if(FALSE){
  freqs <- bind_cols(labs, all_labs_data$molecular) %>% 
    filter(date > make_date(2021, 1, 1)) %>%
    group_by(Laboratorio) %>%
    summarize(freq = sum(total), .groups = "drop") %>% 
    ungroup()
  freqs %>% View()
}

message("Processing lab data.")

labs <- labs %>%
  mutate(Laboratorio = case_when(str_detect(Laboratorio, "toledo") ~ "Toledo",
                                 str_detect(Laboratorio, "bcel") ~ "BCEL",
                                 str_detect(Laboratorio, "landr") ~ "Landrón",
                                 str_detect(Laboratorio, "cmt") ~ "CMT",
                                 str_detect(Laboratorio, "villa ana") ~ "Villa Ana",
                                 str_detect(Laboratorio, "labcorp") ~ "LabCorp",
                                 str_detect(Laboratorio, "quest") ~ "Quest",
                                 str_detect(Laboratorio, "borinquen") ~ "Borinquen",
                                 str_detect(Laboratorio, "coreplus") ~ "CorePlus",
                                 str_detect(Laboratorio, "martin\\s") ~ "Marin",
                                 str_detect(Laboratorio, "noy") ~ "Noy",
                                 str_detect(Laboratorio, "hato rey pathology|hrp") ~ "HRP",
                                 str_detect(Laboratorio, "inno") ~ "Inno Diagnostics",
                                 str_detect(Laboratorio, "immuno reference lab") ~ "Immuno Reference",
                                 str_detect(Laboratorio, "forense") ~ "Ciencias Forense",
                                 #str_detect(Laboratorio, "nichols") ~ "Quest USA",
                                 #str_detect(Laboratorio, "southern pathology services") ~ "Southern Pathology",
                                 TRUE ~ "Otros"))

molecular <- all_labs_data$molecular %>% 
  mutate(testType = "Molecular",
         positives = positives + presumptivePositives,
         negatives = negatives,
         tests = positives + negatives) %>%
  select(testType, positives, tests)
molecular <- bind_cols(labs, molecular) 

serological <-  all_labs_data$serological %>%
  mutate(testType = "Serological",
         positives = positives,
         negatives = negatives,
         tests = positives + negatives) %>%
  select(testType, positives, tests)
serological <- bind_cols(labs, serological) 

antigens <-  all_labs_data$antigens %>%
  mutate(testType = "Antigens",
         positives = positives,
         negatives = negatives,
         tests = positives + negatives) %>%
  select(testType, positives, tests)
antigens <- bind_cols(labs, antigens) 

labs <- bind_rows(molecular, serological, antigens) %>%
  filter(date >= first_day & date <= today()) %>%
  group_by(testType, date, Laboratorio) %>%
  summarize(positives = sum(positives),
            tests = sum(tests),
            missing_city = sum(totalMissingCity),
            missing_phone = sum(totalMissingPhoneNumber),
            .groups = "drop")


lab_positivity <- function(dat){
  day_seq <- seq(first_day + weeks(1), max(labs$date), by = "day")
  map_df(day_seq, function(the_day){
    ret <- dat %>% 
      filter(date > the_day - weeks(1) & date <= the_day) %>%
      summarize(date = the_day, 
                n = sum(tests),
                tests_week_avg  = n / 7, 
                fit = ifelse(n==0, 0, sum(positives) / n),
                lower = qbinom(0.025, n, fit) / n,
                upper = qbinom(0.975, n, fit) / n) %>%
      select(date, fit, lower, upper, tests_week_avg)
  })
}

fits <- labs %>% 
  nest_by(testType, Laboratorio) %>%
  summarize(lab_positivity(data), .groups = "drop") %>%
  group_by(testType, date) %>%
  mutate(prop = tests_week_avg / sum(tests_week_avg)) 

labs <- left_join(fits, labs, by = c("testType", "date", "Laboratorio"))


## For Eddie
lab_tab <- all_labs_data %>%
  select(-molecular, -serological, -antigens) %>%
  rename(Laboratorio = entityName,
         date = collectedDate) %>%
  mutate(date = as_date(date),
         Laboratorio = str_trim(str_remove_all(tolower(Laboratorio), "\t|inc|\\.|\\,")))

## wrange some of the names
lab_tab <- lab_tab %>%
  mutate(Laboratorio = case_when(str_detect(Laboratorio, "toledo") ~ "Toledo",
                                 str_detect(Laboratorio, "bcel") ~ "BCEL",
                                 str_detect(Laboratorio, "landr") ~ "Landrón",
                                 str_detect(Laboratorio, "cmt") ~ "CMT",
                                 str_detect(Laboratorio, "villa ana") ~ "Villa Ana",
                                 str_detect(Laboratorio, "labcorp") ~ "LabCorp",
                                 str_detect(Laboratorio, "quest") ~ "Quest",
                                 str_detect(Laboratorio, "borinquen") ~ "Borinquen",
                                 str_detect(Laboratorio, "coreplus") ~ "CorePlus",
                                 str_detect(Laboratorio, "martin\\s") ~ "Marin",
                                 str_detect(Laboratorio, "noy") ~ "Noy",
                                 str_detect(Laboratorio, "hato rey pathology|hrp") ~ "HRP",
                                 str_detect(Laboratorio, "inno") ~ "Inno Diagnostics",
                                 str_detect(Laboratorio, "immuno reference lab") ~ "Immuno Reference",
                                 str_detect(Laboratorio, "forense") ~ "Ciencias Forense",
                                 str_detect(Laboratorio, "nichols") ~ "Quest USA",
                                 str_detect(Laboratorio, "southern pathology services") ~ "Southern Pathology",
                                 TRUE ~ str_replace(str_to_title(Laboratorio), "Ii", "II")))

molecular <- all_labs_data$molecular %>% 
  mutate(testType = "Molecular",
         tests = positives + presumptivePositives + negatives) %>%
  select(testType, tests)
molecular <- bind_cols(lab_tab, molecular) 

serological <-  all_labs_data$serological %>%
  mutate(testType = "Serological",
         tests = positives + negatives) %>%
  select(testType, tests)
serological <- bind_cols(lab_tab, serological) 

antigens <-  all_labs_data$antigens %>%
  mutate(testType = "Antigens",
         tests = positives + negatives) %>%
  select(testType, tests)
antigens <- bind_cols(lab_tab, antigens) 

lab_tab <- bind_rows(molecular, serological, antigens) %>%
  filter(date >= first_day & date <= today()) %>%
  group_by(testType, date, Laboratorio) %>%
  summarize(tests = sum(tests),.groups = "drop")

lab_tab  <- lab_tab %>% group_by(Laboratorio, testType) %>% 
  mutate(total = sum(tests), .groups = "drop") %>% 
  ungroup() %>%
  group_by(testType) %>%
  mutate(Laboratorio = ifelse(total < 100, "Otros", Laboratorio)) %>%
  group_by(testType, date, Laboratorio) %>%
  summarize(tests = sum(tests),.groups = "drop") 

save(labs, file = file.path(rda_path, "labs.rda"))

save(lab_tab, file = file.path(rda_path, "lab_tab.rda"))
