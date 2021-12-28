# -- Libraries   
library(tidyverse)
library(lubridate)
library(splines)

## if on server, save with full path
## if not on server, save to home directory
if(Sys.info()["nodename"] == "fermat.dfci.harvard.edu"){
  rda_path <- "/homes10/rafa/dashboard/covidpr/rdas"
} else{
  rda_path <- "rdas"
}


# moving average ----------------------------------------------------------

ma7 <- function(d, y, k = 7) 
  tibble(date = d, moving_avg = as.numeric(stats::filter(y, rep(1/k, k), side = 1)))

sum7 <- function(d, y, k = 7) 
  tibble(date = d, moving_sum = as.numeric(stats::filter(y, rep(1, k), side = 1)))


positivity <- function(dat){
  day_seq <- seq(first_day + weeks(1), max(dat$date), by = "day")
  map_df(day_seq, function(the_day){
    dat %>% filter(date > the_day - weeks(1) & date <= the_day) %>%
      mutate(obs = entry_date <= the_day) %>%
      summarize(date = the_day, 
                people_positives_week = n_distinct(patientId[result == "positive"]),
                people_total_week = n_distinct(patientId),
                fit = people_positives_week / people_total_week,
                lower = qbinom(0.025, people_total_week, fit) / people_total_week,
                upper = qbinom(0.975, people_total_week, fit) / people_total_week,
                obs_people_positives_week = n_distinct(patientId[result == "positive" & obs]),
                obs_people_total_week = n_distinct(patientId[obs]),
                obs_fit = obs_people_positives_week / obs_people_total_week,
                obs_lower = qbinom(0.025, obs_people_total_week, obs_fit) / obs_people_total_week,
                obs_upper = qbinom(0.975, obs_people_total_week, obs_fit) / obs_people_total_week) %>%
      select(date, fit, lower, upper, obs_fit, obs_lower, obs_upper, people_positives_week, people_total_week) 
  })
}

# -- Fixed values

first_day <- make_date(2020, 3, 12)


all_tests_with_id <- readRDS(file.path(rda_path, "all_tests_with_id.rds"))
# Compute time it takes tests to come in ----------------------------------

test_types <- c("Molecular", "Antigens", "Molecular+Antigens")
original_test_types <- c("Molecular", "Antigens")

mol_anti <- all_tests_with_id %>%
  filter(date >= first_day & testType %in% c("Molecular", "Antigens") & 
           result %in% c("positive", "negative")) %>%
  mutate(testType = "Molecular+Antigens") 

all_cases <- all_tests_with_id %>%  
  bind_rows(mol_anti) %>%
  filter(date>=first_day & result == "positive" &
           testType %in% test_types) %>%
  arrange(testType, patientId, date) %>%
  group_by(testType, patientId) %>% ##newId takes reinfection into account
  mutate(newId = paste0(patientId, "-", 
                        c(0, floor(diff(as.numeric(date))/90)))) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  arrange(testType, newId, date) %>%
  group_by(testType, newId) %>%
  slice(1) %>% 
  ungroup() %>%
  select(-patientId, -result) %>%
  arrange(testType, date)

message("Computing lag statistics.")

rezago <- all_tests_with_id  %>% 
  filter(result %in% c("positive", "negative") & 
           testType %in% original_test_types &
           resultCreatedAt >= collectedDate) %>% ## based on @midnucas suggestion: can't be added before it's reported
  group_by(testType) %>%
  mutate(diff = (as.numeric(resultCreatedAt) - as.numeric(collectedDate)) / (60 * 60 * 24),
         Resultado = factor(result, labels = c("Negativos", "Positivos"))) %>%
  ungroup %>%
  select(testType, date, Resultado, diff) %>%
  filter(!is.na(diff))

# By Region ---------------------------------------------------------------

## We now create the tests data table but by region
## The code is repetivie becuase this was added after we had the code for the global case

message("Computing by region statistics.")

tests_by_region <- all_tests_with_id %>%
  bind_rows(mol_anti) %>%
  filter(date >= first_day & 
           testType %in% test_types & 
           result %in% c("positive", "negative")) %>%
  group_by(testType, region, date) %>%
  summarize(people_positives = n_distinct(patientId[result == "positive"]),
            people_total = n_distinct(patientId),
            tests_positives = sum(result == "positive"),
            tests_total = n(),
            .groups = "drop") %>%
  mutate(rate = people_positives / people_total)

## run the function on each test type
fits <- all_tests_with_id %>% 
  bind_rows(mol_anti) %>%
  mutate(entry_date = as_date(orderCreatedAt)) %>%
  filter(date >= first_day & testType %in% test_types & 
           result %in% c("positive", "negative")) %>%
  nest_by(testType, region) %>%
  summarize(positivity(data), .groups = "drop")

## add new variable to test data frame
tests_by_region <- left_join(tests_by_region, fits, by = c("testType", "region", "date")) 

## compute weekly totals for positive tests and total tests
tests_by_region <- tests_by_region %>% 
  group_by(testType, region) %>%
  mutate(tests_positives_week = sum7(d = date, y = tests_positives)$moving_sum) %>%
  mutate(tests_total_week = sum7(d = date, y = tests_total)$moving_sum) %>%
  ungroup()

# compute unique cases ------------------------------------------------------------
cases_by_region <- all_cases %>%
  group_by(testType, region, date) %>% 
  summarize(cases = n(), .groups = "drop") 

# Make sure all dates are included
cases_by_region <- left_join(select(tests_by_region, testType, region, date), 
                             cases_by_region, 
                             by = c("testType", "region", "date")) %>%
  replace_na(list(cases = 0)) 

# compute daily weekly average and add to cases data frame
fits <- cases_by_region %>% 
  group_by(testType, region) %>%
  do(ma7(d = .$date, y = .$cases)) %>%
  rename(cases_week_avg = moving_avg)
cases_by_region <- left_join(cases_by_region, fits, by = c("testType", "region", "date"))

## add new cases and weekly average to tests data frame
tests_by_region <- left_join(tests_by_region, cases_by_region, by = c("testType", "region", "date")) %>%
  mutate(cases_plus_negatives = (people_total_week - people_positives_week + cases_week_avg * 7),
         cases_rate = cases_week_avg * 7 / cases_plus_negatives,
         cases_plus_negatives_daily = people_total - people_positives + cases,
         cases_rate_daily = cases / cases_plus_negatives_daily)

## Not using unique negatives so removed them to save time
## Compute unique negatives

# compute unique negative cases ------------------------------------------------------------
negative_cases_by_region <- all_tests_with_id %>%
  bind_rows(mol_anti) %>%
  filter(date>=first_day & result == "negative" &
           testType %in% test_types) %>%
  group_by(testType, region, patientId) %>%
  arrange(date) %>%
  slice(1) %>%
  ungroup() %>%
  select(-patientId, -result) %>%
  arrange(testType, date) %>%
  group_by(testType, region, date) %>%
  summarize(negative_cases = n(), .groups = "drop")

# Make sure all dates are included
negative_cases_by_region <-  select(tests_by_region, testType, region, date) %>%
  left_join(negative_cases_by_region, by = c("testType", "region", "date")) %>%
  replace_na(list(negative_cases = 0))

# compute daily weekly average and add to negative_cases data frame
fits <- negative_cases_by_region %>%
  group_by(testType, region) %>%
  do(ma7(d = .$date, y = .$negative_cases)) %>%
  rename(negative_cases_week_avg = moving_avg)
negative_cases_by_region <- left_join(negative_cases_by_region, fits, by = c("testType", "region", "date"))

## add new cases and weekly average to tests data frame
tests_by_region <- left_join(tests_by_region, negative_cases_by_region, by = c("testType", "region", "date"))

## add regions populations
pop_by_region <- read_csv("https://raw.githubusercontent.com/rafalab/covidpr/main/data/poblacion-region.csv",
                          skip = 1, col_names = c("rn", "region", "poblacion")) %>% 
  select(region, poblacion) %>%
  mutate(region = factor(region, levels = region[order(poblacion, decreasing = TRUE)]))

tests_by_region$region <- factor(as.character(tests_by_region$region), 
                                 levels = c(levels(pop_by_region$region), "No reportada"))

# By age ------------------------------------------------------------------
message("Computing by age statistics.")

age_starts <- c(0, 10, 15, 20, 30, 40, 65, 75)
age_ends <- c(9, 14, 19, 29, 39, 64, 74, Inf)

## compute daily totals
age_levels <- paste(age_starts, age_ends, sep = " a ")
age_levels[length(age_levels)] <- paste0(age_starts[length(age_levels)],"+")

tests_by_age <- all_tests_with_id %>% 
  bind_rows(mol_anti) %>%
  mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
         age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange, levels = age_levels)) %>%
  filter(date >= first_day & 
           testType %in% test_types & 
           result %in% c("positive", "negative")) %>%
  group_by(testType, date, ageRange) %>% 
  summarize(people_positives = n_distinct(patientId[result == "positive"]),
            people_total = n_distinct(patientId),
            tests_positives = sum(result == "positive"),
            tests_total = n(),
            .groups = "drop") %>%
  mutate(rate = people_positives / people_total)

## run the function on each test type
fits <- all_tests_with_id %>% 
  bind_rows(mol_anti) %>%
  mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
         age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange, levels = age_levels)) %>%
  mutate(entry_date = as_date(orderCreatedAt)) %>%
  filter(date >= first_day & testType %in% test_types & 
           result %in% c("positive", "negative")) %>%
  nest_by(testType, ageRange) %>%
  summarize(positivity(data), .groups = "drop")

## add new variable to test data frame
tests_by_age <- left_join(tests_by_age, fits, by = c("testType", "date", "ageRange")) 

## compute weekly totals for positive tests and total tests
tests_by_age <- tests_by_age %>% 
  group_by(testType, ageRange) %>%
  mutate(tests_positives_week = sum7(d = date, y = tests_positives)$moving_sum) %>%
  mutate(tests_total_week = sum7(d = date, y = tests_total)$moving_sum) %>%
  mutate(tests_week_avg = ma7(d = date, y = tests_total)$moving_avg) %>%
  ungroup()

# compute daily new cases
cases_by_age <- all_cases %>%
  mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
         age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange, levels = age_levels)) %>% 
  group_by(testType, date, ageRange) %>% 
  summarize(cases = n(), .groups = "drop") 

# Make sure all dates are included
cases_by_age <-  left_join(select(tests_by_age, testType, date, ageRange), cases_by_age, by = c("testType", "date", "ageRange")) %>%
  replace_na(list(cases = 0)) 

# compute daily weekly average and add to cases data frame
fits <- cases_by_age %>% 
  group_by(testType, ageRange) %>%
  do(ma7(d = .$date, y = .$cases)) %>%
  rename(cases_week_avg = moving_avg)
cases_by_age <- left_join(cases_by_age, fits, by = c("testType", "date", "ageRange"))

## add new cases and weekly average to tests data frame
tests_by_age <- left_join(tests_by_age, cases_by_age, by = c("testType", "date", "ageRange")) %>%
  mutate(cases_plus_negatives = (people_total_week - people_positives_week + cases_week_avg * 7),
         cases_rate = cases_week_avg * 7 / cases_plus_negatives,
         cases_plus_negatives_daily = people_total - people_positives + cases,
         cases_rate_daily = cases / cases_plus_negatives_daily)

# compute unique negative cases ------------------------------------------------------------
negative_cases_by_age <- all_tests_with_id %>%
  bind_rows(mol_anti) %>%
  mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
         age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange, levels = age_levels)) %>%
  filter(date>=first_day & result == "negative" &
           testType %in% test_types) %>%
  group_by(testType, ageRange, patientId) %>%
  arrange(date) %>%
  slice(1) %>%
  ungroup() %>%
  select(-patientId, -result) %>%
  arrange(testType, date) %>%
  group_by(testType, ageRange, date) %>%
  summarize(negative_cases = n(), .groups = "drop")

# Make sure all dates are included
negative_cases_by_age <-  select(tests_by_age, testType, ageRange, date) %>%
  left_join(negative_cases_by_age, by = c("testType", "ageRange", "date")) %>%
  replace_na(list(negative_cases = 0))

# compute daily weekly average and add to negative_cases data frame
fits <- negative_cases_by_age %>%
  group_by(testType, ageRange) %>%
  do(ma7(d = .$date, y = .$negative_cases)) %>%
  rename(negative_cases_week_avg = moving_avg)

negative_cases_by_age <- left_join(negative_cases_by_age, fits, by = c("testType", "ageRange", "date"))

## add new cases and weekly average to tests data frame
tests_by_age <- left_join(tests_by_age, negative_cases_by_age, by = c("testType", "ageRange", "date"))

## add age populations
pop_by_age <- read_csv("https://raw.githubusercontent.com/rafalab/pr-covid/master/dashboard/data/poblacion-por-edad.csv") %>%
  rename(ageRange = agegroup, poblacion = population) %>%
  mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
         age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange, levels = age_levels)) %>% 
  group_by(ageRange) %>%
  summarize(poblacion = sum(poblacion), .groups = "drop") %>%
  mutate(ageRange = str_replace(ageRange, "-", " to ")) %>%
  mutate(ageRange = factor(ageRange, levels = levels(tests_by_age$ageRange)))

## Rezagos muerte
url <- "https://bioportal.salud.gov.pr/api/administration/reports/deaths/summary"

deaths <- jsonlite::fromJSON(url) %>%
  mutate(date = as_date(ymd_hms(deathDate, tz = "America/Puerto_Rico"))) %>%
  mutate(date = if_else(date < first_day | date > today(), 
                        as_date(ymd_hms(reportDate, tz = "America/Puerto_Rico")),
                        date)) %>%
  mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
         age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange, levels = age_levels)) 

rezago_mort <- deaths %>% 
  filter(!is.na(date)) %>%
  mutate(bulletin_date = as_date(ymd_hms(reportDate, tz = "America/Puerto_Rico"))) %>%
  arrange(date, bulletin_date) %>%
  mutate(diff = (as.numeric(bulletin_date) - as.numeric(date))) %>%
  select(date, diff)


# -- Save data

message("Saving data.")

save(rezago, file = file.path(rda_path, "rezago.rda"))

save(rezago_mort, file = file.path(rda_path, "rezago_mort.rda"))

save(tests_by_region, pop_by_region, file = file.path(rda_path, "regions.rda"))

save(deaths_by_age, tests_by_age, pop_by_age, file = file.path(rda_path, "by-age.rda"))
