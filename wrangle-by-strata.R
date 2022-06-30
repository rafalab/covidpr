# -- Libraries   
library(tidyverse)
library(lubridate)
library(data.table)

## if on server, save with full path
## if not on server, save to home directory
if(grepl("fermat|leo", Sys.info()["nodename"])){
  rda_path <- "/homes10/rafa/dashboard/covidpr/rdas"
} else{
  rda_path <- "rdas"
}

# moving average ----------------------------------------------------------

ma7 <- function(y, k = 7) as.numeric(stats::filter(y, rep(1/k, k), side = 1))

sum7 <- function(y, k = 7) as.numeric(stats::filter(y, rep(1, k), side = 1))


first_day <- make_date(2020, 3, 12)

last_complete_day <- today() - 1

age_starts <- c(0, 10, 15, 20, 30, 40, 65, 75)
age_ends <- c(9, 14, 19, 29, 39, 64, 74, Inf)

age_levels <- paste(age_starts, age_ends, sep = " a ")
age_levels[length(age_levels)] <- paste0(age_starts[length(age_levels)],"+")

##load data
all_tests_with_id <- readRDS(file.path(rda_path, "all_tests_with_id.rds"))

all_tests_with_id [, `:=`(
  age_start = as.numeric(str_extract(ageRange, "^\\d+")),
  age_end = as.numeric(str_extract(ageRange, "\\d+$")))]
all_tests_with_id[, ageRange := factor(age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))], levels = age_levels)]


test_types <- c("Molecular", "Antigens", "Molecular+Antigens")
original_test_types <- c("Molecular", "Antigens")


# -- Computing observed positivity rate
## adding a new test type that combines molecular and antigens
mol_anti <- all_tests_with_id[date >= first_day & testType %in% c("Molecular", "Antigens") & result %in% c("positive", "negative")]
mol_anti[, testType := "Molecular+Antigens"]


cases <- rbind(all_tests_with_id, mol_anti)[
  date >= first_day & result == "positive" &
    testType %in% test_types][order(testType, patientId, date)]

cases[, n:= .N,  keyby = c("testType", "patientId")]
cases[, days:=0]
cases[n>1, days := c(0, floor(diff(as.numeric(date))/90)), by = c("testType", "patientId")]
cases[, newId := paste(paste0(patientId, "-", days))]
cases <- cases[order(testType, newId, date)]
cases <- cases[, .SD[1], keyby = c("testType", "newId")] 
cases[, patientId := NULL]
cases[, result:=NULL]

cases <- cases[order(testType, date)]
all_cases <- copy(cases)
rm(cases)

message("Computing lag statistics.")

rezago <- all_tests_with_id[result %in% c("positive", "negative") & 
                              testType %in% original_test_types &
                              resultCreatedAt >= collectedDate]
## based on @midnucas suggestion: can't be added before it's reported
rezago[, `:=`(diff = (as.numeric(resultCreatedAt) - as.numeric(collectedDate)) / (60 * 60 * 24),
              Resultado = factor(result, labels = c("Negativos", "Positivos")))]
rezago <- rezago[!is.na(diff), .(testType, date, Resultado, diff)]


## compute daily totals bt region

all_dates <- CJ(testType = test_types, 
                region = factor(levels(all_tests_with_id$region)),
                date =  seq(first_day, max(all_tests_with_id$date), by = "day"))

tests_by_region <- rbind(all_tests_with_id, mol_anti)[date >= first_day & 
           testType %in% test_types & 
           result %in% c("positive", "negative")]
tests_by_region <- tests_by_region[, .(people_positives = uniqueN(patientId[result == "positive"]),
                   people_total = uniqueN(patientId),
                   tests_positives = sum(result == "positive"),
                   tests_total = .N), keyby = c("testType", "region", "date")]
tests_by_region <- merge(tests_by_region, all_dates, by = c("testType", "region", "date"), all.y = TRUE)
tests_by_region[is.na(tests_by_region)] <- 0
tests_by_region[, rate :=  people_positives / people_total]

## define function to compute weekly distinct cases
## and use this to compute percent of people with positive tests

positivity <- function(dat){
  day_seq <- seq(first_day + weeks(1), max(dat$date), by = "day")
  res <- lapply(day_seq, function(the_day){
    tmp <- dat[date > the_day - weeks(1) & date <= the_day]
    tmp[, obs := entry_date <= the_day]
    return(c(
      uniqueN(tmp$patientId[tmp$pos]),
      uniqueN(tmp$patientId),
      uniqueN(tmp$patientId[tmp$pos & tmp$obs]),
      uniqueN(tmp$patientId[tmp$obs])))
  })
  
  res <- do.call(rbind, res)
  colnames(res) <- c("people_positives_week", "people_total_week",
                     "obs_people_positives_week", "obs_people_total_week")
  res <- as.data.table(res)
  res$date <- day_seq
  res[, fit := people_positives_week / people_total_week]
  res[, obs_fit := obs_people_positives_week / obs_people_total_week]
  res[, lower := qbinom(0.025, people_total_week, fit) / people_total_week]
  res[, upper := qbinom(0.975, people_total_week, fit) / people_total_week]
  res[, obs_lower := qbinom(0.025, obs_people_total_week, obs_fit) / obs_people_total_week]
  res[, obs_upper := qbinom(0.975, obs_people_total_week, obs_fit) / obs_people_total_week]
  return(res[, .(date, fit, lower, upper, obs_fit, obs_lower, obs_upper, people_positives_week, people_total_week)])
}

message("Computing positivity.")

## run the function on each test type
fits <- rbind(all_tests_with_id, mol_anti) 
fits <- fits[date >= first_day & testType %in% test_types & result %in% c("positive", "negative")]
fits[,`:=`(pos = result == "positive", entry_date = as_date(resultCreatedAt))]
fits <- fits[, positivity(.SD), by = c("testType", "region")]  

## add new variable to test data frame
tests_by_region <- merge(tests_by_region, fits, by = c("testType", "region", "date"), all.x = TRUE) 
cols <- c("people_positives_week", "people_total_week")
tests_by_region[, (cols) := lapply(.SD, replace_na, 0), .SDcols = cols]

## compute weekly totals for positive tests and total tests
tests_by_region[, `:=`(tests_positives_week = sum7(tests_positives),
             tests_total_week = sum7(tests_total)), by = c("testType", "region")]

# compute unique cases ------------------------------------------------------------

cases_by_region <- all_cases[, .(cases = .N), keyby = .(testType, region, date)]
  
# Make sure all dates are included
cases_by_region <-  merge(tests_by_region[, .(testType, region, date)], 
                          cases_by_region, 
                by = c("testType", "region", "date"), all.x = TRUE) 
cases_by_region[, cases := replace_na(cases, 0)]
  
# compute daily weekly average and add to cases data frame
cases_by_region[, cases_week_avg := ma7(cases), by = c("testType", "region")]

## add new cases and weekly average to tests data frame
tests_by_region <- merge(tests_by_region, cases_by_region, by = c("testType", "region", "date"), all.x = TRUE)
tests_by_region[, cases_plus_negatives := (people_total_week - people_positives_week + cases_week_avg * 7)]
tests_by_region[, `:=`(cases_rate = cases_week_avg * 7 / cases_plus_negatives,
             cases_plus_negatives_daily = people_total - people_positives + cases)]
tests_by_region[, cases_rate_daily := cases / cases_plus_negatives_daily]
         
## Compute unique negatives

# compute unique negative cases ------------------------------------------------------------
message("Computing unique negatives.")

negative_cases_by_region <- rbind(all_tests_with_id, mol_anti)[
  date>=first_day & result == "negative" &
           testType %in% test_types]
negative_cases_by_region <- negative_cases_by_region[order(date)]
negative_cases_by_region <- negative_cases_by_region[, .SD[1], keyby = .(testType, patientId)]
negative_cases_by_region[, patientId:=NULL]
negative_cases_by_region[, result:=NULL]
negative_cases_by_region <- negative_cases_by_region[, .(negative_cases = .N), keyby = .(testType, region, date)]
negative_cases_by_region <- negative_cases_by_region[order(testType, region, date)]
  
# Compute daily weekly average and add to negative_cases data frame
# but first make sure all dates are included
negative_cases_by_region <-  merge(tests_by_region[, .(testType, region, date)], negative_cases_by_region, by = c("testType", "region", "date"), all.x = TRUE)
negative_cases_by_region[, negative_cases := replace_na(negative_cases, 0)]
  
# compute daily weekly average and add to negative_cases data frame
negative_cases_by_region[, negative_cases_week_avg := ma7(negative_cases), by = c("testType", "region")]

## add new cases and weekly average to tests data frame
tests_by_region <- merge(tests_by_region, negative_cases_by_region, by = c("testType", "region", "date"), all.x=TRUE)

## add regions populations
pop_by_region <- read_csv("https://raw.githubusercontent.com/rafalab/covidpr/main/dashboard/data/poblacion-region.csv",
                          skip = 1, col_names = c("rn", "region", "poblacion")) %>% 
  select(region, poblacion) %>%
  mutate(region = factor(region, levels = region[order(poblacion, decreasing = TRUE)]))

tests_by_region$region <- factor(as.character(tests_by_region$region), 
                                 levels = c(levels(pop_by_region$region), "No reportada"))

## By Age

message("Computing by age statistics.")

all_dates <- CJ(testType = test_types, 
                ageRange = factor(age_levels),
                date =  seq(first_day, max(all_tests_with_id$date), by = "day"))


## compute daily totals
tests_by_age <- rbind(all_tests_with_id, mol_anti)[date >= first_day & 
                                                        testType %in% test_types & 
                                                        result %in% c("positive", "negative")]

tests_by_age <- tests_by_age[, .(people_positives = uniqueN(patientId[result == "positive"]),
                                       people_total = uniqueN(patientId),
                                       tests_positives = sum(result == "positive"),
                                       tests_total = .N), keyby = c("testType", "ageRange", "date")]
tests_by_age <- merge(tests_by_age, all_dates, by = c("testType", "ageRange", "date"), all.y = TRUE)
tests_by_age[is.na(tests_by_age)] <- 0
tests_by_age[, rate :=  people_positives / people_total]


message("Computing positivity.")

## run the function on each test type
fits <- rbind(all_tests_with_id, mol_anti) 
fits <- fits[date >= first_day & testType %in% test_types & result %in% c("positive", "negative")]
fits[,`:=`(pos = result == "positive", entry_date = as_date(resultCreatedAt))]
fits <- fits[, positivity(.SD), by = c("testType", "ageRange")]  

## add new variable to test data frame
tests_by_age <- merge(tests_by_age, fits, by = c("testType", "ageRange", "date"), all.x = TRUE) 
cols <- c("people_positives_week", "people_total_week")
tests_by_age[, (cols) := lapply(.SD, replace_na, 0), .SDcols = cols]

## compute weekly totals for positive tests and total tests
tests_by_age[, `:=`(tests_positives_week = sum7(tests_positives),
                       tests_total_week = sum7(tests_total)), by = c("testType", "ageRange")]

# compute unique cases ------------------------------------------------------------

cases_by_age <- all_cases[, .(cases = .N), keyby = .(testType, ageRange, date)]

# Make sure all dates are included
cases_by_age <-  merge(tests_by_age[, .(testType, ageRange, date)], 
                          cases_by_age, 
                          by = c("testType", "ageRange", "date"), all.x = TRUE) 
cases_by_age[, cases := replace_na(cases, 0)]

# compute daily weekly average and add to cases data frame
cases_by_age[, cases_week_avg := ma7(cases), by = c("testType", "ageRange")]

## add new cases and weekly average to tests data frame
tests_by_age <- merge(tests_by_age, cases_by_age, by = c("testType", "ageRange", "date"), all.x = TRUE)
tests_by_age[, cases_plus_negatives := (people_total_week - people_positives_week + cases_week_avg * 7)]
tests_by_age[, `:=`(cases_rate = cases_week_avg * 7 / cases_plus_negatives,
                       cases_plus_negatives_daily = people_total - people_positives + cases)]
tests_by_age[, cases_rate_daily := cases / cases_plus_negatives_daily]

## Compute unique negatives

# compute unique negative cases ------------------------------------------------------------
message("Computing unique negatives.")

negative_cases_by_age <- rbind(all_tests_with_id, mol_anti)[
  date>=first_day & result == "negative" &
    testType %in% test_types]

negative_cases_by_age <- negative_cases_by_age[order(date)]
negative_cases_by_age <- negative_cases_by_age[, .SD[1], keyby = .(testType, patientId)]
negative_cases_by_age[, patientId:=NULL]
negative_cases_by_age[, result:=NULL]
negative_cases_by_age <- negative_cases_by_age[, .(negative_cases = .N), keyby = .(testType, ageRange, date)]
negative_cases_by_age <- negative_cases_by_age[order(testType, ageRange, date)]

# Compute daily weekly average and add to negative_cases data frame
# but first make sure all dates are included
negative_cases_by_age <-  merge(tests_by_age[, .(testType, ageRange, date)], negative_cases_by_age, by = c("testType", "ageRange", "date"), all.x = TRUE)
negative_cases_by_age[, negative_cases := replace_na(negative_cases, 0)]

# compute daily weekly average and add to negative_cases data frame
negative_cases_by_age[, negative_cases_week_avg := ma7(negative_cases), by = c("testType", "ageRange")]

## add new cases and weekly average to tests data frame
tests_by_age <- merge(tests_by_age, negative_cases_by_age, by = c("testType", "ageRange", "date"), all.x=TRUE)

#

pop_by_age <- read_csv("https://raw.githubusercontent.com/rafalab/covidpr/main/dashboard/data/poblacion-por-edad.csv") %>%
  rename(ageRange = agegroup, poblacion = population) %>%
  mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
         age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange, levels = age_levels)) %>% 
  group_by(ageRange) %>%
  summarize(poblacion = sum(poblacion), .groups = "drop") %>%
  mutate(ageRange = str_replace(ageRange, "-", " to ")) %>%
  mutate(ageRange = factor(ageRange, levels = levels(tests_by_age$ageRange)))

message("Computing deaths by age.")

load(file.path(rda_path, "deaths.rda"))

## Use this to assure all dates are included
all_dates <- data.frame(date = seq(first_day, max(c(last_complete_day, deaths$date)), by = "day"))
deaths_by_age <- deaths %>%
  filter(!is.na(ageRange)) %>%
  group_by(date, ageRange) %>%
  summarize(deaths = n(), .groups = "drop") %>%
  ungroup() %>%
  full_join(all_dates, by = "date") %>%
  complete(date, nesting(ageRange), fill = list(deaths = 0)) %>%
  filter(!is.na(ageRange)) %>%
  group_by(ageRange) %>%
  mutate(deaths_week_avg = ma7(deaths))

# -- Save data

message("Saving data.")

save(rezago, file = file.path(rda_path, "rezago.rda"))

save(tests_by_region, pop_by_region, file = file.path(rda_path, "regions.rda"))

save(deaths_by_age, tests_by_age, pop_by_age, file = file.path(rda_path, "by-age.rda"))


