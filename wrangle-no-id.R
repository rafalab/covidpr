## This is so we have municipio data
# -- Libraries   
library(tidyverse)
library(lubridate)
library(splines)
if(grepl("fermat|leo", Sys.info()["nodename"])){
  rda_path <- "/homes10/rafa/dashboard/covidpr/rdas"
} else{
  rda_path <- "rdas"
}

# moving average ----------------------------------------------------------

ma7 <- function(d, y, k = 7) 
  tibble(date = d, moving_avg = as.numeric(stats::filter(y, rep(1/k, k), side = 1)))

sum7 <- function(d, y, k = 7) 
  tibble(date = d, moving_sum = as.numeric(stats::filter(y, rep(1, k), side = 1)))

# -- Fixed values
first_day <- make_date(2020, 3, 12)

last_day <- now(tz= "America/Puerto_Rico") + days(1)

last_complete_day <- today() - 1

the_years <- seq(2020, year(today()))

age_levels <-  c("0 to 9", "10 to 19", "20 to 29", "30 to 39", "40 to 49", "50 to 59", "60 to 69", 
                 "70 to 79", "80 to 89", "90 to 99", "100 to 109", "110 to 119", "120 to 129")

imputation_delay  <- 2

alpha <- 0.05

## filter by date example: ?createdAtStartDate=2021-09-09T04:00:00Z&createdAtEndDate=2021-09-10T04:00:00Z
#test_url <- "https://bioportal.salud.pr.gov/api/administration/reports/minimal-info-unique-tests"
test_url <- "https://bioportal-apim.salud.pr.gov/bioportal/administration/reports/minimal-info-unique-tests"

the_days <- unique(c(first_day, make_date(the_years[-1], 1, 1), last_day))
the_days <- the_days |> with_tz(tzone = "GMT") |> format("%Y-%m-%dT%H:%M:%SZ")


get_bioportal <- function(url){
  setDT(jsonlite::fromJSON(
    rawToChar(
      httr::GET(url, httr::content_type('application/json'),
                httr::add_headers('Accept-Enconding'="br"))$content)
  ))
}

original_test_types <- c("Molecular", "Antigens")

tmp <- merge(data.frame(test_type = original_test_types), 
             data.frame(start = head(the_days, -1), end =  tail(the_days,-1)))

## Define the queries
queries <- apply(tmp, 1, function(x)
  paste0(test_url,
         "?testType=", x[1],
         "&createdAtStartDate=", x[2], "&createdAtEndDate=", x[3]))


# Reading and wrangling cases data from database ---------------------------
message("Reading case data.")

all_tests <- lapply(queries, get_bioportal)
all_tests <- do.call(rbind, all_tests) |> unique() ##unique removes duplicate rows

message("Processing case data.")

all_tests <- all_tests %>%  
  rename(patientCity = city) %>%
  as_tibble() %>%
  mutate(testType = str_to_title(testType),
         testType = ifelse(testType == "Antigeno", "Antigens", testType),
         collectedDate  = mdy(collectedDate),
         reportedDate   = mdy(reportedDate),
         createdAt      = mdy_hm(createdAt),
         ageRange       = na_if(ageRange, "N/A"),
         ageRange       = factor(ageRange, levels = age_levels),
         patientCity    = ifelse(patientCity == "Loiza", "Loíza", patientCity),
         patientCity    = ifelse(patientCity == "Rio Grande", "Río Grande", patientCity),
         patientCity    = factor(patientCity),
         result         = tolower(result),
         result         = case_when( 
           (grepl("covid", result) | grepl("sars-cov-2", result)) &
             grepl("positive", result) ~ "positive",
           grepl("influenza", result) ~ "other",
           grepl("positive", result) ~ "positive",
           grepl("negative", result) ~ "negative",
           result == "not detected" ~ "negative",
           TRUE ~ "other")) %>%
  arrange(reportedDate, collectedDate) %>%
  filter(testType %in% test_types)

## fixing bad dates: if you want to remove bad dates instead, change FALSE TO TRUE
if(FALSE){
  ## remove bad dates
  all_tests <- all_tests %>% 
    filter(!is.na(collectedDate) & year(collectedDate) %in% the_years & collectedDate <= today()) %>%
    mutate(date = collectedDate) 
} else{
  ## Impute missing dates and remove inconsistent dates
  all_tests <- all_tests %>% 
    mutate(date = if_else(collectedDate > reportedDate, reportedDate, collectedDate)) %>% ## if collectedDate is in the future make it reportedDate
    mutate(date = if_else(is.na(collectedDate), reportedDate - days(imputation_delay),  collectedDate)) %>%
    mutate(date = if_else(!year(date) %in% the_years, reportedDate - days(imputation_delay),  date)) %>%
    filter(year(date) %in% the_years & date <= today()) %>%
    arrange(date, reportedDate)
}

# -- summaries stratified by age group and patientID
mol_anti_2 <-  all_tests %>%
  filter(date >= first_day & testType %in% c("Molecular", "Antigens") & 
           result %in% c("positive", "negative")) %>%
  mutate(testType = "Molecular+Antigens") 

tests_by_strata <- all_tests %>%  
  bind_rows(mol_anti_2) %>%
  filter(date >= first_day & testType %in% test_types & 
           result %in% c("positive", "negative")) %>%
  filter(date>=first_day) %>%
  mutate(patientCity = fct_explicit_na(patientCity, "No reportado")) %>%
  mutate(ageRange = fct_explicit_na(ageRange, "No reportado")) %>%
  group_by(testType, date, patientCity, ageRange, .drop = FALSE) %>%
  summarize(positives = sum(result == "positive"), tests = n(), .groups="drop") %>%
  ungroup()

save(tests_by_strata, file = file.path(rda_path, "tests_by_strata.rda"))

## For backward compatibility keep only molecular
all_tests <- all_tests %>%  filter(testType == "Molecular")
saveRDS(all_tests, file = file.path(rda_path, "all_tests.rds"), compress = "xz")

