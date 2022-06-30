## Run this once
## then wrangle automatically updates it

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

# -- Fixed values
pr_pop <- 3285874 ## population of puerto rico

icu_beds <- 229 #if available beds is missing change to this

first_day <- make_date(2020, 3, 12)

last_day <- make_date(2022, 6, 15)

the_years <- seq(2020, year(today()))

age_levels <-  paste(seq(0, 125, 5), "to", seq(4, 129, 5))

imputation_delay  <- 2

alpha <- 0.05

## Load latest data
message("Downloading dataset.")

## filter by date example: ?createdAtStartDate=2021-09-09T04:00:00Z&createdAtEndDate=2021-09-10T04:00:00Z
cases_url <- "https://bioportal.salud.pr.gov/api/administration/reports/orders/basic"

last_day<- last_day |>
  with_tz(tzone = "GMT") |>
  str_replace(" ", "T") |>
  paste0("Z")


cases_url_molecular <-  paste0(cases_url, 
                               "?testType=Molecular",
                               "&createdAtEndDate=",
                               last_day)

cases_url_antigens <- paste0(cases_url, 
                             "?testType=Antigens",
                             "&createdAtEndDate=",
                             last_day)

get_bioportal <- function(url){
  setDT(jsonlite::fromJSON(
    rawToChar(
      httr::GET(url, httr::content_type('application/json'),
                httr::add_headers('Accept-Enconding'="br"))$content)
  ))
}

original_test_types <- c("Molecular", "Antigens")

# Reading and wrangling cases data from database ---------------------------
message("Reading case data.")

all_tests_with_id_molecular <- get_bioportal(cases_url_molecular)
all_tests_with_id_antigens <- get_bioportal(cases_url_antigens)
all_tests_with_id <- rbind(all_tests_with_id_molecular, all_tests_with_id_antigens)
rm(all_tests_with_id_molecular, all_tests_with_id_antigens); gc(); gc()

all_tests_with_id <- setDT(all_tests_with_id)
all_tests_with_id <- distinct(all_tests_with_id)

message("Processing case data.")
all_tests_with_id[, result := factor(tolower(result))]
all_tests_with_id[, `:=`(testType = fifelse(str_to_title(testType) == "Antigeno", "Antigens", testType),
                         collectedDate = ymd_hms(collectedDate, tz = "America/Puerto_Rico"),
                         reportedDate = ymd_hms(reportedDate, tz = "America/Puerto_Rico"),
                         orderCreatedAt = ymd_hms(orderCreatedAt, tz = "America/Puerto_Rico"),
                         resultCreatedAt = ymd_hms(resultCreatedAt, tz = "America/Puerto_Rico"),
                         ageRange       = factor(na_if(ageRange, "N/A"), levels = age_levels),
                         region = recode(region,`N/A` = "No reportada",
                                         Bayamon = "Bayamón", 
                                         Mayaguez = "Mayagüez", .missing = "No reportada"),
                         result = fcase( 
                           (grepl("covid", result) | grepl("sars-cov-2", result)) &
                             grepl("positive", result),  "positive",
                           grepl("influenza", result),  "other",
                           grepl("positive", result), "positive",
                           grepl("negative", result), "negative",
                           result == "not detected", "negative", default = "other"))]

all_tests_with_id <- all_tests_with_id[testType %in% original_test_types]
all_tests_with_id <- all_tests_with_id[order(reportedDate, collectedDate, patientId)] 


## adding  levels to regions
pop_by_region <- read_csv("https://raw.githubusercontent.com/rafalab/covidpr/main/dashboard/data/poblacion-region.csv",
                          skip = 1, col_names = c("rn", "region", "poblacion")) %>% 
  select(region, poblacion) %>%
  mutate(region = factor(region, levels = region[order(poblacion, decreasing = TRUE)]))

all_tests_with_id$region <- factor(all_tests_with_id$region, levels = c(levels(pop_by_region$region), "No reportada"))

## Impute missing dates
all_tests_with_id[, date := fifelse(collectedDate > reportedDate, reportedDate, collectedDate)] ## if collectedDate is in the future make it reportedDate
all_tests_with_id[, date := fifelse(is.na(collectedDate), reportedDate - days(imputation_delay),  collectedDate)]
all_tests_with_id[, date := as_date(fifelse(!year(date) %in% the_years, reportedDate - days(imputation_delay),  date))]
all_tests_with_id <- all_tests_with_id[year(date) %in% the_years & date <= today()]
all_tests_with_id <- all_tests_with_id[order(date, reportedDate)]

saveRDS(all_tests_with_id , file = file.path(rda_path, "all_tests_with_id.rds"))

