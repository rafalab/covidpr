# -- Libraries   
library(tidyverse)
library(lubridate)
library(splines)

# fit glm spline ----------------------------------------------------------
# no longer used. we now use moving average to match other dashboards
# spline_fit <- function(d, y, n = NULL, 
#                        week_effect = TRUE, 
#                        knots_per_month = 2, 
#                        family = quasibinomial, 
#                        alpha = 0.05){
#   
#   z <- qnorm(1 - alpha/2)
#   
#   x <- as.numeric(d)
#   
#   df  <- round(knots_per_month * length(x) / 30) + 1
#   
#   if(family()$family %in% c("binomial", "quasibinomial")){
#     if(is.null(n)) stop("Must supply n with binomial or quasibinomial")
#     y <- cbind(y, n-y)
#   }
#   
#   if(week_effect){
#     
#     w <- factor(wday(d))
#     contrasts(w) <- contr.sum(length(levels(w)), contrasts = TRUE)
#     w <- model.matrix(~w)[,-1]
#     
#     glm_fit  <- glm(y ~ ns(x, df = df, intercept = TRUE) + w - 1, family = family)
#     
#   } else {
#     
#     glm_fit  <- glm(y ~ ns(x, df = df, intercept = TRUE) - 1, family = family)
#     
#   }
#   
#   glm_pred <- predict(glm_fit, type = "terms", se.fit = TRUE)
#   
#   fit <- family()$linkinv(glm_pred$fit[,1])
#   
#   lower <- family()$linkinv(glm_pred$fit[,1] - z * glm_pred$se.fit[,1])
#   
#   upper <- family()$linkinv(glm_pred$fit[,1] + z * glm_pred$se.fit[,1])
#   
#   return(tibble(date = d, fit = fit, lower = lower, upper = upper))  
# }

## if on server, save with full path
## if not on server, save to home directory
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
pr_pop <- 3285874 ## population of puerto rico

icu_beds <- 229 #if available beds is missing change to this

first_day <- make_date(2020, 3, 12)

last_complete_day <- today() - 1

the_years <- seq(2020, year(today()))

age_levels <-  paste(seq(0, 125, 5), "to", seq(4, 129, 5))

imputation_delay  <- 2

alpha <- 0.05

## Load latest data
message("Loading previous dataset.")

prev_all_tests_with_id <- readRDS(file = file.path(rda_path, "all_tests_with_id.rds"))

molecular_last_download <- filter(prev_all_tests_with_id, testType == "Molecular") %>%
  pull(resultCreatedAt) %>% max() %>% with_tz(tzone = "GMT") %>%
  str_replace(" ", "T") %>%
  paste0("Z")

antigen_last_download <- filter(prev_all_tests_with_id, testType == "Antigens") %>%
  pull(resultCreatedAt) %>% max() %>% with_tz(tzone = "GMT") %>%
  str_replace(" ", "T") %>%
  paste0("Z")

## filter by date example: ?createdAtStartDate=2021-09-09T04:00:00Z&createdAtEndDate=2021-09-10T04:00:00Z
cases_url <- "https://bioportal.salud.pr.gov/api/administration/reports/orders/basic"

cases_url_molecular <-  paste0(cases_url, 
                               "?testType=Molecular", 
                               "&createdAtStartDate=",
                               molecular_last_download)
cases_url_antigens <- paste0(cases_url, 
                             "?testType=Antigens",
                             "&createdAtStartDate=",
                             antigen_last_download)

get_bioportal <- function(url){
  jsonlite::fromJSON(
    rawToChar(
      httr::GET(url, httr::content_type('application/json'),
                httr::add_headers('Accept-Enconding'="br"))$content)
  )
}

# get_bioportal_2 <- function(url){
#   y <- rawToChar(
#     httr::GET(url, httr::content_type('application/json'),
#               httr::add_headers('Accept-Enconding'="br"))$content, multiple = TRUE)
#   starts <- which(y=="{")
#   ends <- c(starts[-1]-2,length(y)-1)
#   y <- sapply(seq_along(starts), function(i) {
#     ind <- (starts[i]):(ends[i])
#     tmp <- str_c(y[ind], collapse="")
#     if(str_detect(tmp, "Serological")) return("") else return(tmp)
#   })
#   jsonlite::fromJSON(paste0("[", str_c(y[y!=""], collapse=","), "]"))
# }
# 
# get_bioportal <- function(url){
#   y <- rawToChar(
#     httr::GET(url, httr::content_type('application/json'),
#               httr::add_headers('Accept-Enconding'="br"))$content, multiple = TRUE)
#   starts <- which(y == "{")
#   ends <- c(starts[-1]-2, length(y) - 1)
#   y <- sapply(seq_along(starts), function(i) {
#     ind <- (starts[i]):(ends[i])
#     return(str_c(y[ind], collapse=""))
#   })
#   jsonlite::fromJSON(paste0("[", str_c(y, collapse=","), "]"))
# }


#test_types <- c("Molecular", "Serological", "Antigens", "Molecular+Antigens")
#original_test_types <- c("Molecular", "Serological", "Antigens")
test_types <- c("Molecular", "Antigens", "Molecular+Antigens")
original_test_types <- c("Molecular", "Antigens")

# Reading and wrangling cases data from database ---------------------------
message("Reading case data.")

all_tests_with_id_molecular <- get_bioportal(cases_url_molecular)
all_tests_with_id_antigens <- get_bioportal(cases_url_antigens)
all_tests_with_id <- rbind(all_tests_with_id_molecular, all_tests_with_id_antigens)
rm(all_tests_with_id_molecular, all_tests_with_id_antigens); gc(); gc()

message("Processing case data.")

all_tests_with_id <- all_tests_with_id %>%  
  as_tibble() %>%
  mutate(testType = str_to_title(testType),
         testType = ifelse(testType == "Antigeno", "Antigens", testType),
         collectedDate = ymd_hms(collectedDate, tz = "America/Puerto_Rico"),
         reportedDate = ymd_hms(reportedDate, tz = "America/Puerto_Rico"),
         orderCreatedAt = ymd_hms(orderCreatedAt, tz = "America/Puerto_Rico"),
         resultCreatedAt = ymd_hms(resultCreatedAt, tz = "America/Puerto_Rico"),
         ageRange       = na_if(ageRange, "N/A"),
         ageRange       = factor(ageRange, levels = age_levels),
         region = na_if(region, "N/A"),
         region = ifelse(region == "Bayamon", "Bayamón", region),
         region = ifelse(region == "Mayaguez", "Mayagüez", region),
         region = replace_na(region, "No reportada"),
         region = factor(region),
         result = tolower(result),
         result = case_when( 
           (grepl("covid", result) | grepl("sars-cov-2", result)) &
             grepl("positive", result) ~ "positive",
           grepl("influenza", result) ~ "other",
                            grepl("positive", result) ~ "positive",
                            grepl("negative", result) ~ "negative",
                            result == "not detected" ~ "negative",
                            TRUE ~ "other")) %>%
  arrange(reportedDate, collectedDate, patientId) %>%
  filter(testType %in% test_types)

## fixing bad dates: if you want to remove bad dates instead, change FALSE TO TRUE
if(FALSE){
  ## remove bad dates
  all_tests_with_id <- all_tests_with_id %>% 
    filter(!is.na(collectedDate) & year(collectedDate) %in% the_years & collectedDate <= today()) %>%
    mutate(date = as_date(collectedDate))
} else{
  ## Impute missing dates
  all_tests_with_id <- all_tests_with_id %>% 
    mutate(date = if_else(collectedDate > reportedDate, reportedDate, collectedDate)) %>% ## if collectedDate is in the future make it reportedDate
    mutate(date = if_else(is.na(collectedDate), reportedDate - days(imputation_delay),  collectedDate)) %>%
    mutate(date = if_else(!year(date) %in% the_years, reportedDate - days(imputation_delay),  date)) %>%  
    mutate(date = as_date(date)) %>%
    filter(year(date) %in% the_years & date <= today()) %>%
    arrange(date, reportedDate)
}

## remove the replicates
all_tests_with_id <- union(prev_all_tests_with_id, all_tests_with_id) %>% arrange(date, reportedDate)

added_records <- pmax(nrow(all_tests_with_id) - nrow(prev_all_tests_with_id), 0)

if(added_records>0){

# -- Computing observed positivity rate
## adding a new test type that combines molecular and antigens
mol_anti <- all_tests_with_id %>%
  filter(date >= first_day & testType %in% c("Molecular", "Antigens") & 
           result %in% c("positive", "negative")) %>%
  mutate(testType = "Molecular+Antigens") 


## compute daily totals
tests <- all_tests_with_id %>%
  bind_rows(mol_anti) %>%
  filter(date >= first_day & 
           testType %in% test_types & 
           result %in% c("positive", "negative")) %>%
  group_by(testType, date) %>%
  summarize(people_positives = n_distinct(patientId[result == "positive"]),
            people_total = n_distinct(patientId),
            tests_positives = sum(result == "positive"),
            tests_total = n(),
            .groups = "drop") %>%
  mutate(rate = people_positives / people_total)

## define function to compute weekly distinct cases
## and use this to compute percent of people with positive tests

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

message("Computing positivity.")

## run the function on each test type
fits <- all_tests_with_id %>% 
  bind_rows(mol_anti) %>%
  mutate(entry_date = as_date(orderCreatedAt)) %>%
  filter(date >= first_day & testType %in% test_types & 
           result %in% c("positive", "negative")) %>%
  nest_by(testType) %>%
  summarize(positivity(data), .groups = "drop")
  
## add new variable to test data frame
tests <- left_join(tests, fits, by = c("testType", "date")) 

## compute weekly totals for positive tests and total tests
tests <- tests %>% 
  group_by(testType) %>%
  mutate(tests_positives_week = sum7(d = date, y = tests_positives)$moving_sum) %>%
  mutate(tests_total_week = sum7(d = date, y = tests_total)$moving_sum) %>%
  ungroup()

# compute unique cases ------------------------------------------------------------
cases <- all_tests_with_id %>%  
  bind_rows(mol_anti) %>%
  filter(date>=first_day & result == "positive" &
           testType %in% test_types) %>%
  arrange(testType, patientId, date) %>%
  group_by(testType, patientId) %>% ##newId takes reinfection into account
  mutate(days = c(0, floor(diff(as.numeric(date))/90))) %>%
  mutate(newId = paste0(patientId, "-", days)) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  arrange(testType, newId, date) %>%
  group_by(testType, newId) %>%
  slice(1) %>% 
  ungroup() %>%
  select(-patientId, -result) %>%
  arrange(testType, date) 

## Do reinfections by age here since instead of 
age_starts <- c(0, 10, 15, 20, 30, 40, 65, 75)
age_ends <- c(9, 14, 19, 29, 39, 64, 74, Inf)

## compute daily totals
age_levels <- paste(age_starts, age_ends, sep = " a ")
age_levels[length(age_levels)] <- paste0(age_starts[length(age_levels)],"+")

reinfections <- cases %>%
  mutate(reinfection = days>0) %>%
  mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
         age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange, levels = age_levels)) %>%
  group_by(testType, ageRange, reinfection, date) %>%
  summarize(cases = n(), .groups = "drop") %>%
  right_join(crossing(date=seq(first_day, today(), by ="day"), 
                      ageRange = factor(age_levels, levels = age_levels),
                      testType = test_types, reinfection = c(TRUE, FALSE)), 
             by = c("testType",  "ageRange", "reinfection", "date")) %>%
  mutate(cases = replace_na(cases, 0)) %>%
  arrange(testType, ageRange, cases, reinfection)

if(FALSE){
  ##check with plot

  reinfections %>% filter(testType ==  "Molecular+Antigens" & reinfection) %>%
    ggplot(aes(date, cases)) + geom_col() + facet_wrap(~ageRange)

  reinfections %>% 
    filter(date>=make_date(2021,7,1)) %>%
    filter(testType ==  "Molecular+Antigens" & reinfection) %>%
    group_by(date) %>%
    summarize(cases=sum(cases), .groups = "drop") %>% 
    ggplot(aes(date, cases)) + 
    geom_col() +
    theme_bw()
}

cases <- cases %>%
  group_by(testType, date) %>% 
  summarize(cases = n(), .groups = "drop") 

# Make sure all dates are included
cases <-  left_join(select(tests, testType, date), cases, by = c("testType", "date")) %>%
  replace_na(list(cases = 0)) 

# compute daily weekly average and add to cases data frame
fits <- cases %>% 
  group_by(testType) %>%
  do(ma7(d = .$date, y = .$cases)) %>%
  rename(cases_week_avg = moving_avg)
cases <- left_join(cases, fits, by = c("testType", "date"))

## add new cases and weekly average to tests data frame
tests <- left_join(tests, cases, by = c("testType", "date")) %>%
  mutate(cases_plus_negatives = (people_total_week - people_positives_week + cases_week_avg * 7),
         cases_rate = cases_week_avg * 7 / cases_plus_negatives,
         cases_plus_negatives_daily = people_total - people_positives + cases,
         cases_rate_daily = cases / cases_plus_negatives_daily)
         
## Compute unique negatives

# compute unique negative cases ------------------------------------------------------------
message("Computing unique negatives.")

negative_cases <- all_tests_with_id %>%  
  bind_rows(mol_anti) %>%
  filter(date>=first_day & result == "negative" &
           testType %in% test_types) %>%
  group_by(testType, patientId) %>%
  arrange(date) %>%
  slice(1) %>% 
  ungroup() %>%
  select(-patientId, -result) %>%
  arrange(testType, date) %>%
  group_by(testType, date) %>% 
  summarize(negative_cases = n(), .groups = "drop")

# Make sure all dates are included
negative_cases <-  select(tests, testType, date) %>% 
  left_join(negative_cases, by = c("testType", "date")) %>%
  replace_na(list(negative_cases = 0))

# compute daily weekly average and add to negative_cases data frame
fits <- negative_cases %>% 
  group_by(testType) %>%
  do(ma7(d = .$date, y = .$negative_cases)) %>%
  rename(negative_cases_week_avg = moving_avg)
negative_cases <- left_join(negative_cases, fits, by = c("testType", "date"))

## add new cases and weekly average to tests data frame
tests <- left_join(tests, negative_cases, by = c("testType", "date"))

## the following are diagnostic plots
if(FALSE){
  library(scales)
  
  source("functions.R")
  lag_to_complete <- 7
  last_day <- today() - days(lag_to_complete)
  
  ## check positivity rate
  
  plot_positivity(tests, first_day, today(), type = "Molecular", show.all = FALSE) +
    geom_smooth(method = "loess", formula = "y~x", span = 0.2, method.args = list(degree = 1, weight = tests$tests), color = "red", lty =2, fill = "pink") 
  
  plot_positivity(tests, first_day, today(), type = "Molecular", show.all = TRUE) 
  ## check test plot
  plot_test(tests, first_day, today())
  plot_test(tests, first_day, today(), type  = "Serological")
  plot_test(tests, first_day, today(), type  = "Antigens")
  plot_test(tests, first_day, today(), type  = "Molecular+Antigens")

  ## check cases plot
  ys <- TRUE
  plot_cases(cases, yscale = ys)
  plot_cases(cases, first_day, today(), type  = "Serological", yscale = ys)
  plot_cases(cases, first_day, today(), type  = "Antigens", yscale = ys)
  plot_cases(cases, first_day, today(), type  = "Molecular+Antigens", yscale = ys)
  
}

} else{ load(file.path(rda_path, "data.rda"))} ## if no new records, tests or cases not created so need to load

# --Mortality and hospitlization
# use old handmade database to fill in the blanks
old_hosp_mort <- read_csv("https://raw.githubusercontent.com/rafalab/covidpr/main/dashboard/data/DatosMortalidad.csv") %>%
  mutate(date = mdy(Fecha)) %>%
  filter(date >= first_day) %>%
  arrange(date) %>%
  select(date, HospitCOV19, CamasICU_disp, CamasICU)
# we started keeping track of available beds on 2020-09-20
# hosp_mort <- hosp_mort %>%
#   replace_na(list(CamasICU_disp = icu_beds))

httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))
url <- "https://covid19datos.salud.gov.pr/estadisticas_v2/download/data/sistemas_salud/completo"
hosp_mort <- try({
  read.csv(text = rawToChar(httr::content(httr::GET(url)))) %>% 
    mutate(date = as_date(FE_REPORTE)) %>%
    filter(date >= first_day) %>%
    full_join(old_hosp_mort, by = "date") %>%
    arrange(date) %>%
    ## add columns to match old table
    mutate(HospitCOV19 = ifelse(is.na(CAMAS_ADULTOS_COVID), HospitCOV19, CAMAS_ADULTOS_COVID),
           CamasICU = ifelse(is.na(CAMAS_ICU_COVID), CamasICU, CAMAS_ICU_COVID),
           CamasICU_disp = ifelse(is.na(CAMAS_ICU_DISP), CamasICU_disp, CAMAS_ICU_DISP))
})

if(class(hosp_mort)[1] == "try-error"){
  load(file.path(rda_path, "hosp_mort.rda"))
  hosp_mort <- hosp_mort %>% 
    select(date, FE_REPORTE, CAMAS_ADULTOS_COVID, CAMAS_ADULTOS_NOCOVID, CAMAS_ADULTOS_OCC, 
           CAMAS_ADULTOS_DISP, CAMAS_ADULTOS_TOTAL, CAMAS_ICU_COVID, CAMAS_ICU_NOCOVID, CAMAS_ICU_OCC, 
           CAMAS_ICU_DISP, CAMAS_ICU_TOTAL, CAMAS_PED_COVID, CAMAS_PED_NOCOVID, CAMAS_PED_OCC, 
           CAMAS_PED_DISP, CAMAS_PED_TOTAL, CAMAS_PICU_COVID, CAMAS_PICU_NOCOVID, CAMAS_PICU_OCC, 
           CAMAS_PICU_DISP, CAMAS_PICU_TOTAL, VENT_ADULTOS_COVID, VENT_ADULTOS_NOCOVID, VENT_ADULTOS_OCC, 
           VENT_ADULTOS_DISP, VENT_ADULTOS_TOTAL, VENT_PED_COVID, VENT_PED_NOCOVID, VENT_PED_OCC, 
           VENT_PED_DISP, VENT_PED_TOTAL, CUARTOS_PRESNEG_OCC, CUARTOS_PRESNEG_DISP, CUARTOS_PRESNEG_TOTAL, 
           VENT_ORD, VENT_REC, VENT_ENTR, CONVALECIENTES, HospitCOV19, CamasICU, CamasICU_disp)
  hosp_mort <- full_join(hosp_mort, old_hosp_mort, by = "date") %>%
    arrange(date) %>%
    mutate(HospitCOV19 = ifelse(is.na(HospitCOV19.x), HospitCOV19.y, HospitCOV19.x),
           CamasICU = ifelse(is.na(CamasICU.x), CamasICU.y, CamasICU.x),
           CamasICU_disp = ifelse(is.na(CamasICU_disp.x), CamasICU_disp.y, CamasICU_disp.x)) %>%
    select(-contains(".x"), -contains(".y"))
}
# -- seven day averages 
# deaths gets added later
# fits <- with(hosp_mort, 
#              ma7(d = date, y = IncMueSalud))
# hosp_mort$mort_week_avg <- fits$moving_avg

fits <- with(hosp_mort, 
             ma7(d = date, y = HospitCOV19))
hosp_mort$hosp_week_avg <- fits$moving_avg

fits <- with(hosp_mort, 
             ma7(d = date, y = CamasICU))
hosp_mort$icu_week_avg <- fits$moving_avg

ind <- which(!is.na(hosp_mort$CAMAS_PED_COVID))
fits <- with(hosp_mort[ind,], 
             ma7(d = date, y = CAMAS_PED_COVID))
hosp_mort$ped_hosp_week_avg <- rep(NA, nrow(hosp_mort))
hosp_mort$ped_hosp_week_avg[ind] <- fits$moving_avg

ind <- which(!is.na(hosp_mort$CAMAS_PICU_COVID))
fits <- with(hosp_mort[ind,], 
             ma7(d = date, y = CAMAS_PICU_COVID))
hosp_mort$picu_week_avg <- rep(NA, nrow(hosp_mort))
hosp_mort$picu_week_avg[ind] <- fits$moving_avg


## Vaccine data
## We will add it to hosp_mort data
url <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv"
vaccines <- read_csv(url) %>% 
  filter(location == "Puerto Rico") %>%
  select(date, total_distributed, total_vaccinations, people_vaccinated, people_fully_vaccinated) %>%
  arrange(date) 

## fill in NAs
for(j in which(names(vaccines)!="date")){
  for(i in 2:nrow(vaccines)){
    if(is.na(vaccines[[j]][i])) vaccines[[j]][i] <- max(vaccines[[j]][1:(i-1)], na.rm=TRUE)
  }
}

## fill in the NAs
hosp_mort <- full_join(hosp_mort, vaccines, by = "date") 

## downloading deaths by age
age_starts <- c(0, 10, 15, 20, 30, 40, 65, 75)
age_ends <- c(9, 14, 19, 29, 39, 64, 74, Inf)

## compute daily totals
age_levels <- paste(age_starts, age_ends, sep = " a ")
age_levels[length(age_levels)] <- paste0(age_starts[length(age_levels)],"+")

message("Computing Deaths")

url <- "https://bioportal.salud.pr.gov/api/administration/reports/deaths/summary"

bioportal_deaths <- jsonlite::fromJSON(url) %>%
  mutate(date = as_date(ymd_hms(deathDate, tz = "America/Puerto_Rico"))) %>%
  mutate(date = if_else(date < first_day | date > today(), 
                        as_date(ymd_hms(reportDate, tz = "America/Puerto_Rico")),
                        date)) %>%
  mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
         age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange, levels = age_levels)) 

## replace the death data with latest from dashboard

url <- "https://covid19datos.salud.gov.pr/estadisticas_v2/download/data/defunciones/completo"
deaths <- try({
  read.csv(text = rawToChar(httr::content(httr::GET(url)))) %>% 
    rename(ageRange = TX_GRUPO_EDAD, date = FE_MUERTE) %>%
    mutate(date = as_date(ymd_hms(date, tz = "America/Puerto_Rico"))) %>%
    mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
           age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
    mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
    mutate(ageRange = factor(ageRange, levels = age_levels)) 
})

if(class(deaths)[1] == "try-error"){
  ## if no dashboard data replace the death data with BioPortal data 
  deaths <- bioportal_deaths
}

hosp_mort <- deaths %>%
  group_by(date) %>%
  summarize(deaths = n(), .groups = "drop") %>%
  full_join(hosp_mort, by = "date") %>%
  arrange(date) %>%
  mutate(deaths = replace_na(deaths,0)) %>%
  mutate(IncMueSalud = deaths,
         mort_week_avg =  ma7(date, deaths)$moving_avg) %>%
  select(-deaths)


## Rezagos muerte

rezago_mort <- bioportal_deaths %>% 
  filter(!is.na(date)) %>%
  mutate(bulletin_date = as_date(ymd_hms(reportDate, tz = "America/Puerto_Rico"))) %>%
  arrange(date, bulletin_date) %>%
  mutate(diff = (as.numeric(bulletin_date) - as.numeric(date))) %>%
  select(date, diff)

## Save results

## define date and time of latest download
the_stamp <- now(tzone="America/Puerto_Rico")
  
save(first_day, last_complete_day, added_records,
     alpha, the_stamp, 
     tests, cases,
     hosp_mort, pr_pop, 
     file = file.path(rda_path, "data.rda"))

save(reinfections,  file = file.path(rda_path, "reinfections.rda"))

     ## save this as backup in case salud dashboard down
save(hosp_mort, file = file.path(rda_path, "hosp_mort.rda"))

save(rezago_mort, file = file.path(rda_path, "rezago_mort.rda"))

## for use in wrangle-by-strata.R
save(deaths, last_complete_day, file = file.path(rda_path, "deaths.rda"))

saveRDS(all_tests_with_id, file = file.path(rda_path, "all_tests_with_id.rds"))



