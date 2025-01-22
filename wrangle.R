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

# -- Fixed values
pr_pop <- 3285874 ## population of puerto rico

icu_beds <- 229 #if available beds is missing change to this

first_day <- make_date(2020, 3, 12)

last_complete_day <- today() - 1

the_years <- seq(2020, year(today()))

age_levels <-  paste(seq(0, 125, 5), "to", seq(4, 129, 5))

imputation_delay  <- 2

alpha <- 0.05

test_types <- c("Molecular", "Antigens", "Molecular+Antigens")

original_test_types <- c("Molecular", "Antigens")

## Load latest data
message("Loading previous dataset.")

prev_all_tests_with_id <- readRDS(file = file.path(rda_path, "all_tests_with_id.rds"))

#molecular_last_download <- max(prev_all_tests_with_id[testType == "Molecular"]$orderCreatedAt) |>
molecular_last_download <- as_datetime(today(tz= "America/Puerto_Rico") - days(10)) |>
  with_tz(tzone = "GMT") |>
  format("%Y-%m-%dT%H:%M:%SZ")

#antigen_last_download <- max(prev_all_tests_with_id[testType == "Antigens"]$orderCreatedAt) |>
antigen_last_download <- as_datetime(today(tz= "America/Puerto_Rico") - days(10)) |>
  with_tz(tzone = "GMT") |>
  format("%Y-%m-%dT%H:%M:%SZ")

#createdAtStartDate=03-12-2020T04:00:00Z&createdAtEndDate=06-15-2022T04:00:00Z

## filter by date example: ?createdAtStartDate=2021-09-09T04:00:00Z&createdAtEndDate=2021-09-10T04:00:00Z
#cases_url <- "https://bioportal.salud.pr.gov/api/administration/reports/orders/basic"
cases_url <- "https://api-bioportal-prod-eastus2-01.azurewebsites.net/administration/reports/orders/basic"

cases_url_molecular <-  paste0(cases_url, 
                               "?testType=Molecular", 
                               "&createdAtStartDate=",
                               molecular_last_download,
                               "&createdAtEndDate=",
                               format(now(tzone="GMT")+days(1), "%Y-%m-%dT%H:%M:%SZ"))

cases_url_antigens <- paste0(cases_url, 
                             "?testType=Antigens",
                             "&createdAtStartDate=",
                             antigen_last_download,
                             "&createdAtEndDate=",
                             format(now(tzone="GMT")+days(1), "%Y-%m-%dT%H:%M:%SZ"))


get_bioportal <- function(url){
  setDT(jsonlite::fromJSON(
    rawToChar(
      httr::GET(url, httr::content_type('application/json'),
                httr::add_headers('Accept-Enconding'="br"))$content
  )))
}

# Reading and wrangling cases data from database ---------------------------
message("Reading case data.")

all_tests_with_id_molecular <- get_bioportal(cases_url_molecular)
all_tests_with_id_antigens <- get_bioportal(cases_url_antigens)
all_tests_with_id <- rbind(all_tests_with_id_molecular, all_tests_with_id_antigens)
rm(all_tests_with_id_molecular, all_tests_with_id_antigens); gc(); gc() 

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

all_tests_with_id$region <- factor(all_tests_with_id$region, levels = levels(prev_all_tests_with_id$region))

## Impute missing dates
all_tests_with_id[, date := fifelse(collectedDate > reportedDate, reportedDate, collectedDate)] ## if collectedDate is in the future make it reportedDate
all_tests_with_id[, date := fifelse(is.na(collectedDate), reportedDate - days(imputation_delay),  collectedDate)]
all_tests_with_id[, date := as_date(fifelse(!year(date) %in% the_years, reportedDate - days(imputation_delay),  date))]
all_tests_with_id <- all_tests_with_id[year(date) %in% the_years & date <= today()]
all_tests_with_id <- all_tests_with_id[order(date, reportedDate)]

## remove the replicates
all_tests_with_id <- funion(prev_all_tests_with_id, all_tests_with_id)[order(date, reportedDate)]

added_records <- pmax(nrow(all_tests_with_id) - nrow(prev_all_tests_with_id), 0)

if(added_records > 0){
  
  all_dates <- CJ(testType = test_types, date =  seq(first_day, max(all_tests_with_id$date), by = "day"))
  # -- Computing observed positivity rate
  ## adding a new test type that combines molecular and antigens
  mol_anti <- all_tests_with_id[date >= first_day & testType %in% c("Molecular", "Antigens") & result %in% c("positive", "negative")]
  mol_anti[, testType := "Molecular+Antigens"]
  
  
  ## compute daily totals
  tests <- rbind(all_tests_with_id, mol_anti)[date >= first_day & 
             testType %in% test_types & 
             result %in% c("positive", "negative")]
  tests <- tests[, .(people_positives = uniqueN(patientId[result == "positive"]),
                     people_total = uniqueN(patientId),
                     tests_positives = sum(result == "positive"),
                     tests_total = .N), by = c("testType", "date")]
  tests <- merge(tests, all_dates, by = c("testType", "date"), all.y = TRUE)
  tests[is.na(tests)] <- 0
  tests[, rate :=  people_positives / people_total]
  
  ## define function to compute weekly distinct cases
  ## and use this to compute percent of people with positive tests
  
  positivity <- function(dat){
    dat<-copy(dat)
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
  fits <- fits[, positivity(.SD), by = "testType"]
  
  ## add new variable to test data frame
  tests <- merge(tests, fits, by = c("testType", "date"), all.x = TRUE) 
  cols <- c("people_positives_week", "people_total_week")
  tests[, (cols) := lapply(.SD, replace_na, 0), .SDcols = cols]
  
  ## compute weekly totals for positive tests and total tests
  tests[, `:=`(tests_positives_week = sum7(tests_positives),
               tests_total_week = sum7(tests_total)), by = "testType"]
  
  # compute unique cases ------------------------------------------------------------
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
  
  ## Do reinfections by age here since instead of 
  age_starts <- c(0, 10, 15, 20, 30, 40, 65, 75)
  age_ends <- c(9, 14, 19, 29, 39, 64, 74, Inf)
  
  ## compute daily totals
  age_levels <- paste(age_starts, age_ends, sep = " a ")
  age_levels[length(age_levels)] <- paste0(age_starts[length(age_levels)],"+")
  
  reinfections <- copy(cases)
  reinfections[, `:=`(
    reinfection = days>0,
    age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
    age_end = as.numeric(str_extract(ageRange, "\\d+$")))] 
  reinfections[, ageRange := factor(age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))], levels = age_levels)]
  reinfections <- reinfections[, .(cases = .N), keyby =  .(testType, ageRange, reinfection, date)]
  reinfections$ageRange <- forcats::fct_na_value_to_level(reinfections$ageRange, "No reportada")
   
  reinfections <- merge(CJ(date=seq(first_day, today(), by ="day"), 
                           ageRange = factor(levels(reinfections$ageRange), levels = levels(reinfections$ageRange)),
                           testType = test_types, reinfection = c(TRUE, FALSE)),
                        reinfections, all.x = TRUE, 
                        by = c("testType",  "ageRange", "reinfection", "date"))
  reinfections[, cases := replace_na(cases, 0)]
  reinfections[, reinfection := fifelse(reinfection, "reinfection", "new")]
  reinfections <- dcast(reinfections, ...  ~ reinfection, value.var = "cases", fill = 0)
  reinfections <- reinfections[order(testType, date, ageRange)]
  reinfections[, new_week_avg := ma7(new), keyby = c("testType", "ageRange")]
  reinfections[, reinfection_week_avg := ma7(reinfection), keyby = c("testType", "ageRange")]
  
  if(FALSE){
    ##check with plot
  
    reinfections %>% filter(testType ==  "Molecular+Antigens") %>%
      ggplot(aes(date, reinfection)) + geom_col() + geom_line(aes(y=reinfection_week_avg)) + facet_wrap(~ageRange)
  
    reinfections %>% filter(testType ==  "Molecular+Antigens") %>%
      ggplot(aes(date, reinfection/(new+reinfection))) + geom_col() + 
      geom_line(aes(y=reinfection_week_avg/(new_week_avg+reinfection_week_avg))) +
                  facet_wrap(~ageRange)
    
    reinfections %>% 
      filter(date>=make_date(2021,7,1)) %>%
      filter(testType ==  "Molecular+Antigens") %>%
      group_by(date) %>%
      summarize(cases=sum(reinfection), .groups = "drop") %>% 
      ggplot(aes(date, cases)) + 
      geom_col() +
      theme_bw()
  }
  
  cases <- cases[, .(cases = .N), keyby = .(testType, date)]
    
  # Make sure all dates are included
  cases <-  merge(tests[, .(testType, date)], cases, by = c("testType", "date"), all.x = TRUE) 
  cases[, cases := replace_na(cases, 0)]
    
  # compute daily weekly average and add to cases data frame
  cases[, cases_week_avg := ma7(cases), by = "testType"]
  
  ## add new cases and weekly average to tests data frame
  tests <- merge(tests, cases, by = c("testType", "date"), all.x = TRUE)
  tests[, cases_plus_negatives := (people_total_week - people_positives_week + cases_week_avg * 7)]
  tests[, `:=`(cases_rate = cases_week_avg * 7 / cases_plus_negatives,
               cases_plus_negatives_daily = people_total - people_positives + cases)]
  tests[, cases_rate_daily := cases / cases_plus_negatives_daily]
           
  ## Compute unique negatives
  
  # compute unique negative cases ------------------------------------------------------------
  message("Computing unique negatives.")
  
  negative_cases <- rbind(all_tests_with_id, mol_anti)[
    date>=first_day & result == "negative" &
             testType %in% test_types]
  negative_cases <- negative_cases[order(date)]
  negative_cases <- negative_cases[, .SD[1], keyby = .(testType, patientId)]
  negative_cases[, patientId:=NULL]
  negative_cases[, result:=NULL]
  negative_cases <- negative_cases[, .(negative_cases = .N), keyby = .(testType, date)]
  negative_cases <- negative_cases[order(testType, date)]
    
  # Compute daily weekly average and add to negative_cases data frame
  # but first make sure all dates are included
  negative_cases <-  merge(tests[, .(testType, date)], negative_cases, by = c("testType", "date"), all.x = TRUE)
  negative_cases[, negative_cases := replace_na(negative_cases, 0)]
    
  # compute daily weekly average and add to negative_cases data frame
  negative_cases[, negative_cases_week_avg := ma7(negative_cases), by = "testType"]
  
  ## add new cases and weekly average to tests data frame
  tests <- merge(tests, negative_cases, by = c("testType", "date"), all.x=TRUE)
  
  ## the following are diagnostic plots
  if(FALSE){
    library(scales)
    
    source("dashboard/functions.R")
    lag_to_complete <- 7
    last_day <- today() - days(lag_to_complete)
    
    ## check positivity rate
    
    plot_positivity(tests, first_day, today(), type = "Molecular") +
      geom_smooth(method = "loess", formula = "y~x", span = 0.05, method.args = list(degree = 1, weight = tests$tests), color = "red", lty =2, fill = "pink") 
    
    plot_positivity(tests, first_day, today(), type = "Molecular") 
    ## check test plot
    plot_test(tests, first_day, today())
    plot_test(tests, first_day, today(), type  = "Antigens")
    plot_test(tests, first_day, today(), type  = "Molecular+Antigens")
  
    ## check cases plot
    ys <- TRUE
    plot_cases(cases, yscale = ys)
    plot_cases(cases, first_day, today(), type  = "Antigens", yscale = ys)
    plot_cases(cases, first_day, today(), type  = "Molecular+Antigens", yscale = ys)
    
  }

} else{ 
  load(file.path(rda_path, "data.rda"))
} ## if no new records, tests or cases not created so need to load

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
url <- "https://covid19datos.salud.pr.gov/estadisticas_v2/download/data/sistemas_salud/completo"
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
           VENT_ORD, VENT_REC, VENT_ENTR, HospitCOV19, CamasICU, CamasICU_disp)
  hosp_mort <- full_join(hosp_mort, old_hosp_mort, by = "date") %>%
    arrange(date) %>%
    mutate(HospitCOV19 = ifelse(is.na(HospitCOV19.x), HospitCOV19.y, HospitCOV19.x),
           CamasICU = ifelse(is.na(CamasICU.x), CamasICU.y, CamasICU.x),
           CamasICU_disp = ifelse(is.na(CamasICU_disp.x), CamasICU_disp.y, CamasICU_disp.x)) %>%
    select(-contains(".x"), -contains(".y"))
}
# -- seven day averages 
hosp_mort$hosp_week_avg <- ma7(hosp_mort$HospitCOV19)

hosp_mort$icu_week_avg <- ma7(hosp_mort$CamasICU)
  
ind <- which(!is.na(hosp_mort$CAMAS_PED_COVID))
hosp_mort$ped_hosp_week_avg <- rep(NA, nrow(hosp_mort))
hosp_mort$ped_hosp_week_avg[ind] <- ma7(hosp_mort$CAMAS_PED_COVID[ind])

ind <- which(!is.na(hosp_mort$CAMAS_PICU_COVID))
hosp_mort$picu_week_avg <- rep(NA, nrow(hosp_mort))
hosp_mort$picu_week_avg[ind] <- ma7(hosp_mort$CAMAS_PICU_COVID[ind])


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

#url <- "https://bioportal.salud.pr.gov/api/administration/reports/deaths/summary"
#url <- "https://api-bioportal-prod-eastus2-01.azurewebsites.net/administration/reports/deaths/summary"
url <- "https://biostatistics.salud.pr.gov/deaths/covid-19/minimal"
deaths <- try({
  get_bioportal(url) %>% 
    mutate(date = ymd(deathDate)) %>%
     mutate(date = if_else(date < first_day | date > today(), 
                           ymd(deathReportDate),
                           date)) %>%
     mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
            age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
     mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
     mutate(ageRange = factor(ageRange, levels = age_levels)) %>%
     rename(region = physicalRegion, reportDate = deathReportDate)
       
  
  ##old code for https://api-bioportal-prod-eastus2-01.azurewebsites.net/administration/reports/deaths/summary
  # jsonlite::fromJSON(url) %>%
  #   mutate(date = as_date(ymd_hms(deathDate, tz = "America/Puerto_Rico"))) %>%
  #   mutate(date = if_else(date < first_day | date > today(), 
  #                         as_date(ymd_hms(reportDate, tz = "America/Puerto_Rico")),
  #                         date)) %>%
  #   mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
  #          age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  #   mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  #   mutate(ageRange = factor(ageRange, levels = age_levels)) 
  #   get_bioportal(url) %>%
  # mutate(date = as_date(ymd_hms(deathDate, tz = "America/Puerto_Rico"))) %>%
})

url <- "https://covid19datos.salud.pr.gov/estadisticas_v2/download/data/defunciones/completo"
dashboard_deaths <- try({
  read.csv(text = rawToChar(httr::content(httr::GET(url)))) %>% 
    rename(ageRange = TX_GRUPO_EDAD, date = FE_MUERTE) %>%
    mutate(date = as_date(ymd_hms(date, tz = "America/Puerto_Rico"))) %>%
    mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
           age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
    mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
    mutate(ageRange = factor(ageRange, levels = age_levels)) 
})

if(class(deaths)[1] == "try-error"){
  if(class(dashboard_deaths)[1] == "try-error"){
	load(file.path(rda_path, "deaths.rda"))
  } else{
	## if no dashboard data replace the death data with BioPortal data 
  	deaths <- dashboard_deaths
  }
}

## make sure all dates appear, even if hosp are not being reported
all_dates <- data.table(date = seq(first_day, max(last_complete_day, max(hosp_mort$date)), by = "day"))
hosp_mort <- left_join(all_dates, hosp_mort, by = "date")
hosp_mort <- deaths %>%
  group_by(date) %>%
  summarize(deaths = n(), .groups = "drop") %>%
  full_join(hosp_mort, by = "date") %>%
  arrange(date) %>%
  mutate(deaths = replace_na(deaths,0)) %>%
  mutate(IncMueSalud = deaths,
         mort_week_avg =  ma7(deaths)) %>%
  select(-deaths)
rm(all_dates)

## Rezagos muerte

if(class(deaths)[1] != "try-error"){
  rezago_mort <- deaths %>% 
  	filter(!is.na(date)) %>%
  	mutate(bulletin_date = ymd(reportDate)) %>%
  	arrange(date, bulletin_date) %>%
  	mutate(diff = (as.numeric(bulletin_date) - as.numeric(date))) %>%
  	select(date, diff)
} else{
  load(file.path(rda_path, "rezago_mort.rda"))
}

## Save results
## define date and time of latest download
the_stamp <- now(tzone="America/Puerto_Rico")

tests <- as.data.frame(tests)
cases <- as.data.frame(cases)
hosp_mort <- as.data.frame(hosp_mort)

save(first_day, last_complete_day, added_records,
     alpha, the_stamp, 
     tests, cases,
     hosp_mort, pr_pop, 
     file = file.path(rda_path, "data.rda"))

if(added_records>0){
  reinfections <- as.data.frame(reinfections)
  save(reinfections,  file = file.path(rda_path, "reinfections.rda"))
}

## save this as backup in case salud dashboard down
save(hosp_mort, file = file.path(rda_path, "hosp_mort.rda"))

rezago_mort <- as.data.frame(rezago_mort)
save(rezago_mort, file = file.path(rda_path, "rezago_mort.rda"))

## for use in wrangle-by-strata.R
deaths <- as.data.frame(deaths)
save(deaths, last_complete_day, file = file.path(rda_path, "deaths.rda"))

#if(added_records>0) saveRDS(all_tests_with_id, file = file.path(rda_path, "all_tests_with_id.rds"))

