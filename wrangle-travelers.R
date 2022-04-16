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

# -- Fixed values

first_day <- make_date(2020, 3, 12)

the_years <- seq(2020, year(today()))

alpha <- 0.05

## add passanger data

message("Computing traveler statistics.")

url <- "https://BioPortal.salud.pr.gov/api/administration/reports/travels/total-forms-by-reported-arrival-date"

travelers <- jsonlite::fromJSON(url) %>%
  mutate(date = mdy(date)) %>%
  filter(year(date)>=2020) %>%
  filter(date <= today()) %>%
  mutate(residents_week_avg =  ma7(date, residents)$moving_avg,
         perc_residents = percentageResidentsArrivedWithNegativePcrResults/100,
         perc_residents_week_avg = ma7(date, perc_residents*residents)$moving_avg/residents_week_avg,
         
         short = nonResidentsStayingLessThan5Days,
         short_week_avg =  ma7(date, short)$moving_avg,
         perc_short = percentageNonResidentsStayingLessThan5DaysArrivedWithNegativePcrResults/100,
         perc_short_week_avg  = ma7(date, perc_short*short)$moving_avg/short_week_avg,
         
         long = nonResidentsStaying5DaysOrMore,
         long_week_avg =  ma7(date, nonResidentsStaying5DaysOrMore)$moving_avg,
         perc_long = percentageResidentsStaying5DaysOrMoreArrivedWithNegativePcrResults/100,
         perc_long_week_avg  = ma7(date, perc_long*long)$moving_avg/long_week_avg,
         
         vaccine = selfReportedTotalArrivedVaccinatedBothDoses / total,
         vaccine_week_avg = ma7(date, selfReportedTotalArrivedVaccinatedBothDoses)$moving_avg / 
           ma7(date, total)$moving_avg ) %>%
  select(date, 
         residents, perc_residents, short, perc_short, long, perc_long,
         residents_week_avg, perc_residents_week_avg, short_week_avg, perc_short_week_avg, long_week_avg, perc_long_week_avg, 
         vaccine, vaccine_week_avg)

save(travelers, file = file.path(rda_path, "travelers.rda"))

