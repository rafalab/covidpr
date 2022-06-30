library(tidyverse)
library(data.table)
## if on server, save with full path
## if not on server, save to home directory
if(grepl("fermat|leo", Sys.info()["nodename"])){
  rda_path <- "/homes10/rafa/dashboard/covidpr/rdas"
} else{
  rda_path <- "rdas"
}
prev_all_tests_with_id <- readRDS(file = file.path(rda_path, "all_tests_with_id.rds"))
prev_all_tests_with_id <- distinct(prev_all_tests_with_id)
prev_all_tests_with_id <- setDT(prev_all_tests_with_id)

saveRDS(prev_all_tests_with_id , file = file.path(rda_path, "all_tests_with_id.rds"))

