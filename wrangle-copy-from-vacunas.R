if(grepl("fermat|ask2me-phys", Sys.info()["nodename"])){
  rda_path <- "/homes10/rafa/dashboard/covidpr/rdas"
} else{
  rda_path <- "rdas"
}
load(url("https://github.com/rafalab/vacunaspr/raw/main/rdas/tabs.rda"))
vacunas_summary_tab <- summary_tab
save(vacunas_summary_tab, file = file.path(rda_path, "vacunas_summary_tab.rda"))
