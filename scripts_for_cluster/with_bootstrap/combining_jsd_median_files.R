

jsd_median_files <- dir(pattern = "JSD_results_median",recursive = T,full.names = T)

message("Combining these files:\n")
message(paste0(jsd_median_files,"\n"))

all_median_jsd <- purrr::map_dfr(jsd_median_files, readRDS)

saveRDS(object = all_median_jsd,file = "combined_JSD_medians.rds")

