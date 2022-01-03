# combine fig2 data


fig2_files <- dir(recursive=T,pattern="fig2")

all_fig2_data <- NULL

for(file_ind in seq_along(fig2_files)){
  
  this_fig2 <- readRDS(file = fig2_files[file_ind])
  
  all_fig2_data <- rbind(all_fig2_data,this_fig2)
  
}


saveRDS(object = all_fig2_data, file = "all_data_for_figure2.rds")

