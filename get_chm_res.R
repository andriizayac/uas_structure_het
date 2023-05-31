pkgs <- c("terra", "tidyverse", "sf")
sapply(pkgs, require, character.only=TRUE, quietly=TRUE)

f <- list.files("~/../../Volumes/az_drive/uas_2022/", pattern = "_chm_", full.names = FALSE, recursive = TRUE) |>
  lapply(function(x){str_split(pattern = "/",x)[[1]][1]}) |> unlist() 

resdf <- data.frame(site = NULL, scale = NULL, res = NULL, chm.res = NULL)
for(i in f){
  # ======================
  site <- i
  # =====================================================
  chmres <- list.files("~/../../Volumes/az_drive/uas_2022/", pattern = paste0(i,".*_chm_"), full.names = TRUE, recursive = TRUE) |>
    rast() |> res() |> first()
  
  rl <- list.files("~/../../Volumes/az_drive/wave_chms/", pattern = site, full.names = TRUE)
  rlist <- lapply(rl, function(x) {
    rast(x) |> res() |> first() 
  }) |> unlist()
  out <- data.frame(site = i, scale = 1:9, res = rlist, chm.res = chmres)
  resdf <- bind_rows(resdf, out)
}
write.csv(resdf, "data/chm_resolutions.csv", row.names = FALSE)
