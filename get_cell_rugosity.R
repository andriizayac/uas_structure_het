pkgs <- c("terra", "tidyverse", "sf")
sapply(pkgs, require, character.only=TRUE, quietly=TRUE)

# === This script extracts canopy top rugosity from the CHMs
# refs: Atkins et al. 2018 & https://www.neonscience.org/resources/learning-hub/tutorials/structural-diversity-discrete-return
cells <- st_read("data/cells.geojson") |> 
  filter(site != "DuncanSaddle") |> st_buffer(5) |>
  mutate(rugosity = NA)

# - crop templates, 5m resolution 
l0 <- list.files("~/../../Volumes/az_drive/uas_2022_5m/", full.names = TRUE)
rnull <- lapply(l0[-grep("DuncanSaddle", l0)], function(x) {
  tmp <- rast(x) 
  set.crs(tmp, "epsg:32611")
  return(tmp)
})
# - rasters: chms
lchm <- list.files("~/../../Volumes/az_drive/uas_2022/", pattern = "_chm_", full.names = TRUE, recursive = TRUE) 
chms <- lapply(1:10, function(x) {
  rast(lchm[-grep("DuncanSaddle", lchm)][x]) |> crop(ext(rnull[[x]]) + 5) })

for(i in 1:length(unique(cells$site))) {
  j <- unique(cells$site)[i]
  idx <- which(cells$site == j)
  
  rug <- extract(chms[[i]], cells[idx,], fun = sd, method = "simple", weights = FALSE, exact = FALSE)
  cells$rugosity[idx] <- rug[,2]
}

cells |> 
  as.data.frame() |> select(cid, site, rugosity) |> 
  write.csv("data/cell_rugosity.csv", row.names = FALSE)

