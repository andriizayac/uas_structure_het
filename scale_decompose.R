pkgs <- c("terra", "tidyverse", "wavethresh", "sf", "exactextractr")
sapply(pkgs, require, character.only=TRUE, quietly=TRUE)

# 1. CHM wavelet decomposition
# 2. Integrate raster heterogeneity with cell+plant data

# === 1. CHM wavelet decomposition
# --- wavelet decomposition -----------------
w.filter.number = 1; w.family = "DaubExPhase"
# -------------------------------------------

# ===
sites <- unique(st_read("data/cells.geojson") |> pull(site))

# === find CHM's and ROI's for each site
f <- list.files("~/../../Volumes/az_drive/uas_2022/", pattern = "_chm_", full.names = TRUE, recursive = TRUE)
rois <- list.files("~/../../Volumes/az_drive/FieldData//", pattern = "_5m_.*tif$", full.names = TRUE, recursive = TRUE)

# === from here until 2 could be iterated over sites or processed manually one site at a time
# ====================
site <- sites[1]
# ====================

file <- f[grepl(site, f)]
roi <- rast(rois[grepl(site, rois)])
chm <- rast(file) |> crop(ext(roi) + 10)
# chm <- classify(chm, rcl = matrix(c(-Inf, 0, 0, 
#                                     3, Inf, 0), nc = 3, byrow = TRUE))
# chm <- subst(chm, NA, 0)
# === Decompose the raster/matrix into tiles stored in a list
cmat <- as.matrix(chm, wide = TRUE)

n <- 512 # the size of each tile: n x n
xtiles <- floor(ncol(chm)/n)
ytiles <- floor(nrow(chm)/n)

xN <- n * xtiles/2 
yN <- n * ytiles/2
xmp <- floor(ncol(chm)/2)
ymp <- floor(nrow(chm)/2)

idl <- xmp - xN; idr <- xmp + (xN-1) # left and right boundaries
idb <- ymp + (yN); idt <- ymp - (yN-1) # bottom and top boundaries

e <- as.matrix(chm, wide = TRUE)[idt:idb, idl:idr]
e.rast <- chm[idt:idb, idl:idr, drop = FALSE]

# === Break down the CHM in n-size tiles and store them in a nested list
# raster index starts in NW corner
# list[[level1]][[level2]]  
#    |_____level 1 - stores rows 
#    |_____level 2 - stores columns
# ie columns are stored within rows
lout <- list()

xidx <- 1 # column index
for(i in 1:xtiles) {
  ltmp <- list()
  yidx <- 1; # row index
  for(j in 1:ytiles) {
    tmp <- e[yidx:c(j*n), xidx:c(i*n)]
    yidx <- yidx + n
    ltmp[[j]] <- tmp
    # === optionally, can store each tile as a csv file
    # out <- paste0("~/Downloads/tmp_rast/tile_", i, ".csv")
    # write.table(tmp, out, row.names = FALSE, col.names = FALSE, sep=',')
  }
  lout[[i]] <- ltmp
  xidx <- xidx + n
}

# === The following chunk can be used to verify the indxing was correct by re-assembling the tiles
# === Bring the tiles back together into a larger image: 
# output should be same as the cropped CHM
# --- ensures the loops and indexes are correct
# out <- NULL;
# for(i in 1:xtiles){
#   tmp <- NULL;
#   for(j in 1:ytiles){
#     tmp = rbind(tmp, lout[[i]][[j]])
#   }
#   out <- cbind(out, tmp)
# }
# plot(e.rast)
# plot(rast(out))
# ==== end reassembling

# === Apply wavelet transformation and store the output in a nested list
# --- Two options:
# (1) Total = 1: quantifies the total amount of energy per scale (ie x^2 + y^2) - used in the paper
# (2) Total = 0: adds raw difference coefficients (ie x + y)
total <- 1 
# list storing wavelet coefs - quantifies spatial variation
ldiff <- list()
# list storing full transform output - optionally, for data records
lro <- list()

for(i in 1:length(lout)){ # rows
  lcol <- list()
  ldiffc <- list()
  for(j in 1:length(lout[[i]])){ # columns
    dat <- lout[[i]][[j]]
    
    wn <- imwd(dat, filter.number = w.filter.number, 
               family = w.family, 
               type = "wavelet", bc = "symmetric")
    
    lcol[[j]] <- wn
    
    D <- list()
    for(k in (wn$nlevels-1):0) { # decomposition levels
      HH <- wn[[paste0("w", k,"L1")]] # horizontal
      hh <- matrix(HH, nc = sqrt(length(HH)))

      VV <- wn[[paste0("w", k,"L2")]] # vertical
      vv <- matrix(VV, nc = sqrt(length(VV)))
      
      DD <- wn[[paste0("w", k,"L3")]] # diagonal
      dd <- matrix(DD, nc = sqrt(length(DD)))
      
      SS <- wn[[paste0("w", k,"L4")]] # averages
      ss <- matrix(SS, nc = sqrt(length(SS)))
      
      if(total){
        D[[k+1]] <- hh^2 + vv^2 + dd^2 
      } else {
        D[[k+1]] <- hh + vv + dd 
      }
      
    }
    ldiffc[[j]] <- D
  }
  lro[[i]] <- lcol # imwd.objects
  ldiff[[i]] <- ldiffc # difference coefs for nlevel-1:8
}
# === reassemble wavelet coefs back into #nlevels CHM images
rlist <- list()

par(mfrow = c(2,4))
plot(rast(e), legend = FALSE, col=topo.colors(64))
for(k in 1:(log2(n))){
  level <- k
  
  out <- NULL;
  for(i in 1:xtiles){
    tmp <- NULL;
    for(j in 1:ytiles){
      tmp = rbind(tmp, ldiff[[i]][[j]][[level]]) # lout[[i]][[j]])# 
    }
    out <- cbind(out, tmp)
  }
  
  plot(rast(out), col=topo.colors(64), legend = FALSE)
  rlist[[k]] <- rast(out)
}

# === end assembling difference coefs tiles

# --- project transformed rasters
for(i in 1:length(rlist)){
  tmp <- rlist[[i]]
  ext(tmp) <- ext(e.rast)
  crs(tmp) <- "epsg:32611"
  rlist[[i]] <- tmp
}

# === export projected, transformed rasters for futher analysis
sapply(1:length(rlist), function(x){
  writeRaster(rlist[[x]], paste0("~/../../Volumes/az_drive/wave_chms/", site, "_scale_wtransform_L", x, "_", n, ".tif"))
})
# ======================



# 2. Integrate raster hetergen with cell+plant data
cells <- st_read("data/cells.geojson")|> 
  st_transform(32611)
plants <- st_read("../plant_data.geojson") |>
  st_transform(32611)
blines <- st_read("data/burn_lines.geojson") |>
  st_transform(32611)

# ======================
site <- sites[1]
# =====================================================
rl <- list.files("~/../../Volumes/az_drive/wave_chms/", pattern = site, full.names = TRUE)
rlist <- lapply(rl, rast)
# =====================================================
plants0 <- filter(plants, site == !!site) 
blines0 <- filter(blines, site == !!site) 
cells0 <- filter(cells, site == !!site) |> 
  mutate(burnt = lengths(st_intersects(geometry, filter(blines0, burn == 1), sparse = TRUE)))

# === scale-dependent heterogeneity

for(i in 1:length(rlist)){
  var <- paste0("Scale_", i)
  cells0 |>
    mutate(!!var := exact_extract(rlist[[i]], st_buffer(cells0, 5), fun = "sum")) -> cells0
}

plants0 |>
  filter(Class == "Shrub") |>
  as.data.frame() |>
  dplyr::select(CID, Ht_gt_25) |>
  group_by(CID) |>
  summarize(N = n(), juv = sum(Ht_gt_25), mat = sum(Ht_gt_25 == 0)) -> plantdf
  
cells0 |>
  left_join(plantdf, by = c("cid" = "CID")) |>
  mutate(across(.cols = c(N, juv, mat), ~replace_na(.x, 0))) |>
  # mutate(N = NA, juv = NA, mat = NA) |> # for DuncanSaddle
  # dplyr::select(-path, -layer) |>
  pivot_longer(cols = starts_with("Scale")) |> 
  rename(scale = name, heterogen = value) |> 
  as.data.frame() |> dplyr::select(-geometry, -path, -layer) -> celldf

write.csv(celldf,
          paste0("~/../../Volumes/az_drive/wave_chms/cell_heterogen/", site, "_cell_heterogen.csv"),
          row.names = FALSE)
# === end
