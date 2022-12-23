pkgs <- c("terra", "tidyverse", "wavethresh", "sf")
sapply(pkgs, require, character.only=TRUE, quietly=TRUE)

# --- wavelet decomposition -----------------
w.filter.number = 1; w.family = "DaubExPhase"
# -------------------------------------------


f <- list.files("~/../../Volumes/az_drive/temp/CorralsTrail/", pattern = "_chm_", full.names = TRUE)#list.files("~/Downloads/tmp", pattern = "_chm_", full.names = TRUE)
roi <- rast("~/Desktop/UAV_demography/FieldData/2022_campaign_field_data/cells_data/CorrlasTrail_5m_raster.tif")
chm <- rast(f[2]) |> crop(ext(roi) + 10)
chm <- classify(chm, rcl = matrix(c(-Inf, 0, 0, 
                                    3, Inf, 0), nc = 3, byrow = TRUE))

# === Decompose the matrix into tiles stored in a list
cmat <- as.matrix(chm, wide = TRUE)

n <- 256
xtiles <- floor(ncol(chm)/n)
ytiles <- floor(nrow(chm)/n)

xN <- n * xtiles/2 
yN <- n * ytiles/2
xmp <- floor(ncol(chm)/2)
ymp <- floor(nrow(chm)/2)

idl <- xmp - xN; idr <- xmp + (xN-1) # left and right boundaries
idb <- ymp + (yN); idt <- ymp - (yN-1) # top and bottom boundaries

e <- cmat[idt:idb, idl:idr]
e_rast <- chm[idt:idb, idl:idr, drop = FALSE]

# raster index starts in NW corner
# list - level 1 - rows 
#   |___ - level 2 - columns
# ie columns are stored within rows
lout <- list()

xidx <- 1; 
for(i in 1:xtiles) {
  ltmp <- list()
  yidx <- 1;
  for(j in 1:ytiles) {
    tmp <- e[yidx:c(j*n), xidx:c(i*n)]
    yidx <- yidx + n
    ltmp[[j]] <- tmp
    # out <- paste0("~/Downloads/tmp_rast/tile_", i, ".csv")
    # write.table(tmp, out, row.names = FALSE, col.names = FALSE, sep=',')
  }
  lout[[i]] <- ltmp
  xidx <- xidx + n
}

# === Bring the tiles back together into a larger image
out <- NULL;
for(i in 1:xtiles){
  tmp <- NULL;
  for(j in 1:ytiles){
    tmp = rbind(tmp, lout[[i]][[j]])
  }
  out <- cbind(out, tmp)
}
plot(rast(e))
plot(rast(out))
# ==== end reassembling



# === apply wavelet transformation and store the output
total <- 1
# list storing wavelet coefs
ldiff <- list()
# list storing full transform output
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
        # hh <- matrix(approx(HH, n = n^2)$y, nc = n) # nc = sqrt(length(HH)))
      VV <- wn[[paste0("w", k,"L2")]] # vertical
        vv <- matrix(VV, nc = sqrt(length(VV)))
        # vv <- matrix(approx(VV, n = n^2)$y, nc = n) #nc = sqrt(length(VV)))
      DD <- wn[[paste0("w", k,"L3")]] # diagonal
        dd <- matrix(DD, nc = sqrt(length(DD)))
        # dd <- matrix(approx(DD, n = n^2)$y, nc = n) #nc = sqrt(length(DD)))
      SS <- wn[[paste0("w", k,"L4")]] # averages
        ss <- matrix(SS, nc = sqrt(length(SS)))
        # ss <- matrix(approx(SS, n = n^2)$y, nc = n) #nc = sqrt(length(SS)))
      if(total){
        D[[k+1]] <- hh^2 + vv^2 + dd^2 # + abs(ss)
      } else {
        D[[k+1]] <- hh + vv + dd #+ ss
      }
      
    }
    ldiffc[[j]] <- D
  }
  lro[[i]] <- lcol # imwd.objects
  ldiff[[i]] <- ldiffc # difference coefs interpolate for nlevel-1:1
}

# === reassemble wavelet coefs back into #nlevels CHM images
rlist <- list()

par(mfrow = c(2,3))
# plot(rast(e), legend = FALSE, col=topo.colors(64))
for(k in 1:(log2(n)-1 - 1)){
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

# === relationship between cover and structural heterogeneity
em <- rast(e) |> classify(rcl = matrix(c(-Inf, .25, 0,
                                         .25, Inf, 1), nc = 3, byrow = TRUE)) |>
  as.matrix(wide = TRUE)

r <- rast(ncols=4, nrows=4, xmin = n, xmax = n*5, ymin = n, ymax = n*5)
values(r) <- 1:ncell(r)
plot(rast(em))
plot(r, add = TRUE)

rpts <- as.points(r) |> st_as_sf() |>
  mutate(x = st_coordinates(geometry)[,1], 
         y = st_coordinates(geometry)[,2]) |>
  rename(id = lyr.1) |>
  st_buffer(50) |>
  dplyr::select(-geometry)


plot(rast(em))
plot(rpts$geometry, add = TRUE)

plot(rlist[[4]])


points(rpts$x/vec[4], rpts$y/vec[4], pch = 19)

# --- calculate response - cover in 1m^2, and predictors at 7 scales
em |> rast() |>
  terra::extract(y = rpts, touches = TRUE) |>
  group_by(ID) |>
  summarize(mean = mean(lyr.1)) |>
  ungroup() -> dfy

rdf <- rpts |>
  as.data.frame() |> 
  dplyr::select(-geometry)
vec <- rev(2^(2:8))
xout <- list()

for(k in 1:(log2(n)-1)) {
  rdf |>
    mutate(x = x/vec[k], y = y/vec[k]) |>
    st_as_sf(coords = c("x", "y")) |>
    st_buffer(50/vec[k])-> tmp_in # within 1m^2 only
  rdf |>
    mutate(x = x/vec[k], y = y/vec[k]) |>
    st_as_sf(coords = c("x", "y")) |> 
    st_buffer(n/vec[k])-> tmp_out # within 1m^2 + neighborhood
  
  
  rlist[[k]] |>
    terra::extract(y = tmp_out, touches = TRUE) |>
    group_by(ID) |>
    summarize(sum = sum(lyr.1)) |>
    ungroup() |> dplyr::select(-ID) |> 
    mutate(var = sum - rlist[[k]] |>
             terra::extract(y = tmp_in, touches = TRUE) |>
             group_by(ID) |>
             summarize(sum = sum(lyr.1)) |>
             ungroup() |> dplyr::select(-ID) |> pull(sum)) -> dfx_tmp
  
  xout[[k]] <- dfx_tmp[,"var"]
}

dfx <- do.call(cbind, xout) 
names(dfx) <- paste0("scale_", 2:8)

# --- combine outputs
df <- bind_cols(dfy, dfx)

cs <- cor(dfy$mean, dfx)

data.frame(corr = t(cs)) |> 
  mutate(scale = 2:8) |>
  pivot_longer(corr) |>
  ggplot() +
  geom_line(aes(x = scale, y = value), lwd = 2) +
  labs(x = "Scale", y = "Correlation") +
  theme_bw()


# === Field data

# --- project transformed rasters
for(i in 1:length(rlist)){
  tmp <- rlist[[i]]
  ext(tmp) <- ext(e_rast)
  rlist[[i]] <- tmp
}
# ---
cells <- st_read("~/Desktop/UAV_demography/FieldData/SurveyPts_shapes_2022/CorralsTrail/CorralsTrail_pts_wgs84.shp") |> 
  st_transform(32611)
fd <- read.csv("~/Desktop/UAV_demography/FieldData/2022_campaign_field_data/plant_data/CorralsTrail_2022_plants.csv") |> 
  filter(!is.na(ID))
ppts <- read.csv("~/Desktop/UAV_demography/FieldData/2022_campaign_field_data/gps/CorralsTrail/CorralsTrail_20220516_20_plants.csv") 
fd |>
  filter(ID != 0) |>
  left_join(ppts, by = c("ID" = "label")) |>
  st_as_sf(coords = c("long", "lat"), crs = 4326) |>
  st_transform(32611) -> plant_df
# ---
burnt <- st_read("~/Desktop/UAV_demography/FieldData/2022_campaign_field_data/burn_lines/CorralsTrail_burn_line_wgs84.geojson") |>
  st_transform(32611) |> filter(burnt == 1)

plant_df |>
  as.data.frame() |> 
  dplyr::select(CID, Class, Species, Ht_gt_25) |>
  group_by(CID) |>
  summarize(pHt25 = sum(Ht_gt_25)/n(), N = n(), 
            Shrubs = sum(Class == "Shrub"), Grasses = sum(Class == "Grass"), 
            ARTR = sum(Species == "ARTR"), PUTR = sum(Species == "PUTR")) -> df


# = scale-dependent heterogeneity
sout <- matrix(NA, nc = length(rlist), nr = nrow(cells))
names()
for(i in 1:ncol(sout)){
  tmp <- rlist[[i]] |>
    terra::extract(y = st_buffer(cells, 5), touches = TRUE, ID = FALSE) |>
    group_by(ID) |>
    summarize(sum = sum(lyr.1)) |>
    ungroup() |> pull(sum)
  
  sout[,i] <- tmp
}
  
cells |> 
  mutate(burnt = lengths(st_intersects(geometry, burnt, sparse = TRUE))) |>
  bind_cols(sout) |> 
  as.data.frame() |> 
  dplyr::select(-geometry) |> 
  right_join(df, by = c("cid" = "CID")) -> out

plot(rast(cor(out |> dplyr::select(-cid))))

out |>
  mutate(burnt = as.factor(burnt)) |> 
  pivot_longer(cols = starts_with("...")) |>
  mutate(scale = gsub("...", "Scale_", name)) |>
  dplyr::select(-name) |> rename(heterogen = value) |>
  mutate(juv = N*pHt25) |>
  ggplot(aes(y = juv, x = log(heterogen), colour = burnt)) + 
  geom_point(size = 3) +
  facet_wrap(~scale) + 
  theme_bw()

f1 <- glm(ARTR ~ 1 + (log(...4)|burnt), family = poisson, data = out)
summary(f1)
as.data.frame(VarCorr(f1))
