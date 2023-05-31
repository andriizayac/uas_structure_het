pkgs <- c("tidyverse", "terra", "sf", "brms", "ggplot2", "ggridges", "patchwork", "wavethresh")
sapply(pkgs, require, character.only=TRUE, quietly=TRUE)

# === scale resolutions: used for labels - an average from across sites
ress <- read.csv("data/chm_resolutions.csv") |>
  filter(site != "DuncanSaddle") |>
  dplyr::select(scale, res) |>
  group_by(scale) |> summarize_all(mean) |> pull(res) |> round(2)
# ===

# === colour palette
colfun <- colorRampPalette(c("#25858EFF", "#FDE725FF"))
cbcols <- colfun(64)

# === Modeling results

# === model 1.0 - stand-level wildfire effects
m1.0 <- readRDS("../wavelet_models/m1_0_burn_scale_models.rds")
map(m1.0, as.data.frame) |>
  reduce(bind_rows) |>
  mutate(scale = paste0("Scale_", sort(rep(1:9, 4000)) )) |> 
  dplyr::select(b_burnt, scale) |>
  group_by(scale) |>
  summarize(mu = mean(b_burnt), sd = sd(b_burnt), 
            l05 = quantile(b_burnt, 0.05), 
            l95 = quantile(b_burnt, 0.95)) |> ungroup() |>
  ggplot(aes(x = scale, y = mu, group = 1)) + geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = l05, ymax = l95), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#101010") +
  scale_x_discrete(limits = rev, labels = rev(ress)) + 
  labs(x = "Scale, m", 
       y = expression(paste("Fire effect, ", beta))) +
  theme_bw()
ggsave("figures/fig_3a.pdf", width = 5.3, height = 3.3)

 # === model 2.0 - sparse model
m2.0 <- readRDS("../wavelet_models/m2_0.rds")

m2.0 |> as.data.frame() |> 
  dplyr::select(contains("Scale")) |> 
  pivot_longer(cols = everything()) |> 
  group_by(name) |> summarize(mean = mean(value), 
                              lo = quantile(value, probs = 0.025), 
                              hi = quantile(value, probs = 0.975)) |> mutate(res = ress) |> 
  ungroup() |> 
  mutate(col = factor(as.numeric(sign(lo) == sign(hi))) ) |>
  ggplot(aes(x = name, y = mean, group = 1)) + geom_point(size = 2) + 
  geom_line() +
  geom_errorbar(aes(ymin = lo, ymax = hi, width = .1)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#101010") +
  scale_x_discrete(limits = rev, labels = rev(ress)) +
  scale_color_manual(values = c("gray", "black")) +
  labs(x = "Scale, m", y = expression(paste("Effect on recruit abundance, ", beta))) +
  theme_bw()
ggsave("figures/fig_3b.pdf", width = 5.3, height = 3.3)

# === model 2.1 and 2.2 - predictive power
acc2.1 <- read.csv("../wavelet_models/accuracy2_1.csv") |>
  mutate(across(c(mae, rmse), function(x) x/(pi*5^2)))
acc2.1tot <- read.csv("../wavelet_models/accuracy2_tot_het.csv") |> 
  mutate(across(c(mae, rmse), function(x) x/(pi*5^2)))

acc2.1 |> 
  dplyr::select(-k) |>
  group_by(scale) |>
  summarize_all(mean) |> ungroup() |> 
  ggplot(aes(x = scale)) + 
  geom_point(aes(y = r2)) + geom_line(aes(y = r2)) +
  geom_hline(yintercept = mean(acc2.1tot$r2), linewidth = .5) +
  geom_line(aes(y = r2l), linetype = "dotted", linewidth = .5) +
  geom_line(aes(y = r2u), linetype = "dotted", linewidth = .5)  +
  geom_line(aes(y = mae), colour = "#666666", linetype = "dashed", linewidth = 2) +
  geom_hline(yintercept = mean(acc2.1tot$mae), colour = "#666666", linetype = "dashed", linewidth = .5) +
  scale_y_continuous(name = expression(paste(R^2)),
                     sec.axis = sec_axis(~.*1, 
                                         name = expression(paste("Mean absulute error, N", ~m^{-2})) )) +
  scale_x_reverse(name = "Scale, m", breaks = 1:9, labels = ress) + 
  annotate("text", x = c(8.5), y = c(0.37), label = "Total", alpha = .7) +
  theme_bw() +
  theme(axis.title.y.right = element_text(colour = "#666666", face = "bold")) 
ggsave("figures/fig_4a.pdf", width = 5.5, height = 3.5)

# === model 2.2 - predictive power
acc2.2 <- read.csv("../wavelet_models/accuracy2_2.csv") |> 
  mutate(across(c(mae, rmse), function(x) x/(pi*5^2)))
acc2.2tot <- read.csv("../wavelet_models/accuracy2_tot_het.csv") |> 
  mutate(across(c(mae, rmse), function(x) x/(pi*5^2)))

acc2.2 |> 
  dplyr::select(-k) |>
  group_by(scale) |>
  summarize_all(mean) |> ungroup() |> 
  ggplot(aes(x = scale )) + 
  # reference line
  geom_hline(yintercept = mean(acc2.2tot$r2)) +
  geom_point(aes(y = r2)) + geom_line(aes(y = r2)) +
  geom_line(aes(y = r2l), linetype = "dotted", linewidth = .5) +
  geom_line(aes(y = r2u), linetype = "dotted", linewidth = .5)  +
  # secondary axis with MAE
  geom_line(aes(y = mae/1.2), colour = "#666666", linetype = "dashed", linewidth = 2) +
  geom_hline(yintercept = mean(acc2.2tot$mae)/1.2, colour = "#666666", linetype = "dashed", linewidth = .5) +
  scale_y_continuous(name = expression(paste(R^2)),
                     sec.axis = sec_axis(~.*1.2, 
                                         name = expression(paste("Mean absulute error, N", ~m^{-2})) )) +
  scale_x_reverse(name = "Scale, m", breaks = 1:9, labels = ress) + 
  annotate("text", x = c(1.3), y = c(0.41), label = "Total", alpha = .7) +
  theme_bw() +
  theme(axis.title.y.right = element_text(colour = "#666666", face = "bold"))
ggsave("figures/fig_4b.png", width = 5.5, height = 3.5)

# === model 3.0 - Inferential part
# === load data
dft <- read.csv("../wavelet_models/data_demo.csv")
# ===

m3.0 <- readRDS("../wavelet_models/m3_1.rds")
pars <- m3.0 |> as.data.frame() |> dplyr::select(-b_Scale_2sd)

# --- coefficient plot
f4a <- pars |>
  dplyr::select(contains("Scale"),contains("logMat")) |> 
  pivot_longer(cols = everything()) |> 
  ggplot(aes(x = value, y = name, height = ..density..)) + 
  geom_density_ridges(scale = 1, stat = "density", panel_scaling = TRUE, rel_min_height = 0.025, alpha = .3) + 
  labs(y = "", x = expression(paste("Coefficient estimate, ", beta) )) +
  scale_y_discrete(labels = c("logMature", "logMature:Heterogeneity", "Heterogneity")) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#101010") +
  theme_bw()

# --- marginal effect plot stratified by 
# individuals: log(1) and log(25) (density=[1/25, 1] m^2)
het <- seq(min(dft$Scale_1sd), max(dft$Scale_1sd), l = nrow(dft))
yhat.lo <- apply(pars, 1, function(x){ yhat <- (x[1] + x[2]*log(3) + x[3]*het + x[4]*log(1)*het) }) 
yhat.hi <- apply(pars, 1, function(x){ yhat <- (x[1] + x[2]*log(12) + x[3]*het + x[4]*log(12)*het) }) 

df <- data.frame(het = het) |>
  mutate(yhat.lo.mu = apply(yhat.lo, 1, mean),
         yhat.lo.lo = apply(yhat.lo, 1, quantile, probs = 0.05),
         yhat.lo.hi = apply(yhat.lo, 1, quantile, probs = 0.95)) |>
  mutate(yhat.hi.mu = apply(yhat.hi, 1, mean),
         yhat.hi.lo = apply(yhat.hi, 1, quantile, probs = 0.05),
         yhat.hi.hi = apply(yhat.hi, 1, quantile, probs = 0.95))
  
cols <- c('low' = "#482173FF", 'high' = "#51C56AFF")
f4b <- df|> 
  ggplot(aes(x = het)) +
  geom_line(aes(y = yhat.lo.mu, colour = "#482173FF"), size = 1) +
  geom_line(aes(y = yhat.hi.mu, colour = "#51C56AFF"),  size = 1) +
  geom_point(data = dft |> 
               filter(burnt == 1), 
             aes(x = Scale_1sd, y = log(juv)), inherit.aes = FALSE, alpha = .1) +
  geom_ribbon(aes(ymin = yhat.lo.lo, ymax = yhat.lo.hi), fill = "#482173FF", alpha = .33) +
  geom_ribbon(aes(ymin = yhat.hi.lo, ymax = yhat.hi.hi), fill = "#51C56AFF", alpha = .33) +
  geom_line(aes(y = yhat.lo.mu), colour = "#482173FF", linewidth = 1) +
  geom_line(aes(y = yhat.hi.mu), colour = "#51C56AFF",  linewidth = 1) +
  scale_colour_manual(name = "N, mature", values =  cols) +
  labs(x = expression(paste("log[ Heterogeneity, ", Sigma~bold(D)[1] ^2,"(",bold(x),",",bold(y),") ]" )), 
       y = "log[ N, juvenile ]") +
  theme_bw()
f4a + f4b
ggsave("figures/fig_4_1.png", width = 10, height = 4)

# === rasters: canopy height and heterogeneity across scales 
library(wavethresh)
# purpose: demonstrate an example
# excerpt from scale_decompose.R
# InitialPointBurnt
f <- list.files("~/../../Volumes/az_drive/uas_2022/", pattern = "_chm_", full.names = TRUE, recursive = TRUE)
rois <- list.files("~/../../Volumes/az_drive/FieldData//", pattern = "_5m_.*tif$", full.names = TRUE, recursive = TRUE)
# ====================
site <- "InitialPointBurn"
# ====================
f.plant <- st_read("~/../../Volumes/az_drive/FieldData/FieldData_2022/plant_data.geojson") |> 
  filter(ID == !!paste0(site, "_1782"))


file <- f[grepl(site, f)]
roi <- rast(rois[grepl(site, rois)])
chm <- rast(file) |> crop(ext(st_buffer(f.plant, 5))) 

n <- 512
e <- as.matrix(chm, wide = TRUE)[0:n, 0:n]
e.rast <- chm[0:n, 0:n, drop = FALSE]
# --- wavelet decomposition -----------------
w.filter.number = 1; w.family = "DaubExPhase"
# -------------------------------------------
wn <- imwd(e, filter.number = w.filter.number, 
           family = w.family, 
           type = "wavelet", bc = "symmetric")
wnl <- list()
wnh <- list()
for(i in (wn$nlevels-1):0 ) {
  
  wnl[[i+1]] <- imwr( nullevels(wn, levelstonull = c(wn$nlevels-1):i) ) |> rast()
  
  HH <- wn[[paste0("w", i,"L1")]] # horizontal
  hh <- matrix(HH, nc = sqrt(length(HH)))

  VV <- wn[[paste0("w", i,"L2")]] # vertical
  vv <- matrix(VV, nc = sqrt(length(VV)))

  DD <- wn[[paste0("w", i,"L3")]] # diagonal
  dd <- matrix(DD, nc = sqrt(length(DD)))
  
  tmp <- rast(hh^2 + vv^2 + dd^2)
  ext(tmp) <- ext(e.rast)
  crs(tmp) <- "epsg:32611"
  
  wnh[[i+1]] <- tmp
}

par(mfrow = c(4,2), mai = c(0,0,0,0))
for(i in 2:length(wnh)){ plot(wnh[[i]], col = cbcols)}

png('figures/archive/fig3b4.png', height = 320, width = 390, units = "px")
# export chms: fig2b[1:3], at aggregate levels (4, 16, 128)
# plot(aggregate(e.rast, fact = 128), col = cbcols, axes = FALSE, mar = c(1.1, .5, 1.1, 4.75), plg = list(cex = 1.5) )
# export transforms: fig2b[4:6] with whn[[3, 6, 8]]
plot(wnh[[8]], col = cbcols, axes = FALSE, mar = c(1.1, .5, 1.1, 4.75), plg = list(cex = 1.5))
dev.off()




# imwr( nullevels(wn, levelstonull = c(2) )) |> rast() |> plot()

# === extract a single plot as inlet
cellid <- 35
cells <- st_read("data/cells.geojson") |> 
  filter(site == "Cold", cid == cellid) |> 
  st_buffer(5.1) |> st_transform(32611) 
plants <- st_read("../plant_data.geojson") |>
  filter(site == "Cold", CID == cellid, Class == "Shrub", Ht_gt_25 == 0) |> 
  st_transform(32611) |> 
  st_intersection(cells)

chm <- rast("~/../../Volumes/az_drive/uas_2022/Cold/Cold_20220610_chm_wgs84utm11n.tif") |> 
  crop(st_as_sfc(st_bbox(cells)) ) |> mask(cells)
  

# pdf("~/Downloads/temp_fig.pdf")
par(mai = c(0,0,0,0))
plot(chm, col = cbcols, legend = FALSE, axes = FALSE)
plot(plants$geometry, pch = 19, add = TRUE)
plot(cells$geometry, border = "#7030A0", lwd = 5, add = TRUE)
# plot(st_buffer(cells, .01)$geometry, add = TRUE, lwd = 2, border = "white")
dev.off()
# ===
