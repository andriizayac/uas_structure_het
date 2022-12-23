pkgs <- c("tidyverse", "terra", "sf", "brms", "ggplot2", "ggridges", "patchwork", "wavethresh")
sapply(pkgs, require, character.only=TRUE, quietly=TRUE)

# === scale resolutions: used for labels - an average from across sites
ress <- read.csv("data/chm_resolutions.csv") |>
  filter(site != "DuncanSaddle") |>
  dplyr::select(scale, res) |>
  group_by(scale) |> summarize_all(mean) |> pull(res) |> round(2)
# ===

# === load data
dft <- read.csv("data/data_demo.csv")
# ===


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
  scale_x_discrete(limits = rev, labels = rev(ress)) + 
  labs(x = "Scale resolution, m", 
       y = expression(paste("Fire effect, ", beta))) +
  theme_bw()
ggsave("figures/fig_1.png", width = 6, height = 4)

# === model 2.0 - sparse model
m2.0 <- readRDS("../wavelet_models/m2_0.rds")

m2.0 |> as.data.frame() |> 
  dplyr::select(contains("Scale")) |> 
  pivot_longer(cols = everything()) |> 
  group_by(name) |> summarize(mean = mean(value), 
                              lo = quantile(value, probs = 0.025), 
                              hi = quantile(value, probs = 0.975)) |> mutate(res = ress) |> 
  ggplot(aes(y = name)) + geom_point(aes(x = mean), size = 2) + 
  geom_errorbar(aes(xmin = lo, xmax = hi, width = .1)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#101010") +
  scale_y_discrete(labels = ress ) +
  labs(y = "Scale resolution, m", x = expression(paste("Coefficient estimate, ", beta))) +
  theme_bw()
# ggsave("figures/fig_2.png", width = 4, height = 6)

# === model 2.1 and 2.2 - predictive power
acc2.1 <- read.csv("../wavelet_models/accuracy2_1.csv")

acc2.1 |> 
  dplyr::select(-k) |>
  group_by(scale) |>
  summarize_all(mean) |> ungroup() |> 
  ggplot(aes(x = scale)) + 
  geom_point(aes(y = r2)) + geom_line(aes(y = r2)) +
  geom_line(aes(y = r2l), linetype = "dotted", size = .5) +
  geom_line(aes(y = r2u), linetype = "dotted", size = .5)  +
  geom_line(aes(y = mae/100), colour = "#666666", linetype = "dashed", size = 2) +
  scale_y_continuous(name = expression(paste(R^2)),
                     sec.axis = sec_axis(~.*100, name = "Mean absulute error, N")) +
  scale_x_reverse(name = "Scale resolution, m", breaks = 1:9, labels = ress) + 
  theme_bw() +
  theme(axis.title.y.right = element_text(colour = "#666666", face = "bold"))
# ggsave("figures/fig_3a.png", width = 6, height = 4)

# === model 2.2 - predictive power
acc2.2 <- read.csv("../wavelet_models/accuracy2_2.csv")

acc2.2 |> 
  dplyr::select(-k) |>
  group_by(scale) |>
  summarize_all(mean) |> ungroup() |> 
  ggplot(aes(x = scale)) + 
  geom_point(aes(y = r2)) + geom_line(aes(y = r2)) +
  geom_line(aes(y = r2l), linetype = "dotted", size = .5) +
  geom_line(aes(y = r2u), linetype = "dotted", size = .5)  +
  geom_line(aes(y = mae/100), colour = "#666666", linetype = "dashed", size = 2) +
  scale_y_continuous(name = expression(paste(R^2)),
                     sec.axis = sec_axis(~.*100, name = "Mean absulute error, N")) +
  scale_x_reverse(name = "Scale resolution, m", breaks = 1:9, labels = ress) + 
  theme_bw() +
  theme(axis.title.y.right = element_text(colour = "#666666", face = "bold"))
# ggsave("figures/fig_3b.png", width = 6, height = 4)

# === model 3.0 - inference
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
  geom_line(aes(y = yhat.lo.mu), colour = "#482173FF", size = 1) +
  geom_line(aes(y = yhat.hi.mu), colour = "#51C56AFF",  size = 1) +
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
for(i in 2:length(wnh)){ plot(wnh[[i]])}

png('figures/archive/fig2b6.png', height = 320, width = 360, units = "px")
plot(wnh[[3]], axes = FALSE, mar = c(1.1, 1.1, 1.1, 4.1))

dev.off()


imwr( nullevels(wn, levelstonull = c(2) )) |> rast() |> plot()

