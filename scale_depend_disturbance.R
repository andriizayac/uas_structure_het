pkgs <- c("terra", "tidyverse", "sf")
sapply(pkgs, require, character.only=TRUE, quietly=TRUE)

# === Quantify scale-dependent variation between burnt/unburnt sites

# === Load datasets
# ============
ress <- read.csv("data/chm_resolutions.csv") |>
  filter(site != "DuncanSaddle") |>
  dplyr::select(scale, res) |>
  group_by(scale) |> summarize_all(mean) |> pull(res) |> round(2)
# ===
# ===
# ===========
cells <- st_read("~/../../Volumes/az_drive/FieldData/FieldData_2022/cells.geojson")
plants <- st_read("~/../../Volumes/az_drive/FieldData/FieldData_2022/plant_data.geojson")
blines <- st_read("~/../../Volumes/az_drive/FieldData/FieldData_2022/burn_lines.geojson") |>
  st_transform(32611) |> filter(burn == 1)
# ===========
# ===
dat <- list.files("~/../../Volumes/az_drive/wave_chms/cell_heterogen/", pattern = ".csv", 
                  full.names = TRUE)
dflist <- lapply(dat, read.csv)
df <- do.call(rbind, dflist) |>
  left_join(blines |> as.data.frame() |>
              select(site, year, elevation), c("site" = "site")) |>
  mutate(year = as.numeric(year)) |>
  filter(site != "LowerDryCreek")
# === 

df |>
  filter(is.na(heterogen) == FALSE) |>
  ggplot(aes(x = as.factor(scale), y = (heterogen), fill = as.factor(burnt))) + 
  geom_boxplot() + 
  labs(x = "Scale", y = "Heterogeneity", fill = "Burn status") +
  facet_wrap(.~ site, scales = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df |>
  # group_by(site) |>
  # mutate(heterogen = heterogen/sum(heterogen)) |> ungroup() |>
  dplyr::select(cid, site, burnt, scale, year, elevation, heterogen) |>
  pivot_wider(names_from = scale, values_from = heterogen) |>
  rowwise() |>
  mutate(tot_heterogen = sum(across(starts_with("Scale_")))) -> het

# === difference model for each scale
k <- 9
mlist <- list()

for(i in 1:k) {
  resp <- paste0("Scale_", i)
  f <- formula(paste0(resp, "~ (1|site) + burnt"))
  
  mlist[[i]] <- brms::brm(f, data = het, cores = 4)
}
# saveRDS(mlist, "../wavelet_models/m1_0_burn_scale_models.rds")
mlist <- readRDS("../wavelet_models/m1_0_burn_scale_models.rds")

map(mlist, as.data.frame) |>
  reduce(bind_rows) |> 
  mutate(scale = paste0("Scale_", sort(rep(1:9, 4000)) )) |>
  select(b_burnt, scale) |>
  group_by(scale) |>
  summarize(mu = mean(b_burnt), sd = sd(b_burnt), 
            l05 = quantile(b_burnt, 0.05), 
            l95 = quantile(b_burnt, 0.95)) |> ungroup() |>
  ggplot(aes(x = scale, y = mu, group = 1)) + geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = l05, ymax = l95), width = 0.1) +
  scale_x_discrete(limits = rev, labels = rev(ress) ) + 
  labs(x = "Scale (m)", 
       y = expression(paste("Fire effect ( ", psi,"(x*,y*) )"))) +
  theme_bw()


# --- 
library(brms)
dat <- dif |>
  mutate(across(c(year, elevation), 
                function(x) {(x-mean(x))/(2*sd(x))} ))
mb <- brm(diff ~ (1|scale) + year + elevation,
          data = dat, family = Gamma(link = "log"))
mcmc_plot(mb, variable = c("b_year", "b_elevation")) + 
  labs(x = "Effect size") +
  scale_y_discrete(labels = c("Year", "Elevation")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12))


# ==========================================================
# === PCA and NMDS
# ==========================================================
library(ggord)
library(vegan)
df |>
  dplyr::select(cid, site, burnt, scale, heterogen) |>
  #mutate(heterogen = log(heterogen)) |> 
  # group_by(site, burnt, scale) |>
  # summarize(heterogen = sum(heterogen)) |>
  ungroup() |> 
  pivot_wider(names_from = scale, values_from = heterogen) |>
  dplyr::select(starts_with("Scale")) -> out

# ---
pca <- princomp(~ ., data = out)
biplot(pca, choices = c(1,2))
ggord(pca, as.factor(out$burnt))
# ---
nms <- metaMDS(out, k = 2, trymax = 100)
ordiplot(nms, type = "n")
ordihull(nms, groups = paste0("Burnt_",out$burnt), draw="polygon", 
         col="grey90", label = TRUE)
orditorp(nms, display = "species", col = "red", air = 0.01)
#orditorp(nms, display = "sites", cex = 1.25, air = 0.01)

# ---
# === export ggplot
ggsave("figures/Scales_fig_4.pdf", width = 6, height = 3)
   





# ============ Archive ==============================
df |>
  dplyr::select(site, burnt, scale, year, elevation, heterogen, N) |>
  group_by(site) |>
  mutate(heterogen = heterogen/sum(heterogen)) |> # (heterogen-mean(heterogen))/(2*sd(heterogen)) ) |>
  ungroup() |>
  group_by(site, burnt, scale) |>
  summarize(het = mean(heterogen),
            sd_het = sd(heterogen),
            year = mean(year), elevation = mean(elevation),
            N = mean(N, na.rm = TRUE)) |>
  ungroup() -> het

set.seed(123)
het |>
  mutate(scale = paste0(scale, ": ", ress, "m")) |>
  ggplot(aes(x = as.factor(burnt), y = (het), colour = scale, group = scale)) +
  geom_point(size = 1, position=position_dodge(0.5)) + geom_line() +
  geom_errorbar(aes(ymin = het-sd_het,
                    ymax = het+sd_het), width = .5,
                position = position_dodge(0.5)) +
  scale_color_viridis_d() +
  labs(x = "Burn status", 
       y = expression(paste("Structural heterogeneity, ", Sigma~bold(D)[s]^2,"(", bold(x),",",bold(y), ")"))) +
  theme_bw() +
  facet_wrap(.~ site, scales = "free_y", ncol = 2)
ggsave("figures/fig_S2.png", width = 6, height = 6)
# 
# 
# het |>
#   dplyr::select(-sd_het) |>
#   arrange(site, desc(burnt), scale) |>
#   group_by(site, scale) |>
#   mutate(diff = diff(het)) |>
#   filter(burnt == 1) |> ungroup() -> dif
# dif |>
#   # mutate(scale = as.numeric(rep(round(ress[1,], 2), 4)) ) |>
#   ggplot(aes(y = diff, x = log(N), colour = scale, group = scale)) +
#   geom_jitter(size = 2.5) + geom_line() +
#   #scale_x_discrete(labels = round(ress[1,], 2)) +
#   #labs(x = "Spatial scale (m)", y = str_wrap("Difference b/w burnt and unburnt", 30)) +
#   theme_bw()

# # ---
# set.seed(123)
# dif |>
#   # filter(scale == "Scale_1") |>
#   # pivot_wider(names_from = scale, values_from = heterogen) |>
#   # group_by(site, burnt) |>
#   # summarize(diff = mean(diff), year = mean(year), elevation = mean(elevation)) |>
#   # ungroup() |> 
#   ggplot(aes(x = as.factor(year), y = (diff), fill = site)) + 
#   geom_boxplot() +
#   geom_jitter(alpha = .5, width = 0.1) + 
#   labs(x = "Year burnt", y = str_wrap("Difference b/w burnt and unburnt", 30)) +
#   theme_bw() 
