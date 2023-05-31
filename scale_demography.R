pkgs <- c("tidyverse", "terra", "sf")
sapply(pkgs, require, character.only=TRUE, quietly=TRUE)

# === Demographic processes and recruitment
# ===========
cells <- st_read("data/cells.geojson")
plants <- st_read("../plant_data.geojson")
blines <- st_read("data/burn_lines.geojson") |>
  st_transform(32611) |> filter(burn == 1)
# ===========
rugosity <- read.csv('data/cell_rugosity.csv') |> pull(rugosity)
# ===
dat <- list.files("data/cell_heterogen/", pattern = ".csv", full.names = TRUE)
dflist <- lapply(dat, read.csv)
df <- do.call(bind_rows, dflist) |>
  left_join(blines |> as.data.frame() |>
              dplyr::select(site, year, elevation), c("site" = "site")) |>
  mutate(year = as.numeric(year), heterogen = log(heterogen)) |> 
  filter(site != "DuncanSaddle")
  # group_by(site, scale) |>
  # mutate(heterogen = scale(heterogen))|>
  # ungroup()
# === 

# --- wide format 
df |>
  mutate(site = as.factor(site)) |>
  pivot_wider(names_from = scale, values_from = heterogen) |>
  rowwise() |>
  mutate(tot_heterogen = sum(across(starts_with("Scale_"))) ) |>
  ungroup() -> dft # for glm models

# ---

# --- various figures
df |>
  group_by(site, burnt, cid) |>
  summarize(N = mean(N, na.rm = TRUE), juv = mean(juv, na.rm = TRUE), mat = mean(mat, na.rm = TRUE), 
            heterogen = sum(heterogen), 
            year = mean(year), elevation = mean(elevation)) |>
  ungroup() |>
  # filter(scale == "Scale_5") |>
  ggplot(aes(x = (heterogen), y = N, colour = burnt)) + 
  geom_point() + 
  labs(x = "log_Heterogeneity", y = "Shrub abundance") +
  theme_bw() +
  facet_wrap(.~site, scales = "free") 

 
df |>
  filter(burnt == 1) |>
  pivot_wider(names_from = scale, values_from = heterogen) |>
  ggplot(aes(x = as.factor(year), y = (juv), fill = site)) +
  geom_boxplot() + 
  geom_jitter(width = .15, alpha = .5) + 
  labs(x = "Year", y = "N, juvenile") +
  theme_bw()

dft |>
  ggplot(aes(x = Scale_6, y = log(juv+1), colour = site)) + 
  labs(x = "log( Heterogeneity )", y = "N, juvenile") +
  geom_point() +
  theme_bw() +
  facet_wrap(.~burnt, scales = "free")

# ---
dft |>
  filter(burnt == 1) |> dplyr::select(-site, -burnt, -cid) |>
  cor() |>corrplot::corrplot()

# === DAG schemes
# DAG for m2.0 is obsoloete. Not necessary since the model is not causal in any sense.
# dag2.0 <- dagitty::dagitty('dag{
#                           Y[pos="0,0"]
#                           J0[pos="0,1"]
#                           M[pos="1,2"]
#                           N[pos="1,1"]
#                           H[pos="2,1"]
#                           J1[pos="2,0"]
#                           U[pos="1,0.5"]
#                           
#                           U <- Y -> J0 -> M -> N
#                           J0 -> N -> H
#                           H -> J1
#                           U -> H
# }')
# plot(dag2.0)
# dagitty::impliedConditionalIndependencies(dag2.0)
# 
# dagitty::paths(dag2.0, "J1", "H")
# === 
dag3.0 <- dagitty::dagitty('dag{
                          J[pos="0,.5"]
                          H2[pos=".5,.66"]
                          H1[pos=".5,.88"]
                          M[pos="1,.5"]
                          Y[pos="1,0"]
                          Y[pos="1,0"]
                          S[pos="1.25,.25"]
                          E[pos="1.5,.5"]

                          H1 -> J <- H2 <- M <- E
                          H2 <- H1
                          H1 <- M -> J <- Y -> M
                          Y <- S -> E
                          
}')
plot(dag3.0)
dagitty::impliedConditionalIndependencies(dag3.0)
dagitty::paths(dag3.0, from = c("M"), to = "J", Z = c("H2", "S", "H1"))
# ===


# === Models 2 and 3
library(brms)

# --- variable transformations
dft |> 
  mutate(across(starts_with("Scale_"), ~ .x/ sd(.x), 
                .names = "{.col}sd"),
         tot_heterogen_sd = tot_heterogen/sd(tot_heterogen)) |>
  mutate(logN = log(N+1), 
              logMat = log(mat+1) / 2*sd(log(mat+1)),
              yearsince = abs(year - 2022) / 2*sd(abs(year - 2022)), 
              elevsc = (elevation - mean(elevation))/(sd(elevation)),
         rugosity = rugosity, 
         rugosity_sc = rugosity/sd(rugosity)) -> dft
# ---

# === Model 2: sparse model for variable (scale) selection
# --- formula
x1 <- "1 + "
x2 <- paste0("Scale_", 1:9, "sd", collapse = " + ")

f2 <- brms::brmsformula(paste("juv ~ ", x1, x2))

m2.0 <- brm(f2, data = dft, family = negbinomial(), 
            prior = set_prior(horseshoe(df = 1, df_global = 2, par_ratio = 0.01)), 
            control = list(adapt_delta = 0.95), 
            cores = 4, iter = 2000)
# saveRDS(m2.0, "../wavelet_models/m2_0.rds")
m2.0 <- readRDS("../wavelet_models/m2_0.rds")

mcmc_plot(m2.0, pars = "^b_")
mcmc_plot(m2.0, pars = "^b_Scale_7")

# === Predictive power of m2.0
# --- test the predictive power (r2 and MAE) using two ways:
# - 1. successive removal of higher-order heterogeneity
# - 2. using the selected variables from the range of scales, one at a time

# - 1. successive removal using k-fold/site validation
#m2.1 <- list() # level 1 = variable; level 2 = site
#mpred2.1 <- list()
for(i in 1:9){
  x2.1a <- paste0("1 + elevsc + ")
  x2.1b <- paste0("Scale_", 1:i, "sd", collapse = " + ")
  f2.1 <- brms::brmsformula(paste("juv ~ ", x2.1a, x2.1b))
  l <- list(); l2 <- list()
  for(j in 1:nlevels(dft$site)){
    dft |> 
      mutate(site = as.numeric(site)) |>
      filter(site != j) -> dat
  
    l[[j]] <- brm(f2.1, data = dat, family = negbinomial(), 
              prior = set_prior(horseshoe(df = 1, par_ratio = 0.1)), 
              control = list(adapt_delta = 0.95), 
              cores = 4, iter = 2000)
    l2[[j]] <- predict(l[[j]])
  }
  m2.1[[i]] <- l
  mpred2.1[[i]] <- l2
}
# saveRDS(m2.1, "../wavelet_models/m2_1.rds")
m2.1 <- readRDS("../wavelet_models/m2_1.rds")

acclist <- list()
for(j in 1:9){
  print(paste0("j = ", j))
  acc2.1 <- data.frame(scale = j, k = 1:10, 
                       r2 = NA, r2l = NA, r2u = NA, 
                       elpd = NA, looic = NA, 
                       mae = NA, rmse = NA)
  for(i in 1:10){
    print(paste0("i = ", i))
    acc2.1[i,c('r2', 'r2l', 'r2u')] <- bayes_R2(m2.1[[j]][[i]])[c(1,3,4)]
    acc2.1[i,c('elpd', 'looic')] <- loo(m2.1[[j]][[i]])[[1]][c(1,3),1]
    
    dft |> 
      mutate(site = as.numeric(site)) |>
      filter(site == i) -> ndat
    
    yhat <- predict(m2.1[[j]][[i]], newdata = ndat)
    acc2.1[i, c('mae', 'rmse')] <- c( mean(abs(ndat$juv-yhat)), sqrt(mean((ndat$juv-yhat)^2)) )
  }
  acclist[[j]] <- acc2.1
}
acc2.1 <- do.call(bind_rows, acclist)
write.csv(acc2.1, "../wavelet_models/accuracy2_1.csv", row.names = FALSE)


# --- 2. predicting using optimal resolution (opt resolution for prediction)
m2.2 <- list()
mpred2.2 <- list()
for(i in 1:9){
  x2.2 <- paste0("1 + elevsc + Scale_", i,"sd")
  # x2.2 <- paste0("1 + elevsc + tot_heterogen_sd")
  f2.2 <- brms::brmsformula(paste("juv ~ ", x2.2))
  
  l <- list(); l2 <- list()
  for(j in 1:nlevels(dft$site)){
    dft |> 
      mutate(site = as.numeric(site)) |>
      filter(site != j) -> dat
    
  l[[j]] <- brm(f2.2, data = dat, family = negbinomial(), 
                   cores = 4, iter = 2000)
  l2[[j]] <- predict(l[[j]])
  }
  m2.2[[i]] <- l
  mpred2.2[[i]] <- l2
}
# saveRDS(m2.2, "../wavelet_models/m2_2.rds")
# m2.2 <- readRDS("../wavelet_models/m2_2.rds")

acclist <- list()
for(j in 1:1){#9
  print(paste0("j = ", j))
  acc2.2 <- data.frame(scale = j, k = 1:10, 
                     r2 = NA, r2l = NA, r2u = NA, 
                     elpd = NA, looic = NA, 
                     mae = NA, rmse = NA)
  for(i in 1:10){
    print(paste0("i = ", i))
  acc2.2[i,c('r2', 'r2l', 'r2u')] <- bayes_R2(m2.2[[j]][[i]])[c(1,3,4)]
  acc2.2[i,c('elpd', 'looic')] <- loo(m2.2[[j]][[i]])[[1]][c(1,3),1]
  
  dft |> 
    mutate(site = as.numeric(site)) |>
    filter(site == i) -> ndat
  
  yhat <- predict(m2.2[[j]][[i]], newdata = ndat)
  acc2.2[i, c('mae', 'rmse')] <- c( mean(abs(ndat$juv-yhat)), sqrt(mean((ndat$juv-yhat)^2)) )
  }
  acclist[[j]] <- acc2.2
}
acc2.2 <- do.call(bind_rows, acclist)
write.csv(acc2.2, "../wavelet_models/accuracy2_2.csv", row.names = FALSE)


# === Model 3: using m2.0 (above) find if SH amplifies the recovery (DAG 3.0)
f3 <- brms::brmsformula("juv ~ 1 + (1|site)  + logMat + Scale_1sd + logMat*Scale_1sd + Scale_2sd")

m3.2 <- brm(f3, 
            data = dft |> 
              filter(burnt == 1, mat > 0) |> #mat > 0,
              mutate(logMat = log(mat+1)/sd(log(mat+1)), 
                     Scale_1sd = Scale_1/sd(Scale_1)), 
            family = negbinomial(), 
            control = list(adapt_delta = 0.95), 
            cores = 4, iter = 2000)
saveRDS(m3.2, "../wavelet_models/m3_2.rds")

mcmc_plot(m3.2, pars = "^b_")
plot(conditional_effects(m3.1, effects = c("Scale_1sd:logMat")), points = TRUE )
# ====================================
# --- 

# === export ggplot
ggsave("figures/demo_fig_2.pdf", width = 7, height = 4)
