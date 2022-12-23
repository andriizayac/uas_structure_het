library(wavethresh)

source("00_helper_functions.R")

# Test on DSM profiles in 1D (NS or EW)
td <- read.table("data/dem/res_64/matrix/crp_regen1_2021_dem_64_DOWN010.csv", header = FALSE, sep = ",") |> 
  as.matrix() |> 
  rev_raster()

image(td)
N <- nrow(td)
ns <- td[, N/2] # north-south
ew <- td[N/2, ] # east-west

n1 <- ns

# visualize an example of a DSM profile 
par(mfrow = c(3, 1))
plot(n1, type = "l", col = "magenta", lwd = 4, main = "DSM")
wn1 <- wd(n1, filter.number = 1, family = "DaubExPhase")
# yhat <- wr(threshold(wn1, levels = 2:3, type = "hard", 
#                      policy = "universal", by.level = TRUE))
# lines(yhat, type = "l", lty = "dashed")
plot(accessC(wn1, level = 2), type = "o", pch = 19, cex = 1, main = "Smoothing coefficients, level = 2")
plot(accessD(wn1, level = 2), type = "o", pch = 19, cex = 1, main = "Difference coefficients, level = 2")
dev.off()

# visualize an example of what will become a predictor matrix:
# smoothing coefs at level = 0
# difference coefs up to n.level-1
par(mfrow = c(8, 1), mar = c(0.5, 3, .5, 0))
plot(n1, type = "l", col = "magenta", lwd = 1, xaxt = "n")
plot(accessC(wn1, level = 0), type = "o", pch = 19, cex = 1, xaxt = "n")

En <- rep(0, 7)
En[1] <- 0 # accessC(wn1, level = 0)^2

for(i in 1:6) {
  plot( accessD(wn1, level = i), type = "o", pch = 19, cex = 1, xaxt = "n")
  En[i+1] <- sum( accessD(wn1, level = i)^2 )
}
dev.off()

# construct a 1D predictor (stack coefficients) for all observations
N <- length(n1)
n <- length(ew)
X <- matrix(NA, nr = n, nc = N-1)
for(m in 1:n){
  # use Haar wavelet because it handles boundaries well and sparsity is not of interest
  wn <- wd(ns[[m]], filter.number = 1, family = "DaubExPhase")
  
  D <- NULL
  for(i in (log2(N)-1):1) {
    D <- c(D, accessD(wn, level = i))
    # print(length(accessD(wn, level = i)))
  }
  X[m, 1:(N-2)] <- D
  X[m,(N-1)] <- accessC(wn, level = 0) 
}

# verify that the smoothing (level = 0, ie constants) correspond to the elevation gradient 
yh <- sapply(ew, mean)

plot(pts$elev ~ X[,N-1])

# check the agreement between the smoothing constants and elevation
summary(lm(pts$elev ~ X[,N-1]))

dat <- bind_cols(growth = pts$growth, X) %>% 
  filter(!is.na(growth))

# run model
library(brms)

mod_ns <- brm(growth ~ ., data = dat, 
           prior = set_prior(lasso(df = 1, scale = 10)))


res <- as.data.frame(mod) %>% 
  dplyr::select(contains("..."))

plot(resmu, pch = 19, cex = .5)

resmu <- apply(res, 2, mean)
reshpdi <- apply(res, 2, quantile, probs = c(0.33, 0.68))

sig <- apply(reshpdi, 2, function(x) { diff(sign(x)) }) + 1

plot(resmu, pch = 19, cex = .5, col = sig)




#### extra

# - perform wavelet transform and store prediction matrices
N <- 128
f <- list.files("data/DSM", pattern = paste0("dsm_2d_", N), full.names = TRUE)
n1 <- readRDS(f)$KOE003[N/2,]

dev.off()
par(mfrow = c(3, 1))
plot(n1, type = "l", col = "magenta", lwd = 4, main = "DSM")
wn1 <- wd(n1, filter.number = 1, family = "DaubExPhase")
# yhat <- wr(threshold(wn1, levels = 2:3, type = "hard", 
#                       policy = "universal", by.level = TRUE))
# lines(yhat, type = "l", lty = "dashed")
plot(accessC(wn1, level = 2), type = "o", pch = 19, cex = 1, main = "Smoothing coefficients, level = 2")
abline(h = 0, lty = "dashed")
plot(accessD(wn1, level = 3), type = "o", pch = 19, cex = 1, main = "Difference coefficients, level = 2")
abline(h = 0, lty = "dashed")

# difference coefs up to n.level-1
dev.off()
par(mfrow = c(wn1$nlevels+1, 1), mar = c(0.5, 4, .5, 0))
plot(n1, type = "l", col = "magenta", lwd = 1, xaxt = "n", ylab = "DSM")
plot(accessC(wn1, level = 0), type = "o", pch = 8, cex = 1, xaxt = "n", ylab = "Smooth_0")

for(i in 1:(wn1$nlevels-1)) {
  plot(accessD(wn1, level = i), type = "o", pch = 19, cex = .25, 
       xaxt = "n", ylab = paste0("Diff. l.=", i))
}

