# This script is a tutorial of the Intro to Wavelets: 
# https://www.bioconductor.org/help/course-materials/2003/Milan/PDF/Lab9.pdf 

library(wavethresh)
library(Rwave)


# Display wavelet functions (example)
draw(filter.number = 5, family="DaubExPhase",  col = "blue", xlab = "", main = NULL) +
  abline(h = 0, v = 0, lty = 3, lwd = 1, col = "gray")

# Explore a multiresolution analysis
test.data1 <- example.1()$y

# conduct the DWT with the selected wavelet function 
# this fn by default decomposes the signal completely (ie up to a constant, signal average)
wds1 <- wd(test.data1, filter.number = 2, family="DaubExPhase")

# access and plot the coefficients
# accessC - extracts smoothing coefficients at a specified level
# accessD - extracts detail coefficients at a specified level
plot(accessC(wds1, level = 5), type = "l")

plot(test.data1, type = "l", lwd = 5)
points(accessC(wds1, level = 9), pch = 19, cex = .2, col = "purple")

# More illustrative example with the chirp data
test.data2 <- simchirp(512)$y

wds2 <- wd(test.data2)

par(mfrow = c(2, 2))
plot(accessC(wds2, level = 8), type = "l")
plot(accessD(wds2, level = 8), type = "l")
plot(accessC(wds2, level = 7), type = "l")
plot(accessD(wds2, level = 7), type = "l")
dev.off()

# An example with the CWT
# color intensity corresponds to the coefficient magnitude
par(mfrow = c(2, 1))
plot(test.data2, type = "l", lwd = 2)
cwtchrip <- cwt(test.data2, 5, 12)

# Explore the sparsity of wavelets (0 count at each level)
# experiment with sparsity threshold and wavelet functions
wds1N <- wd(test.data2, filter.number = 1, family="DaubExPhase")

levels <- 3:(wds1N$nlevels - 1)
nthresh <- length(levels)
d <- NULL
dz <- 0
for (i in 1:nthresh) {
       d <- accessD(wds1N, level = levels[i])
       d[abs(d) <= 0.001] = 0
       dz = dz + sum(d == 0)
       cat("Level: ", levels[i], " there are ", sum(d == 0), " zeroes\n")
   wds1Nd <- putD(wds1N, level = levels[i], v = d)
}

# inverse transform and the reconstruction error after "sparse'ing" the detail coeffiecients
par(mfrow = c(2, 1))
ts.plot(wr(wds1Nd))
ts.plot(test.data2 - wr(wds1Nd))

# Translation invariance example
# simulate data
k1 <- 32
k2 <- 42
par(mfrow = c(3, 1), mar = c(1.5, 1.5, 1.5, 0), mgp = c(5, 0.4,0))
x <- (1:128)/128
unshft <- c(rep(1, 64), rep(-1, 64))
noise <- rnorm(128)
y <- 2 * unshft + noise
plot(x, y, type = "l")
yest <- wr(threshold(wd(y, filter.number = 1, family = "DaubExPhase"), levels = 0:6))
lines(x, yest)

if (k1 == 0) y1 <- y else y1 <- c(y[(1 + k1):128], y[1:k1])
if (k1 == 0) shft1 <- unshft else shft2 <- c(unshft[(1 + k1):128],unshft[1:k1])
plot(x, y1, type = "l")
y1est <- wr(threshold(wd(y1, filter.number = 1, family = "DaubExPhase"), levels = 0:6))
lines(x, y1est)

if (k2 == 0) y2 <- y else y2 <- c(y[(1 + k2):128], y[1:k2])
if (k2 == 0) shft2 <- unshft else shft2 <- c(unshft[(1 + k2):128],unshft[1:k2])
plot(x, y2, type = "l")
y2est <- wr(threshold(wd(y2, filter.number = 1, family = "DaubExPhase"), levels = 0:6))
lines(x, y2est)

# remove coefficients for higher levels 
# for(i in 3:(wn1$nlevels - 1)) {
#   d <- rep(0, 2^i)
#   wn1 <- putD(wn1, level = i, v = d)
# }
# yhat0 <- wr(wn1)
# plot(yhat0, type = "l")
=
