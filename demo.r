########################################################################
source("simu_2drarma.r")
source("fit_2drarma.r")

# fit.2drarma(y, ar = NA, ma = NA, ne = 1)
# y is the evaluated image
# ar and ma are the number of autoregressive and moving average parameters
# n is the number of rows
# k is the number of columns 
# ne is related to the employed neighborhood 
# ne = 1: strongly causal neighborhood (used in Palm, Bayer, and Cintra ?) 
# ne = 2: non causal region
# ne = 3: semi-causal

############################### ne = 1 #################################

alpha =  0.3569 # intercept
phi = c(0.2155, 0.2032, 0.15) # AR parameters
theta = c(0.1529, 0.1744, 0.1998) # MA parameters

y = simu.2drarma(alpha, phi, theta, p = 1, q = 1, n = 10 , k = 10, ne = 1) # ARMA simulation

image(y,col=gray((0:255)/255)) # image plot (tem uma rotacao estranha)
library(raster)
plot(raster(y), col = gray.colors(255)) # better plot

fit.2drarma(y, 1, 1, ne = 1) # adjust
fit.2drarma(y, 2, 2, ne = 1) # adjust
fit.2drarma(y, 1, 2, ne = 1) # adjust

y = simu.2drarma(alpha, NA, theta, p = NA, q = 1, n = 10 , k = 10, ne = 1) # MA simu
fit.2drarma(y, 1, NA, ne = 1)
fit.2drarma(y, 2, NA, ne = 1)

y = simu.2drarma(alpha, phi, NA, p = 1, q = NA, n = 10 , k = 10, ne = 1) # AR simu
fit.2drarma(y, NA, 1, ne = 1)
fit.2drarma(y, NA, 2, ne = 1)

############################### ne = 2 #################################

alpha = 0.3
phi = c(0.21, 0.015, 0.15, 0.5, 0.07, 0.2, 0.011, 0.21)
theta = c(0.1, 1.17, 0.02, 0.12, 0.34, 0.67, 0.042, 0.21)

simu.2drarma(alpha, phi, theta, p = 1, q = 1, n = 20 , k = 20, ne = 2)
fit.2drarma(y, 1, 1, ne = 2)
fit.2drarma(y, NA, 1, ne = 2)

############################### ne = 3 #################################

alpha = 0.3
phi = c(0.2155, 0.2032, 0.15, 0.2, 0.34)
theta = c(0.1529, 0.1744, 0.1998, 0.2, 0.34)

y= simu.2drarma(alpha, phi, theta, p = 1, q = 1, n = 10 , k = 10, ne = 3)

fit.2drarma(y, 1, 1, ne = 3)
fit.2drarma(y, 1, NA, ne = 3)
fit.2drarma(y, NA, 2, ne = 3)

########################################################################

source("fit_2drarma.r")

library(R.matlab)
library(raster)

im_1 = readMat("m_1_p_1.mat") # read the image
im_1 = im_1$im
im = raster(im_1)
plot(im, col = gray.colors(255))

y = im_1[588:(588+49),604:(604+49)] # selected region
plot(raster(y), col = gray.colors(255))

## 2DRARMA(1,1)
fit1 = fit.2drarma(y, 1, 1, ne = 1)
plot(raster(fit1$fitted), col = gray.colors(255))
fit2 = fit.2drarma(y, 1, 1, ne = 2)
plot(raster(fit2$fitted), col = gray.colors(255))
fit3 = fit.2drarma(y, 1, 1, ne = 3)
plot(raster(fit3$fitted), col = gray.colors(255))

## 2DRARMA(1,0)
fit4 = fit.2drarma(y, 1, NA, ne = 1)
plot(raster(fit4$fitted), col = gray.colors(255))
fit5 = fit.2drarma(y, 1, NA, ne = 2)
plot(raster(fit5$fitted), col = gray.colors(255))
fit6 = fit.2drarma(y, 1, NA, ne = 3)
plot(raster(fit6$fitted), col = gray.colors(255))

## 2DRARMA(0,1)
fit7 = fit.2drarma(y, NA, 1, ne = 1)
plot(raster(fit7$fitted), col = gray.colors(255))
fit8 = fit.2drarma(y, NA, 1, ne = 2)
plot(raster(fit8$fitted), col = gray.colors(255))
fit9 = fit.2drarma(y, NA, 1, ne = 3)
plot(raster(fit9$fitted), col = gray.colors(255))
