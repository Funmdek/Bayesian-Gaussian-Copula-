library(rstan)
library(ggplot2)
library(ggpubr)
library(readxl)
library(R2BayesX)
library(tidyverse)
library(Matrix)
setwd("C:\\Users\\14483065\\Documents\\Spatial")
secd = read.csv('main_extract.csv', header=TRUE)
##secd$state
y1 <- secd$mobidity
y2 <- secd$mortality
N=length(y1)
W <- read.csv("https://raw.githubusercontent.com/eosafu/Model/master/Nig_adj_Matrix.txt",sep=",",header = F) # Nigeria adjacency matrix
colnames(W) <- NULL
W <- as(as.matrix(W[,-1]),"sparseMatrix")
image(W)
# Compute spatial effect precision matrix
P <- as(diag(rowSums(W)),"sparseMatrix")
q <- nrow(P)
rho <- 0.99
D <- P-rho*W
D_W_inv = solve(D) # enter on stan as data
K=secd$state
y=c(y1, y2)
x<-as.matrix(data.frame(x1=secd$Male, X2=secd$prim, X3=secd$sec, X4=secd$gm, X5=secd$high, X6=secd$Poorer, x7=secd$Middle,
                        x8=secd$Richer, x9=secd$Richest, x10=secd$still_breast, x11=secd$irregular_breast, x12=secd$priv_hospital, 
                        x13=secd$govt_hospital, x14=secd$home, x15=secd$hausa, x16=secd$igbo, x17=secd$yoruba, x18=secd$islam, 
                        x19=secd$islam, x20=secd$christian, x21=secd$underweight, x22=secd$normalweight, x23=secd$overweight))
z<-as.matrix(data.frame(x1=secd$Male, X2=secd$prim, X3=secd$sec, X4=secd$gm, X5=secd$high, X6=secd$Poorer, x7=secd$Middle,
                        x8=secd$Richer, x9=secd$Richest, x10=secd$still_breast, x11=secd$irregular_breast, x12=secd$priv_hospital, 
                        x13=secd$govt_hospital, x14=secd$home, x15=secd$hausa, x16=secd$igbo, x17=secd$yoruba, x18=secd$islam, 
                        x19=secd$islam, x20=secd$christian, x21=secd$underweight, x22=secd$normalweight, x23=secd$overweight))
data_list <- list(N=N, y=y, x=x, z=z, K=K, mu=rep(0, 37), D_W_inv=as.matrix(D_W_inv))
stan_model <- stan(file = "Main_code.stan", 
                   data = data_list,
                   iter = 5000,
                   warmup = 2500,
                   thin = 1,
                   cores= 1,
                   chains = 4)
print(summary(stan_model))
####Access posterior samples
posterior_samples <- extract(stan_model)

# Diagnose convergence
stan_diag <- monitor(stan_model)

# Plot trace and density plots
traceplot(stan_model)
plot(stan_model)
traceplot(stan_model, pars = c("beta1", "beta10"), inc_warmup = TRUE, nrow = 2)
sampler_params <- get_sampler_params(stan_model, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)
lapply(sampler_params, summary, digits = 2)
control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE)
pairs(stan_model, pars = c("beta1", "beta10", "lp__"), las = 1)
posterior_samples <- as.array(stan_model)
sigma_est <- mean(posterior_samples[, , "sigma"])
