library(rstan)
library(ggplot2)
library(ggpubr)
library(readxl)
library(R2BayesX)
library(tidyverse)
library(Matrix)
setwd("C:\\Users\\Dawodu oluwaegun\\Desktop\\Spatial")
secd = read.csv('main_extract.csv', header=TRUE)
##secd$state
y1 <- secd$mobidity
y2 <- secd$mortality
N=length(y1)
state <- c("Abia.1",        "Abuja.1" , "Adamawa.1",     "Akwa Ibom.1" ,  "Anambra.1"  ,   "Bauchi.1",      "Bayelsa.1",    
           "Benue.1",       "Borno.1" ,  "Cross River.1", "Delta.1",       "Ebonyi.1",      "Edo.1",         "Ekiti.1",      
           "Enugu.1",       "Gombe.1"  ,     "Imo.1",         "Jigawa.1",      "Kaduna.1",    "Kano.1" ,       "Katsina.1",    
           "Kebbi.1",       "Kogi.1"  ,      "Kwara.1",     "Lagos.1" ,    "Nassarawa.1" , 
           "Niger.1",       "Ogun.1" ,       "Ondo.1",       "Osun.1" ,       "Oyo.1"    ,     "Plateau.1"   ,  "Rivers.1",     
           "Sokoto.1",      "Taraba.1",      "Yobe.1",        "Zamfara.1")

numBUGS= c(22,37,8,1,2,3,31,5,6,7,23,32,4,33,24,34,9,25,10,11,12,26,27,13,14,35,15,16,17,28,18,19,20,21,29,30,36)
##sstateid = c(22)
nstates = 37 #length(state)

##### fff = https://github.com/eosafu/Poly--GammaInference/blob/main/publishedcode.R
###k=length(numBUGS)

W <- read.csv("https://raw.githubusercontent.com/eosafu/Model/master/Nig_adj_Matrix.txt",sep=",",header = F) # Nigeria adjacency matrix
colnames(W) <- NULL
W <- as(as.matrix(W[,-1]),"sparseMatrix")
image(W)
# Compute spatial effect precision matrix
P <- as(diag(rowSums(W)),"sparseMatrix")
q <- nrow(P)
rho <- 0.9
D <- P-rho*W
D_W_inv = solve(D)
W = matrix(0, nrow = nstates, ncol = nstates)

neighbour_list = list(c(22, 20, 7, 1, 9, 2, 24, 32), # needs to be at the same order as numBUGS 
                      c(37, 15, 10, 35, 27), c(8, 29, 34, 6), c(1, 20, 22, 7), c(2, 27, 24, 22, 9, 23, 4, 20), 
                      c(3, 25, 11, 19, 29, 34, 10, 30), c(31, 20, 23), c(5, 35, 29, 27, 24, 32, 7), c(6, 30, 34, 8), 
                      c(7, 5, 32, 22, 1), c(23, 17, 4, 20, 2, 31), c(32, 24, 5, 7, 22), c(4, 17, 23, 27, 2), 
                      c(33, 28, 27, 13, 27), c(24, 27, 5, 32, 2, 22), c(34, 3, 30, 6, 8, 29), c(9, 20, 22, 2), 
                      c(25, 11, 3, 30, 12), c(10, 15, 36, 12, 11, 3, 19, 35, 37), c(11, 25, 12, 3, 10), 
                      c(12, 36, 10, 11, 25), c(26, 15, 36, 21), c(27, 13, 15, 37, 35, 5, 24, 2, 33, 17, 4), 
                      c(13, 15, 18, 28, 33, 27), c(14, 16), c(35, 10, 37, 19, 29, 5, 11), c(15, 13, 27, 37, 10, 36, 26), 
                      c(16, 14, 28, 18, 17), c(17, 16, 28, 33, 4, 23, 27), c(28, 18, 16, 17, 33, 13), c(18, 13, 28, 16), 
                      c(19, 10, 3, 35, 29), c(20, 23, 31, 9, 22, 2, 1), c(21, 26, 36), c(29, 5, 35, 19, 3, 34, 8), 
                      c(30, 3, 34, 6), c(36, 26, 15, 21, 12, 10)) 

for(i in 1:nstates){
  where = NULL
  num_neighbours = length(neighbour_list[[i]][-1])
  
  for(j in 1:nstates){
    for(k in 1:num_neighbours){
      if(numBUGS[j] == neighbour_list[[i]][k + 1]){
        where[k] = j
      }
    }
  }
  W[i,where] = 1
}
W=forceSymmetric(W)
D = diag(rowSums(W))
D_W_inv = solve(D-0.9*W) # enter on stan as data
K=secd$state
###sort(unique(K))
###K=secd$sstate
##sort(unique (k))
##phi=D_W_inv
y=c(y1, y2)
X1<-secd$prim
X2<-secd$urban
X3=secd$sec
x<-as.matrix(data.frame(X1=secd$prim, X2=secd$urban, X3=secd$sec))
z<-as.matrix(data.frame(X1=secd$prim, X2=secd$urban, X3=secd$sec))
#K=ncol(x)
#z<-as.matrix(data.frame(X1=secd$prim, X2=secd$urban))
colnames(x) = NULL
#data_list <- list(N=N, y=y, x=secd$prim)
#3data_list <- list(N=N, y=y, x=x, z=z, K=2, D=2)
data_list <- list(N=N, y=y, x=x, z=z, K=K, mu=rep(0, 37), D_W_inv=as.matrix(D_W_inv))
stan_model <- stan(file = "Main_code.stan", 
                   data = data_list,
                   iter = 10,
                   warmup = 1,
                   thin = 1,
                   chains = 1)
print(summary(stan_model))
plot(stan_model)
traceplot(stan_model, pars = c("beta1", "beta10"), inc_warmup = TRUE, nrow = 2)
sampler_params <- get_sampler_params(stan_model, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)
lapply(sampler_params, summary, digits = 2)
control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE)
pairs(stan_model, pars = c("beta1", "beta10", "lp__"), las = 1)
posterior_samples <- as.array(stan_model)
sigma_est <- mean(posterior_samples[, , "sigma"])
