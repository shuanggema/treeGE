library(MASS)
library(rare)

source('Functions.R')
tt=1
k=c(20,40,60,80,100)[tt]

p=200
q=5
n=200
 
Data<-data_generate(n = 200, p =p, q=q, k =k ,sparsity =0.4, snr = 5, rr = 50, tau=0.005,corr_str='AR')

lambda=0.05
alpha=0.75
    
Betatrue=matrix(Data$beta,nc=1)
Xitrue=Data$eta
alphatrue=Data$alpha

re=lasso_rare(Data$Z,Data$X,Data$Y,lambda,alpha,Data$HC,0.05) 
beta_our=re$beta
xi_our=re$eta
alpha_our=re$alpha
 
  
            
          