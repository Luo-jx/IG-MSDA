library(MASS)
library(mvtnorm)
library(Matrix)
library(huge)
library(ggplot2)
library(msda)
library(sparseLDA)
library(penalizedLDA)
library(Rcpp)
library(pROC)
library(CVXR)

source("IG-MSDA.R")

#################################################
## case 1 
G = 3
p = 150
n1 = 50
n2 = 50
n3 = 50
n = n1 + n2 + n3
b = 30    
pi <- c(1/3,1/3,1/3)


mu_star = matrix(0,p,G)
mu_star[c(5,10,15,20,25),2] = -1
mu_star[,3] = -mu_star[,2]


A_B <- 0.3*diag(1,5,5) + 0.7*matrix(rep(1,25),5,5)
Omega_star <- as.matrix(bdiag(bdiag(rep(list(A_B),b)),diag(p-b*5)))
Sigma_star <- solve(Omega_star)
delta_star <- mu_star[,-1] - mu_star[,1]
theta_star <- Omega_star %*% delta_star

X1 = mvrnorm(n1,mu_star[,1],Sigma_star)
X2 = mvrnorm(n2,mu_star[,2],Sigma_star)
X3 = mvrnorm(n3,mu_star[,3],Sigma_star)
X = rbind(X1,X2,X3)
y1 = rep(1,n1)
y2 = rep(2,n2)
y3 = rep(3,n3)
data = rbind(cbind(X1,y1),cbind(X2,y2),cbind(X3,y3))
class = data[,p+1]


n1_test = 300
n2_test = 300
n3_test = 300
n_test = n1_test + n2_test + n3_test
X1_test = mvrnorm(n1_test,mu_star[,1],Sigma_star)
X2_test = mvrnorm(n2_test,mu_star[,2],Sigma_star)
X3_test = mvrnorm(n3_test,mu_star[,3],Sigma_star)
X_test = rbind(X1_test,X2_test,X3_test)
y1_test = rep(1,n1_test)
y2_test = rep(2,n2_test)
y3_test = rep(3,n3_test)
data_test = rbind(cbind(X1_test,y1_test),cbind(X2_test,y2_test),cbind(X3_test,y3_test))
class_test = data_test[,p+1]

f = matrix(0, n_test, G)
for (jj in 1:G){
  f[,jj] = dmvnorm(X_test,mu_star[,jj],Sigma_star)
}

obj_IG_MSDA = cv.IG_MSDA(X,class,nlambda = 2,lambda = c(0.52,0.55),X_test = X_test)
MCR_IG_MSDA = sum(obj_IG_MSDA$class_pred_opt != class_test)/n_test
MS_IG_MSDA = MS(obj_IG_MSDA$theta_hat_opt)
TP_IG_MSDA = TP(obj_IG_MSDA$theta_hat_opt,theta_star)
FP_IG_MSDA = FP(obj_IG_MSDA$theta_hat_opt,theta_star)

#####################################################
## Case 2

G = 3
p = 150
n1 = 50
n2 = 50
n3 = 50
n = n1 + n2 + n3
pi <- c(1/3,1/3,1/3)
iter = 100              

mu_star = matrix(0,p,G)
mu_star[c(1:5),2] = 2.4
mu_star[c(6:10),2] = 1.2
mu_star[,3] = -mu_star[,2]


Omega_star = matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    if(abs(i-j) == 1)
      Omega_star[i,j] = -2/3
  }
} 

Omega_star = Omega_star + diag(5/3,p,p)
Omega_star[1,1] = 4/3
Omega_star[p,p] = 4/3
set.seed(100)
huge.plot(Omega_star)

Sigma_star <- solve(Omega_star)
delta_star <- mu_star[,-1] - mu_star[,1]
theta_star <- Omega_star %*% delta_star

X1 = mvrnorm(n1,mu_star[,1],Sigma_star)
X2 = mvrnorm(n2,mu_star[,2],Sigma_star)
X3 = mvrnorm(n3,mu_star[,3],Sigma_star)
X = rbind(X1,X2,X3)
y1 = rep(1,n1)
y2 = rep(2,n2)
y3 = rep(3,n3)
data = rbind(cbind(X1,y1),cbind(X2,y2),cbind(X3,y3))
class = data[,p+1]

n1_test = 300
n2_test = 300
n3_test = 300
n_test = n1_test + n2_test + n3_test
X1_test = mvrnorm(n1_test,mu_star[,1],Sigma_star)
X2_test = mvrnorm(n2_test,mu_star[,2],Sigma_star)
X3_test = mvrnorm(n3_test,mu_star[,3],Sigma_star)
X_test = rbind(X1_test,X2_test,X3_test)
y1_test = rep(1,n1_test)
y2_test = rep(2,n2_test)
y3_test = rep(3,n3_test)
data_test = rbind(cbind(X1_test,y1_test),cbind(X2_test,y2_test),cbind(X3_test,y3_test))
class_test = data_test[,p+1]

f = matrix(0, n_test, G)
for (jj in 1:G){
  f[,jj] = dmvnorm(X_test,mu_star[,jj],Sigma_star)
}

obj_IG_MSDA = cv.IG_MSDA(X,class,nlambda = 5,lambda = c(0.5,0.65),X_test = X_test)
MCR_IG_MSDA = sum(obj_IG_MSDA$class_pred_opt != class_test)/n_test
MS_IG_MSDA = MS(obj_IG_MSDA$theta_hat_opt)
TP_IG_MSDA = TP(obj_IG_MSDA$theta_hat_opt,theta_star)
FP_IG_MSDA = FP(obj_IG_MSDA$theta_hat_opt,theta_star)

#####################################################
## Case 3
G = 3
p = 150
n1 = 50
n2 = 50
n3 = 50
n = n1 + n2 + n3
pi <- c(1/3,1/3,1/3)
iter = 100


Omega_star = matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    if(abs(i-j) != 0)
      Omega_star[i,j] = -0.5*rbinom(1,1,0.03)
  }
}
Omega_star = Omega_star + 3.5*diag(p)
Omega_star = t(Omega_star)/2 + Omega_star/2

set.seed(100)
huge.plot(Omega_star)

mu_star = matrix(0,p,G)
mu_star[order(apply(Omega_star,2,sum),decreasing = T)[1:5],2] = 0.8
mu_star[,3] = -mu_star[,2]

Sigma_star <- solve(Omega_star)
delta_star <- mu_star[,-1] - mu_star[,1]
theta_star <- Omega_star %*% delta_star

X1 = mvrnorm(n1,mu_star[,1],Sigma_star)
X2 = mvrnorm(n2,mu_star[,2],Sigma_star)
X3 = mvrnorm(n3,mu_star[,3],Sigma_star)
X = rbind(X1,X2,X3)
y1 = rep(1,n1)
y2 = rep(2,n2)
y3 = rep(3,n3)
data = rbind(cbind(X1,y1),cbind(X2,y2),cbind(X3,y3))
class = data[,p+1]


n1_test = 300
n2_test = 300
n3_test = 300
n_test = n1_test + n2_test + n3_test
X1_test = mvrnorm(n1_test,mu_star[,1],Sigma_star)
X2_test = mvrnorm(n2_test,mu_star[,2],Sigma_star)
X3_test = mvrnorm(n3_test,mu_star[,3],Sigma_star)
X_test = rbind(X1_test,X2_test,X3_test)
y1_test = rep(1,n1_test)
y2_test = rep(2,n2_test)
y3_test = rep(3,n3_test)
data_test = rbind(cbind(X1_test,y1_test),cbind(X2_test,y2_test),cbind(X3_test,y3_test))
class_test = data_test[,p+1]

f = matrix(0, n_test, G)
for (jj in 1:G){
  f[,jj] = dmvnorm(X_test,mu_star[,jj],Sigma_star)
}

obj_IG_MSDA = cv.IG_MSDA(X,class,nlambda = 5,lambda = c(0.1,0.2),X_test = X_test)
MCR_IG_MSDA = sum(obj_IG_MSDA$class_pred_opt != class_test)/n_test
MS_IG_MSDA = MS(obj_IG_MSDA$theta_hat_opt)
TP_IG_MSDA = TP(obj_IG_MSDA$theta_hat_opt,theta_star)
FP_IG_MSDA = FP(obj_IG_MSDA$theta_hat_opt,theta_star)
