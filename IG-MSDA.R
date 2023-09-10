
## Model Size
MS = function(theta){
  p = nrow(theta)
  MS = 0
  for(j in 1:p){
    if (norm(theta[j,],'2') != 0){
      MS=MS+1
    }
  }
  return(MS)
}

## True positives
TP = function(theta_hat,theta_star){
  p = nrow(theta_hat)
  TP = 0
  for (j in 1:p){
    if (norm(theta_star[j,],'2')!=0 && norm(theta_hat[j,],'2')!=0){
      TP = TP + 1
    }
  }
  return(TP)
}

## False positives
FP = function(theta_hat,theta_star){
  p = nrow(theta_hat)
  FP = 0
  for (j in 1:p){
    if (norm(theta_star[j,],'2')==0 && norm(theta_hat[j,],'2')!=0){
      FP = FP + 1
    }
  }
  return(FP)
}

## False negatives
FN = function(theta_hat,theta_star){
  p = nrow(theta_hat)
  FN = 0
  for (j in 1:p){
    if (norm(theta_star[j,],'2')!=0 && norm(theta_hat[j,],'2')==0){
      FN = FN + 1
    }
  }
  return(FN)
}


cv.IG_MSDA = function(X,class,nlambda = 6,lambda = c(0.1,0.5),X_test){
  
  result_all = IG_MSDA(X,class,nlambda = nlambda,lambda = lambda,X_test)
  class_error = vector(length = nlambda)
  theta_error = vector(length = nlambda)
  
  for (i in 1:nlambda){
    class_error[i] = sum(result_all$class_pred[[i]] != class_test)/n_test
    theta_error[i] = norm(result_all$theta_hat[[i]] - theta_star,'F')
  }
  id.min = which(class_error == min(class_error))
  if (length(id.min) > 0){
    fn = vector(length = 0)
    for (ii in 1:(length(id.min))){
      fn[ii] = FN(result_all$theta_hat[[id.min[ii]]],theta_star)
    }
    id.min = id.min[which.min(fn)]
    if (length(id.min) > 0){
      id.min = id.min[which.min(theta_error[id.min])]
    }
  }
  
  cvoutput_theta = result_all$theta_hat[[id.min]]
  cvoutput_class_pred = result_all$class_pred[[id.min]]
  cvoutput_id.min = id.min
  
  cvoutput = list(result_all = result_all, id.min = id.min, 
                  theta_hat_opt = cvoutput_theta, class_pred_opt = cvoutput_class_pred)
  return(cvoutput)
}


#################################################################
IG_MSDA = function(X,class,nlambda = 6,lambda = c(0.1,0.5),X_test){
  
  G <- length(unique(class))
  p <- ncol(X)
  n <- nrow(X)
  n_test <- nrow(X_test)
  
  lambda = seq(lambda[1],lambda[2],length.out = nlambda)
  
  #------------------------------------------------------------------------
  ## Estimate pi
  pi_hat <- c(rep(0,G))
  for(k in 1:G){
    pi_hat[k] <- sum(class == k)/n
  }
  
  #------------------------------------------------------------------------
  ## Estimate mu
  # mu_hat p*G
  mu_hat <- matrix(0,p,G)
  for (k in 1:G){
    mu_hat[,k] <- apply(X[which(class == k),],2,mean)
  }
  
  #mu_hat[,k]:use all variables
  
  #------------------------------------------------------------------------
  # Estimate  Sigma
  Sig_hat <- matrix(0,p,p)
  for (i in 1:n){
    Sig_hat <- Sig_hat + (X[i,] - mu_hat[,class[i]]) %*% t(X[i,] - mu_hat[,class[i]])
  }
  Sig_hat <- Sig_hat/n
  Sig_hat <- Sig_hat + 1e-04*diag(p)
  
  # cholesky for Sig_hat
  val_Sig = eigen(Sig_hat)$value
  
  L = chol(Sig_hat)
  
  #------------------------------------------------------------------------
  # Estimate  delta
  delta_hat <- mu_hat[,-1] - mu_hat[,1]             
  
  #------------------------------------------------------------------------
  #                    
  neighb <- huge(X,lambda = huge.select(huge(X),criterion = "stars")$opt.lambda)$path[[1]]
  
  index <- function(x){which(x == 1)}
  neighbour <- apply(neighb,2,index)
  
  zeroindex = list()
  for (j in 1:p){
    zeroindex[[j]] = which(neighb[,j] == 0)
  }
  #------------------------------------------------------------------------
  # 
  phi = vector(length = p)
  for (j in 1:p){
    phi[j] <- sqrt(length(neighbour[[j]]))/norm(delta_hat[j,],'2')
  }
  
  #------------------------------------------------------------------------
  v_mat = Variable(p*(G-1),p)
  g = 2
  v_sum_j = sum_entries(v_mat[c(((g-2)*p+1):((g-1)*p)),],axis = 1,keepdims = T)
  L_n = 1/2 * p_norm(L %*% v_sum_j,2)^2 - (t(delta_hat[,g-1]) %*% v_sum_j)
  for (g in 3:G){
    v_sum_j = sum_entries(v_mat[c(((g-2)*p+1):((g-1)*p)),],axis = 1,keepdims = T)
    L_n = L_n + 1/2 * p_norm(L %*% v_sum_j,2)^2 - (t(delta_hat[,g-1]) %*% v_sum_j)
  }
  
  penalty = function(v_mat, lambda, G, p, phi){
    sum = 0
    for (j in 1:p){
      sum = sum + phi[j] * cvxr_norm(v_mat[,j],2)
    }
    lambda * sum
  }
  
  output_theta = list()
  output_class_pred = list()
  
  constr = list()
  for (j in 1:p){
    for(g in 2:G){
      constr = append(constr,v_mat[zeroindex[[j]]+(g-2)*p,j] == 0)
    }
  }
  
  for (il in 1:nlambda){
    obj = L_n + penalty(v_mat, lambda[il], G, p, phi)
    prob = Problem(Minimize(obj),constr)
    result = solve(prob)
    
    v_hat = round(result$getValue(v_mat),4)
    theta_hat = matrix(0,p,G-1)
    for (g in 2:G){
      for(j in 1:p){
        theta_hat[,g-1] = theta_hat[,g-1] + v_hat[c((1+(g-2)*p):((g-1)*p)),j]
      }
    }
    
    D_hat = matrix(0,n_test,G-1)
    for(i in 1:n_test){
      for(g in 2:G){
        D_hat[i,g-1] = t(X_test[i,] - (mu_hat[,1]+mu_hat[,g])/2) %*% theta_hat[,g-1] + log(pi_hat[g]/pi_hat[1])
      }
    }
    
    class_pred = vector(length = n_test)
    for (i in 1:n_test){
      if (all(D_hat[i,] < 0)){
        class_pred[i] = 1
      }else{
        class_pred[i] = which.max(D_hat[i,]) + 1
      }
    }
    
    output_theta[[il]] = theta_hat
    output_class_pred[[il]] = class_pred
  }
  
  output = list(lambda = lambda,theta_hat = output_theta,class_pred = output_class_pred,prior = pi_hat)
  return(output)
}
