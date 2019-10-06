library(Matrix)
library(CVXR)
library(rrpack)

# The Tien's generation
p1 = p2 =  100
p = p1 + p2 
q = 24
r1 = 4
r2 = 6
n = 90

#correlation among the drugs, i.e. the noises are correlation w.r.t q
Omega = matrix(0.9, ncol = q, nrow = q); diag(Omega) = 1
# generate covariance matrix
Sigma = qr.solve(Omega)
out = eigen(Sigma, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
# generate noise
noise = matrix(rnorm(q*n), nrow = n, ncol = q) %*% S.sqrt



rho = 0.9
#correlation with each X_k
Omega1 = matrix(rho, ncol = p1, nrow = p1)
Omega2 = matrix(rho, ncol = p2, nrow = p2)
diag(Omega1) = 1
diag(Omega2) = 1
# generate covariance matrix
Sigma1 = qr.solve(Omega1)
Sigma2 = qr.solve(Omega2)
# generate covariance
out1 = eigen(Sigma1, symmetric = TRUE)
S.sqrt1 = out1$vectors %*% diag(out1$values^0.5) %*% t(out1$vectors)
out2 = eigen(Sigma2, symmetric = TRUE)
S.sqrt2 = out2$vectors %*% diag(out2$values^0.5) %*% t(out2$vectors)


msee = msee.rrr = mspe = mspe.rrr = enet.mse = enet.msp  = c()
for (ss in 1:30) {
  #Xlist = list( # need to be centered
  #  X1 = scale(matrix(rnorm(n*p1), nr = n, nc = p1) , scale=F) ,
  #  X2 = scale(matrix(rnorm(n*p2), nr = n, nc = p2) , scale=F) )
  
  # X is correlated
  Xlist = list( # need to be centered
    X1 = scale(matrix(rnorm(n*p1), nr = n, nc = p1) %*% S.sqrt1 , scale=F) ,
    X2 = scale(matrix(rnorm(n*p2), nr = n, nc = p2) %*% S.sqrt2 , scale=F) )
  
  # Setting S1: local low-rank
  #b = list(
  #  beta1 = matrix(rnorm(q*r1), nr = q, nc = r1) %*% matrix(rnorm(p1*r1), nr = r1,nc = p1),
  #  beta2 = matrix(rnorm(q*r2), nr = q, nc = r2) %*% matrix(rnorm(p2*r2), nr = r2,nc = p2)
  #)
  #beta0 = do.call(cbind, b)
  
  # Setting S2: low-rank + sparsity
  #b = list(
  # beta1 = matrix(rnorm(q*r1), nr = q, nc = r1) %*% matrix(rnorm(p1*r1), nr = r1,nc = p1),
  # beta2 = matrix(rbinom(q*p2, 1, 0.5) , nr = q ,nc = p2)*rnorm(q*p2,0)
  #)
  #beta0 = do.call(cbind, b)
  
  # Setting S3: GLOBAL LOW-RANK
  #r0 = 2
  #beta0 = t(matrix(rnorm(p*r0), nr = p, nc = r0) %*% matrix(rnorm(q*r0), nr = r0,nc = q))
  
  # Setting S4: sparsity
  beta0 = matrix(rbinom(q*p, 1, 0.2) , nr = q, nc = p)* matrix( rnorm(p*q) ,nr = q,nc = p)
  
  X = do.call(cbind, Xlist)
  Y = X%*%t(beta0) +  rnorm(q*n) #   # + noise
  
  # Calculate weights
  ww = rep(NA,length(Xlist))
  for (i in 1:length(Xlist)){
    ww[i] = svd(Xlist[[i]])$d[1]*(sqrt(q)
                                  + sqrt(as.numeric(rankMatrix(Xlist[[i]]))))/n
  }
  
  # Reshape X as list
  fit = iRRR(Y, Xlist, lam1= .001,
             paramstruct = list(maxrho = 1,
                                rho = .00001,
                                Niter = 2000,
                                weight = ww,
                                Tol = 1e-7))
  ## checking the Mean squared errors    ############
  ha = do.call(rbind,fit$A)
  msee[ss] = mean((beta0 - t(ha))^2) 
  mspe[ss] = mean((Y - X%*%ha)^2) 
  
  # global low-rank
  fit.rr = cv.rrr(Y, X, nfold = 5)
  be.rrr = fit.rr$coef
  mspe.rrr[ss] = mean((Y - X%*%be.rrr)^2);  
  msee.rrr[ss] = mean((t(beta0) - be.rrr)^2);
  
  y_star <- as.vector(Y)
  x_star <- kronecker(Diagonal(q), X)  #cbind(rep(1,n), X)
  library(glmnet)
  fit.enet <- cv.glmnet(x_star, y_star, 
                   intercept = F, 
                   alpha = 0.2, 
                   standardize.response = F)
  enet.pred =  predict(fit.enet, newx = x_star,s="lambda.min")
  enet.msp[ss] = mean((y_star - enet.pred)^2)
  enet.b = matrix(coef(fit.enet, s="lambda.min")[-1], nrow = q, ncol = p)
  enet.mse[ss] = mean((beta0 - enet.b)^2)
}
mean(mspe)
mean(mspe.rrr)
mean(enet.msp)
mean(msee)
mean(msee.rrr)
mean(enet.mse)









library(CVXR)
beta1 <- Variable(rows = p1,cols = q); beta2 <- Variable(rows = p2,cols = q)
loss <- sum_squares(Y - Xlist$X1%*%beta1 - Xlist$X2%*%beta2 ) /n
lamda = .01
obj <- Minimize(loss + lamda*cvxr_norm(beta1,'nuc')*ww[1] 
                + lamda*cvxr_norm(beta2,'nuc')*ww[2] )
prob <- Problem(obj);   a <- solve(prob)
b1 = a$getValue(beta1); b2 = a$getValue(beta2)
beta.h = rbind(b1,b2)
# mean squared error  # can NOT computed w real data
# mean squared prediction error # can be computed w real data
mean((Y - X%*%beta.h)^2);mean((t(beta0) - beta.h)^2)






# test with nuclear norm
library(CVXR)
teta <- Variable(rows = p,cols = q)
loss <- sum_squares(Y - X %*% teta)/n 
lamda = .01
obj.nuclear <- Minimize(loss + lamda*cvxr_norm(teta,'nuc'))
prob.nuclear <- Problem(obj.nuclear)
aa <- solve(prob.nuclear)
beta.hh = aa$getValue(teta)
# mean squared error  # can NOT computed w real data
mean((Y - X%*%beta.hh)^2);  mean((t(beta0) - beta.hh)^2);








#correlation among the drugs, i.e. the noises are correlation
Omega = matrix(0.3, ncol = q, nrow = q)
diag(Omega) = 1
# generate covariance matrix
Sigma = qr.solve(Omega)
# generate noise
Z = matrix(rnorm(q*n), nrow = n, ncol = q)
out = eigen(Sigma, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
noise = Z %*% S.sqrt

Y = X%*%t(beta0) + noise