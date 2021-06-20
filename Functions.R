#############################################
#####    codes for data generating
##################################################
corr_setting<-function(p,q,corr_str,rho,G_id=1){
  if (G_id==1){
    if (corr_str=='AR'){
      sigmaG<-matrix(0,nrow=p,ncol=p)
      diag(sigmaG)<-rep(1,p)
      for(i in 1:p){
        for(j in 1:(i-1)){
          sigmaG[i,j]<-sigmaG[j,i]<-rho^(i-j)
        }
      }
    } else if (corr_str=='Band1'){
      sigmaG=matrix(0,p,p)
      diag(sigmaG)=0.5
      for (i in 1:(p-1)){
        sigmaG[i,i+1]=0.3
      }
      sigmaG=sigmaG+t(sigmaG)
    } else if (corr_str=='Band2'){
      sigmaG=matrix(0,p,p)
      diag(sigmaG)=0.5
      for (i in 1:(p-2)){
        sigmaG[i,i+1]=0.5
        sigmaG[i,i+2]=0.3
      }
      sigmaG=sigmaG+t(sigmaG)
    }
  } else {
    sigmaG=0
  }
  
  sigmaE<-matrix(0,nrow=q,ncol=q)
  diag(sigmaE)<-rep(1,q)
  for(i in 1:q){
    for(j in 1:(i-1)){
      sigmaE[i,j]<-sigmaE[j,i]<-0.3^(i-j)
    }
  }
  
  return(list(sigmaG=sigmaG,sigmaE=sigmaE))
}
data_generate <- function(n, p, q, k, sparsity, snr, rr, tau,corr_str) {
   centroids <- 1/(k:1) # 1-dimensional centroids with varing pairwise distances
   dist_centroids_min <- abs(diff(centroids)) 
   dist_centroids_min <- c(dist_centroids_min, dist_centroids_min[k-1]) # min pairwise distance for each centroid
   groups_assign <- c(rep(seq(1, k/2), each = ceiling(p/(2*k))), rep(seq(k/2 + 1, k), each = ceiling(3*p/(2*k))))
   if (p/(2*k) != ceiling(p/(2*k))) {
        p2 = length(groups_assign)
        groups_assign <- groups_assign[-match((k - (p2-p) + 1):k, groups_assign)]
   }
   membership_assign <- matrix(0, nrow = p, ncol = k)
   for (i in 1:k) membership_assign[which(groups_assign == i), i] <- 1
    # Generate hier tree from p latent vectors
   ZZ <- membership_assign %*% centroids + tau * (membership_assign %*% dist_centroids_min) * rnorm(p)
   hc <- hclust(dist(ZZ))
   hc$height <- sqrt(hc$height)
   sigma_GE=corr_setting(p,q,corr_str,0.3)
   X<-mvrnorm(n,rep(0,q),sigma_GE$sigmaE)
   X[X[,1]>0,1]=1
   X[X[,1]<=0,1]=0
   X[X[,4]>0,4]=1
   X[X[,4]<=0,4]=0
   Z <-mvrnorm(n,rep(0,p),sigma_GE$sigmaG)
   for (j in 1:p){
        temp1=quantile(Z[,j],0.98)
        temp2=quantile(Z[,j],0.995)
        temp3=Z[,j]
        Z[temp3<=temp1,j]=0
        Z[(temp3>temp1)&(temp3<=temp2),j]=1
        Z[temp3>temp2,j]=2
    }
    Z.pool <- Z %*% membership_assign
    group.support <- rep(1, len = k)
      if (sparsity > 0) {
        group.support[(k/2-k*sparsity/2+1):(k/2)] <- 0 # Zero out k*sparsity/2 smaller groups
        group.support[(k/2+1):(k/2+k*sparsity/2)] <- 0 # Zero out k*sparsity/2 larger groups
      }
      beta.pool <- group.support * runif(k, 0.8, 1.5) * rep(c(0, 1,0,0,0), k/5)
      ###############################################################
      eta_1.pool <- group.support * runif(k,0.8, 1.5) 
      eta_2.pool <- group.support * runif(k, 0.8, 1.5) 
      eta_3.pool <- group.support * runif(k, 0.8, 1.5) 
      eta_4.pool <- group.support * runif(k, 0, 0) 
      eta_5.pool <- group.support * runif(k, 0, 0) 
      
      ################################################################  
      alpha_true=matrix(runif(q,0.8,1.2),q,1)
      
      b_true=matrix(0,q+1,k)
      b_true[1,]=beta.pool
      
      xi_1= eta_1.pool*beta.pool
      xi_2= eta_2.pool*beta.pool
      xi_3= eta_3.pool*beta.pool
      xi_4= eta_4.pool*beta.pool
      xi_5= eta_5.pool*beta.pool
      
      b_true[2,]=xi_1
      b_true[3,]=xi_2
      b_true[4,]=xi_3
      b_true[5,]=xi_4
      b_true[6,]=xi_5
      b_vector=matrix(b_true,(q+1)*k,1)#######Z1 X1Z1 X2Z1......
      eta=cbind(eta_1.pool,eta_2.pool,eta_3.pool,eta_4.pool,eta_5.pool)
       ##############################################################  generate y
      pp=k*(q+1)
      W=matrix(0,n,pp)
      W[,seq(from=1,to=k*(q+1),by=(q+1))]=Z.pool
      for (i in 1:n){
        temp3=matrix(X[i,],q,1)%*%Z.pool[i,]
        W[i,setdiff(seq(from=1,to=k*(q+1),by=1),seq(from=1,to=k*(q+1),by=q+1))]=matrix(temp3,k*q,1)
      }
      WW=X%*%alpha_true+W%*%b_vector
       y=WW+matrix(rnorm( n), n, 1)
        beta_true=t(membership_assign%*%beta.pool)
      eta_true=t(membership_assign%*%t(b_true[-1,]))
      return(list(Y=y, X=X,Z=Z,Z.pool=Z.pool, HC=hc,beta=beta_true , eta=eta_true ,alpha=alpha_true,member=membership_assign)) # make each col its own list element
    }
    
#############################################
#####    codes for a tree-based gene-environment interaction analysis
##################################################
    S <- function(a,k){
      s <- a
      s[which(abs(s) <= k)] <- 0
      s[which(s > k)] <- s[which(s > k)] - k
      s[which(s < -k)] <- s[which(s < -k)] + k
      return(s)
    }
#######################################################################our approach
lasso_rare<- function(Z,X,y,lambda,alpha,hc,rho){
      n<-dim(X)[1]
      p<-dim(Z)[2]
      q<-dim(X)[2]
      W=array(0,dim=c(n,q,p))
      for (i in 1:n){                       ########################
        temp3=matrix(X[i,],q,1)%*%Z[i,]
        W[i,,]=temp3
      }
      beta0=matrix(0,p,1)
      theta0=matrix(0,q,p)
      alpha_c0=solve(t(X)%*%X+0.001*diag(1,q))%*%t(X)%*%y
      rr=y-X%*%alpha_c0
      #######################################################above is step 1
      loop_time=1  
      diff_v=1
      objective=mean(rr^2)/2
      beta1=beta0
      theta1=theta0
      loop_time=1  
      c=matrix(0,p,1)
      A0 <- tree.matrix(hc)
      nnodes <- ncol(A0)
      v_1=v_2=v_3=beta_1=beta_2=beta_3=beta0=rep(0,p)
      u_1=u_2=gamma_1=gamma_2=gamma=rep(0,nnodes)
      gamma_k= matrix(rep(rep(0,nnodes),q),nr=q)
      v_1k=v_2k=v_3k=xi_1=xi_2=xi_3=xi=rep(0,p)
      u_1k=u_2k=gamma_1k=gamma_2k=gammak=rep(0,nnodes)
      while ((diff_v>1e-2) && (loop_time<500)){###############CD +ADMM
        eta=matrix(0,n,p)
        for (j in 1:p){
          eta[,j]=Z[,j]+ rowSums((matrix(1,n,1)%*%t(theta0[,j]))*W[,,j])
        }
        y1=rr+eta%*%matrix(beta0,nc=1)
        fit <- rarefit(y = y1, X = eta, hc = hc, intercept = F, lambda=lambda,alpha = alpha, rho = 0.05, eps1 = 1e1, eps2 = 1e1, maxite = 1e1)
        beta1 <- fit$beta[[1]]
        gamma<- fit$gamma[[1]]
        rr=rr-eta%*%matrix((beta1-beta0),nc=1)
        beta0=beta1
        active_id=which(beta0!=0)
        tide_W= W
        for(k in 1:q){
          for(j in 1:p){
            tide_W[,k,j]= rowSums((matrix(1,n,1)%*%t(beta1[j]))*W[,k,j])
          }
        }
         for (k in 1:q){
          W_t= tide_W[,k,]
          xi_q= theta0
          xi_q[k,]=0
          w_yt=matrix(0,nr=n)
           for(kin in 1:q){
            w_yt0= tide_W[,kin,]%*%matrix(xi_q[kin,],nr=p)
            w_yt=w_yt+w_yt0
          }
           eta=matrix(0,n,p)
           for (j in active_id){
            eta[,j]=beta0[j]*W[,k,j]
          }
          y2=rr+eta%*%matrix(theta0[k,],nc=1)
          fit <- rarefit(y = y2, X = eta, hc = hc, intercept = F, lambda=lambda,alpha = alpha, rho = 0.05, eps1 =  1e1, eps2 = 1e1, maxite = 1e1)
          theta1[k,active_id] <- fit$beta[[1]][active_id]
          gamma_k[k,active_id]<- fit$gamma[[1]][active_id]
         rr=rr-eta%*%matrix((theta1[k,]-theta0[k,]),nc=1)
          theta0=theta1
        }
         ##################################################################above is step3
        temp=rr+X%*%alpha_c0
        alpha_c1=solve(t(X)%*%X+0.001*diag(1,q))%*%t(X)%*%temp
         rr=rr-X%*%(alpha_c1-alpha_c0)
        ##########################################################above is step4 
        alpha_c0=alpha_c1
        objective1=mean(rr^2)/2+lambda*(1-alpha)*(sum(abs(theta1))+sum(abs(beta1)))+lambda*alpha*(sum(abs(gamma_k))+sum(abs(gamma)))
         diff_v=abs(objective1-objective)/abs(objective)
        objective=objective1
         loop_time=loop_time+1  
       }
      ##########################################################above is step5 
      RSS=sum(rr^2)
       beta1[abs(beta1)<1e-2]=0
      beta_temp=(matrix(1,q,1)%*%t(beta1))*theta1
      beta_temp[abs(beta_temp)<1e-2]=0  
      df=sum(beta_temp!=0)+sum(beta1!=0)
      result=list(alpha=alpha_c1,beta=beta1,eta=beta_temp,RSS=RSS,df=df,diff_v=diff_v)
       return(result)
     }
  
  


