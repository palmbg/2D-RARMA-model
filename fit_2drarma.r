fit.2drarma = function (y, ar = NA, ma = NA, ne = 1)
{
  # ne is related to the employed neighborhood 
  # ne = 1: strongly causal neighborhood (used in Palm, Bayer, and Cintra ?) 
  # ne = 1: 2D-RARMA(p,q), 2D-RARMA(p,q), 2D-RARMA(p,q)
  # ne = 2: non causal region
  # ne = 2: 2D-RARMA(p,q), 2D-RARMA(p,q), 2D-RARMA(p,q)
  # ne = 3: semi-causal
  # ne = 3: 2D-RARMA(p,q), 2D-RARMA(p,q), 2D-RARMA(p,q)
  
  
  maxit1 = 1000
  
  p = max(ar) # AR order
  q = max(ma) # MA order
  n = dim(y)[1] # n is the number of rows
  k = dim(y)[2] # k is the number of columns 
  m = max(p,q,na.rm=T)
  
  if(ne == 1)
  {
    ynew = log(y)

  # inicializacao dos parametros alpha e phi (beta)
  if(any(is.na(ar)==F)) # se nao tem componente ar
  {
    p1 = (p+1)^2-1
    
    XX = c()
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx1 = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
        XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
      }
    }
    
    P = XX[,2:dim(XX)[2]]
    Y = as.matrix(XX[,1])

    Z <- cbind(rep(1,(n-m)*(k-m)),P)
    
    x = as.matrix(Z)
    ajuste = lm.fit(x, Y)
    mqo = c(ajuste$coef)
    
  }else{
    ynew = log(y)
    Z <- as.matrix(rep(1,(n-m)*(k-m)))
    
    Y = as.vector(ynew[(m+1):n,(m+1):k])
    
    x = as.matrix(Z)
    ajuste = lm.fit(x, Y)
    mqo = c(ajuste$coef)
  }
  
  
  ############
  
  ######### ARMA model
  if(any(is.na(ar)==F) && any(is.na(ma)==F))
  { 
    q1 = (q+1)^2-1
   
    reg = c(mqo, rep(0,q1)) # initializing the parameter values
    
    loglik = function(z) 
    {
      alpha = z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      
      eta = error = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          
          yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = as.matrix(yy[(length(yy)-1):1])
          
          eta[i,j]  = alpha + phi%*%y_new1 + theta%*%error_new
          error[i,j] = ynew[i,j]-eta[i,j]
        }
        
      }
      
      mu = as.vector(exp((t(eta[(m+1):n,(m+1):k]))))
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))
      
      ll = suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2))))
      
      sum(ll)
  
    } 
    
    escore = function(z) 
    {
      alpha = z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      
      eta = error = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          
          yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = as.matrix(yy[(length(yy)-1):1])
          
          eta[i,j]  = alpha + phi%*%y_new1 + theta%*%error_new
          error[i,j] = ynew[i,j]-eta[i,j]
        }
        
      }
      
      mu = as.vector(exp(t(eta[(m+1):n,(m+1):k])))
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))
  
      dmu = as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu)))
      
      mT = diag(mu)
      
      ###################################
      
      # dalpha
      deta.dalpha = matrix(0, nrow = n, ncol=k)
      
      for(i in (m+1):n)
      {
        for(j in (m+1):k)
        {
          xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
          deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
          deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
        }
        
      }
      
      a = as.vector(deta.dalpha[(m+1):n,(m+1):k])
      
      # dphi 
      P1 = rbind(rep(0,p1),P)
      deta.dphi = matrix(0, ncol=p1,nrow=dim(P1)[1])
      dsum.phi = matrix(0, ncol=p1,nrow=q1)
    
      for(i in m:dim(P1)[1])
      {
        if( i == m)
        {
          deta.dphi[i,]<- P1[i,] 
        }else{
          
          for(j in 1:q1)
          {
            dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
          }
          deta.dphi[i,]<- P1[i,] - apply(dsum.phi,2,mean)
        }
      }
      
      rP <- deta.dphi[-m,]
      
      #dtheta
      XX = c()
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx1 = as.vector(t(error[(i-ma):i,(j-ma):j]))
          XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        }
      }
      
      R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
      deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
      dsum.theta = matrix(0, ncol=q1,nrow=q1)
      
        for(i in m:dim(R)[1])
        {
          if( i == m)
          {
            deta.dtheta[i,]<- R[i,] 
          }else{
            for(j in 1:q1)
            {
              dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
            }
            deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
          }
        }
      rR <- deta.dtheta[-m,]
      
      ###################################
      
      Ualpha =  - a %*% mT %*% dmu
      Uphi =    - t(rP) %*% mT %*% dmu
      Utheta =  - t(rR) %*% mT %*% dmu

      rval = c(Ualpha,Uphi,Utheta)
    }
    
    if(any(is.na(ar))==F) names_phi = c(paste("phi",1:p1,sep=""))
    
    if(any(is.na(ma))==F) names_theta = (paste("theta",1:q1,sep=""))
    
    names_par <- c("alpha",names_phi,names_theta)
    
    opt = optim(reg, loglik, 
                 escore,
                 method = "BFGS", 
                 hessian = TRUE,
                 control = list(maxit = 400, abstol = 1e-6,
                                factr = 1e20))
    
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z = c()
    z$conv = opt$conv
    coef = opt$par
    names(coef) = names_par
    z$coeff = coef

    alpha = coef[1]
    phi = coef[2:(p1+1)]
    theta = coef[(p1+2):(p1+q1+1)]


    z$alpha = alpha
    z$phi = phi
    z$theta = theta
    
    etahat = errorhat = matrix(0, ncol=k,nrow=n)
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
        y_new1 = as.matrix(xx[(length(xx)-1):1])
        
        yy = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        errorhat_new = as.matrix(yy[(length(yy)-1):1])
        
        etahat[i,j]  = alpha + phi%*%y_new1 + theta%*%errorhat_new
        errorhat[i,j] = ynew[i,j]-etahat[i,j]
      }
      
    }

    z$fitted = exp(etahat[(m+1):n,(m+1):k])
    z$etahat = etahat[(m+1):n,(m+1):k]
    z$errorhat = errorhat[(m+1):n,(m+1):k]
  
    y1 = as.vector(t(y[(m+1):n,(m+1):k]))
    
    ###################################
    
    # dalpha
    deta.dalpha = matrix(0, nrow = n, ncol=k)
    
    for(i in (m+1):n)
    {
      for(j in (m+1):k)
      {
        xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
        deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
        deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
      }
      
    }
    
    a = as.vector(deta.dalpha[(m+1):n,(m+1):k])
    
    # dphi 
    P1 = rbind(rep(0,p1),P)
    deta.dphi = matrix(0, ncol=p1,nrow=dim(P1)[1])
    dsum.phi = matrix(0, ncol=p1,nrow=q1)
    
    for(i in m:dim(P1)[1])
    {
      if( i == m)
      {
        deta.dphi[i,]<- P1[i,] 
      }else{
        
        for(j in 1:q1)
        {
          dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
        }
        deta.dphi[i,]<- P1[i,] - apply(dsum.phi,2,mean)
      }
    }
    
    rP <- deta.dphi[-m,]
    
    ##dtheta
    XX = c()
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx1 = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
      }
    }
    
    R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
    deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
    dsum.theta = matrix(0, ncol=q1,nrow=q1)
    
    for(i in m:dim(R)[1])
    {
      if( i == m)
      {
        deta.dtheta[i,]<- R[i,] 
      }else{
        for(j in 1:q1)
        {
          dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
        }
        deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
      }
    }
    rR <- deta.dtheta[-m,]
    
    ###################################
    
    muhat2 = as.vector(exp(t(etahat[(m+1):n,(m+1):k])))
    W = diag(((4)/(muhat2^2))*(muhat2^2))

    Kaa = t(a) %*% W %*% a
    Kpa = t(rP) %*% W %*% a
    Kap = t(Kpa)
    Kta = t(rR) %*% W %*% a
    Kat = t(Kta)
    Kpp = t(rP) %*% W %*% rP
    Kpt = t(rP) %*% W %*% rR
    Ktp = t(Kpt)
    Ktt = t(rR) %*% W %*% rR

    K = rbind(
      cbind(Kaa,Kap,Kat),
      cbind(Kpa,Kpp,Kpt),
      cbind(Kta,Ktp,Ktt)
    )
    

  }  
  
  
  ##### only AR model
  if(any(is.na(ar)==F) && any(is.na(ma)==T))
  {
    
    p1 = (p+1)^2-1
    q1 = 0
    reg = c(mqo) # initializing the parameter values
 
    
    loglik = function(z) 
    {
      alpha = z[1]
      phi = z[2:(p1+1)] 
      
      eta = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          eta[i,j]  = alpha + phi%*%y_new1 
        }
        
      }
      
      mu = as.vector(exp((t(eta[(m+1):n,(m+1):k]))))
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))
      
      ll = suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2))))
      
      sum(ll)
      
    } 
    
    escore <- function(z) 
    {
      alpha = z[1]
      phi = z[2:(p1+1)] 
      
      eta = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          eta[i,j]  = alpha + phi%*%y_new1 
        }
        
      }
      
      mu = as.vector(exp(t(eta[(m+1):n,(m+1):k])))
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))
      
      dmu = as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu)))
      
      mT = diag(mu)
      a = as.vector(rep(1,(n-m)*(k-m)))
      
      Ualpha =  - t(a) %*% mT %*% dmu
      Uphi =    - t(P) %*% mT %*% dmu
      
      rval = c(Ualpha,Uphi)
    }
    
    if(any(is.na(ar))==F) names_phi = c(paste("phi",1:p1,sep=""))
    
    names_par <- c("alpha",names_phi)
    
    opt = optim(reg, loglik, 
                escore,
                method = "BFGS", 
                hessian = TRUE,
                control = list(maxit = 400, abstol = 1e-6,
                               factr = 1e20))
    
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z = c()
    z$conv = opt$conv
    coef = (opt$par)
    names(coef) = names_par
    z$coeff = coef
    
    alpha = coef[1]
    phi = coef[2:(p1+1)]

    z$alpha = alpha
    z$phi = phi
    
    etahat = matrix(0, ncol=k,nrow=n)
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
        y_new1 = as.matrix(xx[(length(xx)-1):1])
        etahat[i,j]  = alpha + phi%*%y_new1 
  
      }
      
    }
    
    
    z$fitted = exp(etahat[(m+1):n,(m+1):k])
    z$etahat = etahat[(m+1):n,(m+1):k]
    
    y1 = as.vector(t(y[(m+1):n,(m+1):k]))
    
    muhat2 = as.vector(exp(t(etahat[(m+1):n,(m+1):k])))
    W = diag(((4)/(muhat2^2))*(muhat2^2))
    
    a = as.vector(rep(1,(n-m)*(k-m)))
    
    Kaa = t(a) %*% W %*% a
    Kpa = t(P) %*% W %*% a
    Kap = t(Kpa)
    Kpp = t(P) %*% W %*% P
    
    K = rbind(
      cbind(Kaa,Kap),
      cbind(Kpa,Kpp)
    )

    
  }
  
  ######### MA model
  if(any(is.na(ar)==T) && any(is.na(ma)==F))
  { 
    p1 = 0
    q1 = (q+1)^2-1
    reg = c(mqo[1],rep(0,q1)) # initializing the parameter values
    
    loglik = function(z) 
    {
      alpha = z[1]
      theta = z[2:(q1+1)]
      
      eta = error = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = as.matrix(yy[(length(yy)-1):1])
          
          eta[i,j]  = alpha + theta%*%error_new
          error[i,j] = ynew[i,j]-eta[i,j]
        }
        
      }
      
      mu = as.vector(exp((t(eta[(m+1):n,(m+1):k]))))
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))
      
      ll = suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2))))
      
      sum(ll)
      
    } 
    
    escore = function(z) 
    {
      alpha = z[1]
      theta = z[2:(q1+1)]
      
      eta = error = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          
          yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = as.matrix(yy[(length(yy)-1):1])
          
          eta[i,j]  = alpha + theta%*%error_new
          error[i,j] = ynew[i,j]-eta[i,j]
        }
        
      }
      
      mu = as.vector(exp(t(eta[(m+1):n,(m+1):k])))
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))
      
      dmu = as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu)))
      
      mT = diag(mu)
      
      ###################################
      
      # dalpha
      deta.dalpha = matrix(0, nrow = n, ncol=k)
      
      for(i in (m+1):n)
      {
        for(j in (m+1):k)
        {
          xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
          deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
          deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
        }
        
      }
      
      a = as.vector(deta.dalpha[(m+1):n,(m+1):k])

      #dtheta
      XX = c()
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx1 = as.vector(t(error[(i-ma):i,(j-ma):j]))
          XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        }
      }
      
      R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
      deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
      dsum.theta = matrix(0, ncol=q1,nrow=q1)
      
      for(i in m:dim(R)[1])
      {
        if( i == m)
        {
          deta.dtheta[i,]<- R[i,] 
        }else{
          for(j in 1:q1)
          {
            dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
          }
          
          deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
        }
      }
      rR <- deta.dtheta[-m,]
      
      ###################################
      
      Ualpha =  - a %*% mT %*% dmu
      Utheta =  - t(rR) %*% mT %*% dmu
      
      rval = c(Ualpha,Utheta)
    }
    
    if(any(is.na(ma))==F) names_theta = (paste("theta",1:q1,sep=""))
    
    names_par <- c("alpha",names_theta)
    
    opt = optim(reg, loglik, 
                escore,
                method = "BFGS", 
                hessian = TRUE,
                control = list(maxit = 400, abstol = 1e-6,
                               factr = 1e20))
    
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z = c()
    z$conv = opt$conv
    coef = opt$par
    names(coef) = names_par
    z$coeff = coef
    
    alpha = coef[1]
    theta = coef[2:(q1+1)]
    
    z$alpha = alpha
    z$theta = theta
    
    etahat = errorhat = matrix(0, ncol=k,nrow=n)
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
      
        yy = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        errorhat_new = as.matrix(yy[(length(yy)-1):1])
        
        etahat[i,j]  = alpha + theta%*%errorhat_new
        errorhat[i,j] = ynew[i,j]-etahat[i,j]
      }
      
    }
    
    
    z$fitted = exp(etahat[(m+1):n,(m+1):k])
    z$etahat = etahat[(m+1):n,(m+1):k]
    z$errorhat = errorhat[(m+1):n,(m+1):k]
    
    y1 = as.vector(t(y[(m+1):n,(m+1):k]))
    
    ###################################
    
    # dalpha
    deta.dalpha = matrix(0, nrow = n, ncol=k)
    
    for(i in (m+1):n)
    {
      for(j in (m+1):k)
      {
        xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
        deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
        deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
      }
      
    }
    
    a = as.vector(deta.dalpha[(m+1):n,(m+1):k])
    
    #dtheta
    XX = c()
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx1 = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
      }
    }
    
    R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
    deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
    dsum.theta = matrix(0, ncol=q1,nrow=q1)
    
    for(i in m:dim(R)[1])
    {
      if( i == m)
      {
        deta.dtheta[i,]<- R[i,] 
      }else{
        for(j in 1:q1)
        {
          dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
        }
        deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
      }
    }
    rR <- deta.dtheta[-m,]

    
    ###################################
    
    muhat2 = as.vector(exp(t(etahat[(m+1):n,(m+1):k])))
    W = diag(((4)/(muhat2^2))*(muhat2^2))
    
    Kaa = t(a) %*% W %*% a
    Kta = t(rR) %*% W %*% a
    Kat = t(Kta)
    Ktt = t(rR) %*% W %*% rR
    
    K = rbind(
      cbind(Kaa,Kat),
      cbind(Kta,Ktt)
    )
    
  }
  
  z$image = y
  y1 = y[(m+1):n,(m+1):k]
  z$fitted = exp(z$etahat)

  ##############################################

  # residuals

  z$resid = qnorm(pr(y1,z$fitted)) # quantile residuals
  # 
  
  ############################################

  vcov = chol2inv(chol(K))
  z$vcov = vcov

  stderror = sqrt(diag(vcov))
  z$stderror = stderror

  z$zstat = abs(z$coef/stderror)
  z$pvalues = 2*(1 - pnorm(z$zstat) )

  model_presentation = cbind(round(z$coef,4),round(z$stderror,4),
                             round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)=c("Estimate","Std. Error","z value","Pr(>|z|)")

  z$model = model_presentation

  # 
   return(z)
  
  }
  
  if(ne == 2)
  {
    ynew = log(y)
    # inicializacao dos parametros alpha e phi (beta)
    if(any(is.na(ar)==F)) # se nao tem componente ar
    {
      p1 = (2*p+1)^2-1
      
      XX = c()
      
      for (i in (m+1):(n-m))
      {
        for (j in (m+1):(k-m))
        {
          xx1 = as.vector(t(ynew[(i-ar):(i+ar),(j-ar):(j+ar)]))
          XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        }
      }
      
      P = XX[,2:dim(XX)[2]]
      Y = as.matrix(XX[,1])
      
      Z <- cbind(rep(1,(n-(2*m))*(k-(2*m))),P)
      
      x = as.matrix(Z)
      ajuste = lm.fit(x, Y)
      mqo = c(ajuste$coef)
      
    }else{
      ynew = log(y)
      Z <- as.matrix(rep(1,(n-(2*m))*(k-(2*m))))
      
      Y = as.vector(ynew[(m+1):(n-m),(m+1):(k-m)])
      
      x = as.matrix(Z)
      ajuste = lm.fit(x, Y)
      mqo = c(ajuste$coef)
    }
    
    
    ############
    
    ######### ARMA model
    if(any(is.na(ar)==F) && any(is.na(ma)==F))
    { 
      q1 = (2*q+1)^2-1
      
      reg = c(mqo, rep(0,q1)) # initializing the parameter values
      
      loglik = function(z) 
      {
        alpha = z[1]
        phi = z[2:(p1+1)] 
        theta = z[(p1+2):(p1+q1+1)]
        
        eta = error = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):(n-m))
        {
          for (j in (m+1):(k-m))
          {
            xx = as.vector(t(ynew[(i-ar):(i+ar),(j-ar):(j+ar)]))
            y_new1 = as.matrix(xx[(length(xx)-1):1])
            
            yy = as.vector(t(error[(i-ma):(i+ma),(j-ma):(j+ma)]))
            error_new = as.matrix(yy[(length(yy)-1):1])
            
            eta[i,j]  = alpha + phi%*%y_new1 + theta%*%error_new
            error[i,j] = ynew[i,j]-eta[i,j]
          }
          
        }
        
        mu = as.vector(exp((t(eta[(m+1):(n-m),(m+1):(k-m)]))))
        y1 = as.vector(t(y[(m+1):(n-m),(m+1):(k-m)]))
        
        ll = suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2))))
        
        sum(ll)
        
      } 
      
      escore = function(z) 
      {
        alpha = z[1]
        phi = z[2:(p1+1)] 
        theta = z[(p1+2):(p1+q1+1)]
        
        eta = error = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):(n-m))
        {
          for (j in (m+1):(k-m))
          {
            xx = as.vector(t(ynew[(i-ar):(i+ar),(j-ar):(j+ar)]))
            y_new1 = as.matrix(xx[(length(xx)-1):1])
            
            yy = as.vector(t(error[(i-ma):(i+ma),(j-ma):(j+ma)]))
            error_new = as.matrix(yy[(length(yy)-1):1])
            
            eta[i,j]  = alpha + phi%*%y_new1 + theta%*%error_new
            error[i,j] = ynew[i,j]-eta[i,j]
          }
          
        }
        
        mu = as.vector(exp(t(eta[(m+1):(n-m),(m+1):(k-m)])))
        y1 = as.vector(t(y[(m+1):(n-m),(m+1):(k-m)]))
        
        dmu = as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu)))
        
        mT = diag(mu)
        
        ###################################
        
        # dalpha
        deta.dalpha = matrix(0, nrow = n, ncol=k)
        
        for(i in (m+1):(n-m))
        {
          for(j in (m+1):(k-m))
          {
            xx = as.vector(t(deta.dalpha[(i-ma):(i+ma),(j-ma):(j+ma)]))
            deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
            deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
          }
          
        }
        
        a = as.vector(deta.dalpha[(m+1):(n-m),(m+1):(k-m)])
        
        # dphi 
        P1 = rbind(rep(0,p1),P)
        deta.dphi = matrix(0, ncol=p1,nrow=dim(P1)[1])
        dsum.phi = matrix(0, ncol=p1,nrow=q1)
        
        for(i in m:dim(P1)[1])
        {
          if( i == m)
          {
            deta.dphi[i,]<- P1[i,] 
          }else{
            
            for(j in 1:q1)
            {
              dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
            }
            deta.dphi[i,]<- P1[i,] - apply(dsum.phi,2,mean)
          }
        }
        
        rP <- deta.dphi[-m,]
        
        #dtheta
        XX = c()
        
        for (i in (m+1):(n-m))
        {
          for (j in (m+1):(k-m))
          {
            xx1 = as.vector(t(error[(i-ma):(i+ma),(j-ma):(j+ma)]))
            XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
          }
        }
        
        R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
        deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
        dsum.theta = matrix(0, ncol=q1,nrow=q1)
        
        for(i in m:dim(R)[1])
        {
          if( i == m)
          {
            deta.dtheta[i,]<- R[i,] 
          }else{
            for(j in 1:q1)
            {
              dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
            }
            deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
          }
        }
        rR <- deta.dtheta[-m,]
        
        ###################################
        
        Ualpha =  - a %*% mT %*% dmu
        Uphi =    - t(rP) %*% mT %*% dmu
        Utheta =  - t(rR) %*% mT %*% dmu
        
        rval = c(Ualpha,Uphi,Utheta)
      }
      
      if(any(is.na(ar))==F) names_phi = c(paste("phi",1:p1,sep=""))
      
      if(any(is.na(ma))==F) names_theta = (paste("theta",1:q1,sep=""))
      
      names_par <- c("alpha",names_phi,names_theta)
      
      opt = optim(reg, loglik, 
                  escore,
                  method = "BFGS", 
                  hessian = TRUE,
                  control = list(maxit = 400, abstol = 1e-6,
                                 factr = 1e20))
      
      
      if (opt$conv != 0)
        warning("FUNCTION DID NOT CONVERGE!")
      
      z = c()
      z$conv = opt$conv
      coef = opt$par
      names(coef) = names_par
      z$coeff = coef
      
      alpha = coef[1]
      phi = coef[2:(p1+1)]
      theta = coef[(p1+2):(p1+q1+1)]
      
      
      z$alpha = alpha
      z$phi = phi
      z$theta = theta
      
      etahat = errorhat = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):(n-m))
      {
        for (j in (m+1):(k-m))
        {
          xx = as.vector(t(ynew[(i-ar):(i+ar),(j-ar):(j+ar)]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          
          yy = as.vector(t(errorhat[(i-ma):(i+ma),(j-ma):(j+ma)]))
          errorhat_new = as.matrix(yy[(length(yy)-1):1])
          
          etahat[i,j]  = alpha + phi%*%y_new1 + theta%*%errorhat_new
          errorhat[i,j] = ynew[i,j]-etahat[i,j]
        }
        
      }
      
      z$fitted = exp(etahat[(m+1):(n-m),(m+1):(k-m)])
      z$etahat = etahat[(m+1):(n-m),(m+1):(k-m)]
      z$errorhat = errorhat[(m+1):n,(m+1):(k-m)]
      
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))
      
      ###################################
      
      # dalpha
      deta.dalpha = matrix(0, nrow = n, ncol=k)
      
      for(i in (m+1):(n-m))
      {
        for(j in (m+1):(k-m))
        {
          xx = as.vector(t(deta.dalpha[(i-ma):(i+ma),(j-ma):(j+ma)]))
          deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
          deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
        }
        
      }
      
      a = as.vector(deta.dalpha[(m+1):(n-m),(m+1):(k-m)])
      
      # dphi 
      P1 = rbind(rep(0,p1),P)
      deta.dphi = matrix(0, ncol=p1,nrow=dim(P1)[1])
      dsum.phi = matrix(0, ncol=p1,nrow=q1)
      
      for(i in m:dim(P1)[1])
      {
        if( i == m)
        {
          deta.dphi[i,]<- P1[i,] 
        }else{
          
          for(j in 1:q1)
          {
            dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
          }
          deta.dphi[i,]<- P1[i,] - apply(dsum.phi,2,mean)
        }
      }
      
      rP <- deta.dphi[-m,]
      
      ##dtheta
      XX = c()
      
      for (i in (m+1):(n-m))
      {
        for (j in (m+1):(k-m))
        {
          xx1 = as.vector(t(errorhat[(i-ma):(i+ma),(j-ma):(j+ma)]))
          XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        }
      }
      
      R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
      deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
      dsum.theta = matrix(0, ncol=q1,nrow=q1)
      
      for(i in m:dim(R)[1])
      {
        if( i == m)
        {
          deta.dtheta[i,]<- R[i,] 
        }else{
          for(j in 1:q1)
          {
            dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
          }
          deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
        }
      }
      rR <- deta.dtheta[-m,]
      
      ###################################
      
      muhat2 = as.vector(exp(t(etahat[(m+1):(n-m),(m+1):(k-m)])))
      W = diag(((4)/(muhat2^2))*(muhat2^2))
      
      Kaa = t(a) %*% W %*% a
      Kpa = t(rP) %*% W %*% a
      Kap = t(Kpa)
      Kta = t(rR) %*% W %*% a
      Kat = t(Kta)
      Kpp = t(rP) %*% W %*% rP
      Kpt = t(rP) %*% W %*% rR
      Ktp = t(Kpt)
      Ktt = t(rR) %*% W %*% rR
      
      K = rbind(
        cbind(Kaa,Kap,Kat),
        cbind(Kpa,Kpp,Kpt),
        cbind(Kta,Ktp,Ktt)
      )
      
      
    }  
    
    
    ##### only AR model
    if(any(is.na(ar)==F) && any(is.na(ma)==T))
    {
      
      p1 = (2*p+1)^2-1
      q1 = 0
      reg = c(mqo) # initializing the parameter values
      
      
      loglik = function(z) 
      {
        alpha = z[1]
        phi = z[2:(p1+1)] 
        
        eta = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):(n-m))
        {
          for (j in (m+1):(k-m))
          {
            xx = as.vector(t(ynew[(i-ar):(i+ar),(j-ar):(j+ar)]))
            y_new1 = as.matrix(xx[(length(xx)-1):1])
            
            eta[i,j]  = alpha + phi%*%y_new1 
          }
          
        }
        
        mu = as.vector(exp((t(eta[(m+1):(n-m),(m+1):(k-m)]))))
        y1 = as.vector(t(y[(m+1):(n-m),(m+1):(k-m)]))
        
        ll = suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2))))
        
        sum(ll)
        
      } 
      
      escore <- function(z) 
      {
        alpha = z[1]
        phi = z[2:(p1+1)] 
        
        eta = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):(n-m))
        {
          for (j in (m+1):(k-m))
          {
            xx = as.vector(t(ynew[(i-ar):(i+ar),(j-ar):(j+ar)]))
            y_new1 = as.matrix(xx[(length(xx)-1):1])
            
            eta[i,j]  = alpha + phi%*%y_new1 
          }
          
        }
        
        mu = as.vector(exp(t(eta[(m+1):(n-m),(m+1):(k-m)])))
        y1 = as.vector(t(y[(m+1):(n-m),(m+1):(k-m)]))
        
        dmu = as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu)))
        
        mT = diag(mu)
        a = as.vector(rep(1,(n-2*m)*(k-2*m)))
        
        Ualpha =  - t(a) %*% mT %*% dmu
        Uphi =    - t(P) %*% mT %*% dmu
        
        rval = c(Ualpha,Uphi)
      }
      
      if(any(is.na(ar))==F) names_phi = c(paste("phi",1:p1,sep=""))
      
      names_par <- c("alpha",names_phi)
      
      opt = optim(reg, loglik, 
                  escore,
                  method = "BFGS", 
                  hessian = TRUE,
                  control = list(maxit = 400, abstol = 1e-6,
                                 factr = 1e20))
      
      
      if (opt$conv != 0)
        warning("FUNCTION DID NOT CONVERGE!")
      
      z = c()
      z$conv = opt$conv
      coef = (opt$par)
      names(coef) = names_par
      z$coeff = coef
      
      alpha = coef[1]
      phi = coef[2:(p1+1)]
      
      z$alpha = alpha
      z$phi = phi
      
      etahat = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):(n-m))
      {
        for (j in (m+1):(k-m))
        {
          xx = as.vector(t(ynew[(i-ar):(i+ar),(j-ar):(j+ar)]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          
          etahat[i,j]  = alpha + phi%*%y_new1 
          
        }
        
      }
      
      
      z$fitted = exp(etahat[(m+1):(n-m),(m+1):(k-m)])
      z$etahat = etahat[(m+1):(n-m),(m+1):(k-m)]
      
      y1 = as.vector(t(y[(m+1):(n-m),(m+1):(k-m)]))
      
      muhat2 = as.vector(exp(t(etahat[(m+1):(n-m),(m+1):(k-m)])))
      W = diag(((4)/(muhat2^2))*(muhat2^2))
      
      a = as.vector(rep(1,(n-2*m)*(k-2*m)))
      
      Kaa = t(a) %*% W %*% a
      Kpa = t(P) %*% W %*% a
      Kap = t(Kpa)
      Kpp = t(P) %*% W %*% P
      
      K = rbind(
        cbind(Kaa,Kap),
        cbind(Kpa,Kpp)
      )
      
      
    }
    
    ######### MA model
    if(any(is.na(ar)==T) && any(is.na(ma)==F))
    { 
      p1 = 0
      q1 = (2*q+1)^2-1
      reg = c(mqo[1],rep(0,q1)) # initializing the parameter values
      
      loglik = function(z) 
      {
        alpha = z[1]
        theta = z[2:(q1+1)]
        
        eta = error = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):(n-m))
        {
          for (j in (m+1):(k-m))
          {
            yy = as.vector(t(error[(i-ma):(i+ma),(j-ma):(j+ma)]))
            error_new = as.matrix(yy[(length(yy)-1):1])
            
            eta[i,j]  = alpha + theta%*%error_new
            error[i,j] = ynew[i,j]-eta[i,j]
          }
          
        }
        
        mu = as.vector(exp((t(eta[(m+1):(n-m),(m+1):(k-m)]))))
        y1 = as.vector(t(y[(m+1):(n-m),(m+1):(k-m)]))
        
        ll = suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2))))
        
        sum(ll)
        
      } 
      
      escore = function(z) 
      {
        alpha = z[1]
        theta = z[2:(q1+1)]
        
        eta = error = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):(n-m))
        {
          for (j in (m+1):(k-m))
          {
            
            yy = as.vector(t(error[(i-ma):(i+ma),(j-ma):(j+ma)]))
            error_new = as.matrix(yy[(length(yy)-1):1])
            
            eta[i,j]  = alpha + theta%*%error_new
            error[i,j] = ynew[i,j]-eta[i,j]
          }
          
        }
        
        mu = as.vector(exp(t(eta[(m+1):(n-m),(m+1):(k-m)])))
        y1 = as.vector(t(y[(m+1):(n-m),(m+1):(k-m)]))
        
        dmu = as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu)))
        
        mT = diag(mu)
        
        ###################################
        
        # dalpha
        deta.dalpha = matrix(0, nrow = n, ncol=k)
        
        for(i in (m+1):(n-m))
        {
          for(j in (m+1):(k-m))
          {
            xx = as.vector(t(deta.dalpha[(i-ma):(i+ma),(j-ma):(j+ma)]))
            deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
            deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
          }
          
        }
        
        a = as.vector(deta.dalpha[(m+1):(n-m),(m+1):(k-m)])
        
        #dtheta
        XX = c()
        
        for (i in (m+1):(n-m))
        {
          for (j in (m+1):(k-m))
          {
            xx1 = as.vector(t(error[(i-ma):(i+ma),(j-ma):(j+ma)]))
            XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
          }
        }
        
        R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
        deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
        dsum.theta = matrix(0, ncol=q1,nrow=q1)
        
        for(i in m:dim(R)[1])
        {
          if( i == m)
          {
            deta.dtheta[i,]<- R[i,] 
          }else{
            for(j in 1:q1)
            {
              dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
            }
            
            deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
          }
        }
        rR <- deta.dtheta[-m,]
        
        ###################################
        
        Ualpha =  - a %*% mT %*% dmu
        Utheta =  - t(rR) %*% mT %*% dmu
        
        rval = c(Ualpha,Utheta)
      }
      
      if(any(is.na(ma))==F) names_theta = (paste("theta",1:q1,sep=""))
      
      names_par <- c("alpha",names_theta)
      
      opt = optim(reg, loglik, 
                  escore,
                  method = "BFGS", 
                  hessian = TRUE,
                  control = list(maxit = 400, abstol = 1e-6,
                                 factr = 1e20))
      
      if (opt$conv != 0)
        warning("FUNCTION DID NOT CONVERGE!")
      
      z = c()
      z$conv = opt$conv
      coef = opt$par
      names(coef) = names_par
      z$coeff = coef
      
      alpha = coef[1]
      theta = coef[2:(q1+1)]
      
      z$alpha = alpha
      z$theta = theta
      
      etahat = errorhat = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):(n-m))
      {
        for (j in (m+1):(k-m))
        {
          
          yy = as.vector(t(errorhat[(i-ma):(i+ma),(j-ma):(j+ma)]))
          errorhat_new = as.matrix(yy[(length(yy)-1):1])
          
          etahat[i,j]  = alpha + theta%*%errorhat_new
          errorhat[i,j] = ynew[i,j]-etahat[i,j]
        }
        
      }
      
      
      z$fitted = exp(etahat[(m+1):(n-m),(m+1):(k-m)])
      z$etahat = etahat[(m+1):(n-m),(m+1):(k-m)]
      z$errorhat = errorhat[(m+1):(n-m),(m+1):(k-m)]
      
      y1 = as.vector(t(y[(m+1):(n-m),(m+1):(k-m)]))
      
      ###################################
      
      # dalpha
      deta.dalpha = matrix(0, nrow = n, ncol=k)
      
      for(i in (m+1):(n-m))
      {
        for(j in (m+1):(k-m))
        {
          xx = as.vector(t(deta.dalpha[(i-ma):(i+ma),(j-ma):(j+ma)]))
          deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
          deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
        }
        
      }
      
      a = as.vector(deta.dalpha[(m+1):(n-m),(m+1):(k-m)])
      
      #dtheta
      XX = c()
      
      for (i in (m+1):(n-m))
      {
        for (j in (m+1):(k-m))
        {
          xx1 = as.vector(t(errorhat[(i-ma):(i+ma),(j-ma):(j+ma)]))
          XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        }
      }
      
      R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
      deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
      dsum.theta = matrix(0, ncol=q1,nrow=q1)
      
      for(i in m:dim(R)[1])
      {
        if( i == m)
        {
          deta.dtheta[i,]<- R[i,] 
        }else{
          for(j in 1:q1)
          {
            dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
          }
          deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
        }
      }
      rR <- deta.dtheta[-m,]
      
      
      ###################################
      
      muhat2 = as.vector(exp(t(etahat[(m+1):(n-m),(m+1):(k-m)])))
      W = diag(((4)/(muhat2^2))*(muhat2^2))
      
      Kaa = t(a) %*% W %*% a
      Kta = t(rR) %*% W %*% a
      Kat = t(Kta)
      Ktt = t(rR) %*% W %*% rR
      
      K = rbind(
        cbind(Kaa,Kat),
        cbind(Kta,Ktt)
      )
      
    }
    
    z$image = y
    y1 = y[(m+1):(n-m),(m+1):(k-m)]
    z$fitted = exp(z$etahat)
    
    ##############################################
    
    # residuals
    
    z$resid = qnorm(pr(y1,z$fitted)) # quantile residuals
    # 
    
    ############################################
    
    vcov = chol2inv(chol(K))
    z$vcov = vcov
    
    stderror = sqrt(diag(vcov))
    z$stderror = stderror
    
    z$zstat = abs(z$coef/stderror)
    z$pvalues = 2*(1 - pnorm(z$zstat) )
    
    model_presentation = cbind(round(z$coef,4),round(z$stderror,4),
                               round(z$zstat,4),round(z$pvalues,4))
    colnames(model_presentation)=c("Estimate","Std. Error","z value","Pr(>|z|)")
    
    z$model = model_presentation
    
    # 
    return(z)
    
  }
  
  if(ne == 3)
  {
    ynew = log(y)
    # inicializacao dos parametros alpha e phi (beta)
    if(any(is.na(ar)==F)) # se nao tem componente ar
    {
      p1 = 2*(p+1)^2-2-ar

      XX = c()
      
      for (i in (m+1):n)
      {
        for (j in (m+1):(k-m))
        {
          xx1 = as.vector(t(ynew[(i-ar):i,(j-ar):(j+ar)]))
          XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        }
      }
      
      P = XX[,2:dim(XX)[2]]
      Y = as.matrix(XX[,1])
      
      Z <- cbind(rep(1,(n-m)*(k-2*m)),P)
      
      x = as.matrix(Z)
      ajuste = lm.fit(x, Y)
      mqo = c(ajuste$coef)
      
    }else{
      ynew = log(y)
      Z <- as.matrix(rep(1,(n-m)*(k-2*m)))
      
      Y = as.vector(ynew[(m+1):n,(m+1):(k-m)])
      
      x = as.matrix(Z)
      ajuste = lm.fit(x, Y)
      mqo = c(ajuste$coef)
    }
    
    
    ############
    
    ######### ARMA model
    if(any(is.na(ar)==F) && any(is.na(ma)==F))
    { 
      q1 = 2*(q+1)^2-2-ma
      
      reg = c(mqo, rep(0,q1)) # initializing the parameter values
      
      loglik = function(z) 
      {
        alpha = z[1]
        phi = z[2:(p1+1)] 
        theta = z[(p1+2):(p1+q1+1)]
        
        eta = error = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):n)
        {
          for (j in (m+1):(k-m))
          {
            xx = as.vector(t(ynew[(i-ar):i,(j-ar):(j+ar)]))
            y_new1 = as.matrix(xx[(length(xx)-1):1])
            
            yy = as.vector(t(error[(i-ma):i,(j-ma):(j+ma)]))
            error_new = as.matrix(yy[(length(yy)-1):1])
            
            eta[i,j]  = alpha + phi%*%y_new1 + theta%*%error_new
            error[i,j] = ynew[i,j]-eta[i,j]
          }
          
        }
        
        mu = as.vector(exp((t(eta[(m+1):n,(m+1):(k-m)]))))
        y1 = as.vector(t(y[(m+1):n,(m+1):(k-m)]))
        
        ll = suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2))))
        
        sum(ll)
        
      } 
      
      escore = function(z) 
      {
        alpha = z[1]
        phi = z[2:(p1+1)] 
        theta = z[(p1+2):(p1+q1+1)]
        
        eta = error = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):n)
        {
          for (j in (m+1):(k-m))
          {
            xx = as.vector(t(ynew[(i-ar):i,(j-ar):(j+ar)]))
            y_new1 = as.matrix(xx[(length(xx)-1):1])
            
            yy = as.vector(t(error[(i-ma):i,(j-ma):(j+ma)]))
            error_new = as.matrix(yy[(length(yy)-1):1])
            
            eta[i,j]  = alpha + phi%*%y_new1 + theta%*%error_new
            error[i,j] = ynew[i,j]-eta[i,j]
          }
          
        }
        
        mu = as.vector(exp(t(eta[(m+1):n,(m+1):(k-m)])))
        y1 = as.vector(t(y[(m+1):n,(m+1):(k-m)]))
        
        dmu = as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu)))
        
        mT = diag(mu)
        
        ###################################
        
        # dalpha
        deta.dalpha = matrix(0, nrow = n, ncol=k)
        
        for(i in (m+1):n)
        {
          for(j in (m+1):(k-m))
          {
            xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):(j+ma)]))
            deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
            deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
          }
          
        }
        
        a = as.vector(deta.dalpha[(m+1):n,(m+1):(k-m)])
        
        # dphi 
        P1 = rbind(rep(0,p1),P)
        deta.dphi = matrix(0, ncol=p1,nrow=dim(P1)[1])
        dsum.phi = matrix(0, ncol=p1,nrow=q1)
        
        for(i in m:dim(P1)[1])
        {
          if( i == m)
          {
            deta.dphi[i,]<- P1[i,] 
          }else{
            
            for(j in 1:q1)
            {
              dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
            }
            deta.dphi[i,]<- P1[i,] - apply(dsum.phi,2,mean)
          }
        }
        
        rP <- deta.dphi[-m,]
        
        #dtheta
        XX = c()
        
        for (i in (m+1):n)
        {
          for (j in (m+1):(k-m))
          {
            xx1 = as.vector(t(error[(i-ma):i,(j-ma):(j+ma)]))
            XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
          }
        }
        
        R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
        deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
        dsum.theta = matrix(0, ncol=q1,nrow=q1)
        
        for(i in m:dim(R)[1])
        {
          if( i == m)
          {
            deta.dtheta[i,]<- R[i,] 
          }else{
            for(j in 1:q1)
            {
              dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
            }
            deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
          }
        }
        rR <- deta.dtheta[-m,]
        
        ###################################
        
        Ualpha =  - a %*% mT %*% dmu
        Uphi =    - t(rP) %*% mT %*% dmu
        Utheta =  - t(rR) %*% mT %*% dmu
        
        rval = c(Ualpha,Uphi,Utheta)
      }
      
      if(any(is.na(ar))==F) names_phi = c(paste("phi",1:p1,sep=""))
      
      if(any(is.na(ma))==F) names_theta = (paste("theta",1:q1,sep=""))
      
      names_par <- c("alpha",names_phi,names_theta)
      
      opt = optim(reg, loglik, 
                  escore,
                  method = "BFGS", 
                  hessian = TRUE,
                  control = list(maxit = 400, abstol = 1e-6,
                                 factr = 1e20))
      
      
      if (opt$conv != 0)
        warning("FUNCTION DID NOT CONVERGE!")
      
      z = c()
      z$conv = opt$conv
      coef = opt$par
      names(coef) = names_par
      z$coeff = coef
      
      alpha = coef[1]
      phi = coef[2:(p1+1)]
      theta = coef[(p1+2):(p1+q1+1)]
      
      
      z$alpha = alpha
      z$phi = phi
      z$theta = theta
      
      etahat = errorhat = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):(k-m))
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):(j+ar)]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          
          yy = as.vector(t(errorhat[(i-ma):i,(j-ma):(j+ma)]))
          errorhat_new = as.matrix(yy[(length(yy)-1):1])
          
          etahat[i,j]  = alpha + phi%*%y_new1 + theta%*%errorhat_new
          errorhat[i,j] = ynew[i,j]-etahat[i,j]
        }
        
      }
      
      z$fitted = exp(etahat[(m+1):n,(m+1):(k-m)])
      z$etahat = etahat[(m+1):n,(m+1):(k-m)]
      z$errorhat = errorhat[(m+1):n,(m+1):(k-m)]
      
      y1 = as.vector(t(y[(m+1):n,(m+1):(k-m)]))
      
      ###################################
      
      # dalpha
      deta.dalpha = matrix(0, nrow = n, ncol=k)
      
      for(i in (m+1):n)
      {
        for(j in (m+1):(k-m))
        {
          xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):(j+ma)]))
          deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
          deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
        }
        
      }
      
      a = as.vector(deta.dalpha[(m+1):n,(m+1):(k-m)])
      
      # dphi 
      P1 = rbind(rep(0,p1),P)
      deta.dphi = matrix(0, ncol=p1,nrow=dim(P1)[1])
      dsum.phi = matrix(0, ncol=p1,nrow=q1)
      
      for(i in m:dim(P1)[1])
      {
        if( i == m)
        {
          deta.dphi[i,]<- P1[i,] 
        }else{
          
          for(j in 1:q1)
          {
            dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
          }
          deta.dphi[i,]<- P1[i,] - apply(dsum.phi,2,mean)
        }
      }
      
      rP <- deta.dphi[-m,]
      
      ##dtheta
      XX = c()
      
      for (i in (m+1):n)
      {
        for (j in (m+1):(k-m))
        {
          xx1 = as.vector(t(errorhat[(i-ma):i,(j-ma):(j+ma)]))
          XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        }
      }
      
      R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
      deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
      dsum.theta = matrix(0, ncol=q1,nrow=q1)
      
      for(i in m:dim(R)[1])
      {
        if( i == m)
        {
          deta.dtheta[i,]<- R[i,] 
        }else{
          for(j in 1:q1)
          {
            dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
          }
          deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
        }
      }
      rR <- deta.dtheta[-m,]
      
      ###################################
      
      muhat2 = as.vector(exp(t(etahat[(m+1):n,(m+1):(k-m)])))
      W = diag(((4)/(muhat2^2))*(muhat2^2))
      
      Kaa = t(a) %*% W %*% a
      Kpa = t(rP) %*% W %*% a
      Kap = t(Kpa)
      Kta = t(rR) %*% W %*% a
      Kat = t(Kta)
      Kpp = t(rP) %*% W %*% rP
      Kpt = t(rP) %*% W %*% rR
      Ktp = t(Kpt)
      Ktt = t(rR) %*% W %*% rR
      
      K = rbind(
        cbind(Kaa,Kap,Kat),
        cbind(Kpa,Kpp,Kpt),
        cbind(Kta,Ktp,Ktt)
      )
      
      
    }  
    
    
    ##### only AR model
    if(any(is.na(ar)==F) && any(is.na(ma)==T))
    {
      
      p1 = 2*(p+1)^2-2-ar
      q1 = 0
      reg = c(mqo) # initializing the parameter values
      
      
      loglik = function(z) 
      {
        alpha = z[1]
        phi = z[2:(p1+1)] 
        
        eta = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):n)
        {
          for (j in (m+1):(k-m))
          {
            xx = as.vector(t(ynew[(i-ar):i,(j-ar):(j+ar)]))
            y_new1 = as.matrix(xx[(length(xx)-1):1])
            
            eta[i,j]  = alpha + phi%*%y_new1 
          }
          
        }
        
        mu = as.vector(exp((t(eta[(m+1):n,(m+1):(k-m)]))))
        y1 = as.vector(t(y[(m+1):n,(m+1):k]))
        
        ll = suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2))))
        
        sum(ll)
        
      } 
      
      escore <- function(z) 
      {
        alpha = z[1]
        phi = z[2:(p1+1)] 
        
        eta = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):n)
        {
          for (j in (m+1):(k-m))
          {
            xx = as.vector(t(ynew[(i-ar):i,(j-ar):(j+ar)]))
            y_new1 = as.matrix(xx[(length(xx)-1):1])
            
            eta[i,j]  = alpha + phi%*%y_new1 
          }
          
        }
        
        mu = as.vector(exp(t(eta[(m+1):n,(m+1):(k-m)])))
        y1 = as.vector(t(y[(m+1):n,(m+1):(k-m)]))
        
        dmu = as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu)))
        
        mT = diag(mu)
        a = as.vector(rep(1,(n-m)*(k-2*m)))
        
        Ualpha =  - t(a) %*% mT %*% dmu
        Uphi =    - t(P) %*% mT %*% dmu
        
        rval = c(Ualpha,Uphi)
      }
      
      if(any(is.na(ar))==F) names_phi = c(paste("phi",1:p1,sep=""))
      
      names_par <- c("alpha",names_phi)
      
      opt = optim(reg, loglik, 
                  escore,
                  method = "BFGS", 
                  hessian = TRUE,
                  control = list(maxit = 400, abstol = 1e-6,
                                 factr = 1e20))
      
      
      if (opt$conv != 0)
        warning("FUNCTION DID NOT CONVERGE!")
      
      z = c()
      z$conv = opt$conv
      coef = (opt$par)
      names(coef) = names_par
      z$coeff = coef
      
      alpha = coef[1]
      phi = coef[2:(p1+1)]
      
      z$alpha = alpha
      z$phi = phi
      
      etahat = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):(k-m))
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):(j+ar)]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          
          etahat[i,j]  = alpha + phi%*%y_new1 
          
        }
        
      }
      
      
      z$fitted = exp(etahat[(m+1):n,(m+1):(k-m)])
      z$etahat = etahat[(m+1):n,(m+1):(k-m)]
      
      y1 = as.vector(t(y[(m+1):n,(m+1):(k-m)]))
      
      muhat2 = as.vector(exp(t(etahat[(m+1):n,(m+1):(k-m)])))
      W = diag(((4)/(muhat2^2))*(muhat2^2))
      
      a = as.vector(rep(1,(n-m)*(k-2*m)))
      
      Kaa = t(a) %*% W %*% a
      Kpa = t(P) %*% W %*% a
      Kap = t(Kpa)
      Kpp = t(P) %*% W %*% P
      
      K = rbind(
        cbind(Kaa,Kap),
        cbind(Kpa,Kpp)
      )
      
      
    }
    
    ######### MA model
    if(any(is.na(ar)==T) && any(is.na(ma)==F))
    { 
      p1 = 0
      q1 = 2*(q+1)^2-2-ma
      reg = c(mqo[1],rep(0,q1)) # initializing the parameter values
      
      loglik = function(z) 
      {
        alpha = z[1]
        theta = z[2:(q1+1)]
        
        eta = error = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):n)
        {
          for (j in (m+1):(k-m))
          {
            yy = as.vector(t(error[(i-ma):i,(j-ma):(j+ma)]))
            error_new = as.matrix(yy[(length(yy)-1):1])
            
            eta[i,j]  = alpha + theta%*%error_new
            error[i,j] = ynew[i,j]-eta[i,j]
          }
          
        }
        
        mu = as.vector(exp((t(eta[(m+1):n,(m+1):(k-m)]))))
        y1 = as.vector(t(y[(m+1):n,(m+1):(k-m)]))
        
        ll = suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2))))
        
        sum(ll)
        
      } 
      
      escore = function(z) 
      {
        alpha = z[1]
        theta = z[2:(q1+1)]
        
        eta = error = matrix(0, ncol=k,nrow=n)
        
        for (i in (m+1):n)
        {
          for (j in (m+1):(k-m))
          {
            
            yy = as.vector(t(error[(i-ma):i,(j-ma):(j+ma)]))
            error_new = as.matrix(yy[(length(yy)-1):1])
            
            eta[i,j]  = alpha + theta%*%error_new
            error[i,j] = ynew[i,j]-eta[i,j]
          }
          
        }
        
        mu = as.vector(exp(t(eta[(m+1):n,(m+1):(k-m)])))
        y1 = as.vector(t(y[(m+1):n,(m+1):(k-m)]))
        
        dmu = as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu)))
        
        mT = diag(mu)
        
        ###################################
        
        # dalpha
        deta.dalpha = matrix(0, nrow = n, ncol=k)
        
        for(i in (m+1):n)
        {
          for(j in (m+1):(k-m))
          {
            xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):(j+ma)]))
            deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
            deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
          }
          
        }
        
        a = as.vector(deta.dalpha[(m+1):n,(m+1):(k-m)])
        
        #dtheta
        XX = c()
        
        for (i in (m+1):n)
        {
          for (j in (m+1):(k-m))
          {
            xx1 = as.vector(t(error[(i-ma):i,(j-ma):(j+ma)]))
            XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
          }
        }
        
        R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
        deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
        dsum.theta = matrix(0, ncol=q1,nrow=q1)
        
        for(i in m:dim(R)[1])
        {
          if( i == m)
          {
            deta.dtheta[i,]<- R[i,] 
          }else{
            for(j in 1:q1)
            {
              dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
            }
            
            deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
          }
        }
        rR <- deta.dtheta[-m,]
        
        ###################################
        
        Ualpha =  - a %*% mT %*% dmu
        Utheta =  - t(rR) %*% mT %*% dmu
        
        rval = c(Ualpha,Utheta)
      }
      
      if(any(is.na(ma))==F) names_theta = (paste("theta",1:q1,sep=""))
      
      names_par <- c("alpha",names_theta)
      
      opt = optim(reg, loglik, 
                  escore,
                  method = "BFGS", 
                  hessian = TRUE,
                  control = list(maxit = 400, abstol = 1e-6,
                                 factr = 1e20))
      
      
      if (opt$conv != 0)
        warning("FUNCTION DID NOT CONVERGE!")
      
      z = c()
      z$conv = opt$conv
      coef = opt$par
      names(coef) = names_par
      z$coeff = coef
      
      alpha = coef[1]
      theta = coef[2:(q1+1)]
      
      z$alpha = alpha
      z$theta = theta
      
      etahat = errorhat = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):(k-m))
        {
          
          yy = as.vector(t(errorhat[(i-ma):i,(j-ma):(j+ma)]))
          errorhat_new = as.matrix(yy[(length(yy)-1):1])
          
          etahat[i,j]  = alpha + theta%*%errorhat_new
          errorhat[i,j] = ynew[i,j]-etahat[i,j]
        }
        
      }
      
      
      z$fitted = exp(etahat[(m+1):n,(m+1):(k-m)])
      z$etahat = etahat[(m+1):n,(m+1):(k-m)]
      z$errorhat = errorhat[(m+1):n,(m+1):(k-m)]
      
      y1 = as.vector(t(y[(m+1):n,(m+1):(k-m)]))
      
      ###################################
      
      # dalpha
      deta.dalpha = matrix(0, nrow = n, ncol=k)
      
      for(i in (m+1):n)
      {
        for(j in (m+1):(k-m))
        {
          xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):(j+ma)]))
          deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
          deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
        }
        
      }
      
      a = as.vector(deta.dalpha[(m+1):n,(m+1):(k-m)])
      
      #dtheta
      XX = c()
      
      for (i in (m+1):n)
      {
        for (j in (m+1):(k-m))
        {
          xx1 = as.vector(t(errorhat[(i-ma):i,(j-ma):(j+ma)]))
          XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        }
      }
      
      R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
      deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
      dsum.theta = matrix(0, ncol=q1,nrow=q1)
      
      for(i in m:dim(R)[1])
      {
        if( i == m)
        {
          deta.dtheta[i,]<- R[i,] 
        }else{
          for(j in 1:q1)
          {
            dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
          }
          deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
        }
      }
      rR <- deta.dtheta[-m,]
      
      
      ###################################
      
      muhat2 = as.vector(exp(t(etahat[(m+1):n,(m+1):(k-m)])))
      W = diag(((4)/(muhat2^2))*(muhat2^2))
      
      Kaa = t(a) %*% W %*% a
      Kta = t(rR) %*% W %*% a
      Kat = t(Kta)
      Ktt = t(rR) %*% W %*% rR
      
      K = rbind(
        cbind(Kaa,Kat),
        cbind(Kta,Ktt)
      )
      
    }
    
    z$image = y
    y1 = y[(m+1):n,(m+1):(k-m)]
    z$fitted = exp(z$etahat)
    
    ##############################################
    
    # residuals
    
    z$resid = qnorm(pr(y1,z$fitted)) # quantile residuals
    # 
    
    ############################################
    
    vcov = chol2inv(chol(K))
    z$vcov = vcov
    
    stderror = sqrt(diag(vcov))
    z$stderror = stderror
    
    z$zstat = abs(z$coef/stderror)
    z$pvalues = 2*(1 - pnorm(z$zstat) )
    
    model_presentation = cbind(round(z$coef,4),round(z$stderror,4),
                               round(z$zstat,4),round(z$pvalues,4))
    colnames(model_presentation)=c("Estimate","Std. Error","z value","Pr(>|z|)")
    
    z$model = model_presentation
    
    # 
    return(z)
    
  }
}

