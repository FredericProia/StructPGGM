library(lbfgs)
library(ks)


# vectorialisation symetrique avec elements diagonaux en tete
vechd = function(X, q){
  VX = vech(X)
  if (q > 1){
    L = cumsum(c(1, q:2))
    R = c(VX[L], VX[-L])
  } else{
    R = VX[1]
  }
  return(R)
}


# operation inverse de la vectorialisation precedente
invvechd = function(VX, q){
  if (q > 1){
  L = cumsum(c(1, q:2))
  Lb = setdiff(1:length(VX), L)
  R = VX
  R[L] = VX[1:q]
  R[Lb] = VX[(q+1):length(VX)]
  } else{
    R = VX[1]
  }
  return(invvech(R))
}


# fonction a minimiser pour l'optimisation de Oyy
VraisOyy = function(x, EO, StEO, LtEO, EOStEO, EOLtEO, Syy, eta, L, b, q){
  X = invvechd(x, q)
  X = as.matrix(X)
  iX = solve(X)
  if (sum(is.nan(x)) == 0){
    if (min(eigen(X)$value) <= 0){
      Obj = +Inf
    } else{
      Obj = -log(det(X)) + sum(diag(Syy%*%X)) + sum(diag(StEO%*%iX%*%EO)) + eta*(sum(diag(LtEO%*%iX%*%EO)))^b
    }
  } else{
    Obj = +Inf
  }
  return(Obj)
}


# gradient de la fonction a minimiser pour l'optimisation de Oyy
GradOyy = function(x, EO, StEO, LtEO, EOStEO, EOLtEO, Syy, L, eta, b, q){
  X = invvechd(x, q)
  iX = solve(X)
  Grad = -iX + Syy - iX%*%EOStEO%*%iX - eta*b*(iX%*%EOLtEO%*%iX)*(sum(diag(LtEO%*%iX%*%EO)))^(b-1)
  return(vechd(Grad, q))
}


# fonction a minimiser pour l'optimisation de Oyx
VraisOyx = function(x, iEstOyy, Syx, Sxx, L, eta, b, q, p){
  X = invvec(x, nrow=q, ncol=p)
  #iOyy = solve(EstOyy)
  if (sum(is.nan(x)) == 0){
    Obj = sum(diag(Sxx%*%t(X)%*%iEstOyy%*%X)) + 2*sum(diag(t(Syx)%*%X)) + eta*(sum(diag(L%*%t(X)%*%iEstOyy%*%X)))^b
  } else{
    Obj = +Inf
  }
  return(Obj)
}


# gradient de la fonction a minimiser pour l'optimisation de Oyx
GradOyx = function(x, iEstOyy, Syx, Sxx, eta, L, b, q, p){
  X = invvec(x, nrow=q, ncol=p)
  #iOyy = solve(EstOyy)
  Grad = 2*iEstOyy%*%X%*%Sxx + 2*Syx + 2*eta*b*(iEstOyy%*%X%*%L)*(sum(diag(L%*%t(X)%*%iEstOyy%*%X)))^(b-1)
  return(vec(Grad))
}


# procedure d'estimation : alternance d'optimisations Oyy/Oyx jusqu'a convergence
Estim = function(XSample, YSample, lam, mu, eta, b, L, itmax, VarOracle, cv, affich){
  n = nrow(XSample)
  p = ncol(XSample)
  q = ncol(YSample)
  
  Sxx = t(XSample)%*%XSample/n
  Syy = t(YSample)%*%YSample/n
  Syx = t(YSample)%*%XSample/n
  SSample = matrix(0, nrow=q+p, ncol=q+p)
  SSample[1:q, 1:q] = Syy
  SSample[1:q, (q+1):(q+p)] = Syx
  SSample[(q+1):(q+p), 1:q] = t(Syx)
  SSample[(q+1):(q+p), (q+1):(q+p)] = Sxx
  OSample = solve(SSample+diag(p+q)/10)
  if (n >= 15){
    InitOyy = OSample[1:q, 1:q]
    InitOyx = OSample[1:q, (q+1):(q+p)]
  } else{
    InitOyy = diag(q)
    InitOyx = matrix(0, nrow=q, ncol=p)
  }
  EstOyy = InitOyy
  EstOyx = InitOyx
  
  contgen = TRUE
  cptgen = 0
  while( (contgen) & (cptgen < 3*itmax) ){
    cptgen = cptgen+1
    if (cptgen <= 3*itmax/10){
      eps = 5*10^(-3)
    } else if (cptgen <= 2*3*itmax/10){
      eps = 10^(-2)
    } else if (cptgen <= 5*3*itmax/10){
      eps = 5*10^(-2)
    } else{
      eps = 10^(-1)
    }
    
    if (is.null(VarOracle)){
      EstOyyp = EstOyy
      
      if (q == 1){
        EO = vec(EstOyx)
        tEO = EO
      } else{
        EO = EstOyx
        tEO = t(EstOyx)
      }
      StEO = Sxx%*%tEO
      LtEO = L%*%tEO
      EOStEO = EO%*%Sxx%*%tEO
      EOLtEO = EO%*%L%*%tEO
        
      cont = TRUE
      cptyy = 0
      while( (cont) & (cptyy < itmax) ){
        cptyy = cptyy+1
        precyy = 10^(-5 + ceiling(cptyy/round(itmax/4)))
        #mlsyy = 1000
        mlsyy = 100
        maxityy = 100
        #maxityy = 1000
        if (cv){
          precyy = precyy*5
          mlsyy = mlsyy/10
          maxityy = 500
        }
        #Init = InitOyy + ifelse(q > 1, diag(runif(q, -0.1, 0.1)), runif(1, -0.1, 0.1))
        Init = EstOyyp
        Lbfgsyy = lbfgs(call_eval = VraisOyy, call_grad = GradOyy, vars=vechd(Init, q), EO = EO, StEO = StEO, LtEO = LtEO, EOStEO = EOStEO, EOLtEO = EOLtEO, Syy = Syy, L = L, eta = eta, b = b, q = q, orthantwise_c=lam, orthantwise_start=q, gtol=0.1, epsilon=precyy, max_linesearch=mlsyy, max_iterations=maxityy, invisible=1)
        EstOyy = invvechd(Lbfgsyy$par, q)
        #if ( (Lbfgsyy$convergence != 0) | (Lbfgsyy$value == +Inf) | (is.nan(Lbfgsyy$value)) ){
        if ( (Lbfgsyy$value == +Inf) | (is.nan(Lbfgsyy$value)) ){
          cont = TRUE
        } else{
          cont = FALSE
        }
      }
      varyy = norm(as.matrix(EstOyy - EstOyyp))
      #if (varyy <= eps*max(1, norm(as.matrix(EstOyyp)))){
        #contgen = FALSE
      #}
    } else{
        EstOyyp = VarOracle
        EstOyy = VarOracle
        varyy = 0
        cptyy = 0
        precyy = 0
    }

    EstOyxp = EstOyx
    cont = TRUE
    cptyx = 0
    Init = InitOyx
    iEstOyy = solve(EstOyy)
    while( (cont) & (cptyx < itmax) ){
      cptyx = cptyx+1
      precyx = 10^(-5 + ceiling(cptyx/round(itmax/4)))
      #mlsyx = 1000
      mlsyx = 100
      maxityx = 100
      #maxityx = 1000
      if (cv){
        precyx = precyx*5
        mlsyx = mlsyx/10
        maxityx = 500
      }
      #Init = InitOyx + matrix(runif(p*q, -0.1, 0.1), nrow=q, ncol=p)
      Init = EstOyxp
      Lbfgsyx = lbfgs(call_eval = VraisOyx, call_grad = GradOyx, vars=vec(Init), iEstOyy = iEstOyy, Syx = Syx, Sxx = Sxx, eta = eta, L = L, b = b, q = q, p = p, orthantwise_c=mu, gtol=0.1, epsilon=precyx, max_linesearch=mlsyx, max_iterations=maxityx, invisible=1)
      EstOyx = invvec(Lbfgsyx$par, nrow=q, ncol=p)
      #if ( (Lbfgsyx$convergence != 0) | (Lbfgsyx$value == +Inf) | (is.nan(Lbfgsyx$value)) ){
      if ( (Lbfgsyx$value == +Inf) | (is.nan(Lbfgsyx$value)) ){
        cont = TRUE
      } else{
        cont = FALSE
      }
    }
    varyx = norm(as.matrix(EstOyx - EstOyxp))
    
    #eps = 10^(-4 + ceiling(cptgen/round(5*itmax/3)))
    if (affich){
      cat(paste("Itération ", cptgen, " - (Oyy : it ", cptyy, " [prec ", precyy, "] / Oyx : it ", cptyx, " [prec ", precyx, "]) - err ", round(varyy/max(1, norm(as.matrix(EstOyyp))), 7), "/", round(varyx/max(1, norm(as.matrix(EstOyxp))), 7), " [stop ", eps, "]\n", sep = ""))
    }
    
    if ( (varyy <= eps*max(1, norm(as.matrix(EstOyyp)))) & (varyx <= eps*max(1, norm(as.matrix(EstOyxp)))) ){
      contgen = FALSE
    }
  }
  
  return(list(EstOyy = EstOyy, EstOyx = EstOyx))
}

