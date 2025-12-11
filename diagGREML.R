diagGREML = function(y, A, X, varcmp=NULL, is_EM=F, itermax=100, constraint=F, magic0316=F){
    # estimate g(s) and e variance compoents using only the diagonal elements of GRM(s) (DiagGREML)
    # y: phenotypes (n*1 vector), where n is the sample size
    # A: a list of diagonal element vectors (n*1) of r(>=1) GRMs, where r is the number of GRM components
    # X: covariates (n*k matrix), where k is the number of covariates. It should include a column of ones, espeicially when y is not centered.
    # varcmp: the r+1 initial variance values (sum of this vector should be equal to Vp);
    # is_EM: whether the initial variance values will be updated by EM algoritm
    # itermax: the maximum iteration for diagGREML
    # constraint: constraint the variance components in the non-negative value space
    # magic0316: the magic value used in GCTA-GREML to make the change of varaince components in each step smaller. It can bypass some covergence issues.
    # assume: no missing values in GRMs; the individual in y, A, X should be in the same order
    # missing data: remove the individuals with missing data for y or X
    Vp = var(y); prev_logL = -1e20; prevarcmp=varcmp
    is_e_only = sum(sapply(A, function(x) length(which(x!=0))))==0
    Vlist = c(A, list(error=rep(1,length(y))))
    for(iter in 0:itermax){
        ## coverge or next iteration
        if(iter==1 & is_EM){
            prevarcmp = varcmp
            prev_logL = logL
            varcmp = as.numeric((varcmp^2 * yPVPy + varcmp*length(y) - varcmp^2*trPV)/length(y))
        }else if(iter != 0){
            ### decide whehter it's converged
            if(sum((varcmp-prevarcmp)^2)/sum(varcmp^2)<1e-8 & (abs(dlogL)<1e-4 | (abs(dlogL)<1e-2 & dlogL<0))){
              message("converged")
              break
            }else if(is.na(dlogL)){
              stop("dlogL NA")
            }

            prevarcmp = varcmp
            prev_logL = logL
            if(magic0316 & dlogL>1){varcmp = varcmp+delta*0.316
            }else{varcmp = varcmp+delta}
            if(constraint){varcmp = pmax(var(y)*10e-6, varcmp)}
        }

        ## calcu_V (n*n matrix)
        V = Reduce("+", lapply(1:length(Vlist), function(x) Vlist[[x]]*varcmp[x]))
        logdet_V = sum(log(V))
        
        ## logL
        Vi = 1/V
        Vi_X = do.call("cbind", lapply(1:ncol(X),function(j) X[,j]*Vi))   ### n*k matrix
        Xt_Vi_X = crossprod(X, Vi_X)   ### k*k matrix
        Xt_Vi_X_i = solve(Xt_Vi_X, tol=1e-50)
        Vi_X_Xt_Vi_X_i = Vi_X %*% Xt_Vi_X_i   ### n*k matrix
        logdet_Xt_Vi_X = as.numeric(determinant(Xt_Vi_X)$modulus)
        ### P is a dense n*n matrix, so try to avoid calculate P
        Py = (Vi*y) - Vi_X_Xt_Vi_X_i %*% crossprod(Vi_X, y)
        logL = -0.5 * as.numeric(logdet_V + logdet_Xt_Vi_X + crossprod(y, Py))
        dlogL = logL - prev_logL

        ## AI matrix: (r+1)*(r+1) symmetry matrix
        VPy = lapply(Vlist, function(x) x*Py)    ## a list with length r+1 of n*1 vector
        PVPy = lapply(VPy, function(x) Vi*x - Vi_X_Xt_Vi_X_i %*% crossprod(Vi_X, x))
        H = diag(length(varcmp))
        for(i in 1:length(varcmp)){
            for(j in 1:i){
                H[i,j]=H[j,i] = as.numeric(crossprod(VPy[[i]], PVPy[[j]]))/2
        }}
        
        ### DL matrix (r+1)*1 vector
        ### tv(PA)=tr(V^-1A)-tr(V^-1X(X'V^-1X)^-1X'V^-1A)
        ### tr(V^-1X(X'V^-1X)^-1X'V^-1A) = tr(X'V^-1AV^-1X(X'V^-1X)^-1), the former is a n*n matrix, the latter is a k*k matrix
        trPV = sapply(1:length(Vlist), function(i)
            sum(Vi*Vlist[[i]])-sum(diag(crossprod(do.call("cbind", lapply(1:ncol(X), function(j) Vi_X[,j]*Vlist[[i]])), Vi_X_Xt_Vi_X_i))))
        yPVPy = sapply(VPy, function(x) as.numeric(crossprod(Py,x)))
        R = -1/2*(trPV - yPVPy)

        if(is_e_only){Hi=ginv(H)
        }else{Hi=solve(H)}
        delta = Hi %*% R         ### (r+1)*1 vector

        message(paste(date(), "; Iteration: ", iter, "; varcmp: ", paste(round(varcmp,3), collapse=", "),
               "; logL:", round(logL,4), "; dlogL:", round(dlogL, 4), sep=""))
    }
    ### Hi (the inverse of AI matrix is the sampling variance/covariance matrix of the estimates of variance components)
    beta_hat = Xt_Vi_X_i %*% crossprod(Vi_X, y)
    beta_var = Xt_Vi_X_i
    return(list(varcmp=varcmp, logL=logL, Vlistnames=names(Vlist), Xcolnames=colnames(X), Hi=Hi, beta_hat=beta_hat, beta_var=beta_var))
}

fit_diagGREML = function(y, A, X, varcmp=NULL, is_EM=F, itermax=100, constraint=F, magic0316=F, fullonly=F){
    ## y: vector, A: a list; X: a matrix
    library(Matrix);library(MASS)
    message(paste("Start the DiagGREML analysis at ", date(), sep="")); t1 <- Sys.time()
    message(paste(length(y), " individuals in the phenotype vector y", sep=""))
    message(paste(length(A), " GRM components in the GRM list A", sep=""))
    message(paste(dim(X)[1], " individuals for ", dim(X)[2], " covariates in the covariate matrix X", sep=""))
    message("")
    
    ## missing data for phenotypes, covariates, and GRMs
    missing_y = is.na(y); missing_X = apply(X,1,function(x) sum(is.na(x))) != 0
    missing_A = Reduce("+", lapply(A, function(x) is.na(x))) != 0
    message(paste("missing data: ", length(which(missing_y)), " for phenotypes; ",
        length(which(missing_X)), " for covariate matrix; ", 
        length(which(missing_A)), " for GRM components",  sep=""))
    non_missing_index = which(!missing_y & !missing_X & !missing_A)
    y = y[non_missing_index]
    A = lapply(A, function(x) x[non_missing_index])
    X = X[non_missing_index,,drop=F]
    message(paste("after removing the missing data in phenotypes/covariates/GRMs, ", length(non_missing_index), " individuals remained", sep=""))
    message("")

    ### varcmp initation
    if(is.null(varcmp)){varcmp = rep(var(y)/(length(A)+1), length(A)+1)
    }else if(length(varcmp)!=(length(A)+1)){stop("ERROR: the length of varcmp should be the number of GRMs + 1!")}

    ## drop zero-variance covariate(s) and digonal GRM(s)
    if(ncol(X) == 1 & length(which(X!=1))==0){
        message("no covariates are included and the input covariate matrix is a design matrix of 1")
    }else{
        covar_var0 = which(apply(X, 2, function(x) var(x))==0 & !(apply(X, 2, function(x) sum(x!=1))==0))
        if(length(covar_var0)!=0){X = X[,-covar_var0,drop=F]}
        message(paste(ncol(X), " covariates remained, dropping ", ifelse(length(covar_var0)==0, "zero",
            paste(paste(covar_var0, collapse=", ", sep=""), "th", sep="")), " covariates", sep=""))
    }
    A_var0 = which(sapply(A, var)==0)
    if(length(A_var0)!=0){A = A[-A_var0]; varcmp = varcmp[-A_var0]}
    message(paste(length(A), " GRMs remained, dropping ", ifelse(length(A_var0)==0, "zero",
        paste(paste(A_var0, collapse=", ", sep=""), "th", sep="")), " GRMs", sep=""))
    message("")

    full = diagGREML(y=y, A=A, X=X, varcmp=varcmp, is_EM=is_EM, itermax=itermax, constraint=constraint, magic0316=magic0316)
    if(!fullonly){
    reduce = diagGREML(y=y, A=list(A1=rep(0, length(y))), X=X, varcmp=c(0, sum(varcmp)), is_EM=F, itermax=itermax, constraint=constraint, magic0316=F)
    LRT = (full$logL-reduce$logL)*2
    p_LRT = pchisq(LRT, df=length(A), lower.tail=F)
    logL0=reduce$logL
    }else{
    LRT=NA; p_LRT=NA;logL0=NA
    }
    Vp = sum(full$varcmp)
    VarVp = sum(full$Hi)
    h2 = full$varcmp[1:length(A)]/Vp
    Varh2 = sapply(1:length(A), function(x) h2[x]^2*(diag(full$Hi)[x]/full$varcmp[x]^2+VarVp/Vp^2-2*sum(full$Hi[,x])/(full$varcmp[x]*Vp)))
    p_chisq = pchisq(full$varcmp^2/diag(full$Hi), df=1, lower.tail=F)
    p_chisq_beta = pchisq(full$beta_hat^2/diag(full$beta_var), df=1, lower.tail=F)
    npar = ncol(X)+length(A)
    n = nrow(X)
    AIC = 2*npar - 2*full$logL
    BIC = log(n)*npar - 2*full$logL
    message("")
    message(paste("End the DiagGREML analysis at ", date(), sep="")); t2 <- Sys.time()
    library(hms)
    message(paste("overall time: ", hms(hours=as.numeric(difftime(t2,t1, units="hour"))), sep=""))
    return(list(varcmp=full$varcmp, Vlistnames=full$Vlistnames, p_chisq=p_chisq,
                logL = full$logL, logL0=logL0, LRT = LRT, df=length(A), p_LRT=p_LRT, n=n, 
                h2=h2, se_h2 = sqrt(Varh2), Hi=full$Hi,
                Xcolnames=full$Xcolnames, beta_hat=full$beta_hat, beta_var=full$beta_var, p_chisq_beta=p_chisq_beta,
                npar=npar, AIC=AIC, BIC=BIC))
}
