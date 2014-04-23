
slope.com <- function( y, x, groups, method="SMA", alpha=0.05, data=NULL, intercept=TRUE, V=array( 0, c( 2,2,length(unique(groups)) ) ), group.names=sort(unique(groups)), ci=TRUE, bs=TRUE )
{
    if ( nargs() < 3 )
    {
        stop('Sorry, no can do without three arguments -- Y, X, GROUPS')
    }

    if ( is.null(data)==FALSE )
    {
        attach(data)
    }

    dat    <- cbind(y, x)
    g      <- length(group.names)

    # Find sample size, variances for each group:
    n      <- matrix( 0, g, 1 )
    res.df <- matrix( 0, g, 1 )
    z      <- matrix( 0, g, 3 )
    do.bs  <- bs
    bs     <- matrix( NA, 3, g, dimnames=list(c("slope","lower.CI.lim","upper.CI.lim"),group.names) )
    for (i in 1:g)
    {
        iref   <- ( groups==group.names[i] )
        iref   <- iref & ( is.na(x+y) == FALSE )
        n[i]   <- sum(iref)
        if ( intercept==FALSE )
        {
           xi <- t(dat[iref, ]) %*% dat[iref, ] / n[i] - V[, , i]
        }
        else 
        {
           if (n[i]>1)
               { xi <- cov(dat[iref, ]) - V[, , i] }
           else if (n[i]==1)
               { xi <- matrix(0,2,2) } #leave as zero for n[i]=1
        }
        z[i,]     <- c( xi[1,1], xi[2,2], xi[1,2] )
        if (do.bs==TRUE & n[i]>1)
            {
            slopei    <- slope.test(y[iref],x[iref],method=method, alpha=alpha, V=V[,,i], intercept=intercept )
            bs[,i]    <- c(slopei$b, slopei$ci)
            }
    }
    if (intercept==FALSE)
        { res.df <- n-1 }
    else
        { res.df <- n-2 }

    if ( is.null(data)==FALSE )
    {
        detach(data)
    }

    # Find common slope:
    lambda <- 1 #only actually used for the major axis.
    res    <- b.com.est( z, n, method, lambda, res.df=res.df )

    # Calculate LR:
    dets <- z[,1]*z[,2] - z[,3]^2 #This is l1*l2 under Halt.
    arguments <- list( l1=dets, l2=1, z=z, n=n, method=method, crit=0, lambda=lambda, res.df=res.df)
    LR     <- lr.b.com(res$b, arguments) 
    # if lambda is being estimated, check endpoint LR values:
    if ( (method==3) | ( method=='lamest' ) )
    {
        res0     <- b.com.est( z, n, 2, lambda=10^-9, res.df ) # to find est when lambda=0
        arguments <- list( l1=dets, l2=1, z=z, n=n, method=method, crit=0, lambda=10^-9, res.df=res.df)
        LR0      <- lr.b.com(res0$b, arguments) 
        resInf   <- b.com.est( z, n, 2, 10^9, res.df ) # to find est when lambda=inf
        arguments <- list( l1=dets, l2=1, z=z, n=n, method=method, crit=0, lambda=10^9, res.df=res.df)
        LRinf    <- lr.b.com(resInf$b, arguments) 
        LR       <- min(LR,LR0,LRinf)
        if ( LR==LR0 )    { res <- res0 }
        if ( LR==LRinf )  { res <- resInf }
    }
    
    # Record values for arguments separately
    b      <- res$b
    bi     <- res$bi
    l1     <- res$l1
    l2     <- res$l2
    lambda <- res$lambda
     
    # Calculate P-value:
    Pvalue <- 1 - pchisq( LR, g - 1 - sum(n==1) ) #don't count any singleton groups in df

    # Calculate a CI for common slope
    if ( (method==1) | (method=='SMA') )
    {
       # Getting variance of common slope
       varBs <- ( z[,1] - (z[,3]^2)/z[,2] ) / z[,2]
       varBs <- varBs / res.df
    } else
    if ( (method==2) | (method=='MA') )
    {
       varBs <- 1 / ( l2/l1 + l1/l2 - 2)*( lambda + b^2)^2
       varBs <- varBs / res.df
    
    }
    if ( (method==3) | (method=='lamest') )
    #Still work to be done to calculate CI for lamest.
    {
       varBs <- NA
    }
    varB <- 1 / sum( 1 / varBs, na.rm=TRUE )
 #   print(c(varB, lambda))
    crit <- qchisq( 1 - alpha, 1 )
    if ( (method==3) | (method=='lamest') )
    {
       ci = FALSE
    }
    bCI=NA
    if ( ci == TRUE )
    {
       bCI  <- com.ci( b, varB, crit, z, n, l1, l2, method, lambda, res.df )
    }
    if (lambda==10^-9) { lambda=0 }
    if (lambda==10^9) { lambda=Inf }
    list( LR=LR, p=Pvalue, b=b, ci=bCI, varb=varB, lambda=lambda, bs=bs )
}

