
shift.com <- function( y, x, groups, data=NULL, method="SMA", intercept=TRUE,  V=array( 0, c( 2,2,length(unique(groups)) ) ), group.names=sort(unique(groups)) )
{
    if ( is.null(data)==FALSE )
    {
        attach(data)
    }

    y <- as.matrix(y)
    x <- as.matrix(x)
    groups <- as.matrix(groups)

    nobs <- length(groups)
    g    <- length(group.names)

    res  <- slope.com( y, x, groups, method, data=data, intercept=intercept, V=V, ci=FALSE, bs=FALSE )
   #  res  <- slope.com( y, x, groups, method, V=V, ci=FALSE, bs=FALSE )
    
    #   print(res$lambda)
    lr   <- res$lr
    p    <- res$p
    b    <- res$b
    varb <- res$varb
    
  #  print(res)

    n        <- matrix( 0, g, 1 )
    varAxis  <- n
    as       <- n
    means    <- matrix( 0, g, 2 )

    if ( (method=="SMA") | method==1 )
    {
        axis       <- y + b*x
        coefV1     <- 1 #The coef of V[1,1,:] in var(axis).
        coefV2     <- b^2 #The coef of V[2,2,:] in var(axis).
        mean.ref   <- 2 #Ref for the column of means to use as coef of var(b)
    }
    if ( (method=="MA") | method==2 )
    {
        axis       <- b*y + x
        coefV1     <- b^2 #The coef of V[1,1,:] in var(axis).
        coefV2     <- 1
        mean.ref   <- 1 #Ref for the column of means to use as coef of var(b)
    }
 
    for ( i in 1:g )
    {
        iref       <- ( groups==group.names[i] )
        iref       <- iref & ( is.na(x+y) == FALSE )
        n[i]       <- sum( iref )
        means[i,1] <- mean( y[iref] ) 
        means[i,2] <- mean( x[iref] )
        as[i]      <- mean( axis[iref] )
        varAxis[i] <- var( axis[iref] )
    }
    varAxis    <- varAxis - coefV1*V[1,1,] - coefV2*V[2,2,]
    varAxis    <- varAxis * (n-1) / (n-2)
    mean.for.b <- means[,mean.ref]

    varAs <- diag( array(varAxis/n) ) + varb*mean.for.b%*%t(mean.for.b)
    
    
    print(cbind(mean.for.b, c(varAxis/n), as)) 
    print(varb)
    print(mean.for.b%*%t(mean.for.b))
    print(varAs)
    
    varAs[n==1,] <- 0 #For singleton groups
    varAs[,n==1] <- 0
    df     <- g - 1 - sum(n==1)
    L      <- matrix(0,df,g)
    L[,n>1] <- cbind( matrix( 1, df, 1), diag( array( -1, df), nrow=df ) )
    stat  <- t(L%*%as)%*%solve(L%*%varAs%*%t(L), tol=1.0e-050 )%*%(L%*%as)

    # list( varAxis=varAxis, n=n, varb=varb, means=means )

    pvalue <- 1 - pchisq( stat, df )    # remove during sims

    if ( is.null(data)==FALSE )
    {
        detach(data)
    }

    list( stat=stat, p=pvalue, f.mean=as.vector(as) )

}
