
elev.com <- function( y, x, groups, data=NULL, method="SMA", alpha=0.05, V=array( 0, c( 2,2,length(unique(groups)) ) ), group.names=sort(unique(groups)) )
{
    if ( is.null(data)==FALSE )
    {
        attach(data)
    }

    x      <- as.matrix( x )
    y      <- as.matrix( y )
    groups <- as.matrix( groups )

    nobs <- length( groups )
    g    <- length( group.names )

    res  <- slope.com( y, x, groups, method=method, V=V, bs=FALSE, ci=FALSE )
    lr   <- res$lr
    p    <- res$p
    b    <- res$b
    varb <- res$varb

    n      <- matrix( 0, g, 1 )
    varres <- matrix( 0, g, 1 )
    means  <- matrix( 0, g, 2 )
    res    <- y - b*x

    for ( i in 1:g )
    {
        iref       <- ( groups==group.names[i] )
        iref       <- iref & ( is.na(x+y) == FALSE )
        n[i]       <- sum( iref )
        means[i,1] <- mean( y[iref] ) 
        means[i,2] <- mean( x[iref] )
        varres[i]  <- ( var( res[iref] ) - V[1,1,i] - b^2*V[2,2,i] )
    }
    varres <- varres *( n - 1 )/( n - 2 )

    as     <- means[,1] - b*means[,2]
    names(as) <- group.names
    varas  <- diag( array( varres/n ) ) + varb*means[,2]%*%t( means[,2] )
    
 #   print(cbind(means[,2], varres/n,as))
 #   print(c(varb, n))
 #   print(means[,2]%*%t( means[,2] ))
 #   print(varas)
    
    varas[n==1,] <- 0 #For singleton groups
    varas[,n==1] <- 0
    df     <- g - 1 - sum(n==1)
    L      <- matrix(0,df,g)
    L[,n>1] <- cbind( matrix( 1, df, 1), diag( array( -1, df), nrow=df ) )
    print(L)    
    stat   <- t(L%*%as)%*%( solve( L%*%varas%*%t(L) ) )%*%(L%*%as)
    
    pvalue <- 1 - pchisq( stat, df ) # remove during sims
    sinv=matrix(0,g,g)
    sinv[n>1,n>1]   <- solve( varas[n>1,n>1] )
    a      <- (matrix(1,1,g)%*%sinv%*%as)/ sum( sum( sinv) )
    vara   <- 1 / sum( sum( sinv ) )
    crit   <- qchisq( 1 - alpha, 1 )
    a.ci   <- c( a - sqrt( crit*vara), a + sqrt( crit*vara ) )

    if ( is.null(data)==FALSE )
    {
        detach(data)
    }

    list( stat=stat, p=pvalue , a=a, ci=a.ci, as=as )
}
