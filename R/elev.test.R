
elev.test <- function( y, x, test.value=0, data=NULL, alpha=0.05, method="SMA", V=matrix(0,2,2) )
{
    if ( is.null(data)==FALSE )
    {
        attach(data)
    }

    iref <- ( is.na(x+y) == FALSE ) #to remove NA cases
    n    <- sum(iref)
    res.df <- n - 2
    fcrit  <- qf( 1-alpha, 1, res.df )
    dat    <- cbind( y[iref], x[iref] )
    vr     <- ( var( dat ) - V )*( n - 1 )
    r      <- vr[1,2]/( ( vr[1,1]*vr[2,2] )^0.5 )

    if ( (method==0) | (method=="OLS") )
    {
        b       <- vr[1,2] / vr[2,2]
        var.res <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] ) / res.df
        var.b   <- var.res / vr[2,2]
    }
    else if ( (method==1) | (method=="SMA") )
    {
        b       <- sign( vr[1,2] )*sqrt( vr[1,1] / vr[2,2] )
        var.res <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] ) / res.df
        var.b   <- ( vr[1,1] - (vr[1,2]^2)/vr[2,2] ) / res.df / vr[2,2]
    }
    else if ( (method==2) | (method=="MA") )
    {
        fac     <- vr[1,1] - vr[2,2]
        b       <- ( fac + sqrt( fac^2 + 4*vr[1,2]^2) ) / 2 / vr[1,2]
    
    
    #elev versions
        var.res <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] ) / res.df
        var.fit <- ( b^2*vr[1,1] + 2*b*vr[1,2] + vr[2,2] ) / res.df / ( 1 + b^2 )
        var.b   <- 1 / ( var.res/var.fit + var.fit/var.res - 2)*( 1 + b^2 )^2 / res.df    # Use Fisher info      
    }


    means    <- apply(dat,2,mean)
    a        <- means[1] - b*means[2]
    var.a    <- var.res/n + var.b*means[2]^2
    t      <- (a - test.value)/sqrt(var.a)
    pvalue <- 2*pt( -abs(t), res.df )

    if ( is.null(data)==FALSE )
    {
        detach(data)
    }

print(c(b, var.b, var.res, 0, res.df))   
    list( v=var.a, t=t, a=a, p=pvalue, a.ci=c( a-sqrt(var.a*fcrit), a+sqrt(var.a*fcrit) ), test.value=test.value )
}
