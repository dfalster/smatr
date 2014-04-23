
line_df.cis <- function( y, x, alpha=0.05, data=NULL, method="SMA", intercept=TRUE, V=matrix(0,2,2), f.crit=0 )
{

    if ( is.null(data)==FALSE )
    {
        attach(data)
    }
    iref <- ( is.na(x+y) == FALSE ) #to remove NA cases
    n    <- sum(iref)

    # if the line is forced through the origin, df are n-1 not n-2
    if ( intercept == TRUE )
    {
        res.df <- n-2
    }
    else
    { 
        res.df <- n-1 
    }

    if ( f.crit==0 )
    {
        f.crit <- qf( 1 - alpha, 1, res.df )
    }

    dat  <- data.frame( y, x )
    datm <- as.matrix( dat[iref,] )
    #if the line is forced through the origin, SS are estimated without centring the data.
    if ( intercept == TRUE ) 
    {
        vr <- ( var(dat[iref,]) - V )*(n-1) 
    }
    else 
    {
        vr <- t(datm[iref,])%*%datm[iref,] - V*n
    }

    r   <- vr[1,2] / sqrt( vr[1,1]*vr[2,2] )
    cis <- matrix( 0, 2, 2)

    if ( (method==0) | (method=="OLS") )
    {
        lab      <- "coef(reg)"
        b        <- vr[1,2] / vr[2,2]
        var.res  <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] ) / res.df
        var.b    <- var.res / vr[2,2]
        cis[2,1] <- b - sqrt(var.b)*sqrt(f.crit)
        cis[2,2] <- b + sqrt(var.b)*sqrt(f.crit)
    }
    if ( (method==1) | (method=="SMA") )
    {
        lab      <- "coef(SMA)"
        b        <- sign( vr[1,2] ) * sqrt( vr[1,1] / vr[2,2] )
        bigb     <- f.crit * ( 1 - r^2 ) / res.df
        cis[2,1] <- b*( sqrt(bigb+1) - sqrt(bigb) )
        cis[2,2] <- b*( sqrt(bigb+1) + sqrt(bigb) )
        var.res  <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] ) / res.df
        var.b    <- ( vr[1,1] - vr[1,2]^2/vr[2,2] ) / res.df/vr[2,2]
    }
    if ( (method==2) | (method=="MA") )
    {
        lab      <- "coef(MA)"
        fac      <- vr[1,1] - vr[2,2]
        b        <- ( fac + sqrt( fac^2 + 4*vr[1,2]^2 ) ) / 2 / vr[1,2]
        q        <- f.crit*( vr[1,1]*vr[2,2] - vr[1,2]^2 ) / res.df
        cis[2,1] <- (fac + sqrt( fac^2 + 4*vr[1,2]^2 - 4*q ) ) / 2 / ( vr[1,2] + sqrt(q) )
        cis[2,2] <- (fac + sqrt( fac^2 + 4*vr[1,2]^2 - 4*q ) ) / 2 / ( vr[1,2] - sqrt(q) )
        if ( (fac^2 + 4*vr[1,2]^2 - 4*q ) < 0 )
        {
            cis[2,1] <- -Inf
            cis[2,2] <-  Inf
        }
        var.res  <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] ) / res.df
        var.fit  <- ( b^2*vr[1,1] + 2*b*vr[1,2] + vr[2,2] ) / res.df
        var.b    <- 1 / ( var.res/var.fit + var.fit/var.res - 2 )*( 1 + b^2 )^2 / res.df    # Use Fisher info
     
    }

    if (intercept == TRUE)
    {
	  means    <- apply(datm,2,mean)
        a        <- means[1] - b*means[2]
	  var.a    <- var.res/n + var.b*means[2]^2
	  
    
        cis[1,1] <- a - sqrt(var.a)*sqrt(f.crit)
        cis[1,2] <- a + sqrt(var.a)*sqrt(f.crit)
    }
    else
    {
        a        <- 0
        cis[1,]  <- NA
    }
   # print("syy, sxx, Sxy, Var b, a, resid")
   # print(c(vr[1,1], vr[2,2], vr[1,2], var.b, var.a, var.res) )
    
    coeff           <- rbind( a, b )
    coef.names      <- c( "elevation", "slope" )
    coeff           <- data.frame( coeff, cis )
    names(coeff)    <- c( lab, "lower limit", "upper limit" )
    rownames(coeff) <- coef.names

    if ( is.null(data)==FALSE )
    {
        detach(data)
    }
  
  
   print(c(b, var.b, var.res, 0, res.df)) 
    return(coeff)
}
