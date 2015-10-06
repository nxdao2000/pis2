## this file defines methods for the 'pmif' and 'pmifList' classes

## extract the estimated log likelihood
setMethod('logLik','pmif',function(object,...)object@loglik)

## pmifList class
setClass(
         'pmifList',
         contains='list',
         validity=function (object) {
           if (!all(sapply(object,is,'pmif'))) {
             retval <- paste0(
                              "error in ",sQuote("c"),
                              ": dissimilar objects cannot be combined"
                              )
             return(retval)
           }
           d <- sapply(object,function(x)dim(x@conv.rec))
           if (!all(apply(d,1,diff)==0)) {
             retval <- paste0(
                              "error in ",sQuote("c"),
                              ": to be combined, ",sQuote("pmif"),
                              " objects must have chains of equal length"
                              )
             return(retval)
           }
           TRUE
         }
         )

setMethod(
          'c',
          signature=signature(x='pmif'),
          definition=function (x, ...) {
            y <- list(...)
            if (length(y)==0) {
              new("pmifList",list(x))
            } else {
              p <- sapply(y,is,'pmif')
              pl <- sapply(y,is,'pmifList')
              if (any(!(p||pl)))
                stop("cannot mix ",sQuote("pmif"),
                     " and non-",sQuote("pmif")," objects")
              y[p] <- lapply(y[p],list)
              y[pl] <- lapply(y[pl],as,"list")
              new("pmifList",c(list(x),y,recursive=TRUE))
            }
          }
          )

setMethod(
          'c',
          signature=signature(x='pmifList'),
          definition=function (x, ...) {
            y <- list(...)
            if (length(y)==0) {
              x
            } else {
              p <- sapply(y,is,'pmif')
              pl <- sapply(y,is,'pmifList')
              if (any(!(p||pl)))
                stop("cannot mix ",sQuote("pmif"),
                     " and non-",sQuote("pmif")," objects")
              y[p] <- lapply(y[p],list)
              y[pl] <- lapply(y[pl],as,"list")
              new("pmifList",c(as(x,"list"),y,recursive=TRUE))
            }
          }
          )

setMethod(
          "[",
          signature=signature(x="pmifList"),
          definition=function(x, i, ...) {
            new('pmifList',as(x,"list")[i])
          }
          )

## extract the convergence record as a coda::mcmc object
setMethod(
          'conv.rec',
          signature=signature(object='pmif'),
          function (object, pars, ...) {
            if (missing(pars)) pars <- colnames(object@conv.rec)
            coda::mcmc(object@conv.rec[,pars,drop=FALSE])
          }
          )

## extract the convergence records as a coda::mcmc.list object
setMethod(
          'conv.rec',
          signature=signature(object='pmifList'),
          definition=function (object, ...) {
            f <- selectMethod("conv.rec","pmif")
            coda::mcmc.list(lapply(object,f,...))
          }
          )

## plot pmif object
setMethod(
          "plot",
          signature=signature(x='pmif'),
          function (x, y, ...) {
            if (!missing(y)) {
              y <- substitute(y)
              warning(sQuote(y)," is ignored")
            }
            pmif.diagnostics(list(x))
          }
          )


setMethod(
          "plot",
          signature=signature(x='pmifList'),
          definition=function (x, y, ...) {
            if (!missing(y)) {
              y <- substitute(y)
              warning(sQuote(y)," is ignored")
            }
            pmif.diagnostics(x)
          }
          )

pmif.diagnostics <- function (z) {
  ## assumes that x is a list of pmifs with identical structure
  mar.multi <- c(0,5.1,0,2.1)
  oma.multi <- c(6,0,5,0)
  xx <- z[[1]]
  estnames <- xx@pars

  ## plot pmif convergence diagnostics
  other.diagnostics <- c("loglik", "log.prior","nfail")
  plotnames <- c(other.diagnostics,estnames)
  nplots <- length(plotnames)
  n.per.page <- min(nplots,10)
  nc <- if (n.per.page<=4) 1 else 2
  nr <- ceiling(n.per.page/nc)
  oldpar <- par(mar=mar.multi,oma=oma.multi,mfcol=c(nr,nc))
  on.exit(par(oldpar)) 
  low <- 1
  hi <- 0
  iteration <- seq(0,xx@Nmcmc)
  while (hi<nplots) {
    hi <- min(low+n.per.page-1,nplots)
    for (i in seq(from=low,to=hi,by=1)) {
      n <- i-low+1
      dat <- sapply(z,conv.rec,pars=plotnames[i])
      matplot(
              y=dat, 
              x=iteration,
              axes = FALSE,
              xlab = "",
              ylab = "",
              type = "l"
              )
      box()
      y.side <- 2
      axis(y.side,xpd=NA)
      mtext(plotnames[i],y.side,line=3)
      do.xax <- (n%%nr==0||n==n.per.page)
      if (do.xax) axis(1,xpd=NA)
      if (do.xax) mtext("pmif iteration",side=1,line=3)
    }  
    low <- hi+1
    mtext("pmif convergence diagnostics",3,line=2,outer=TRUE)
  }
  invisible(NULL)
}

compare.pmif <- function (z) {
  if (!is.list(z)) z <- list(z)
  if (!all(sapply(z,function(x)is(x,'pmif'))))
    stop("compare.pmif error: ",sQuote("z"),
         " must be a pmif object or a list of pmif objects",call.=FALSE)
  warning(sQuote("compare.pmif")," is deprecated.\n",
          "Use ",sQuote("diagnostics")," instead.",call.=FALSE)
  pmif.diagnostics(z)
}
