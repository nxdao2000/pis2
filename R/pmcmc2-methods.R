## this file defines methods for the 'pmcmc2' and 'pmcmc2List' classes

## extract the estimated log likelihood
setMethod('logLik','pmcmc2',function(object,...)object@loglik)

## pmcmc2List class
setClass(
         'pmcmc2List',
         contains='list',
         validity=function (object) {
           if (!all(sapply(object,is,'pmcmc2'))) {
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
                              ": to be combined, ",sQuote("pmcmc2"),
                              " objects must have chains of equal length"
                              )
             return(retval)
           }
           TRUE
         }
         )

setMethod(
          'c',
          signature=signature(x='pmcmc2'),
          definition=function (x, ...) {
            y <- list(...)
            if (length(y)==0) {
              new("pmcmc2List",list(x))
            } else {
              p <- sapply(y,is,'pmcmc2')
              pl <- sapply(y,is,'pmcmc2List')
              if (any(!(p||pl)))
                stop("cannot mix ",sQuote("pmcmc2"),
                     " and non-",sQuote("pmcmc2")," objects")
              y[p] <- lapply(y[p],list)
              y[pl] <- lapply(y[pl],as,"list")
              new("pmcmc2List",c(list(x),y,recursive=TRUE))
            }
          }
          )

setMethod(
          'c',
          signature=signature(x='pmcmc2List'),
          definition=function (x, ...) {
            y <- list(...)
            if (length(y)==0) {
              x
            } else {
              p <- sapply(y,is,'pmcmc2')
              pl <- sapply(y,is,'pmcmc2List')
              if (any(!(p||pl)))
                stop("cannot mix ",sQuote("pmcmc2"),
                     " and non-",sQuote("pmcmc2")," objects")
              y[p] <- lapply(y[p],list)
              y[pl] <- lapply(y[pl],as,"list")
              new("pmcmc2List",c(as(x,"list"),y,recursive=TRUE))
            }
          }
          )

setMethod(
          "[",
          signature=signature(x="pmcmc2List"),
          definition=function(x, i, ...) {
            new('pmcmc2List',as(x,"list")[i])
          }
          )

## extract the convergence record as a coda::mcmc object
setMethod(
          'conv.rec',
          signature=signature(object='pmcmc2'),
          function (object, pars, ...) {
            if (missing(pars)) pars <- colnames(object@conv.rec)
            coda::mcmc(object@conv.rec[,pars,drop=FALSE])
          }
          )

## extract the convergence records as a coda::mcmc.list object
setMethod(
          'conv.rec',
          signature=signature(object='pmcmc2List'),
          definition=function (object, ...) {
            f <- selectMethod("conv.rec","pmcmc2")
            coda::mcmc.list(lapply(object,f,...))
          }
          )

## plot pmcmc2 object
setMethod(
          "plot",
          signature=signature(x='pmcmc2'),
          function (x, y, ...) {
            if (!missing(y)) {
              y <- substitute(y)
              warning(sQuote(y)," is ignored")
            }
            pmcmc2.diagnostics(list(x))
          }
          )


setMethod(
          "plot",
          signature=signature(x='pmcmc2List'),
          definition=function (x, y, ...) {
            if (!missing(y)) {
              y <- substitute(y)
              warning(sQuote(y)," is ignored")
            }
            pmcmc2.diagnostics(x)
          }
          )

pmcmc2.diagnostics <- function (z) {
  ## assumes that x is a list of pmcmc2s with identical structure
  mar.multi <- c(0,5.1,0,2.1)
  oma.multi <- c(6,0,5,0)
  xx <- z[[1]]
  estnames <- xx@pars

  ## plot pmcmc2 convergence diagnostics
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
      if (do.xax) mtext("pmcmc2 iteration",side=1,line=3)
    }  
    low <- hi+1
    mtext("pmcmc2 convergence diagnostics",3,line=2,outer=TRUE)
  }
  invisible(NULL)
}

compare.pmcmc2 <- function (z) {
  if (!is.list(z)) z <- list(z)
  if (!all(sapply(z,function(x)is(x,'pmcmc2'))))
    stop("compare.pmcmc2 error: ",sQuote("z"),
         " must be a pmcmc2 object or a list of pmcmc2 objects",call.=FALSE)
  warning(sQuote("compare.pmcmc2")," is deprecated.\n",
          "Use ",sQuote("diagnostics")," instead.",call.=FALSE)
  pmcmc2.diagnostics(z)
}
