## define the pmcmc2 class
setClass(
         'pmcmc2',
         contains='pfilterd2.pomp',
         slots=c(
           pars = 'character',
           Nmcmc = 'integer',
           random.walk.sd = 'numeric',
           conv.rec = 'matrix',
           log.prior = 'numeric'
  
           )
         )

pmcmc2.internal <- function (object, Nmcmc,
                            start, pars,
                            rw.sd, Np, Nacceptance=0,
                            tol, max.fail,
                            verbose,
                            .ndone = 0L,
                            .prev.pfp = NULL, .prev.log.prior = NULL,
                            .getnativesymbolinfo = TRUE) {

  object <- as(object,"pomp")
  gnsi <- as.logical(.getnativesymbolinfo)
  .ndone <- as.integer(.ndone)

  if (missing(start))
    stop(sQuote("start")," must be specified",call.=FALSE)
  if (length(start)==0)
    stop(
         sQuote("start")," must be specified if ",
         sQuote("coef(object)")," is NULL",
         call.=FALSE
         )
  start.names <- names(start)
  if (is.null(start.names))
    stop("pmcmc2 error: ",sQuote("start")," must be a named vector",call.=FALSE)

  if (missing(rw.sd))
    stop("pmcmc2 error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
  rw.names <- names(rw.sd)
  if (is.null(rw.names) || any(rw.sd<0))
    stop("pmcmc2 error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
  if (!all(rw.names%in%start.names))
    stop("pmcmc2 error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
  rw.names <- names(rw.sd[rw.sd>0])
  if (length(rw.names) == 0)
    stop("pmcmc2 error: ",sQuote("rw.sd")," must have one positive entry for each parameter to be estimated",call.=FALSE)

  if (missing(pars))
    stop("pmcmc2 error: ",sQuote("pars")," must be specified",call.=FALSE)
  if (length(pars)==0)
    stop("pmcmc2 error: at least one parameter must be estimated",call.=FALSE)
  if (
      !is.character(pars) ||
      !all(pars%in%start.names) ||
      !all(pars%in%rw.names)
      )
    stop(
         "pmcmc2 error: ",
         sQuote("pars"),
         " must be a mutually disjoint subset of ",
         sQuote("names(start)"),
         " and must correspond to positive random-walk SDs specified in ",
         sQuote("rw.sd"),
         call.=FALSE
         )

  if (!all(rw.names%in%pars)) {
    extra.rws <- rw.names[!(rw.names%in%pars)]
    warning(
            "pmcmc2 warning: the variable(s) ",
            paste(extra.rws,collapse=", "),
            " have positive random-walk SDs specified, but are not included in ",
            sQuote("pars"),
            ". These random walk SDs are ignored.",
            call.=FALSE
            )
  }
  rw.sd <- rw.sd[pars]
  rw.names <- names(rw.sd)

  ntimes <- length(time(object))
  if (missing(Np))
    stop("pmcmc2 error: ",sQuote("Np")," must be specified",call.=FALSE)
  if (is.function(Np)) {
    Np <- try(
              vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
              silent=FALSE
              )
    if (inherits(Np,"try-error"))
      stop("if ",sQuote("Np")," is a function, it must return a single positive integer")
  }
  if (length(Np)==1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np)!=(ntimes+1))
    stop(sQuote("Np")," must have length 1 or length ",ntimes+1)
  if (any(Np<=0))
    stop("number of particles, ",sQuote("Np"),", must always be positive")
  if (!is.numeric(Np))
    stop(sQuote("Np")," must be a number, a vector of numbers, or a function")
  Np <- as.integer(Np)

  if (missing(Nacceptance))
    Nacceptance=0
  if (missing(Nmcmc))
    stop("pmcmc2 error: ",sQuote("Nmcmc")," must be specified",call.=FALSE)
  Nmcmc <- as.integer(Nmcmc)
  if (Nmcmc<0)
    stop("pmcmc2 error: ",sQuote("Nmcmc")," must be a positive integer",call.=FALSE)

  if (verbose) {
    cat("performing",Nmcmc,"pmcmc2 iteration(s) to estimate parameter(s)",
        paste(pars,collapse=", "))
    cat(" using random-walk with SD\n")
    print(rw.sd)
    cat("using",Np,"particles\n")
  }

  theta <- start

  conv.rec <- matrix(
                     data=NA,
                     nrow=Nmcmc+1,
                     ncol=length(theta)+3,
                     dimnames=list(
                       rownames=seq(from=0,to=Nmcmc,by=1),
                       colnames=c('loglik','log.prior','nfail',names(theta))
                       )
                     )

  if (!all(is.finite(theta[pars]))) {
    stop(
         sQuote("pmcmc2"),
         " error: cannot estimate non-finite parameters: ",
         paste(
               pars[!is.finite(theta[pars])],
               collapse=","
               ),
         call.=FALSE
         )
  }

  if (.ndone==0L) { ## compute prior and likelihood on initial parameter vector
    pfp <- try(
               pfilter2.internal(
                                object=object,
                                params=theta,
                                Np=Np,
                                tol=tol,
                                max.fail=max.fail,
                                pred.mean=FALSE,
                                pred.var=FALSE,
                                filter.mean=TRUE,
                                save.states=FALSE,
                                save.params=FALSE,
                                .transform=FALSE,
                                verbose=verbose,
                                .getnativesymbolinfo=gnsi
                                ),
               silent=FALSE
               )
    if (inherits(pfp,'try-error'))
      stop("pmcmc2 error: error in ",sQuote("pfilter2"),call.=FALSE)
    log.prior <- dprior(object,params=theta,log=TRUE,.getnativesymbolinfo=gnsi)
    gnsi <- FALSE
  } else { ## has been computed previously
    pfp <- .prev.pfp
    log.prior <- .prev.log.prior
  }
  conv.rec[1,names(theta)] <- theta
  conv.rec[1,c(1,2,3)] <- c(pfp@loglik,log.prior,pfp@nfail)

  for (n in seq_len(Nmcmc)) { # main loop

    theta.prop <- theta
    theta.prop[pars] <- rnorm(n=length(pars),mean=theta[pars],sd=rw.sd)

    ## run the particle filter on the proposed new parameter values
    pfp.prop <- try(
                    pfilter2.internal(
                                     object=pfp,
                                     params=theta.prop,
                                     Np=Np,
                                     tol=tol,
                                     max.fail=max.fail,
                                     pred.mean=FALSE,
                                     pred.var=FALSE,
                                     filter.mean=TRUE,
                                     save.states=FALSE,
                                     save.params=FALSE,
                                     .transform=FALSE,
                                     verbose=verbose,
                                     .getnativesymbolinfo=gnsi
                                     ),
                    silent=FALSE
                    )
    if (inherits(pfp.prop,'try-error'))
      stop("pmcmc2 error: error in ",sQuote("pfilter2"),call.=FALSE)
    log.prior.prop <- dprior(object,params=theta.prop,log=TRUE,.getnativesymbolinfo=gnsi)
    gnsi <- FALSE

    ## pmcmc2 update rule (OK because proposal is symmetric)
    if (runif(1) < exp(pfp.prop@loglik+log.prior.prop-pfp@loglik-log.prior)) {
      pfp <- pfp.prop
      theta <- theta.prop
      log.prior <- log.prior.prop
      Nacceptance <- Nacceptance + 1
    }

    ## store a record of this iteration
    conv.rec[n+1,names(theta)] <- theta
    conv.rec[n+1,c(1,2,3)] <- c(pfp@loglik,log.prior,pfp@nfail)

    if (verbose) cat("pmcmc2 iteration ",n," of ",Nmcmc," completed\n")

  }
  cat("Acceptance rate:",Nacceptance/Nmcmc,"\n")
  new(
      "pmcmc2",
      pfp,
      params=theta,
      Nmcmc=Nmcmc,
      pars=pars,
      random.walk.sd=rw.sd,
      Np=Np,
      tol=tol,
      conv.rec=conv.rec,
      log.prior=log.prior
      )
}

setMethod(
          "pmcmc2",
          signature=signature(object="pomp"),
          function (object, Nmcmc = 1,
                    start, pars, rw.sd, Np, 
                    tol = 1e-17, max.fail = 0,
                    verbose = getOption("verbose"),
                    ...) {
            
            if (missing(start)) start <- coef(object)
  	    if (missing(rw.sd))
              stop("pmcmc2 error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
            if (missing(pars)) pars <- names(rw.sd)[rw.sd>0]
            if (missing(Np))
              stop("pmcmc2 error: ",sQuote("Np")," must be specified",call.=FALSE)
              
            pmcmc2.internal(
                           object=object,
                           Nmcmc=Nmcmc,
                           start=start,
                           pars=pars,
                           rw.sd=rw.sd,
                           Np=Np,
			   Nacceptance= 0,
                           tol=tol,
                           max.fail=max.fail,
                           verbose=verbose,
                           ...
                           )
          }
          )

setMethod(
          "pmcmc2",
          signature=signature(object="pfilterd2.pomp"),
          function (object, Nmcmc = 1, Np, tol, ...) {

            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol
            pmcmc2(
                  object=as(object,"pomp"),
                  Nmcmc=Nmcmc,
                  Np=Np,
		  tol=tol,
                  ...
                  )
          }
          )

setMethod(
          "pmcmc2",
          signature=signature(object="pmcmc2"),
          function (object, Nmcmc,
                    start, pars, rw.sd,
                    Np, tol, max.fail = 0,
                    verbose = getOption("verbose"),
                    ...) {

            if (missing(Nmcmc)) Nmcmc <- object@Nmcmc
            if (missing(start)) start <- coef(object)
            if (missing(pars)) pars <- object@pars
            if (missing(rw.sd)) rw.sd <- object@random.walk.sd
            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol
	    
            pmcmc2(
                  object=as(object,"pomp"),
                  Nmcmc=Nmcmc,
                  start=start,
                  pars=pars,
                  rw.sd=rw.sd,
                  Np=Np,
		  tol=tol,
                  max.fail=max.fail,
                  verbose=verbose,
                  ...
                  )
          }
          )

setMethod(
          'continue',
          signature=signature(object='pmcmc2'),
          function (object, Nmcmc = 1, ...) {

            ndone <- object@Nmcmc

            obj <- pmcmc2(
                         object=object,
                         Nmcmc=Nmcmc,
			 ...,
                         .ndone=ndone,
                         .prev.pfp=as(object,"pfilterd2.pomp"),
                         .prev.log.prior=object@log.prior
                         )
            
            obj@conv.rec <- rbind(
                                  object@conv.rec[,colnames(obj@conv.rec)],
                                  obj@conv.rec[-1,]
                                  )
            obj@Nmcmc <- as.integer(ndone+Nmcmc)
            obj
          }
          )
