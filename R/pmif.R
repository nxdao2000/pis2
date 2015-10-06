## define the pmif class
setClass(
         'pmif',
         contains='pfilterd2.pomp',
         slots=c(
           pars = 'character',
           Nmcmc = 'integer',
           cooling.type='character',
           cooling.fraction='numeric',
           method='character',           
		   random.walk.sd = 'numeric',
           conv.rec = 'matrix',
           log.prior = 'numeric',
           post="array"
           )
         )


pmif.internal <- function (object, Nmcmc,
                            start, pars,
                            rw.sd, Np, Nacceptance=0,
                            lag=0,
                            cooling.type, cooling.fraction, cooling.factor,
                            method,
                            tol, max.fail,
                            verbose, 
                            .ndone = 0L,
							paramMatrix = NULL,
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
  
  if (missing(lag) || lag<0)
	lag=0
  start.names <- names(start)
  if (is.null(start.names))
    stop("pmif error: ",sQuote("start")," must be a named vector",call.=FALSE)

  if (missing(rw.sd))
    stop("pmif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
  rw.names <- names(rw.sd)
  rwsdMat <-rw.sd
  if (is.null(rw.names) || any(rw.sd<0))
    stop("pmif error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
  if (!all(rw.names%in%start.names))
    stop("pmif error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
  rw.names <- names(rw.sd[rw.sd>0])
  if (length(rw.names) == 0)
    stop("pmif error: ",sQuote("rw.sd")," must have one positive entry for each parameter to be estimated",call.=FALSE)

  if (missing(pars))
    stop("pmif error: ",sQuote("pars")," must be specified",call.=FALSE)
  if (length(pars)==0)
    stop("pmif error: at least one parameter must be estimated",call.=FALSE)
  if (
      !is.character(pars) ||
      !all(pars%in%start.names) ||
      !all(pars%in%rw.names)
      )
    stop(
         "pmif error: ",
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
            "pmif warning: the variable(s) ",
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
    stop("pmif error: ",sQuote("Np")," must be specified",call.=FALSE)
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
  if (missing(cooling.fraction))
    cooling.fraction<-0.2
  #  stop("mif error: ",sQuote("cooling.fraction")," must be specified",call.=FALSE)
  cooling.fraction <- as.numeric(cooling.fraction)
  if ((length(cooling.fraction)!=1)||(cooling.fraction<0)||(cooling.fraction>1))
    stop("mif error: ",sQuote("cooling.fraction")," must be a number between 0 and 1",call.=FALSE)
  
  cooling <- cooling.function(
                              type=cooling.type,
                              perobs=(method=="pmif2")||(method=="mif3")||(method=="mif4"),
                              fraction=cooling.fraction,
                              ntimes=ntimes
                              )

  if (missing(Nacceptance))
    Nacceptance=0
  if (missing(Nmcmc))
    stop("pmif error: ",sQuote("Nmcmc")," must be specified",call.=FALSE)
  Nmcmc <- as.integer(Nmcmc)
  if (Nmcmc<0)
    stop("pmif error: ",sQuote("Nmcmc")," must be a positive integer",call.=FALSE)
    if ((method=="mif2")&&(Np[1L]!=Np[ntimes+1]))
    stop("the first and last values of ",sQuote("Np")," must agree when method = ",sQuote("mif2"))
  
  if (verbose) {
    cat("performing",Nmcmc,"pmif iteration(s) to estimate parameter(s)",
        paste(pars,collapse=", "))
    cat(" using random-walk with SD\n")
    print(rw.sd)
    cat("using",Np,"particles\n")
  }
  
  theta <- start

    
  sigma <- rep(0,length(start))
  names(sigma) <- start.names
  
  rw.names <- names(rw.sd)
  
  sigma[rw.names] <- rw.sd

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
         sQuote("pmif"),
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
                                lag=0,
                                verbose=verbose,
                                .getnativesymbolinfo=gnsi
                                ),
               silent=FALSE
               )
    if (inherits(pfp,'try-error'))
      stop("pmif error: error in ",sQuote("pfilter2"),call.=FALSE)
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
    cool.sched <- cooling(nt=1,m=.ndone+n)
    sigma.n <- sigma*cool.sched$alpha
    names(sigma.n)<-names(theta)
    
    theta.prop[pars] <- rnorm(n=length(pars),mean=theta[pars],sd=sigma.n)

    ## run the particle filter on the proposed new parameter values
    pfp.prop <- try(
                    pfilter2.internal(
                      object=object,
                      params=theta.prop, 
                      Np=Np,
                      tol=tol,
                      max.fail=max.fail,
                      pred.mean=(n==Nmcmc),
                      pred.var=TRUE,
                      filter.mean=TRUE,
                      cooling=cooling,
                      cooling.m=.ndone+n,
                      .mif2=(method=="mif2"),
                      .corr=(method=="mif5"),
                      .wn =(method=="mif4")||(method=="mif5"),
                      .rw.sd=sigma.n[pars],
                      .transform=FALSE,
                      save.states=FALSE, 
                      save.params=FALSE,
                      lag=0,
                      verbose=verbose,
                      .getnativesymbolinfo=gnsi),
                    silent=FALSE
                    )
     
    if (inherits(pfp.prop,'try-error'))
      stop("pmif error: error in ",sQuote("pfilter2"),call.=FALSE)
    ntimes <- length(time(object))
    theta.prop[pars] <- pfp.prop@filter.mean[pars,ntimes,drop=FALSE]
      #rowMeans(pfp.prop@filter.mean[pars,,drop=FALSE])
    log.prior.prop <- dprior(object,params=theta.prop,log=TRUE,.getnativesymbolinfo=gnsi)
    gnsi <- FALSE

    ## pmif update rule (OK because proposal is symmetric)
    if (runif(1) < exp(pfp.prop@loglik+log.prior.prop-pfp@loglik-log.prior)) {
      pfp <- pfp.prop
      theta <- theta.prop
      paramMatrix <- pfp.prop@paramMatrix
      log.prior <- log.prior.prop
      Nacceptance <- Nacceptance + 1
    }

    ## store a record of this iteration
    conv.rec[n+1,names(theta)] <- theta
    conv.rec[n+1,c(1,2,3)] <- c(pfp@loglik,log.prior,pfp@nfail)

    if (verbose) cat("pmif iteration ",n," of ",Nmcmc," completed\n")

  }
  cat("Acceptance rate:",Nacceptance/Nmcmc,"\n")
  new(
      "pmif",
      pfp,
      params=theta,
      Nmcmc=Nmcmc,
      pars=pars,
      random.walk.sd=rw.sd,
      Np=Np,
      tol=tol,
      conv.rec=conv.rec,
      method=method,
      cooling.type=cooling.type,
      cooling.fraction=cooling.fraction,
      post = paramMatrix,
      paramMatrix=if (method!="mif") paramMatrix else array(data=numeric(0),dim=c(0,0)),
      lag=0,
      
      log.prior=log.prior
      )
}

setMethod(
          "pmif",
          signature=signature(object="pomp"),
          function (object, Nmcmc = 1,
                    start, pars, rw.sd, Np,
                    cooling.type = c("geometric","hyperbolic"),
                    cooling.fraction, cooling.factor,
                    method = c("pmif","punweighted","pfp","pmif2","pmif3","pmif4"),
                    lag=0, 
                    tol = 1e-17, max.fail = 0,
                    verbose = getOption("verbose"),
                    ...) {
            
            method <- match.arg(method)
            if(missing(lag)) lag<-0            
            if (missing(start)) start <- coef(object)
  	    if (missing(rw.sd))
              stop("pmif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
            if (missing(pars)) pars <- names(rw.sd)[rw.sd>0]
            if (missing(Np))
              stop("pmif error: ",sQuote("Np")," must be specified",call.=FALSE)
            cooling.type <- match.arg(cooling.type)  
            pmif.internal(
                           object=object,
                           Nmcmc=Nmcmc,
                           start=start,
                           pars=pars,
                           rw.sd=rw.sd,
                           Np=Np,
			   Nacceptance= 0,
                           tol=tol,
                           cooling.type=cooling.type,
                           cooling.factor=cooling.factor,
                           cooling.fraction=cooling.fraction,
                           lag=0,
                           method=method,
                           max.fail=max.fail,
                           verbose=verbose,
                           ...
                           )
          }
          )

setMethod(
          "pmif",
          signature=signature(object="pfilterd2.pomp"),
          function (object, Nmcmc = 1, Np, tol, ...) {

            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol
            pmif(
                  object=as(object,"pomp"),
                  Nmcmc=Nmcmc,
                  Np=Np,
		  tol=tol,
			lag=0,
                  ...
                  )
          }
          )

setMethod(
          "pmif",
          signature=signature(object="pmif"),
          function (object, Nmcmc,
                    start, pars, rw.sd,
                    Np, cooling.type, cooling.fraction, lag=0,
                    method,tol, max.fail = 0,
                    verbose = getOption("verbose"),
                    
                    ...) {

            if (missing(Nmcmc)) Nmcmc <- object@Nmcmc
            if (missing(start)) start <- coef(object)
            if (missing(pars)) pars <- object@pars
            if (missing(cooling.type)) cooling.type <- object@cooling.type
            if (missing(cooling.fraction)) cooling.fraction <- object@cooling.fraction
            if (missing(method)) method <- object@method
            if (missing(rw.sd)) rw.sd <- object@random.walk.sd
            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol
            if (missing(lag)) lag <- 0
			
            pmif(
                  object=as(object,"pomp"),
                  Nmcmc=Nmcmc,
                  start=start,
                  pars=pars,
                  rw.sd=rw.sd,
                  Np=Np,
                  cooling.type=cooling.type,
                  cooling.fraction=cooling.fraction,
                  lag=0,
                  method=method,
		          tol=tol,
                  max.fail=max.fail,
                  verbose=verbose,
                  ...
                  )
          }
          )

setMethod(
          'continue',
          signature=signature(object='pmif'),
          function (object, Nmcmc = 1, ...) {

            ndone <- object@Nmcmc

            obj <- pmif(
                         object=object,
                         Nmcmc=Nmcmc,
                         lag =0,
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
