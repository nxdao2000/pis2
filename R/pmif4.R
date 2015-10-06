## define the pmif4 class
setClass(
         'pmif4',
         contains='pfilterd2.pomp',
         slots=c(
           transform = "logical",
           ivps = 'character',
           pars = 'character',
           Nmcmc = 'integer',
           particles = 'function',
           var.factor='numeric',
           ic.lag='integer',
           cooling.type='character',
           cooling.fraction='numeric',
           method='character',
           random.walk.sd = 'numeric',
           conv.rec = 'matrix',
           log.prior = 'numeric',
           post="array"
           )
         )


pmif4.internal <- function (object, Nmcmc,
                            start, pars, ivps,
                            particles,
                            rw.sd, Np, Nacceptance=0, var.factor, ic.lag,
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
  
  
  start.names <- names(start)
  if (is.null(start.names))
    stop("pmif4 error: ",sQuote("start")," must be a named vector",call.=FALSE)

  if (missing(rw.sd))
    stop("pmif4 error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
  rw.names <- names(rw.sd)
  rwsdMat <-rw.sd
  if (is.null(rw.names) || any(rw.sd<0))
    stop("pmif4 error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
  if (!all(rw.names%in%start.names))
    stop("pmif4 error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
  rw.names <- names(rw.sd[rw.sd>0])
  if (length(rw.names) == 0)
    stop("pmif4 error: ",sQuote("rw.sd")," must have one positive entry for each parameter to be estimated",call.=FALSE)

  if (missing(pars))
    stop("pmif4 error: ",sQuote("pars")," must be specified",call.=FALSE)
  if (length(pars)==0)
    stop("pmif4 error: at least one parameter must be estimated",call.=FALSE)
  if (
      !is.character(pars) ||
      !is.character(ivps) ||
      !all(pars%in%start.names) ||
      !all(ivps%in%start.names) ||
      any(pars%in%ivps) ||
      any(ivps%in%pars) ||
      !all(pars%in%rw.names) ||
      !all(ivps%in%rw.names)
      )
    stop(
         "pmif4 error: ",
         sQuote("pars")," and ",sQuote("ivps"),
         " must be mutually disjoint subsets of ",
         sQuote("names(start)"),
         " and must have a positive random-walk SDs specified in ",
         sQuote("rw.sd"),
         call.=FALSE
         )
  
  if (!all(rw.names%in%c(pars,ivps))) {
    extra.rws <- rw.names[!(rw.names%in%c(pars,ivps))]
    warning(
            ngettext(length(extra.rws),"mif warning: the variable ",
                     "mif warning: the variables "),
            paste(sQuote(extra.rws),collapse=", "),
            ngettext(length(extra.rws)," has positive random-walk SD specified, but is included in neither ",
                     " have positive random-walk SDs specified, but are included in neither "),
            sQuote("pars")," nor ",sQuote("ivps"),
            ngettext(length(extra.rws),". This random walk SD will be ignored.",
                     ". These random walk SDs will be ignored."),
            call.=FALSE
            )
  }
  rw.sd <- rw.sd[c(pars,ivps)]
  rw.names <- names(rw.sd)
  
  ntimes <- length(time(object))
  if (is.null(Np)) stop("mif error: ",sQuote("Np")," must be specified",call.=FALSE)
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
  
  ic.lag <- as.integer(ic.lag)
  if ((length(ic.lag)!=1)||(ic.lag<1))
    stop("mif error: ",sQuote("ic.lag")," must be a positive integer",call.=FALSE)
  if (ic.lag>ntimes) {
    warning(
            "mif warning: ",sQuote("ic.lag")," = ",ic.lag," > ",ntimes,
            " = length(time(",sQuote("object"),"))",
            " is nonsensical.  Setting ",sQuote("ic.lag")," = ",ntimes,".",
            call.=FALSE
            )
    ic.lag <- length(time(object))
  }
  if ((length(pars)==0)&&(ic.lag<length(time(object)))) {
    warning(
            "mif warning: only IVPs are to be estimated, yet ",sQuote("ic.lag")," = ",ic.lag,
            " < ",ntimes," = length(time(",sQuote("object"),")),",
            " so unnecessary work is to be done.",
            call.=FALSE
            )
  }
  
  ## the following deals with the deprecated option 'cooling.factor'
  if (!missing(cooling.factor)) {
    warning(sQuote("cooling.factor")," is deprecated.\n",
            "See ",sQuote("?mif")," for instructions on specifying the cooling schedule.",
            call.=FALSE)
    cooling.factor <- as.numeric(cooling.factor)
    if ((length(cooling.factor)!=1)||(cooling.factor<0)||(cooling.factor>1))
      stop("mif error: ",sQuote("cooling.factor")," must be a number between 0 and 1",call.=FALSE)
    if (missing(cooling.fraction)) {
      cooling.fraction <- cooling.factor^50
    } else {
      warning("specification of ",sQuote("cooling.factor"),
              " is overridden by that of ",sQuote("cooling.fraction"),
              call.=FALSE)
    }
  }

  if (missing(cooling.fraction))
    stop("mif error: ",sQuote("cooling.fraction")," must be specified",call.=FALSE)
  cooling.fraction <- as.numeric(cooling.fraction)
  if ((length(cooling.fraction)!=1)||(cooling.fraction<0)||(cooling.fraction>1))
    stop("mif error: ",sQuote("cooling.fraction")," must be a number between 0 and 1",call.=FALSE)
  cooling.fraction <- as.numeric(cooling.fraction)
  if ((length(cooling.fraction)!=1)||(cooling.fraction<0)||(cooling.fraction>1))
    stop("mif error: ",sQuote("cooling.fraction")," must be a number between 0 and 1",call.=FALSE)
  
  cooling <- cooling.function(
                              type=cooling.type,
                              perobs=(method=="pmif4")||(method=="mif3")||(method=="mif4"),
                              fraction=cooling.fraction,
                              ntimes=ntimes
                              )

  if (missing(Nacceptance))
    Nacceptance=0
  if (missing(Nmcmc))
    stop("pmif4 error: ",sQuote("Nmcmc")," must be specified",call.=FALSE)
  Nmcmc <- as.integer(Nmcmc)
  if (Nmcmc<0)
    stop("pmif4 error: ",sQuote("Nmcmc")," must be a positive integer",call.=FALSE)
    if ((method=="mif2")&&(Np[1L]!=Np[ntimes+1]))
    stop("the first and last values of ",sQuote("Np")," must agree when method = ",sQuote("mif2"))
  
  if (verbose) {
    cat("performing",Nmcmc,"pmif4 iteration(s) to estimate parameter(s)",
        paste(pars,collapse=", "))
    cat(" using random-walk with SD\n")
    print(rw.sd)
    cat("using",Np,"particles\n")
  }
  
  theta <- start

  sigma <- rep(0,length(start))
  names(sigma) <- start.names
  
  rw.sd <- rw.sd[c(pars,ivps)]
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
  
  if (!all(is.finite(theta[c(pars,ivps)]))) {
    stop(
      sQuote("pmif2"),
      " error: cannot estimate non-finite parameters: ",
      paste(
        pars[!is.finite(theta[pars])],
        collapse=","
      ),
      call.=FALSE
    )
  }
  
  obj <- as(object,"pomp")
  
 if (Nmcmc>0) {
    tmp.pmif <- new("pmif4",object,particles=particles,Np=Np[1L])
  } else {
    pfp <- obj
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
      stop("pmcmc error: error in ",sQuote("pfilter2"),call.=FALSE)
    log.prior <- dprior(object,params=theta,log=TRUE,.getnativesymbolinfo=gnsi)
    gnsi <- FALSE
  } else { ## has been computed previously
    pfp <- .prev.pfp
    log.prior <- .prev.log.prior
  }
  conv.rec[1,names(theta)] <- theta
  conv.rec[1,c(1,2,3)] <- c(pfp@loglik,log.prior,pfp@nfail)
  
  have.parmat <- !(is.null(paramMatrix) || length(paramMatrix)==0)
  newtheta<-theta
  theta.prop <- theta
  for (n in seq_len(Nmcmc)) { # main loop

    ## get the intensity of artificial noise from the cooling schedule
    cool.sched <- cooling(nt=1,m=.ndone+n)
    sigma.n <- sigma*cool.sched$alpha
    names(sigma.n)<-names(theta)
    
    ## initialize the parameter portions of the particles
    P <- try(
             particles(
                       tmp.mif,
                       Np=Np[1L],
                       center=theta,
                       sd=sigma.n*var.factor
                       ),
             silent = FALSE
             )
    if (inherits(P,"try-error")) 
      stop("mif error: error in ",sQuote("particles"),call.=FALSE)

    
    pfp <- try(
               pfilter2.internal(
                                object=obj,
                                params=P, 
                                Np=Np,
                                tol=tol,
                                max.fail=max.fail,
                                pred.mean=(n==Nmcmc),
                                pred.var=((method=="pmif4")||(TRUE)),
                                filter.mean=TRUE,
                                cooling=cooling,
                                cooling.m=.ndone+n,
                                .corr=(method=="mif5"),
                                .wn =TRUE,
                                .rw.sd=sigma.n[pars],
                                .transform=FALSE,
                                save.states=FALSE, 
                                save.params=FALSE,
                                lag=lag,
                                verbose=verbose,
                                .getnativesymbolinfo=gnsi
                                ),
               silent=FALSE
               )
    if (inherits(pfp,"try-error")) 
      stop("mif error: error in ",sQuote("pfilter2"),call.=FALSE)

    gnsi <- FALSE
    
               oldtheta1<-theta
               oldtheta<-theta
               npars<-length(pars)
               names(oldtheta)<-names(theta)
               newtheta[c(pars,ivps)]<-0
               npars<-length(theta)
               Hessian<-array(0,dim=c(npars,npars))
               colnames(Hessian)<-names(theta)
               rownames(Hessian)<-names(theta)
               if (lag>0){
                 
                 phat<-pfp@phats
                 names(phat)<-names(theta)
                 covhat<-pfp@covhats
                 colnames(covhat)<-names(theta)
                 rownames(covhat)<-names(theta)
                 Hessian[c(pars,ivps),c(pars,ivps)]<-covhat[c(pars,ivps),c(pars,ivps)]-ntimes*diag(sigma.n[c(pars,ivps)]^2)
                 newtheta[c(pars,ivps)] <- phat[c(pars,ivps)]-ntimes*oldtheta[c(pars,ivps)]  
                 
                 Hessian[c(pars,ivps),c(pars,ivps)]<-0.5*(Hessian[c(pars,ivps),c(pars,ivps)]+t(Hessian[c(pars,ivps),c(pars,ivps)]))
                 v1 <- cool.sched$gamma*rwsdMat[c(pars,ivps)]^2  
                 newtheta[c(pars,ivps)]<- solve(Hessian[c(pars,ivps),c(pars,ivps)])%*%newtheta[c(pars,ivps)]*v1  
                 theta.prop[c(pars,ivps)]  <-  oldtheta1[c(pars,ivps)]-newtheta[c(pars,ivps)]
                 
               }

    ## run the particle filter on the proposed new parameter values
    pfp.prop <- try(
      pfilter2.internal(
        object=object,
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
      stop("pmif4 error: error in ",sQuote("pfilter2"),call.=FALSE)
    log.prior.prop <- dprior(object,params=theta.prop,log=TRUE,.getnativesymbolinfo=gnsi)
    gnsi <- FALSE

    ## pmif4 update rule (OK because proposal is symmetric)
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

    if (verbose) cat("pmif4 iteration ",n," of ",Nmcmc," completed\n")

  }
  cat("Acceptance rate:",Nacceptance/Nmcmc,"\n")
  new(
      "pmif4",
      pfp,
      params=theta,
      Nmcmc=Nmcmc,
      ivps=ivps,
      pars=pars,
      particles=particles,
      var.factor=var.factor,
      ic.lag=ic.lag,
      random.walk.sd=sigma[rw.names],
      Np=Np,
      tol=tol,
      conv.rec=conv.rec,
      method=method,
      cooling.type=cooling.type,
      cooling.fraction=cooling.fraction,
      post = paramMatrix,
      paramMatrix= paramMatrix,
      lag=lag,
      
      log.prior=log.prior
      )
}

setMethod(
          "pmif4",
          signature=signature(object="pomp"),
          function (object, Nmcmc = 1,
                    start, pars, ivps = character(0),
                    particles, rw.sd, Np, ic.lag, var.factor,
                    cooling.type = c("geometric","hyperbolic"),
                    cooling.fraction, cooling.factor,
                    method = c("pmif4","punweighted","pfp","pmif4","pmif4","pmif4"),
                    lag=0, 
                    tol = 1e-17, max.fail = 0,
                    verbose = getOption("verbose"),
                    ...) {
            
            method <- match.arg(method)
            if(missing(lag)) lag<-0            
            if (missing(start)) start <- coef(object)
  	    if (missing(rw.sd))
              stop("pmif4 error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
            if (missing(pars)) pars <- names(rw.sd)[rw.sd>0]
            if (missing(Np))
              stop("pmif4 error: ",sQuote("Np")," must be specified",call.=FALSE)
            cooling.type <- match.arg(cooling.type)
            if (missing(particles)) { # use default: normal distribution
              particles <- default.pomp.particles.fun
            } else {
              particles <- match.fun(particles)
              if (!all(c('Np','center','sd','...')%in%names(formals(particles))))
                stop(
                     "mif error: ",
                     sQuote("particles"),
                     " must be a function of prototype ",
                     sQuote("particles(Np,center,sd,...)"),
                     call.=FALSE
                     )
            }
              
            pmif4.internal(
                           object=object,
                           Nmcmc=Nmcmc,
                           start=start,
                           pars=pars,
                           ivps=ivps,
                           particles=particles,
                           rw.sd=rw.sd,
                           Np=Np,
			   Nacceptance= 0,
                           tol=tol,
                           cooling.type=cooling.type,
                           cooling.factor=cooling.factor,
                           cooling.fraction=cooling.fraction,
						   var.factor=var.factor,
                           ic.lag=ic.lag,
                           lag=lag,
                           method=method,
                           max.fail=max.fail,
                           verbose=verbose,
                           ...
                           )
          }
          )

setMethod(
          "pmif4",
          signature=signature(object="pfilterd2.pomp"),
          function (object, Nmcmc = 1, Np, tol, ...) {

            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol
            pmif4(
                  object=as(object,"pomp"),
                  Nmcmc=Nmcmc,
                  Np=Np,
		  tol=tol,
                  ...
                  )
          }
          )

setMethod(
          "pmif4",
          signature=signature(object="pmif4"),
          function (object, Nmcmc,
                    start, pars, ivps,
                    particles, rw.sd,
                    Np, ic.lag, lag, var.factor, cooling.type, cooling.fraction, 
                    method,tol, max.fail = 0,
                    verbose = getOption("verbose"),
                    
                    ...) {

            if (missing(Nmcmc)) Nmcmc <- object@Nmcmc
            if (missing(start)) start <- coef(object)
            if (missing(pars)) pars <- object@pars
            if (missing(ivps)) ivps <- object@ivps
            if (missing(particles)) particles <- object@particles
            if (missing(rw.sd)) rw.sd <- object@random.walk.sd
            if (missing(ic.lag)) ic.lag <- object@ic.lag
            if (missing(var.factor)) var.factor <- object@var.factor
            if (missing(cooling.type)) cooling.type <- object@cooling.type
            if (missing(cooling.fraction)) cooling.fraction <- object@cooling.fraction
            if (missing(method)) method <- object@method
            if (missing(rw.sd)) rw.sd <- object@random.walk.sd
            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol
            if (missing(lag)) lag <- object@lag
            pmif4(
                  object=as(object,"pomp"),
                  Nmcmc=Nmcmc,
                  start=start,
                  pars=pars,
                  ivps=ivps,
                  particles=particles,
                  rw.sd=rw.sd,
                  Np=Np,
                  cooling.type=cooling.type,
                  cooling.fraction=cooling.fraction,
                  var.factor=var.factor,
                  ic.lag=ic.lag,
                  lag=lag,
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
          signature=signature(object='pmif4'),
          function (object, Nmcmc = 1, ...) {

            ndone <- object@Nmcmc

            obj <- pmif4(
                         object=object,
                         Nmcmc=Nmcmc,
                         lag =lag,
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


