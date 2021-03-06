library(pis2)
library(coda)
library(lattice)
pompExample(ou2)
Nmcmc<-100
NP<-1000


dprior.ou2 <- function (params, log, ...) {
  f <- sum(dunif(params,min=coef(ou2)-1,max=coef(ou2)+1,log=TRUE))
  if (log) f else exp(f)
}
p.truth <- coef(ou2)
guess2 <- guess1 <- p.truth
guess1[c('x1.0','x2.0','alpha.2','alpha.3')] <- 0.8*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]

set.seed(123)

f1 <- pmcmc(
            pomp(ou2,dprior=dprior.ou2),
            start=guess1,
            Nmcmc=Nmcmc,
            rw.sd=c(alpha.2=0.02,alpha.3=0.02),
            Np=NP,
            max.fail=100, 
            verbose=FALSE
            )
set.seed(123)

f2 <- pmif(
  pomp(ou2,dprior=dprior.ou2),
  start=guess1,
  Nmcmc=Nmcmc,
  rw.sd=c(alpha.2=0.02,alpha.3=0.02),
  Np=NP,
  cooling.type="geometric",
  cooling.factor=1,
  method="pmif",
  lag=0,
  max.fail=100, 
  verbose=FALSE
)

set.seed(123)
f3 <- pmif2(
  pomp(ou2,dprior=dprior.ou2),
  start=guess1,
  Nmcmc=Nmcmc,
  rw.sd=c(alpha.2=0.02,alpha.3=0.02),
  Np=NP,
  cooling.type="geometric",
  cooling.factor=1,
  method="pmif2",
  ic.lag=100,
  var.factor=1,
  lag=0,
  max.fail=100, 
  verbose=FALSE
)

set.seed(123)
f4 <- pmif3(
  pomp(ou2,dprior=dprior.ou2),
  start=guess1,
  Nmcmc=Nmcmc,
  rw.sd=c(alpha.2=0.02,alpha.3=0.02),
  Np=NP,
  cooling.type="geometric",
  cooling.factor=1,
  method="pmif3",
  ic.lag=100,
  var.factor=1,
  lag=1,
  max.fail=100, 
  verbose=FALSE
)
set.seed(123)
f5 <- pmif4(
  pomp(ou2,dprior=dprior.ou2),
  start=guess1,
  Nmcmc=Nmcmc,
  rw.sd=c(alpha.2=0.02,alpha.3=0.02),
  Np=NP,
  cooling.type="geometric",
  cooling.factor=1,
  method="pmif4",
  ic.lag=100,
  var.factor=1,
  lag=1,
  max.fail=100, 
  verbose=FALSE
)



pdf(file="ou2-pmcmc.pdf")
result<-mcmc.list(coda::mcmc(conv.rec(f1))[,c("alpha.3","alpha.2","loglik")],coda::mcmc(conv.rec(f2))[,c("alpha.3","alpha.2","loglik")],coda::mcmc(conv.rec(f3))[,c("alpha.3","alpha.2","loglik")],coda::mcmc(conv.rec(f4))[,c("alpha.3","alpha.2","loglik")],coda::mcmc(conv.rec(f5))[,c("alpha.3","alpha.2","loglik")])
save(
     result,
     file="ou2-10000.rda",compress='xz')



plot(conv.rec(f1,c("alpha.2","alpha.3","loglik")))
plot(conv.rec(f2,c("alpha.2","alpha.3","loglik")))
plot(conv.rec(f3,c("alpha.2","alpha.3","loglik")))
plot(conv.rec(f4,c("alpha.2","alpha.3","loglik")))
plot(conv.rec(f5,c("alpha.2","alpha.3","loglik")))

dev.off()

pdf(file="ou2-pmcmc1.pdf")
result<-mcmc.list(coda::mcmc(conv.rec(f1))[,c("alpha.3","alpha.2","loglik")],coda::mcmc(conv.rec(f2))[,c("alpha.3","alpha.2","loglik")],coda::mcmc(conv.rec(f3))[,c("alpha.3","alpha.2","loglik")],coda::mcmc(conv.rec(f4))[,c("alpha.3","alpha.2","loglik")],coda::mcmc(conv.rec(f5))[,c("alpha.3","alpha.2","loglik")])

xyplot(result)
densityplot(result)
qqmath(result, start = (Nmcmc/2))
acfplot(result)
plot(result, trace = TRUE, density = TRUE, smooth = FALSE, auto.layout = TRUE)
dev.off()



