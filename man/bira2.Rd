\name{bira2}
\alias{bira2}
\alias{bira2r1}
\alias{mdes.bira2}
\alias{power.bira2}
\alias{mrss.bira2}
\alias{mdes.bira2r1}
\alias{power.bira2r1}
\alias{mrss.bira2r1}
\title{Two-Level Blocked Individual-level Random Assignment Design}

\description{
For two-level randomized block designs (treatment at level 1, with random effects across level 2 blocks), use \code{mdes.bira2()} to calculate the minimum detectable effect size, \code{power.bira2()} to calculate the statistical power, and \code{mrss.bira2()} to calculate the minimum required sample size.
}

\usage{
mdes.bira2(power=.80, alpha=.05, two.tailed=TRUE,
           rho2,  omega2, p=.50, g2=0, r21=0, r2t2=0,
           n, J)

power.bira2(es=.25, alpha=.05, two.tailed=TRUE,
            rho2,  omega2, g2=0, p=.50, r21=0, r2t2=0,
            n, J)

mrss.bira2(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
           n, J0=10, tol=.10,
           rho2, omega2, g2=0, p=.50, r21=0, r2t2=0)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{es}{effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome between level 2 units (unconditional ICC2).}
  \item{omega2}{treatment effect heterogeneity as ratio of treatment effect variance among level 2 units to the residual variance at level 2.}
  \item{p}{average proportion of level 1 units randomly assigned to treatment within level 2 units.}
  \item{g2}{number of covariates at level 2.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r2t2}{proportion of treatment effect variance among level 2 units explained by level 2 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{level 2 sample size.}
  \item{J0}{starting value for \code{J}.}
  \item{tol}{tolerance to end iterative process for finding \code{J}.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in power calculation.}
  \item{df}{degrees of freedom.}
  \item{ncp}{noncentrality parameter.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{mdes}{minimum detectable effect size.}
  \item{J}{number of level 2 units.}
}

\examples{
# cross-checks
mdes.bira2(rho2=.17, omega2=.50, n=15, J=20)
power.bira2(es=.366, rho2=.17, omega2=.50, n=15, J=20)
mrss.bira2(es=.366, rho2=.17, omega2=.50, n=15)
}
