\name{ira}
\alias{ira}
\alias{ira1r1}
\alias{bira2f1}
\alias{bira2c1}
\alias{mdes.ira}
\alias{power.ira}
\alias{mrss.ira}
\alias{mdes.ira1r1}
\alias{power.ira1r1}
\alias{mrss.ira1r1}
\alias{mdes.bira2f1}
\alias{power.bira2f1}
\alias{mrss.bira2f1}
\alias{mdes.bira2c1}
\alias{power.bira2c1}
\alias{mrss.bira2c1}
\title{Individual-level Random Assignment Design}

\description{
For an individual-level random assignment design (simple randomized controlled trial) use \code{mdes.ira()} to calculate the minimum detectable effect size, \code{power.ira()} to calculate the statistical power, and \code{mrss.ira()} to calculate the minimum required sample size (number of subjects).

If fixed block effects exist (randomized block designs), that is, treatment effect varies from one block to another but it is not random, use \code{mdes.bira2f1()} to calculate the minimum detectable effect size, \code{power.bira2f1()} to calculate the statistical power, and \code{mrss.bira2f1()} to calculate the minimum required sample size (number of subjects per block).

If block effect is constant, that is, treatment effect does not change from one block to another, use \code{mdes.bira2c1()} to calculate the minimum detectable effect size, \code{power.bira2c1()} to calculate the statistical power, and \code{mrss.bira2c1()} to calculate the minimum required sample size (number of subjects per block).
}

\usage{
mdes.ira(power=.80, alpha=.05, two.tailed=TRUE,
         p=.50, g1=0, r21=0, n)

power.ira(es=.25, alpha=.05, two.tailed=TRUE,
          p=.50, g1=0, r21=0, n)

mrss.ira(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
         n0=10, tol=.10,
         p=.50, g1=0, r21=0)

mdes.bira2f1(power=.80, alpha=.05, two.tailed=TRUE,
             p=.50, g1=0, r21=0, n, J)

power.bira2f1(es=.25, alpha=.05, two.tailed=TRUE,
              p=.50, g1=0, r21=0, n, J)

mrss.bira2f1(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
             J, n0=10, tol=.10,
             p=.50, g1=0, r21=0)

mdes.bira2c1(power=.80, alpha=.05, two.tailed=TRUE,
             p=.50, g1=0, r21=0,
             n, J)

power.bira2c1(es=.25, alpha=.05, two.tailed=TRUE,
              p=.50, g1=0, r21=0,
              n, J)

mrss.bira2c1(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
             J, n0=10, tol=.10,
             p=.50, g1=0, r21=0)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{es}{effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{p}{proportion of units randomly assigned to the treatment condition.}
  \item{g1}{number of covariates.}
  \item{r21}{proportion of variance in the outcome explained by covariates.}
  \item{n}{sample size (number of subjects per block, if exist, calculated using simple average or harmonic mean)}
  \item{J}{level 2 sample size (number of blocks).}
  \item{n0}{starting value for \code{n}.}
  \item{tol}{tolerance to end iterative process for finding \code{n} or \code{J.}}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in power calculation.}
  \item{df}{degrees of freedom.}
  \item{ncp}{noncentrality parameter.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{mdes}{minimum detectable effect size.}
  \item{n}{sample size (number of subjects).}
}

\examples{
# cross-checks
mdes.ira(n=250)
power.ira(es=.356, n=250)
mrss.ira(es=.356)
# cross-checks
mdes.bira2c1(n=15, J=20)
power.bira2c1(es=.325, n=15, J=20)
mrss.bira2c1(es=.325, J=20)
# cross-checks
mdes.bira2f1(n=15, J=20)
power.bira2f1(es=.325, n=15, J=20)
mrss.bira2f1(es=.325, J=20)
}


