\name{cra2}
\alias{cra2}
\alias{cra2r2}
\alias{bcra3f2}
\alias{mdes.cra2}
\alias{power.cra2}
\alias{mrss.cra2}
\alias{mdes.cra2r2}
\alias{mdes.mod221}
\alias{mdes.mod222}
\alias{mdesd.mod221}
\alias{mdesd.mod222}
\alias{power.cra2r2}
\alias{power.mod221}
\alias{power.mod222}
\alias{power.med211}
\alias{power.med221}
\alias{mrss.cra2r2}
\alias{mrss.mod221}
\alias{mrss.mod222}
\alias{mdes.bcra3f2}
\alias{power.bcra3f2}
\alias{mrss.bcra3f2}

\title{Two-level Cluster-randomized Trials to Detect Main, Moderation and Mediation Effects}

\description{
For main treatment effects, use \code{mdes.cra2()} to calculate the minimum detectable effect size, \code{power.cra2()} to calculate the statistical power, \code{mrss.cra2()} to calculate the minimum required sample size (number of clusters).

For moderator at level 1, use \code{mdes.mod221()} to calculate the minimum detectable effect size, \code{power.mod221()} to calculate the statistical power, \code{mrss.mod221()} to calculate the minimum required sample size (number of clusters).

For moderator at level 2, use \code{mdes.mod222()} to calculate the minimum detectable effect size, \code{power.mod222()} to calculate the statistical power, \code{mrss.mod222()} to calculate the minimum required sample size (number of clusters).

For mediator at level 1 and level 2, use \code{power.med211()} to calculate the statistical power for the 2-1-1 mediation, and \code{power.med221()} for the 2-2-1 mediation.

For cluster-randomized block designs (treatment at level 2, with fixed effects across level 3 blocks), use \code{mdes.bcra3f2()} to calculate the minimum detectable effect size, \code{power.bcra3f2()} to calculate the statistical power, and \code{mrss.bcra3f2()} to calculate the minimum required sample size (number of clusters per block).
}

\usage{
mdes.cra2(power=.80, alpha=.05, two.tailed=TRUE,
          rho2, p=.50, g2=0, r21=0, r22=0,
          n, J)

mdesd.mod221(power=.80, alpha=.05, two.tailed=TRUE,
            rho2, omegam2, g1=0, r21=0, r2m2=0,
            p=.50, q=NULL, n, J)

mdesd.mod222(power=.80, alpha=.05, two.tailed=TRUE,
            rho2, g2=0, r21=0, r22=0,
            p=.50, q=NULL, n, J)

power.cra2(es=.25, alpha=.05, two.tailed=TRUE,
           rho2, g2=0, p=.50, r21=0, r22=0,
           n, J)

power.mod221(es=.25, alpha=.05, two.tailed=TRUE,
             rho2, omegam2, g1=0, r21=0, r2m2=0,
             p=.50, q=NULL, n, J)

power.mod222(es=.25, alpha=.05, two.tailed=TRUE,
             rho2, g2=0, r21=0, r22=0,
             p=.50, q=NULL, n, J)

power.med211(esa, esb1, esB, escp, two.tailed = TRUE, alpha = .05,
             mc = FALSE, nsims = 1000, ndraws = 1000,
             rhom2, rho2, r21, r22, r2m1, r2m2,
             p, n, J)

power.med221(esa, esb, escp, two.tailed = TRUE, alpha = .05,
             mc = FALSE, nsims = 1000, ndraws = 1000,
             rho2, r22, r21, r2m2,
             p = .50, n, J)

mrss.cra2(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
          n, J0=10, tol=.10,
          rho2, g2=0, p=.50, r21=0, r22=0)

mrss.mod221(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
            n, J0=10, tol=.10, rho2, omegam2, g1=0, r21=0, r2m2=0,
            p=.50, q=NULL)

mrss.mod222(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
            n, J0=10, tol=.10, rho2, g2=0, r21=0, r22=0,
            p=.50, q=NULL)

mrss.mod222(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
            n, J0=10, tol=.10, rho2, g2=0, r21=0, r22=0,
            p=.50, q=NULL)

mdes.bcra3f2(power=.80, alpha=.05, two.tailed=TRUE,
             rho2, p=.50, g2=0, r21=0, r22=0,
             n, J, K)

power.bcra3f2(es=.25, alpha=.05, two.tailed=TRUE,
              rho2, p=.50, g2=0, r21=0, r22=0,
              n, J, K)

mrss.bcra3f2(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
             n, K, J0=10, tol=.10,
             rho2, p=.50, g2=0, r21=0, r22=0)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{es, esa, esb, esb1, esB, escp}{effect size for main/moderator effects, or for path coefficients a (treatment - mediator), b (level 2 mediator - outcome), b1 (level 1 mediator - outcome), B (overall mediator - outcome) or cp (direct treatment - outcome) in the mediation model.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{FALSE} for one-tailed hypothesis testing.}
  \item{rho2}{proportion of variance in the outcome between level 2 units (unconditional ICC2).}
  \item{rhom2}{proportion of variance in the mediator between level 2 units.}
  \item{omegam2}{ratio of the unconditional variance in the moderator effect that is between level 2 units to the residual variance between level 2 units in the null model.}
  \item{p}{proportion of level 2 units randomly assigned to treatment.}
  \item{q}{proportion of level 1 or level 2 units in the moderator subgroup.}
  \item{g1}{number of covariates at level 1.}
  \item{g2}{number of covariates at level 2.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{r2m1}{proportion of mediator variance at level 1 explained by level 1 covariates.}
  \item{r2m2}{proportion of variance in the moderator effect that is explained by level 2 predictors. For the mediation model, proportion of variance in the mediator explained by level 2 predictors.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{level 2 sample size.}
  \item{K}{number of level 3 units (blocks).}
  \item{J0}{starting value for \code{J}.}
  \item{tol}{tolerance to end iterative process for finding \code{J}.}
  \item{mc}{logical; \code{TRUE} for monte carlo simulation based power.}
  \item{nsims}{number of replications, if \code{mc = TRUE}.}
  \item{ndraws}{number of draws from the distribution of the path coefficients for each replication, if \code{mc = TRUE}.}
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
# cross-checks for the main effect
mdes.cra2(rho2=.17, n=15, J=20)
power.cra2(es=.629, rho2=.17, n=15, J=20)
mrss.cra2(es=.629, rho2=.17, n=15)

# cross-checks for the randomly varying cont. L1 moderator effect
mdesd.mod221(rho2=.17, omegam2=.10, n=15, J=20)
power.mod221(es=.3563, rho2=.17, omegam2 =.10, n=15, J=20)
mrss.mod221(es=.3563, rho2=.17, omegam2 =.10, n=15)

# cross-checks for the non-randomly varying cont. L1 moderator effect
mdesd.mod221(rho2=.17, omegam2=0, n=15, J=20)
power.mod221(es=0.2957, rho2=.17, omegam2 =0, n=15, J=20)
mrss.mod221(es=0.2957, rho2=.17, omegam2 =0, n=15)

# cross-checks for the randomly varying bin. L1 moderator effect
mdesd.mod221(rho2=.17, omegam2=.10, q=.50, n=15, J=20)
power.mod221(es=.647, rho2=.17, omegam2 =.10, q=.50, n=15, J=20)
mrss.mod221(es=.647, rho2=.17, omegam2 =.10, q=.50, n=15)

# cross-checks for the non-randomly varying bin. L1 moderator effect
mdesd.mod221(rho2=.17, omegam2=0, q=.50, n=15, J=20)
power.mod221(es=0.5915, rho2=.17, omegam2 =0, q=.50, n=15, J=20)
mrss.mod221(es=0.5915, rho2=.17, omegam2 =0, q=.50, n=15)

# cross-checks for the cont. L2 moderator effect
mdesd.mod222(rho2=.17, n=15, J=100)
power.mod222(es=0.2742, rho2=.17, n=15, J=100)
mrss.mod222(es=0.2742, rho2=.17, n=15)

# cross-checks for the bin. L2 moderator effect
mdesd.mod222(rho2=.17, q=.50, n=15, J=100)
power.mod222(es=0.5485, rho2=.17, q=.50, n=15, J=100)
mrss.mod222(es=0.5485, rho2=.17, q=.50, n=15)

# 2-2-1 mediation
power.med221(esa=0.6596, esb=0.1891, escp=.1,
             rho2=.15, r22=.52, r21=.40, r2m2=.50,
             n=100, J=40, p=.5)

# 2-1-1 mediation
power.med211(esa=0.4135, esb1=0.0670, esB=0.3595, escp=.1,
            rhom2=.3, rho2=.3, r22=.6, r21=.6, r2m2=.6, r2m1=.6,
            n=30, J=80, p=.1)

# cross-checks for cluster-randomized block design
# treatment at level 2, with fixed effects across level 3 blocks
mdes.bcra3f2(rho2=.10, n=20, J=44, K=5)
power.bcra3f2(es = .145, rho2=.10, n=20, J=44, K=5)
mrss.bcra3f2(es = .145, rho2=.10, n=20, K=5)
}
