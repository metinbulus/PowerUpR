\name{replication}
\alias{rep}
\alias{mdh}
\alias{q.test}
\alias{power.rep}
\alias{mdh.rep}
\alias{replication}
\title{Unambiguous Test of Replication for Ensemble of Studies}

\description{
Use \code{power.rep()} for the statistical power, and \code{mdh.rep()} for the minimum detectable heterogeneity. Functions implement methods designed to conduct unambiguous test of replication for ensemble of studies (Hedges & Schauer, 2019).
}

\usage{
power.rep(k = 2L, mdh = 1/4, mdh.null = 0, alpha = .05)

mdh.rep(k = 2L, mdh.max = 15, alpha = .05, power = 0.80,
        mdh.null = 0, step = .001, plot = FALSE)

q.test(es, sse, mdh.null = 0, alpha = .05)
}

\arguments{
  \item{k}{number of replications.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{alpha}{probability of type I error.}
  \item{mdh}{minimum detectable heterogeneity (MDH).}
  \item{mdh.null}{MDH for null hypothesis.}
  \item{mdh.max}{maximum of possible MDH values for grid search.}
  \item{step}{step size to generate possible MDH values.}
  \item{plot}{logical; if TRUE plots MDH - power curve.}
  \item{es}{vector; effect sizes for k studies.}
  \item{sse}{vector; standardized standard errors for k studies.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in the calculation.}
  \item{df}{degrees of freedom.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{mdh}{minimum detectable heterogeneity (MDH).}
  \item{Q}{Q statistics.}
  \item{df}{degrees of freedom.}
  \item{p.value}{p values.}
}

\examples{
# cross-checks
power.rep(k=20L, mdh = 0.50)
mdh.rep(k=20L, power = .39)

# observed power
# random sample of 10 studies from Hedges & Schauer (2019) example
ES <- c(0.52, 0.62, 0.62, 0.42, 0.96, 0.11, 0.34, 1.11, 0.79, 0.35)
SSE <- c(0.22, 0.07, 0.22, 0.21, 0.20, 0.22, 0.23, 0.24, 0.21, 0.19)
q <- q.test(es = ES, sse = SSE)
power.rep(mdh = q$Q / q$df) # exact
power.rep(mdh = q$Q / q$df, mdh.null = 1/4) # approximate
}

\references{
Hedges, L. V., & Schauer, J. (2018). Statistical analyses for studying replication: Meta-analytic perspectives. Psychological Methods, 24(5), 557-570. http://dx.doi.org/10.1037/met0000189
}
