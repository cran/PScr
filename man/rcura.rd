\name{rcura}
\alias{rcura}
\title{
Simulation of Right-Censored Data from a Power Series Cure Rate Model
}
\description{
This function generates right-censored survival data from a power series cure rate model. The number of latent competing causes follows a power series distribution (Poisson, Bernoulli, negative binomial, geometric or polylogarithm) and the time-to-event of each cause follows a Weibull, gamma or Birnbaum-Saunders distribution. The generated data can be analysed with \code{\link{EM.PScr}}.
}
\usage{
rcura(n, beta, z = NULL, q = NULL, alpha, sigma, C, model, dist)
}
\arguments{
\item{n}{
sample size (a positive integer).}
\item{beta}{
numeric vector of regression coefficients (of length \eqn{r}{r}) related to the power parameter \eqn{\theta}{theta} of the distribution of the latent causes. Its length must be equal to the number of rows of \code{z}.}
\item{z}{
matrix of covariates with \eqn{r}{r} rows (one per covariate, the first one usually being the intercept) and \code{n} columns (one per individual). If \code{NULL} (default), an intercept row of ones is used and, if \code{length(beta) > 1}, the remaining covariates are generated as independent Bernoulli(0.5) variables.}
\item{q}{
positive numeric parameter of the distribution of the latent causes. It is required for \code{model = "binneg"} (size parameter of the negative binomial distribution) and for \code{model = "polilog"} (parameter of the polylogarithm distribution). It is ignored for the other models. Default is \code{NULL}.}
\item{alpha}{
positive shape parameter of the distribution of the time-to-event of each latent cause.}
\item{sigma}{
positive scale parameter of the distribution of the time-to-event of each latent cause.}
\item{C}{
positive censoring time, common to all individuals.}
\item{model}{
character string with the distribution of the latent causes: \code{"poisson"}, \code{"bernoulli"}, \code{"binneg"} (negative binomial), \code{"geometrica"} (geometric) or \code{"polilog"} (polylogarithm).}
\item{dist}{
character string with the distribution of the time-to-event of each latent cause: \code{"weibull"}, \code{"gamma"} or \code{"bs"} (Birnbaum-Saunders).}
}
\details{
Let \eqn{M_i}{M_i} be the number of latent causes of individual \eqn{i}{i}, \eqn{i = 1, \ldots, n}{i = 1, ..., n}, with power parameter \eqn{\theta_i = g(z_i^{\top}\beta)}{theta_i = g(z_i' beta)}. Given \eqn{M_i > 0}{M_i > 0}, the time-to-event is \eqn{T_i = \min(W_{i1}, \ldots, W_{iM_i})}{T_i = min(W_i1, ..., W_iM_i)}, where the \eqn{W_{ij}}{W_ij} are independent and identically distributed with the distribution selected in \code{dist}. Individuals with \eqn{M_i = 0}{M_i = 0} are cured and have \eqn{T_i = \infty}{T_i = Inf}. The observed time is \eqn{t_i = \min(T_i, C)}{t_i = min(T_i, C)} and the failure indicator is \eqn{\delta_i = I(T_i \leq C)}{delta_i = I(T_i <= C)}.

The distributions available for the latent causes, the link functions and the corresponding cure fractions \eqn{p_0 = P(M = 0)}{p0 = P(M = 0)} are:

\tabular{llll}{
\code{model} \tab distribution of \eqn{M}{M} \tab link \eqn{g}{g} \tab cure fraction \eqn{p_0}{p0} \cr
\code{"poisson"} \tab Poisson(\eqn{\theta}{theta}) \tab exponential \tab \eqn{\exp(-\theta)}{exp(-theta)} \cr
\code{"bernoulli"} \tab Bernoulli(\eqn{\theta/(1+\theta)}{theta/(1 + theta)}) \tab exponential \tab \eqn{1/(1+\theta)}{1/(1 + theta)} \cr
\code{"binneg"} \tab negative binomial(\eqn{q}{q}, \eqn{1-\theta}{1 - theta}) \tab logit \tab \eqn{(1-\theta)^q}{(1 - theta)^q} \cr
\code{"geometrica"} \tab geometric (\code{"binneg"} with \eqn{q = 1}{q = 1}) \tab logit \tab \eqn{1-\theta}{1 - theta} \cr
\code{"polilog"} \tab \eqn{X - 1}{X - 1}, with \eqn{P(X = m) = \theta^m / (m^q \mathrm{Li}_q(\theta))}{P(X = m) = theta^m / (m^q Li_q(theta))}, \eqn{m \geq 1}{m >= 1} \tab logit \tab \eqn{\theta / \mathrm{Li}_q(\theta)}{theta / Li_q(theta)} \cr
}

The time-to-event distributions are parameterized by a shape \code{alpha} and a scale \code{sigma}: \code{"weibull"} uses \code{\link[stats]{rweibull}}, \code{"gamma"} uses \code{\link[stats]{rgamma}} and \code{"bs"} uses \code{\link[VGAM]{rbisa}}. This is the same parameterization used by \code{\link{EM.PScr}}, so the simulated data can be fitted with the matching \code{dist} code.

The cure fraction and the percentage of censoring of a simulated sample depend on \code{beta}, \code{alpha}, \code{sigma} and \code{C}. If the time-to-event distribution has a heavy tail relative to \code{C}, a fraction of the susceptible individuals is censored before the event occurs, so the Kaplan-Meier curve does not level off at \eqn{p_0}{p0} and the cure fraction is difficult to identify. The empirical cure fraction and the proportion of censored observations of each sample are returned as attributes.
}
\value{
A \code{data.frame} with \code{n} rows and two columns:
\item{t}{observed times, \eqn{\min(T, C)}{min(T, C)}.}
\item{delta}{failure indicators: 1 if the event was observed and 0 if the observation was censored.}

The data frame has the following attributes, which can be extracted with \code{\link[base]{attr}}:
\item{z}{the matrix of covariates used to generate the data.}
\item{theta}{numeric vector with the power parameter \eqn{\theta_i}{theta_i} of each individual.}
\item{M}{integer vector with the simulated number of latent causes of each individual.}
\item{p0_empirico}{empirical cure fraction, that is, the proportion of individuals with \code{M = 0}.}
\item{cens_pct}{proportion of censored observations.}
}
\references{
Birnbaum ZW and Saunders SC. (1969). A new family of life distributions. Journal of Applied Probability 6 (2), 319-327.

Cancho VG, Louzada F and Ortega EMM. (2013). The power series cure rate model: an application to a cutaneous melanoma data. Communications in Statistics - Simulation and Computation 42 (3), 586-602. \doi{10.1080/03610918.2011.639971}

Gallardo DI, Gomez YM and De Castro M. (2018). A flexible cure rate model based on the polylogarithm distribution. Journal of Statistical Computation and Simulation 88 (11), 2137-2149.

Gallardo DI, Romeo JS and Meyer R. (2017). A simplified estimation procedure based on the EM algorithm for the power series cure rate model. Communications in Statistics - Simulation and Computation 46 (8), 6342-6359. \doi{10.1080/03610918.2016.1202276}

Yakovlev AY and Tsodikov AD. (1996). Stochastic Models of Tumor Latency and Their Biostatistical Applications. World Scientific, Singapore.
}
\author{
Daniel Jana, Yolanda Gomez and Diego Gallardo
}
\note{
The function uses the random number generator, so results are reproducible only after calling \code{\link[base]{set.seed}}. The computing time grows linearly with \code{n}.

Correspondence with the codes of \code{\link{EM.PScr}}: \code{model = "poisson"} corresponds to \code{model = 1}, \code{"binneg"} and \code{"geometrica"} to \code{model = 3} (with \eqn{q = 1}{q = 1} in the geometric case), \code{"bernoulli"} to \code{model = 4} and \code{"polilog"} to \code{model = 5}; \code{dist = "weibull"}, \code{"gamma"} and \code{"bs"} correspond to \code{dist = 2}, \code{3} and \code{4}, respectively. The logarithmic and Flory-Schulz models and the slash half-normal and log-normal time-to-event distributions available in \code{\link{EM.PScr}} are not implemented in \code{rcura}.
}
\seealso{
\code{\link{EM.PScr}}
}
\examples{
set.seed(2026)

## Poisson cure rate model with Weibull times and no covariates.
## With theta = -log(0.7) the cure fraction is p0 = exp(-theta) = 0.7
theta <- -log(0.7)
dat <- rcura(n = 500, beta = log(theta), alpha = 2, sigma = 3,
             C = 8, model = "poisson", dist = "weibull")
head(dat[dat$delta == 1, ])  # first observed events
attr(dat, "p0_empirico")     # empirical cure fraction (approximately 0.7)
attr(dat, "cens_pct")        # proportion of censored observations

## The Kaplan-Meier curve levels off near the cure fraction
fit <- survival::survfit(survival::Surv(t, delta) ~ 1, data = dat)
plot(fit, conf.int = FALSE, xlab = "Time", ylab = "Survival probability")
abline(h = 0.7, lty = 2)

## Negative binomial model (q = 2) with gamma times and cure fraction 0.4.
## Since p0 = (1 - theta)^q, we have theta = 1 - p0^(1 / q)
q <- 2
theta <- 1 - 0.4^(1 / q)
dat2 <- rcura(n = 200, beta = qlogis(theta), q = q, alpha = 2, sigma = 1.4,
              C = 8, model = "binneg", dist = "gamma")
mean(dat2$delta == 0)

## Intercept and a binary covariate (covariates in rows, individuals in columns)
n <- 200
z <- rbind(1, rbinom(n, 1, 0.5))
dat3 <- rcura(n = n, beta = c(-0.5, 0.7), z = z, alpha = 0.5, sigma = 2.5,
              C = 8, model = "poisson", dist = "bs")
table(dat3$delta)

\donttest{
## Fit the first simulated data set with EM.PScr (model = 1: Poisson, dist = 2: Weibull)
fit.em <- EM.PScr(t = dat$t, delta = dat$delta,
                  z = matrix(1, nrow = 1, ncol = nrow(dat)), model = 1, dist = 2)
fit.em$estimate

## Polylogarithm model
dat4 <- rcura(n = 200, beta = 0.5, q = 2, alpha = 2, sigma = 3,
              C = 8, model = "polilog", dist = "weibull")
attr(dat4, "p0_empirico")
}
}
\keyword{datagen}
\keyword{survival}
