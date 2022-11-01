\name{EM.PScr}
\alias{EM.PScr}
\title{
Maximum likelihood estimation based on EM algorithm for the Power Series cure rate model
}
\description{
This function provides the maximum likelihood estimation based on the EM algorithm for the
Power Series cure rate model
}
\usage{
EM.PScr(t, delta, z, model = 1, dist = 1, max.iter = 1000, 
           prec = 1e-04)
}
\arguments{
  \item{t}{
observed times}
  \item{delta}{
failure indicators}
  \item{z}{
matrix of covariates (with n rows and r columns)}
  \item{model}{
distribution to be used for the concurrent causes: 
1 for Poisson, 
2 for logarithmic,
3 for negative binomial, 
4 for bernoulli and 
5 for polylogarithm (Gallardo et al. 2018).
}
  \item{dist}{
distribution to be used for the time-to-event: 
1 for slash half-normal (Gallardo et al., 2022), 
2 for Weibull,
3 for gamma and 
4 for Birnbaum-Saunders.}
  \item{max.iter}{
maximum number of iterations to be used by the algorithm}
  \item{prec}{
precision (in absolute value) for the parameters to stop the algorithm.}
}
\details{
The EM algorithm for the model is implemented as in Gallardo et al. (2017).
}
\value{
\item{estimate}{a matrix containing the estimated parameters and their standard error}
\item{loglike}{the estimated log-likelihood function evaluated in the maximum likelihood
	estimators}
\item{AIC}{the Akaike information criterion}
\item{BIC}{the Bayesian (also known as Schwarz) information criterion}
}
\references{
Gallardo DI, Gomez YM and De Castro M. (2018). A flexible cure rate model based on 
the polylogarithm distribution. Journal of Statistical Computation and Simulation 
88 (11), 2137-2149

Gallardo DI, Gomez YM, Gomez HJ, Gallardo-Nelson MJ, Bourguignon M. (2022) The slash 
half-normal distribution applied to a cure rate model with application to bone marrow 
transplantation. Mathematics, Submitted.

Gallardo DI, Romeo JS and Meyer R. (2017). A simplified estimation procedure based 
on the EM algorithm for the power series cure rate model. Communications in 
Statistics-Simulation and Computation 46 (8), 6342-6359.
}
\author{
Diego I. Gallardo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\donttest{
require(mstate)
data(ebmt4)
attach(ebmt4)
t = srv / 365.25 # Time in years
delta=srv.s
prophy=as.factor(proph)
year2=ifelse(year=="1985-1989",0,1)
z=t(model.matrix(~proph-1))
#Computes the estimation for Poisson-Slash half-normal cure rate model
EM.PScr(t, delta, z, model=1, dist=1)
}
}

