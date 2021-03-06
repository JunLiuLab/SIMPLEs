% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_ztrunc.R
\name{ztruncnorm}
\alias{ztruncnorm}
\title{Fit a zero-inflated censored normal distribution.}
\usage{
ztruncnorm(
  dat,
  cutoff = 0.1,
  iter = 200,
  s_upper = 4,
  s_lower = 0.1,
  p_min = 0.6,
  p_max = 1
)
}
\arguments{
\item{dat}{Data matrix.}

\item{cutoff}{The value below cutoff is treated as no expression. Default =
0.1.}

\item{iter}{Max number of iteration for optimization. Default = 100.}

\item{s_upper}{upper bound of standard deviation. Default = 4.}

\item{s_lower}{lower bound of standard deviation. Default = 0.1.}

\item{p_min}{lower bound of the probability of the normal component. Default
= 0.6.}

\item{p_max}{upper bound of the probability of the normal component. Can be either a scaler or a vector. If it is a vector, each element is the upper bound for a gene. Default = 1.}
}
\value{
\code{ztruncnorm} returns a list of results in the following order.
\enumerate{
  \item{param}{A matrix with mean and standard deviation for each gene. }
  \item{pg}{A vector for the probablity of the normal component.}
}
}
\description{
Fit each gene expression by a zero-inflated censored normal distribution by
direct maximize the log-likelihood and return MLE of mean, variance and
probability of mixture component.
}
\details{
It will return NA if number of nonzero is less than 4 indicating that the gene has no expression
or if the optimization is not convergence or error
which will be filled with sample mean and std, p_max afterwards.
}
\author{
Zhirui Hu, \email{zhiruihu@g.harvard.edu}

Songpeng Zu, \email{songpengzu@g.harvard.edu}
}
