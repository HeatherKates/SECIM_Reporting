% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/construct_database.R
\name{construct_database}
\alias{construct_database}
\title{Construct in-house or public MS2 database for metid.}
\usage{
construct_database(
  path = ".",
  version = "0.0.1",
  metabolite.info.name = "metabolite.info.csv",
  source = "Michael Snyder Lab",
  link = "http://snyderlab.stanford.edu/",
  creater = "Xiaotao Shen",
  email = "shenxt1990@163.com",
  rt = TRUE,
  mz.tol = 15,
  rt.tol = 30,
  threads = 3
)
}
\arguments{
\item{path}{Work directory.}

\item{version}{The version of you database. Default is 0.0.1.}

\item{metabolite.info.name}{The metabolite information table name, it must be csv format.
The demo data can be got from the `demoData` package.
Please see \url{https://tidymass.github.io/metid/articles/metid.html}}

\item{source}{The source of your database.}

\item{link}{Website link of the source.}

\item{creater}{Creater name. For example, Xiaotao Shen.}

\item{email}{email address.}

\item{rt}{Do the metabolites have RT information or not?. If not, set it as FALSE.}

\item{mz.tol}{m/z tolerance for the match between metabolites and precursor m/z of MS2 spectra.}

\item{rt.tol}{RT tolerance for the match between metabolites and precursor m/z of MS2 spectra.}

\item{threads}{The number of threads}
}
\value{
A databaseClass object.
}
\description{
Construct MS2 spectra database according to mzXML data and compound information table (csv format).
}
\examples{
\dontrun{
database0.0.1 <- construct_database(
path = ".",
version = "0.0.1",
metabolite.info.name = "metabolite.info.csv",
creater = "dumine",
email = "dumine@zju.edu.cn",
rt = FALSE,
mz.tol = 15,
rt.tol = 30,
threads = 5
)
}
}
\seealso{
The example and demo data of this function can be found
\url{https://tidymass.github.io/metid/articles/metid.html}
}
\author{
Xiaotao Shen
\email{shenxt1990@outlook.com}
}
