% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{annotate_symbol}
\alias{annotate_symbol}
\title{Add transcript symbol column from ensembl id}
\usage{
annotate_symbol(.data, .ensembl, action = "add")
}
\arguments{
\item{.data}{A `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> |}

\item{.ensembl}{A character string. The column that is represents ensembl gene id}

\item{action}{A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).}
}
\value{
A `tbl` object including additional columns for transcript symbol
}
\description{
annotate_symbol() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with the the additional transcript symbol column
}
\details{
\lifecycle{experimental}

This is useful since different resources use ensembl IDs while others use gene symbol IDs.
}
\examples{
\donttest{



	counts_ensembl \%>\% annotate_symbol(ens)

}

}
