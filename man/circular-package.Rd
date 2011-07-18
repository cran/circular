\name{Circular}
\alias{Circular}
\alias{circular-package}
\title{Package `circular': summary information}
\description{
  The package `circular' provides functions for the statistical analysis
  and graphics representation of circular data (observations which are
  angles). It originally started as a porting from S-plus to R of
  functions developed for the book: Circular Statistics, from "Topics in
  circular Statistics" (2001) S. Rao Jammalamadaka and A. SenGupta,
  World Scientific. Now, it has an S3 implementation and several new
  functions and datasets.  
}

\section{Version}{
The version level of the package is given by the command 
\code{packageDescription("circular")}.
The most recent version of the package can be obtained from
the R-Forge repository at \url{https://r-forge.r-project.org/projects/circular/}
}

\section{Author}{
Claudio Agostinelli, Department of Environmental Sciences, Informatics
and Statistics, Ca' Foscari University, Venice, Italy
(\url{http://www.dst.unive.it/~claudio})

Ulric Lund, Department of Statistics, California Polytechnic State
University, San Luis Obispo, California, USA
(\url{http://statweb.calpoly.edu/ulund/})
}

\section{Licence}{
This package and its documentation are usable under the terms of the 
"GNU General Public License", a copy of which is distributed with the 
package. While the software is freely usable, it would be appreciated
if a reference is inserted in publications or other work
which makes use of it; for this purpose, see the command
\code{citation("circular")}. 
}
\section{Acknowledgements}{
The package has evolved through several versions, developed over some years.

Many thanks to all that points out bugs, provide suggestions and comments.

The functions \code{watson.wiliams.test} and \code{wallraff.test} are
developed by Jean-Olivier Irisson (\url{http://jo.irisson.com/})

The functions \code{dcarthwrite}, \code{dgenvonmises},
\code{(d,r)katojones}, \code{djonespewsey} are developed by Federico
Rotolo (\url{http://www.stat.unipd.it/cv/federico_rotolo})

The function \code{rose.diag} has contribution by Hiroyoshi Arai (\url{mailto:h_arai@aoni.waseda.jp})

The function \code{windrose} is developed by Matthew Pocernich \url{http://www.rap.ucar.edu/staff/pocernic/}

Dataset \code{swallows} is kindly provided by Dimitri Giunchi \url{http://unimap.unipi.it/cercapersone/dettaglio.php?ri=2504&template=dettaglio.tpl}

If I miss to report your contribution please let me know by email at \url{mailto:claudio@unive.it}
}

\keyword{circle}
\keyword{circular}
\keyword{distribution}
\keyword{univariate}
\keyword{hplot}
\keyword{htest}
