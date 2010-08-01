\name{heatmap.circular}
\alias{heatmap.circular}

\title{
 Draw a Heat Map for circular data
}
\description{
  A heat map is a false color image (basically
  \code{\link{image}(t(x))}) with a dendrogram added to the left side
  and to the top.  Typically, reordering of the rows and columns
  according to some set of values (row or column means) within the
  restrictions imposed by the dendrogram is carried out. See also \code{\link[stats]{heatmap}}.
}
\usage{
heatmap.circular(x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
distfun = dist.circular, hclustfun = hclust, 
reorderfun = function(d, w) reorder(d, w), add.expr, symm = FALSE, 
revC = identical(Colv, "Rowv"), na.rm = TRUE, margins = c(5, 5), 
lwid = c(1, 4), lhei = c(1, 4), ColSideColors, RowSideColors, 
NAColors = "black", cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), 
labRow = NULL, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, 
keep.dendro = FALSE, annotate.expr, annotate = rep(NA, 4), 
verbose = getOption("verbose"), ...)
}
\arguments{
  \item{x}{numeric matrix of class \code{\link{circular}} of the values to be plotted.}
  \item{Rowv}{determines if and how the \emph{row} dendrogram should be
    computed and reordered.  Either a \code{\link{dendrogram}} or a
    vector of values used to reorder the row dendrogram or
    \code{\link{NA}} to suppress any row dendrogram (and reordering) or
    by default, \code{\link{NULL}}, see \sQuote{Details} below.}
  \item{Colv}{determines if and how the \emph{column} dendrogram should be
    reordered.  Has the same options as the \code{Rowv} argument above and
    \emph{additionally} when \code{x} is a square matrix, \code{Colv =
      "Rowv"} means that columns should be treated identically to the
    rows (and so if there is to be no row dendrogram there will not be a
    column one either).}
  \item{distfun}{function used to compute the distance (dissimilarity)
    between both rows and columns.  Defaults to \code{\link{dist.circular}}.}
  \item{hclustfun}{function used to compute the hierarchical clustering
    when \code{Rowv} or \code{Colv} are not dendrograms.  Defaults to
    \code{\link{hclust}}. Should take as argument a result of \code{distfun}
    and return an object to which \code{\link{as.dendrogram}} can be applied.}
  \item{reorderfun}{function(d,w) of dendrogram and weights for
    reordering the row and column dendrograms.  The default uses
    \code{\link{reorder.dendrogram}}.}
  \item{add.expr}{expression that will be evaluated after the call to
    \code{image}.  Can be used to add components to the plot.}
  \item{symm}{logical indicating if \code{x} should be treated
    \bold{symm}etrically; can only be true when \code{x} is a square matrix.}
  \item{revC}{logical indicating if the column order should be
    \code{\link{rev}}ersed for plotting, such that e.g., for the
    symmetric case, the symmetry axis is as usual.}
  \item{na.rm}{logical indicating whether \code{NA}'s should be removed.}
  \item{margins}{numeric vector of length 2 containing the margins
    (see \code{\link{par}(mar= *)}) for column and row names, respectively.}
  \item{lwid}{a vector of values for the widths of columns on the device.
          Relative widths are specified with numeric values. Absolute
          widths (in centimetres) are specified with the \code{\link{lcm}()}
          function (see \code{\link{layout}}).}
  \item{lhei}{a vector of values for the heights of rows on the device.
          Relative and absolute heights can be specified, see \code{lwid} above.}
  \item{ColSideColors}{(optional) character vector of length \code{ncol(x)}
    containing the color names for a horizontal side bar that may be used to
    annotate the columns of \code{x}.}
  \item{RowSideColors}{(optional) character vector of length \code{nrow(x)}
    containing the color names for a vertical side bar that may be used to
    annotate the rows of \code{x}.}
  \item{NAColors}{
    the color used to plot missing values.
  }
  \item{cexRow, cexCol}{positive numbers, used as \code{cex.axis} in
    for the row or column axis labeling.  The defaults currently only
    use number of rows or columns, respectively.}
  \item{labRow, labCol}{character vectors with row and column labels to
    use; these default to \code{rownames(x)} or \code{colnames(x)},
    respectively.}
  \item{main, xlab, ylab}{main, x- and y-axis titles; defaults to none.}
  \item{keep.dendro}{logical indicating if the dendrogram(s) should be
    kept as part of the result (when \code{Rowv} and/or \code{Colv} are
    not NA).}
  \item{annotate}{annotation in the four external side of the figure. A positive value in a position means you want annotate something in that position (1=bottom, 2=left, 3=top, 4=right). For instance, \code{annotate=c(0.1, NA, NA, 1, 1)} means you want to annotate one thing on the bottom with dimension \code{0.1} and two things on right each with dimension \code{1}.}
  \item{annotate.expr}{must be a list of expressions with the same length as \code{annotate}. For instance for \code{annotate=c(0.1, NA, NA, 1, 1)} you must have something as \code{annotate.expr=list(expr1, NA, NA, expr2, expr2)} where \code{expr1} etc. must be a valid \code{R} expression able to produce a plot.}
  \item{verbose}{logical indicating if information should be printed.}
  \item{\dots}{additional arguments passed on to \code{\link{image}},
    e.g., \code{col} specifying the colors.}}
\details{
  If either \code{Rowv} or \code{Colv} are dendrograms they are honored
  (and not reordered).  Otherwise, dendrograms are computed as
  \code{dd <- as.dendrogram(hclustfun(distfun(X)))} where \code{X} is
  either \code{x} or \code{t(x)}.

  If either is a vector (of \sQuote{weights}) then the appropriate
  dendrogram is reordered according to the supplied values subject to
  the constraints imposed by the dendrogram, by \code{\link{reorder}(dd,
    Rowv)}, in the row case.
  If either is missing, as by default, then the ordering of the
  corresponding dendrogram is by the mean direction value of the rows/columns,
  i.e., in the case of rows, \code{Rowv <- rowMeans(x, na.rm=na.rm)}.
  If either is \code{\link{NULL}}, \emph{no reordering} will be done for
  the corresponding side.

  Unless \code{Rowv = NA} (or \code{Colw = NA}), the original rows and
  columns are reordered \emph{in any case} to match the dendrogram,
  e.g., the rows by \code{\link{order.dendrogram}(Rowv)} where
  \code{Rowv} is the (possibly \code{\link{reorder}()}ed) row
  dendrogram.

  \code{heatmap()} uses \code{\link{layout}} and draws the
  \code{\link{image}} in the lower right corner of a 2x2 layout.
  Consequentially, it can \bold{not} be used in a multi column/row
  layout, i.e., when \code{\link{par}(mfrow= *)} or \code{(mfcol= *)}
  has been called.
}

\value{
\code{\link{par}(mfrow= *)} or \code{(mfcol= *)}
  has been called.
}
\author{
Claudio Agostinelli using the code from \code{\link[stats]{heatmap}}.
}

\seealso{\code{\link{dist.circular}}, \code{\link{heatmap}}, \code{\link{image}}, \code{\link{hclust}}}

%\examples{
%
%}

\keyword{hplot}
