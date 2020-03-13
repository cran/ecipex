\name{ecipex}
\alias{ecipex}
\title{Calculate fine structure isotope pattern}
\usage{
  ecipex(formulas, isoinfo = ecipex::nistiso, limit = 1e-12,
    id = FALSE, sortby = "abundance", gross = FALSE, groupby = "mass")
}
\arguments{
  \item{formulas}{a character vector specifying the
  chemical formulas whose isotope patterns are to be
  calculated. The elements specified must be present in
  \code{isoinfo} and must be compatible with
  \code{\link{count.elements}}.}

  \item{isoinfo}{a data frame that specifies the masses and
  isotopic abundances of the elements to be used in the
  calculations. Must include the variables \code{element}
  (the chemical symbol), \code{mass} and \code{abundance}
  so that each row specifies a unique isotope. Defaults to
  \code{\link{nistiso}}.}

  \item{limit}{isotopologues with abundances below this
  value are ignored.}

  \item{id}{determines whether the full isotopic
  composition of each isotopic variant should be specified
  in the output. Setting this to \code{TRUE} can increase
  memory requirements substantially. This argument has no
  effect if \code{gross} is set to \code{TRUE}.}

  \item{sortby}{should be one of \code{"abundance"} or
  \code{"mass"}, depending on whether the output should be
  sorted by abundance (high to low) or mass (low to high).}

  \item{gross}{if \code{TRUE} the gross isotope
  patterns will be returned instead of the fine structure ones.
  }

  \item{groupby}{determines how the fine structure
  isotope patterns are grouped together when the gross structure
  isotope patterns are calculated. If it is equal to \code{"mass"}
  then isotopologues whose masses get rounded to the same integer number are 
  grouped together. If it is equal to \code{"nucleons"} then 
  isotopologues with the same number of nucleons get grouped 
  together.  This argument is only used if \code{gross} 
  is set to \code{TRUE}.

  }
}
\value{
  A list of data frames containing the fine structure
  isotope patterns of each formula in \code{formulas}. The
  list names are determined by \code{formulas}, so that the
  isotope pattern of any given formula is easily extracted.
  If \code{id} is set to \code{TRUE}, the output will
  include additional columns listing the counts of each
  distinct isotope of each element. If \code{gross}
  is set to \code{TRUE}, the output will instead contain the
  gross structure isotope patterns of each formula in \code{formulas}.
}
\description{
  This function can calculate either the fine structure isotope
  patterns or the gross structure isotope patterns for a given set 
  of chemical formulas. It returns the isotope patterns as a list 
  of data frames which may be sorted by mass or abundance.
}
\details{
  The fine structure isotope pattern of each formula is calculated by
  applying the multi-dimensional fast Fourier transform to
  a simplex-based representation of each element's isotopic
  abundance. The algorithm is most efficient when the atom
  counts of the same elements in different formulas are
  roughly similar. Performance can also be improved by
  increasing \code{limit} although this should be done with
  care. It is generally not advisable to reduce
  \code{limit} below its default value, as this can
  increase memory requirements significantly and is
  unlikely to provide information of any value, since the
  natural isotopic abundances can be quite variable. If \code{gross} 
  is set to \code{TRUE} then the gross structure isotope patterns 
  are calculated directly from the fine structure isotope patterns.
  Note that for centroids with extremely low abundances (say less than
  one billionth of the total) the centroided masses
  can be somewhat inaccurate due to floating point errors.
  
}
\examples{
# a simple molecule
iso_H2O <- ecipex("H2O")[[1]]
iso_H2O

# reduce limit for larger molecule and sort output by mass
iso_C254H338N65O75S6 <- ecipex("C254H338N65O75S6", limit=1e-8, sortby="mass")[[1]]
head(iso_C254H338N65O75S6)

# check that sum of all abundances is still close to 1
sum(iso_C254H338N65O75S6$abundance)

# inspect the full isotope pattern, the fine structure, and the full pattern on a log scale
par(mfrow=c(1,3))
plot(iso_C254H338N65O75S6, t="h")
plot(iso_C254H338N65O75S6, t="h", xlim=c(5691.29, 5691.31))
plot(iso_C254H338N65O75S6, t="h", log="y")

# calculate isotopic abundances with enriched Carbon-13
modifiediso <- nistiso 
modifiediso[modifiediso$element=="C",3] <- c(0.9, 0.1)
ecipex("C2", isoinfo=modifiediso)

# the isotope pattern can be calculated quickly if the elements only have 2 stable isotopes
system.time(iso_C10000H10000 <- ecipex("C10000H10000", limit=1e-8)[[1]])

# this is typically a more demanding calculation because S has 4 stable isotopes
system.time(iso_S50 <- ecipex("S50", limit=1e-8)[[1]])

# if the limit is greater than the most abundant isotopologue the output is uninformative
iso_C10000H10000_useless <- ecipex("C10000H10000", limit=0.015)

# calculate the isotope patterns of multiple formulas, and include the detailed isotopic composition
multisopatterns <- ecipex(c("H2O", "CO2", "O2", "C8H18", "C60"), sortby="mass", id=TRUE)

# inspect C8H18 in particular
multisopatterns$C8H18

# make sure all abundances are close to 1
sapply(multisopatterns, function(x){sum(x$abundance)})

# due to floating point errors, the following are not identical
iso_C60_almostComplete <- ecipex("C60", limit= 0)[[1]]
iso_C60_reallyComplete <- ecipex("C60", limit= -1)[[1]]
# the latter includes negative isotopic abundances because the floating point errors are orders of 
# magnitude greater than the "true" abundances. The variations in natural isotopic abundances will
# typically be much greater than floating point errors.

# calculate the gross structure isotope pattern, grouping the fine structure isotopologues by mass
ecipex("C6H14N4O2", sortby="mass", gross=TRUE, groupby="mass")
}
\references{
Ipsen, A., Efficient Calculation of Exact Fine Structure Isotope Patterns via the Multidimensional Fourier Transform, Anal. Chem., 2014, 86 (11), pp 5316-5322

\url{http://pubs.acs.org/doi/abs/10.1021/ac500108n}
}
\seealso{
  \code{\link{nistiso}}
}
\keyword{fft}
\keyword{fine}
\keyword{isotope}
\keyword{pattern}
\keyword{structure}

