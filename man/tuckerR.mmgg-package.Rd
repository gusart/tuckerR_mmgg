\name{tuckerR.mmgg-package}
\alias{tuckerR.mmgg-package}
\alias{tuckerR.mmgg}
\docType{package}
\title{
Three-Mode Principal Components: Tuckers  Models
}
\description{
  This package performs Three-Mode Principal Components  using Tuckers  Models and  plot  interactive  Biplot. The Three-Mode Principal Component Analysis, provides both useful analytic and graphic tools to study and characterize phytogenetic resources, especially when the influence of environmental factors are possible.
}
\details{
\tabular{ll}{
Package: \tab tuckerR.mmgg\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-10-27\cr
License: \tab GPL-3\cr
}
}
\author{
Marta Marticorena, Gustavo Gimenez, Cecilia Gonzalez, Sergio Bramardi
}
\references{
 Marticorena, M.; Bramardi, S.; Defacio, R. (2010). Characterization of
  maize populations in different environmental conditions by means of
  Three-Mode Principal Components Analysis. Revista Ciencia e
  Investigacion   Agraria. 37(3): 93-105.
  
  Timmerman, M.E., and H. Kiers. 2000. Three-mode principal components
 analysis. Choosing numbers of components and sensitivity to local
 optima. The British Journal of the Mahematical and Statistical Psychology 53: 1-16.
}
\keyword{ package }
\examples{
data(maize_pop)
prueba1 <- tucker2R(maize_pop, amb=2, stand=TRUE, nc1=3, nc2=3)
plotear_tucker2(prueba1)
}
