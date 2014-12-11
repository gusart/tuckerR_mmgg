\name{plotear_tucker2}
\alias{plotear_tucker2}
\title{Plot an interactive Biplot
}
\description{The interactive Biplot consists of combining two of the modes, obtaining markers for the individuals  and vectors for the
variables that were concatenated with the conditions. To plot the interactive Biplot this function need the output for the tucker2R function.
}
\usage{
plotear_tucker2(saltuck)
}
\arguments{
  \item{saltuck}{is a list with the results of the algorithm to plot the
    biplot, where the names of the conditions are well kept.
  }
}
\details{The interactive Biplot allows for the visualization of the
  inter-structure of the differents data tables.
}
\references{
 Marticorena, M.; Bramardi, S.; Defacio, R. (2010). Characterization of
  maize populations in different environmental conditions by means of
  Three-Mode Principal Components Analysis. Revista Ciencia e
  Investigacion   Agraria. 37(3): 93-105.

  Gabriel, K.R. 1971. The biplot graphic display of matrices with appications to principal components analysis. Biometrika. 58, 453-467.
}
\author{Marta Marticorena, Gustavo Gimenez, Cecilia Gonzalez, Sergio Bramardi
}
\examples{
data(maize_pop)
prueba1 <- tucker2R(maize_pop, amb=2, stand=TRUE, nc1=3, nc2=3)
plotear_tucker2(prueba1)
}
\encoding{UTF-8}