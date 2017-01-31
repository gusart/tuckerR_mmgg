#'@method print marta
print.marta <- function(mostrar, ...){
  cat("Means by variables","\n")
  print(as.matrix(mostrar$Resultado$MEANS))
  cat("Environment:",mostrar$Ambientes,"\n")
  cat("Cases:",mostrar$Resultados$individuos,"\n")
  cat("The 'G' matrix:","\n")
  print(mostrar$matrizG)
  cat("Total Variability:",mostrar$Resultados$vartot,"\n")
  cat("iteration number:",mostrar$Resultados$iteraciones,"\n")
  cat("Explain variability:",mostrar$Resultados$SCExplicada,"\n")
}
