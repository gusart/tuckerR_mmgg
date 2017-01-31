#'@method print diff
print.diff <- function(mues_diff, ...){
  cat("All posible models","\n")
  print(as.matrix(mues_diff$models))
  cat("Critic Value",mues_diff$critic_value,"\n")
}
