prob_bus <- function(mins, prob){
  return((1 - prob)^mins)
}

prob_bday <- function(people){
  list <- c(365:(365-people+1))
  return(1 - (1/365)^people * prod(list))
}

expLikelihood <- function(tauList, waitTimes){
    return((1/tauList^length(waitTimes) * exp(-1/tauList * sum(waitTimes))))
}