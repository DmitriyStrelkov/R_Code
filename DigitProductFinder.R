#

Euler40 = function(N){
  
  num_list <- c(1:N)
  number <- paste(num_list, collapse = '')
  number <- substr(number, 1, N)
  number <- as.numeric(strsplit(number, "")[[1]])
  
  answer = 1
  power = 0
  while(10^power < N){
    answer = answer * number[10^power]
    power = power + 1
  }
  return(answer)
}