gcd = function(a,b) {
  remainder = 0
  while (b != 0){
    remainder = b
    b= a %% b
    a = remainder
  }
  return(a)
}

lcm = function(a,b){
  abs((a*b)/(gcd(a,b)))
}

Euler5 = function(maxVal = 20) {
  mult_list = c(1:maxVal)
  lcm_new = 1
  for (mult in mult_list){
    lcm_new = lcm(lcm_new, mult) #lcm(a,b,c) = lcm(a, lcm(b,c))
  }
  return(lcm_new)
}