#Written in R by Dmitriy Strelkov

Euler2 = function(maxVal = 4000000){
  #We create 3 variable
  #The first counts our fibonacci list index and allows us to add new values to the list
  #The second is the beginning of our list of fibonacci numbers
  #The third is where we will be adding our fibonacci numbers
  i = 3
  fib_list = c(1,1)
  fib_sum = 0
  
  #This while loop creates a new list element that is the sum of the previous two elements
  #As long as the last value created is lower than the max value
  #This will create a list of all fibonacci numbers below the max value and one above
  while(fib_list[i-1]<= maxVal) {
    fib_list[i] = fib_list[i-1] + fib_list[i-2]
    i = i + 1
  }
  
  #Then, for every item in the list we check to make sure that the value is below the max value
  #And we also check to make sure that the value is divisible by 2(even)
  #If both are true, then we add that value to our sum
  for(n in fib_list) {
    if((n %% 2) == 0 && n < maxVal){
      fib_sum = fib_sum + n
    }
  }
  
  return(fib_sum)
  #Answer checks out in Project Euler!
}
