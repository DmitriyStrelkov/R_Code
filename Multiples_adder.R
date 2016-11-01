#Written in R by Dmitriy Strelkov
#Find the sum of all numbers that are multiples of 3 or 5 (no duplicates) and are less than 1000

Euler1 = function(maxN = 1000) {
  Threes = seq(0,(maxN-1), by=3)
  Fives = seq(0,(maxN-1), by=5)
  Duplicates = seq(0,-(maxN-1), by=-15)
  
  #Three lists are created: One for threes, one for fives, and one for the duplicates(multiples of 15)
  #We only go up to 999 because the questions asks for the sum of numbers less than 1000
  #The duplicates are negative to simplify every operation into a summation
  #Next, we sum the elements of each vector, collapsing each vector into an integer
  
  Threes_total = sum(Threes)
  Fives_total = sum(Fives)
  Duplicates_total = sum(Duplicates)
  
  #Finally we sum all of the resultant integers together to get our answer
  
  Answer = Threes_total + Fives_total + Duplicates_total
  return(Answer)
  
  #Answer checks out in ProjectEuler!
}