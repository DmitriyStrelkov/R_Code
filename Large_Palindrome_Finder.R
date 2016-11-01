Euler4 = function(digits = 3) {
  #First, we need to create a couple of functions that allow us to check if the number is a palindrome
  #These two functions are taken from the following source http://bit.ly/2cOiumf
  #The first function, strrev, takes a value, converts it to a string and reverses it
  strrev = function(x){
    paste(substring(x, nchar(x):1, nchar(x):1),
          collapse = "")
  }
  
  #The second function, palindrome, takes a value, converts it to a string, then compares it to its reverse
  palindrome = function(x){
    x[1] == strrev(x[1])
  }
  
  #Create an empty list of palindromes
  dromelist = c(0)
  
  #Then, we multiply every three digit value by every other three digit value
  #and check if the product is a palindrome
  #This method is more computationally expensive than other methods, but it creates a list of all
  #palindromes that fit the criteria
  for(x in 100:999){
    for(y in 100:999){
      if(palindrome(x * y) & (x*y) > dromelist[1]){
        dromelist = c(x*y)
      }
    }
  }
  #We sort the list of palindromes in ascending order and then take the last value in the list
  return(dromelist)
  
  #Answer checks out in Project Euler!
}