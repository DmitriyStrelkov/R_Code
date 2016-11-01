Euler3 = function(number = 600851475143){
  #Code directly taken from http://pastebin.com/M0qk0Kn8
  #First two variables are defined: One for the lowest prime factor and the number in question
 factor = 2
 num = number
 
 #Next, we will use a while loop to begin dividing our number in question by its factors
 #The loop finds the smallest factor of the number and divides the number by that factor, resulting in a new number
 #Eventually, the number is its only factor and that is the greatest prime factor of the original number
 while(num > factor){
   if (num %% factor == 0){
     num = num/factor
     factor = 2
   } else{
      factor = factor + 1
   }
 }
 return(factor)
 #Answer checks out in Project Euler!
}
