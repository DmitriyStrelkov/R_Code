#Metropolis algorithm
# 1. Pick left or right with fair coin
# 2. If new pop is bigger then move to next island
# 3. If new pop is smaller then move with probability of P_next/P_current
# The output is a histogram of frequency of island visitation, where the visitation
# frequency is directly proportional to the island population.
# Each island has a population equal to its index(island 1 has 1 person, island 2 has 2 people, etc.)

#Create a function to choose our direction
flip = function(){
  sample( x=c(-1,1), prob=c(.5,.5), size=1, replace=TRUE)
}


island_sampler = function(island_amount, steps){
  #Make islands
  #Need at least 4 islands. Why? Also, 4 islands has a lot of bugs. In general,
  #if there is an error in the output, try the same input again. Error usually resolves itself.
  #Also, 1000+ steps should be used for a representative sample.
  
  #Create islands and the populations of those islands.
  islands = c(0:(island_amount + 1))
  P_islands = c(0,islands - 1 ,0)
  
  
  #Choose a starting location
  current_island = sample(1:island_amount, 1)
  island_frequency = c()

  for(i in 1:steps){
    move = flip()
    #This next part of code is buggy and returns error "argument of length zero" occasionally.
    
    #If the population of the next island is higher than the population of the current island, move.
    if(P_islands[move + current_island] > P_islands[current_island]){
      current_island = move + current_island
    } 
    #If the population of the next island is lower than the population of the current island, move
    #with a probability equal to the ratio of the population of the new island to the current island.
    else if(P_islands[move + current_island] < P_islands[current_island]){
      P_move = P_islands[move + current_island]/P_islands[current_island]
      move_roll = sample( x=c(-1,0), 
                          prob=c(P_move, 1- P_move),
                          size = 1,
                          replace = FALSE)
      current_island = move_roll + current_island
    } 
    #Create a memory of each island visit.
    island_frequency = c(island_frequency, current_island)
  }
  #For some reason, subtract 3 to get the histogram to display the proper island number. This has something
  #to do with the fact that the minimum islands allowed are 4. I don't know why this happens.
  island_frequency = island_frequency - 3
  return(hist(island_frequency))
}

#Another common error is the "NA in probability vector" error, I don't know what causes this. Code runs 
#well most of the time.

#Conclusion is that this a working Metropolis algorithm, but there are errors that occur as a result
#of using dummy islands with 0 population to create a zero percent chance of moving to those islands.