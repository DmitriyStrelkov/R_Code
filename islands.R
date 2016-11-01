#Metropolis algorithm
# 1. Pick left or right with fair coin
# 2. If new pop is bigger then move to next island
# 3. If new pop is smaller then move with probability of P_next/P_current

#Create a function to choose our direction
flip = function(){
  sample( x=c(-1,1), prob=c(.5,.5), size=1, replace=TRUE)
}


island_sampler = function(island_amount, steps){
  #Make islands
  islands = c(0:(island_amount + 1))
  P_islands = c(0,islands - 1 ,0)
  
  
  #Choose a starting location
  current_island = sample(1:island_amount, 1)
  island_frequency = c()
  
  for(i in 1:steps){
    move = flip()
    if(P_islands[move + current_island] > P_islands[current_island]){
      current_island = move + current_island
    } else if(P_islands[move + current_island] < P_islands[current_island]){
      P_move = P_islands[move + current_island]/P_islands[current_island]
      move_roll = sample( x=c(-1,0), 
                          prob=c(P_move, 1- P_move),
                          size = 1,
                          replace = TRUE)
      current_island = move_roll + current_island
    }
    island_frequency = c(island_frequency, current_island)
  }
  return(hist(island_frequency))
}