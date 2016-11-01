#Code originally taken from section 4.5 of "Doing Bayesian Analysis 2nd Edition" by John Kruschke
#Code edited by Dmitriy Strelkov in R.

flipCoin = function(N,pHeads){ # Specify the total number of flips, denoted N.
  #Create a function that flips the coin N times with a heads bias pHeads
  #N is a positive integer, pHeads must be between 0 and 1

# Generate a random sample of N flips (heads=1, tails=0):

flipSequence = sample( x=c(0,1), prob=c(1-pHeads,pHeads), size=N, replace=TRUE)
#Samples from the vector 0 or 1(tails or heads), with probabilities determined by pHeads,
#creating a sample of size N. Replace is TRUE to refresh the vector on every sample.

#Compute the running proportion of heads:

r = cumsum( flipSequence ) # Cumulative sum: Number of heads at each step.
#How many times heads was flipped

n = 1:N # Number of flips at each step.
#How many times the coin was flipped

runProp = r / n # Component by component division.
#The simulated probability of flipping heads

# Graph the running proportion:

plot( n , runProp , type="o" , log="x" , col="skyblue" ,
      
      xlim=c(1,N) , ylim=c(0.0,1.0) , cex.axis=1.5 ,
      
      xlab="Flip Number" , ylab="Proportion Heads" , cex.lab=1.5 ,
      
      main="Running Proportion of Heads" , cex.main=1.5 )

#Plots the proportion of heads as the coin is flipped.

# Plot a dotted horizontal reference line:

abline( h=pHeads , lty="dotted" )

# Display the beginning of the flip sequence:

flipLetters = paste( c("T","H")[flipSequence[1:10]+1] , collapse="" )

displayString = paste0( "Flip Sequence=" , flipLetters , "." )

text( N , .9 , displayString , adj=c(1,0.5) , cex=1.3 )

# Display the relative frequency at the end of the sequence.

text( N , .8 , paste("End Proportion =",runProp[N]) , adj=c(1,0.5) , cex=1.3 )
}

#Task: Change the probability of heads to .8 and correct the reference line to reflect that change.
#The reference line moves automatically with the probability.
#All we have to do is call our newly made function with a desired number of flips and the desired probability.

flipCoin(500,.8)

#Done!