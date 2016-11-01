#Written in R by Dmitriy Strelkov
#Collection of different functions necessary to flip a coin, produce a prior distribution of coin bias,
#produce a list of bias values, and to produce a normalized posterior plot.


#This function flips a coin and returns a vector of T/F(heads/tails) values as its result.
flipCoin <- function(N,pHeads){
  
  flipSequence <- sample( x=c(0,1), prob=c(1-pHeads,pHeads), size=N, replace=TRUE) == 1
  
  return(flipSequence)
}

#This function takes our distribution function(in this case, the function was obtained
#from Kruschke's "Doing Bayesian Data Analysis" Chapter 5) and calculates a probability for each
#possible bias.
getLikelihood <- function(data, thetalist){
  
  likelihood <- c()
  headsCount <- sum(data)
  tailsCount <- length(data) - sum(data)
  
  likelihood <- (thetalist ^ headsCount) * ((1-thetalist) ^ (tailsCount))
  return(likelihood)
}

#This function gets the posterior probability distribution, normalizes it, and plots it.
getPosterior <- function(prior, likelihood){
  posterior <- (prior * likelihood)
  posterior <- posterior/sum(posterior) #normalization
  return(plot(thetalist, posterior, type = "l"))
}

#Generates a prior probability distribution that is not biased.
get_fair_prior <- function(thetalist){
  fprior <- getLikelihood(c(0,1), thetalist)
  return(fprior)
}

#Generates a prior probability distribution that favors heads.
get_heads_prior <- function(thetalist){
  hprior <- getLikelihood(c(0,1,1), thetalist)
  return(hprior)
}

#Generates a prior probability distribution that favors tails.
get_tails_prior <- function(thetalist){
  tprior <- getLikelihood(c(0,0,1), thetalist)
  return(tprior)
}

#Creates a list of biases to be tested.
get_thetalist = function(){
  step <- readline("How many values of theta do you want to test?")
  step <- as.numeric(step)
  list_of_thetas <- seq(0, 1, by = 1/step)
  return(list_of_thetas)
}

#This function generates our posterior given a prior
posterior_generator = function(){
  
  #get the type of prior you want
  prior_type <- readline("What bias exists in your prior? f for fair, h for heads, t for tails.")
  
  if(prior_type == "f"){
    prior <- get_fair_prior(thetalist)
  } else if(prior_type == "h"){
    prior <- get_heads_prior(thetalist)
  } else{
    prior <- get_tails_prior(thetalist)
  }
  #ask for data
  need_data <- readline("Do you have data already? y/n") == "n"
  
  if(need_data){
    #flip a coin with some bias and some amount of flips
    flips <- as.numeric(readline("How many flips?"))
    bias <- as.numeric(readline("What bias?"))
    data <- flipCoin(flips, bias)
  } else{
    data <- readline("Please enter your data now.")
    data <- as.numeric(strsplit(as.character(data), "")[[1]]) == 1
  }
  
  likelihood <- getLikelihood(data, thetalist)
  getPosterior(prior, likelihood)
}

#Now lets run our program!
thetalist <- get_thetalist()
posterior_generator()

#data must be entered as one long string of 0's for tails and 1's for heads.