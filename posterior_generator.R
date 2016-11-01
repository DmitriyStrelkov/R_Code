source("flip_likelihood_posterior.R")

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
