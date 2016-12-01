#Written in R by Dmitriy Strelkov
#Problem 1: Determine whether a given data set is exponential or pareto.

#The data set in question was provided by Eric Corwin, called ExponentialOrPowerLaw.dat

#To read the data set into R, enter the following commands:

# data <- read.delim("ExponentialOrPowerLaw.dat", sep = "\n")
# data <- t(data)

#From now on I will be using 'data' to refer to this data.

#In order to fully determine whether or not the data is exponential or pareto, I will do
#N things:

# 1. Make a function returns the probability of an exponential distribution.
# 2. Make a function returns the probability of a pareto distribution.
# 3. Make a function that returns the parameters associated with the exponential distribution.
# 4. Make a function that returns the parameters associated with the pareto distribution.
# 5. Make a function that returns the probability distribution of the parameters of the exponential distribution.
# 6. Make a function that returns the probability distribution of the parameters of the pareto distribution.
# 7. Make a function that determines which model is more likely.
# 8. Make a function that displays the parameter probability distribution for the likelier model.

# 1. Make a function returns the probability of an exponential distribution.

expAnswer = function(taulist, data){
  
  # Create a list of log likelihoods to determine the most likely tau(time constant) for the
  # exponential distribution.
  logLikelihood = length(data) * log(1/taulist) + (-1/taulist)*sum(data)
  
  # Empty list for our likelihood.
  likelihood_list <- c()
  
  # The index at which the log likelihood is max is the most likely value of tau.
  tau_true <- taulist[which.max(exp(logLikelihood - max(logLikelihood)))]

  # This first sub-function, f, returns the log likelihood of a given tau value.
  f <- function(tau){length(data) * log(1/tau) + (-1/tau) * sum(data)}
  
  # This second sub-function, g, returns the likelihood of any tau value, scaled such that the
  # most likely tau value is represented as 1.
  g <- function(tau){exp(f(tau) - f(tau_true))}
  
  # Generate a list of our likelihoods.
  for(i in 1:length(taulist)){
    likelihood_list[i] <- g(taulist[i])
  }
  
  # When normalized, this value is the probability that this model represents the data tested.
  # Since we scaled the likelihoods by a factor of the exponential of f(tau_true), we will
  # keep track of this scaling factor and use it later when comparing to the pareto distribution.
  scaled_probability <- (sum(likelihood_list) * (taulist[2] - taulist[1]))
  
  #We return the scaled probability and the associated exponential scaling factor.
  return(c(scaled_probability, f(tau_true)))
}

# 2. Make a function returns the probability of an pareto distribution.

paretoAnswer = function(alist, data){
  
  # The approach is essentially identical to above. Here we create an initial log likelihood.
  # For the b value in the pareto distribution, we simply take the minimum of the data and call it
  # our cutoff point.
  b <- min(data)
  oldlogLikelihood <- sum(-(alist[1] + 1)*log(data) + alist[1]*log(b) + log(alist[1]))
  
  
  # Empty variable for to hold the final a value.
  a_f <- numeric()
  
  # We iterate through all values of a, and if the new log likelihood is higher than the old
  # log likelihood, we replace the log likelihood and its associated a value.
  for(i in 1:length(alist)){
    a <- alist[i]
    newlogLikelihood <- sum(-(a + 1)*log(data) + a*log(b) + log(a))
    
    if(newlogLikelihood > oldlogLikelihood){
      oldlogLikelihood <- newlogLikelihood
      a_f <- a
    }
  } 
  
  
  # Same f subfunction as before, returning the log likelihood of any a value.
  f <- function(a){sum(log(a) + a*log(b) - (a + 1)*log(data))}
  
  # Same g subfunction as before, returning the likelihood scaled by an exponential factor of 
  # f(a_f)
  g <- function(a){exp(f(a) - f(a_f))}
  
  likelihood_list <- numeric()
  
  # Generate a list of scaled likelihoods from all a values.
  for(i in 1:length(alist)){
    likelihood_list[i] <- g(alist[i])
  }
  
  # As above, we determine the scaled probability that the data in question supports a pareto distribution.
  scaled_probability <- (sum(likelihood_list) * (alist[2] - alist[1]))
  
  # Finally we return the scaled probability and the exponential scaling factor.
  return(c(scaled_probability, f(a_f)))
}

# 3. Make a function that returns the parameters associated with the exponential distribution.

# This is the first part of the first function, except this one stops when it finds the true tau
# value.
tau_finder <- function(taulist, data){
  logLikelihood <- length(data) * log(1/taulist) + (-1/taulist)*sum(data)
  
  tau_true <- taulist[which.max(exp(logLikelihood - max(logLikelihood)))]
  
  return(tau_true)
}

# 4. Make a function that returns the parameters associated with the pareto distribution.

# This is the pareto analog to the above function.
a_finder <- function(alist, data){
  oldlogLikelihood <- sum(-(alist[1] + 1)*log(data) + alist[1]*log(min(data)) + log(alist[1]))
  
  a_f <- numeric()
  
  for(i in 1:length(alist)){
    a <- alist[i]
    newlogLikelihood <- sum(-(a + 1)*log(data) + a*log(min(data)) + log(a))
    
    if(newlogLikelihood > oldlogLikelihood){
      oldlogLikelihood <- newlogLikelihood
      a_f <- a
    }
  } 
  return(a_f)
}

# 5. Make a function that returns the probability distribution of the parameters of the exponential distribution.

exponential_p_distro <- function(taulist, data){
  
  # Begin by shortening the list of graphed taus to only contain the region of interest.
  
  tau_true <- tau_finder(taulist, data)
  tau_low <- .95 * tau_true
  tau_high <- 1.053 * tau_true
  
  graphed_tau <- seq(tau_low, tau_high, by =  tau_true/1000)
  
  # Same subfunctions as before.
  f <- function(tau){length(data) * log(1/tau) + (-1/tau) * sum(data)}
  g <- function(tau){exp(f(tau) - f(tau_true))}
  
  loglikelihood_list <- numeric()
  
  # We create a list of log likelihoods for each tau of interest.
  for(i in 1:length(graphed_tau)){
    loglikelihood_list[i] <- f(graphed_tau[i])
  }
  
  # Now we generate a normalized probability of any graphed tau value.
  normalized_probability <- 1/sum(exp(loglikelihood_list - max(loglikelihood_list))) * 
    exp(loglikelihood_list - max(loglikelihood_list))
  
  # Plot of the taus of interest and their probabilities.
  plot(graphed_tau, normalized_probability, xlab = 'Tau', ylab = 'Probability')
  
  # Output an 95% confidence interval
  i <- 1
  j <- length(graphed_tau)
  
  while(sum(normalized_probability[i:j]) > .96){
    i <- i + 1
    j <- j - 1
  }
  
  print(paste(c("The probability of tau being between", graphed_tau[i], "and", graphed_tau[j], "is", 
                round(sum(normalized_probability[i:j]), 3)), collapse = ' '))
  
  # Output the most likely value of tau.
  print(paste(c("The most likely value of tau is", graphed_tau[which.max(normalized_probability)]), 
               collapse = ' '))
  
  # Return the expected tau constant.
  return(tau_true)
}

# 6. Make a function that returns the probability distribution of the parameters of the pareto distribution.

pareto_p_distro <- function(alist, data){
  
  # This function is the pareto analog to the one above it.
  
  # Begin by selecting a values in the region of interest as before. Recall that a must be greater than 0.
  a_f <- a_finder(alist, data)
  a_l <- .95 * a_f
  a_h <- 1.053 * a_f
  b <- min(data)
  
  graphed_a <- seq(a_l, a_h, by = a_f/1000)
  
  # Same subfunctions as before.
  f <- function(a){sum(log(a) + a*log(b) - (a + 1)*log(data))}
  g <- function(a){exp(f(a) - f(a_f))}
  
  loglikelihood_list <- numeric()
  
  # Create a list of the log likelihoods.
  for(i in 1:length(graphed_a)){
    loglikelihood_list[i] <- f(graphed_a[i])
  }
  
  # Generate a normalized probability distribution for the values of a.
  normalized_probability <- 1/sum(exp(loglikelihood_list - max(loglikelihood_list))) * 
    exp(loglikelihood_list - max(loglikelihood_list))
  
  # Plot the distribution.
  plot(graphed_a, normalized_probability, xlab = 'a', ylab = 'Probability')
  
  # Print out the 95% confidence interval.
  i <- 1
  j <- length(graphed_a)
  
  while(sum(normalized_probability[i:j]) > .96){
    i <- i + 1
    j <- j - 1
  }
  
  print(paste(c("The probability of a being between", graphed_a[i], "and", graphed_a[j], "is", 
                round(sum(normalized_probability[i:j]), 3)), collapse = ' '))
  
  # Print the most likely value of a.
  print(paste(c("The most likely value of a is", graphed_a[which.max(normalized_probability)]), 
               collapse = ' '))
  
  # Return the most likely value of a.
  return(a_f)
}

# 7. Make a function that determines which model is more likely.

exp_or_pareto_tester <- function(taulist, alist, data){
  
  # Model 1 will be exponential, Model 2 will be pareto. Get the scaled probabilities and exponential
  # scaling factors for each model.
  exponential_model_results <- expAnswer(taulist, data)
  pareto_model_results <- paretoAnswer(alist, data)
  
  # Since the scaling factors are exponential scaling factors, the total Bayes factor is equal to
  # the ratio of the exponentials of the scaling factors multiplied by the ratio of the 
  # respective scaled model probabilities.
  
  # We use the logarithm of the Bayes factor to determine which model is more likely. If the BF
  # is greater than 0, then model 1 is more likely. If it is less than 0, model 2 is more likely.
  scaling_factor_powers <- c(exponential_model_results[2], pareto_model_results[2])
  scaled_model_probabilities <- c(exponential_model_results[1], pareto_model_results[1])
  
  BF_log <- (scaling_factor_powers[1] - scaling_factor_powers[2]) + 
    log(scaled_model_probabilities[1]/scaled_model_probabilities[2])
  
  # Return which distribution is controlling the data and how much more likely it is than the other model.
  # Then, return the model number.
  if((scaling_factor_powers[1] - scaling_factor_powers[2]) + 
     log(scaled_model_probabilities[1]/scaled_model_probabilities[2]) > 0){
    print("The distribution is exponential. The probability and relative likelihood is displayed below.")
    print("Model one is more likely than model two by the exponential of:")
    print(BF_log)
    
    percentage <- (exp(scaling_factor_powers[1]) * scaled_model_probabilities[1])/
      sum(exp(scaling_factor_powers) * scaled_model_probabilities)
    
    print("The probability of model one is:")
    if(is.nan(percentage)){
      print(1)
      
      #Return the model #
      return(1)
    } else{
      print(percentage)
      return(1)
    }
  }else if((scaling_factor_powers[1] - scaling_factor_powers[2]) + 
           log(scaled_model_probabilities[1]/scaled_model_probabilities[2]) < 0){
    print("The distribution is pareto. The probability and relative likelihood is displayed below.")
    print("Model two is more likely than model one by the exponential of:")
    print(-BF_log)
    
    percentage <- (exp(scaling_factor_powers[1]) * scaled_model_probabilities[1])/
      sum(exp(scaling_factor_powers) * scaled_model_probabilities)
    
    print("The probability of model two is:")
    if(is.nan(percentage)){
      print(1)
      return(2)
    } else{
      print(percentage)
      return(2)
    }
  }
  #This function can identify whether or not a distribution is exponential/pareto very well.
}

# 8. Make a function that displays the parameter probability distribution for the likelier model.


Resultant_distribution_finder <- function(taulist, alist, data){
  # Once again, model 1 is exponential, model 2 is pareto.
  # We find the model responsible for the data, then call the appropriate probability distribution generator.
  m <- exp_or_pareto_tester(taulist, alist, data)
  
  if(m == 1){
    exponential_p_distro(taulist, data)
  } else if(m == 2){
    pareto_p_distro(alist, data)
  }
}

# Now we shall test!
library('extraDistr')
exp_data <- rexp(10000, 1/19) # tau = 19
pareto_data <- rpareto(10000, 1.2, 5) # a = 1.2, b = min(data) = 5

real_data <- read.delim("ExponentialOrPowerLaw.dat", sep = "\n")
real_data <- t(real_data)
real_data <- c(real_data)

taulist <- seq(17, 23, by = .1)
alist <- seq(.7, 1.4, by = .01)

# Proving it works by trying data distributed according to known values.
par(mfrow = c(3,1))
Resultant_distribution_finder(taulist, alist, exp_data)
Resultant_distribution_finder(taulist, alist, pareto_data)

# Now for the moment of truth!
Resultant_distribution_finder(taulist, alist, real_data)

# Looks like an exponential distribution with a tau of 20.04, distributed similarly to
# rexp(10000, 1/20.04)