# Written in R by Dmitriy Strelkov.
# Problem 2: Write a program to distinguish between additional models of data.

# I chose to write a program that can distinguish between a normal distribution and a student-t distribution.
# They look very similar to the naked eye, and so it can be difficult to tell what type of distribution
# is really governing our data.

# My approach to this problem will be essentially identical to my approach to the previous problem.

# 1. Make a function returns the probability of a normal distribution.
# 2. Make a function returns the probability of a student-t distribution.
# 3. Make a function that returns the parameters associated with the normal distribution.
# 4. Make a function that returns the parameters associated with the student-t distribution.
# 5. Make a function that returns the probability distribution of the parameters of the normal distribution.
# 6. Make a function that returns the probability distribution of the parameters of the student-t distribution.
# 7. Make a function that determines which model is more likely.
# 8. Make a function that displays the parameter probability distribution for the likelier model.

# 1. Make a function returns the probability of a normal distribution.
normal_answer <- function(sigmalist, data){
  
  # This function returns the scaled probability and exponential scaling factor for the normally distributed
  # model.
  mean <- mean(data)
  loglikelihood <- numeric()
  
  # Compute the log likelihood for each value of sigma.
  i <- 1
  for(sigma in sigmalist){
    loglikelihood[i] <- sum(log((1/sqrt(2*sigma^2*pi))*exp(-(data - mean)^2/(2*sigma^2))))
    i <- i + 1
  }
  # This is the index of the most likely value of sigma.
  index <- which.max(loglikelihood)
  
  # The f and g subfunctions again. F returns the log likelihood of the input, g returns the 
  # exponentially scaled likelihood.
  f <- function(sigma){sum(log((1/sqrt(2*sigma^2*pi))*exp(-(data - mean)^2/(2*sigma^2))))}
  g <- function(sigma){exp(f(sigma) - f(sigmalist[index]))}
  
  # Generate the list of scaled likelihoods.
  likelihood_list <- numeric()
  for(i in 1:length(sigmalist)){
    likelihood_list[i] <- g(sigmalist[i])
  }
  
  # Find the scaled probability.
  scaled_probability <- (sum(likelihood_list) * (sigmalist[2] - sigmalist[1]))
  
  # Return the scaled probability and the exponential scaling factor.
  return(c(scaled_probability, f(sigmalist[index])))
}

# 2. Make a function returns the probability of a student-t distribution.

student_t_answer <- function(nulist, data){
  
  # Same as above, compute the log likelihood to find the most likely value of nu(degrees of freedom).
  logLikelihood <- numeric()
  i <- 1
  
  for(nu in nulist){
    logLikelihood[i] <- sum(log((gamma((nu + 1)/2)/(sqrt(nu*pi)*gamma(nu/2))) * 
                                  (1 + data^2/nu)^(-(nu + 1)/2)))
    i <- i + 1
  }
  # Find the max likelihood which gives the most likely value of nu.
  index <- which.max(logLikelihood)
  nu_true <- nulist[index]
  
  # Same subfunctions as before.
  f <- function(nu){sum(log((gamma((nu + 1)/2)/(sqrt(nu*pi)*gamma(nu/2))) * 
                              (1 + data^2/nu)^(-(nu + 1)/2)))}
  g <- function(nu){exp(f(nu) - f(nu_true))}
  
  # Generate scaled likelihoods.
  likelihood_list <- numeric()
  for(i in 1:length(nulist)){
    likelihood_list[i] <- g(nulist[i])
  }
  
  # Generate scaled probability.
  scaled_probability <- (sum(likelihood_list) * (nulist[2] - nulist[1]))
  
  # Return scaled probability and exponential scaling factor.
  return(c(scaled_probability, f(nulist[index])))
}

# 3. Make a function that returns the parameters associated with the normal distribution.
# 4. Make a function that returns the parameters associated with the student-t distribution.

# In the interest of saving time, I embedded these features into the following functions.

# 5. Make a function that returns the probability distribution of the parameters of the normal distribution.

normal_p_distro <- function(sigmalist, data){
  
  # Find the most likely sigma.
  mean <- mean(data)
  loglikelihood <- numeric()
  
  i <- 1
  for(sigma in sigmalist){
    loglikelihood[i] <- sum(log((1/sqrt(2*sigma^2*pi))*exp(-(data - mean)^2/(2*sigma^2))))
    i <- i + 1
  }
  index <- which.max(loglikelihood)
  
  # Create a small region around the most likely sigma as the region of interest.
  sigma_true <- sigmalist[index]
  sigma_low <- .95 * sigma_true
  sigma_high <- 1.053 * sigma_true
  
  graphed_sigma <- seq(sigma_low, sigma_high, by =  sigma_true/1000)
  
  # Same subfunctions as before.
  f <- function(sigma){sum(log((1/sqrt(2*sigma^2*pi))*exp(-(data - mean)^2/(2*sigma^2))))}
  g <- function(sigma){exp(f(sigma) - f(sigma_true))}
  
  loglikelihood_list <- numeric()
  
  # Generate a list of log likelihoods.
  for(i in 1:length(graphed_sigma)){
    loglikelihood_list[i] <- f(graphed_sigma[i])
  }
  
  # Create a normalized probability distribution of the most likely values of sigma.
  normalized_probability <- 1/sum(exp(loglikelihood_list - max(loglikelihood_list))) * 
    exp(loglikelihood_list - max(loglikelihood_list))
  
  # Plot the resultant distribution.
  plot(graphed_sigma, normalized_probability, xlab = 'Sigma', ylab = 'Probability')
  
  # Print a 95% confidence interval.
  i <- 1
  j <- length(graphed_sigma)
  
  while(sum(normalized_probability[i:j]) > .96){
    i <- i + 1
    j <- j - 1
  }
  
  print(paste(c("The probability of sigma being between", graphed_sigma[i], "and", graphed_sigma[j], "is", 
                round(sum(normalized_probability[i:j]), 3)), collapse = ' '))
  
  # Print the most likely value of sigma.
  print(paste(c("The most likely value of sigma is", graphed_sigma[which.max(normalized_probability)]), 
              collapse = ' '))
  
  # Return the expected sigma.
  return(sigmalist[index])
}

# 6. Make a function that returns the probability distribution of the parameters of the student-t distribution.

student_t_p_distro <- function(nulist, data){
  
  # First find the most likely value of nu to identify a region of interest.
  # It is important to note that the values of nu must be whole numbers because they represent degrees of freedom.
  # They also cannot be negative.
  logLikelihood <- numeric()
  i <- 1
  
  for(nu in nulist){
    logLikelihood[i] <- sum(log((gamma((nu + 1)/2)/(sqrt(nu*pi)*gamma(nu/2))) * 
                                  (1 + data^2/nu)^(-(nu + 1)/2)))
    i <- i + 1
  }
  index <- which.max(logLikelihood)
  
  # Create region of interest.
  nu_true <- nulist[index]
  nu_low <-  nu_true - 1
  nu_high <- nu_true + 1
  
  graphed_nu <- seq(nu_low, nu_high, by =  nu_true/1000)
  
  # Same subfunctions as before: f for generating the log likelihood and g for generating the scaled probability.
  f <- function(nu){sum(log((gamma((nu + 1)/2)/(sqrt(nu*pi)*gamma(nu/2))) * 
                              (1 + data^2/nu)^(-(nu + 1)/2)))}
  g <- function(nu){exp(f(nu) - f(nu_true))}
  
  loglikelihood_list <- numeric()
  
  # Generate the list of log likelihoods.
  for(i in 1:length(graphed_nu)){
    loglikelihood_list[i] <- f(graphed_nu[i])
  }
  
  # Generate the normalized probability.
  normalized_probability <- 1/sum(exp(loglikelihood_list - max(loglikelihood_list))) * 
    exp(loglikelihood_list - max(loglikelihood_list))
  
  # Plot the distribution.
  plot(graphed_nu, normalized_probability, xlab = 'Nu', ylab = 'Probability')
  
  # Print the 95% confidence interval.
  i <- 1
  j <- length(graphed_nu)
  
  while(sum(normalized_probability[i:j]) > .96){
    i <- i + 1
    j <- j - 1
  }
  
  print(paste(c("The probability of nu being between", graphed_nu[i], "and", graphed_nu[j], "is", 
                round(sum(normalized_probability[i:j]), 3)), collapse = ' '))
  
  # Print the most likely value of nu.
  print(paste(c("The most likely value of nu is", graphed_nu[which.max(normalized_probability)]), 
              collapse = ' '))
  
  # Return the expected amount of degrees of freedom.
  return(nulist[index])
}

# 7. Make a function that determines which model is more likely.
norm_or_t_tester <- function(sigmalist, nulist, data){
  
  # This is literally the same program as in Problem 1 except with different variable names to reflect
  # the different models used. I will refrain from reiterating the comments on this program.
  
  # Normal will be Model 1, student-t will be Model 2.
  normal_model_results <- normal_answer(sigmalist, data)
  t_model_results <- student_t_answer(nulist, data)
  
  scaling_factor_powers <- c(normal_model_results[2], t_model_results[2])
  scaled_model_probabilities <- c(normal_model_results[1], t_model_results[1])
  
  BF_log <- (scaling_factor_powers[1] - scaling_factor_powers[2]) + 
    log(scaled_model_probabilities[1]/scaled_model_probabilities[2])
  
  if((scaling_factor_powers[1] - scaling_factor_powers[2]) + 
     log(scaled_model_probabilities[1]/scaled_model_probabilities[2]) > 0){
    
    print("The distribution is normal. The probability and relative likelihood is displayed below.")
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
    
    print("The distribution is student-t. The probability and relative likelihood is displayed below.")
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
  #This function can identify whether or not a distribution is normal/student_t very well.
}

# 8. Make a function that displays the parameter probability distribution for the likelier model.

Resultant_distribution_finder_2 <- function(sigmalist, nulist, data){
  
  # This is also a copy of the function I used in problem 1.
  # Once again, model 1 is normal, model 2 is student_t.
  
  m <- norm_or_t_tester(sigmalist, nulist, data)
  
  if(m == 1){
    normal_p_distro(sigmalist, data)
  } else if(m == 2){
    student_t_p_distro(nulist, data)
  }
}

# Let's test it on some data!

norm_data <- rnorm(10000, 0, 1.2) # sigma = 1.2
t_data <- rt(10000, 3) # nu = 3

sigmalist <- seq(.7, 1.4, by = .01)
nulist <- seq(1, 10, by = 1)

# Drum roll please...
par(mfrow = c(2,1))
Resultant_distribution_finder_2(sigmalist, nulist, norm_data)
Resultant_distribution_finder_2(sigmalist, nulist, t_data)

# Looks like we return accurate results!