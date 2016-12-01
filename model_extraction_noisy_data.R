# Written in R by Dmitriy Strelkov
# Problem 3: Write a program to estimate the model parameters of an exponential function that is
# affected by Gaussian noise, and estimate the model parameters of the noise as well.

# Library required: geophys

#install.packages('geophys')
library('geophys')

exponential_with_normal_noise <- function(lambdalist, sigmalist, data){
  
  # These three variables will help us construct our set of log likelihoods. Current_like is used
  # as an intermediary between the data and the log likelihood because, due to the error function,
  # some of our likelihood values will be zero, and we need to get rid of those values before we take
  # the log.
  current_like <- c()
  sigma_loglikelihood <- c()
  lambda_loglikelihood <- c()
  
  # Convenient indices.
  i <- 1
  j <- 1
  
  # We will approach the determination of sigma and lambda through the constraining of variables. First
  # we constrain lambda to a reasonable value, 1 over the mean. Then, we find the most likely value of sigma.
  # Then we constrain sigma to that most likely value and repeat the same process for finding lambda.
  # This cuts down on computation time and also produces accurate results.
  
  lambda_init <- 1/mean(data)
  
  for(sigma in sigmalist){
    
    # Generate the likelihood of all data points for a given sigma.
    current_like <- (lambda_init/2)*exp((sigma^2 * lambda_init^2)/(2)) * exp(-data*lambda_init)*
          (1 + erf((data - lambda_init*sigma^2)/(sigma*sqrt(2))))
    
    # Keep all numbers greater than zero.
    current_like <- current_like[current_like > 0]
    
    # Create a list of the total log likelihoods for all sigma.
    sigma_loglikelihood[i] <- sum(log(current_like))
    
    i <- i + 1
  }
  
  # Determine the maximum log likelihood and the corresponding sigma value.
  sigma_index <- which(sigma_loglikelihood == max(sigma_loglikelihood))
  sigma_true <- sigmalist[sigma_index]
  
  # Repeat the process for all lambda, fixing sigma as sigma_true.
  for(lambda in lambdalist){
    current_like <- (lambda/2)*exp((sigma_true^2 * lambda^2)/(2)) * exp(-data*lambda)*
      (1 + erf((data - lambda*sigma_true^2)/(sigma_true*sqrt(2))))
    current_like <- current_like[current_like > 0]
      
    lambda_loglikelihood[j] <- sum(log(current_like))
      
    j <- j + 1
  }
  
  # Determine the maximum log likelihood and the corresponding lambda (and tau) value.
  lambda_index <- which(lambda_loglikelihood == max(lambda_loglikelihood))
  
  lambda_true <- lambdalist[lambda_index]
  tau_true <- 1/lambda_true
  
  # We will construct a matrix of likelihoods so we can visualize the 2D probability distribution
  # of lambda and sigma.
  lambda_matrix <- matrix(exp(lambda_loglikelihood - max(lambda_loglikelihood)),
                          nrow = length(lambda_loglikelihood), ncol = 1)
  sigma_matrix <- matrix(exp(sigma_loglikelihood - max(sigma_loglikelihood)),
                         nrow = 1, ncol = length(sigma_loglikelihood))
  
  loglikelihood <- lambda_matrix %*% sigma_matrix

  # Print the most probable values of tau, lambda, sigma.
  print(paste(c('The most probable tau value for the underlying exponential distribution is:', tau_true)))
  print(paste(c('The most probable lambda(1/tau) value for the underlying exponential distribution is:', 
                lambda_true)))
  print(paste(c('The most probable sigma value for the noise is:', sigma_true)))
  
  # Return the 2D probability distribution.
  return(persp(lambdalist, sigmalist, loglikelihood, phi = 37, theta = 60))
}

# Let's test it out on some data!

exp_data_normal_noise_data <- rexp(10000, 1/21) + rnorm(10000, 0, .8)

lambdalist <- 1/(rev(seq(18, 24, by = .1)))
sigmalist <- seq(.3, 1.5, by = .01)

# IT'S HAPPENING!!!
exponential_with_normal_noise(lambdalist, sigmalist, exp_data_normal_noise_data)

# This makes me very happy :)