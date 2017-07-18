#Written in R by Dmitriy Strelkov

setwd('C:\\Users\\Dmitriy\\DoE\\Final\\Data')
data <- read.csv('Tea.csv')

install.packages('plyr')
library('plyr')

### IMPORTANT PLEASE READ ###
# This code is designed with pauses in it and is meant to be looked at in tandem with the attached report.
# Proceed through the code slowly so that things do not get out of order; lest the story be spoiled.
# There are a vast amount of comments with varying degrees of detail in this code and they are meant to explain
# my thought process as I was writing the code. The results obtained here are essentially summarized and
# organized for clarity in the report, and as such the best way to fully understand what I did and why I did it
# is to consider the two documents together rather than separately.

### THANK YOU FOR READING ### 

# We are going to investigate the tea data from the other DoE class.

# The data in the file uses these 4 factors to classify scores into categories. Each column has a classifier
# that includes one of each of the following:

# G or B (green/black)
# 6 or 3 (6min/3min)
# S or N (sweetener/no sweetener)
# H or L (high leaf mass/low leaf mass)

#In looking at this data, we will try to determine the following things:
# Which of the following four different factors has the most significant effect on perceived quality:
#  - Green vs. Black tea
#  - 6 minute vs. 3 minute brew time
#  - Sweetened vs. Unsweetened tea
#  - 3 grams of leaves vs 1.5 grams of leaves in the brew


# How significant is the prior preference of tea in the determination of perceived quality? Is it more or
# less significant than brewing factors?

# How do the rating distributions vary between preferences for each type of tea?

# What are the probability distributions of tea scores for each taste preference?

total_results_aggregate <- c()
par(mfrow = c(1,1))

for(n in 2:16){
  total_results_aggregate <- c(total_results_aggregate, data[,n])
}

hist(total_results_aggregate, 100)
readline()

#Looks like a normal distribution for the sum of all data. Let's check using the following algorithm.

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
  sigma_low <- .5 * sigma_true
  sigma_high <- 1.5 * sigma_true
  
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

sigmalist <- seq(1.5, 2.5, by = .01)
normal_p_distro(sigmalist, total_results_aggregate) -> sigmatest
readline()

# Now we'll compare a normal distribution with the resultant parameters to the data's entirity.
par(mfrow = c(1,2))

rand_norm_like_data <- rnorm(180, mean(total_results_aggregate), sigmatest)

hist(rand_norm_like_data, 10)
hist(total_results_aggregate, 10)
readline()

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

print(normal_answer(sigmalist, rand_norm_like_data))
normal_p_distro(sigmalist, rand_norm_like_data)

print(normal_answer(sigmalist, total_results_aggregate))
normal_p_distro(sigmalist, total_results_aggregate)
readline()


#The actual data is about equally likely to be a normal distribution as the normally
#distributed random sample.


# Let's compare the differences between the data for green tea and black tea.

green_tea_data <- c()
for(n in 2:8){
  green_tea_data <- c(green_tea_data, data[,n])
}

black_tea_data <- c()
for(n in 9:16){
  black_tea_data <- c(black_tea_data, data[,n])
}

hist(green_tea_data, 100); text(8, 12, labels = paste('Mean = ', round(mean(green_tea_data), 4), sep = '\n'))
hist(black_tea_data, 100); text(3, 18, labels = paste('Mean = ', round(mean(black_tea_data), 4), sep = '\n'))
readline()

# Looks like on average, people like green tea better than they like black tea. In fact, there are no '1' ratings
# for green tea. Additionally, there are two '10' ratings for green tea and there aren't any for black!
# Let's take a look at the distribution of tea preferences and see if that has anything to do with it...

par(mfrow = c(1,1))
preferences <- data[,17]

# Since the preferences have levels, use barplot.

barplot(table(preferences))
readline()

# We will use the following distributions of scores as baselines for each of the different taste preferences:
#  - Sweet preference all tea
#  - Bitter preference all tea
#  - No preference all tea

# Looks like most people have no preference, and sweet preferences win out slightly over bitter. This trend
# is reflected in the earlier look at mean ratings. Let's look now at how preference changes these distributions.

# We will compare the baseline distributions the following distributions:
#  - Sweet preference green tea
#  - Sweet preference black tea
#  - Bitter preference green tea
#  - Bitter preference black tea
#  - No preference green tea
#  - No preference black tea


sweet_pref_all <- c()
for(n in 1:12){
  if(data$Preference[n] == "Sweet"){
    sweet_pref_all <- c(sweet_pref_all, as.numeric(data[n, 2:16]))
  }
}

bitter_pref_all <- c()
for(n in 1:12){
  if(data$Preference[n] == "Bitter"){
    bitter_pref_all <- c(bitter_pref_all, as.numeric(data[n, 2:16]))
  }
}

no_pref_all <- c()
for(n in 1:12){
  if(data$Preference[n] == "No Preference"){
    no_pref_all <- c(no_pref_all, as.numeric(data[n, 2:16]))
  }
}

par(mfrow = c(1,3))

hist(no_pref_all, 100); text(2, 20, labels = paste('Mean = ', round(mean(no_pref_all), 4), sep = '\n'))
hist(bitter_pref_all, 100) ; text(8, 6, labels = paste('Mean = ', round(mean(bitter_pref_all), 4), sep = '\n'))
hist(sweet_pref_all, 100); text(3, 7, labels = paste('Mean = ', round(mean(sweet_pref_all), 4), sep = '\n'))

readline()

# Looks like the people who have no preference seem to, on average, like tea better than either group with
# a taste preference. Of those who have a taste preference, lower scores are seen in those with sweet
# taste preferences.

no_pref_green <- c()
for(n in 1:12){
  if(data$Preference[n] == "No Preference"){
    no_pref_green <- c(no_pref_green, as.numeric(data[n, 2:8]))
  }
}

no_pref_black <- c()
for(n in 1:12){
  if(data$Preference[n] == "No Preference"){
    no_pref_black <- c(no_pref_black, as.numeric(data[n, 9:16]))
  }
}

par(mfrow = c(1,3))

hist(no_pref_all, 100); text(2, 20, labels = paste('Mean = ', round(mean(no_pref_all), 4), sep = '\n'))
hist(no_pref_green, 100); text(4, 12, labels = paste('Mean = ', round(mean(no_pref_green), 4), sep = '\n'))
hist(no_pref_black, 100) ; text(3, 12, labels = paste('Mean = ', round(mean(no_pref_black), 4), sep = '\n'))
readline()

# Once again, we see that green tea is favored over black tea. Not only is green tea the more popular among
# those with no preference, this group/tea combination has the highest average rating so far. Black tea, meanwhile,
# is more popular in those with no taste preference than it is with all observers in total, but its rating
# is still below the average total tea score.

bitter_pref_green <- c()
for(n in 1:12){
  if(data$Preference[n] == "Bitter"){
    bitter_pref_green <- c(bitter_pref_green, as.numeric(data[n, 2:8]))
  }
}

bitter_pref_black <- c()
for(n in 1:12){
  if(data$Preference[n] == "Bitter"){
    bitter_pref_black <- c(bitter_pref_black, as.numeric(data[n, 9:16]))
  }
}

hist(bitter_pref_all, 100) ; text(8, 6, labels = paste('Mean = ', round(mean(bitter_pref_all), 4), sep = '\n'))
hist(bitter_pref_green, 100); text(6, 4, labels = paste('Mean = ', round(mean(bitter_pref_green), 4), sep = '\n'))
hist(bitter_pref_black, 100); text(4, 3, labels = paste('Mean = ', round(mean(bitter_pref_black), 4), sep = '\n'))
readline()

# For the first time, we see distributions that don't look like the normal distribution. Interestingly
# enough, this group exhibits a stronger preference for green tea and a stronger dislike of black tea than
# the group with no preferences in taste. This would seem to suggest that green tea is more bitter than
# black tea. We will examine the distribution of the data later. For now we will continue and just get a visual
# on the distribution of preferences for those who prefer sweet tea.

sweet_pref_green <- c()
for(n in 1:12){
  if(data$Preference[n] == "Sweet"){
    sweet_pref_green <- c(sweet_pref_green, as.numeric(data[n, 2:8]))
  }
}

sweet_pref_black <- c()
for(n in 1:12){
  if(data$Preference[n] == "Sweet"){
    sweet_pref_black <- c(sweet_pref_black, as.numeric(data[n, 9:16]))
  }
}

hist(sweet_pref_all, 100); text(3, 7, labels = paste('Mean = ', round(mean(sweet_pref_all), 4), sep = '\n'))
hist(sweet_pref_green, 100); text(6, 4, labels = paste('Mean = ', round(mean(sweet_pref_green), 4), sep = '\n'))
hist(sweet_pref_black, 100); text(3, 5, labels = paste('Mean = ', round(mean(sweet_pref_black), 4), sep = '\n'))
readline()

# The people who prefer sweet tea have an overall disdain for all kinds of tea. In fact, this group has the lowest
# score for black tea and green tea. The distribution of green scores seems to be somewhat uniform, and the
# distribution of black scores seems to resemble a normal distribution.

# Finally, let's compare all of the histograms.

par(mfrow = c(3,3))

hist(no_pref_all, 100); text(2, 20, labels = paste('Mean = ', round(mean(no_pref_all), 4), sep = '\n'))
hist(no_pref_green, 100); text(4, 12, labels = paste('Mean = ', round(mean(no_pref_green), 4), sep = '\n'))
hist(no_pref_black, 100) ; text(3, 12, labels = paste('Mean = ', round(mean(no_pref_black), 4), sep = '\n'))
hist(bitter_pref_all, 100) ; text(8, 6, labels = paste('Mean = ', round(mean(bitter_pref_all), 4), sep = '\n'))
hist(bitter_pref_green, 100); text(6, 4, labels = paste('Mean = ', round(mean(bitter_pref_green), 4), sep = '\n'))
hist(bitter_pref_black, 100); text(4, 3, labels = paste('Mean = ', round(mean(bitter_pref_black), 4), sep = '\n'))
hist(sweet_pref_all, 100); text(3, 7, labels = paste('Mean = ', round(mean(sweet_pref_all), 4), sep = '\n'))
hist(sweet_pref_green, 100); text(6, 4, labels = paste('Mean = ', round(mean(sweet_pref_green), 4), sep = '\n'))
hist(sweet_pref_black, 100); text(3, 5, labels = paste('Mean = ', round(mean(sweet_pref_black), 4), sep = '\n'))
readline()

# We can see that black tea is less popular among all groups, and that those with sweet preferences 
# have less appreciation for tea than the other groups.

# I will perform the same analysis on the rest of the data in a similar fashion. Since experience and 
# taste preference have a nearly 1:1 correlation, I will omit experience data, as 'experience' with tea seems
# to be an inherently poor preference parameter. The one and only interesting thing about tea experience is that
# individuals with more tea experience have stronger preferences for bitter tea. We will talk about why this
# could be in the discussion.

# Cut out extraneous experience data from frame.
data <- data[, 1:17]

# The next factor most significant factor I hypothesize will be sweetener vs no sweetener.

# Make data of tea ratings that are sweetened and unsweetened.
sweetener_added_data <- c(data[,2], data[,5], data[,6], data[,7], data[,9], data[,10], data[,15], data[,16])
nosweetener_data <- c(data[,3], data[,4], data[,8], data[,11], data[,12], data[,13], data[,14])

# Make vectors to indicate column positions; easier to manipulate.
sweetener_vector <- c(2,5,6,7,9,10,15,16)
nosweetener_vector <-c(3,4,8,11,12,13,14)

no_pref_sweetener <- c()
for(n in 1:12){
  if(data$Preference[n] == "No Preference"){
    no_pref_sweetener <- c(no_pref_sweetener, as.numeric(data[n, sweetener_vector]))
  }
}

no_pref_nosweetener <- c()
for(n in 1:12){
  if(data$Preference[n] == "No Preference"){
    no_pref_nosweetener <- c(no_pref_nosweetener, as.numeric(data[n, nosweetener_vector]))
  }
}

par(mfrow = c(1,3))

hist(no_pref_all, 100); text(2, 20, labels = paste('Mean = ', round(mean(no_pref_all), 4), sep = '\n'))
hist(no_pref_sweetener, 100); text(4, 12, labels = paste('Mean = ', round(mean(no_pref_sweetener), 4), sep = '\n'))
hist(no_pref_nosweetener, 100) ; text(3, 12, labels = paste('Mean = ', round(mean(no_pref_nosweetener), 4), sep = '\n'))
readline()

# Looks like the data here is suggesting that the difference between sweetener and no sweetener is indiscernable
# between those who have no taste preference, which is what we would expect. Both distributions seem normal.

bitter_pref_sweetener <- c()
for(n in 1:12){
  if(data$Preference[n] == "Bitter"){
    bitter_pref_sweetener <- c(bitter_pref_sweetener, as.numeric(data[n, sweetener_vector]))
  }
}

bitter_pref_nosweetener <- c()
for(n in 1:12){
  if(data$Preference[n] == "Bitter"){
    bitter_pref_nosweetener <- c(bitter_pref_nosweetener, as.numeric(data[n, nosweetener_vector]))
  }
}

hist(bitter_pref_all, 100) ; text(8, 6, labels = paste('Mean = ', round(mean(bitter_pref_all), 4), sep = '\n'))
hist(bitter_pref_sweetener, 100); text(6, 4, labels = paste('Mean = ', round(mean(bitter_pref_sweetener), 4), sep = '\n'))
hist(bitter_pref_nosweetener, 100); text(5, 2.5, labels = paste('Mean = ', round(mean(bitter_pref_nosweetener), 4), sep = '\n'))
readline()

# The data is as expected. The mean is significantly lower for sweetened tea, as expected in the group
# who prefer bitter tastes. Taking into account the low sample sizes, the distributions seem fairly normal.

sweet_pref_sweetener <- c()
for(n in 1:12){
  if(data$Preference[n] == "Sweet"){
    sweet_pref_sweetener <- c(sweet_pref_sweetener, as.numeric(data[n, sweetener_vector]))
  }
}

sweet_pref_nosweetener <- c()
for(n in 1:12){
  if(data$Preference[n] == "Sweet"){
    sweet_pref_nosweetener <- c(sweet_pref_nosweetener, as.numeric(data[n, nosweetener_vector]))
  }
}

hist(sweet_pref_all, 100); text(3, 7, labels = paste('Mean = ', round(mean(sweet_pref_all), 4), sep = '\n'))
hist(sweet_pref_sweetener, 100); text(6, 7, labels = paste('Mean = ', round(mean(sweet_pref_sweetener), 4), sep = '\n'))
hist(sweet_pref_nosweetener, 100); text(3, 5, labels = paste('Mean = ', round(mean(sweet_pref_nosweetener), 4), sep = '\n'))
readline()

# And once again, the data is as expected. This group prefers tea that is sweetened, although their enjoyment
# of tea in general is much lower than the other two groups. The non-sweetened tea seems to have a sinc
# function-esque distribution.

# Finally, let's compare all of the histograms.

par(mfrow = c(3,3))

hist(no_pref_all, 100); text(2, 20, labels = paste('Mean = ', round(mean(no_pref_all), 4), sep = '\n'))
hist(no_pref_sweetener, 100); text(4, 12, labels = paste('Mean = ', round(mean(no_pref_sweetener), 4), sep = '\n'))
hist(no_pref_nosweetener, 100) ; text(3, 12, labels = paste('Mean = ', round(mean(no_pref_nosweetener), 4), sep = '\n'))
hist(bitter_pref_all, 100) ; text(8, 6, labels = paste('Mean = ', round(mean(bitter_pref_all), 4), sep = '\n'))
hist(bitter_pref_sweetener, 100); text(6, 4, labels = paste('Mean = ', round(mean(bitter_pref_sweetener), 4), sep = '\n'))
hist(bitter_pref_nosweetener, 100); text(5, 2.5, labels = paste('Mean = ', round(mean(bitter_pref_nosweetener), 4), sep = '\n'))
hist(sweet_pref_all, 100); text(3, 7, labels = paste('Mean = ', round(mean(sweet_pref_all), 4), sep = '\n'))
hist(sweet_pref_sweetener, 100); text(6, 7, labels = paste('Mean = ', round(mean(sweet_pref_sweetener), 4), sep = '\n'))
hist(sweet_pref_nosweetener, 100); text(3, 5, labels = paste('Mean = ', round(mean(sweet_pref_nosweetener), 4), sep = '\n'))
readline()

# The group averages act in accord with their taste preferences. In the bitter taste preference group, this
# factor was more significant than whether or not the tea is green/black. In the no taste preference, this
# factor was insignificant. In the sweet taste preference, this factor was about as significant as the
# choice between green and black tea.



# The next factor to examine is brewing time (6 vs 3 minutes).

# Make data of tea ratings for long and short brew times.

long_brew_data <- c(data[,2], data[,3], data[,4], data[,9], data[,10], data[,11], data[,12])
short_brew_data <- c(data[,5], data[,6], data[,7], data[,8], data[,13], data[,14], data[,15], data[,16])

# Make vectors to indicate column positions; easier to manipulate.
long_brew_vector <- c(2,3,4,9,10,11,12)
short_brew_vector <-c(5,6,7,8,13,14,15,16)

no_pref_long_brew <- c()
for(n in 1:12){
  if(data$Preference[n] == "No Preference"){
    no_pref_long_brew <- c(no_pref_long_brew, as.numeric(data[n, long_brew_vector]))
  }
}

no_pref_short_brew <- c()
for(n in 1:12){
  if(data$Preference[n] == "No Preference"){
    no_pref_short_brew <- c(no_pref_short_brew, as.numeric(data[n, short_brew_vector]))
  }
}

par(mfrow = c(1,3))

hist(no_pref_all, 100); text(2, 20, labels = paste('Mean = ', round(mean(no_pref_all), 4), sep = '\n'))
hist(no_pref_long_brew, 100); text(4, 10, labels = paste('Mean = ', round(mean(no_pref_long_brew), 4), sep = '\n'))
hist(no_pref_short_brew, 100) ; text(3, 12, labels = paste('Mean = ', round(mean(no_pref_short_brew), 4), sep = '\n'))
readline()

# Looks like the brewing time has almost no effect on the people who have no taste preference.

bitter_pref_long_brew <- c()
for(n in 1:12){
  if(data$Preference[n] == "Bitter"){
    bitter_pref_long_brew <- c(bitter_pref_long_brew, as.numeric(data[n, long_brew_vector]))
  }
}

bitter_pref_short_brew <- c()
for(n in 1:12){
  if(data$Preference[n] == "Bitter"){
    bitter_pref_short_brew <- c(bitter_pref_short_brew, as.numeric(data[n, short_brew_vector]))
  }
}


hist(bitter_pref_all, 100) ; text(8, 6, labels = paste('Mean = ', round(mean(bitter_pref_all), 4), sep = '\n'))
hist(bitter_pref_long_brew, 100); text(9, 2, labels = paste('Mean = ', round(mean(bitter_pref_long_brew), 4), sep = '\n'))
hist(bitter_pref_short_brew, 100); text(8, 4, labels = paste('Mean = ', round(mean(bitter_pref_short_brew), 4), sep = '\n'))
readline()

# In the group of people who have bitter taste preferences, brew time makes a huge difference in quality.

sweet_pref_long_brew <- c()
for(n in 1:12){
  if(data$Preference[n] == "Sweet"){
    sweet_pref_long_brew <- c(sweet_pref_long_brew, as.numeric(data[n, long_brew_vector]))
  }
}

sweet_pref_short_brew <- c()
for(n in 1:12){
  if(data$Preference[n] == "Sweet"){
    sweet_pref_short_brew <- c(sweet_pref_short_brew, as.numeric(data[n, short_brew_vector]))
  }
}

hist(sweet_pref_all, 100); text(3, 7, labels = paste('Mean = ', round(mean(sweet_pref_all), 4), sep = '\n'))
hist(sweet_pref_long_brew, 100); text(4, 4, labels = paste('Mean = ', round(mean(sweet_pref_long_brew), 4), sep = '\n'))
hist(sweet_pref_short_brew, 100); text(3, 5, labels = paste('Mean = ', round(mean(sweet_pref_short_brew), 4), sep = '\n'))
readline()

# In those who have sweet taste preferences, brew time is not significant.

# Finally, let's compare all of the histograms.


par(mfrow = c(3,3))

hist(no_pref_all, 100); text(2, 20, labels = paste('Mean = ', round(mean(no_pref_all), 4), sep = '\n'))
hist(no_pref_long_brew, 100); text(4, 10, labels = paste('Mean = ', round(mean(no_pref_long_brew), 4), sep = '\n'))
hist(no_pref_short_brew, 100) ; text(3, 12, labels = paste('Mean = ', round(mean(no_pref_short_brew), 4), sep = '\n'))
hist(bitter_pref_all, 100) ; text(8, 6, labels = paste('Mean = ', round(mean(bitter_pref_all), 4), sep = '\n'))
hist(bitter_pref_long_brew, 100); text(9, 2, labels = paste('Mean = ', round(mean(bitter_pref_long_brew), 4), sep = '\n'))
hist(bitter_pref_short_brew, 100); text(8, 4, labels = paste('Mean = ', round(mean(bitter_pref_short_brew), 4), sep = '\n'))
hist(sweet_pref_all, 100); text(3, 7, labels = paste('Mean = ', round(mean(sweet_pref_all), 4), sep = '\n'))
hist(sweet_pref_long_brew, 100); text(4, 4, labels = paste('Mean = ', round(mean(sweet_pref_long_brew), 4), sep = '\n'))
hist(sweet_pref_short_brew, 100); text(3, 5, labels = paste('Mean = ', round(mean(sweet_pref_short_brew), 4), sep = '\n'))
readline()

# The only group who exhibits a significant reaction to brewing time is the group who prefers bitter taste.
# This is probably because a longer brewing time results in a denser flavor profile, and since the bitterness
# in tea comes from the brewing of the leaves, this makes sense.



# The last factor to investigate is the mass of leaves used in the brew. This can be either 3 or 1.5 grams,
# high or low.

# Make data of tea ratings for high and low leaf mass.

high_mass_data <- c(data[,2], data[,4], data[,5], data[,8], data[,9], data[,11], data[,12], data[,13], 
                    data[,15], data[,16])
low_mass_data <- c(data[,3], data[,6], data[,7], data[,10], data[,14])

# Note, the amount of data of each factor is skewed more heavily in this case: high mass brew is represented
# significantly better than low mass brew.

# I get lazy here and skip to the final plot with all 9 histograms rather than looking at each one 1 by 1. Sorry.

# Make vectors to indicate column positions; easier to manipulate.
high_mass_vector <- c(2,4,5,8,9,11,12,13,15,16)
low_mass_vector <-c(3,6,7,10,14)

no_pref_high_mass <- c()
for(n in 1:12){
  if(data$Preference[n] == "No Preference"){
    no_pref_high_mass <- c(no_pref_high_mass, as.numeric(data[n, high_mass_vector]))
  }
}

no_pref_low_mass <- c()
for(n in 1:12){
  if(data$Preference[n] == "No Preference"){
    no_pref_low_mass <- c(no_pref_low_mass, as.numeric(data[n, low_mass_vector]))
  }
}

par(mfrow = c(1,3))

hist(no_pref_all, 100); text(2, 20, labels = paste('Mean = ', round(mean(no_pref_all), 4), sep = '\n'))
hist(no_pref_high_mass, 100); text(4, 14, labels = paste('Mean = ', round(mean(no_pref_high_mass), 4), sep = '\n'))
hist(no_pref_low_mass, 100) ; text(3, 7, labels = paste('Mean = ', round(mean(no_pref_low_mass), 4), sep = '\n'))

#readline()

# 

bitter_pref_high_mass <- c()
for(n in 1:12){
  if(data$Preference[n] == "Bitter"){
    bitter_pref_high_mass <- c(bitter_pref_high_mass, as.numeric(data[n, high_mass_vector]))
  }
}

bitter_pref_low_mass <- c()
for(n in 1:12){
  if(data$Preference[n] == "Bitter"){
    bitter_pref_low_mass <- c(bitter_pref_low_mass, as.numeric(data[n, low_mass_vector]))
  }
}


hist(bitter_pref_all, 100) ; text(8, 5, labels = paste('Mean = ', round(mean(bitter_pref_all), 4), sep = '\n'))
hist(bitter_pref_high_mass, 100); text(8, 4, labels = paste('Mean = ', round(mean(bitter_pref_high_mass), 4), sep = '\n'))
hist(bitter_pref_low_mass, 100); text(8, 3, labels = paste('Mean = ', round(mean(bitter_pref_low_mass), 4), sep = '\n'))

#readline()

# 

sweet_pref_high_mass <- c()
for(n in 1:12){
  if(data$Preference[n] == "Sweet"){
    sweet_pref_high_mass <- c(sweet_pref_high_mass, as.numeric(data[n, high_mass_vector]))
  }
}

sweet_pref_low_mass <- c()
for(n in 1:12){
  if(data$Preference[n] == "Sweet"){
    sweet_pref_low_mass <- c(sweet_pref_low_mass, as.numeric(data[n, low_mass_vector]))
  }
}

hist(sweet_pref_all, 100); text(3, 7, labels = paste('Mean = ', round(mean(sweet_pref_all), 4), sep = '\n'))
hist(sweet_pref_high_mass, 100); text(4, 6, labels = paste('Mean = ', round(mean(sweet_pref_high_mass), 4), sep = '\n'))
hist(sweet_pref_low_mass, 100); text(3, 4, labels = paste('Mean = ', round(mean(sweet_pref_low_mass), 4), sep = '\n'))

#readline()

# 

# Finally, let's compare all of the histograms.

par(mfrow = c(3,3))

hist(no_pref_all, 100); text(2, 20, labels = paste('Mean = ', round(mean(no_pref_all), 4), sep = '\n'))
hist(no_pref_high_mass, 100); text(4, 14, labels = paste('Mean = ', round(mean(no_pref_high_mass), 4), sep = '\n'))
hist(no_pref_low_mass, 100) ; text(3, 7, labels = paste('Mean = ', round(mean(no_pref_low_mass), 4), sep = '\n'))
hist(bitter_pref_all, 100) ; text(8, 5, labels = paste('Mean = ', round(mean(bitter_pref_all), 4), sep = '\n'))
hist(bitter_pref_high_mass, 100); text(8, 4, labels = paste('Mean = ', round(mean(bitter_pref_high_mass), 4), sep = '\n'))
hist(bitter_pref_low_mass, 100); text(8, 3, labels = paste('Mean = ', round(mean(bitter_pref_low_mass), 4), sep = '\n'))
hist(sweet_pref_all, 100); text(3, 7, labels = paste('Mean = ', round(mean(sweet_pref_all), 4), sep = '\n'))
hist(sweet_pref_high_mass, 100); text(4, 6, labels = paste('Mean = ', round(mean(sweet_pref_high_mass), 4), sep = '\n'))
hist(sweet_pref_low_mass, 100); text(3, 4, labels = paste('Mean = ', round(mean(sweet_pref_low_mass), 4), sep = '\n'))
readline()

# No significant effects. However, the second most significant factor affecting the no preference group.

# 650 lines of code later, we have determined the following:

#  - The mass of the tea leaves is essentially insignificant when compared to the other factors, but in the
#    no preferences group it is the second most significant factor.
#  - The brew time of the tea is significant for bitter taste preferences and insignificant otherwise.
#  - Whether or not sweetener is added is the most significant factor affecting bitter taste preferences,
#    and is slightly less significant than type of tea for sweet taste preferences. It is insignificant
#    for no taste preference.
#  - The type of tea is the most significant factor affecting tea scores across all three groups.
#  - In general, black tea is less popular than green tea.



# Now that we have these facts assembled, and lists of various data, it is time to combine the data in a way
# that lets us create probability distributions of tea scores for each taste preference. We will begin with
# the no preferences group, then the sweet preference group, then the bitter preference group.



# For the first group we will omit sweetener presence and brew time from our considerations for simplicity.

# No Preference - Probability distribution of scores - Tea type + Leaf mass

# We will create a contour plot of the count frequencies due to each parameter. Then, we will extract the
# most probable tea score by taking the indices where the maximum occurs, multiplying these indices
# together, and taking their nth root (n is the number of indices that specify maxima) to obtain the most
# probable tea score given these data and parameters.


# Count the occurences of each score in each parameter.
no_pref_green_freq <- count(no_pref_green)

# Fill in scores that got 0 counts to allow for contour plots.
if(min(no_pref_green_freq$x) != 1){
  for(n in (min(no_pref_green_freq$x) - 1):1){
    no_pref_green_freq <- rbind(c(n, 0), no_pref_green_freq)
  }
}
if(max(no_pref_green_freq$x) != 10){
  for(n in (max(no_pref_green_freq$x) + 1):10){
    no_pref_green_freq <- rbind(no_pref_green_freq, c(n, 0))
  }
}

no_pref_high_mass_freq <- count(no_pref_high_mass)
if(min(no_pref_high_mass_freq$x) != 1){
  for(n in (min(no_pref_high_mass_freq$x) - 1):1){
    no_pref_high_mass_freq <- rbind(c(n, 0), no_pref_high_mass_freq)
  }
}
if(max(no_pref_high_mass_freq$x) != 10){
  for(n in (max(no_pref_high_mass_freq$x) + 1):10){
    no_pref_high_mass_freq <- rbind(no_pref_high_mass_freq, c(n, 0))
  }
}

no_pref_low_mass_freq <- count(no_pref_low_mass)
if(min(no_pref_low_mass_freq$x) != 1){
  for(n in (min(no_pref_low_mass_freq$x) - 1):1){
    no_pref_low_mass_freq <- rbind(c(n, 0), no_pref_low_mass_freq)
  }
}
if(max(no_pref_low_mass_freq$x) != 10){
  for(n in (max(no_pref_low_mass_freq$x) + 1):10){
    no_pref_low_mass_freq <- rbind(no_pref_low_mass_freq, c(n, 0))
  }
}

# Separate the counts into separate vectors.
npg_counts <- no_pref_green_freq$freq
nphm_counts <- no_pref_high_mass_freq$freq
nplm_counts <- no_pref_low_mass_freq$freq

# Create a 1xn matrix and an nx1 matrix to allow for contour plotting.
npg_counts <- t(t(npg_counts))
nphm_counts <- matrix(nphm_counts, nrow = 1, ncol = length(nphm_counts))
nplm_counts <- matrix(nplm_counts, nrow = 1, ncol = length(nplm_counts))

# Create the z dimension.
GH_count_density <- npg_counts %*% nphm_counts
GL_count_density <- npg_counts %*% nplm_counts

# Normalize.
GH_count_density <- GH_count_density/sum(GH_count_density)
GL_count_density <- GL_count_density/sum(GL_count_density)

# Plot.
par(mfrow = c(2,2))
persp(no_pref_green_freq$x, no_pref_high_mass_freq$x, GH_count_density, theta = 37.5, phi = 37.5,
      xlab = 'Tea type score', ylab = 'Leaf mass score', zlab = 'Count density', main = "Green Tea
      High Leaf Mass")

GH_indices <- as.numeric(which(GH_count_density == max(GH_count_density), arr.ind = TRUE))
GH_result <- round((prod(GH_indices)^(1/length(GH_indices))), 3)

text(.10679, .29275, GH_result)

persp(no_pref_green_freq$x, no_pref_low_mass_freq$x, GL_count_density, theta = 37.5, phi = 37.5,
      xlab = 'Tea type score', ylab = 'Leaf mass score', zlab = 'Count density', main = "Green Tea
      Low Leaf Mass")

GL_indices <- as.numeric(which(GL_count_density == max(GL_count_density), arr.ind = TRUE))
GL_result <- round((prod(GL_indices)^(1/length(GL_indices))), 3)

text(.1472, .2945, GL_result)

# The most probable score is displayed on the persp peak.

# Now we will do the same for black tea.

no_pref_black_freq <- count(no_pref_black)
if(min(no_pref_black_freq$x) != 1){
  for(n in (min(no_pref_black_freq$x) - 1):1){
    no_pref_black_freq <- rbind(c(n, 0), no_pref_black_freq)
  }
}
if(max(no_pref_black_freq$x) != 10){
  for(n in (max(no_pref_black_freq$x) + 1):10){
    no_pref_black_freq <- rbind(no_pref_black_freq, c(n, 0))
  }
}

npb_counts <- no_pref_black_freq$freq

npb_counts <- t(t(npb_counts))

BH_count_density <- npb_counts %*% nphm_counts
BL_count_density <- npb_counts %*% nplm_counts

# Normalize
BH_count_density <- BH_count_density/sum(BH_count_density)
BL_count_density <- BL_count_density/sum(BL_count_density)

persp(no_pref_black_freq$x, no_pref_high_mass_freq$x, BH_count_density, theta = 37.5, phi = 37.5,
      xlab = 'Tea type score', ylab = 'Leaf mass score', zlab = 'Count density', main = "Black Tea
      High Leaf Mass")

BH_indices <- as.numeric(which(BH_count_density == max(BH_count_density), arr.ind = TRUE))
BH_result <- round((prod(BH_indices)^(1/length(BH_indices))), 3)

text(.10679, .29275, BH_result)

persp(no_pref_black_freq$x, no_pref_low_mass_freq$x, BL_count_density, theta = 37.5, phi = 37.5,
      xlab = 'Tea type score', ylab = 'Leaf mass score', zlab = 'Count density', main = "Black Tea
      Low Leaf Mass")

BL_indices <- as.numeric(which(BL_count_density == max(BL_count_density), arr.ind = TRUE))
BL_result <- round((prod(BL_indices)^(1/length(BL_indices))), 3)

text(.1472, .2945, BL_result)
readline()

# Interesting results! These data indicate that while the mean preference for high leaf mass is
# higher than the mean preference for low leaf mass, when taking into account the additional factor of green
# tea, the most likely score is higher in the low leaf mass. Why does this happen? There are fewer counts
# in the low leaf mass data, so each count has more weight. The few low scores have a strong effect on the
# mean of the lower leaf data. At the same time, a larger proportion of the data lies at high scores. In the
# high mass data, the maximum score represents 21% of the data, in the low mass data, the two maximum scores
# represent 45% of the data.

# As expected, the black tea max score is significantly lower than the green tea, but the relationship between
# the high and low leaf mass remains. Among people with no preference, low leaf mass is more likely to yield
# a higher score.

# Next we will examine the group who has sweet taste preferences.



# Sweet taste preferences - Probability distribution of scores - Tea type + sweetener

# Same approach as last time.

# Count the occurences of each score in each parameter.
sweet_pref_green_freq <- count(sweet_pref_green)

# Fill in scores that got 0 counts to allow for contour plots.
if(min(sweet_pref_green_freq$x) != 1){
  for(n in (min(sweet_pref_green_freq$x) - 1):1){
    sweet_pref_green_freq <- rbind(c(n, 0), sweet_pref_green_freq)
  }
}
if(max(sweet_pref_green_freq$x) != 10){
  for(n in (max(sweet_pref_green_freq$x) + 1):10){
    sweet_pref_green_freq <- rbind(sweet_pref_green_freq, c(n, 0))
  }
}

sweet_pref_sweetener_freq <- count(sweet_pref_sweetener)
if(min(sweet_pref_sweetener_freq$x) != 1){
  for(n in (min(sweet_pref_sweetener_freq$x) - 1):1){
    sweet_pref_sweetener_freq <- rbind(c(n, 0), sweet_pref_sweetener_freq)
  }
}
if(max(sweet_pref_sweetener_freq$x) != 10){
  for(n in (max(sweet_pref_sweetener_freq$x) + 1):10){
    sweet_pref_sweetener_freq <- rbind(sweet_pref_sweetener_freq, c(n, 0))
  }
}

sweet_pref_nosweetener_freq <- count(sweet_pref_nosweetener)
if(min(sweet_pref_nosweetener_freq$x) != 1){
  for(n in (min(sweet_pref_nosweetener_freq$x) - 1):1){
    sweet_pref_nosweetener_freq <- rbind(c(n, 0), sweet_pref_nosweetener_freq)
  }
}
if(max(sweet_pref_nosweetener_freq$x) != 10){
  for(n in (max(sweet_pref_nosweetener_freq$x) + 1):10){
    sweet_pref_nosweetener_freq <- rbind(sweet_pref_nosweetener_freq, c(n, 0))
  }
}

# Separate the counts into separate vectors.
spg_counts <- sweet_pref_green_freq$freq
sps_counts <- sweet_pref_sweetener_freq$freq
spn_counts <- sweet_pref_nosweetener_freq$freq

# Create a 1xn matrix and an nx1 matrix to allow for contour plotting.
spg_counts <- t(t(spg_counts))
sps_counts <- matrix(sps_counts, nrow = 1, ncol = length(sps_counts))
spn_counts <- matrix(spn_counts, nrow = 1, ncol = length(spn_counts))

# Create the z dimension.
GS_count_density <- spg_counts %*% sps_counts
GN_count_density <- spg_counts %*% spn_counts

# Normalize.
GS_count_density <- GS_count_density/sum(GS_count_density)
GN_count_density <- GN_count_density/sum(GN_count_density)

# Plot.
par(mfrow = c(2,2))
persp(sweet_pref_green_freq$x, sweet_pref_sweetener_freq$x, GS_count_density, theta = 37.5, phi = 37.5,
      xlab = 'Tea type score', ylab = 'Sweetener presence score', zlab = 'Count density', main = "Green Tea
      Sweetened")

GS_indices <- as.numeric(which(GS_count_density == max(GS_count_density), arr.ind = TRUE))
GS_result <- round((prod(GS_indices)^(1/length(GS_indices))), 3)

text(.10679, .29275, GS_result)

persp(sweet_pref_green_freq$x, sweet_pref_nosweetener_freq$x, GN_count_density, theta = 37.5, phi = 37.5,
      xlab = 'Tea type score', ylab = 'Sweetener presence score', zlab = 'Count density', main = "Green Tea
      Not Sweetened")

GN_indices <- as.numeric(which(GN_count_density == max(GN_count_density), arr.ind = TRUE))
GN_result <- round((prod(GN_indices)^(1/length(GN_indices))), 3)

text(.1472, .2945, GN_result)


sweet_pref_black_freq <- count(sweet_pref_black)
if(min(sweet_pref_black_freq$x) != 1){
  for(n in (min(sweet_pref_black_freq$x) - 1):1){
    sweet_pref_black_freq <- rbind(c(n, 0), sweet_pref_black_freq)
  }
}
if(max(sweet_pref_black_freq$x) != 10){
  for(n in (max(sweet_pref_black_freq$x) + 1):10){
    sweet_pref_black_freq <- rbind(sweet_pref_black_freq, c(n, 0))
  }
}

spb_counts <- sweet_pref_black_freq$freq

spb_counts <- t(t(spb_counts))

BS_count_density <- spb_counts %*% sps_counts
BN_count_density <- spb_counts %*% spn_counts

# Normalize
BS_count_density <- BS_count_density/sum(BS_count_density)
BN_count_density <- BN_count_density/sum(BN_count_density)

persp(sweet_pref_black_freq$x, sweet_pref_sweetener_freq$x, BS_count_density, theta = 37.5, phi = 37.5,
      xlab = 'Tea type score', ylab = 'Sweetner presence score', zlab = 'Count density', main = "Black Tea
      Sweetened")

BS_indices <- as.numeric(which(BS_count_density == max(BS_count_density), arr.ind = TRUE))
BS_result <- round((prod(BS_indices)^(1/length(BS_indices))), 3)

text(.10679, .29275, BS_result)

persp(sweet_pref_black_freq$x, sweet_pref_nosweetener_freq$x, BN_count_density, theta = 37.5, phi = 37.5,
      xlab = 'Tea type score', ylab = 'Sweetner presence score', zlab = 'Count density', main = "Black Tea
      Not Sweetened")

BN_indices <- as.numeric(which(BN_count_density == max(BN_count_density), arr.ind = TRUE))
BN_result <- round((prod(BN_indices)^(1/length(BN_indices))), 3)

text(.1472, .2945, BN_result)
readline()

# Once again, really interesting results. The group that says that they prefer sweet tastes have a much lower
# projected score for sweetened green tea than they do for unsweetened black tea. Also interesting is that the
# most likely score for both unsweetened green tea and sweetened black tea is the same.



# Finally, we have the bitter taste preference group. We have to analyze three different parameters:
#  - Tea type
#  - Brew time
#  - Sweetener presence

# Good luck.

# Count the occurences of each score in each parameter.


# Fill in scores that got 0 counts to allow for contour plots.

# Tea type
bitter_pref_green_freq <- count(bitter_pref_green)
if(min(bitter_pref_green_freq$x) != 1){
  for(n in (min(bitter_pref_green_freq$x) - 1):1){
    bitter_pref_green_freq <- rbind(c(n, 0), bitter_pref_green_freq)
  }
}

# hacking my way through the middle of the data frame. this fixes a bug
bitter_pref_green_freq <- rbind(bitter_pref_green_freq[1:5,], c(6,0), bitter_pref_green_freq[6:8,])
bitter_pref_green_freq <- rbind(bitter_pref_green_freq[1:7,], c(8,0), bitter_pref_green_freq[8:9,])

if(max(bitter_pref_green_freq$x) != 10){
  for(n in (max(bitter_pref_green_freq$x) + 1):10){
    bitter_pref_green_freq <- rbind(bitter_pref_green_freq, c(n, 0))
  }
}

bitter_pref_black_freq <- count(bitter_pref_black)
if(min(bitter_pref_black_freq$x) != 1){
  for(n in (min(bitter_pref_black_freq$x) - 1):1){
    bitter_pref_black_freq <- rbind(c(n, 0), bitter_pref_black_freq)
  }
}
if(max(bitter_pref_black_freq$x) != 10){
  for(n in (max(bitter_pref_black_freq$x) + 1):10){
    bitter_pref_black_freq <- rbind(bitter_pref_black_freq, c(n, 0))
  }
}

# Sweetener presence
bitter_pref_sweetener_freq <- count(bitter_pref_sweetener)
if(min(bitter_pref_sweetener_freq$x) != 1){
  for(n in (min(bitter_pref_sweetener_freq$x) - 1):1){
    bitter_pref_sweetener_freq <- rbind(c(n, 0), bitter_pref_sweetener_freq)
  }
}
if(max(bitter_pref_sweetener_freq$x) != 10){
  for(n in (max(bitter_pref_sweetener_freq$x) + 1):10){
    bitter_pref_sweetener_freq <- rbind(bitter_pref_sweetener_freq, c(n, 0))
  }
}

bitter_pref_nosweetener_freq <- count(bitter_pref_nosweetener)
if(min(bitter_pref_nosweetener_freq$x) != 1){
  for(n in (min(bitter_pref_nosweetener_freq$x) - 1):1){
    bitter_pref_nosweetener_freq <- rbind(c(n, 0), bitter_pref_nosweetener_freq)
  }
}
if(max(bitter_pref_nosweetener_freq$x) != 10){
  for(n in (max(bitter_pref_nosweetener_freq$x) + 1):10){
    bitter_pref_nosweetener_freq <- rbind(bitter_pref_nosweetener_freq, c(n, 0))
  }
}

# Brew time
bitter_pref_long_brew_freq <- count(bitter_pref_long_brew)
if(min(bitter_pref_long_brew_freq$x) != 1){
  for(n in (min(bitter_pref_long_brew_freq$x) - 1):1){
    bitter_pref_long_brew_freq <- rbind(c(n, 0), bitter_pref_long_brew_freq)
  }
}
if(max(bitter_pref_long_brew_freq$x) != 10){
  for(n in (max(bitter_pref_long_brew_freq$x) + 1):10){
    bitter_pref_long_brew_freq <- rbind(bitter_pref_long_brew_freq, c(n, 0))
  }
}

bitter_pref_short_brew_freq <- count(bitter_pref_short_brew)
if(min(bitter_pref_short_brew_freq$x) != 1){
  for(n in (min(bitter_pref_short_brew_freq$x) - 1):1){
    bitter_pref_short_brew_freq <- rbind(c(n, 0), bitter_pref_short_brew_freq)
  }
}
# this fixes another bug
bitter_pref_short_brew_freq <- rbind(bitter_pref_short_brew_freq[1:7,], c(8,0), bitter_pref_short_brew_freq[8:9,])

if(max(bitter_pref_short_brew_freq$x) != 10){
  for(n in (max(bitter_pref_short_brew_freq$x) + 1):10){
    bitter_pref_short_brew_freq <- rbind(bitter_pref_short_brew_freq, c(n, 0))
  }
}


# Separate the counts into separate vectors.
bpg_counts <- bitter_pref_green_freq$freq
bpb_counts <- bitter_pref_black_freq$freq
bps_counts <- bitter_pref_sweetener_freq$freq
bpn_counts <- bitter_pref_nosweetener_freq$freq
bplb_counts <- bitter_pref_long_brew_freq$freq
bpsb_counts <- bitter_pref_short_brew_freq$freq

# Choose one of each parameter to be the x, y, and z axis.
# Create a 3D matrix that contains our counting density.

# Lets start with green tea, sweet, long brew.
par(mfrow = c(2,2))

GSLB <- array(0, dim = c(10, 10, 10))
for(i in 1:10){
  for(j in 1:10){
    for(k in 1:10){
      GSLB[i,j,k] <- bpg_counts[i] * bps_counts[j] * bplb_counts[k]
    }
  }
}

# Normalize
GSLB <- GSLB/sum(GSLB)

GSLB_indices <- as.numeric(which(GSLB == max(GSLB), arr.ind = TRUE))
GSLB_result <- round((prod(GSLB_indices)^(1/length(GSLB_indices))), 3)

persp(bitter_pref_green_freq$x, bitter_pref_sweetener_freq$x, GSLB[,,GSLB_indices[3]], theta = 37.5,
      phi = 37.5, xlab = 'Tea type score', ylab = 'Sweetener presence score', 
      zlab = 'Count density', main = "Green Tea Sweetened \n Long Brew \n (shown at most probable long brew score)")

text(-.149, .2885, label = GSLB_result)

# So since we can't really get a true 4D plot, we make a 3D matrix of the count density, find the indices of
# maximum likelihood, and then show the count density at the slice that represents the most probable score
# due to a long brew.

# Now we will do the rest of the green tea combinations: sweet/short brew, non-sweet/long brew, 
# non-sweet/short brew.


# Green, sweet, short brew
GSSB <- array(0, dim = c(10, 10, 10))
for(i in 1:10){
  for(j in 1:10){
    for(k in 1:10){
      GSSB[i,j,k] <- bpg_counts[i] * bps_counts[j] * bpsb_counts[k]
    }
  }
}

GSSB <- GSSB/sum(GSSB)

GSSB_indices <- as.numeric(which(GSSB == max(GSSB), arr.ind = TRUE))
GSSB_result <- round((prod(GSSB_indices)^(1/length(GSSB_indices))), 3)

persp(bitter_pref_green_freq$x, bitter_pref_sweetener_freq$x, GSSB[,,GSSB_indices[3]], theta = 37.5,
      phi = 37.5, xlab = 'Tea type score', ylab = 'Sweetener presence score', 
      zlab = 'Count density', main = "Green Tea Sweetened \n Short Brew \n (shown at most probable short brew score)")

text(-.149, .2885, label = GSSB_result)


# Green, non-sweet, long brew
GNLB <- array(0, dim = c(10, 10, 10))
for(i in 1:10){
  for(j in 1:10){
    for(k in 1:10){
      GNLB[i,j,k] <- bpg_counts[i] * bpn_counts[j] * bplb_counts[k]
    }
  }
}

GNLB <- GNLB/sum(GNLB)

GNLB_indices <- as.numeric(which(GNLB == max(GNLB), arr.ind = TRUE))
GNLB_result <- round((prod(GNLB_indices)^(1/length(GNLB_indices))), 3)

persp(bitter_pref_green_freq$x, bitter_pref_nosweetener_freq$x, GNLB[,,GNLB_indices[3]], theta = 37.5,
      phi = 37.5, xlab = 'Tea type score', ylab = 'Sweetener presence score', 
      zlab = 'Count density', main = "Green Tea Not Sweetened \n Long Brew \n (shown at most probable long brew score)")

text(-.149, .2885, label = GNLB_result)

# Green, non-sweet, short brew
GNSB <- array(0, dim = c(10, 10, 10))
for(i in 1:10){
  for(j in 1:10){
    for(k in 1:10){
      GNSB[i,j,k] <- bpg_counts[i] * bpn_counts[j] * bpsb_counts[k]
    }
  }
}

GNSB <- GNSB/sum(GNSB)

GNSB_indices <- as.numeric(which(GNSB == max(GNSB), arr.ind = TRUE))
GNSB_result <- round((prod(GNSB_indices)^(1/length(GNSB_indices))), 3)

persp(bitter_pref_green_freq$x, bitter_pref_nosweetener_freq$x, GNSB[,,GNSB_indices[3]], theta = 37.5,
      phi = 37.5, xlab = 'Tea type score', ylab = 'Sweetener presence score', 
      zlab = 'Count density', main = "Green Tea Not Sweetened \n Short Brew \n (shown at most probable short brew score)")

text(-.149, .2885, label = GNSB_result)
readline()

# In the group preferring bitter tastes, unsweetened, long brew green tea seems to be the best choice out of
# the varieties of green tea. As expected, with a bitter taste preference, the sweetened teas are projected
# to score less.

# Now for black tea...

# Black tea, sweet, long brew
BSLB <- array(0, dim = c(10, 10, 10))
for(i in 1:10){
  for(j in 1:10){
    for(k in 1:10){
      BSLB[i,j,k] <- bpb_counts[i] * bps_counts[j] * bplb_counts[k]
    }
  }
}

BSLB <- BSLB/sum(BSLB)

BSLB_indices <- as.numeric(which(BSLB == max(BSLB), arr.ind = TRUE))
BSLB_result <- round((prod(BSLB_indices)^(1/length(BSLB_indices))), 3)

persp(bitter_pref_black_freq$x, bitter_pref_sweetener_freq$x, BSLB[,,BSLB_indices[3]], theta = 37.5,
      phi = 37.5, xlab = 'Tea type score', ylab = 'Sweetener presence score', 
      zlab = 'Count density', main = "Black Tea Sweetened \n Long Brew \n (shown at most probable long brew score)")

text(-.149, .2885, label = BSLB_result)

# Black tea, sweet, short brew
BSSB <- array(0, dim = c(10, 10, 10))
for(i in 1:10){
  for(j in 1:10){
    for(k in 1:10){
      BSSB[i,j,k] <- bpb_counts[i] * bps_counts[j] * bpsb_counts[k]
    }
  }
}

BSSB <- BSSB/sum(BSSB)

BSSB_indices <- as.numeric(which(BSSB == max(BSSB), arr.ind = TRUE))
BSSB_result <- round((prod(BSSB_indices)^(1/length(BSSB_indices))), 3)

persp(bitter_pref_black_freq$x, bitter_pref_sweetener_freq$x, BSSB[,,BSSB_indices[3]], theta = 37.5,
      phi = 37.5, xlab = 'Tea type score', ylab = 'Sweetener presence score', 
      zlab = 'Count density', main = "Black Tea Sweetened \n Short Brew \n (shown at most probable short brew score)")

text(-.149, .2885, label = BSSB_result)

# Black tea, non-sweet, long brew
BNLB <- array(0, dim = c(10, 10, 10))
for(i in 1:10){
  for(j in 1:10){
    for(k in 1:10){
      BNLB[i,j,k] <- bpb_counts[i] * bpn_counts[j] * bplb_counts[k]
    }
  }
}

BNLB <- BNLB/sum(BNLB)

BNLB_indices <- as.numeric(which(BNLB == max(BNLB), arr.ind = TRUE))
BNLB_result <- round((prod(BNLB_indices)^(1/length(BNLB_indices))), 3)

persp(bitter_pref_black_freq$x, bitter_pref_nosweetener_freq$x, BNLB[,,BNLB_indices[3]], theta = 37.5,
      phi = 37.5, xlab = 'Tea type score', ylab = 'Sweetener presence score', 
      zlab = 'Count density', main = "Black Tea Not Sweetened \n Long Brew \n (shown at most probable long brew score)")

text(-.149, .2885, label = BNLB_result)

# Black tea, non-sweet, short brew
BNSB <- array(0, dim = c(10, 10, 10))
for(i in 1:10){
  for(j in 1:10){
    for(k in 1:10){
      BNSB[i,j,k] <- bpb_counts[i] * bpn_counts[j] * bpsb_counts[k]
    }
  }
}

BNSB <- BNSB/sum(BNSB)

BNSB_indices <- as.numeric(which(BNSB == max(BNSB), arr.ind = TRUE))
BNSB_result <- round((prod(BNSB_indices)^(1/length(BNSB_indices))), 3)

persp(bitter_pref_black_freq$x, bitter_pref_nosweetener_freq$x, BNSB[,,BNSB_indices[3]], theta = 37.5,
      phi = 37.5, xlab = 'Tea type score', ylab = 'Sweetener presence score', 
      zlab = 'Count density', main = "Black Tea Not Sweetened \n Short Brew \n (shown at most probable short brew score)")

text(-.149, .2885, label = BNSB_result)

# Similar conclusion as above. In fact, the most popular tea among those preferring bitter tastes is projected
# to be black tea with no sweetener and long brewed. This makes sense, as out of any of the tea options this would
# result in the most bitter taste profile.

# And thus concludes our foray into the world of tea tasting. We have looked at the probability distributions
# of tea scores based on certain parameters, and we have answered many of the questions we set out to answer.

# Future directions:
#  - Figure out how to plot the 3D matrix, maybe some kind of time evolution or heatmap addition.
#  - More trials; but this is not something I could control.

# This was a really fun project, I wish there was a bit more structure to it, but I feel like I have learned
# a lot and have the ability to apply what I know to new and exciting situations. This was the best
# class I have taken in my life so far.

# Happy Holidays!