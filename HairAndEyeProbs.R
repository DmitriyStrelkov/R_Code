show( HairEyeColor ) # Show data
#Shows the builtin object HairEyeColor

EyeHairFreq <- apply( HairEyeColor, c("Eye","Hair"), sum ) # Sum across sex
#Sums the Male and Female versions of the distribution to produce one array

EyeHairProp <- EyeHairFreq / sum( EyeHairFreq ) # joint proportions,Table 4.1
#Takes each element in EyeHairFreq and divides it by the total elements in EyeHairFreq to produce
#an array that represents the proportionality of each possible result

show( round( EyeHairProp , 2 ) )
#Rounds the proportions to the nearest two decimal places

HairFreq <- apply( HairEyeColor , c("Hair") , sum ) # Sum across sex and eye
#Sums the occurences of every hair color across both genders and all eye colors

HairProp <- HairFreq / sum( HairFreq ) # marginal proportions,Table 4.1
#Takes each hair color and their occurences and divides by the total hair occurences
#to yield the marginal proportions of each hair color

show( round( HairProp , 2 ) )
#Rounds the marginal proportions of each hair color to two decimal places

EyeFreq <- apply( HairEyeColor , c("Eye") , sum ) # Sum across sex and eye
#Sums the occurences of every eye color across both genders and all hair colors

EyeProp <- EyeFreq / sum( EyeFreq ) # marginal proportions,Table 4.1
#Takes each eye color and their occurences and divides by the total eye occurences
#to yield the marginal proportions of each eye color

show( round( EyeProp , 2 ) )
#Rounds the marginal proportions of each eye color to two decimal places

EyeHairProp["Blue",] / EyeProp["Blue"] # conditional prob,Table 4.2
#Takes each value from the row "Blue" in EyeHairProp, and divides it by the value of "Blue" in EyeProp
#This creates a table of the conditional probability of any given hair color if eyes are blue.
#This happens because EyeHairProp["Blue"] calls the joint probability of each hair color and blue eyes,
#and EyeProp["Blue"] calls the marginal probability of blue eyes.
#Dividing these two yields the conditional probability.

#Task: Compute the conditional probability of all hair color given brown eyes and the conditional probability
#of all eye colors given brown hair.

CondProbGivenBrownEyes <- EyeHairProp["Brown",]/EyeProp["Brown"]
#For the conditional probability of any hair color given brown eyes, we take the joint probability of
#all hair colors and brown eyes and divide it by the marginal probability of brown eyes.

CondProbGivenBrownHair <- EyeHairProp[,"Brown"]/HairProp["Brown"]
#For the conditional probability of any eye color given brown hair, we take the joint probability of
#all eye colors and brown hair and divide it by the marginal probability of brown hair.
#By putting the comma before "Brown" in EyeHairProp[,"Brown"], we call the column(hair color) "Brown"
#rather than the row.
