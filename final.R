source("posterior_generator.R")
thetalist <- get_thetalist()
posterior_generator()

#data must be entered as one long string of 0's for tails and 1's for heads.