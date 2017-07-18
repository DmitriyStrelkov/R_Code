# MIGS: Multifunctional Integrated Gaussian beam Simulator
# 2D Gaussian Beam and Flow Cytometer Flow Cell Simulator
# Written in R by Dmitriy Strelkov

# All beam dimensions should be entered as diameters not as radii. They are labeled capital R.

# Definition of Gaussian beam. P is power in watts, R is 1/e^2 beam diameter in any unit, r is the chosen distance.
# If a list is input for r, a Gaussian beam intensity graph will be generated for the given beam parameters.
Gaussian <- function(P, R, r){
  Intensity <- (P/(pi * (.5 * R)^2 * 1/2)) * exp((-2*(r^2)/((.5 * R)^2)))
  return(Intensity)
}

# This function returns the beam size at any point in the beam given the wavelength, beam width(min spot size), and
# position. All distances should be in microns.
Beam_Propagation <- function(wavelength, beam_diameter, z){
  beam_waist <- .5 * beam_diameter
  Rayleigh_range <- (pi * beam_waist^2)/wavelength
  w_z <- beam_waist * sqrt(1 + (z/Rayleigh_range)^2)
  #print(paste('The Rayleigh Range is', round(Rayleigh_range * 10^-3, 5), "mm", sep = ' '))
  return(2 * w_z)
}


Real_Intensity_Matrix_Maker <- function(P, R_x, R_y, range){

    x_values <- seq(-R_x, R_x, by = R_x/100)
    y_values <- seq(-R_y, R_y, by = R_y/100)
    #range_values <- seq(-range, range, by = range/100)

  # Add zeros to values.
  if(!(is.element(0, x_values))){
    x_zero_index <- which(x_values == min(x_values[x_values > 0])) - 1
    
    x_values <- append(x_values, 0, x_zero_index)
    y_values <- append(y_values, max(y_values) + R_x/100, length(y_values))
    #range_values <- append(range_values, max(range_values) + R_x/100, length(range_values))
  }
  if(!(is.element(0, y_values))){
    y_zero_index <- which(y_values == min(y_values[y_values > 0])) - 1
    
    y_values <- append(y_values, 0, y_zero_index)
    x_values <- append(x_values, max(x_values) + R_y/100, length(x_values))
    #range_values <- append(range_values, max(range_values) + R_y/100, length(range_values))
  }
    # if(!(is.element(0, range_values))){
    #   range_zero_index <- which(range_values == min(range_values[range_values > 0])) - 1
    #
    #   range_values <- append(range_values, 0, range_zero_index)
    #   x_values <- append(x_values, max(x_values) + range/100, length(x_values))
    #   y_values <- append(y_values, max(y_values) + range/100, length(y_values))
    # }
    #
  x_intensity <- Gaussian(P, R_x, x_values)
  y_intensity <- Gaussian(P, R_y, x_values)

  I_0 <- (2/(pi * (R_x/2) * (R_y/2))) * P
  I_xyz <- matrix(0, ncol = length(x_values), nrow = length(y_values))
  # if(R_x > R_y){
  #   for(i in 1:length(x_values)){
  #     for(j in 1:length(y_values)){
  #       radius <- sqrt(x_values[i]^2 + y_values[j]^2)
  #       I_xyz[i,j] <- 
  #     }
  #   }
  # }

  x_start <- min(which(x_values > -range))
  x_end <- max(which(x_values < range))
  y_start <- min(which(y_values > -range))
  y_end <- max(which(y_values < range))


  if(R_x > R_y){
    for(i in x_start:x_end){
      for(j in y_start:y_end){
        I_xyz[i,j] <- I_0 * (exp((-2 * x_values[i]^2)/(.5 * R_x)^2) *
                               exp((-2 * x_values[j]^2)/(.5 * R_y)^2))
      }
    }
  } else if(R_y >= R_x){
    for(i in x_start:x_end){
      for(j in y_start:y_end){
        I_xyz[i,j] <- I_0 * (exp((-2 * y_values[i]^2)/(.5 * R_x)^2) *
                               exp((-2 * y_values[j]^2)/(.5 * R_y)^2))
      }
    }
  }
  return(I_xyz)
}


# This function generates an intensity matrix given power and beam dimensions then normalizes the matrix 
# such that the sum of all intensities is 1. This will help in determining what % of power is captured in a 
# certain area of the beam.
Sum_1_Intensity_Matrix_Maker <- function(P, R_x, R_y){
  x_values <- seq(-R_x, R_x, by = R_x/50)
  y_values <- seq(-R_y, R_y, by = R_y/50)
  
  # Add zeros to values.
  if(!(is.element(0, x_values))){
    x_zero_index <- which(x_values == min(x_values[x_values > 0])) - 1
    
    x_values <- append(x_values, 0, x_zero_index)
    y_values <- append(y_values, max(y_values) + R_y/50, length(y_values))
  }
  if(!(is.element(0, y_values))){
    y_zero_index <- which(y_values == min(y_values[y_values > 0])) - 1
    
    y_values <- append(y_values, 0, y_zero_index)
    x_values <- append(x_values, max(x_values) + R_x/50, length(x_values))
  }
  
  if(R_x > R_y){
    x_intensity <- Gaussian(P, R_x, x_values)
    y_intensity <- Gaussian(P, R_y, x_values)
    
    z_beam <- x_intensity %*% t(y_intensity)
    
  } else{
    x_intensity <- Gaussian(P, R_x, y_values)
    y_intensity <- Gaussian(P, R_y, y_values)
    
    z_beam <- x_intensity %*% t(y_intensity)
  }
  z_beam_sum_1 <- z_beam/sum(z_beam)
  return(z_beam_sum_1)
}

# This function generates an intensity matrix given power and beam dimensions, then 
# converts every value into the fraction of max power. This is useful for finding the % intensity at a given value.
Max_1_Intensity_Matrix_Maker <- function(P, R_x, R_y){
  Max_1_Matrix <- Real_Intensity_Matrix_Maker(P, R_x, R_y)
  Max_1_Matrix <- Max_1_Matrix/max(Max_1_Matrix)
  return(Max_1_Matrix)
}

# This function creates a matrix that represents the radius of any given xy position.
Radial_Matrix_Maker <- function(x, y, a, b){
  Radii <- matrix(0, nrow = length(x), ncol = length(y))
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      Radii[i,j] <- sqrt((x[i]^2)/a^2 + (y[j]^2)/b^2)
    }
  }
  return(Radii)
}

# This function plots a Gaussian beam profile of a 2D beam given the intensity matrix and the 1/e^2 beam diameter.
Beam_Plotter <- function(P, R_x, R_y){
  Intensity_Matrix <- Real_Intensity_Matrix_Maker(P, R_x, R_y)
  x_values <- seq(-R_x, R_x, by = R_x/50)
  y_values <- seq(-R_y, R_y, by = R_y/50)
  
  # Add zeros to values.
  if(!(is.element(0, x_values))){
    x_zero_index <- which(x_values == min(x_values[x_values > 0])) - 1
    
    x_values <- append(x_values, 0, x_zero_index)
    y_values <- append(y_values, max(y_values) + R_y/50, length(y_values))
  }
  if(!(is.element(0, y_values))){
    y_zero_index <- which(y_values == min(y_values[y_values > 0])) - 1
    
    y_values <- append(y_values, 0, y_zero_index)
    x_values <- append(x_values, max(x_values) + R_x/50, length(x_values))
  }
  
  if(R_x > R_y){
    persp(x_values, x_values, Intensity_Matrix, theta = 35, phi = 35, xlab = 'X position', ylab = "Y position",
          zlab = 'Intensity', main = 'Intensity Profile at Beam Waist')
  } else{
    persp(y_values, y_values, Intensity_Matrix, theta = 35, phi = 35, xlab = 'X position', ylab = "Y position",
                 zlab = 'Intensity', main = 'Intensity Profile at Beam Waist')
  }
}

# This function finds the peak-valley variation of the long axis in the beam. Inputs are the power,
# the beam dimensions, and the x dimension and y dimension of interrogation. Output is the ratio of the max and min.
Peak_Valley_Width_Interrogator <- function(P, R_x, R_y, x_dim, y_dim){
  Peak_Valley_Width_Ratio <- numeric()
  
  if(R_x > R_y){
    Peak_Valley_Width_Ratio <- Gaussian(P, R_x, 0)/Gaussian(P, R_x, (.5 * x_dim))
  } else if(R_y >= R_x){
    Peak_Valley_Width_Ratio <- Gaussian(P, R_y, 0)/Gaussian(P, R_y, (.5 * y_dim))
  }
  
  Peak_Valley_Percentage <- (Peak_Valley_Width_Ratio - 1) * 100
  
  long_axis_result <- paste('The peak to valley ratio of the long axis at the short axis maximum is', 
                            round(Peak_Valley_Percentage, 1), '%', sep = ' ')
  print(long_axis_result)
  return(Peak_Valley_Width_Ratio)
}

# This function shows the contained power of a given beam in a window.
# Core must be have an area that is greater than 4% of the area of the beam.
# +- 0.5 um resolution.
Contained_Power_Finder <- function(P, R_x, R_y, x_dim, y_dim){
  if(x_dim > y_dim){
    Window_Intensity_Matrix <- Real_Intensity_Matrix_Maker(P, R_x, R_y, x_dim)
  } else if(y_dim >= x_dim){
    Window_Intensity_Matrix <- Real_Intensity_Matrix_Maker(P, R_x, R_y, y_dim)
  }
  if(R_x > R_y){
    Real_Intensity_Matrix <- Real_Intensity_Matrix_Maker(P, R_x, R_y, R_x)
  } else if(R_y >= R_x){
    Real_Intensity_Matrix <- Real_Intensity_Matrix_Maker(P, R_x, R_y, R_y)
    }
  x_values <- seq(-R_x, R_x, by = R_x/100)
  y_values <- seq(-R_y, R_y, by = R_y/100)

  # Add zeros to values.
  # if(!(is.element(0, x_values))){
  #   x_zero_index <- which(x_values == min(x_values[x_values > 0])) - 1
  #   
  #   x_values <- append(x_values, 0, x_zero_index)
  #   y_values <- append(y_values, max(y_values) + R_y/101, length(y_values))
  # }
  # if(!(is.element(0, y_values))){
  #   y_zero_index <- which(y_values == min(y_values[y_values > 0])) - 1
  #   
  #   y_values <- append(y_values, 0, y_zero_index)
  #   x_values <- append(x_values, max(x_values) + R_x/101, length(x_values))
  # }
  
  radius_x_range <- seq(-x_dim, x_dim, by = x_dim/100)
  radius_y_range <- seq(-y_dim, y_dim, by = y_dim/100)

  if(R_x > R_y){
    Radii <- Radial_Matrix_Maker(radius_x_range, radius_x_range, (x_dim/2), (y_dim/2))
  } else if(R_y >= R_x){
    Radii <- Radial_Matrix_Maker(radius_y_range, radius_y_range, (x_dim/2), (y_dim/2))
  }
  Radii <- Radii * (R_x/x_dim)
  Window_Area_Matrix <- (Radii <= 1) * 1
  
  Contained_Intensity_Matrix <- Window_Intensity_Matrix * Window_Area_Matrix
  contained_power_fraction <- sum(Contained_Intensity_Matrix)/sum(Real_Intensity_Matrix)
  
  # Square boxes
  # x_dim_index <- which(abs(x_values - round(.5 * x_dim))==min(abs(x_values - round(.5 * x_dim))))
  # y_dim_index <- which(abs(y_values - round(.5 * y_dim))==min(abs(y_values - round(.5 * y_dim))))
  # 
  # max_index <- which(Real_Intensity_Matrix == max(Real_Intensity_Matrix), arr.ind = TRUE)[1]
  # 
  # x_dim_range <- abs(x_dim_index - max_index)
  # y_dim_range <- abs(y_dim_index - max_index)
  # 
  # contained_power_fraction <- sum(Real_Intensity_Matrix[(max_index - x_dim_range):(max_index + x_dim_range),
  #                                         (max_index - y_dim_range):(max_index + y_dim_range)])/sum(Real_Intensity_Matrix)
  
  contained_power_value <- contained_power_fraction * P
  on_axis_intensity_max <- (2 * P)/((.5 * R_x * 10^-4) * (.5 * R_y * 10^-4) * pi)
  message(paste('The fraction of total power contained in the window is', contained_power_fraction, sep = ' '))
  message(paste('The total power contained in the window is', contained_power_value, 'Watts', sep = ' '))
  message(paste('The maximum on axis intensity is', on_axis_intensity_max, 'W/cm^2',
              sep = ' '))
  return(contained_power_fraction)
}

# This function finds the ratio between the power of the beam in the center of the flow channel and at the back side
# of the flow channel. Inputs are wavelength, power, beam dimensions at focus, back to center distance of the flow channel,
# and window dimensions. Output is the ratio between the power in the window at the center of the  channel and the 
# power in the window at the back of the channel. All lengths should be in microns.

# This function returns the coefficient of variation.
CV_Finder <- function(P, R_x, R_y){
  Intensity_Matrix <- Real_Intensity_Matrix_Maker(P, R_x, R_y, R_x)
  CV <- sd(Intensity_Matrix)/mean(Intensity_Matrix)
  result <- paste("The coefficient of variation is", CV, sep = ' ')
  return(result)
}

# This function returns the changing intensity profile as the beam propagates along the optic axis. The intensity
# is shown along the x and z dimensions; the y dimension is assumed to be at the intensity maximum (y = 0).

Propagation_Intensity_Plotter <- function(P, wavelength, R_x, R_y, z, window){
  x_values <- seq(-R_x, R_x, by = R_x/101)
  z_values <- seq(-z, z, by = z/101)
  w_z <- 0
  
  # Add zeros to values.
  if(!(is.element(0, x_values))){
    x_zero_index <- which(x_values == min(x_values[x_values > 0])) - 1
    
    x_values <- append(x_values, 0, x_zero_index)
    z_values <- append(z_values, max(z_values) + z/101, length(z_values))
  }
  
  x_intensity <- Gaussian(P, R_x, x_values)
  y_intensity <- Gaussian(P, R_y, 0)
  
  
  Intensity_Matrix <- matrix(0, nrow= length(z_values), ncol = length(x_values))
  I_0 <- max(x_intensity) + y_intensity

  for(i in 1:length(z_values)){
    for(j in 1:length(x_values)){
      w_z <- Beam_Propagation(wavelength, R_x, z_values[i])
      Intensity_Matrix[i,j] <- I_0 * (R_x/w_z)^2 * (exp((-2 * x_values[j]^2)/(.5 * R_x)^2))
    }
  }
  Intensity_Matrix <- Intensity_Matrix/max(Intensity_Matrix)
  filled.contour(z_values, x_values, Intensity_Matrix, color = topo.colors, 
                 main = 'Intensity Heatmap Through Flow Cell',
                 xlab = 'Position in Flow Channel', ylab = 'Horizontal Beam Location')

  # Zoom in on window.
  window_values <- seq(-(.5 * window), (.5 * window), by = window/194)
  if(!(is.element(0, window_values))){
    window_zero_index <- which(window_values == min(window_values[window_values > 0])) - 1
    
    window_values <- append(window_values, 0, window_zero_index)
  }
  window_start <- max(which(x_values < min(window_values)))
  window_end <- min(which(x_values > max(window_values)))
  
  filled.contour(z_values, x_values[window_start:window_end], Intensity_Matrix[,window_start:window_end],
                 color = topo.colors, 
                 main = 'Intensity Heatmap Through Flow Cell',
                 xlab = 'Position in Flow Channel', ylab = 'Horizontal Beam Location')
}

# This function integrates the beam plotter, the peak to valley finder, the coefficient of 
# variation finder, the contained power finder, and the propagation intensity plotter. Inputs 
# are power(Watts), e^-2 beam dimensions(microns), the interrogation window dimensions(microns),
# wavelength(microns), and z distance from center to edge of flow cell(microns).
Integrated_Function <- function(P, wavelength, R_x, R_y, x_dim, y_dim, z){
  Real_Matrix <- Real_Intensity_Matrix_Maker(P, R_x, R_y)
  window <- sqrt(x_dim * y_dim)
  
  Beam_Plotter(Real_Matrix, R_x, R_y)
  print(Peak_Valley_Width_Interrogator(P, R_x, R_y, x_dim, y_dim))
  print(CV_Finder(Real_Matrix))
  Contained_Power_Finder(P, R_x, R_y, x_dim, y_dim)
  Propagation_Intensity_Plotter(P, wavelength, R_x, R_y, z, window)
}

# Everything in microns
Beam_Profile_Through_Lens <- function(wavelength, M2, f, D){
  min_spot_size <- (4 * wavelength * M2 * f)/(pi * D)
  return(min_spot_size)
}

# This function calculates the photon flux at which a particular fluorophore saturates. MESF is fluors present on the 
# particle, mol_abs_coeff is the molar absorbtivity coefficient (epsilon), lifetime is the fluorescence lifetime of the 
# fluor in ns, QY is quantum yield of the fluor, P is power in Watts, R_x is the beam diameter in the x axis, R_y is the 
# beam diameter in the y axis, and particle_diameter is the particle diameter in the same units as the beam profile. 
# Wavelength is in microns. Flow rate in meters/sec.

# Function is tested with AF488 as the standard. The variable values are as follows:
# (270000, 73000, 4.1, .92, .1, Contained_Power_Finder(.1, 40, 10, 7.65, 7.65), 40, 10, 7.65, .488, 4)
Photon_Flux_State_Population <- function(
  MESF, mol_abs_coeff, lifetime, QY, P, Contained_Power, R_x, R_y, particle_diameter, wavelength, flow_rate){
  
  # Only one side of the particle is seen by the laser. For a spherical particle, we assume that only half of the 
  # fluors attached to a particle can see the laser.
  # UPDATE: is this valid? No, the MESF value already accounts for this. Corrected to remove this.
  
  # To determine how many photons are absorbed by the particle, we need to know the scattering cross section of 
  # the fluor. Units are cm squared.
  sigma <- log(10) * ((10^3)/(6.023 * 10^23)) * mol_abs_coeff 
  
  # Now we have to figure out the maximum photon flux (particle center illuminated by center of beam).
  power_on_target_fraction <- Contained_Power
  power_on_target_true <- P * power_on_target_fraction
  
  # Now we can know how many photons are incident per second.
  h <- 6.626 * 10^-34 # Joule * sec
  c <- 3 * 10^8 # meters/sec
  lambda <- wavelength * 10^-6 # meters
  
  photon_energy <- (h * c)/lambda
  
  photons_per_second_on_target <- power_on_target_true/photon_energy
  
  # Nominally, the flow rate of the NxT is 4 m/s. Figure out how long the particle spends in the beam.
  
  t_illuminated <- R_y/(flow_rate * 10^6)
  
  # Now, assume that all of the fluors are in the ground state before illumination. The time it takes for a 
  # population of excited molecules to decay to 1/e is the lifetime of the fluor. The excitation/decay cycle will happen 
  # every lifetime. How many lifetimes are in the illumination time?
  
  decay_cycle_count <- floor((t_illuminated)/(lifetime * 10^-9))
  
  # Now we simulate our 'leaky barrel' that is our fluor population. The amount of fluors that are absorbed become the 
  # excited state population at the beginning of the cycle. We lose 63.2% of our excited state population every fluorescence 
  # lifetime. The quantum yield of the fluor determines how much of this decay results in the emission of a photon.
  
  # We will set our populations now. We do care about ground state depletion, but we care more about the rate of change 
  # of the excited state. When the promotion to the excited state and decay to the ground state reach an equilibrium, 
  # we can say that a certain level of saturation has been achieved.
  
  # Set ground/excited state pop to zero. Set initial photons absorbed/emitted to zero.
  N_1_current_cycle <- MESF
  N_2_current_cycle <- 0
  photons_absorbed_record <- c()
  photons_emitted_record <- c()
  N_1_record_ex <- c()
  N_2_record_ex <- c()
  N_1_record_decay <- c()
  N_2_record_decay <- c()
  
  # Scale down the amount of lifetimes counted to save computing time. Save the scaling factor to later compute 
  # power out.
  lifetime_scaling <- numeric()
  if(decay_cycle_count > 500){
    lifetime_scaling <- decay_cycle_count/500
    decay_cycle_count <- 500
  }
  for(cycle in 1:decay_cycle_count){
    # We are going to calculate and record the rate of change of the ground/excited state population.
    
    # Calculate how many photons get absorbed. Photon count times ground state pop per square cm times the cross 
    # section (chance of absorption).
    
    photons_absorbed_current_lifetime  <- round(photons_per_second_on_target * 
      (N_1_current_cycle/(pi * (.5 * particle_diameter * 10^-4)^2)) * sigma * lifetime * 10^-9)

    photons_absorbed_record <- c(photons_absorbed_record, photons_absorbed_current_lifetime)
    
    # Calculate new ground/excited state population.
    N_1_current_cycle <- N_1_current_cycle - round(photons_absorbed_current_lifetime)
    N_2_current_cycle <- N_2_current_cycle + round(photons_absorbed_current_lifetime)

    # Determine how many photons are emitted, then determine the population of the excited state at 
    # the end of the fluorescence lifetime. This vector will keep track of the photons emitted at each cycle.
    photons_emitted_current_cycle <- round((1 - exp(-1)) * N_2_current_cycle * QY)
    photons_emitted_record <- c(photons_emitted_record, photons_emitted_current_cycle)
    
    # Append the ground/excited state population to the record. These numbers represent the population of the energy 
    # states before decay.
    N_2_record_ex <- c(N_2_record_ex, N_2_current_cycle)
    N_1_record_ex <- c(N_1_record_ex, N_1_current_cycle)
    
    
    # Determine how the excited state population at the end of the cycle. Round again.
    N_1_current_cycle <- round(N_1_current_cycle + (N_2_current_cycle * (1 - exp(-1))))
    N_2_current_cycle <- round(N_2_current_cycle * exp(-1))

    
    # Append the ground/excited state population to the record. These numbers represent the population of the energy 
    # states after decay.
    N_2_record_decay <- c(N_2_record_decay, N_2_current_cycle)
    N_1_record_decay <- c(N_1_record_decay, N_1_current_cycle)
    
  }
  # Make some plots of the cycles.
  # plot(c(1:decay_cycle_count), photons_absorbed_record, xlab = 'Cycle', ylab = 'Photons absorbed')
  # plot(c(1:decay_cycle_count), photons_emitted_record, xlab = 'Cycle', ylab = 'Photons emitted')
  # plot(c(1:decay_cycle_count), N_2_record_ex, xlab = 'Cycle', ylab = 'Excited state population before decay')
  # plot(c(1:decay_cycle_count), N_1_record_ex, xlab = 'Cycle', ylab = 'Ground state population before decay')
  # plot(c(1:decay_cycle_count), N_2_record_decay, xlab = 'Cycle', ylab = 'Excited state population after decay')
  # plot(c(1:decay_cycle_count), N_1_record_decay, xlab = 'Cycle', ylab = 'Ground state population after decay')
  
  # Power out in watts.
  P_out <- (sum(photons_emitted_record)* lifetime_scaling * photon_energy)/(t_illuminated)
  
  # Average number of ground/excited state photons. Since the population dynamics quickly reach steady state 
  # (within a few lifetimes), we take the mean of the maximum and minimum populations of each state. That is, we 
  # take the mean of the population at the 'start' of the lifetime and the 'end' of the lifetime in order to 
  # arrive at the average population of that state.
  avg_ex_state_pop <- floor((max(N_2_record_ex) + max(N_2_record_decay))/2)
  avg_g_state_pop <- ceiling((min(N_1_record_ex) + min(N_1_record_decay))/2)
  if(avg_g_state_pop < 0){
    stop("Power too high.")
  }
  avg_photons_absorbed_per_lifetime <- round(photons_per_second_on_target * 
                                               (avg_g_state_pop/(pi * (.5 * particle_diameter * 10^-4)^2)) *
                                               sigma * lifetime * 10^-9)
  avg_photons_emitted_per_lifetime <- round(avg_photons_absorbed_per_lifetime * QY)
  
  # Ratio of ground state and excited state populations.
  g_ex_pop_ratio <- (avg_g_state_pop/avg_ex_state_pop)
  
  # Ng_1/P_1 = s * Ng_2/P_2
  
  # We will use the following value to show ground state population with excitation intensity.
  photon_flux_cm_squared <- (photons_per_second_on_target * 10^8)/power_on_target_fraction
  
  # Going to return the results that are relevant to us.
  results <- c(photon_flux_cm_squared,
               avg_photons_absorbed_per_lifetime, avg_photons_emitted_per_lifetime, avg_ex_state_pop,
               avg_g_state_pop, g_ex_pop_ratio, P_out) # 7 params
  # max(N_2_record_ex), min(N_1_record_ex),
  # max(N_2_record_decay), min(N_1_record_decay)
  return(results)
}

# This function tells us the fractional relationship between signal gain and power increase. 
Signal_Power_Gain_Max_Set <- function(MESF, mol_abs_coeff, lifetime, QY, P, Contained_Power,
                              R_x, R_y, particle_diameter, wavelength, flow_rate){
  # We will use the ratio of the state populations multiplied by the average photons emitted to create a 
  # measure of return on power investment.
  
  # First we have to find the power where avg photons emitted multiplied by state pop ratio is at a max. This 
  # will set our bar for 100% return on power. We'll use a metropolis algorithm to do so. 1000 iterations should 
  # suffice.
  P_current <- P
  P_test <- .1 * P
  efficiency_value_record <- c()
  efficiency_value_current <- suppressMessages(
    prod(Photon_Flux_State_Population(MESF, mol_abs_coeff, lifetime, QY, P_current, Contained_Power, 
                                      R_x, R_y, particle_diameter, wavelength, flow_rate)[c(3,6)]))
  
  # Select the photons emitted and the state ratios, and find the maximum efficiency value.
  for(i in 1:500){
    
    efficiency_value_test <- suppressMessages(
      prod(Photon_Flux_State_Population(MESF, mol_abs_coeff, lifetime, QY, P_test, Contained_Power,
                                        R_x, R_y, particle_diameter, wavelength, flow_rate)[c(3,6)]))
  
      if(efficiency_value_test > efficiency_value_current){
        
        efficiency_value_current <- efficiency_value_test
        P_current <- P_test
        P_test <- rnorm(1, mean = 1, sd = .05) * P_current
        efficiency_value_record <- c(efficiency_value_record, efficiency_value_current)
        
      } else if(efficiency_value_current >= efficiency_value_test){
        
        move_probability <- efficiency_value_test/(efficiency_value_current + efficiency_value_test)
        roll <- sample(c(0,1), 1, prob = c((1 - move_probability), move_probability))
        
        if(roll == 1){
          
          efficiency_value_current <- efficiency_value_test
          P_current <- P_test
          P_test <- rnorm(1, mean = 1, sd = .05) * P_current
          efficiency_value_record <- c(efficiency_value_record, efficiency_value_current)
          
        } else if(roll == 0){
          
          P_test <- rnorm(1, mean = 1, sd = .05) * P_current
          efficiency_value_record <- c(efficiency_value_record, efficiency_value_current)
        }
      }
    }
  efficiency_value_record <- hist(efficiency_value_record, 500)
  most_freq_value_index <- which(max(efficiency_value_record$counts) == efficiency_value_record$counts, arr.ind = T)[1]
  max_efficiency_value <- efficiency_value_record$breaks[most_freq_value_index]
  return(max_efficiency_value)
}

Signal_Power_Gain <- function(MESF, mol_abs_coeff, lifetime, QY, P, Contained_Power,
                              R_x, R_y, particle_diameter, wavelength, flow_rate, max_efficiency){
  efficiency <- (suppressMessages(
    prod(Photon_Flux_State_Population(MESF, mol_abs_coeff, lifetime, QY, P, Contained_Power, 
                                      R_x, R_y, particle_diameter, wavelength, flow_rate)[c(3,6)])))/max_efficiency
  
  return(efficiency)
}

Fluor_Sat_Prediction_Plotter <- function(MESF, mol_abs_coeff, lifetime, QY, P_max, R_x, R_y, particle_diameter,
                                         wavelength, flow_rate){
  P <- seq(.001, (P_max + .001), by = .05)
  cp <- suppressMessages(Contained_Power_Finder(.1, R_x, R_y, particle_diameter, particle_diameter))
  max_eff <- Signal_Power_Gain_Max_Set(MESF, mol_abs_coeff, lifetime, QY, .1, cp, R_x, R_y, particle_diameter,
                                       wavelength, flow_rate)
  
  while(max_eff > MESF){
    max_eff <- Signal_Power_Gain_Max_Set(MESF, mol_abs_coeff, lifetime, QY, .1, cp, R_x, R_y, particle_diameter,
                                         wavelength, flow_rate)
  }
  
  plot(1, xlim = c(0, (P_max + .1 * P_max)), ylim = c(0, 1), type = 'n', xlab = 'Power (W)', ylab = 'Return on Power Input',
       main = 'Fluor Power Input Efficiency')
  RoPI_record <- c()
  
  for(i in P){
    efficiency <- Signal_Power_Gain(MESF, mol_abs_coeff, lifetime, QY, i, cp, R_x, R_y, particle_diameter,
                                    wavelength, flow_rate, max_eff)
    RoP <- efficiency
    RoPI_record <- c(RoPI_record, RoP)
    points(i, RoP, pch = 19)
  }
  
  results <- data.frame(P, RoPI_record)
  return(results)
}

# Because the signal will be scaled up to fit user data, the MESF value will not actually affect the plotted signal.
Fluor_Signal_Prediction_Plotter <- function(MESF, mol_abs_coeff, lifetime, QY, P_max, R_x, R_y, particle_diameter,
                                            wavelength, flow_rate){
  P <- seq(.001, (P_max + .001), by = .05)
  cp <- suppressMessages(Contained_Power_Finder(.1, R_x, R_y, particle_diameter, particle_diameter))
  signal <- c(0)
  signal_total <- c()
  RoPI <- Fluor_Sat_Prediction_Plotter(MESF, mol_abs_coeff, lifetime, QY, P_max, R_x, R_y, particle_diameter,
                                       wavelength, flow_rate)$RoPI_record
  
  plot(1, xlim = c(0, (P_max + .1 * P_max)), ylim = c(0, MESF), type = 'n', xlab = 'Power (W)', ylab = 'Signal',
       main = 'Predicted Fluor Signal')
  j <- 1
  for(i in P){
    signal <- c(signal,
                (Photon_Flux_State_Population(MESF, mol_abs_coeff, lifetime, QY, i, cp, R_x, R_y, particle_diameter,
                                              wavelength, flow_rate)[3] - 
                   signal[j])* RoPI[j])
    signal_total <- c(signal_total, sum(signal))
    j <- j + 1
  }
  
  true_signal_power_index <- (as.numeric(
    readline('Enter the power (in mW) of the signal you will use as a reference (multiple of 50).'))/50) + 1
  true_signal <- as.numeric(readline('Enter the signal at that laser power.'))
  signal_scaling_factor <- true_signal/signal_total[true_signal_power_index]
  
  i <- 2
  while(RoPI[i - 1] >= .5){
    signal_total[i] <- signal_total[i] * signal_scaling_factor * RoPI[i - 1]
    i <- i + 1
  }
  for(i in i:length(signal_total)){
    signal_total[i] <- signal_total[i - 1] * (((P[i]/P[i - 1] - 1) * RoPI[i]) + 1)
  }
  
  growth <- c(0, 0)
  for(j in 2:(length(signal_total) - 1)){
    growth[j + 1] <- signal_total[j + 1]/signal_total[j]
  }
  
  saturation_fraction <- signal_total/max(signal_total)
  
  P_signal_table <- data.frame(P, signal_total, growth, saturation_fraction)
  points(P_signal_table$P, P_signal_table$signal_total, pch = 19)
  return(P_signal_table)
}


# Fluor specific plots.

AF488_Sat_Plotter <- function(){
  P <- seq(.001, 1.051, by = .05)
  cp <- suppressMessages(Contained_Power_Finder(.1, 40, 10, 7.65, 7.65))
  max_eff <- Signal_Power_Gain_Max_Set(270000, 73000, 4.1, .92, .1, cp, 40, 10, 7.65, .488, 4)
  
  while(max_eff > 270000){
    max_eff <- Signal_Power_Gain_Max_Set(270000, 73000, 4.1, .92, .1, cp, 40, 10, 7.65, .488, 4)
  }
  
  plot(1, xlim = c(0, 1), ylim = c(0, 1), type = 'n', xlab = 'Power (W)', ylab = 'Return on Power Input',
       main = 'AF488')
  RoPI_record <- c()
  
  for(i in P){
    efficiency <- Signal_Power_Gain(270000, 73000, 4.1, .92, i, cp, 40, 10, 7.65, .488, 4, max_eff)
    RoP <- efficiency
    RoPI_record <- c(RoPI_record, RoP)
    points(i, RoP, pch = 19)
  }
  
  results <- data.frame(P, RoPI_record)
  return(results)
}

# Because the signal will be scaled up to fit user data, the MESF value will not actually affect the plotted signal.
AF488_Signal_Plotter <- function(){
  P <- seq(.001, 1.051, by = .05)
  cp <- suppressMessages(Contained_Power_Finder(.1, 40, 10, 7.65, 7.65))
  signal <- c(0)
  signal_total <- c()
  RoPI <- AF488_Sat_Plotter()$RoPI_record
  
  plot(1, xlim = c(0, 1), ylim = c(0, 270000), type = 'n', xlab = 'Power (W)', ylab = 'Signal',
       main = 'AF488')
  j <- 1
  for(i in P){
    signal <- c(signal,
                (Photon_Flux_State_Population(270000, 73000, 4.1, .92, i, cp, 40, 10, 7.65, .488, 4)[3] - 
                   signal[j])* RoPI[j])
    signal_total <- c(signal_total, sum(signal))
    j <- j + 1
  }
  
  # true_signal <- as.numeric(readline('Enter the fluorescent signal value at 50 mW.'))
  # signal_scaling_factor <- true_signal/signal_total[2]
  signal_scaling_factor <- 1

  i <- 2
  while(RoPI[i - 1] >= .5){
    signal_total[i] <- signal_total[i] * signal_scaling_factor * RoPI[i - 1]
    i <- i + 1
  }
  for(i in i:length(signal_total)){
    signal_total[i] <- signal_total[i - 1] * (((P[i]/P[i - 1] - 1) * RoPI[i]) + 1)
  }
  
  growth <- c(0, 0)
  for(j in 2:(length(signal_total) - 1)){
    growth[j + 1] <- signal_total[j + 1]/signal_total[j]
  }
  
  saturation_fraction <- signal_total/max(signal_total)

  P_signal_table <- data.frame(P, signal_total, growth, saturation_fraction)
  points(P_signal_table$P, P_signal_table$signal_total, pch = 19)
  return(P_signal_table)
}



# AF 532 - see if it works! Same MESF, QY = .61, lifetime = 2.5 ns, wavelength = .532, epsilon = 81000
AF532_Sat_Prediction_Plotter <- function(){
  P <- seq(.001, 1.351, by = .05)
  cp <- suppressMessages(Contained_Power_Finder(.1, 40, 10, 7.65, 7.65))
  max_eff <- Signal_Power_Gain_Max_Set(270000, 81000, 2.5, .61, .1, cp, 40, 10, 7.65, .532, 4)
  
  while(max_eff > 270000){
    max_eff <- Signal_Power_Gain_Max_Set(270000, 81000, 2.5, .61, .1, cp, 40, 10, 7.65, .532, 4)
  }
  
  plot(1, xlim = c(0, 1.4), ylim = c(0, 1), type = 'n', xlab = 'Power (W)', ylab = 'Return on Power Input',
       main = 'AF532')
  RoPI_record <- c()
  for(i in P){
    efficiency <- Signal_Power_Gain(270000, 81000, 2.5, .61, i, cp, 40, 10, 7.65, .532, 4, max_eff)
    RoP <- efficiency
    RoPI_record <- c(RoPI_record, RoP)
    points(i, RoP, pch = 19)
  }
  results <- data.frame(P, RoPI_record)
  return(results)
}

AF532_Signal_Prediction_Plotter <- function(){
  P <- seq(.001, 1.351, by = .05)
  cp <- suppressMessages(Contained_Power_Finder(.1, 40, 10, 7.65, 7.65))
  signal <- c(0)
  signal_total <- c()
  RoPI <- AF532_Sat_Prediction_Plotter()$RoPI_record
  plot(1, xlim = c(0, 1.4), ylim = c(0, 270000), type = 'n', xlab = 'Power (W)', ylab = 'Signal',
       main = 'AF532')
  j <- 1
  for(i in P){
    signal <- c(signal,
                (Photon_Flux_State_Population(270000, 81000, 2.5, .61, i, cp, 40, 10, 7.65, .532, 4)[3] - 
                   signal[j])* RoPI[j])
    signal_total <- c(signal_total, sum(signal))
    j <- j + 1
  }
  
  # Since we don't have any data, we will use a signal scaling factor of 1
  # true_signal <- as.numeric(readline('Enter the fluorescent signal value at 50 mW.'))
  # signal_scaling_factor <- true_signal/signal_total[2]
  
  signal_scaling_factor <- 1
  i <- 2
  while(RoPI[i - 1] >= .5){
    signal_total[i] <- signal_total[i] * signal_scaling_factor * RoPI[i - 1]
    i <- i + 1
  }
  for(i in i:length(signal_total)){
    signal_total[i] <- signal_total[i - 1] * (((P[i]/P[i - 1] - 1) * RoPI[i]) + 1)
  }
  
  growth <- c(0, 0)
  for(j in 2:(length(signal_total) - 1)){
    growth[j + 1] <- signal_total[j + 1]/signal_total[j]
  }
  
  saturation_fraction <- signal_total/max(signal_total)
  
  P_signal_table <- data.frame(P, signal_total, growth, saturation_fraction)
  points(P_signal_table$P, P_signal_table$signal_total, pch = 19)
  return(P_signal_table)
}

# AF 647 - see if it works! Same MESF, QY = .33, lifetime = 1.0 ns, wavelength = .647, epsilon = 270,000
AF647_Sat_Prediction_Plotter <- function(){
  P <- seq(.001, .851, by = .05)
  cp <- suppressMessages(Contained_Power_Finder(.1, 40, 10, 7.65, 7.65))
  max_eff <- Signal_Power_Gain_Max_Set(270000, 270000, 1.0, .33, .1, cp, 40, 10, 7.65, .647, 4)
  
  while(max_eff > 270000){
    max_eff <- Signal_Power_Gain_Max_Set(270000, 270000, 1.0, .33, .1, cp, 40, 10, 7.65, .647, 4)
  }
  
  plot(1, xlim = c(0, .9), ylim = c(0, 1), type = 'n', xlab = 'Power (W)', ylab = 'Return on Power Input',
       main = 'AF647')
  RoPI_record <- c()
  for(i in P){
    efficiency <- Signal_Power_Gain(270000, 270000, 1.0, .33, i, cp, 40, 10, 7.65, .647, 4, max_eff)
    RoP <- efficiency
    RoPI_record <- c(RoPI_record, RoP)
    points(i, RoP, pch = 19)
  }
  results <- data.frame(P, RoPI_record)
  return(results)
}

AF647_Signal_Prediction_Plotter <- function(){
  P <- seq(.001, .851, by = .05)
  cp <- suppressMessages(Contained_Power_Finder(.1, 40, 10, 7.65, 7.65))
  signal <- c(0)
  signal_total <- c()
  RoPI <- AF647_Sat_Prediction_Plotter()$RoPI_record
  plot(1, xlim = c(0, .9), ylim = c(0, 270000), type = 'n', xlab = 'Power (W)', ylab = 'Signal',
       main = 'AF647')
  j <- 1
  for(i in P){
    signal <- c(signal,
                (Photon_Flux_State_Population(270000, 270000, 1.0, .33, i, cp, 40, 10, 7.65, .647, 4)[3] - 
                   signal[j])* RoPI[j])
    signal_total <- c(signal_total, sum(signal))
    j <- j + 1
  }
  
  # Since we don't have any data, we will use a signal scaling factor of 1
  # true_signal <- as.numeric(readline('Enter the fluorescent signal value at 50 mW.'))
  # signal_scaling_factor <- true_signal/signal_total[2]
  
  signal_scaling_factor <- 1
  i <- 2
  while(RoPI[i - 1] >= .5){
    signal_total[i] <- signal_total[i] * signal_scaling_factor * RoPI[i - 1]
    i <- i + 1
  }
  for(i in i:length(signal_total)){
    signal_total[i] <- signal_total[i - 1] * (((P[i]/P[i - 1] - 1) * RoPI[i]) + 1)
  }
  
  growth <- c(0, 0)
  for(j in 2:(length(signal_total) - 1)){
    growth[j + 1] <- signal_total[j + 1]/signal_total[j]
  }
  
  saturation_fraction <- signal_total/max(signal_total)
  
  P_signal_table <- data.frame(P, signal_total, growth, saturation_fraction)
  points(P_signal_table$P, P_signal_table$signal_total, pch = 19)
  return(P_signal_table)
}

