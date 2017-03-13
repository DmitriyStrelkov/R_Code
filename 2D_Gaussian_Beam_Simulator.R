# 2D Gaussian Beam Simulator
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


Real_Intensity_Matrix_Maker <- function(P, R_x, R_y){
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
  
  x_intensity <- Gaussian(P, R_x, x_values)
  y_intensity <- Gaussian(P, R_y, x_values)
  
  I_0 <- sqrt(max(x_intensity)^2 + max(y_intensity)^2)
  I_xyz <- matrix(0, ncol = length(x_values), nrow = length(y_values))
 
  if(R_x > R_y){
    for(i in 1:length(x_values)){
      for(j in 1:length(y_values)){
        I_xyz[i,j] <- I_0 * (exp((-2 * x_values[i]^2)/(.5 * R_x)^2) * exp((-2 * x_values[j]^2)/(.5 * R_y)^2))
      }
    }
  } else{
    for(i in 1:length(x_values)){
      for(j in 1:length(y_values)){
        I_xyz[i,j] <- I_0 * (exp((-2 * y_values[i]^2)/(.5 * R_x)^2) * exp((-2 * y_values[j]^2)/(.5 * R_y)^2))
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
      Radii[j,i] <- sqrt((x[i]^2)/a^2 + (y[j]^2)/b^2)
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
  
  if(R_x > R_y){
    Peak_Valley_Width_Ratio <- Gaussian(P, R_x, 0)/Gaussian(P, R_x, (.5 * x_dim))
  } else if(R_y > R_x){
    Peak_Valley_Width_Ratio <- Gaussian(P, R_y, 0)/Gaussian(P, R_y, (.5 * y_dim))
  }
  long_axis_result <- paste('The peak to valley ratio of the long axis at the short axis maximum is', 
                            Peak_Valley_Width_Ratio, sep = ' ')
  print(long_axis_result)
  return(Peak_Valley_Width_Ratio)
}

# This function shows the contained power of a given beam in a window.

Contained_Power_Finder <- function(P, R_x, R_y, x_dim, y_dim){
  Real_Intensity_Matrix <- Real_Intensity_Matrix_Maker(P, R_x, R_y)
  
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
    Radii <- Radial_Matrix_Maker(x_values, x_values, (R_x/2), (R_y/2))
  } else if(R_y <= R_x){
    Radii <- Radial_Matrix_Maker(y_values, y_values, (R_x/2), (R_y/2))
  }
  Window_Area_Matrix <- (Radii <= 1) * 1
  Contained_Intensity_Matrix <- Real_Intensity_Matrix * Window_Area_Matrix
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
  contained_s_power_density <- contained_power_value/((.5 * x_dim * 10^-4) * (.5 * y_dim * 10^-4) * pi) * 2
  print(paste('The fraction of total power contained in the window is', contained_power_fraction, sep = ' '))
  print(paste('The total power contained in the window is', contained_power_value, 'Watts', sep = ' '))
  print(paste('The max surface power density in the window at the focus is', contained_s_power_density, 'W/cm^2',
              sep = ' '))
  return(contained_power_fraction)
}

# This function finds the ratio between the power of the beam in the center of the flow channel and at the back side
# of the flow channel. Inputs are wavelength, power, beam dimensions at focus, back to center distance of the flow channel,
# and window dimensions. Output is the ratio between the power in the window at the center of the  channel and the 
# power in the window at the back of the channel. All lengths should be in microns.

# This function returns the coefficient of variation.
CV_Finder <- function(Intensity_Matrix){
  CV <- sd(Intensity_Matrix)/mean(Intensity_Matrix)
  result <- paste("The coefficient of variation is", CV, sep = ' ')
  return(result)
}

# This function returns the changing intensity profile as the beam propagates along the optic axis. The intensity
# is shown along the x and z dimensions; the y dimension is assumed to be at the intensity maximum (y = 0).

Propagation_Intensity_Plotter <- function(P, wavelength, R_x, R_y, z, window){
  x_values <- seq(-R_x, R_x, by = R_x/50)
  z_values <- seq(-z, z, by = z/50)
  w_z <- 0
  
  # Add zeros to values.
  if(!(is.element(0, x_values))){
    x_zero_index <- which(x_values == min(x_values[x_values > 0])) - 1
    
    x_values <- append(x_values, 0, x_zero_index)
    z_values <- append(z_values, max(z_values) + z/50, length(z_values))
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
  window_values <- seq(-(.5 * window), (.5 * window), by = window/100)
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
