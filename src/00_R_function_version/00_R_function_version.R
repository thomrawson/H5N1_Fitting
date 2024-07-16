# This is an orderly script - edit it to suit your needs. You might include
#
# * orderly2::orderly_parameters():
#       declare parameters that your report accepts
# * orderly2::orderly_description():
#       describe your report with friendly names, descriptions and metadata
# * orderly2::orderly_resource():
#       declare files in your source tree that are inputs
# * orderly2::orderly_shared_resource():
#       use files from the root directory's 'shared/' directory
# * orderly2::orderly_dependency():
#       use files from a previously-run packet
# * orderly2::orderly_artefact():
#       declare files that you promise to produce, and describe them
# * orderly2::orderly_strict_mode():
#       enable some optional checks

## At what time do we start testing cattle that want to cross state borders?
orderly2::orderly_parameters(time_test = 136)  #Testing was effective Monday, April 29, 2024, if we assume starting Dec 15th, that's 136 days
orderly2::orderly_parameters(no_iterations = 100)
orderly2::orderly_parameters(beta = 1.6)
orderly2::orderly_parameters(gamma = 0.07)
orderly2::orderly_parameters(sigma = 0.14)
orderly2::orderly_parameters(alpha = 0.0002)
orderly2::orderly_parameters(num_cows_tested = 30)
orderly2::orderly_parameters(N_time_steps_days = 250)
orderly2::orderly_parameters(dt = 0.25)

# Write parameters to a .txt file
list(time_test = time_test, no_iterations = no_iterations, beta = beta, gamma = gamma,
     sigma = sigma, alpha = alpha, num_cows_tested = num_cows_tested,
     N_time_steps_days = N_time_steps_days, dt = dt) -> parameters
writeLines(
  paste(names(parameters), parameters, sep = " = "),
  con = "parameters.txt"
)
################################################################################################################
##This is the same as script 2, except that it runs it for multiple iterations and plots the 95% CI output
#INPUTS
#Initial conditions
#State_start_index - a vector of length N_states detailing the starting index for herds in state i
#State_end_index - a vector of length N_states detailing the ending index for herds in state i
#S - a vector containing the number of susceptibles for each herd. 
#E - a vector containing the number of susceptibles for each herd. 
#I - a vector containing the number of susceptibles for each herd. 
#R - a vector containing the number of susceptibles for each herd. 
#Model parameters
#beta - transmission within-herd
#alpha - transmission within-state
#gamma - recovery rate
#sigma - exposed to infected
#Herds_in_State[N_states] - How many herds are in each state?
#p_state_export[N_states] - The probability that a herd in a state exports a shipment of cattle
#p_cow_exported[N_states] - The probability that a cow in an exporting state is exported, this is calculated from mean shipment sizes
#State_to_State_shipment_probability[N_states, N_states] - A matrix where row i tells you the probabilities of which state a shipment from state i will go to. Hence, rows will sum to 1
#state_travel_allowed - a logical TRUE/FALSE, is travel between states for infected cattle allowed?
##################################################################################################################
library(tidyverse)

##Load the model function:
source("function.R")

##Load some pre-generated cow data
load("Data/Mean_Cattle_movement_data.RData")
load("Data/State_Movement_Data.RData")

##Wyoming and Rhode Island have such small dairy cow populations that we don't know how many herds they have.
## Let's just assume they both have 2 for now.
Herds_in_State <- State_Herd_Numbers$no_dairy_herds_2023
Herds_in_State[is.na(Herds_in_State)] <- 2

N_herd <- sum(Herds_in_State)
N_state <- length(Herds_in_State)

N_time_steps <- N_time_steps_days/dt

##Average herd size is # of cows divided by # of herds
State_Herd_Numbers$mean_herd_size <- (State_Herd_Numbers$milk_cows_Q4_2023_thousands/Herds_in_State)*1000

#Probability a state exports is its annual exports, divided by 365, divided by total number of herds
p_state_export <- rep(0,48)
for(i in 1:length(p_state_export)){
  p_state_export[i] <- sum(Mean_Cattle_export_shipments[i,])/(365*Herds_in_State[i])
}

#Probability an individual cow in the herd is exported. mean number shipped/mean size of herd.
p_cow_exported <- rep(0,48)
for(i in 1:length(p_cow_exported)){
  p_cow_exported[i] <- Mean_Cattle_mean_shipment_size[i,1]/State_Herd_Numbers$mean_herd_size[i]
}
##TODO: I need to go back and alter the probability for Rhode Island, no time, just set it same as Wyoming for now
p_cow_exported[is.na(p_cow_exported)] <- p_cow_exported[48]

## Scale each row of annual exports matrix to sum to 1. You don't actually have to do this for rmultinom in R, but good practice.
State_to_State_shipment_probability <- Mean_Cattle_export_shipments
for(i in 1:48){
  State_to_State_shipment_probability[i,] <- State_to_State_shipment_probability[i,]/sum(State_to_State_shipment_probability[i,])
}  

## Some placeholder model param values, eventually we'll fit to these
# beta <- 0.5
# gamma <- 1/10
# sigma <- 1/7
# alpha <- 0

## Initial Conditions
S_ini <- rep(0, N_herd)
E_ini <- rep(0, N_herd)
I_ini <- rep(0, N_herd)
R_ini <- rep(0, N_herd)

State_start_index <- rep(0, N_state)
State_end_index <- rep(0, N_state)

State_start_index[1] <- 1
State_end_index[1] <- Herds_in_State[1]
for(i in 2:length(State_start_index)){
  State_start_index[i] <- State_end_index[(i-1)]+1
  State_end_index[i] <- State_end_index[(i-1)]+Herds_in_State[i]
}

for(i in 1:N_state){
  index_start <- State_start_index[i]
  index_end <- State_end_index[i]
  index_seq <- c(index_start:index_end)
  
  S_ini[index_seq] <- ceiling(State_Herd_Numbers$mean_herd_size[i])
}

## And we're gonna start off with 5 cases in Texas
Texas_index <- which(All_states == "Texas")
I_ini[State_start_index[Texas_index]] <- 5
S_ini[State_start_index[Texas_index]] <- S_ini[State_start_index[Texas_index]] - I_ini[State_start_index[Texas_index]]


#N_time_steps <- 1000
S_tot <- array(0, dim = c(N_time_steps, no_iterations))
I_tot <- array(0, dim = c(N_time_steps, no_iterations))
R_tot <- array(0, dim = c(N_time_steps, no_iterations))

#Let's track state totals too
S_States <- array(0, dim = c(48, N_time_steps, no_iterations))
I_States <- array(0, dim = c(48, N_time_steps, no_iterations))
R_States <- array(0, dim = c(48, N_time_steps, no_iterations))


#################################################################################
##RUNNING THE MODEL
for(k in 1:no_iterations){
  
  if(k%%5 == 0){
    print(k)
  }
  
  track_time <- 0
  S_input <- S_ini
  E_input <- E_ini
  I_input <- I_ini
  R_input <- R_ini
  
for(i in 1:N_time_steps){
  
  ## At a certain time, we want to start testing cows for H5N1
  track_time <- track_time + dt
  if(track_time < time_test){
    state_travel_allowed <- TRUE
  }else{
    state_travel_allowed <- FALSE
  }
  
  ##Run the next time step
  output <- cattle_epi_model(dt = dt, #The size of the time step
                             S = S_input, E = E_input, I = I_input, R = R_input, #Initial conditions
                             State_start_index = State_start_index,
                             State_end_index = State_end_index,
                             Herds_in_State = Herds_in_State,
                             p_state_export = p_state_export*dt, 
                             p_cow_exported = p_cow_exported*dt,
                             num_cows_tested = num_cows_tested,
                             State_to_State_shipment_probability = State_to_State_shipment_probability,
                             state_travel_allowed = state_travel_allowed,
                             beta = beta, alpha = alpha, gamma = gamma, sigma = sigma #Model parameters
  )
  S_input <- output[[1]]
  E_input <- output[[2]]
  I_input <- output[[3]]
  R_input <- output[[4]]
  
  S_tot[i,k] <- sum(S_input)
  I_tot[i,k] <- sum(I_input)
  R_tot[i,k] <- sum(R_input)
  
  for(j in 1:48){
    index_start <- State_start_index[j]
    index_end <- State_end_index[j]
    index_seq <- c(index_start:index_end)
    S_States[j,i,k] <- sum(S_input[index_seq])
    I_States[j,i,k] <- sum(I_input[index_seq])
    R_States[j,i,k] <- sum(R_input[index_seq])
  }
  
}
}


# Initialize lists to store results
S_mean_values <- apply(S_States, c(1, 2), mean)
S_ci_values <- apply(S_States, c(1, 2), function(x) quantile(x, c(0.025, 0.975)))
I_mean_values <- apply(I_States, c(1, 2), mean)
I_ci_values <- apply(I_States, c(1, 2), function(x) quantile(x, c(0.025, 0.975)))
R_mean_values <- apply(R_States, c(1, 2), mean)
R_ci_values <- apply(R_States, c(1, 2), function(x) quantile(x, c(0.025, 0.975)))

S_tot_mean_values <- apply(S_tot, 1, mean)
S_tot_ci_values <- apply(S_tot, c(1), function(x) quantile(x, c(0.025, 0.975)))
I_tot_mean_values <- apply(I_tot, 1, mean)
I_tot_ci_values <- apply(I_tot, c(1), function(x) quantile(x, c(0.025, 0.975)))
R_tot_mean_values <- apply(R_tot, 1, mean)
R_tot_ci_values <- apply(R_tot, c(1), function(x) quantile(x, c(0.025, 0.975)))

no_data <- length(S_mean_values)
S_df <- data.frame(Time = rep(0,no_data),
                       US_State = rep("A", no_data),
                       Compartment = rep("S", no_data),
                       mean_value = rep(0,no_data),
                       lower_ci = rep(0,no_data),
                       upper_ci = rep(0,no_data))
I_df <- data.frame(Time = rep(0,no_data),
                   US_State = rep("A", no_data),
                   Compartment = rep("I", no_data),
                   mean_value = rep(0,no_data),
                   lower_ci = rep(0,no_data),
                   upper_ci = rep(0,no_data))
R_df <- data.frame(Time = rep(0,no_data),
                   US_State = rep("A", no_data),
                   Compartment = rep("R", no_data),
                   mean_value = rep(0,no_data),
                   lower_ci = rep(0,no_data),
                   upper_ci = rep(0,no_data))

index <- 0
for(i in 1:N_time_steps){
  for(j in 1:N_state){
   index <- index + 1
   S_df$Time[index] <- i*dt
   S_df$US_State[index] <- All_states[j]
   S_df$mean_value[index] <- S_mean_values[j,i]
   S_df$lower_ci[index] <- S_ci_values[1,j,i]
   S_df$upper_ci[index] <- S_ci_values[2,j,i]
   
   I_df$Time[index] <- i*dt
   I_df$US_State[index] <- All_states[j]
   I_df$mean_value[index] <- I_mean_values[j,i]
   I_df$lower_ci[index] <- I_ci_values[1,j,i]
   I_df$upper_ci[index] <- I_ci_values[2,j,i]
   
   R_df$Time[index] <- i*dt
   R_df$US_State[index] <- All_states[j]
   R_df$mean_value[index] <- R_mean_values[j,i]
   R_df$lower_ci[index] <- R_ci_values[1,j,i]
   R_df$upper_ci[index] <- R_ci_values[2,j,i]
     
  }
}



# Combine all into a master data frame
master_df <- bind_rows(S_df, I_df, R_df)
# Ensure the Compartment column is a factor with levels in the desired order
master_df$Compartment <- factor(master_df$Compartment, levels = c("S", "I", "R"))



S_tot_df <- data.frame(Time = seq(from = 0.25, by = 0.25, length.out = N_time_steps),
                   Compartment = rep("S", N_time_steps),
                   mean_value = S_tot_mean_values,
                   lower_ci = S_tot_ci_values[1,],
                   upper_ci = S_tot_ci_values[2,])
I_tot_df <- data.frame(Time = seq(from = 0.25, by = 0.25, length.out = N_time_steps),
                   Compartment = rep("I", N_time_steps),
                   mean_value = I_tot_mean_values,
                   lower_ci = I_tot_ci_values[1,],
                   upper_ci = I_tot_ci_values[2,])
R_tot_df <- data.frame(Time = seq(from = 0.25, by = 0.25, length.out = N_time_steps),
                   Compartment = rep("R", N_time_steps),
                   mean_value = R_tot_mean_values,
                   lower_ci = R_tot_ci_values[1,],
                   upper_ci = R_tot_ci_values[2,])
# Combine all into a master data frame
master_tot_df <- bind_rows(S_tot_df, I_tot_df, R_tot_df)
# Ensure the Compartment column is a factor with levels in the desired order
master_tot_df$Compartment <- factor(master_tot_df$Compartment, levels = c("S", "I", "R"))


dir.create("State_SIR_trajectories")

ggplot(master_tot_df) +
  geom_line(aes(x = Time, y = mean_value, color = Compartment)) +
  geom_ribbon(aes(x = Time, ymin = lower_ci, ymax = upper_ci,
                  fill = Compartment), alpha = 0.2) +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.3))) +
  scale_color_manual(values = c("S" = "#4B59B4", "I" = "#b44b59", "R" = "#59B44B")) +
  scale_fill_manual(values = c("S" = "#4B59B4", "I" = "#b44b59", "R" = "#59B44B")) +
  ggtitle("Total US Dairy Cattle") +
  theme_classic() -> tot_plot

ggsave(filename = "00_Total_US_Cases.png",
       path = 'State_SIR_trajectories', plot = tot_plot,
       dpi=300, height=6, width=6, units="in")



for(i in 1:48){
  state_hold <- All_states[i]
  data_hold <- filter(master_df, US_State == state_hold)
  ggplot(data_hold) +
    geom_line(aes(x = Time, y = mean_value, color = Compartment)) +
    geom_ribbon(aes(x = Time, ymin = lower_ci, ymax = upper_ci,
                    fill = Compartment), alpha = 0.2) +
    ggtitle(sprintf("%s", state_hold)) +
    theme(axis.text = element_text(size = rel(1.2)),
          axis.title = element_text(size = rel(1.3)),
          legend.text = element_text(size = rel(1.2)),
          legend.title = element_text(size = rel(1.3))) +
    scale_color_manual(values = c("S" = "#4B59B4", "I" = "#b44b59", "R" = "#59B44B")) +
    scale_fill_manual(values = c("S" = "#4B59B4", "I" = "#b44b59", "R" = "#59B44B")) +
    theme_classic() -> plot_hold
  
  ggsave(filename = sprintf("%s.png", state_hold),
         path = 'State_SIR_trajectories', plot = plot_hold,
         dpi=300, height=6, width=6, units="in")
  
}



