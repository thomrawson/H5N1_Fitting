##Build a simple-to-parse function explaining each forward time-step of my model.
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


cattle_epi_model <- function(dt = 0.25, #The size of the time step
                             S, E, I, R, #Initial conditions
                             State_start_index,
                             State_end_index,
                             Herds_in_State,
                             p_state_export, p_cow_exported,
                             num_cows_tested,
                             State_to_State_shipment_probability,
                             state_travel_allowed,
                             beta, alpha, gamma, sigma #Model parameters
){
  
  
  
  N_herds <- length(S)
  N_states <- length(State_start_index)
  
  #Calculate the total herd sizes, N:
  N <- rep(0,N_herds)
  
  N <- S + E + I + R
  
  
  #First, we deal with the progression between disease compartments. Exports come after.
  #Calculate the force of infection, lambda, for each herd in each state:
  lambda <- rep(0, N_herds)
  p_SE <- rep(0, N_herds)
  
  ## Other probabilities of transition:
  p_EI <- 1 - exp(-sigma * dt) # E to I
  p_IR <- 1 - exp(-gamma * dt) # I to R
  
  ## Calculate disease progressions:
  n_SE <- rep(0, N_herds)
  n_EI <- rep(0, N_herds)
  n_IR <- rep(0, N_herds)
  
  new_S <- rep(0, N_herds)
  new_E <- rep(0, N_herds)
  new_I <- rep(0, N_herds)
  new_R <- rep(0, N_herds)
  
  for(i in 1:N_states){
    start_index <- State_start_index[i]
    end_index <- State_end_index[i]
    index_seq <- c(start_index:end_index)
    
    lambda[index_seq] <- ifelse(
      N[index_seq] == 0,
      0,
      beta *((I[index_seq]/N[index_seq]) + (alpha)*((sum(I[index_seq]) - I[index_seq])/(sum(N[index_seq]) - N[index_seq])) )
    )
    
    #Calculate the probability of transition too:
    p_SE[index_seq] <- 1 - exp(-lambda[index_seq] * dt) #S to E
    
    n_SE[index_seq] <- rbinom(length(index_seq),S[index_seq], p_SE[index_seq])
    
  }
  
  n_EI <- rbinom(length(n_EI),E, p_EI)
  n_IR <- rbinom(length(n_IR),I, p_IR)
  
  
  
  #And change the populations (we do this BEFORE calculating import/exports)
  new_S <- S - n_SE
  new_E <- E + n_SE - n_EI
  new_I <- I + n_EI - n_IR
  new_R <- R + n_IR
  
  
  ## Now calculate the exports from each herd
  Exports_S <- rep(0, N_herds)
  Exports_E <- rep(0, N_herds)
  Exports_I <- rep(0, N_herds)
  Exports_R <- rep(0, N_herds)
  Exports_N <- rep(0, N_herds)
  
  for(i in 1:N_states){
    start_index <- State_start_index[i]
    end_index <- State_end_index[i]
    index_seq <- c(start_index:end_index)
    Export_switch <- rbinom(length(index_seq),1,p_state_export[i])
    
    # <- (Will it export?)*(How many exported?)
    Exports_S[index_seq] <- Export_switch*rbinom(length(index_seq),new_S[index_seq], p_cow_exported[i])
    Exports_E[index_seq] <- Export_switch*rbinom(length(index_seq),new_E[index_seq], p_cow_exported[i])
    Exports_I[index_seq] <- Export_switch*rbinom(length(index_seq),new_I[index_seq], p_cow_exported[i])
    Exports_R[index_seq] <- Export_switch*rbinom(length(index_seq),new_R[index_seq], p_cow_exported[i])
    Exports_N[index_seq] <- Exports_S[index_seq] + Exports_E[index_seq] + Exports_I[index_seq] + Exports_R[index_seq]
  }
  
  
  ##Now we calculate where that export is going to.
  Imports_S <- rep(0, N_herds)
  Imports_E <- rep(0, N_herds)
  Imports_I <- rep(0, N_herds)
  Imports_R <- rep(0, N_herds)
  
  for(i in 1:N_states){
    start_index <- State_start_index[i]
    end_index <- State_end_index[i]
    index_seq <- c(start_index:end_index)
    
    for(j in index_seq){
      if(Exports_N[j] == 0){} else{
        #Can't really vectorise out this loop without removing the ability for two herds to export to the same herd.
        #Could do this if the speed up is desired, but don't think it'll work.
        Destination_State <- which(rmultinom(1,1,State_to_State_shipment_probability[i,]) == 1)
        Destination_Herd <- State_start_index[Destination_State] + ceiling(runif(1, 0, Herds_in_State[Destination_State])) - 1
        
        
        ## Eventually, the US implemented testing up to 30 cows in a shipment out-of-state for H5N1. Any positives, and the shipment couldn't happen.
        if(state_travel_allowed){
          Imports_S[Destination_Herd] <- Imports_S[Destination_Herd] + Exports_S[j]
          Imports_E[Destination_Herd] <- Imports_E[Destination_Herd] + Exports_E[j]
          Imports_I[Destination_Herd] <- Imports_I[Destination_Herd] + Exports_I[j]
          Imports_R[Destination_Herd] <- Imports_R[Destination_Herd] + Exports_R[j]
        } else{
          Imports_S_hold <-  Exports_S[j]
          Imports_E_hold <-  Exports_E[j]
          Imports_I_hold <-  Exports_I[j]
          Imports_R_hold <-  Exports_R[j]
          Imports_N_hold <- Imports_S_hold + Imports_E_hold + Imports_I_hold + Imports_R_hold
          
          if(Destination_State != i){
            infection_detected <- rhyper(1,m = Imports_I_hold, n = Imports_N_hold, k = min(Imports_N_hold,num_cows_tested)) ##Max 30 sampled
            if(infection_detected == 0){
              Imports_S[Destination_Herd] <- Imports_S[Destination_Herd] + Imports_S_hold
              Imports_E[Destination_Herd] <- Imports_E[Destination_Herd] + Imports_E_hold
              Imports_I[Destination_Herd] <- Imports_I[Destination_Herd] + Imports_I_hold
              Imports_R[Destination_Herd] <- Imports_R[Destination_Herd] + Imports_R_hold
            }
            
          }else{
            ## If it's staying within-state
            Imports_S[Destination_Herd] <- Imports_S[Destination_Herd] + Imports_S_hold
            Imports_E[Destination_Herd] <- Imports_E[Destination_Herd] + Imports_E_hold
            Imports_I[Destination_Herd] <- Imports_I[Destination_Herd] + Imports_I_hold
            Imports_R[Destination_Herd] <- Imports_R[Destination_Herd] + Imports_R_hold
          }
          
        }
        
      }
      
    }
  }
  
  ##Then add/remove all the imports/exports
  new_S <- new_S - Exports_S + Imports_S
  new_E <- new_E - Exports_E + Imports_E
  new_I <- new_I - Exports_I + Imports_I
  new_R <- new_R - Exports_R + Imports_R
  
  output_list <- list(new_S, new_E, new_I, new_R)
  return(output_list)
}