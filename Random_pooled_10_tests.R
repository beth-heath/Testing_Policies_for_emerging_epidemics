#Random pooled function for 10 tests per day

## extracted from hitchings
infect_contacts <- function(potential_contacts,beta_value){
  num_neighbours_susc <- length(potential_contacts)
  # Sample from each group of neighbours in turn
  # First choose how many neighbours each node infects
  num_contacts_inf <- rbinom(1,num_neighbours_susc,1-exp(-beta_value))
  # Then sample from the neighbours
  # If one node gets picked twice by different nodes, just discard the duplicate.
  infectees_susc <- c()
  if(num_contacts_inf>0){
    sample_indices <- ceiling((runif(1) + 1:num_contacts_inf-1)*num_neighbours_susc/num_contacts_inf)
    infectees_temp <- potential_contacts[sample_indices] 
    infectees_susc <- funique(infectees_temp)
  }
  return(infectees_susc)
}

infect_from_source <- function( s_nodes, v_nodes, e_nodes_info, direct_VE,incperiod_shape, incperiod_rate,rate_from_source){
  # identify susceptibles
  infectees_susc <- c()
  s_hr_l <- s_nodes==1
  ##!! not trial ppts for now
  s_hr <- g_name[s_hr_l & v_nodes==0 & t_nodes==0]
  # infect
  if(length(s_hr)>0)
    infectees_susc <- funique(infect_contacts(s_hr,beta_value=rate_from_source))
  newinfected <- infectees_susc
  # infect vaccinated
  if(sum(v_nodes)>0){
    infectees_susc <- c()
    beta_v <- rate_from_source*(1-direct_VE)
    s_hr <- g_name[s_hr_l & v_nodes==1 & t_nodes==0]
    if(length(s_hr)>0)
      infectees_susc <- funique(infect_contacts(s_hr,beta_value=beta_v))
    newinfected <- c(newinfected,infectees_susc)
  }
  
  if (length(newinfected)>0) {
    # Give each newly exposed node an incubation/latent period
    inc_periods <- incperiod_const + rgamma(length(newinfected),shape=incperiod_shape,rate=incperiod_rate)
    # Add them to e_nodes and remove from s_nodes and v_nodes
    for(i in 1:length(newinfected))
      e_nodes_info <- rbind(e_nodes_info,c(newinfected[i],0,inc_periods[i]))
  }
  return(e_nodes_info)
}

spread <- function( s_nodes, v_nodes, e_nodes_info, current_infectious, direct_VE,incperiod_shape, incperiod_rate,susc_list=contact_list,beta_scalar=1){
  # Spread will create new infected nodes from two sources: infectious nodes within the the study
  # population, and external pressure from the source population
  # Inputs:
  # s_nodes, e_nodes and i_nodes are susceptible, exposed and infected nodes
  # beta_base is the hazard of infection for one contact
  # incperiod_shape and rate are used to assign each newly exposed node a latent/incubation period
  # length, currently drawn from a gamma distribution
  scaled_beta <- beta_scalar * beta_base
  # Process: go through list of i_nodes, and choose a random number of its susceptible neighbours to
  # be infected, according to beta_base and choose a random number of its susceptible vaccinated neighbours to
  # be infected, according to beta_base and direct_VE
  infectees_susc <- c()
  # Get a list of all neighbours of all infected nodes
  potential_contacts <- c()
  #potential_hr_contacts <- c()
  #potential_neighbours <- c()
  for(i in current_infectious) {
    #potential_hr_contacts <- c(potential_hr_contacts,high_risk_list[[i]],household_list[[i]])
    #potential_neighbours <- c(potential_neighbours,contact_of_contact_list[[i]])
    potential_contacts <- c(potential_contacts,susc_list[[i]])
  }
  # remove duplications
  #potential_contacts <- potential_contacts[!potential_contacts%in%potential_hr_contacts]
  # infect high risk
  #s_hr_l <- s_nodes[potential_hr_contacts]==1
  #v_hr_l <- v_nodes[potential_hr_contacts]
  #s_hr <- potential_hr_contacts[s_hr_l & v_hr_l==0]
  #if(length(s_hr)>0)
  #  infectees_hr_susc <- funique(infect_contacts(s_hr,beta_value=scaled_beta*high_risk_scalar))
  # infect neighbours
  #s_nb_l <-  s_nodes[potential_neighbours]==1
  #v_nb_l <-  v_nodes[potential_neighbours]
  #s_nb <- potential_neighbours[s_nb_l & v_nb_l==0]
  #if(length(s_nb)>0)
  #  infectees_n_susc <- funique(infect_contacts(s_nb,beta_value=scaled_beta*neighbour_scalar))
  # infect other contacts
  if(length(potential_contacts)>0){
    s_contacts_l <- s_nodes[potential_contacts]==1
    v_contacts_l <- v_nodes[potential_contacts]
    s_contacts <- potential_contacts[s_contacts_l & v_contacts_l==0]
    if(length(s_contacts)>0)
      infectees_susc <- funique(infect_contacts(s_contacts,beta_value=scaled_beta))
  }
  newinfected <- funique(infectees_susc) #funique(c(infectees_susc,infectees_hr_susc,infectees_n_susc))
  # infect vaccinated
  if(sum(v_nodes)>0){
    infectees_susc <- c()
    beta_v <- scaled_beta*(1-direct_VE)
    # infect high risk
    #s_hr <- potential_hr_contacts[s_hr_l & v_hr_l==1]
    #if(length(s_hr)>0)
    #  infectees_hr_susc <- funique(infect_contacts(s_hr,beta_value=beta_v*high_risk_scalar))
    # infect neighbours
    #s_nb <- potential_neighbours[s_nb_l & v_nb_l==1]
    #if(length(s_nb)>0)
    #  infectees_n_susc <- funique(infect_contacts(s_nb,beta_value=beta_v*neighbour_scalar))
    # infect other contacts
    if(length(potential_contacts)>0){
      s_contacts <- potential_contacts[s_contacts_l & v_contacts_l==1]
      if(length(s_contacts)>0)
        infectees_susc <- funique(infect_contacts(s_contacts,beta_value=beta_v))
    }
    newinfected <- c(newinfected,funique(infectees_susc))
  }
  
  if (length(newinfected)>0) {
    # Give each newly exposed node an incubation/latent period
    inc_periods <- incperiod_const + rgamma(length(newinfected),shape=incperiod_shape,rate=incperiod_rate)
    # Add them to e_nodes and remove from s_nodes and v_nodes
    for(i in 1:length(newinfected))
      e_nodes_info <- rbind(e_nodes_info,c(newinfected[i],0,inc_periods[i]))
  }
  return(e_nodes_info)
}

recover <- function(e_nodes_info,i_nodes_info, infperiod_shape,infperiod_rate,cluster_people_index,time_diff=NULL) {
  # Input is a list of the exposed nodes, 
  # with number of days since infection and total incubation/latent
  # period, and equivalently for the infectious nodes.
  # For each of these nodes, we will add it to newinfectious if the number of days exposed has
  # reached the total length of the incubation period, and equivalently for the infectious nodes.
  
  # Advance infectious nodes first, otherwise you will doubly advance any nodes switching from
  # exposed to infectious at this time step
  indices_to_remove <- i_nodes_info[,2]>=i_nodes_info[,3]
  newremoved <- as.vector(i_nodes_info[indices_to_remove,1])
  
  # Add one day to the length of each infected individual's time infected
  i_nodes_info[,2] <- i_nodes_info[,2]+1
  
  #Have an index for the length of time an individual has spent in isolation
  #q_nodes_info[,2] <- q_nodes_info[,2]+1
  
  # Remove any recovered from i_nodes and add to r_nodes
  i_nodes_info <- i_nodes_info[!indices_to_remove,,drop=FALSE]
  
  # Now advance exposed nodes
  indices_to_remove <- e_nodes_info[,2]>=e_nodes_info[,3]
  newinfectious <- as.vector(e_nodes_info[indices_to_remove,1])
  incubation_days <- as.vector(e_nodes_info[indices_to_remove,2])+1
  
  # Add one day to the length of each infected individual's time infected
  e_nodes_info[,2] <- e_nodes_info[,2]+1
  
  # Remove any progressing from e_nodes_info and add to i_nodes
  e_nodes_info <- e_nodes_info[!indices_to_remove,,drop=FALSE]
  if(length(newinfectious)>0){
    inf_periods <- infperiod_const + rgamma(length(newinfectious),infperiod_shape,rate=infperiod_rate)
    #hosp_time <- rgamma(length(newinfectious),shape=hosp_shape,rate=hosp_rate)
    if(!is.null(time_diff)){
      hosp_time <- c()
      for(i in 1:length(newinfectious))
        if(cluster_people_index[newinfectious[i]]==1){
          hosp_time[i] <- rtruncnorm(1,a=0,mean=hosp_mean,sd=hosp_sd)
        }else{
          hosp_time[i] <- rtruncnorm(1,a=0,mean=hosp_mean_index,sd=hosp_sd_index)
        }
      # person is infectious until the minimal time that they are hospitalised or removed otherwise
      for(i in 1:length(inf_periods)) if(hosp_time[i] < inf_periods[i]) inf_periods[i] <- hosp_time[i] 
      #min_time <- pmin(inf_periods,hosp_time)
      # if person is enrolled soon, they will be hospitalised then
      if(0<=time_diff) for(i in 1:length(inf_periods)) if(time_diff < inf_periods[i]) inf_periods[i] <- time_diff 
    }
    #i_nodes <- rbind(i_nodes,cbind(newinfectious,rep(0,length(newinfectious)),min_time,incubation_days))
    for(i in 1:length(newinfectious)) 
      i_nodes_info <- rbind(i_nodes_info,c(newinfectious[i],0,inf_periods[i],incubation_days[i],runif(1)<observed))
  }
  infected_persons <- i_nodes_info[,1][i_nodes_info[,2] > 0 ]
  # testing_people <- sample(g_name, 100, replace=F)
  # names_of <- infected_persons[infected_persons %in% testing_people]
  list(e_nodes_info, i_nodes_info, newremoved, newinfectious, infected_persons)
}


simulate_contact_network <- function(first_infected,individual_recruitment_times=F,end_time=31,start_day=0,from_source=0,cluster_flag=0,allocation_ratio=0.5,
                                     direct_VE=0,base_rate=0,spread_wrapper=covid_spread_wrapper, spread_wrapper_2 = covid_spread_wrapper_2, 
                                     prob_false_neg, pool_size,non_compliance_prob, tests, pools_every_time){
  # set up info to store
  `%ni%` <- Negate(`%in%`)
  non_compliers <- sample(g_name, length(g_name)*non_compliance_prob, replace = T )
  trajectories <- list()
  over_point<-0
  over_point_1<-0
  trajectories$S <- length(vertices) - 1
  trajectories$E <- 0
  trajectories$I <- 0
  trajectories$R <- 0
  e_nodes_info <- matrix(nrow=0,ncol=3)
  i_nodes_info <- matrix(nrow=0,ncol=5)
  e_nodes <- rep(0,length(g_name))
  i_nodes <- rep(0,length(g_name))
  v_nodes <- rep(0,length(g_name))
  c_nodes <- rep(0,length(g_name))
  s_nodes <- rep(1,length(g_name))
  r_nodes <- rep(0,length(g_name))
  t_nodes <- rep(0,length(g_name))
  z_nodes <- rep(0, length(g_name))
  q_nodes <- rep(0,length(g_name))
  p_nodes <- rep(0, length(g_name))
  a_nodes<-rep(0,length(g_name))
  q_nodes_info <- rep(0, length(g_name))
  w_nodes_info <-rep(0, length(g_name))
  
  
  # generate info for index case
  inc_time <- incperiod_const + rgamma(length(first_infected),shape=incperiod_shape,rate=incperiod_rate)
  ceil_inc_time <- ceiling(inc_time)
  #i_nodes_info <- rbind(i_nodes_info,c(first_infected,rep(0,length(first_infected)),inf_time,inc_time))
  for (i in 1:length(first_infected)) 
    e_nodes_info <- rbind(e_nodes_info,c(first_infected[i],0,inc_time[i]))
  
  s_nodes[first_infected] <- 0
  e_nodes[first_infected] <- 1
  
  #recruitment_time <- round(rgamma(1,shape=recruit_shape,rate=recruit_rate))
  recruitment_time <- ceiling(rtruncnorm(1,a=0,mean=recruit_mean,sd=recruit_sd))
  results <- matrix(nrow=0,ncol=5)#c(first_infected,0,-inc_time,NA),nrow=1)
  numinfectious <- 1

  
  # identify contacts of index case
  contacts_1 <- contact_list[first_infected[1]]
  contacts_2 <-contact_list[first_infected[2]]
  contacts_3 <-contact_list[first_infected[3]]
  contacts_4 <-contact_list[first_infected[4]]
  contacts_5 <-contact_list[first_infected[5]]
  contacts_6 <-contact_list[first_infected[6]]
  contacts_7 <-contact_list[first_infected[7]]
  contacts_8 <-contact_list[first_infected[8]]
  contacts_9 <-contact_list[first_infected[9]]
  contacts_10 <-contact_list[first_infected[10]]
  contacts_11 <- contact_list[first_infected[11]]
  contacts_12 <-contact_list[first_infected[12]]
  contacts_13 <-contact_list[first_infected[13]]
  contacts_14 <-contact_list[first_infected[14]]
  contacts_15 <-contact_list[first_infected[15]]
  contacts_16 <-contact_list[first_infected[16]]
  contacts_17 <-contact_list[first_infected[17]]
  contacts_18 <-contact_list[first_infected[18]]
  contacts_19 <-contact_list[first_infected[19]]
  contacts_20 <-contact_list[first_infected[20]]

  
  
  contacts <-funique(c(contacts_1, contacts_2, contacts_3, contacts_4, contacts_5, contacts_6, contacts_7, contacts_8, contacts_9, contacts_10, contacts_11, contacts_12, contacts_13, contacts_14, contacts_15, contacts_16, contacts_17, contacts_18, contacts_19, contacts_20))

  order_infected <- first_infected
  ## identify high-risk people
  ##!! all household members are high risk. 
  #if(!exists('high_risk_list')) high_risk_list <- lapply(g_name,function(x)c())
  #high_risk <- high_risk_list[[first_infected]]
  ## contacts of contacts
  # identify contacts of index case
  contacts_of_contacts_1 <- contact_of_contact_list[first_infected[1]]
  contacts_of_contacts_2 <- contact_of_contact_list[first_infected[2]]
  contacts_of_contacts_3 <- contact_of_contact_list[first_infected[3]]
  contacts_of_contacts_4 <- contact_of_contact_list[first_infected[4]]
  contacts_of_contacts_5 <- contact_of_contact_list[first_infected[5]]
  contacts_of_contacts_6 <- contact_of_contact_list[first_infected[6]]
  contacts_of_contacts_7 <- contact_of_contact_list[first_infected[7]]
  contacts_of_contacts_8 <- contact_of_contact_list[first_infected[8]]
  contacts_of_contacts_9 <- contact_of_contact_list[first_infected[9]]
  contacts_of_contacts_10 <- contact_of_contact_list[first_infected[10]]
  contacts_of_contacts_11 <- contact_of_contact_list[first_infected[11]]
  contacts_of_contacts_12 <- contact_of_contact_list[first_infected[12]]
  contacts_of_contacts_13 <- contact_of_contact_list[first_infected[13]]
  contacts_of_contacts_14 <- contact_of_contact_list[first_infected[14]]
  contacts_of_contacts_15 <- contact_of_contact_list[first_infected[15]]
  contacts_of_contacts_16 <- contact_of_contact_list[first_infected[16]]
  contacts_of_contacts_17 <- contact_of_contact_list[first_infected[17]]
  contacts_of_contacts_18 <- contact_of_contact_list[first_infected[18]]
  contacts_of_contacts_19 <- contact_of_contact_list[first_infected[19]]
  contacts_of_contacts_20 <- contact_of_contact_list[first_infected[20]]
  
  
  ## add households of high-risk contacts to contacts of contacts
  #if(length(high_risk)>0) 
  # for(hr in high_risk)
  #  contacts_of_contacts <- c(contacts_of_contacts,household_list[[hr]])
  #high_risk <- c(high_risk,household_list[[first_infected]])

  contacts_of_contacts <-funique(c(contacts_of_contacts_1, contacts_of_contacts_2, contacts_of_contacts_3, contacts_of_contacts_4, contacts_of_contacts_5, contacts_of_contacts_6, contacts_of_contacts_7, contacts_of_contacts_8, contacts_of_contacts_9, contacts_of_contacts_10, contacts_of_contacts_11, contacts_of_contacts_12, contacts_of_contacts_13, contacts_of_contacts_14, contacts_of_contacts_15, contacts_of_contacts_16, contacts_of_contacts_17, contacts_of_contacts_18, contacts_of_contacts_19, contacts_of_contacts_20))

  
  cluster_people <- unlist(funique(c(contacts,contacts_of_contacts)))
  cluster_people_index <- g_name%in%cluster_people
  
  
  
  # enroll trial participants
  n_trial_participants <- rbinom(1,length(cluster_people),enrollment_rate)
  trial_participants <- sample(cluster_people,n_trial_participants,replace=T)
  t_nodes[trial_participants] <- 1
  vaccinees <- c()
  t_nodes <<- t_nodes
  if(cluster_flag==0){
    nvacc <- round(length(trial_participants)*allocation_ratio)
    vaccinees <- trial_participants[sort(sample.int(length(trial_participants),nvacc,replace=T))]
  }else{
    if(runif(1)<allocation_ratio)
      vaccinees <- trial_participants
  }
  
  vaccine_incubation_times <- 0
  if(length(vaccinees)>0)
    vaccine_incubation_times <- rgamma(length(vaccinees),shape=vacc_shape,rate=vacc_rate)
  if(individual_recruitment_times==F){
    recruitment_times <- rep(recruitment_time,n_trial_participants) + ceil_inc_time
  }else{
    recruitment_times <- sample(1:recruitment_time,n_trial_participants,replace=T) + ceil_inc_time
  }
  # roll epidemic forward one day at a time
  #sim_time <- recruitment_time + end_time + ceil_inc_time
  sim_time <- 500
  report <- rep(0,500)
  
  positive = matrix(nrow = sim_time, ncol = length(g_name))
  
  for(time_step in 1:sim_time){
    ## vaccination given time to develop immunity
    if(length(vaccinees)>0) {
      developed <- vaccine_incubation_times<=time_step-recruitment_times[trial_participants%in%vaccinees]
      v_nodes[vaccinees[developed]] <- 1
    }
    
    
    
    
    # update everyone's internal clock
    newinfectious <- newremoved <- c()
    if ((nrow(e_nodes_info)>0)||(nrow(i_nodes_info)>0)) {
      time_diff <- NULL
      if(!individual_recruitment_times) time_diff <- recruitment_time-time_step
      rec_list <- recover(e_nodes_info,i_nodes_info,infperiod_shape,infperiod_rate,cluster_people_index=cluster_people_index,time_diff=time_diff)
      e_nodes_info <- rec_list[[1]]
      i_nodes_info <- rec_list[[2]]
      newremoved <- rec_list[[3]]
      i_nodes[newremoved] <- 0
      r_nodes[newremoved] <- 1 
      newinfectious <- rec_list[[4]]
      e_nodes[newinfectious] <- 0
      i_nodes[newinfectious] <- 1
      infected_persons <-g_name[i_nodes==1]
      exposed <- g_name[e_nodes ==1]
      infected <- g_name[i_nodes==1]
      all <- c(exposed, infected)

      
      total_inf <- length(infected_persons)
      infection_period <- 0
      report[time_step] <- total_inf
      if (total_inf ==! 0){
        infection_period = infection_period + 1
      } else {
        infection_period = infection_period + 1
      }

      non_isol <- g_name[q_nodes == 0]
      non_isolating <- non_isol[non_isol %ni% non_compliers]
      #organising the non-isolated population into pools
      if (length(non_isolating) > pool_size){
        pool_1 <- sample(non_isolating, pool_size, replace=T)
        g_name_1 <- non_isolating[non_isolating %ni% pool_1]
      }
      if (length(g_name_1) > pool_size){
        pool_2 <- sample(g_name_1, pool_size, replace=T)
        g_name_2 <- g_name_1[g_name_1 %ni% pool_2]
      }
      if (length(g_name_2) > pool_size){
        pool_3 <- sample(g_name_2, pool_size, replace=T)
        g_name_3 <- g_name_2[g_name_2 %ni% pool_3]
      }
      if (length(g_name_3) > pool_size){
        pool_4 <- sample(g_name_3, pool_size, replace=T)
        g_name_4 <- g_name_3[g_name_3 %ni% pool_4]
      }
      if (length(g_name_4) > pool_size){
        pool_5 <- sample(g_name_4, pool_size, replace=T)
        g_name_5 <- g_name_4[g_name_4 %ni% pool_5]
      }
      if (length(g_name_5) > pool_size){
        pool_6 <- sample(g_name_5, pool_size, replace=T)
        g_name_6 <- g_name_5[g_name_5 %ni% pool_6]
      }
      if (length(g_name_6) > pool_size){
        pool_7 <- sample(g_name_6, pool_size, replace=T)
        g_name_7 <- g_name_6[g_name_6 %ni% pool_7]
      }
      if (length(g_name_7) > pool_size){
        pool_8 <- sample(g_name_7, pool_size, replace=T)
        g_name_8 <- g_name_7[g_name_7 %ni% pool_8]
      }
      if (length(g_name_8) > pool_size){
        pool_9 <- sample(g_name_8, pool_size, replace=T)
        g_name_9 <- g_name_8[g_name_8 %ni% pool_9]
      }
      if (length(g_name_9) > pool_size){
        pool_10 <- sample(g_name_9, pool_size, replace=T)
        g_name_10 <- g_name_9[g_name_9 %ni% pool_10]
      }
      if (length(g_name_9) > pool_size){
        pool_10 <- sample(g_name_9, pool_size, replace=T)
        g_name_10 <- g_name_9[g_name_9 %ni% pool_10]
      }
      if (length(g_name_10) > pool_size){
        pool_11 <- sample(g_name_10, pool_size, replace=T)
        g_name_11 <- g_name_10[g_name_10 %ni% pool_11]
      }
      if (length(g_name_11) > pool_size){
        pool_12 <- sample(g_name_11, pool_size, replace=T)
        g_name_12 <- g_name_11[g_name_11 %ni% pool_12]
      }
      if (length(g_name_12) > pool_size){
        pool_13 <- sample(g_name_12, pool_size, replace=T)
        g_name_13 <- g_name_12[g_name_12 %ni% pool_13]
      }
      if (length(g_name_13) > pool_size){
        pool_14 <- sample(g_name_13, pool_size, replace=T)
        g_name_14 <- g_name_13[g_name_13 %ni% pool_14]
      }
      if (length(g_name_14) > pool_size){
        pool_15 <- sample(g_name_14, pool_size, replace=T)
        g_name_15 <- g_name_14[g_name_14 %ni% pool_15]
      }
      if (length(g_name_15) > pool_size){
        pool_16  <- sample(g_name_15, pool_size, replace=T)
        g_name_16 <- g_name_15[g_name_15 %ni% pool_16]
      }
      if (length(g_name_16) > pool_size){
        pool_17 <- sample(g_name_16, pool_size, replace=T)
        g_name_17 <- g_name_16[g_name_16 %ni% pool_17]
      }
      if (length(g_name_17) > pool_size){
        pool_18 <- sample(g_name_17, pool_size, replace=T)
        g_name_18 <- g_name_17[g_name_17 %ni% pool_18]
      }
      if (length(g_name_18) > pool_size){
        pool_19 <- sample(g_name_18, pool_size, replace=T)
        g_name_19 <- g_name_18[g_name_18 %ni% pool_19]
      }
      if (length(g_name_19) > pool_size){
        pool_20 <- sample(g_name_19, pool_size, replace=T)
        g_name_20 <- g_name_19[g_name_19 %ni% pool_20]
      }
      if (length(g_name_20) > pool_size){
        pool_21 <- sample(g_name_20, pool_size, replace=T)
        g_name_21 <- g_name_20[g_name_20 %ni% pool_21]
      }
      if (length(g_name_21) > pool_size){
        pool_22 <- sample(g_name_21, pool_size, replace=T)
        g_name_22 <- g_name_21[g_name_21 %ni% pool_22]
      }
      if (length(g_name_22) > pool_size){
        pool_23 <- sample(g_name_22, pool_size, replace=T)
        g_name_23 <- g_name_22[g_name_22 %ni% pool_23]
      }
      if (length(g_name_23) > pool_size){
        pool_24 <- sample(g_name_23, pool_size, replace=T)
        g_name_24 <- g_name_23[g_name_23 %ni% pool_24]
      }
      if (length(g_name_24) > pool_size){
        pool_25 <- sample(g_name_24, pool_size, replace=T)
        g_name_25 <- g_name_24[g_name_24 %ni% pool_25]
      }
      if (length(g_name_25) > pool_size){
        pool_26 <- sample(g_name_25, pool_size, replace=T)
        g_name_26 <- g_name_25[g_name_25 %ni% pool_26]
      }
      if (length(g_name_26) > pool_size){
        pool_27 <- sample(g_name_26, pool_size, replace=T)
        g_name_27 <- g_name_26[g_name_26 %ni% pool_27]
      }
      if (length(g_name_27) > pool_size){
        pool_28 <- sample(g_name_27, pool_size, replace=T)
        g_name_28 <- g_name_27[g_name_27 %ni% pool_28]
      }
      if (length(g_name_28) > pool_size){
        pool_29 <- sample(g_name_28, pool_size, replace=T)
        g_name_29 <- g_name_28[g_name_28 %ni% pool_29]
      }
      if (length(g_name_29) > pool_size){
        pool_30 <- sample(g_name_29, pool_size, replace=T)
        g_name_30 <- g_name_29[g_name_29 %ni% pool_30]
      }
      
      
      #finding those who were in positive pools
      people_who_were_in_positive_pools <- g_name[p_nodes == 1]
      #finding individuals who are isolating
      Individuals_who_are_isolating <- g_name[q_nodes == 1]
      #finding those who were in positive pools who are isolating
      pooled_positive <- people_who_were_in_positive_pools[people_who_were_in_positive_pools %in% Individuals_who_are_isolating]
      
      
      infected_pool_1 <- infected_persons[infected_persons %in% pool_1]
      if (length(infected_pool_1) != 0){
        if (runif(1,0,1) > prob_false_neg){
          p_nodes[pool_1] <- 1
          q_nodes[pool_1] <- 1
        }
      }
      infected_pool_2 <- infected_persons[infected_persons %in% pool_2]
      if (length(infected_pool_2) != 0){
        if (runif(1,0,1) > prob_false_neg){
          p_nodes[pool_2] <- 1
          q_nodes[pool_2] <- 1
        }
      }
      infected_pool_3 <- infected_persons[infected_persons %in% pool_3]
      if (length(infected_pool_3) != 0){
        if (runif(1,0,1) > prob_false_neg){
          p_nodes[pool_3] <- 1
          q_nodes[pool_3] <- 1
        }
      }
       infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
       if (length(infected_pool_4) != 0){
        if (runif(1,0,1) > prob_false_neg){
           p_nodes[pool_4] <- 1
           q_nodes[pool_4] <- 1
         }
       }
       infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
       if (length(infected_pool_5) != 0){
         if (runif(1,0,1) > prob_false_neg){
           p_nodes[pool_5] <- 1
           q_nodes[pool_5] <- 1
         }
       }

      
      test_pools <- tests - pools_every_time
      infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
      infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
      infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
      infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
      infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
      infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
      infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
      infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
      infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
      infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
      infected_pool_16 <- infected_persons[infected_persons %in% pool_16]
      infected_pool_17 <- infected_persons[infected_persons %in% pool_17]
      infected_pool_18 <- infected_persons[infected_persons %in% pool_18]
      infected_pool_19 <- infected_persons[infected_persons %in% pool_19]
      infected_pool_20 <- infected_persons[infected_persons %in% pool_20]
      infected_pool_21 <- infected_persons[infected_persons %in% pool_21]
      infected_pool_22 <- infected_persons[infected_persons %in% pool_22]
      infected_pool_23 <- infected_persons[infected_persons %in% pool_23]
      infected_pool_24 <- infected_persons[infected_persons %in% pool_24]
      infected_pool_25 <- infected_persons[infected_persons %in% pool_25]
      
      #testing those who were in positive pools plus extra pools depending on the number of individuals who were in positive pools waiting for testing
      if (length(pooled_positive) > test_pools){
        testing_pool <- sample(pooled_positive, test_pools, replace=F)
        p_nodes[testing_pool] <- 0
        positive_testing <- infected_persons[infected_persons %in% testing_pool]
        positive_testing_wrong <-c()
        for (i in 1:length(positive_testing))
          if (runif(1,0,1) < prob_false_neg){
            positive_testing_wrong[i] <- positive_testing[i]
          }
        pos_test_right <- positive_testing[positive_testing %ni% positive_testing_wrong]
        negative_testing <- testing_pool[testing_pool %ni% pos_test_right]
        q_nodes[negative_testing] <- 0 
      } else if (length(pooled_positive) == (tests - pools_every_time)){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
      } else if ((tests-length(pooled_positive) - pools_every_time)==1){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
      }else if ((tests-length(pooled_positive) - pools_every_time) == 2){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
      }else if ((tests-length(pooled_positive)- pools_every_time)==3){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
      }else if ((tests-length(pooled_positive)- pools_every_time)==4){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
      }else if ((tests-length(pooled_positive)- pools_every_time)==5){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
      } else if ((tests-length(pooled_positive)- pools_every_time)==6){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==7){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==8){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==9){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==10){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
      }else if ((tests-length(pooled_positive) - pools_every_time)==11){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
        if (length(infected_pool_16) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_16] <- 1
            q_nodes[pool_16] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==12){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
        infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
        if (length(infected_pool_16) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_16] <- 1
            q_nodes[pool_16] <- 1
          }
        }
        if (length(infected_pool_17) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_17] <- 1
            q_nodes[pool_17] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==13){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
        infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
        infected_pool_16 <- infected_persons[infected_persons %in% pool_16]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
        if (length(infected_pool_16) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_16] <- 1
            q_nodes[pool_16] <- 1
          }
        }
        if (length(infected_pool_17) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_17] <- 1
            q_nodes[pool_17] <- 1
          }
        }
        if (length(infected_pool_18) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_18] <- 1
            q_nodes[pool_18] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==14){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
        infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
        infected_pool_16 <- infected_persons[infected_persons %in% pool_16]
        infected_pool_17 <- infected_persons[infected_persons %in% pool_17]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
        if (length(infected_pool_16) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_16] <- 1
            q_nodes[pool_16] <- 1
          }
        }
        if (length(infected_pool_17) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_17] <- 1
            q_nodes[pool_17] <- 1
          }
        }
        if (length(infected_pool_18) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_18] <- 1
            q_nodes[pool_18] <- 1
          }
        }
        if (length(infected_pool_19) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_19] <- 1
            q_nodes[pool_19] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==15){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
        infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
        infected_pool_16 <- infected_persons[infected_persons %in% pool_16]
        infected_pool_17 <- infected_persons[infected_persons %in% pool_17]
        infected_pool_18 <- infected_persons[infected_persons %in% pool_18]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
        if (length(infected_pool_16) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_16] <- 1
            q_nodes[pool_16] <- 1
          }
        }
        if (length(infected_pool_17) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_17] <- 1
            q_nodes[pool_17] <- 1
          }
        }
        if (length(infected_pool_18) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_18] <- 1
            q_nodes[pool_18] <- 1
          }
        }
        if (length(infected_pool_19) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_19] <- 1
            q_nodes[pool_19] <- 1
          }
        }
        if (length(infected_pool_20) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_20] <- 1
            q_nodes[pool_20] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==16){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
        infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
        infected_pool_16 <- infected_persons[infected_persons %in% pool_16]
        infected_pool_17 <- infected_persons[infected_persons %in% pool_17]
        infected_pool_18 <- infected_persons[infected_persons %in% pool_18]
        infected_pool_19 <- infected_persons[infected_persons %in% pool_19]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
        if (length(infected_pool_16) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_16] <- 1
            q_nodes[pool_16] <- 1
          }
        }
        if (length(infected_pool_17) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_17] <- 1
            q_nodes[pool_17] <- 1
          }
        }
        if (length(infected_pool_18) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_18] <- 1
            q_nodes[pool_18] <- 1
          }
        }
        if (length(infected_pool_19) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_19] <- 1
            q_nodes[pool_19] <- 1
          }
        }
        if (length(infected_pool_20) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_20] <- 1
            q_nodes[pool_20] <- 1
          }
        }
        if (length(infected_pool_21) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_21] <- 1
            q_nodes[pool_21] <- 1
          }
        }
      }else if ((tests-length(pooled_positive) - pools_every_time)==17){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
        infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
        infected_pool_16 <- infected_persons[infected_persons %in% pool_16]
        infected_pool_17 <- infected_persons[infected_persons %in% pool_17]
        infected_pool_18 <- infected_persons[infected_persons %in% pool_18]
        infected_pool_19 <- infected_persons[infected_persons %in% pool_19]
        infected_pool_20 <- infected_persons[infected_persons %in% pool_20]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
        if (length(infected_pool_16) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_16] <- 1
            q_nodes[pool_16] <- 1
          }
        }
        if (length(infected_pool_17) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_17] <- 1
            q_nodes[pool_17] <- 1
          }
        }
        if (length(infected_pool_18) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_18] <- 1
            q_nodes[pool_18] <- 1
          }
        }
        if (length(infected_pool_19) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_19] <- 1
            q_nodes[pool_19] <- 1
          }
        }
        if (length(infected_pool_20) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_20] <- 1
            q_nodes[pool_20] <- 1
          }
        }
        if (length(infected_pool_21) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_21] <- 1
            q_nodes[pool_21] <- 1
          }
        }
        if (length(infected_pool_22) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_22] <- 1
            q_nodes[pool_22] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==18){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
        infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
        infected_pool_16 <- infected_persons[infected_persons %in% pool_16]
        infected_pool_17 <- infected_persons[infected_persons %in% pool_17]
        infected_pool_18 <- infected_persons[infected_persons %in% pool_18]
        infected_pool_19 <- infected_persons[infected_persons %in% pool_19]
        infected_pool_20 <- infected_persons[infected_persons %in% pool_20]
        infected_pool_21 <- infected_persons[infected_persons %in% pool_21]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
        if (length(infected_pool_16) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_16] <- 1
            q_nodes[pool_16] <- 1
          }
        }
        if (length(infected_pool_17) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_17] <- 1
            q_nodes[pool_17] <- 1
          }
        }
        if (length(infected_pool_18) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_18] <- 1
            q_nodes[pool_18] <- 1
          }
        }
        if (length(infected_pool_19) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_19] <- 1
            q_nodes[pool_19] <- 1
          }
        }
        if (length(infected_pool_20) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_20] <- 1
            q_nodes[pool_20] <- 1
          }
        }
        if (length(infected_pool_21) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_21] <- 1
            q_nodes[pool_21] <- 1
          }
        }
        if (length(infected_pool_22) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_22] <- 1
            q_nodes[pool_22] <- 1
          }
        }
        if (length(infected_pool_23) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_23] <- 1
            q_nodes[pool_23] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==19){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
        infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
        infected_pool_16 <- infected_persons[infected_persons %in% pool_16]
        infected_pool_17 <- infected_persons[infected_persons %in% pool_17]
        infected_pool_18 <- infected_persons[infected_persons %in% pool_18]
        infected_pool_19 <- infected_persons[infected_persons %in% pool_19]
        infected_pool_20 <- infected_persons[infected_persons %in% pool_20]
        infected_pool_21 <- infected_persons[infected_persons %in% pool_21]
        infected_pool_22 <- infected_persons[infected_persons %in% pool_22]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
        if (length(infected_pool_16) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_16] <- 1
            q_nodes[pool_16] <- 1
          }
        }
        if (length(infected_pool_17) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_17] <- 1
            q_nodes[pool_17] <- 1
          }
        }
        if (length(infected_pool_18) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_18] <- 1
            q_nodes[pool_18] <- 1
          }
        }
        if (length(infected_pool_19) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_19] <- 1
            q_nodes[pool_19] <- 1
          }
        }
        if (length(infected_pool_20) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_20] <- 1
            q_nodes[pool_20] <- 1
          }
        }
        if (length(infected_pool_21) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_21] <- 1
            q_nodes[pool_21] <- 1
          }
        }
        if (length(infected_pool_22) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_22] <- 1
            q_nodes[pool_22] <- 1
          }
        }
        if (length(infected_pool_23) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_23] <- 1
            q_nodes[pool_23] <- 1
          }
        }
        if (length(infected_pool_24) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_24] <- 1
            q_nodes[pool_24] <- 1
          }
        }
      } else if ((tests-length(pooled_positive) - pools_every_time)==20){
        p_nodes[pooled_positive] <- 0 
        pos_test <- infected_persons[infected_persons %in% pooled_positive]
        pos_testing_wrong <-c()
        for (i in 1:length(pos_test))
          if (runif(1,0,1) < prob_false_neg){
            pos_testing_wrong[i] <- pos_test[i]
          }
        pos_test_right <- pos_test[pos_test %ni% pos_testing_wrong]
        neg_test <- pooled_positive[pooled_positive %ni% pos_test_right]
        q_nodes[neg_test] <- 0
        infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
        infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
        infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
        infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
        infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
        infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
        infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
        infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
        infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
        infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
        infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
        infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
        infected_pool_16 <- infected_persons[infected_persons %in% pool_16]
        infected_pool_17 <- infected_persons[infected_persons %in% pool_17]
        infected_pool_18 <- infected_persons[infected_persons %in% pool_18]
        infected_pool_19 <- infected_persons[infected_persons %in% pool_19]
        infected_pool_20 <- infected_persons[infected_persons %in% pool_20]
        infected_pool_21 <- infected_persons[infected_persons %in% pool_21]
        infected_pool_22 <- infected_persons[infected_persons %in% pool_22]
        infected_pool_23 <- infected_persons[infected_persons %in% pool_23]
        if (length(infected_pool_6) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_6] <- 1
            q_nodes[pool_6] <- 1
          }
        }
        if (length(infected_pool_7) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_7] <- 1
            q_nodes[pool_7] <- 1
          }
        }
        if (length(infected_pool_8) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_8] <- 1
            q_nodes[pool_8] <- 1
          }
        }
        if (length(infected_pool_9) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_9] <- 1
            q_nodes[pool_9] <- 1
          }
        }
        if (length(infected_pool_10) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_10] <- 1
            q_nodes[pool_10] <- 1
          }
        }
        if (length(infected_pool_11) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_11] <- 1
            q_nodes[pool_11] <- 1
          }
        }
        if (length(infected_pool_12) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_12] <- 1
            q_nodes[pool_12] <- 1
          }
        }
        if (length(infected_pool_13) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_13] <- 1
            q_nodes[pool_13] <- 1
          }
        }
        if (length(infected_pool_14) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_14] <- 1
            q_nodes[pool_14] <- 1
          }
        }
        if (length(infected_pool_15) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_15] <- 1
            q_nodes[pool_15] <- 1
          }
        }
        if (length(infected_pool_16) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_16] <- 1
            q_nodes[pool_16] <- 1
          }
        }
        if (length(infected_pool_17) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_17] <- 1
            q_nodes[pool_17] <- 1
          }
        }
        if (length(infected_pool_18) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_18] <- 1
            q_nodes[pool_18] <- 1
          }
        }
        if (length(infected_pool_19) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_19] <- 1
            q_nodes[pool_19] <- 1
          }
        }
        if (length(infected_pool_20) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_20] <- 1
            q_nodes[pool_20] <- 1
          }
        }
        if (length(infected_pool_21) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_21] <- 1
            q_nodes[pool_21] <- 1
          }
        }
        if (length(infected_pool_22) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_22] <- 1
            q_nodes[pool_22] <- 1
          }
        }
        if (length(infected_pool_23) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_23] <- 1
            q_nodes[pool_23] <- 1
          }
        }
        if (length(infected_pool_24) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_24] <- 1
            q_nodes[pool_24] <- 1
          }
        }
        if (length(infected_pool_25) != 0){
          if (runif(1,0,1) > prob_false_neg){
            p_nodes[pool_25] <- 1
            q_nodes[pool_25] <- 1
          }
        }
      }
      
      
      


      z_nodes[infected_persons] <- 1
      routine <- g_name[z_nodes == 1]
      
      pos_peeps <- g_name[q_nodes == 1]
      
      
      for (i in 1:length(g_name))
        if (p_nodes[i] == 1 & z_nodes[i]==0 & q_nodes[i]==1){
          w_nodes_info[i] <- w_nodes_info[i] + 1
        }
      

      
      #Isolating individuals for a maximum of 14 days
      for (i in 1:length(g_name))
        if (q_nodes[i] == 1){
          q_nodes_info[i] <- q_nodes_info[i] + 1
          if (q_nodes_info[i] >= 14){
            q_nodes[i] <- 2
          }
        } else {
          q_nodes_info[i] <- q_nodes_info[i]
        }
      
      for (i in 1:length(g_name))
        if (q_nodes[i] == 2 & p_nodes[i]==1){
          q_nodes[i] ==0
        }
      
      
      #if (g_name[q_nodes == 1]){
      # g_name[q_nodes_info] <- g_name[q_nodes_info] + 1
      #}
      
      # if(length(names_of >0)){
      #   positive[time_step,1:length(names_of)] <- names_of  
      # }
    }
    # store new cases of infectiousness
    numnewinfectious <- length(newinfectious)
    if (numnewinfectious >5){
      over_point <-over_point + 1 
    }
    if (numnewinfectious >10){
      over_point_1 <-over_point_1 + 1 
    }
    if (numnewinfectious>0) {
      # Update results
      ord <- match(newinfectious,i_nodes_info[,1])
      results <- rbind(results,cbind(newinfectious,time_step,time_step-as.numeric(i_nodes_info[ord,4]),NA,i_nodes_info[ord,5]))
      numinfectious <- numinfectious+numnewinfectious
    }
    if(length(newremoved)>0){
      results[results[,1]%in%newremoved,4] <- time_step
    }
    
    ## spread infection - 
    if (time_step < 15){
      e_nodes_info <- spread_wrapper_2(i_nodes_info,s_nodes,v_nodes,e_nodes_info, pos_peeps, direct_VE)
      s_nodes[e_nodes_info[,1]] <- 0
      e_nodes[e_nodes_info[,1]] <- 1
      order_infected_2 <- c(order_infected,e_nodes_info[,1])
    } else {
      e_nodes_info <- spread_wrapper(i_nodes_info,s_nodes,v_nodes,e_nodes_info, pos_peeps, direct_VE)
      s_nodes[e_nodes_info[,1]] <- 0
      e_nodes[e_nodes_info[,1]] <- 1
      order_infected_2 <- c(order_infected,e_nodes_info[,1])
    }
    
    # infect from source
    rate_from_source <- max((start_day + time_step)*from_source + base_rate, 0)
    if(rate_from_source>0){
      e_nodes_info <- infect_from_source(s_nodes,v_nodes,e_nodes_info,direct_VE,incperiod_shape,incperiod_rate,rate_from_source)
      s_nodes[e_nodes_info[,1]] <- 0
      e_nodes[e_nodes_info[,1]] <- 1
      order_infected <- c(order_infected,e_nodes_info[,1])
    }
    
    # store information
    trajectories$S[time_step+1] <- sum(s_nodes)# + sum(v_nodes) + sum(c_nodes)
    trajectories$E[time_step+1] <- trajectories$S[time_step] - trajectories$S[time_step+1]
    trajectories$I[time_step+1] <- numnewinfectious
    trajectories$R[time_step+1] <- sum(r_nodes)-sum(trajectories$R)
    
  }
  
  
  # store information and format to return
  results <- as.data.frame(results)
  colnames(results) <- c('InfectedNode', 'DayInfectious', 'DayInfected', 'DayRemoved','Observed')
  results$inCluster <- results$InfectedNode%in%cluster_people
  results$contact <- results$InfectedNode%in%contacts
  #results$highrisk <- results$InfectedNode%in%high_risk
  results$inTrial <- results$InfectedNode%in%trial_participants
  results$vaccinated <- results$InfectedNode%in%vaccinees
  results$RecruitmentDay <- recruitment_times[match(results$InfectedNode,trial_participants)]
  peaks <- rep(0, 500) 
  
  for (i in 11:488){
    if (report[i]>report[i+1]){
      if (report[i] >report[i+2]){
        if (report[i] > report[i+3]){
          if (report[i] > report[i+4]){
            if (report[i] > report[i+5]){
              if (report[i] > report[i+6]){
                if (report[i] > report[i+7]){
                  if (report[i] > report[i+8]){
                    if (report[i] > report[i+9]){
                      if (report[i] > report[i+10]){
                        if (report[i] > report[i-10]){
                          if (report[i] > report[i-9]){
                            if (report[i] > report[i-8]){
                              if (report[i] > report[i-7]){
                                if (report[i] > report[i-6]){
                                  if (report[i] > report[i-5]){
                                    if (report[i] > report[i-4]){
                                      if (report[i] > report[i-3]){
                                        if (report[i] > report[i-2]){
                                          if (report[i] > report[i-1]){
                                            peaks[i] <- report[i]
                                            
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }else{
      peaks[i] <- 0
    }
    
  } 
  not_infect <-sum(report == 0)
  not_peaks <- sum(peaks == 0)
  
  resu <- as.data.frame(cbind(length(routine), length(pos_peeps)))
  infections <- length(routine)
  isolations <- sum(q_nodes_info)
  peak <-max(report)
  infection_period<- 500 - not_infect
  dif_peaks <- 500 - not_peaks
  unneccesary_isolations<-sum(w_nodes_info)
  
  
  
  
  
  return(list(results,length(cluster_people),recruitment_times,length(vaccinees),
              length(trial_participants),vaccinees,trial_participants,order_infected,
              vaccine_incubation_times+recruitment_times[trial_participants%in%vaccinees],
              peak, infections, isolations, infection_period, dif_peaks,over_point, over_point_1, unneccesary_isolations))
  
}

