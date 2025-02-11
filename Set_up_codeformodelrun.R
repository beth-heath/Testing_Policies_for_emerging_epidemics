#setwd('')
setwd('/Users/bethanyheath/OneDrive/bethany/PhD/Code for submission')
library(igraph)
library(truncnorm)
library(infotheo)
library(xtable)
library(RColorBrewer)
library(plotrix)
library(profvis)
library(funique)
library(doParallel)
library(foreach)
library(survival)
library(coxme)
library(pracma)
library(dplyr)
library(readxl)
#library(rethinking)
library(rBeta2009)



## build network ###############################################################
#source('Network_set.R')
`%ni%` <- Negate(`%in%`)

## functions #######################################################

source('Random_pooled_10_tests.R')

ref_recruit_day <<- 30
eval_day <<- 25
target_weight <<- 24
ve_est_threshold <<- 0.005
break_threshold <<- 5

#Spread code in the less social interactions case

covid_spread_wrapper <- function(i_nodes_info,s_nodes,v_nodes,e_nodes_info,isolation_individuals, direct_VE){
  # to contacts
  # e infects house and work and anyone - only enodes infected one day ago or more, and only enodes with one day left incubating
  ##!! a subset of i_nodes are nonsymptomatic and therefore continue to infect contacts. these should be a fixed list, not sampled randomly every time.
  current_infectious <- c(i_nodes_info[,1],e_nodes_info[e_nodes_info[,2]>=e_nodes_info[,3],1])
  available_infectors <- current_infectious[current_infectious %ni% isolation_individuals]
  #non_isolating <- current_infectious %ni% isolation_individuals
  #available_infectors <- current_infectious[non_isolating]
  if(length(available_infectors)>0){
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,available_infectors,direct_VE,incperiod_shape,incperiod_rate,susc_list=contact_list,beta_scalar=nonrandom_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,available_infectors,direct_VE,incperiod_shape,incperiod_rate,susc_list=random_list,beta_scalar=random_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
  }
  #i infects house - changing code so that oninfected_personsly those that have tested positive - currently have that testing occurs every day
  infected_and_isolated <- current_infectious[current_infectious %in% isolation_individuals]
  if(length(isolation_individuals)>0){
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info, infected_and_isolated,direct_VE,incperiod_shape,incperiod_rate,susc_list=household_list,beta_scalar=nonrandom_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
  }
  return(e_nodes_info)
  
}

#Spread in the case with more social interactions
covid_spread_wrapper_2 <- function(i_nodes_info,s_nodes,v_nodes,e_nodes_info,isolation_individuals, direct_VE){
  # to contacts
  # e infects house and work and anyone - only enodes infected one day ago or more, and only enodes with one day left incubating
  ##!! a subset of i_nodes are nonsymptomatic and therefore continue to infect contacts. these should be a fixed list, not sampled randomly every time.
  current_infectious <- c(i_nodes_info[,1],e_nodes_info[e_nodes_info[,2]>=e_nodes_info[,3],1])
  available_infectors <- current_infectious[current_infectious %ni% isolation_individuals]
  if(length(current_infectious)>0){
 #   e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=contact_list,beta_scalar=nonrandom_scalar)
#    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=random_list,beta_scalar=random_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=social_list,beta_scalar=nonrandom_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
  }
  #i infects house - changing code so that only those that have tested positive - currently have that testing occurs every day
  infected_and_isolated <- current_infectious[current_infectious %in% isolation_individuals]
  if(length(infected_and_isolated)>0){
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,isolation_individuals,direct_VE,incperiod_shape,incperiod_rate,susc_list=household_list,beta_scalar=0.2)
  }
  return(e_nodes_info)
}



## set up #######################################################

# Per-time-step hazard of infection for a susceptible nodes from an infectious
# neighbour
beta_base <<- 0.01
# Gamma-distribution parameters of incubation and infectious period and wait times
# hist(rgamma(1000,shape=1.43,rate=0.549)+2)
infperiod_shape <<- 1.43
infperiod_rate <<- 0.549
infperiod_const <<- 2
## assume there is no difference between infectious with and without symptoms - all I
#mn <- 6.5 # =shape/rate # 5.2
#sd <- 2.6 # =sqrt(shape/rate^2) # 2.8
# hist(rgamma(1000,shape=13.3,rate=4.16)+2)
incperiod_rate <<- 4.16 # mn/sd^2 # 
incperiod_shape <<- 13.3 # incperiod_rate*mn # 
incperiod_const <<- 2
#hosp_shape_index <<- 2
#hosp_rate_index <<- 0.5
#hosp_shape <<- 2
#hosp_rate <<- 0.5
recruit_shape <<- 5.4
recruit_rate <<- 0.47
#hosp_mean_index <<- 3.85
#hosp_sd_index <<- 2.76
#hosp_mean <<- 3.85
#hosp_sd <<- 2.76
recruit_mean <<- 10.32
recruit_sd <<- 4.79
enrollment_rate <<- 0
nonrandom_scalar <<- 1
random_scalar <<- 1/10
length(E(new_g))
length(E(random_g))
direct_VE <- 0.0

g <<- new_g

g_name <<- as.numeric(as.vector(V(g)$name))
vertices <- V(g)
cluster_size <- hosp_times <- recruit_times <- c()
results_list <- list()


###########################################################################


set_variables_from_gamma_distributions <- function(){
  vacc_shape <<- 3
  vacc_rate <<- 1
  
  infperiod_scale <- 1/infperiod_rate
  incperiod_scale <- 1/incperiod_rate
  vacc_scale <- 1/vacc_rate
  
  mu <- incperiod_shape*incperiod_scale + vacc_shape*vacc_scale
  sig2 <- incperiod_shape*incperiod_scale^2 + vacc_shape*vacc_scale^2
  alpha <- mu^2/sig2
  beta <- sig2/mu
  inc_plus_vacc_shape <<- alpha
  inc_plus_vacc_rate <<- 1/beta
  
  zero <<- 50
  pgamma_vector <<- pgamma(1:100-incperiod_const,shape=inc_plus_vacc_shape,rate=inc_plus_vacc_rate)
  dgamma_vector <<- dgamma(1:100-incperiod_const,shape=inc_plus_vacc_shape,rate=inc_plus_vacc_rate)
  pgamma_inc_vector <<- c(rep(0,zero),pgamma(1:100-incperiod_const,shape=incperiod_shape,rate=incperiod_rate))
  pgamma_vacc_vector <<- c(rep(0,zero),pgamma(1:100,shape=vacc_shape,rate=vacc_rate))
  dgamma_inc_vector <<- c(rep(0,zero),dgamma(1:100-incperiod_const,shape=incperiod_shape,rate=incperiod_rate))
  
  recruitment_time <<- 30
}
set_variables_from_gamma_distributions()

#### Additional code ####
nIter <- 100
e_order <- list()
infected_1<-infected_1<-length_1<-no_peaks_1<- peak_1<-isolated_1<-unneccesary_infections1<-threshold_11 <-threshold_21 <- matrix(nrow = 15, ncol=1)
sdinfected_1<-sdinfected_1<-sdlength_1<-sdno_peaks_1<- sdpeak_1<-sdisolated_1<-sdunneccesary_infections1<-sdthreshold_11 <-sdthreshold_21 <-matrix(nrow =15, ncol=1)
#data_collected = matrix (nrow = length(nIter), ncol = 2)
infected_2<- isolated_2<- peak_2<-length_2<-no_peaks_2 <-unneccesary_infections <- threshold_2<-threshold_1<-c()
#profvis(

#delay_prob_range <-c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

delay_prob_range <-c(2,4,6,8,10,12,14,16,20,22, 24, 26, 28, 30, 32)
other_symps_range <-c(1)



for (i in delay_prob_range){
 for (j in other_symps_range){
   for (iter in 1:nIter){
     #observed <<- delay_prob_range
     observed <<- i
# select random person to start
     first_infected <- sample(g_name, 20)
     inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
     netwk <- simulate_contact_network(first_infected,start_day=iter,from_source=0,cluster_flag=0,individual_recruitment_times=T,spread_wrapper=covid_spread_wrapper,
                                       prob_false_neg = 0.279, tests=10,
                                       non_compliance_prob =0, pools_every_time=5, pool_size=i


     )

     results_list[[iter]] <- netwk[[1]]
     cluster_size[iter] <- netwk[[2]]
     peak_2[iter] <- netwk[[10]]
     infected_2[iter] <- netwk[[11]]
     isolated_2[iter] <- netwk[[12]]
     length_2[iter] <- netwk[[13]]
     no_peaks_2[iter] <- netwk[[14]]
     unneccesary_infections[iter]<-netwk[[17]]
     threshold_1[iter] <-netwk[[15]]
     threshold_2[iter] <-netwk[[16]]
#recruit_times[iter] <- netwk[[3]]
#    e_order[[iter]] <- netwk[[8]][!duplicated(netwk[[8]])]
#data_collected[iter, 1] = netwk[[11]]
#data_collected[iter, 2] = netwk[[12]]
   }

    infected_1[which(delay_prob_range==i),which(other_symps_range==j)]<- mean(infected_2)
    sdinfected_1[which(delay_prob_range==i),which(other_symps_range==j)] <-sd(infected_2)
    length_1[which(delay_prob_range==i),which(other_symps_range==j)]<-mean(length_2)
    sdlength_1[which(delay_prob_range==i),which(other_symps_range==j)]<- sd(length_2)
    no_peaks_1[which(delay_prob_range==i),which(other_symps_range==j)]<- mean(no_peaks_2)
    sdno_peaks_1[which(delay_prob_range==i),which(other_symps_range==j)]<-sd(no_peaks_2)
    peak_1[which(delay_prob_range==i),which(other_symps_range==j)] <- mean(peak_2)
    sdpeak_1[which(delay_prob_range==i),which(other_symps_range==j)]<- sd(peak_2)
    isolated_1[which(delay_prob_range==i),which(other_symps_range==j)]<- mean(isolated_2)
    sdisolated_1[which(delay_prob_range==i),which(other_symps_range==j)] <- sd(isolated_2)
    unneccesary_infections1[which(delay_prob_range==i),which(other_symps_range==j)] <- mean(unneccesary_infections)
    sdunneccesary_infections1[which(delay_prob_range==i),which(other_symps_range==j)]<- sd(unneccesary_infections)
    threshold_11[which(delay_prob_range==i),which(other_symps_range==j)] <- mean(threshold_1)
    sdthreshold_11[which(delay_prob_range==i),which(other_symps_range==j)] <- sd(threshold_1)
    threshold_21[which(delay_prob_range==i),which(other_symps_range==j)] <- mean(threshold_2)
    sdthreshold_21[which(delay_prob_range==i),which(other_symps_range==j)]<- sd(threshold_2)
    print(infected_1)
  }
}



