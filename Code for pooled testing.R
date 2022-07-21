#### 5 tests ####

if (length(pooled_positive) > 5){
  testing_pool <- sample(pooled_positive, 5, replace = F)
  p_nodes[testing_pool] <- 0 
  positive_testing <- infected_persons[infected_persons %in% testing_pool]
  negative_testing <- testing_pool[testing_pool %ni% positive_testing]
  q_nodes[negative_testing] <- 0 
} else if (length(pooled_positive) == 5){
  p_nodes[pooled_positive] <- 0 
  pos_test <- infected_persons[infected_persons %in% pooled_positive]
  neg_test <- pooled_positive[pooled_positive %ni% pos_test]
  q_nodes[neg_test] <- 0
}else if (length(pooled_positive) == 0 ){
  infected_pool_1 <- infected_persons[infected_persons %in% pool_1]
  if (length(infected_pool_1) != 0){
    p_nodes[pool_1] <- 1
    q_nodes[pool_1] <- 1
  }
  
  infected_pool_2 <- infected_persons[infected_persons %in% pool_2]
  if (length(infected_pool_2) != 0){
    p_nodes[pool_2] <- 1
    q_nodes[pool_2] <- 1
  }
  
  infected_pool_3 <- infected_persons[infected_persons %in% pool_3]
  if (length(infected_pool_3) != 0){
    p_nodes[pool_3] <- 1
    q_nodes[pool_3] <- 1
  }
  
  infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
  if (length(infected_pool_4) != 0){
    p_nodes[pool_4] <- 1
    q_nodes[pool_4] <- 1
  }
  infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
  if (length(infected_pool_5) != 0){
    p_nodes[pool_5] <- 1
    q_nodes[pool_5] <- 1
  }
}


#### 10 tests #####

if (length(pooled_positive) > 10){
  testing_pool <- sample(pooled_positive, 10, replace = F)
  p_nodes[testing_pool] <- 0 
  positive_testing <- infected_persons[infected_persons %in% testing_pool]
  negative_testing <- testing_pool[testing_pool %ni% positive_testing]
  q_nodes[negative_testing] <- 0 
} else if (length(pooled_positive) == 10){
  p_nodes[pooled_positive] <- 0 
  pos_test <- infected_persons[infected_persons %in% pooled_positive]
  neg_test <- pooled_positive[pooled_positive %ni% pos_test]
  q_nodes[neg_test] <- 0
}else if (length(pooled_positive) == 0 ){
  infected_pool_1 <- infected_persons[infected_persons %in% pool_1]
  if (length(infected_pool_1) != 0){
    p_nodes[pool_1] <- 1
    q_nodes[pool_1] <- 1
  }
  
  infected_pool_2 <- infected_persons[infected_persons %in% pool_2]
  if (length(infected_pool_2) != 0){
    p_nodes[pool_2] <- 1
    q_nodes[pool_2] <- 1
  }
  
  infected_pool_3 <- infected_persons[infected_persons %in% pool_3]
  if (length(infected_pool_3) != 0){
    p_nodes[pool_3] <- 1
    q_nodes[pool_3] <- 1
  }
  
  infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
  if (length(infected_pool_4) != 0){
    p_nodes[pool_4] <- 1
    q_nodes[pool_4] <- 1
  }
  infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
  if (length(infected_pool_5) != 0){
    p_nodes[pool_5] <- 1
    q_nodes[pool_5] <- 1
  }
  infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
  if (length(infected_pool_6) != 0){
    p_nodes[pool_6] <- 1
    q_nodes[pool_6] <- 1
  }
  infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
  if (length(infected_pool_7) != 0){
    p_nodes[pool_7] <- 1
    q_nodes[pool_7] <- 1
  }
  infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
  if (length(infected_pool_8) != 0){
    p_nodes[pool_8] <- 1
    q_nodes[pool_8] <- 1
  }
  infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
  if (length(infected_pool_9) != 0){
    p_nodes[pool_9] <- 1
    q_nodes[pool_9] <- 1
  }
  infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
  if (length(infected_pool_10) != 0){
    p_nodes[pool_10] <- 1
    q_nodes[pool_10] <- 1
  }
}


###### 15 tests ######

if (length(pooled_positive) > 15){
  testing_pool <- sample(pooled_positive, 15, replace = F)
  p_nodes[testing_pool] <- 0 
  positive_testing <- infected_persons[infected_persons %in% testing_pool]
  negative_testing <- testing_pool[testing_pool %ni% positive_testing]
  q_nodes[negative_testing] <- 0 
} else if (length(pooled_positive) == 15){
  p_nodes[pooled_positive] <- 0 
  pos_test <- infected_persons[infected_persons %in% pooled_positive]
  neg_test <- pooled_positive[pooled_positive %ni% pos_test]
  q_nodes[neg_test] <- 0
}else if (length(pooled_positive) == 5){
  p_nodes[pooled_positive] <- 0 
  pos_test <- infected_persons[infected_persons %in% pooled_positive]
  neg_test <- pooled_positive[pooled_positive %ni% pos_test]
  q_nodes[neg_test] <- 0
  infected_pool_1 <- infected_persons[infected_persons %in% pool_1]
  if (length(infected_pool_1) != 0){
    p_nodes[pool_1] <- 1
    q_nodes[pool_1] <- 1
  }
  
  infected_pool_2 <- infected_persons[infected_persons %in% pool_2]
  if (length(infected_pool_2) != 0){
    p_nodes[pool_2] <- 1
    q_nodes[pool_2] <- 1
  }
  
  infected_pool_3 <- infected_persons[infected_persons %in% pool_3]
  if (length(infected_pool_3) != 0){
    p_nodes[pool_3] <- 1
    q_nodes[pool_3] <- 1
  }
  
  infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
  if (length(infected_pool_4) != 0){
    p_nodes[pool_4] <- 1
    q_nodes[pool_4] <- 1
  }
  infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
  if (length(infected_pool_5) != 0){
    p_nodes[pool_5] <- 1
    q_nodes[pool_5] <- 1
  }
  infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
  if (length(infected_pool_6) != 0){
    p_nodes[pool_6] <- 1
    q_nodes[pool_6] <- 1
  }
  infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
  if (length(infected_pool_7) != 0){
    p_nodes[pool_7] <- 1
    q_nodes[pool_7] <- 1
  }
  infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
  if (length(infected_pool_8) != 0){
    p_nodes[pool_8] <- 1
    q_nodes[pool_8] <- 1
  }
  infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
  if (length(infected_pool_9) != 0){
    p_nodes[pool_9] <- 1
    q_nodes[pool_9] <- 1
  }
  infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
  if (length(infected_pool_10) != 0){
    p_nodes[pool_10] <- 1
    q_nodes[pool_10] <- 1
  }
}else if (length(pooled_positive) == 0 ){
  infected_pool_1 <- infected_persons[infected_persons %in% pool_1]
  if (length(infected_pool_1) != 0){
    p_nodes[pool_1] <- 1
    q_nodes[pool_1] <- 1
  }
  
  infected_pool_2 <- infected_persons[infected_persons %in% pool_2]
  if (length(infected_pool_2) != 0){
    p_nodes[pool_2] <- 1
    q_nodes[pool_2] <- 1
  }
  
  infected_pool_3 <- infected_persons[infected_persons %in% pool_3]
  if (length(infected_pool_3) != 0){
    p_nodes[pool_3] <- 1
    q_nodes[pool_3] <- 1
  }
  
  infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
  if (length(infected_pool_4) != 0){
    p_nodes[pool_4] <- 1
    q_nodes[pool_4] <- 1
  }
  infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
  if (length(infected_pool_5) != 0){
    p_nodes[pool_5] <- 1
    q_nodes[pool_5] <- 1
  }
  infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
  if (length(infected_pool_6) != 0){
    p_nodes[pool_6] <- 1
    q_nodes[pool_6] <- 1
  }
  
  infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
  if (length(infected_pool_7) != 0){
    p_nodes[pool_7] <- 1
    q_nodes[pool_7] <- 1
  }
  
  infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
  if (length(infected_pool_8) != 0){
    p_nodes[pool_8] <- 1
    q_nodes[pool_8] <- 1
  }
  
  infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
  if (length(infected_pool_9) != 0){
    p_nodes[pool_9] <- 1
    q_nodes[pool_9] <- 1
  }
  infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
  if (length(infected_pool_10) != 0){
    p_nodes[pool_10] <- 1
    q_nodes[pool_10] <- 1
  }
  infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
  if (length(infected_pool_11) != 0){
    p_nodes[pool_11] <- 1
    q_nodes[pool_11] <- 1
  }
  infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
  if (length(infected_pool_12) != 0){
    p_nodes[pool_12] <- 1
    q_nodes[pool_12] <- 1
  }
  infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
  if (length(infected_pool_13) != 0){
    p_nodes[pool_13] <- 1
    q_nodes[pool_13] <- 1
  }
  infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
  if (length(infected_pool_14) != 0){
    p_nodes[pool_14] <- 1
    q_nodes[pool_14] <- 1
  }
  infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
  if (length(infected_pool_15) != 0){
    p_nodes[pool_15] <- 1
    q_nodes[pool_15] <- 1
  }
}

#### 20 tests ####

if (length(pooled_positive) > 20){
  testing_pool <- sample(pooled_positive, 20, replace = F)
  p_nodes[testing_pool] <- 0 
  positive_testing <- infected_persons[infected_persons %in% testing_pool]
  negative_testing <- testing_pool[testing_pool %ni% positive_testing]
  q_nodes[negative_testing] <- 0 
} else if (length(pooled_positive) == 20){
  p_nodes[pooled_positive] <- 0 
  pos_test <- infected_persons[infected_persons %in% pooled_positive]
  neg_test <- pooled_positive[pooled_positive %ni% pos_test]
  q_nodes[neg_test] <- 0
}else if (length(pooled_positive) == 10){
  p_nodes[pooled_positive] <- 0 
  pos_test <- infected_persons[infected_persons %in% pooled_positive]
  neg_test <- pooled_positive[pooled_positive %ni% pos_test]
  q_nodes[neg_test] <- 0
  infected_pool_1 <- infected_persons[infected_persons %in% pool_1]
  if (length(infected_pool_1) != 0){
    p_nodes[pool_1] <- 1
    q_nodes[pool_1] <- 1
  }
  
  infected_pool_2 <- infected_persons[infected_persons %in% pool_2]
  if (length(infected_pool_2) != 0){
    p_nodes[pool_2] <- 1
    q_nodes[pool_2] <- 1
  }
  
  infected_pool_3 <- infected_persons[infected_persons %in% pool_3]
  if (length(infected_pool_3) != 0){
    p_nodes[pool_3] <- 1
    q_nodes[pool_3] <- 1
  }
  
  infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
  if (length(infected_pool_4) != 0){
    p_nodes[pool_4] <- 1
    q_nodes[pool_4] <- 1
  }
  infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
  if (length(infected_pool_5) != 0){
    p_nodes[pool_5] <- 1
    q_nodes[pool_5] <- 1
  }
  infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
  if (length(infected_pool_6) != 0){
    p_nodes[pool_6] <- 1
    q_nodes[pool_6] <- 1
  }
  infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
  if (length(infected_pool_7) != 0){
    p_nodes[pool_7] <- 1
    q_nodes[pool_7] <- 1
  }
  infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
  if (length(infected_pool_8) != 0){
    p_nodes[pool_8] <- 1
    q_nodes[pool_8] <- 1
  }
  infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
  if (length(infected_pool_9) != 0){
    p_nodes[pool_9] <- 1
    q_nodes[pool_9] <- 1
  }
  infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
  if (length(infected_pool_10) != 0){
    p_nodes[pool_10] <- 1
    q_nodes[pool_10] <- 1
  }
}else if (length(pooled_positive) == 0 ){
  infected_pool_1 <- infected_persons[infected_persons %in% pool_1]
  if (length(infected_pool_1) != 0){
    p_nodes[pool_1] <- 1
    q_nodes[pool_1] <- 1
  }
  
  infected_pool_2 <- infected_persons[infected_persons %in% pool_2]
  if (length(infected_pool_2) != 0){
    p_nodes[pool_2] <- 1
    q_nodes[pool_2] <- 1
  }
  
  infected_pool_3 <- infected_persons[infected_persons %in% pool_3]
  if (length(infected_pool_3) != 0){
    p_nodes[pool_3] <- 1
    q_nodes[pool_3] <- 1
  }
  
  infected_pool_4 <- infected_persons[infected_persons %in% pool_4]
  if (length(infected_pool_4) != 0){
    p_nodes[pool_4] <- 1
    q_nodes[pool_4] <- 1
  }
  infected_pool_5 <- infected_persons[infected_persons %in% pool_5]
  if (length(infected_pool_5) != 0){
    p_nodes[pool_5] <- 1
    q_nodes[pool_5] <- 1
  }
  infected_pool_6 <- infected_persons[infected_persons %in% pool_6]
  if (length(infected_pool_6) != 0){
    p_nodes[pool_6] <- 1
    q_nodes[pool_6] <- 1
  }
  
  infected_pool_7 <- infected_persons[infected_persons %in% pool_7]
  if (length(infected_pool_7) != 0){
    p_nodes[pool_7] <- 1
    q_nodes[pool_7] <- 1
  }
  
  infected_pool_8 <- infected_persons[infected_persons %in% pool_8]
  if (length(infected_pool_8) != 0){
    p_nodes[pool_8] <- 1
    q_nodes[pool_8] <- 1
  }
  
  infected_pool_9 <- infected_persons[infected_persons %in% pool_9]
  if (length(infected_pool_9) != 0){
    p_nodes[pool_9] <- 1
    q_nodes[pool_9] <- 1
  }
  infected_pool_10 <- infected_persons[infected_persons %in% pool_10]
  if (length(infected_pool_10) != 0){
    p_nodes[pool_10] <- 1
    q_nodes[pool_10] <- 1
  }
  infected_pool_11 <- infected_persons[infected_persons %in% pool_11]
  if (length(infected_pool_11) != 0){
    p_nodes[pool_11] <- 1
    q_nodes[pool_11] <- 1
  }
  infected_pool_12 <- infected_persons[infected_persons %in% pool_12]
  if (length(infected_pool_12) != 0){
    p_nodes[pool_12] <- 1
    q_nodes[pool_12] <- 1
  }
  infected_pool_13 <- infected_persons[infected_persons %in% pool_13]
  if (length(infected_pool_13) != 0){
    p_nodes[pool_13] <- 1
    q_nodes[pool_13] <- 1
  }
  infected_pool_14 <- infected_persons[infected_persons %in% pool_14]
  if (length(infected_pool_14) != 0){
    p_nodes[pool_14] <- 1
    q_nodes[pool_14] <- 1
  }
  infected_pool_15 <- infected_persons[infected_persons %in% pool_15]
  if (length(infected_pool_15) != 0){
    p_nodes[pool_15] <- 1
    q_nodes[pool_15] <- 1

  }
}

