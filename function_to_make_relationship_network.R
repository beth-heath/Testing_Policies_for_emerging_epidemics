`%ni%` <- Negate(`%in%`)
set.seed(123)
#loading in the data
occupancy_data <- read_xls('/Users/bethanyheath/Downloads/Dat.xls',sheet=2)
#Naming the columns after the variables
colnames(occupancy_data) <- occupancy_data[9,]
#Removing the bit at the top of the document that was not needed
occupancy_data <- occupancy_data[-c(1:9),]
#renaming column name
occupancy_data$`Household Size`[occupancy_data$`Household Size`=="11 or more people in household"] <- '11 people or more'
#sums up the number of people in the household
occupancy_data$number <- as.numeric(sapply(occupancy_data$`Household Size`,function(x)strsplit(x,' p')[[1]][1]))
#getting rid of the blank NAs role from inbetween the i and i+1 sections
occupancy_data <- subset(occupancy_data,!is.na(number))
#The 11 plus section does not have the totals written in certain places so this is including them in
occupancy_data[286,5] <- as.character(6432)
occupancy_data[286,7] <- as.character(occupancy_data[286,8])
#Applying across the columns it takes the subset of whose number equals a certain value such as 1 takes their total and sums it. Essentially it is working out the number of households that have a certain number of inidividuals.
households <- sapply(1:10,function(y)sum(as.numeric(subset(occupancy_data,number==y)$Total)))
#This is doing a similar job to the previous question but is multiplying the household by the number of people.
people <- sapply(1:10,function(y)sum(as.numeric(subset(occupancy_data,number==y)$Total))*y)
#Just assigning a number to each of the cases.
occupancy_data$type <- 1:nrow(occupancy_data)

communal <- read_xls('/Users/bethanyheath/Downloads/communal_living.xls',sheet=1)
economic <- read_xls('/Users/bethanyheath/Downloads/communal_living_economic.xls',sheet=1)

#Previous section was about loading in the data that will be used to build the model. 

#Next stage is about creating the network that will be used to build the model.

#set the number of households to 500
number_of_households <- 500
#set up the type for these 500 household by choosing based on the probability taken from the data already loaded in
household_types <- sample(occupancy_data$type,number_of_households,replace=T,prob=occupancy_data$Total)
#renaming variables
colnames(occupancy_data)[1:4] <- c('description','children','adults','elderly')
#setting the size of the household based on what type of household was chosen
household_sizes <- occupancy_data$number[household_types]

label_start <- 0
#hh is a list with the number of households there are then says to make a  full graph then set vertex attribution to the initial start plus household sizes. 
hh <- list()
for(i in 1:number_of_households) {
  hh[[i]] <- make_full_graph(household_sizes[i]) %>%
    set_vertex_attr("name", value = label_start+1:household_sizes[i])
  label_start <- label_start + household_sizes[i]
}
# extract household data frames
attrs <- do.call(rbind,lapply(hh,function(x)igraph::as_data_frame(x,'vertices')))
# combine all
el <- do.call(rbind,lapply(hh,function(x)igraph::as_data_frame(x)))
# convert to network
new_g <- graph_from_data_frame(el, directed = FALSE, vertices = attrs)
# save layout for plotting
pts <- 80
save_layout <- layout_nicely(induced.subgraph(new_g,1:pts))
# add household labels
hh_labels <- rep(1:number_of_households,household_sizes)
#setting the vertex attributes to what has been worked out 
new_g <- set_vertex_attr(new_g,'hh',value=hh_labels)
#ego finds the vertices not further from the neighbourhood of the vertices
household_list <<- lapply(V(new_g),function(x) {cs <- as.vector(unlist(ego(new_g,order=1,nodes=x))); cs[cs!=x]})
#combining the houses together including the children adults and eldery together.
house_makeup <- lapply(2:nrow(occupancy_data)-1,function(x)as.numeric(c(occupancy_data$children[occupancy_data$type==x],
                                                                        occupancy_data$adults[occupancy_data$type==x],
                                                                        occupancy_data$elderly[occupancy_data$type==x])))
#for the residents they have find the number of the residents for each of the number of households. Setting up measures given which household are together.                                                                         occupancy_data$elderly[occupancy_data$type==x])))
demographic_index <- rep(0,length(V(new_g)))
for(i in 1:number_of_households){
  residents <- which(hh_labels==i)
  occupants <- house_makeup[[household_types[i]]]
  labels <- rep(1:3,times=occupants)
  demographic_index[residents] <- labels
}
demographic_index <<- demographic_index

#code to assign ages to the children, adults and eldery sections
min_age <- c(0,10,20,30,40,50,60,70,80)
max_age <- c(9,19,29,39,49,59,69,79,150)
cfr <- c(0.002,0.006,0.03,0.08,0.15,0.6,2.2,5.1,9.3)/100
grouped_cfr <<- c(mean(cfr[1:2]),mean(cfr[3:7]),mean(cfr[8:9]))
tune <<- grouped_cfr/max(grouped_cfr)

#states the number of indivduals in the samples that are children, adults and elderly.
paste0(sum(demographic_index==1),' children')
paste0(sum(demographic_index==2),' adults')
paste0(sum(demographic_index==3),' elderly')

# Next section is about creating work connections between individuals. In the current code this does not take into account: children and eldery not working, communal/institutional living, social

#Next section is counting how many individuals are working
#Create a new measure that will be used to count this up
worker_index <- rep(0,length(demographic_index))
#Next each adult is considered to be in the working population. This is easily adjustable by adding a probability into this.
worker_index[demographic_index==2] <- 1
#Elders are considered in the population but only those that would be below a certain age therefore to take this into account, we get:
worker_index[demographic_index==3&runif(length(worker_index))<0.2] <- 1
#Calculating the number working from summing this measure calculated
n_adults <- sum(worker_index)
# number of workplaces calculated using a rpois distribution using information from https://www.hse.gov.uk/contact/faqs/toilets.htm
number_workplaces <- rpois(1,n_adults/15)
#setting up for the workplace index
workplace_index <- rep(0,length(V(new_g)))
#takes a sample from the number of workplaces for the number of adults who work there are in the sample therefore each working member gets assigned a number for their workplace
workplace_index[worker_index==1] <- sample(1:number_workplaces,n_adults,replace=T)
#add in work connections into the model which works by looking at people who work in the same workplace adds connections between them.
for(i in 1:number_workplaces) {
  workers <- which(workplace_index==i)
  if(length(workers)>1)
    for(j in 2:length(workers))
      for(k in 1:(j-1)){
        new_g <- add_edges(new_g,edges=c(workers[j],workers[k]))
      }
}

## Children in classrooms
young_people <- V(new_g)[demographic_index==1]
Number_classes <- rpois(1, length(young_people)/25)
School_index <- rep(0, length(V(new_g)))
School_index[demographic_index==1] <- sample(1:Number_classes, length(young_people), replace=T)
for(i in 1:Number_classes) {
  pupils <- which(School_index==i)
  if(length(pupils)>1)
    for(j in 2:length(pupils))
      for(k in 1:(j-1)){
        new_g <- add_edges(new_g,edges=c(pupils[j],pupils[k]))
      }
}

par(mar=c(1,1,1,1))
plot(induced.subgraph(new_g,1:pts),layout=save_layout,
     vertex.label=NA,vertex.size=10,vertex.color=c('white',"#660066","#CC99FF")[demographic_index[1:pts]])
legend('topleft',legend=c('<19','19-65','>65'),pt.cex=2,col='black',pch=21, pt.bg=c('white',"#660066","#CC99FF"),bty='n')

#Plotting graph of the connections
#plot.igraph(new_g,vertex.label=NA,vertex.size=1,layout=save_layout)

#Looking at the size of the clusters
cluster_sizes <- sapply(V(new_g),function(x)ego_size(new_g,order=2,nodes=x))
#hist(cluster_sizes,main='',xlab='Cluster size')

#looking at the degree distribution
degreedistribution <- degree.distribution(new_g)*length(E(new_g))
#works out the average contacts per person
average_contacts <- sum(degreedistribution*c(1:length(degreedistribution)-1)/length(E(new_g)))
length(E(new_g))/length(V(new_g))*2

#lists the people that each individual has contact with in the model
contact_list <<- lapply(V(new_g),function(x) {cs <- as.vector(unlist(ego(new_g,order=1,nodes=x))); cs[cs!=x]})

#creates list of those who are contacts of the person's contacts
contact_of_contact_list <<- lapply(V(new_g),function(x) {
  cs <- as.vector(unlist(ego(new_g,order=1,nodes=x))); 
  cofcs <- as.vector(unlist(ego(new_g,order=2,nodes=x))); 
  ccs <- funique(c(cs,cofcs))
  ccs[ccs!=x]
})

## generate random edges network for random transmission
#random graph as according to the Erdos-Renyi model every edge created with the same constant probability in this caser the probabilityu is 10/length of vertices
random_g <- sample_gnp(length(V(new_g)), 10/length(V(new_g)))
#creates a list of all the vertices that they are attached to as a result of these random interactions
random_list <<- lapply(V(random_g),function(x) {cs <- as.vector(unlist(ego(random_g,order=1,nodes=x))); cs[cs!=x]})



                    

## adding in random connections to simulate social connections
social_g <<- new_g
social_g <- add_edges(social_g, edges=c(1,2))

for (i in 1:1731){
  listing_one <-sample(V(social_g),11000, replace = TRUE)
  listing_two <- sample(V(social_g), 11000, replace = TRUE)
  one_pair <-unlist(unname(as.list(listing_one[i])))
  two_pair<-unlist(unname(as.list(listing_two[i])))
  social_g <-add_edges(social_g, edges = c(one_pair, two_pair))
}


social_list <<- lapply(V(social_g),function(x) {cs <- as.vector(unlist(ego(social_g,order=1,nodes=x))); cs[cs!=x]})

  
  
plot(induced.subgraph(social_g,1:pts),layout=save_layout,
     vertex.label=NA,vertex.size=10,vertex.color=c('white',"#660066","#CC99FF")[demographic_index[1:pts]])
legend('topleft',legend=c('<19','19-65','>65'),pt.cex=2,col='black',pch=21, pt.bg=c('white',"#660066","#CC99FF"),bty='n')  
  
  







