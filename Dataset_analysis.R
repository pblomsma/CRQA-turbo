########################################################################################
#
# Dataset_analysis, calculates cross recurrence for multiple timeseries variables. Currenlty the methods are based on the iMap data,
# but hopefully they generalize to all kinds of dataframes.
#
# Author: Peter Blomsma
#
# This research has been funded by a grant No.: PROJ-007246 fromthe European Union, OP Zuid, the Ministry of
# Economic affairs,the Province of Noord-Brabant, and the municipality of Tilburg awarded to Max M. Louwerse.
########################################################################################

calculateRRforDataset <-  function(dataset, grouping_columns, leader_columns, follower_columns = leader_columns, delays = c(-80:80))
{
  require(dplyr)
  
  directory = "Exports/CSV_all"
  
  #NAs are set to zero - should be reported in the paper. 3729 values are NA
  dataset[is.na(dataset)] <- 0
  
  #Divide data in TS per role, session and dialog:
  dataset.groups <- dataset %>% group_by_at(grouping_columns) %>% group_data()
  
  for(leader.index in 4:length(leader_columns))
  {
    for(follower.index in 1:length(follower_columns))
    {
      RR.all_delays <-  calculateRRforVariables(dataset, dataset.groups, leader_columns[leader.index], follower_columns[follower.index], delays, spiketrains = FALSE, simulate_randomness = FALSE, turbo = TRUE)
      RR.all_delays$leader <- leader_columns[leader.index]
      RR.all_delays$follower <- follower_columns[follower.index]
  
      filename = paste(leader.index, follower.index, format(Sys.time(), "DT%Y-%m-%d_%H_%M"), ".csv", sep = "_")
      write.csv(RR.all_delays, paste(directory, filename, sep = "/"))  
    }
  }
  return(TRUE)
}

#Our solution
calculateRRforVariables <- function(
  dataset, 
  dataset.groups, 
  leader, 
  follower, 
  delays= c(-240:240), 
  simulate_randomness = FALSE, 
  only_onsets = FALSE, 
  spiketrains = TRUE,
  turbo = TRUE)
{
  print(paste(leader, follower, sep = " --?-> "))
  
  RR.all_delays <- NULL
  
  for(group.id in 1:nrow(dataset.groups))
  {
    group.rows <- unlist(dataset.groups[group.id,]$.rows)
    
    leader.ts <- dataset[group.rows,] %>% select(leader)
    follower.ts <- dataset[group.rows,] %>% select(follower)
    
    if(spiketrains)
    {
      leader.ts <- toSpikeTrain(leader.ts, only_onsets)
      follower.ts <- toSpikeTrain(follower.ts, only_onsets)
    }
    
    if(turbo)
    {
      RR.pergroup <- calculateRRProfile(leader.ts, follower.ts, delays)
    }
    else
    { 
      #original method
      require(crqa)
      leader.ts[which(leader.ts == 0),] <- -9
      follower.ts[which(follower.ts == 0),] <- -3
      
      ws <-  max(delays)
      RR.pergroup <- (drpdfromts(leader.ts, follower.ts, ws = max(delays), datatype = "categorical", radius = 0.000001))$profile
      
      if(length(RR.pergroup) == 1) #sometimes drpdfromts returns an empty profile.
      {
        RR.pergroup  <- data.frame(RR = replicate((ws * 2) + 1, 0))
      }
    }

    RR.pergroup <- data.frame(RR = RR.pergroup)
    RR.pergroup$group <- group.id #TODO: nice to have, add all info from dataset.groups[group.id,] to the result dataframe.
    RR.pergroup$delay <- delays
    RR.pergroup$RR_chancelevel <- unlist(lapply(delays, calculateRRChancelevel, timeseries1 = leader.ts, timeseries2 = follower.ts, simulations = 50, simulated = simulate_randomness))
    RR.pergroup <- cbind(RR.pergroup, dataset.groups[group.id, 1:(ncol(dataset.groups) - 1)])
    RR.all_delays <- rbind(RR.all_delays, RR.pergroup)
  }
  return(RR.all_delays)
}

