########################################################################################
#
# Methods to cross-recurrence analyse discrete timeseries.
#
# Author: Peter Blomsma
#
# This research has been funded by a grant No.: PROJ-007246 fromthe European Union, OP Zuid, the Ministry of
# Economic affairs,the Province of Noord-Brabant, and the municipality of Tilburgawarded to Max M. Louwerse.
########################################################################################

#Main method, returns the RR for every delay.
calculateRRProfile <-  function(timeseries1, timeseries2, delays = c(1:240))
{
  result <- unlist(lapply(delays, calculateRR, timeseries1, timeseries2))
  names(result) <- delays
  return(result)
}

###########################################################################
# RR_chancelevel / Baseline related function
###########################################################################
calculateRR <- function(delay = 0, timeseries1, timeseries2)
{
  diagonal <- calculateDiagonalLine(timeseries1, timeseries2, delay = delay)
  result <- sum(diagonal)/length(diagonal)
  
  if(is.nan(result))
  {
    result <- 0
  }
  return(result)
}

calculateDiagonalLine <-  function(timeseries1, timeseries2, delay = 0)
{
  # Compares corresponding elements of the timeseries with eachtother. In other words it's constructing a diagonal line
  # of a recurrence plot. On the line, every 1 indicates a match, every 0 a non-match. 
  # 
  # Args: 
  #   timeseries1, is compared from begin to the (end - delay). 
  #   timeseries2, is compared from delay to the end.
  #   delay: delay in element size (i.e. in timeunits).
  #
  # Return:
  # Vector with length of spiketrain - delay, indicating which corresponding elements both contains a spike.
  
  if(is.list(timeseries1))
  {
    timeseries1 <- as.array(timeseries1[,1])
  }
  
  if(is.list(timeseries2))
  {
    timeseries2 <- as.array(timeseries2[,1])
  }
  
  if(abs(delay) >= min(length(timeseries1), length(timeseries2)))
  {
    warning("Delay in calculateDiagonalLine is longer or equal to the length of timeseries, returning empty diagonal")
    return(as.numeric(vector()))
  }
  
  if(delay < 0)
  {
    # A minus delay results in the same calculation as above, however with parameters switched and a positive delay.
    return(calculateDiagonalLine(timeseries2, timeseries1, abs(delay)))
  }
  
  timeseries1 <- timeseries1[(1+delay):length(timeseries1)] 
  timeseries2 <- timeseries2[1:(length(timeseries2) - delay)] 
  
  return(as.numeric(ifelse(timeseries1 == timeseries2, ifelse(timeseries1 == 1, 1, 0), 0)))
}

toSpikeTrain <-  function(timeseries.original, only_onsets = FALSE)
{
  # Identifies at which timepoints the timeseries changed value.
  #
  # Args:
  #  timeseries: timeseries that is encoded in 0s and 1s.
  #  only_onsets: if TRUE it only identifies which timepoints in the series change from 0 to 1 (and not the other way arround)
  #
  #Returns:
  # Timeseries (spiketrain) with the length of the timeseries.original, each elements is either 
  # 0 (= no change happened at that timepoint) or 
  # 1 (= the value of the series changed at that timepoints).
  #
  
  
  if(!is.double(timeseries.original))
  {
    timeseries.original <- as.array(timeseries.original[,1])
  }
  
  if(length(timeseries.original) < 2)
  {
    return(timeseries.original)
  }
    
  change <- c(timeseries.original[1], as.vector(timeseries.original[1:(length(timeseries.original)-1)]))
  
  if(only_onsets)
  {
    return(as.numeric(ifelse(timeseries.original == 1 & change == 0, 1, 0)))
  }
  
  return(as.numeric(ifelse(timeseries.original == change, 0, 1)))
}

calculateRRChancelevel <-  function(delay, timeseries1, timeseries2, simulated, simulations = length(timeseries2))
{
  if(simulated)
  {
    return(mean(replicate(simulations, calculateRR(delay = delay, timeseries1 = sample(timeseries1), timeseries2 = sample(timeseries2)))))
  }
  timeseries1_length <- max(length(timeseries1), nrow(timeseries1))
  timeseries2_length <- max(length(timeseries2), nrow(timeseries2))
  
  return(sum(timeseries1) * (sum(timeseries2)/timeseries1_length) / timeseries2_length)
}