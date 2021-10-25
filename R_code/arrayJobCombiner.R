combineArrayResults <- function(runName, arrayStart, arrayEnd){
  combinedData <- read.table(paste0("../data/runResults/", runName, "_arrayJob", arrayStart, ".txt"), header = TRUE, sep = " ", fill = TRUE)
  incomplete <- ""
  for(i in (arrayStart + 1):arrayEnd){
    filename <- paste0("../data/runResults/", runName, "_arrayJob", i, ".txt")
    if(!file.exists(filename)){
      print(paste0("Skipping file ", filename, " as it does not exist!"))
      next
    }
    newData <- read.table(filename, header = TRUE, sep = " ", fill = TRUE)
    # if(((newData[1, "search"] != "a(0)" | newData[1, "size"] <= 50) & nrow(newData) < 25) | nrow(newData) < 5){
    # #if(((substr(newData[1, "search"], 1, 1) != "u" | newData[1, "size"] <= 50) & nrow(newData) < 100) | nrow(newData) < 20){
    #  print(paste0("File ", filename, " contains unexpected number of rows, only contains ", nrow(newData)))
    #  incomplete <- paste0(incomplete, ",", i)
    # }
    #if(nrow(newData[is.na(newData$cost), ]) > 0){
    #  print(paste0("File ", filename, " contains na values"))
    #  print(newData[is.na(newData$cost), ])
    #}
    combinedData <- rbind(combinedData, newData)
  }
  # print("Incomplete files: ")
  # print(incomplete)
  return(combinedData)
}

averagePerformance <- function(data){
  # Now get averages
  size <- unique(data$size)
  input_type <- unique(data$input_type)
  disaster <- unique(data$disaster)
  search <- unique(data$search)
  xDim <- unique(data$xDim)
  yDim <- unique(data$yDim)
  extraTimeVals <- ("structure_time" %in% colnames(data))
  
  if(extraTimeVals){averagedData <- data.frame(matrix(ncol = 12, nrow = 0))}
  else{averagedData <- data.frame(matrix(ncol = 10, nrow = 0))}
  entries <- 0
  for(size in unique(data$size)){
    for(input_type in unique(data$input_type)){
      for(disaster in unique(data$disaster)){
        for(search in unique(data$search)){
          for(xDim in unique(data$xDim)){
            for(yDim in unique(data$yDim)){
              
              # First check this actually exists
              if(length(data[data$size == size &
                             data$input_type == input_type & 
                             data$disaster == disaster & 
                             data$search == search & 
                             data$xDim == xDim & 
                             data$yDim == yDim &
                             !is.na(data$cost), c("cost")]) == 0){next}
              
              # All good, can add a rows
              entries <- entries + 1
              if(extraTimeVals){
                averagedData[entries, ] <- c(size, input_type, disaster, search, xDim, yDim, 
                                             mean(data[data$size == size &
                                                         data$input_type == input_type & 
                                                         data$disaster == disaster & 
                                                         data$search == search & 
                                                         data$xDim == xDim & 
                                                         data$yDim == yDim &
                                                         !is.na(data$cost), c("cost")]),
                                             sqrt(var(data[data$size == size &
                                                             data$input_type == input_type & 
                                                             data$disaster == disaster & 
                                                             data$search == search & 
                                                             data$xDim == xDim & 
                                                             data$yDim == yDim &
                                                             !is.na(data$cost), c("cost")])),
                                             mean(data[data$size == size &
                                                         data$input_type == input_type & 
                                                         data$disaster == disaster & 
                                                         data$search == search & 
                                                         data$xDim == xDim & 
                                                         data$yDim == yDim &
                                                         !is.na(data$cost), c("augment_time")]),
                                             sqrt(var(data[data$size == size &
                                                             data$input_type == input_type & 
                                                             data$disaster == disaster & 
                                                             data$search == search & 
                                                             data$xDim == xDim & 
                                                             data$yDim == yDim &
                                                             !is.na(data$cost), c("augment_time")])),
                                             mean(data[data$size == size &
                                                         data$input_type == input_type & 
                                                         data$disaster == disaster & 
                                                         data$search == search & 
                                                         data$xDim == xDim & 
                                                         data$yDim == yDim &
                                                         !is.na(data$cost), c("structure_time")]),
                                             mean(data[data$size == size &
                                                         data$input_type == input_type & 
                                                         data$disaster == disaster & 
                                                         data$search == search & 
                                                         data$xDim == xDim & 
                                                         data$yDim == yDim &
                                                         !is.na(data$cost), c("edge_time")]))
              
              }else{
                averagedData[entries, ] <- c(size, input_type, disaster, search, xDim, yDim, 
                                           mean(data[data$size == size &
                                                       data$input_type == input_type & 
                                                       data$disaster == disaster & 
                                                       data$search == search & 
                                                       data$xDim == xDim & 
                                                       data$yDim == yDim &
                                                       !is.na(data$cost), c("cost")]),
                                            sqrt(var(data[data$size == size &
                                                            data$input_type == input_type & 
                                                            data$disaster == disaster & 
                                                            data$search == search & 
                                                            data$xDim == xDim & 
                                                            data$yDim == yDim &
                                                            !is.na(data$cost), c("cost")])),
                                           mean(data[data$size == size &
                                                       data$input_type == input_type & 
                                                       data$disaster == disaster & 
                                                       data$search == search & 
                                                       data$xDim == xDim & 
                                                       data$yDim == yDim &
                                                       !is.na(data$cost), c("augment_time")]),
                                           sqrt(var(data[data$size == size &
                                                           data$input_type == input_type & 
                                                           data$disaster == disaster & 
                                                           data$search == search & 
                                                           data$xDim == xDim & 
                                                           data$yDim == yDim &
                                                           !is.na(data$cost), c("augment_time")])))
              }
            }
          }
        }
      }
    }
  }
  
  if(extraTimeVals){colnames(averagedData) <- c("size", "input_type", "disaster", "search", "xDim", "yDim", "costMean", "costSD", "timeMean", "timeSD", "timeStructMean", "timeEdgeMean")}
  else{colnames(averagedData) <- c("size", "input_type", "disaster", "search", "xDim", "yDim", "costMean", "costSD", "timeMean", "timeSD")}
  return(averagedData)
}

getAveragePerformance <- function(runName, arrayStart, arrayEnd){
  data <- combineArrayResults(runName, arrayStart, arrayEnd)
  if(nrow(data[is.na(data$cost), ]) > 0){
    print("Some of the data looped forever, here it is:")
    print(data[is.na(data$cost), ])
    print("These are the searches that looped forever: ")
    print(unique(data[is.na(data$cost), "search"]))
  }
  if(nrow(data[is.infinite(data$cost), ]) > 0){
    print("Some of the data did not get solved, here it is:")
    print(data[is.infinite(data$cost), ])
    print("These are the searches that did not work sometimes: ")
    print(unique(data[is.infinite(data$cost), "search"]))
    for(search in unique(data[is.infinite(data$cost), "search"])){
      print(paste0("Search ", search, " infeasible ", nrow(data[is.infinite(data$cost) & data$search == search, ])))
    }
  }
  
  averagePerformance(data)
}

addComparisonColumns <- function(runData){
  
  
  for(size in unique(runData$size)){
    for(input_type in unique(runData$input_type)){
      for(disaster in unique(runData$disaster)){
        for(xDim in unique(runData$xDim)){
          for(yDim in unique(runData$yDim)){
            if(length(runData[runData$size == size &
                              runData$input_type == input_type & 
                              runData$disaster == disaster &
                              runData$xDim == xDim & 
                              runData$yDim == yDim &
                              !is.na(runData$cost) & 
                              !is.infinite((runData$cost)), c("cost")]) == 0){next}
            # So far so good, store this
            runData[runData$size == size &
                      runData$input_type == input_type & 
                      runData$disaster == disaster &
                      runData$xDim == xDim & 
                      runData$yDim == yDim, c("bestCost")] <- min(runData[runData$size == size &
                                                                            runData$input_type == input_type & 
                                                                            runData$disaster == disaster &
                                                                            runData$xDim == xDim & 
                                                                            runData$yDim == yDim &
                                                                            !is.na(runData$cost) & 
                                                                            !is.infinite((runData$cost)), c("cost")])
            
            runData[runData$size == size &
                      runData$input_type == input_type & 
                      runData$disaster == disaster &
                      runData$xDim == xDim & 
                      runData$yDim == yDim, c("bestTime")] <- min(runData[runData$size == size &
                                                                            runData$input_type == input_type & 
                                                                            runData$disaster == disaster &
                                                                            runData$xDim == xDim & 
                                                                            runData$yDim == yDim &
                                                                            !is.na(runData$cost) & 
                                                                            !is.infinite((runData$cost)), c("augment_time")])
            
            # Now can do comparison
            for(search in unique(runData$search)){
              if(length(runData[runData$size == size &
                                runData$input_type == input_type & 
                                runData$disaster == disaster &
                                runData$xDim == xDim & 
                                runData$yDim == yDim &
                                runData$search == search &
                                !is.na(runData$cost) & 
                                !is.infinite((runData$cost)), c("cost")]) == 0){next}
              
              
              runData[runData$size == size &
                        runData$input_type == input_type & 
                        runData$disaster == disaster &
                        runData$search == search &
                        runData$xDim == xDim & 
                        runData$yDim == yDim, c("costRatio")] <- runData[runData$size == size &
                                                                           runData$input_type == input_type & 
                                                                           runData$disaster == disaster &
                                                                           runData$search == search &
                                                                           runData$xDim == xDim & 
                                                                           runData$yDim == yDim, c("cost")] / runData[runData$size == size &
                                                                                                                        runData$input_type == input_type & 
                                                                                                                        runData$disaster == disaster &
                                                                                                                        runData$search == search &
                                                                                                                        runData$xDim == xDim & 
                                                                                                                        runData$yDim == yDim, c("bestCost")]
              
              
              runData[runData$size == size &
                        runData$input_type == input_type & 
                        runData$disaster == disaster &
                        runData$search == search &
                        runData$xDim == xDim & 
                        runData$yDim == yDim, c("timeRatio")] <- runData[runData$size == size &
                                                                           runData$input_type == input_type & 
                                                                           runData$disaster == disaster &
                                                                           runData$search == search &
                                                                           runData$xDim == xDim & 
                                                                           runData$yDim == yDim, c("augment_time")] / runData[runData$size == size &
                                                                                                                        runData$input_type == input_type & 
                                                                                                                        runData$disaster == disaster &
                                                                                                                        runData$search == search &
                                                                                                                        runData$xDim == xDim & 
                                                                                                                        runData$yDim == yDim, c("bestTime")]
            }
          }
        }
      }
    }
  }
  return(runData)
}


# Example usage
#combined <- combineArrayResults("smallTest", 1, 48)
#average <- getAveragePerformance("smallTest", 1, 48)
