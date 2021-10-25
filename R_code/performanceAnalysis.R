source("arrayJobCombiner.R")
source("processName.R")

get_data_row <- function(fullData, search, n_start, n_end, n_interval, value_type, type = "mean"){

  data <- c()
  nVals <- c()
  #if(as.character(search) == "p(0)" & value_type == "cost"){print(paste0("Got ", nrow(fullData[fullData$search == search & fullData$size == 100, ]), " rows"))}
  for(n in seq(n_start, n_end, n_interval)){
    total <- 0
    for(disaster in unique(fullData[fullData$size == n & fullData$search == search, "disaster"])){
      #print(paste0("n ", n, " dis ", disaster, " input ", input))
      #print(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster & fullData$input_type == input, value_type])
      if(type == "mean"){total <- total + mean(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster, value_type])}
      if(type == "sd"){total <- total + sqrt(var(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster, value_type]))}
      if(type == "maxPlot"){total <- total + mean(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster, value_type]) + sqrt(var(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster, value_type]))}
    }
    # for(disaster in unique(fullData[fullData$size == n & fullData$search == search, "disaster"])){
    #   for(input in unique(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster, "input_type"])){
    #     #print(paste0("n ", n, " dis ", disaster, " input ", input))
    #     #print(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster & fullData$input_type == input, value_type])
    #     if(type == "mean"){total <- total + mean(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster & fullData$input_type == input, value_type])}
    #     if(type == "sd"){total <- total + sqrt(var(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster & fullData$input_type == input, value_type]))}
    #     if(type == "maxPlot"){total <- total + mean(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster & fullData$input_type == input, value_type]) + sqrt(var(fullData[fullData$size == n & fullData$search == search & fullData$disaster == disaster & fullData$input_type == input, value_type]))}
    #   }
    # }
    if(nrow(unique(fullData[fullData$size == n & fullData$search == search, ])) > 0){
      nVals <- c(nVals, n)
      data <- c(data, mean(total))
    }
    
  }
  #print(data)
  return(list(v1 = data, v2 = nVals))
}

remove_invalid_comparisons <- function(augmentationData){
  count <- 0
  for(size in unique(augmentationData$size)){
    for(seed in unique(augmentationData$seed)){
      for(input_type in unique(augmentationData$input_type)){
        for(disaster in unique(augmentationData$disaster)){
          if(nrow(augmentationData[augmentationData$size == size &
                                   augmentationData$seed == seed &
                                   augmentationData$input_type == input_type &
                                   augmentationData$disaster == disaster &
                                   (is.na(augmentationData$cost) |
                                   is.infinite(augmentationData$cost)), ]) == 0){next}
          count <- count + 1
          augmentationData <- augmentationData[augmentationData$size != size |
                                                 augmentationData$seed != seed |
                                                 augmentationData$input_type != input_type |
                                                 augmentationData$disaster != disaster, ]
        }
      }
    }
  }
  #print(paste0("Found ", count, " to remove"))
  return(augmentationData)
}




augmentation_analysis <- function(preFilename, augmentationData, input_type, disaster, n_start, n_end, n_interval, searches, plotMult = 1.25){

  # Only look at type of input if want to
  if(input_type != "ALL"){
    augmentationData <- augmentationData[augmentationData$input_type == input_type, ]
  }
  # Only look at specific disaster length if want to
  if(disaster != "ALL"){
    augmentationData <- augmentationData[augmentationData$disaster == disaster, ]
  }
  #print(paste0("Start with ", nrow(augmentationData), "rows, got ", nrow(augmentationData[is.na(augmentationData$cost), ]) + nrow(augmentationData[is.infinite(augmentationData$cost), ]), " should remove"))
  #augmentationData <- remove_invalid_comparisons(augmentationData)
  #print(paste0("End with ", nrow(augmentationData), "rows, got ", nrow(augmentationData[is.na(augmentationData$cost), ]) + nrow(augmentationData[is.infinite(augmentationData$cost), ]), " should remove"))
  
  # Note no check for na's as there should not be any by this stage
  augmentationData <- augmentationData[!is.na(augmentationData$cost), ]
  augmentationData <- augmentationData[!is.infinite(augmentationData$cost), ]
  #print(paste0("Second end with ", nrow(augmentationData), "rows, got ", nrow(augmentationData[is.na(augmentationData$cost), ]) + nrow(augmentationData[is.infinite(augmentationData$cost), ]), " should remove"))
  # Got all of the values for the functions, can now plot them
  # Start with the cost 
  x <- seq(n_start, n_end, n_interval)
  
  if(disaster == "ALL"){
    disaster_filename <- "disasterAll"
  }else{
    disaster_filename <- paste0("disaster", disaster)
  }
  
  plot_filename_cost <- paste0("../data/analysis_plots/cost/", preFilename, "-n",n_start,"-",n_interval,"-",n_end, "_input", input_type, "_", disaster_filename, "_augmentation_cost.png")
  range <- c()
  # Get min and max
  for (search in searches){
    range <- c(range, get_data_row(augmentationData, search, n_start, n_end, n_interval, "cost")[[1]])
  }
  range <- c(0.9 * min(range), 1.3 *  max(range))
  
  cex_expansion = 1.5
  cex_points = 2
  
  png(plot_filename_cost, width = plotMult*682, height = plotMult*512)
  plot(get_data_row(augmentationData, searches[[1]], n_start, n_end, n_interval, "cost", "mean")[[2]], get_data_row(augmentationData, searches[[1]], n_start, n_end, n_interval, "cost", "mean")[[1]], 
       #col = "black",
       col = processColour(searches[[1]]),
       type = "b",
       ylab = "Cost (length)",
       xlab = "Graph size",
       ylim = range,
       xaxt = "n",
       #pch = 1,
       pch = processPCH(searches[[1]]),
       cex.lab = cex_expansion,
       cex.axis = cex_expansion,
       cex.main = cex_expansion,
       cex.sub = cex_expansion,
       cex = cex_points)
  pch <- 2
  for(search in searches[-1]){
    #points(get_data_row(augmentationData, search, n_start, n_end, n_interval, "cost", "mean")[[2]], get_data_row(augmentationData, search, n_start, n_end, n_interval, "cost", "mean")[[1]], col = "black", type = "b", pch = pch, cex = cex_points)
    points(get_data_row(augmentationData, search, n_start, n_end, n_interval, "cost", "mean")[[2]], get_data_row(augmentationData, search, n_start, n_end, n_interval, "cost", "mean")[[1]], col = processColour(search), type = "b", pch = processPCH(search), cex = cex_points)
    pch <- pch + 1
  }
  #legend("topleft", legend = processNames(as.character(searches)), pch = 1:length(searches), cex = cex_expansion)
  legend("topleft", legend = processNames(as.character(searches)), pch = processPCHs(searches), col = processColours(searches), cex = cex_expansion)
  axis(1, at = x, cex.axis = cex_expansion)
  out <- dev.off()
  

  plot_filename_time <- paste0("../data/analysis_plots/time/", preFilename, "-n",n_start,"-",n_interval,"-",n_end, "_input", input_type, "_", disaster_filename, "_augmentation_time.png")

  range <- c()
  # Get min and max
  for (search in searches){
    range <- c(range, get_data_row(augmentationData, search, n_start, n_end, n_interval, "augment_time")[[1]])
  }
  range <- c(0.9 * min(range), 1.1 *  max(range))

  png(plot_filename_time, width = plotMult*682, height = plotMult*512)
  plot(get_data_row(augmentationData, searches[[1]], n_start, n_end, n_interval, "augment_time", "mean")[[2]], get_data_row(augmentationData, searches[[1]], n_start, n_end, n_interval, "augment_time", "mean")[[1]],
       col = processColour(searches[[1]]),
       type = "b",
       ylab = "Time (seconds)",
       xlab = "Graph size",
       ylim = range,
       xaxt = "n",
       pch = processPCH(searches[[1]]),
       cex.lab = cex_expansion,
       cex.axis = cex_expansion,
       cex.main = cex_expansion,
       cex.sub = cex_expansion,
       cex = cex_points)
  pch <- 2
  for(search in searches[-1]){
    #points(get_data_row(augmentationData, search, n_start, n_end, n_interval, "augment_time", "mean")[[2]], get_data_row(augmentationData, search, n_start, n_end, n_interval, "augment_time", "mean")[[1]], col = "black", type = "b", pch = pch, cex = cex_points)
    points(get_data_row(augmentationData, search, n_start, n_end, n_interval, "augment_time", "mean")[[2]], get_data_row(augmentationData, search, n_start, n_end, n_interval, "augment_time", "mean")[[1]], col = processColour(search), type = "b", pch = processPCH(search), cex = cex_points)
    pch <- pch + 1
  }
  #legend("topleft", legend = processNames(as.character(searches)), pch = 1:length(searches), cex = cex_expansion)
  legend("topleft", legend = processNames(as.character(searches)), pch = processPCHs(searches), cex = cex_expansion, col = processColours(searches))
  #print(as.character(t(searches$name)))
  axis(1, at = x, cex.axis = cex_expansion)
  out <- dev.off()
  #legend_names <- c()
  #for (row in 1:nrow(searches)){
  #  legend_names <- c(legend_names, searches[row, "name"])
  #}
  #print(legend_names)
  #legend("topleft", legend_names, pch = 1:nrow(searches))
  if("structure_time" %in% colnames(augmentationData)){
    plot_filename_time <- paste0("../data/analysis_plots/time/", preFilename, "-structBuildTimes-n",n_start,"-",n_interval,"-",n_end, "_input", input_type, "_", disaster_filename, "_augmentation_time.png")
    range <- c()
    # Get min and max
    for (search in searches){
      range <- c(range, get_data_row(augmentationData, search, n_start, n_end, n_interval, "structure_time")[[1]])
    }
    range <- c(0.9 * min(range), 1.1 *  max(range))
    png(plot_filename_time, width = plotMult*682, height = plotMult*512)
    plot(get_data_row(augmentationData, searches[[1]], n_start, n_end, n_interval, "structure_time", "mean")[[2]], get_data_row(augmentationData, searches[[1]], n_start, n_end, n_interval, "structure_time", "mean")[[1]],
         col = "black",
         type = "b",
         ylab = "Time (seconds)",
         xlab = "Graph size",
         ylim = range,
         xaxt = "n",
         pch = 1,
         cex.lab = cex_expansion,
         cex.axis = cex_expansion,
         cex.main = cex_expansion,
         cex.sub = cex_expansion,
         cex = cex_points)
    pch <- 2
    for(search in searches[-1]){
      points(get_data_row(augmentationData, search, n_start, n_end, n_interval, "structure_time", "mean")[[2]], get_data_row(augmentationData, search, n_start, n_end, n_interval, "structure_time", "mean")[[1]], col = "black", type = "b", pch = pch, cex = cex_points)
      pch <- pch + 1
    }
    legend("topleft", legend = processNames(as.character(searches)), pch = 1:length(searches), cex = cex_expansion)
    axis(1, at = x, cex.axis = cex_expansion)
    out <- dev.off()
  }
}

full_augmentation_analysis <- function(preFilename, augmentationData, n_start, n_end, n_interval, searches, mult = 1.0){
  augmentation_analysis(preFilename, augmentationData, "ALL", "ALL", n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "ALL", 0.001, n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "ALL", 5, n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "ALL", 100, n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "ALL", 300, n_start, n_end, n_interval, searches, mult)
  
  augmentation_analysis(preFilename, augmentationData, "MST", "ALL", n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "MST", 0.001, n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "MST", 5, n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "MST", 100, n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "MST", 300, n_start, n_end, n_interval, searches, mult)
  
  augmentation_analysis(preFilename, augmentationData, "TopBottom", "ALL", n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "TopBottom", 0.001, n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "TopBottom", 5, n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "TopBottom", 100, n_start, n_end, n_interval, searches, mult)
  augmentation_analysis(preFilename, augmentationData, "TopBottom", 300, n_start, n_end, n_interval, searches, mult)
}



