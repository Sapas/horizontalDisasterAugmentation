source("arrayJobCombiner.R")
source("processName.R")

get_data_row <- function(fullData, search, n_start, n_end, n_interval, value_type, type = "mean"){
  data <- c()
  for(n in seq(n_start, n_end, n_interval)){
    if(type == "mean"){data <- c(data, mean(fullData[fullData$size == n & fullData$search == search, value_type]))}
    if(type == "sd"){data <- c(data, sqrt(var(fullData[fullData$size == n & fullData$search == search, value_type])))}
    if(type == "maxPlot"){data <- c(data, mean(fullData[fullData$size == n & fullData$search == search, value_type]) + sqrt(var(fullData[fullData$size == n, value_type])))}
  }
  return(data)
}




augmentation_analysis <- function(preFilename, augmentationData, input_type, disaster, n_start, n_end, n_interval, searches){

  # Only look at type of input if want to
  if(input_type != "ALL"){
    augmentationData <- augmentationData[augmentationData$input_type == input_type, ]
  }
  # Only look at specific disaster length if want to
  if(disaster != "ALL"){
    augmentationData <- augmentationData[augmentationData$disaster == disaster, ]
  }
  # Note no check for na's as there should not be any by this stage
  augmentationData <- augmentationData[!is.na(augmentationData$cost), ]
  augmentationData <- augmentationData[!is.infinite(augmentationData$cost), ]
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
    range <- c(range, get_data_row(augmentationData, search, n_start, n_end, n_interval, "cost"))
  }
  range <- c(0.9 * min(range), 1.3 *  max(range))
  
  cex_expansion = 1.5
  cex_points = 2
  
  #png(plot_filename_cost, width = 682, height = 512)
  png(plot_filename_cost, width = 2*682, height = 2*512)
  plot(x, get_data_row(augmentationData, searches[[1]], n_start, n_end, n_interval, "cost", "mean"), 
       col = "black", 
       type = "b",
       ylab = "Cost (length)",
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
    points(x, get_data_row(augmentationData, search, n_start, n_end, n_interval, "cost", "mean"), col = "black", type = "b", pch = pch, cex = cex_points)
    pch <- pch + 1
  }
  legend("topleft", legend = processNames(as.character(searches)), pch = 1:length(searches), cex = cex_expansion)
  axis(1, at = x, cex.axis = cex_expansion)
  out <- dev.off()




  plot_filename_time <- paste0("../data/analysis_plots/time/", preFilename, "-n",n_start,"-",n_interval,"-",n_end, "_input", input_type, "_", disaster_filename, "_augmentation_time.png")

  range <- c()
  # Get min and max
  for (search in searches){
    range <- c(range, get_data_row(augmentationData, search, n_start, n_end, n_interval, "augment_time"))
  }
  range <- c(0.9 * min(range), 1.1 *  max(range))

  #png(plot_filename_time, width = 682, height = 512)
  png(plot_filename_time, width = 2*682, height = 2*512)
  plot(x, get_data_row(augmentationData, searches[[1]], n_start, n_end, n_interval, "augment_time", "mean"),
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
    points(x, get_data_row(augmentationData, search, n_start, n_end, n_interval, "augment_time", "mean"), col = "black", type = "b", pch = pch, cex = cex_points)
    pch <- pch + 1
  }
  legend("topleft", legend = processNames(as.character(searches)), pch = 1:length(searches), cex = cex_expansion)
  #print(as.character(t(searches$name)))
  axis(1, at = x, cex.axis = cex_expansion)
  out <- dev.off()
  #legend_names <- c()
  #for (row in 1:nrow(searches)){
  #  legend_names <- c(legend_names, searches[row, "name"])
  #}
  #print(legend_names)
  #legend("topleft", legend_names, pch = 1:nrow(searches))
}

#options(max.print=1000000)
#combined <- combineArrayResults("preliminaryRun", 1, 880)
#average <- getAveragePerformance("preliminaryRun", 1, 880)

#augmentation_analysis("preRun", combined, "ALL", "ALL", 10, 100, 10, as.character(unique(combined$search)))
#augmentation_analysis("preRun", combined, "ALL", "ALL", 10, 100, 10, c("a(0)", "b(0)", "o(2)", "p(0)", "s(0)"))
#augmentation_analysis("preRun", combined, "ALL", 0.001, 10, 100, 10, c("a(0)", "b(0)", "o(2)", "p(0)", "s(0)"))

combinedPre <- combineArrayResults("preliminaryRun", 1, 1040)
averagePre <- getAveragePerformance("preliminaryRun", 1, 1040)

combinedMain <- combineArrayResults("mainRun", 1, 960)
averageMain <- getAveragePerformance("mainRun", 1, 960)

combinedPoly <- combineArrayResults("polynomialRun", 1, 640)
averagePoly <- getAveragePerformance("polynomialRun", 1, 640)

combinedMainAndPoly <- rbind(combinedMain, combinedPoly)

augmentation_analysis("preRun", combinedPre, "ALL", "ALL", 10, 100, 10, c("a(0)", "b(0)", "o(2)", "p(0)", "s(0)"))
augmentation_analysis("preRun", combinedPre, "ALL", 0.001, 10, 100, 10, c("a(0)", "b(0)", "o(2)", "p(0)", "s(0)"))
augmentation_analysis("mainRunInitial", combinedMain, "ALL", "ALL", 10, 100, 10, c("o(2)", "q(2)", "u(2)"))
augmentation_analysis("mainRunInitial", combinedMain, "MST", 300, 10, 100, 10, c("o(2)", "q(2)", "u(2)"))
augmentation_analysis("mainRunWeightWeakConnection", combinedMain, "ALL", "ALL", 10, 100, 10, c("o(1)", "o(2)", "o(3)", "o(4)"))
augmentation_analysis("mainRunWeightWeakConnectionSingleOrigin", combinedMain, "ALL", "ALL", 10, 100, 10, c("q(1)", "q(2)", "q(3)", "q(4)"))
augmentation_analysis("mainRunWeightWeakConnectionStochastic", combinedMain, "ALL", "ALL", 10, 100, 10, c("u(1)", "u(2)", "u(3)", "u(4)"))
augmentation_analysis("mainRunFinalPresentation", combinedMain, "ALL", "ALL", 10, 100, 10, c("o(3)", "q(2)"))
augmentation_analysis("polyComp", combinedMainAndPoly, "ALL", "ALL", 10, 100, 10, c("o(3)", "q(2)", "v(3)", "w(2)"))
augmentation_analysis("polyComp", combinedMainAndPoly, "ALL", 0.001, 10, 100, 10, c("o(3)", "q(2)", "v(3)", "w(2)"))
augmentation_analysis("polyComp", combinedMainAndPoly, "ALL", 300, 10, 100, 10, c("o(3)", "q(2)", "v(3)", "w(2)"))





