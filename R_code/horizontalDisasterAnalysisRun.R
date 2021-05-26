get_data_row <- function(data_frame, search, weight, n_start, n_end, n_interval, value_type){
  data_frame <- data_frame[as.character(data_frame$search) == search, ]
  data_frame <- data_frame[data_frame$search_weight == weight, ]
  
  data <- c()
  for(n in seq(n_start, n_end, n_interval)){
    
    data <- c(data, mean(data_frame[data_frame$size == n, value_type]))
  }
  return(data)
}

augmentation_analysis <- function(filename, input_type, disaster_length, n_start, n_end, n_interval, searches){
	# Get the data and store it
	augmentation <- read.table(filename, header = TRUE)
	# Only look at type of input if want to
	if(input_type != 'ALL'){
		augmentation <- augmentation[augmentation$input_type == input_type, ]
	}
	# Only look at specific disaster length if want to
	if(disaster_length != -1){
		augmentation <- augmentation[augmentation$disaster == disaster_length, ]
	}
	# Quick check to see if any of the remaining ones did not get solved
	augmentation <- augmentation[!is.na(augmentation$cost), ]
	
	# Got all of the values for the functions, can now plot them
	# Start with the cost 
	x <- seq(n_start, n_end, n_interval)
	
	if(disaster_length == -1){
		disaster_filename <- "disasterAll"
	}else{
		disaster_filename <- paste0("disaster", disaster_length)
	}
	plot_filename_cost <- paste0("data/analysis_plots/cost/n",n_start,"-",n_interval,"-",n_end, "_input", input_type, "_", disaster_filename, "_augmentation_cost.png")
	range <- c()
	# Get min and max
	for (row in 1:nrow(searches)){
	  range <- c(range, get_data_row(augmentation, searches[row, "search"], searches[row, "weight"], n_start, n_end, n_interval, "cost"))
	  # print(get_data_row(augmentation, searches[row, "search"], searches[row, "weight"], n_start, n_end, n_interval, "cost"))
	  # print(100*(get_data_row(augmentation, searches[row, "search"], searches[row, "weight"], n_start, n_end, n_interval, "cost") -
	  #         get_data_row(augmentation, 'o', 3, n_start, n_end, n_interval, "cost") ) /
	  #         get_data_row(augmentation, searches[row, "search"], searches[row, "weight"], n_start, n_end, n_interval, "cost"))
	}
	range <- c(0.9 * min(range), 1.3 *  max(range))
	
	cex_expansion = 1.5
	cex_points = 2
	
	png(plot_filename_cost, width = 682, height = 512)
	plot(x, get_data_row(augmentation, searches[1, "search"], searches[1, "weight"], n_start, n_end, n_interval, "cost"), 
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
	for (row in 2:nrow(searches)){
	  points(x, get_data_row(augmentation, searches[row, "search"], searches[row, "weight"], n_start, n_end, n_interval, "cost"), col = "black", type = "b", pch = row, cex = cex_points)
	}
	legend("topleft", legend = as.character(t(searches$name)), pch = 1:nrow(searches), cex = cex_expansion)
	axis(1, at = x, cex.axis = cex_expansion)
	out <- dev.off()
	#legend_names <- c()
	#for (row in 1:nrow(searches)){
	#  legend_names <- c(legend_names, searches[row, "name"])
	#}
	#print("Here")
	#print(legend_names)
	#print(searches$name)
	#legend("topleft", searches$name, pch = 1:nrow(searches))

	
	
	
	plot_filename_time <- paste0("data/analysis_plots/time/n",n_start,"-",n_interval,"-",n_end, "_input", input_type, "_", disaster_filename, "_augmentation_time_clean.png")
	
	range <- c()
	# Get min and max
	for (row in 1:nrow(searches)){
	  range <- c(range, get_data_row(augmentation, searches[row, "search"], searches[row, "weight"], n_start, n_end, n_interval, "augment_time"))
	  # print(get_data_row(augmentation, searches[row, "search"], searches[row, "weight"], n_start, n_end, n_interval, "augment_time"))
	  # print(100*(get_data_row(augmentation, searches[row, "search"], searches[row, "weight"], n_start, n_end, n_interval, "augment_time") -
	  #              get_data_row(augmentation, 'q', 2, n_start, n_end, n_interval, "augment_time") ) /
	  #         get_data_row(augmentation, searches[row, "search"], searches[row, "weight"], n_start, n_end, n_interval, "augment_time"))
	}
	range <- c(0.9 * min(range), 1.1 *  max(range))
	
	png(plot_filename_time, width = 682, height = 512)
	plot(x, get_data_row(augmentation, searches[1, "search"], searches[1, "weight"], n_start, n_end, n_interval, "augment_time"), 
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
	for (row in 2:nrow(searches)){
	  points(x, get_data_row(augmentation, searches[row, "search"], searches[row, "weight"], n_start, n_end, n_interval, "augment_time"), col = "black", type = "b", pch = row, cex = cex_points)
	}
	legend("topleft", legend = as.character(t(searches$name)), pch = 1:nrow(searches), cex = cex_expansion)
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


# The next section defines which schemes are plotted. Uncomment the ones which are desired

# PRE RUN
# filename <- "data/exec_times/pre_run.txt"
# n_start <- 10
# n_end <- 50
# n_interval <- 10

# searches <-data.frame(search=character(), weight=double(), name=character())
# searches<-rbind(searches, data.frame(search="a", weight=0, name="BlockJoin-StrongConnection"))
# searches<-rbind(searches, data.frame(search="b", weight=0, name="BlockJoin-StrongConnection-CutDistance"))
# # searches<-rbind(searches, data.frame(search="c", weight=0, name="LeafOrigin-StrongConnection"))
# # searches<-rbind(searches, data.frame(search="d", weight=0, name="LeafOrigin-StrongConnection-CutDistance"))
# # searches<-rbind(searches, data.frame(search="e", weight=0, name="LeafOrigin-WeakConnection"))
# # searches<-rbind(searches, data.frame(search="i", weight=0, name="LeafJoin-WeakConnection"))
# # searches<-rbind(searches, data.frame(search="j", weight=0, name="LeafJoin-StrongConnection"))
# # searches<-rbind(searches, data.frame(search="k", weight=0, name="LeafJoin-StrongConnection-CutDistance"))
# searches<-rbind(searches, data.frame(search="o", weight=0, name="LeafCombination(2)-WeakConnection"))
# searches<-rbind(searches, data.frame(search="p", weight=0, name="LeafCombination(2)-StrongConnection"))
# searches<-rbind(searches, data.frame(search="s", weight=0, name="LeafCombination(2)-StrongConnection-CutDistance"))


# EXPERIMENTAL RUN
n_start <- 50
n_end <- 100
n_interval <- 10

filename <- "data/exec_times/ExperimentalRun.txt"
searches <-data.frame(search=character(), weight=double(), name=character())
# searches<-rbind(searches, data.frame(search="o", weight=1, name="LeafCombination(1)-WeakConnection"))
# searches<-rbind(searches, data.frame(search="o", weight=2, name="LeafCombination(2)-WeakConnection"))
searches<-rbind(searches, data.frame(search="o", weight=3, name="LeafCombination(3)-WeakConnection"))
# searches<-rbind(searches, data.frame(search="o", weight=4, name="LeafCombination(4)-WeakConnection"))
# searches<-rbind(searches, data.frame(search="q", weight=1, name="LeafCombination(1)-WeakConnection-SingleOrigin"))
searches<-rbind(searches, data.frame(search="q", weight=2, name="LeafCombination(2)-WeakConnection-SingleOrigin"))
# searches<-rbind(searches, data.frame(search="q", weight=3, name="LeafCombination(3)-WeakConnection-SingleOrigin"))
# searches<-rbind(searches, data.frame(search="q", weight=4, name="LeafCombination(4)-WeakConnection-SingleOrigin"))
# searches<-rbind(searches, data.frame(search="u", weight=1, name="LeafCombination(1)-WeakConnection-StochasticSearch(max,5)"))
# searches<-rbind(searches, data.frame(search="u", weight=2, name="LeafCombination(2)-WeakConnection-StochasticSearch(max,5)"))
# searches<-rbind(searches, data.frame(search="u", weight=3, name="LeafCombination(3)-WeakConnection-StochasticSearch(max,5)"))
# searches<-rbind(searches, data.frame(search="u", weight=4, name="LeafCombination(4)-WeakConnection-StochasticSearch(max,5)"))



# This decides the disaster length and the type if graph to analyse. If length is -1, it means average all
augmentation_analysis(filename, 'ALL', -1, n_start, n_end, n_interval, searches)
#augmentation_analysis(filename, 'TopBottom', -1, n_start, n_end, n_interval, searches)
#augmentation_analysis(filename, 'MST', -1, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'ALL', 0.001, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'ALL', 5, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'ALL', 100, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'ALL', 300, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'TopBottom', 0.001, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'MST', 0.001, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'TopBottom', 5, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'MST', 5, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'TopBottom', 100, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'MST', 100, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'TopBottom', 300, n_start, n_end, n_interval, searches)
# augmentation_analysis(filename, 'MST', 300, n_start, n_end, n_interval, searches)


# augmentation <- read.table(filename, header = TRUE)
# nrow(augmentation)
# void <- augmentation[is.na(augmentation$cost), ]
# nrow(void)

# Combine data run
#filename <- "data/exec_times/ExperimentalRun1.txt"
#run1 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun2.txt"
#run2 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun3.txt"
#run3 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun4.txt"
#run4 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun5.txt"
#run5 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun6.txt"
#run6 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun7.txt"
#run7 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun8.txt"
#run8 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun9.txt"
#run9 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun10.txt"
#run10 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun11.txt"
#run11 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun12.txt"
#run12 <- read.table(filename, header = TRUE)
#filename <- "data/exec_times/ExperimentalRun13.txt"
#run13 <- read.table(filename, header = TRUE)
#run <- rbind(run1, run2, run3, run4, run5, run6, run7, run8, run9, run10, run11, run12, run13)
#write.table(run, "data/exec_times/ExperimentalRun.txt", col.names = TRUE, row.names = FALSE)
