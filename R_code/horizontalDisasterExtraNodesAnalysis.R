inputs <- c("MST", "TopBottom")
ns <- seq(10,100,10)
ss <- 1:20
searches <- c("o", "q", "u")
lengths <- c("0.001", "5", "100", "300", "0.001000", "5.000000", "100.000000", "300.000000")

data <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("size", "seed", "input_type", "disaster", "search", "vertices"))

for(input in inputs){
  for(n in ns){
    for(s in ss){
      for(search in searches){
        for(length in lengths){
          filename <- paste0("data/plots/CompletedGraphs/n", n, "_s", s, "_l", length, "_", input, "_search_", search, ".svg")
          if(!file.exists(filename)){next;}
          reader <- file(filename, open = "r")
          readLines(reader, n = 3, warn = FALSE)
          vertices <- 0
          while(substr(readLines(reader, n = 1, warn = FALSE), 3, 8) == "circle" && 
                substr(readLines(reader, n = 1, warn = FALSE), 3, 6) == "text"){
            vertices <- vertices + 1
          }
          data <- rbind(data, data.frame(size = c(n), seed = c(s), input_type = c(input), disaster = c(length), search = c(search), vertices = c(vertices)))
          close(reader)
        }
      }
    }
  }
}

means <- c()
maxs <- c()
mins <- c()

for(n in ns){
  
  means <- c(means, mean(data[data$size == n, "vertices"])/n)
  maxs <- c(maxs, max(data[data$size == n, "vertices"])/n)
  mins <- c(mins, min(data[data$size == n, "vertices"])/n)
  
}


range <- c(0.9 * min(mins), 1.1 *  max(maxs))

png("data/analysis_plots/average_vertices_increase.png", width = 682, height = 512)
plot(ns, maxs, 
     col = "black", 
     type = "b",
     ylab = "Final graph size ratio",
     xlab = "Initial graph size",
     ylim = range,
     xaxt = "n",
     pch = 1)

points(ns, means, col = "black", type = "b", pch = 2)
points(ns, mins, col = "black", type = "b", pch = 3)

legend("topleft", legend = c("Maximum", "Average", "Minimum"), pch = c(1, 2, 3))
axis(1, at = ns)
out <- dev.off()




