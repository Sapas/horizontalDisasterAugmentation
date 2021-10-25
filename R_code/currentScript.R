source("performanceAnalysis.R")


# First call average performance functions. Will not be used but tells you what was not solved
averagePre <- getAveragePerformance("preliminaryRun", 1, 1040)
averageMain <- getAveragePerformance("mainRun", 1, 1600)
# averagePoly <- getAveragePerformance("polynomialRun", 1, 640)
# averageIntermediate <- getAveragePerformance("intermediateWeightsRun", 1, 1680)

# Now grab all of the data
combinedPre <- combineArrayResults("preliminaryRun", 1, 1040)
combinedMain <- combineArrayResults("mainRun", 1, 1600)
combinedPoly <- combineArrayResults("polynomialRun", 1, 640)
combinedIntermediate <- combineArrayResults("intermediateWeightsRun", 1, 1680)

# Prelim analysis
full_augmentation_analysis("preRun", combinedPre, 10, 100, 10, c("a(0)", "b(0)", "o(2)", "p(0)", "s(0)"), 1.0)
#interest <- combinedPre[combinedPre$size == 100 & combinedPre$disaster == 300 & combinedPre$input_type == "TopBottom" & combinedPre$search == "s(0)", ]


# Main analysis
full_augmentation_analysis("mainRunInitial", combinedMain, 10, 100, 10, c("o(2)", "q(2)", "u(2)"))
# Weight analysis
combinedAll <- rbind(combinedMain, combinedPoly)
combinedAll <- rbind(combinedAll, combinedIntermediate)
#full_augmentation_analysis("mainRunWeightWeakConnection", combinedAll, 10, 100, 10, c("o(1)", "o(1.5)", "o(2)", "o(2.5)", "o(3)", "o(3.5)", "o(4)"))
full_augmentation_analysis("mainRunWeightWeakConnection", combinedAll, 10, 100, 10, c("o(1)", "o(2)", "o(3)", "o(4)"))
#full_augmentation_analysis("mainRunWeightWeakConnectionSingleOrigin", combinedAll, 10, 100, 10, c("q(1)", "q(1.5)", "q(2)", "q(2.5)", "q(3)", "q(3.5)", "q(4)"))
full_augmentation_analysis("mainRunWeightWeakConnectionSingleOrigin", combinedAll, 10, 100, 10, c("q(1)", "q(2)", "q(3)", "q(4)"))
#full_augmentation_analysis("mainRunWeightWeakConnectionStochastic", combinedAll, 10, 100, 10, c("u(1)", "u(1.5)", "u(2)", "u(2.5)", "u(3)", "u(3.5)", "u(4)"))
full_augmentation_analysis("mainRunWeightWeakConnectionStochastic", combinedAll, 10, 100, 10, c("u(1)", "u(2)", "u(3)", "u(4)"))
 
# Poly time analysis, compare performance
full_augmentation_analysis("mainRunPoly", combinedAll, 10, 100, 10, c("o(3)", "q(2)", "v(3)", "w(2)"))
# Based on results, might want to show TopBottom 300 in one graph, and the rest in another
allMinusLarge <- combinedAll[combinedAll$disaster != 300 | combinedAll$input_type != "TopBottom", ]
augmentation_analysis("mainRunPolyMinus", allMinusLarge, "ALL", "ALL", 10, 100, 10, c("o(3)", "q(2)", "v(3)", "w(2)"))

# Finally plot of the two best performing algorithms
full_augmentation_analysis("finalGraphs", combinedAll, 10, 100, 10, c("v(3)", "w(2)"))
full_augmentation_analysis("finalGraphs", combinedAll, 10, 100, 10, c("o(3)", "q(2)"))

combinedAll <- combinedAll[!is.na(combinedAll$cost), ]
combinedAll <- combinedAll[!is.infinite(combinedAll$cost), ]

resV3cost <- get_data_row(combinedAll, "v(3)", 10, 100, 10, "cost", "mean")[[1]]
resW2cost <- get_data_row(combinedAll, "w(2)", 10, 100, 10, "cost", "mean")[[1]]
resV3time <- get_data_row(combinedAll, "v(3)", 10, 100, 10, "augment_time", "mean")[[1]]
resW2time <- get_data_row(combinedAll, "w(2)", 10, 100, 10, "augment_time", "mean")[[1]]

(resW2cost - resV3cost) / resW2cost * 100
(resV3time - resW2time) / resV3time * 100

# Working out statistics to report on
combinedPre <- combinedPre[!is.na(combinedPre$cost), ]
combinedPre <- combinedPre[!is.infinite(combinedPre$cost), ]
combinedPre <- combinedPre[combinedPre$disaster == 0.001, ]

resO2 <- get_data_row(combinedPre, "o(2)", 10, 100, 10, "cost", "mean")[[1]]
resP2 <- get_data_row(combinedPre, "p(0)", 10, 100, 10, "cost", "mean")[[1]]
resA <- get_data_row(combinedPre, "a(0)", 10, 100, 10, "cost", "mean")[[1]]
# 
(resP2 - resO2) / resP2 * 100
(resA - resO2) / resA * 100
# 
resA <- get_data_row(combinedPre, "a(0)", 10, 100, 10, "augment_time", "mean")[[1]]
resO2 <- get_data_row(combinedPre, "o(2)", 10, 100, 10, "augment_time", "mean")[[1]]
resB <- get_data_row(combinedPre, "b(0)", 10, 100, 10, "augment_time", "mean")[[1]]
resS <- get_data_row(combinedPre, "s(0)", 10, 100, 10, "augment_time", "mean")[[1]]
resP2 <- get_data_row(combinedPre, "p(0)", 10, 100, 10, "augment_time", "mean")[[1]]
# 
resO2 / resB
resP2 / resO2

# 
resO2
resB
resP2
# 
(resA - resO2) / resA * 100
 
 