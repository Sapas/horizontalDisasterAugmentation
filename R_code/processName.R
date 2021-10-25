processName <- function(name){
  if(name == "a(0)"){return("BlockJoin-StrongConnection")}
  else if(name == "b(0)"){return("BlockJoin-StrongConnection-CutDistance")}
  else if(name == "c(0)"){return("LeafOrigin-StrongConnection")}
  else if(name == "d(0)"){return("LeafOrigin-StrongConnection-CutDistance")}
  else if(name == "e(0)"){return("LeafOrigin-WeakConnection")}
  else if(name == "f(0)"){return("Unkown-f(0)")}
  else if(name == "g(0)"){return("Unkown-g(0)")}
  else if(name == "h(0)"){return("Unkown-h(0)")}
  else if(name == "i(0)"){return("LeafJoin-WeakConnection")}
  else if(name == "j(0)"){return("LeafJoin-StrongConnection")}
  else if(name == "k(0)"){return("LeafJoin-StrongConnection-CutDistance")}
  else if(name == "l(0)"){return("Unkown-l(0)")}
  else if(name == "m(0)"){return("Unkown-m(0)")}
  else if(name == "n(0)"){return("Unkown-n(0)")}
  
  else if(name == "o(1)"){return("LeafCombination(1)-WeakConnection")}
  else if(name == "o(1.5)"){return("LeafCombination(1.5)-WeakConnection")}
  else if(name == "o(2)"){return("LeafCombination(2)-WeakConnection")}
  else if(name == "o(2.5)"){return("LeafCombination(2.5)-WeakConnection")}
  else if(name == "o(3)"){return("LeafCombination(3)-WeakConnection")}
  else if(name == "o(3.5)"){return("LeafCombination(3.5)-WeakConnection")}
  else if(name == "o(4)"){return("LeafCombination(4)-WeakConnection")}
  
  else if(name == "p(0)"){return("LeafCombination(2)-StrongConnection")}
  
  else if(name == "q(1)"){return("LeafCombination(1)-WeakConnection-SingleOrigin")}
  else if(name == "q(1.5)"){return("LeafCombination(1.5)-WeakConnection-SingleOrigin")}
  else if(name == "q(2)"){return("LeafCombination(2)-WeakConnection-SingleOrigin")}
  else if(name == "q(2.5)"){return("LeafCombination(2.5)-WeakConnection-SingleOrigin")}
  else if(name == "q(3)"){return("LeafCombination(3)-WeakConnection-SingleOrigin")}
  else if(name == "q(3.5)"){return("LeafCombination(3.5)-WeakConnection-SingleOrigin")}
  else if(name == "q(4)"){return("LeafCombination(4)-WeakConnection-SingleOrigin")}
  
  else if(name == "r(0)"){return("Unkown-r(0)")}
  else if(name == "s(0)"){return("LeafCombination(2)-StrongConnection-CutDistance")}
  else if(name == "t(0)"){return("Unkown-t(0)")}
  
  else if(name == "u(1)"){return("LeafCombination(1)-WeakConnection-StochasticSearch(max,5)")}
  else if(name == "u(1.5)"){return("LeafCombination(1.5)-WeakConnection-StochasticSearch(max,5)")}
  else if(name == "u(2)"){return("LeafCombination(2)-WeakConnection-StochasticSearch(max,5)")}
  else if(name == "u(2.5)"){return("LeafCombination(2.5)-WeakConnection-StochasticSearch(max,5)")}
  else if(name == "u(3)"){return("LeafCombination(3)-WeakConnection-StochasticSearch(max,5)")}
  else if(name == "u(3.5)"){return("LeafCombination(3.5)-WeakConnection-StochasticSearch(max,5)")}
  else if(name == "u(4)"){return("LeafCombination(4)-WeakConnection-StochasticSearch(max,5)")}
  
  else if(name == "v(1)"){return("FullyPolynomial-LeafCombination(1)-WeakConnection")}
  else if(name == "v(1.5)"){return("FullyPolynomial-LeafCombination(1.5)-WeakConnection")}
  else if(name == "v(2)"){return("FullyPolynomial-LeafCombination(2)-WeakConnection")}
  else if(name == "v(2.5)"){return("FullyPolynomial-LeafCombination(2.5)-WeakConnection")}
  else if(name == "v(3)"){return("FullyPolynomial-LeafCombination(3)-WeakConnection")}
  else if(name == "v(3.5)"){return("FullyPolynomial-LeafCombination(3.5)-WeakConnection")}
  else if(name == "v(4)"){return("FullyPolynomial-LeafCombination(4)-WeakConnection")}
  
  else if(name == "w(1)"){return("FullyPolynomial-LeafCombination(1)-WeakConnection-SingleOrigin")}
  else if(name == "w(1.5)"){return("FullyPolynomial-LeafCombination(1.5)-WeakConnection-SingleOrigin")}
  else if(name == "w(2)"){return("FullyPolynomial-LeafCombination(2)-WeakConnection-SingleOrigin")}
  else if(name == "w(2.5)"){return("FullyPolynomial-LeafCombination(2.5)-WeakConnection-SingleOrigin")}
  else if(name == "w(3)"){return("FullyPolynomial-LeafCombination(3)-WeakConnection-SingleOrigin")}
  else if(name == "w(3.5)"){return("FullyPolynomial-LeafCombination(3.5)-WeakConnection-SingleOrigin")}
  else if(name == "w(4)"){return("FullyPolynomial-LeafCombination(4)-WeakConnection-SingleOrigin")}
  
  else{return(paste0("Unknown", name))}
}

processNames <- function(names){
  processedNames <- c()
  for(name in names){
    processedNames <- c(processedNames, processName(name))
  }
  return(processedNames)
}
  

processPCH <- function(name){
  if(name == "a(0)"){return(1)}
  else if(name == "b(0)"){return(2)}
  else if(name == "c(0)"){return(3)}
  else if(name == "d(0)"){return(4)}
  else if(name == "e(0)"){return(5)}
  else if(name == "f(0)"){return(0)}
  else if(name == "g(0)"){return(0)}
  else if(name == "h(0)"){return(0)}
  else if(name == "i(0)"){return(6)}
  else if(name == "j(0)"){return(7)}
  else if(name == "k(0)"){return(8)}
  else if(name == "l(0)"){return(0)}
  else if(name == "m(0)"){return(0)}
  else if(name == "n(0)"){return(0)}
  
  else if(name == "o(1)"){return(9)}
  else if(name == "o(1.5)"){return(0)}
  else if(name == "o(2)"){return(10)}
  else if(name == "o(2.5)"){return(0)}
  else if(name == "o(3)"){return(11)}
  else if(name == "o(3.5)"){return(0)}
  else if(name == "o(4)"){return(12)}
  
  else if(name == "p(0)"){return(13)}
  
  else if(name == "q(1)"){return(15)}
  else if(name == "q(1.5)"){return(0)}
  else if(name == "q(2)"){return(16)}
  else if(name == "q(2.5)"){return(0)}
  else if(name == "q(3)"){return(17)}
  else if(name == "q(3.5)"){return(0)}
  else if(name == "q(4)"){return(18)}
  
  else if(name == "r(0)"){return(0)}
  else if(name == "s(0)"){return(18)}
  else if(name == "t(0)"){return(0)}
  
  else if(name == "u(1)"){return(22)}
  else if(name == "u(1.5)"){return(0)}
  else if(name == "u(2)"){return(23)}
  else if(name == "u(2.5)"){return(0)}
  else if(name == "u(3)"){return(24)}
  else if(name == "u(3.5)"){return(0)}
  else if(name == "u(4)"){return(25)}
  
  else if(name == "v(1)"){return(0)}
  else if(name == "v(1.5)"){return(0)}
  else if(name == "v(2)"){return(0)}
  else if(name == "v(2.5)"){return(0)}
  else if(name == "v(3)"){return(0)}
  else if(name == "v(3.5)"){return(0)}
  else if(name == "v(4)"){return(0)}
  
  else if(name == "w(1)"){return(0)}
  else if(name == "w(1.5)"){return(0)}
  else if(name == "w(2)"){return(0)}
  else if(name == "w(2.5)"){return(0)}
  else if(name == "w(3)"){return(0)}
  else if(name == "w(3.5)"){return(0)}
  else if(name == "w(4)"){return(0)}
  
  else{return(0)}
}

processPCHs <- function(names){
  processedPCHs <- c()
  for(name in names){
    processedPCHs <- c(processedPCHs, processPCH(name))
  }
  return(processedPCHs)
}

processColour <- function(name){
  if(name == "a(0)"){return("chocolate2")}
  else if(name == "b(0)"){return("chartreuse3")}
  else if(name == "c(0)"){return("cadetblue3")}
  else if(name == "d(0)"){return("burlywood3")}
  else if(name == "e(0)"){return("brown3")}
  else if(name == "f(0)"){return("black")}
  else if(name == "g(0)"){return("black")}
  else if(name == "h(0)"){return("black")}
  else if(name == "i(0)"){return("gold2")}
  else if(name == "j(0)"){return("dodgerblue")}
  else if(name == "k(0)"){return("yellow")}
  else if(name == "l(0)"){return("black")}
  else if(name == "m(0)"){return("black")}
  else if(name == "n(0)"){return("black")}
  
  else if(name == "o(1)"){return("red")}
  else if(name == "o(1.5)"){return("black")}
  else if(name == "o(2)"){return("red2")}
  else if(name == "o(2.5)"){return("black")}
  else if(name == "o(3)"){return("red3")}
  else if(name == "o(3.5)"){return("black")}
  else if(name == "o(4)"){return("red4")}
  
  else if(name == "p(0)"){return("cornflowerblue")}
  
  else if(name == "q(1)"){return("steelblue1")}
  else if(name == "q(1.5)"){return("black")}
  else if(name == "q(2)"){return("steelblue3")}
  else if(name == "q(2.5)"){return("black")}
  else if(name == "q(3)"){return("steelblue4")}
  else if(name == "q(3.5)"){return("black")}
  else if(name == "q(4)"){return("steelblue")}
  
  else if(name == "r(0)"){return("black")}
  else if(name == "s(0)"){return("slategray2")}
  else if(name == "t(0)"){return("black")}
  
  else if(name == "u(1)"){return("sienna1")}
  else if(name == "u(1.5)"){return("black")}
  else if(name == "u(2)"){return("sienna3")}
  else if(name == "u(2.5)"){return("black")}
  else if(name == "u(3)"){return("sienna4")}
  else if(name == "u(3.5)"){return("black")}
  else if(name == "u(4)"){return("sienna")}
  
  else if(name == "v(1)"){return("black")}
  else if(name == "v(1.5)"){return("black")}
  else if(name == "v(2)"){return("black")}
  else if(name == "v(2.5)"){return("black")}
  else if(name == "v(3)"){return("black")}
  else if(name == "v(3.5)"){return("black")}
  else if(name == "v(4)"){return("black")}
  
  else if(name == "w(1)"){return("black")}
  else if(name == "w(1.5)"){return("black")}
  else if(name == "w(2)"){return("black")}
  else if(name == "w(2.5)"){return("black")}
  else if(name == "w(3)"){return("black")}
  else if(name == "w(3.5)"){return("black")}
  else if(name == "w(4)"){return("black")}
  
  else{return("black")}
}

processColours <- function(names){
  processedColours <- c()
  for(name in names){
    processedColours <- c(processedColours, processColour(name))
  }
  return(processedColours)
}


