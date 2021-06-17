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
  else if(name == "o(2)"){return("LeafCombination(2)-WeakConnection")}
  else if(name == "o(3)"){return("LeafCombination(3)-WeakConnection")}
  else if(name == "o(4)"){return("LeafCombination(4)-WeakConnection")}
  else if(name == "p(0)"){return("LeafCombination(2)-StrongConnection")}
  else if(name == "q(1)"){return("LeafCombination(1)-WeakConnection-SingleOrigin")}
  else if(name == "q(2)"){return("LeafCombination(2)-WeakConnection-SingleOrigin")}
  else if(name == "q(3)"){return("LeafCombination(3)-WeakConnection-SingleOrigin")}
  else if(name == "q(4)"){return("LeafCombination(4)-WeakConnection-SingleOrigin")}
  else if(name == "r(0)"){return("Unkown-r(0)")}
  else if(name == "s(0)"){return("LeafCombination(2)-StrongConnection-CutDistance")}
  else if(name == "t(0)"){return("Unkown-t(0)")}
  else if(name == "u(1)"){return("LeafCombination(1)-WeakConnection-StochasticSearch(max,5)")}
  else if(name == "u(2)"){return("LeafCombination(2)-WeakConnection-StochasticSearch(max,5)")}
  else if(name == "u(3)"){return("LeafCombination(3)-WeakConnection-StochasticSearch(max,5)")}
  else if(name == "u(4)"){return("LeafCombination(4)-WeakConnection-StochasticSearch(max,5)")}
  else if(name == "v(1)"){return("FullyPolynomial-LeafCombination(1)-WeakConnection")}
  else if(name == "v(2)"){return("FullyPolynomial-LeafCombination(2)-WeakConnection")}
  else if(name == "v(3)"){return("FullyPolynomial-LeafCombination(3)-WeakConnection")}
  else if(name == "v(4)"){return("FullyPolynomial-LeafCombination(4)-WeakConnection")}
  else if(name == "w(1)"){return("FullyPolynomial-LeafCombination(1)-WeakConnection-SingleOrigin")}
  else if(name == "w(2)"){return("FullyPolynomial-LeafCombination(2)-WeakConnection-SingleOrigin")}
  else if(name == "w(3)"){return("FullyPolynomial-LeafCombination(3)-WeakConnection-SingleOrigin")}
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