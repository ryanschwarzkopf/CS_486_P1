# randomSequence generates a string consisting of ACGT in an even multinomial distribution
# PRE: length of desired string
# POST: returns random string of length with even distribution
randomSequence <- function(length) {
    p <- c("A"=0.25, "C"=0.25, "G"=0.25, "T"=0.25) # probability of A,C,G,T
    x <- rmultinom(n=length, size=1, prob=p) # create the multinomial array
    w <- apply(x, 2, function(y) {names(p)[which(y==1)]}) # name every element of x
    s <- paste(w, collapse='') # concat w into one string
    #cat(s, "\n")
    return(s)
} # end randomSequence

# patternCount finds the number of times pattern occurs in a string.
# PRE: string of text, pattern to be searched
# POST: pattern count
patternCount <- function(pattern, text) {                                                
    count <- 0
    k <- nchar(pattern)
    n <- nchar(text)
    for(i in 1:n-k+1) { # check all substrings of the same length of pattern and compare 
        if(substr(text, i, i+k-1) == pattern) count <- count+1
    }
    return(count)    
} # end patternCount                        

# frequenWords finds the most frequent k-mer(s) in a string
# PRE: string of text, and length of k-mers to be searched
# POST: most frequent k-mer(s)
frequentWords <- function(text, k) {
    freqPatterns <- c()
    n <- nchar(text)
    count <- array(0, c(1, n-k+1))
    len = length(count)
    for(i in 1:len) { # patternCount every kmer in the string, hold in array
        pattern <- substr(text, i, i+k-1)
        count[i] <- patternCount(pattern, text)
    }                
    for(i in 1:len) { # find the max(es) in count and add the kmer values to frequentPatterns
        if(count[i] == max(count)) {
            pattern <- substr(text, i, i+k-1)           
            if(is.null(freqPatterns) || is.na(match(pattern, freqPatterns))){
                freqPatterns <- c(freqPatterns, pattern)
            }
        }        
    }                
    return (freqPatterns)
} # end frequentWords

k_list <- c(3, 6, 9, 12, 15)
printData("output.txt", k_list)

              
# read in .txt file with one string and finds frequentWords for list of k-mers
# print most frequent k-mers
printData <- function(filename, k_list) {
  SARSCoV2 <- read.delim(k_list, "", header = FALSE)
  len <- length(k_list)
  for(i in len) {
    print(k_list[i] + ": ")
    k_list[i] <- frequentWords(SARSCoV2, k_list[i])
    print(k_list[i])
  }
  return(k_list)
}




install.packages("hash")
require(hash)

betterFrequentWords <- function(text, k) {
    freqPatterns = c()
    n = nchar(text)
    count = as.list.hash(c())
    l = n-k+1
    max_count = 0
    
    for(i in 1:l) {
      
      pattern = substr(text, i, i+k-1)
      
      if(is.null(count) || is.na(match(pattern, names(count)))){
        count[pattern] = patternCount(pattern, text)
        
        if(count[pattern] > max_count){
          max_count = unname(unlist(count[pattern]))
        }
      }
    } 
    
    freqPatterns = names(count[count == max_count])
    
    return (freqPatterns)
}

test = randomSequence(10)
print(frequentWords(test, 3))
print(betterFrequentWords(test,3))

install.packages("rbenchmark")
require(rbenchmark)

test = randomSequence(10)
benchmark(frequentWords(test, 3))$elapsed
benchmark(betterFrequentWords(test, 3))$elapsed

compareFrequentWords <- function(l, k){
    text = randomSequence(l)
    r_hash = c()
    r_count = c()
    
    r_count = c(r_count, benchmark(frequentWords(text, k), replications = 1)$elapsed)
    r_hash = c(r_hash, benchmark(betterFrequentWords(text, k), replications = 1)$elapsed)
    
    return(c(r_hash,r_count))
}

compareFrequentWords(10,3)

runAnalysis <- function(){
  
  runtimes = data.frame('length'=numeric(),'k'=numeric(),'time'=numeric(),'method'=character())
  lengths = c(250,500,750,1000,1500)
  k_vals = c(3,6,9,12,15)
  
  for(i in k_vals){
    for(j in lengths){
      
      result = compareFrequentWords(j,i)
      
      add = nrow(runtimes) + 1
      runtimes[add,] = c(j, i, result[1], 'hash')
      
      add = nrow(runtimes) + 1
      runtimes[add,] = c(j, i, result[2], 'count')
    }
  }
  
  runtimes$length = as.numeric(runtimes$length)
  runtimes$time = as.numeric(runtimes$time)
  runtimes$k = as.numeric(runtimes$k)
  return(runtimes)
}

analysis = runAnalysis()

install.packages('ggplot2')
require(ggplot2)

ggplot(data = analysis, aes(x=length,y=time,col=method, shape=factor(k))) + geom_point() + ggtitle('All Runtimes')
ggplot(data = analysis[analysis$k == 3,], aes(x=length,y=time,col=method)) + geom_point() + geom_smooth() + ggtitle('k=3 runtimes')
ggplot(data = analysis[analysis$k == 6,], aes(x=length,y=time,col=method)) + geom_point() + geom_smooth() + ggtitle('k=6 runtimes')
ggplot(data = analysis[analysis$k == 9,], aes(x=length,y=time,col=method)) + geom_point() + geom_smooth() + ggtitle('k=9 runtimes')
ggplot(data = analysis[analysis$k == 12,], aes(x=length,y=time,col=method)) + geom_point() + geom_smooth() + ggtitle('k=12 runtimes')
ggplot(data = analysis[analysis$k == 15,], aes(x=length,y=time,col=method)) + geom_point() + geom_smooth() + ggtitle('k=15 runtimes')
