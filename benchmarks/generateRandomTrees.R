#!/usr/bin/env Rscript

library(stringi)
library(ape)
for(i in seq(from = 5, to = 1005, by=10)){

  tips = stri_rand_strings(i, 3)

  a = rtree(i, rooted = TRUE, tip.label=tips, br=runif,min = 1,max = 10)
  write.tree(a,file = paste("data/tree_",i,".nwk",sep=''))
}
