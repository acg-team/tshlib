#!/usr/bin/env Rscript

library(ggplot2)
library(Rmisc)

source_folder <- "."
files <- list.files(source_folder, pattern='.csv')
header = c("nodes", "moves", "time", "measure", "opt", "tree")

data <- data.frame(nodes=as.numeric(), moves=as.numeric(), time=as.numeric(), measure=as.numeric(), opt=as.character(), tree=as.numeric())
for(file in files){
  k <- gsub(".csv","",file)
  x <- read.csv2(paste(source_folder,'/',file,sep=''), sep=',', header=FALSE, stringsAsFactors = FALSE)
  x$opt <- rep(strsplit(k, "_")[[1]][3], times=dim(x)[1])
  x$tree <- rep(strsplit(k, "_")[[1]][2], times=dim(x)[1])
  data <- rbind(data, x)
}
colnames(data) <- header


tgc <- summarySE(data, measurevar="time", groupvars=c("moves","opt"))
pd <- position_dodge(0.1) # move them .05 to the left and right

theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(tgc, aes(x=moves, y=time, colour=opt, group = opt)) +
  geom_errorbar(aes(ymin=time-se, ymax=time+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.5, fill="white", shape=21) + # 21 is filled circle
  geom_smooth(span = 0.3) +
  labs(title="TSHLIB - Benchmark 1: Graph perturbation using bidirectional node swapping (d>3, n=10)",
       subtitle="Benchmark performed on 2.9 GHz Intel Core i5 with 8 GB 2133 MHz LPDDR3 using a random generated set of binary graphs",
       caption="(C) Lorenzo Gatti & Massimo Maiolo 2017",
       x="Swaps (node swaps + revert swaps)",
       y="Time (microseconds)",
       color=NULL)+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0))+
  scale_colour_hue(name="Optimisation scheme",      # Legend label, use darker colors
                   breaks=unique(tgc$opt),
                   labels=c("CLANG -O0", "ICPC -O3", "CLANG -O3"),
                   l=40)                            # Use darker colors, lightness=40


ggsave('benchmark1.pdf', g, "pdf", width=10, height=10)

tgc <- summarySE(data, measurevar="time", groupvars=c("nodes","opt"))
pd <- position_dodge(0.1) # move them .05 to the left and right

theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(tgc, aes(x=nodes, y=time, colour=opt, group = opt)) +
  geom_errorbar(aes(ymin=time-se, ymax=time+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.5, fill="white", shape=21) + # 21 is filled circle
  geom_smooth(span = 0.3) +
  labs(title="TSHLIB - Benchmark 2: Full tree-space search (d>3, n=10)",
       subtitle="Benchmark performed on 2.9 GHz Intel Core i5 with 8 GB 2133 MHz LPDDR3 using a random generated set of binary graphs",
       caption="(C) Lorenzo Gatti & Massimo Maiolo 2017",
       x="Nodes (internal + terminal)",
       y="Time (microseconds)",
       color=NULL)+
  theme(legend.justification=c(1,0),
      legend.position=c(1,0))+
  scale_colour_hue(name="Optimisation scheme",    # Legend label, use darker colors
                   breaks=unique(tgc$opt),
                   labels=c("CLANG -O0", "ICPC -O3", "CLANG -O3"),
                   l=40)                    # Use darker colors, lightness=40


ggsave('benchmark2.pdf', g, "pdf", width=10, height=10)
