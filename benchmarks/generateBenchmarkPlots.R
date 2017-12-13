#!/usr/bin/env Rscript

library(ggplot2)
library(Rmisc)

log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}



graph_width = 12
graph_height = 9.5
lab_opt = c("debug", "intel", "release", "pll")
label_tests <- c("PLL", "CLANG -O0", "ICPC -O3", "CLANG -O3")


source_folder <- "."
files <- list.files(source_folder, pattern='.csv')
header = c("nodes", "moves", "time", "metric", "opt", "tree")

data <- data.frame(nodes=as.numeric(), moves=as.numeric(), time=as.numeric(), metric=as.numeric(), opt=as.character(), tree=as.numeric())
for(file in files){
  k <- gsub(".csv","",file)
  x <- read.csv2(paste(source_folder,'/',file,sep=''), sep=',', header=FALSE, stringsAsFactors = FALSE)
  x$V4 <- as.double(x$V4)
  x$opt <- rep(strsplit(k, "_")[[1]][3], times=dim(x)[1])
  x$tree <- rep(strsplit(k, "_")[[1]][2], times=dim(x)[1])
  data <- rbind(data, x)
}

#data$metric <- as.numeric(data$metric);
#data$nodes <- as.numeric(data$nodes);
#data$moves <- as.numeric(data$moves);
#data$time <- as.numeric(data$time);

colnames(data) <- header



data_nopll <- data[which(data$opt != 'pll' & data$tree <= 890),]




tgc <- summarySE(data_nopll, measurevar="time", groupvars=c("moves","opt"))
pd <- position_dodge(0.1) # move them .05 to the left and right



theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(tgc, aes(x=moves, y=time, colour=opt, group = opt)) +
  geom_errorbar(aes(ymin=time-se, ymax=time+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.5, aes(fill=opt), shape=21) + # 21 is filled circle
  #geom_smooth(se=FALSE) +
 # geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE) +
  labs(title="TSHLIB - Benchmark 1: Network perturbation using bidirectional node-edge swapping (d>3, n=10)",
       subtitle="Benchmark performed on 2.9 GHz Intel Core i5 with 8 GB 2133 MHz LPDDR3 using a random generated set of hierarchical binary networks",
       caption="(C) Lorenzo Gatti & Massimo Maiolo 2017",
       x="# node-edge swaps (apply + revert)",
       y="computational time / (microseconds)",
       color=NULL)+
  theme(legend.position = "bottom",
        legend.justification = c("center", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6)) + 
  scale_x_log10(minor_breaks=log10_minor_break()) +annotation_logticks(sides = "b")+
  scale_color_manual(values=c("#1976d2", "#689f38", "#e64a19", "#000000"), 
                     name="Optimisation scheme",
                     breaks=as.factor(lab_opt),
                     labels=lab_opt) + 
  
  scale_fill_manual(values=c("#1976d2", "#689f38", "#e64a19", "#000000"), 
                    name="Optimisation scheme",
                    breaks=as.factor(lab_opt),
                    labels=lab_opt) + guides(fill=FALSE)


ggsave('benchmark1.pdf', g, "pdf", width=graph_width, height=graph_height)

tgc <- summarySE(data_nopll, measurevar="time", groupvars=c("nodes","opt"))
pd <- position_dodge(0.1) # move them .05 to the left and right

theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(tgc, aes(x=nodes, y=time, colour=opt, group = opt)) +
  geom_errorbar(aes(ymin=time-se, ymax=time+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.5, aes(fill=opt), shape=21) + # 21 is filled circle
  geom_smooth(se=FALSE) +
  #geom_smooth(method="lm", aes(color="Exp SubstitutionModel"), formula= (y ~ exp(x)), se=FALSE, linetype = 1) +
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE) +
  labs(title="TSHLIB - Benchmark 2: Full network optimisation (d>3, n=10)",
       subtitle="Benchmark performed on 2.9 GHz Intel Core i5 with 8 GB 2133 MHz LPDDR3 using a random generated set of hierarchical binary networks",
       caption="(C) Lorenzo Gatti & Massimo Maiolo 2017",
       x="network size / (nodes)",
       y="computational time / (microseconds)",
       color=NULL)+
  theme(legend.position = "bottom",
        legend.justification = c("center", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6)) +
  scale_color_manual(values=c("#1976d2", "#689f38", "#e64a19", "#000000"), 
                     name="Optimisation scheme",
                     breaks=as.factor(lab_opt),
                     labels=lab_opt) + 
  
  scale_fill_manual(values=c("#1976d2", "#689f38", "#e64a19", "#000000"), 
                    name="Optimisation scheme",
                    breaks=as.factor(lab_opt),
                    labels=lab_opt) + guides(fill=FALSE)

#
ggsave('benchmark2.pdf', g, "pdf", width=graph_width, height=graph_height)


tgc <- summarySE(data_nopll, measurevar="metric", groupvars=c("nodes","opt"))
pd <- position_dodge(0.1) # move them .05 to the left and right

theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(tgc, aes(x=nodes, y=metric, colour=opt, group = opt)) +
  geom_errorbar(aes(ymin=metric-se, ymax=metric+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.5, aes(fill=opt), shape=21) + # 21 is filled circle
  geom_smooth(se=FALSE) +
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE) +
  labs(title="TSHLIB - Benchmark 3: Scaling performances with ingreasing graph size (d>3, n=10)",
       subtitle="Benchmark performed on 2.9 GHz Intel Core i5 with 8 GB 2133 MHz LPDDR3 using a random generated set of hierarchical binary networks",
       caption="(C) Lorenzo Gatti & Massimo Maiolo 2017",
       x="network size / (nodes)",
       y="swap move/computational time (microseconds)",
       color=NULL)+
  theme(legend.position = "bottom",
        legend.justification = c("center", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6)) +
  scale_color_manual(values=c("#1976d2", "#689f38", "#e64a19", "#000000"), 
                     name="Optimisation scheme",
                     breaks=as.factor(lab_opt),
                     labels=lab_opt) + 
  
  scale_fill_manual(values=c("#1976d2", "#689f38", "#e64a19", "#000000"), 
                    name="Optimisation scheme",
                    breaks=as.factor(lab_opt),
                    labels=lab_opt) + guides(fill=FALSE)
 # scale_colour_hue(name="Optimisation scheme",    # Legend label, use darker colors
 #                  breaks=unique(tgc$opt),
 #                  labels=unique(tgc$opt),
 #                  l=40)  +
            # Use darker colors, lightness=40


ggsave('benchmark3.pdf', g, "pdf", width=graph_width, height=graph_height)



#=======================================
data$tree <- as.numeric(data$tree)
data$metric <- as.numeric(data$metric)
data$moves <- as.numeric(data$moves)
data$time <- as.numeric(data$time)
data_treemax75 <- data[which(data$tree <= 75), ]

tgc <- summarySE(data_treemax75, measurevar="metric", groupvars=c("nodes","opt"))
pd <- position_dodge(0.1) # move them .05 to the left and right

theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(tgc, aes(x=nodes, y=metric, colour=opt, group = opt)) +
  geom_errorbar(aes(ymin=metric-se, ymax=metric+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.5, aes(fill=opt), shape=21) + # 21 is filled circle
  #geom_smooth(se=FALSE) +
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE) +
  labs(title="TSHLIB - Benchmark 4: Rearrangement operation performances w.r.t PLL (d>3, n=10)",
       subtitle="Benchmark performed on 2.9 GHz Intel Core i5 with 8 GB 2133 MHz LPDDR3 using a random generated set of hierarchical binary networks",
       caption="(C) Lorenzo Gatti & Massimo Maiolo 2017",
       x="network size / (nodes)",
       y="swap move/computational time (microseconds)",
       color=NULL)+
  scale_x_log10(minor_breaks=log10_minor_break()) +annotation_logticks(sides = "b") +
  theme(legend.position = "bottom",
        legend.justification = c("center", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6)) +
  #scale_colour_brewer(palette = "Set1", 
  #                    name="Optimisation scheme",    # Legend label, use darker colors
  #                    breaks=lab_opt,
  #                    labels=lab_opt) +  guides(fill=FALSE)

  scale_color_manual(values=c("#1976d2", "#689f38", "#e64a19", "#000000"), 
                  name="Optimisation scheme",
                  breaks=as.factor(lab_opt),
                  labels=lab_opt) + 
  
  scale_fill_manual(values=c("#1976d2", "#689f38", "#e64a19", "#000000"), 
                     name="Optimisation scheme",
                     breaks=as.factor(lab_opt),
                     labels=lab_opt) + guides(fill=FALSE)

ggsave('benchmark4.pdf', g, "pdf", width=graph_width, height=graph_height)