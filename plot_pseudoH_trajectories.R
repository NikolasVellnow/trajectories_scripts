#################### positive Selection ####################################

rm(list=ls())
library(lattice)
library(ggplot2)
library("ggpubr")

dat_0.0001 <- read.delim("~/Desktop/Threshold_project/20221008_174406_entropy_delta_entropy_s=0.0001.txt", header=T)
dat_0.0001$selection_coeff <- 0.0001

dat_0.001 <- read.delim("~/Desktop/Threshold_project/20221008_165812_entropy_delta_entropy_s=0.001.txt", header=T)
dat_0.001$selection_coeff <- 0.001

dat_0.01 <- read.delim("~/Desktop/Threshold_project/20221009_101505_entropy_delta_entropy._s=0.01.txt", header=T)
dat_0.01$selection_coeff <- 0.01

dat_0.1 <- read.delim("~/Desktop/Threshold_project/20221009_122739_entropy_delta_entropy_s=0.1.txt", header=T)
dat_0.1$selection_coeff <- 0.1

dat_all <- rbind(dat_0.0001, dat_0.001, dat_0.01, dat_0.1)
dat_all$run <- as.factor(dat_all$run)
dat_all$dominance_coeff <- 0.5

#write.table(dat_all, "~/Desktop/Threshold_project/posSelModel_all_s_entropy_delta_entropy.txt",
#            sep='\t', row.names = F)

mean(dat_all$num_trajectories)

theme_set(
    theme_bw() +
        theme(legend.position = "top")
)

p <- ggplot(dat_all,(aes(x=generation, y=delta_pseudo_entropy, color=run)))+
    geom_line()+
    theme_linedraw()+
    theme(legend.position = "none",
          plot.title=element_text(size=12), strip.text.x = element_text(10))+
    stat_summary(fun=mean, colour="black", geom="line", aes(group = 1), size=1.1)+
    ggtitle("Delta-Entropy (mean # loci: 116) for selection coefficients (columns) & recombination rates (rows)") +
    facet_grid(recombination_rate ~ selection_coeff) +
    theme(strip.text.x = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")),
          strip.text.y = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")))
p

###################### negative selection #########

rm(list=ls())
library(lattice)
library(ggplot2)
library("ggpubr")

dat <- read.delim("~/Desktop/Threshold_project/20221010_115240_entropy_delta_entropy.txt", header=T)

dat$run <- as.factor(dat$run)
dat$dominance_coeff <- 0.5

mean(dat$num_trajectories)

theme_set(
    theme_bw() +
        theme(legend.position = "top")
)

p <- ggplot(dat,(aes(x=generation, y=delta_pseudo_entropy, color=run)))+
    geom_line()+
    theme_linedraw()+
    theme(legend.position = "none",
          plot.title=element_text(size=12), strip.text.x = element_text(10))+
    stat_summary(fun=mean, colour="black", geom="line", aes(group = 1), size=1.1)+
    ggtitle("Delta-Entropy (mean # loci: 86) for selection coefficients (columns) & recombination rates (rows)") +
    facet_grid(recombination_rate ~ selection_coeff) +
    theme(strip.text.x = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")),
          strip.text.y = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")))
p


###################### overdominance selection #########

rm(list=ls())
library(lattice)
library(ggplot2)
library("ggpubr")

dat <- read.delim("~/Desktop/Threshold_project/20221010_172934_entropy_delta_entropy.txt", header=T)

dat$run <- as.factor(dat$run)

mean(dat$num_trajectories)

dat_1 <- dat[dat$dominance_coeff==-0.5,]

theme_set(
    theme_bw() +
        theme(legend.position = "top")
)

p <- ggplot(dat_1,(aes(x=generation, y=delta_pseudo_entropy, color=run)))+
    geom_line()+
    theme_linedraw()+
    theme(legend.position = "none",
          plot.title=element_text(size=12), strip.text.x = element_text(10))+
    stat_summary(fun=mean, colour="black", geom="line", aes(group = 1), size=1.1)+
    ggtitle("Delta-Entropy (mean # loci: 460) for a dominance coefficient of -0.5 & selection coefficients (columns) & recombination rates (rows)") +
    facet_grid(recombination_rate ~ selection_coeff) +
    theme(strip.text.x = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")),
          strip.text.y = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")))
p

###################### negative selection #########

rm(list=ls())
library(lattice)
library(ggplot2)
library("ggpubr")

dat <- read.delim("~/Desktop/Threshold_project/20221010_115240_entropy_delta_entropy.txt", header=T)

dat$run <- as.factor(dat$run)
dat$dominance_coeff <- 0.5

mean(dat$num_trajectories)

theme_set(
    theme_bw() +
        theme(legend.position = "top")
)

p <- ggplot(dat,(aes(x=generation, y=delta_pseudo_entropy, color=run)))+
    geom_line()+
    theme_linedraw()+
    theme(legend.position = "none",
          plot.title=element_text(size=12), strip.text.x = element_text(10))+
    stat_summary(fun=mean, colour="black", geom="line", aes(group = 1), size=1.1)+
    ggtitle("Delta-Entropy (mean # loci: 86) for selection coefficients (columns) & recombination rates (rows)") +
    facet_grid(recombination_rate ~ selection_coeff) +
    theme(strip.text.x = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")),
          strip.text.y = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")))
p


###################### overdominance selection #########

rm(list=ls())
library(lattice)
library(ggplot2)
library("ggpubr")

dat <- read.delim("~/Desktop/Threshold_project/20221010_172934_entropy_delta_entropy.txt", header=T)

dat$run <- as.factor(dat$run)

mean(dat$num_trajectories)

dat_1 <- dat[dat$dominance_coeff==-0.5,]

theme_set(
    theme_bw() +
        theme(legend.position = "top")
)

p <- ggplot(dat_1,(aes(x=generation, y=delta_pseudo_entropy, color=run)))+
    geom_line()+
    theme_linedraw()+
    theme(legend.position = "none",
          plot.title=element_text(size=12), strip.text.x = element_text(10))+
    stat_summary(fun=mean, colour="black", geom="line", aes(group = 1), size=1.1)+
    ggtitle("Delta-Entropy (mean # loci: 460) for a dominance coefficient of -0.5 & selection coefficients (columns) & recombination rates (rows)") +
    facet_grid(recombination_rate ~ selection_coeff) +
    theme(strip.text.x = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")),
          strip.text.y = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")))
p


###################### soft sweep #########

rm(list=ls())
library(lattice)
library(ggplot2)
library("ggpubr")

dat <- read.delim("~/Desktop/Trajectories_project/pseudoEntropySimsSoftSweep_neutral_sel_n400__entropy_delta_entropy.txt", header=T)

dat$run <- as.factor(dat$run)

mean(dat$num_trajectories)

dat_1 <- dat[dat$dominance_coeff==-0.5,]

theme_set(
    theme_bw() +
        theme(legend.position = "top")
)

p <- ggplot(dat_1,(aes(x=generation, y=delta_pseudo_entropy, color=run)))+
    geom_line()+
    theme_linedraw()+
    theme(legend.position = "none",
          plot.title=element_text(size=12), strip.text.x = element_text(10))+
    stat_summary(fun=mean, colour="black", geom="line", aes(group = 1), size=1.1)+
    ggtitle("Delta-Entropy (mean # loci: 460) for a dominance coefficient of -0.5 & selection coefficients (columns) & recombination rates (rows)") +
    facet_grid(recombination_rate ~ selection_coeff) +
    theme(strip.text.x = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")),
          strip.text.y = element_text(size=8, margin=margin(.5,.5,.5,.5, "pt")))
p


