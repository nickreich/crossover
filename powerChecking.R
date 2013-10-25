## Climo et al. 2013 power example

library(clusterPower)
library(lme4)
setwd("~/Documents/work/research/CRXO/manuscripts/LetterToEditor")
source("randomEffectAnalysisFunc.R")
NSIM <- 500

## numbers based on Climo protocol downloaded from NEJM website and paper

monthly.admissions <- c(123.8, 46.3, 51.6, 85.3, 41.8, 111.6, 55.8, 62.3, 72.7)
period.admissions <- round(6*monthly.admissions)
avg.clust.size.per.period <- round(mean(period.admissions))

#monthly.patient.days <- c(692.3, 285.7, 285.7, 425.9, 786.3, 598.8, 299.1, 316.3, 467.1)

mean.length.of.stays <- c(5.6, 6.2, 5.5, 5.0, 18.8, 5.4, 5.4, 5.1, 6.4)
mean.los.overall <- round(mean(mean.length.of.stays))

baselines <- c(8.1, 9.6, 0, 0.4, 5.5, 3.1, 8.5, 2.2, 8.7)/1000
baseline <- mean(baselines)
## from p. 18 section 6.10 of protocol of Climo, they assumed baseline of 0.008
        
delta <- .75 ## based on p.18, section 6.10 of protocol

## qualitatively check bcv against actual baseline numbers
## y ~ pois(\lambda)
## log \lambda = baseline + alpha_c
## alpha_c ~ norm(0, bcv)
clust.means.dist <- exp(rnorm(10000, mean=log(mean(baselines)), sd=bcv^2))
dens <- density(clust.means.dist)
plot(dens, xlim=range(c(dens$x, baselines)))
rug(baselines)

bcv <- .5 ## not given in protocol of Climo et al. 
## BCV chosen to match qualitatively assumed distribution to the baseline numbers. 
## Note: As bcv gets small, crxo and cr results converge.


################################
## single simulation for Climo study itself

crxo.power1 <- power.sim.poisson(n.sim=NSIM, effect.size=log(delta), 
                                n.clusters=12, n.periods=2,
                                cluster.size = avg.clust.size.per.period,
                                at.risk.params = mean.los.overall,
                                period.effect = log(baseline), period.var=0,
                                btw.clust.var=bcv,
                                estimation.function=random.effect.DF,
                                verbose=TRUE)

## n.clusters stays the same, cluster size changes
cr.power1 <- power.sim.poisson(n.sim=NSIM, effect.size=log(delta), 
                              n.clusters=12, n.periods=1, 
                              cluster.size = 5*avg.clust.size.per.period,
                              at.risk.params = mean.los.overall,
                              period.effect = log(baseline), period.var=0,
                              btw.clust.var=bcv,
                              estimation.function=random.effect.DF,
                              verbose=TRUE)

## check profiling
#p <- profr(crxo.power1 <- power.sim.poisson(n.sim=1, effect.size=log(delta), 
#                                 n.clusters=12, n.periods=2,
#                                 cluster.size = avg.clust.size.per.period,
#                                 at.risk.params = mean.los.overall,
#                                 period.effect = log(baseline), period.var=0,
#                                 btw.clust.var=bcv,
#                                 estimation.function=random.effect.DF,
#                                 verbose=TRUE)
#)
#plot(p)

##############################
## simulate power for Climo paper across a range of number of clusters
## NOTE: cluster size doubles for CRTs to have the same total person time as CRXOs

library(doMC)
n.clusters <- seq(10, 100, by=10)
registerDoMC(10)

NSIM=500

## get power for cluster-randomized designs
clust.power <- foreach(n.clust=n.clusters,.combine='c') %dopar% {
        cr.power2 <- power.sim.poisson(n.sim=NSIM, effect.size=log(delta), 
                                       n.clusters=n.clust, n.periods=1, 
                                       cluster.size = 2*avg.clust.size.per.period,
                                       at.risk.params = mean.los.overall,
                                       period.effect = log(baseline), period.var=0,
                                       btw.clust.var=bcv,
                                       estimation.function=random.effect.DF,
                                       verbose=TRUE)
        cr.power2$power
}

## get power for cluster-randomized crossover designs
crxo.power <- foreach(n.clust=n.clusters,.combine='c') %dopar% {
        cr.power2 <- power.sim.poisson(n.sim=NSIM, effect.size=log(delta), 
                                       n.clusters=n.clust, n.periods=2, 
                                       cluster.size = avg.clust.size.per.period,
                                       at.risk.params = mean.los.overall,
                                       period.effect = log(baseline), period.var=0,
                                       btw.clust.var=bcv,
                                       estimation.function=random.effect.DF,
                                       verbose=TRUE)
        cr.power2$power
}

save.image(file="run_20131024_bcv5.rda")

## plot single
qplot(n.clusters, clust.power*100, geom=c("point", "line"), se=FALSE) + 
        ylim(0,100) + ylab("Power (in %)") +
        xlim(0, 100) +
        xlab("Number of clusters") 

single.data <- data.frame(n.clusters=12, power=c(cr.power1$power, crxo.power1$power)*100, 
                          name=c("cluster-randomized", "cluster-randomized crossover"))

## plot both
crxo.dat <- cbind(n.clusters, power=crxo.power*100, name="cluster-randomized crossover")
crxo.dat <- rbind(crxo.dat, c(12, crxo.power1$power*100, "cluster-randomized crossover"))

cr.dat <- cbind(n.clusters, power=clust.power*100, name="cluster-randomized")
cr.dat <- rbind(cr.dat, c(12, cr.power1$power*100, "cluster-randomized"))

dat <- as.data.frame(rbind(crxo.dat, cr.dat))
dat$n.clusters <- as.numeric(as.character(dat$n.clusters))
dat$power <- as.numeric(as.character(dat$power))
dat$name <- relevel(dat$name, ref="cluster-randomized crossover")

qplot(n.clusters, power, geom=c("point", "line"), color=name, data=dat, se=FALSE) + 
        scale_y_continuous(limits = c(0, 101)) + ylab("Power (in %)") +
        scale_x_continuous(limits = c(0, 100), breaks=seq(0, 100, 20), expand=c(.02,.02)) + xlab("Number of clusters") +
        theme(legend.justification=c(1,0), legend.position=c(.9,.4), legend.title=element_blank()) +
        scale_color_brewer(palette="Set1") + 
        #stat_smooth(method="lm", formula = y ~ ns(x, 7), se=FALSE) +
        geom_point(aes(n.clusters, power), shape=1, size=4, data=single.data)
