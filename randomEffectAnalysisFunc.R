## random effect analysis function for clusterPower with appropriate t degrees of freedom
## Nicholas Reich
## October 2013

random.effect.DF <- function(dat, incl.period.effect, outcome.type, alpha)  {
        if (outcome.type == "poisson") {
                offsets <- log(dat[, "at.risk.time"])
        }
        else {
                offsets <- rep(0, nrow(dat))
        }
        if (incl.period.effect == 0) {
                if(outcome.type=="gaussian") {
                        fit <- lmer(y ~ trt + (1 | clust), data = dat, offset = offsets)
                } else {
                        fit <- glmer(y ~ trt + (1 | clust), data = dat, family = outcome.type, 
                                     offset = offsets)
                }
        }
        else {
                if(outcome.type=="gaussian") {
                        fit <- lmer(y ~ trt + per + (1 | clust) - 1, data = dat, 
                                    offset = offsets)
                } else {
                        fit <- glmer(y ~ trt + per + (1 | clust) - 1, data = dat, 
                                     family = outcome.type, offset = offsets)
                }
        }
        n.clust <- length(unique(dat$clust))
        df <- n.clust - 2 ## based on k-2 in Donner & Klar p.118
        t <- qt(1 - alpha/2, df=df) * c(-1, 1)
        est <- coef(summary(fit))["trt", ]
        ci <- est["Estimate"] + t * est["Std. Error"]
        return(c(est["Estimate"], ci))
}