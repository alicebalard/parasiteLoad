devtools::install_github("alicebalard/Parasite_Load")

simpara <- c(k = 2,
             "male:old.alpha" = 1.4,
             "male:young.alpha" = 1.2,
             "male:baby.alpha" = 1.0,
             "female:old.alpha" = 2.0,
             "female:young.alpha" = -1.8,
             "female:baby.alpha" = 1.1,
             "male:old.inter" = 14,
             "male:young.inter" = 12,
             "male:baby.inter" = 10,
             "female:old.inter" = 20,
             "female:young.inter" = 18,
             "female:baby.inter" = 11,
             "male:old.growth" = 2,
             "male:young.growth" = 1,
             "male:baby.growth" = -4,
             "female:old.growth" = 2,
             "female:young.growth" = 0,
             "female:baby.growth" = -1)

################## input end ##################

SimulatedData <- function(param, n){
  gdata <- data.frame(
    group1 = rep(c("male", "female"), each=n/2),
    group2 = sample(c("old", "young", "baby"), n, replace=TRUE)
  )
  gdata$HI<- round(runif(n), 2)
  xloads <- by(
    gdata, 
    gdata$group1:gdata$group2, 
    function (x) {
      pattern <- paste0("^", unique(x$group1), ":", unique(x$group2))
      this.param <- param[grepl(pattern, names(param))]
      loads <- rnbinom(
        n = nrow(x), 
        size = param["k"],
        mu = MeanLoad(
          intercept = this.param[grepl("\\.inter", names(this.param))],
          growth = this.param[grepl("\\.growth", names(this.param))],
          alpha = this.param[grepl("\\.alpha", names(this.param))],
          hybridIndex = x$HI
        )
      )
      cbind(x, loads)
    }
  )
  as.data.frame(do.call("rbind", xloads))
}

set.seed(5)
simdata <- SimulatedData(simpara, 1000)

################## simulation end ##################

################## Test : one discrete group OK ##################

G1 <- glm.hybrid(loads ~ HI * group1, data = simdata)

G2 <- glm.hybrid(loads ~ HI * group1 * group2, data = simdata)

Joelle_data <- read.csv("../examples/Reproduction_WATWM/EvolutionFinalData.csv")
Joelle_data[is.na(Joelle_data)] <- 0
is.na(Joelle_data)

G3 <- glm.hybrid(Trichuris ~ HI * Sex, data = Joelle_data)

glm.hybrid(Aspiculuris.Syphacia ~ HI * Sex, data = Joelle_data)

glm.hybrid(Taenia ~ HI * Sex, data = Joelle_data)

################## Test : two discrete groups OK ##################

G4 <- glm.hybrid(loads ~ HI * group1 * group2, data = simdata)

################## Test : starting parameters OK ##################

# really bad when starting parameters just close to zero
opt.para <- glm.hybrid(formula=loads ~ HI*group1*group2, data=simdata,startValues=simpara)

# ADD AGAIN WHEN / IF ALPHA.START IS ADDED BACK
# glm.h1 <- glm.hybrid(formula=loads~HI*group1*group2, data=simdata, alpha.along = "HI",
#                                  alpha.start=1)
# 
# glm.h1.5 <- glm.hybrid(formula=loads~ HI*group1*group2, data=simdata, alpha.along = "HI",
#                                    alpha.start=1.5)
# 
# glm.h1.9 <- glm.hybrid(formula=loads~ HI*group1*group2, data=simdata, alpha.along = "HI",
#                                    alpha.start=1.9)
# 
# glm.h2.5 <- glm.hybrid(formula=loads~HI*group1*group2, data=simdata, "HI",
#                                    alpha.start=2.5)

# para.table <- cbind(simpara,
#                     opt.sim = opt.para$opt.param[names(simpara)],
#                     opt.nb1 = glm.h1$opt.param[names(simpara)],
#                     opt.nb1.5 = glm.h1.5$opt.param[names(simpara)],
#                     opt.nb1.9 = glm.h1.9$opt.param[names(simpara)],
#                     opt.nb2.5 = glm.h2.5$opt.param[names(simpara)])
# 
# glm.h1$opt.param
# opt.para
# names(simpara)
# para.table
# 
# opt.para$value
# glm.h1$twologlik/2
# glm.h1.5$twologlik/2
# glm.h1.9$twologlik/2
# glm.h2.5$twologlik/2
# ADD AGAIN WHEN / IF ALPHA.START IS ADDED BACK

## replace some of the parameters to not come via glm.nb (start.mod)
## but being entered manually (in a named vector, names have to be
## same as the model would assign)
# glm.try <- glm.hybrid(formula=loads~HI*group1, data=simdata, "HI",
#                                               alpha.start=2.5, start.values=simpara[5:8])
# 
# non.nb <- glm.hybrid(formula=loads~HI*group1, data=simdata, "HI",
#                                  alpha.start=1.5, start.mod = MASS::glm.nb)
# 
# # Compare models
# glm.h0 <- glm.hybrid(formula=loads~HI*group1, data=simdata, alpha.along = "HI",
#                      alpha.start=1)
# 
# anova.hybrid(m1 = glm.h1, m2 = glm.h0)

# ## Plot results for one group
# HI = seq(0,1,0.001)
# df <- data.frame(HI = HI,
#                  ML = MeanLoad(G1$opt.param["female.inter"], G1$opt.param["female.growth"],
#                                G1$opt.param["alpha"], HI))
# 
# ggplot2::ggplot(df, aes(x = HI, y = ML)) +
#   geom_point() +
#   theme_classic()

# ## Correct : 1 k nd alpha per group!!
# ## By row
# LogLik <- function(data, param, group.name, response, alpha.along, whichsign = 1){
#   ## split the name into two
#   gname <- sort(group.name)
#   split.L<- by(data, data[, gname], function(x)  {
#     ## by makes sure we get all levels: get the name of the paramter
#     ## from the values within the by "loop"
#     param.pattern <- unique(interaction(x[, gname], sep=":"))
#     ## construct a regex with it
#     par.regex <- paste0("^k$|alpha|^", param.pattern)
#     ## select from our ugly paramter collection
#     sub.param <- param[grepl(par.regex, names(param))]
#     l.lik <- dnbinom(x[, response],
#                      size=abs(sub.param[names(sub.param) %in% "k"]),
#                      mu=abs(MeanLoad(alpha=sub.param[names(sub.param) %in% "alpha"],
#                                      intercept=sub.param[grepl("inter",
#                                                                names(sub.param))],
#                                      growth=sub.param[grepl("growth", names(sub.param))],
#                                      HI=x[, alpha.along])),
#                      log = TRUE)
#     l.lik
#   })
#   all.l.lik <- unlist(split.L)
#   if(length(all.l.lik)!=nrow(data)){
#     stop("Not all likelihoods considered, group/parameter matching problem")
#   } else{
#     sum(all.l.lik) * whichsign
#   }
# }
