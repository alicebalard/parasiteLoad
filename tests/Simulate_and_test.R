### source or load package once done with packaging
library(devtools)
install_github("alicebalard/Parasite_Load")

simpara <- c(k = 2, alpha = 1.92,
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

################### input end ###################################

#########################    
SimulatedData <- function(param, n){
  gdata <- data.frame(group1 = rep(c("male", "female"), each=n/2),
                      group2 = sample(c("old", "young", "baby"),
                                      n, replace=TRUE))
  gdata$HI<- round(runif(n), 2)
  xloads <- by(gdata, gdata$group1:gdata$group2, function (x) {
    pattern <- paste0("^", unique(x$group1), ":", unique(x$group2))
    this.param <- param[grepl(pattern, names(param))]
    loads <- rnbinom(n = nrow(x), size = param["k"],
                     mu = glm.hybrid:::MeanLoad(intercept=this.param[grepl("\\.inter",
                                                                           names(this.param))],
                                                growth=this.param[grepl("\\.growth",
                                                                        names(this.param))],
                                                alpha=param["alpha"],
                                                HI=x$HI))
    cbind(x, loads)
  })
  as.data.frame(do.call("rbind", xloads))
}

set.seed(5)
simdata <- SimulatedData(simpara, 1000)

glm.hybrid:::LogLik(simdata, simpara, c("group1", "group2"))

### really bad when starting paramters just close to zero
### parameters
opt.para <- glm.hybrid::glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                                   alpha.start=1, start.values=simpara)

glm.h1 <- glm.hybrid::glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                                 alpha.start=1)

glm.h1.5 <- glm.hybrid::glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                                   alpha.start=1.5)

glm.h1.9 <- glm.hybrid::glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                                   alpha.start=1.9)

glm.h2.5 <- glm.hybrid::glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                                   alpha.start=2.5)

## replace some of the parameters to not come via glm.nb (start.mod)
## but being entered manually (in a named vector, names have to be
## same as the model would assign)
glm.try <- glm.hybrid::glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                                  alpha.start=2.5, start.values=simpara[5:8])

non.nb <- glm.hybrid::glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                                 alpha.start=1.5, start.mod=glm.nb)

para.table <- cbind(simpara,
                    opt.sim = opt.para$par[names(simpara)],
                    opt.nb1 = glm.h1$opt.param[names(simpara)],
                    opt.nb1.5 = glm.h1.5$opt.param[names(simpara)],
                    opt.nb1.9 = glm.h1.9$opt.param[names(simpara)],
                    opt.nb2.5 = glm.h2.5$opt.param[names(simpara)])

opt.para$value
glm.h1$twologlik/2
glm.h1.5$twologlik/2
glm.h1.9$twologlik/2
glm.h2.5$twologlik/2

## pairs(para.table)

glm.hybrid::glm.hybrid(formula=loads~group2*HI, data=simdata, "HI")$twologlik/2

glm.hybrid::glm.hybrid(formula=loads~group1*HI, data=simdata, "HI")$twologlik/2

NBglm <- glm.hybrid::glm.nb(formula=loads~group1*group2*HI, data=simdata)
NBglm$twologlik/2

######################################################################################
## ML_bounds tests:
glm.hybrid(formula=loads~group1*HI*group2, data=simdata, "HI")

# Generate random simdata:  (based on Phoung's oocysts counting)
simdata_generator <- function(){
  I <- round(runif(6, 0, 4*10^6), 2)
  S <- round(runif(6, -I, 10^6), 2)
  simparaBS <- c(k = round(abs(runif(1, 1, 8)), 2),
                 alpha = round(runif(1, -2, 2), 2),
                 "male:old.inter" = I[1],
                 "male:young.inter" = I[2],
                 "male:baby.inter" = I[3],
                 "female:old.inter" = I[4],
                 "female:young.inter" = I[5],
                 "female:baby.inter" = I[6],
                 "male:old.growth" = S[1],
                 "male:young.growth" = S[2],
                 "male:baby.growth" = S[3],
                 "female:old.growth" = S[4],
                 "female:young.growth" = S[5],
                 "female:baby.growth" = S[6])
  simdataBS <- SimulatedData(simparaBS, 1000)
  # output:
  list(simparaBS, simdataBS)
}

# Bootstrap the CI:
myBS <- function(){
  simBS <- simdata_generator()
  simpara <- simBS[[1]]
  simdata <- simBS[[2]]
  # Calculate bounds 95%CI:
  bounds <- ML_bounds_Wald(param = simpara, data = simdata,
                           group.name =  c("group1", "group2"))
  
  # Check if the simpara are indeed in the good CI:
  bounds <- as.data.frame(bounds)
  bounds$actualpara <- simpara
  (sum(bounds$LowerBounds <= bounds$actualpara) + sum(bounds$UpperBounds >= bounds$actualpara)) /28*100
}

# Give NA instead of error:
myBS_secured <- function(){
  tryCatch(myBS(), error=function(err) NA)
} 

# Big BS:
result <- replicate(1000, myBS_secured())
mean(na.omit(result))
## 97.32