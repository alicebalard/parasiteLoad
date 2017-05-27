## From AliceTest
## source stuff here  (bugged??)

## Ideas for CI

AliceTest

library(ggplot2)
library(plyr)

## For each UB (and LB) fixed, we want to run a maximum likelihood optimisation:
LogLik.1Pfixed <- function(data, allbut1P, oneP, group.name){
    LogLik(data, c(allbut1P, oneP), group.name)
}

myfn1 <- function(i, param, dataMLbounds, group.name){
    oneP <- dataMLbounds[[i]]$UB.minimum
    names(oneP) <- names(param[i])
    return(c(oneP,
             optim(par = param[-i],
                   fn = LogLik.1Pfixed, ## function to be maximized
                   control = list(fnscale=-1),
                   data = simdata,
                   oneP = oneP,
                   group.name = group.name)$par))
}

## Run over all parameters:
set.of.param <- lapply(1:2, function(i) { ## length(names(simpara))
    myfn1(i,
          param = simpara,
          dataMLbounds = AliceTest,
          group.name = "group1:group2")})

## Calculate the MeanLoad for all our points in dataset,
## for the set.of.parameters:
set.of.param[[1]]
param[["group1:group2"]]

## For 1 set of parameters:
DF <- data.frame(set.of.param[[1]])

DF2 <- data.frame(parameters = sapply(1:nrow(DF), function(i) strsplit(row.names(DF), "\\.")[[i]][1]),
           level = sapply(1:nrow(DF), function(i) strsplit(row.names(DF), "\\.")[[i]][2]),
           value = DF[[1]])

by(DF2, level, function(x){DF2})


myfn2 <- function(data, group.name, param){
    ## split the name into two
    gname <- strsplit(group.name, ":")[[1]]
    ## by makes sure we get all levels: get the name of the paramter
    split.L <- ddply(data, data[, gname], function(x)  {
        data[x]})
}

myfn2(data = simdata,
      group.name = "group1:group2",
      param = set.of.param[[1]])

split.L <- ldply(data, data[, gname], function(x)  {
        unique(interaction(x[, gname], sep=":"))
        split.L
    

        ## from the values within the by "loop"
        param.pattern <- unique(interaction(x[, gname], sep=":"))
        param.pattern
        ## construct a regex with it
        par.regex <- paste0("^k$|alpha|\\.", param.pattern)
        ## select from our ugly paramter collection
        sub.param <- param[grepl(par.regex, names(param))]
        MeanLoad(intercept = sub.param[[grep("int", x = names(sub.param))]],
                 growth = sub.param[[grep("growth", x = names(sub.param))]],
                 alpha = sub.param[["alpha"]],
                 HI = seq(0, 1, 0.2))
    })
    split.L




######








head(simdata)
set.of.param
MeanLoad

myfn2 <- function(data, group.name, param){
  ## split the name into two
    gname <- strsplit(group.name, ":")[[1]]
    split.L<- by(data, data[, gname], function(x)  {
        ## by makes sure we get all levels: get the name of the paramter
        ## from the values within the by "loop"
        param.pattern <- unique(interaction(x[, gname], sep=":"))
        ## construct a regex with it
        par.regex <- paste0("^k$|alpha|\\.", param.pattern)
        ## select from our ugly paramter collection
        sub.param <- param[grepl(par.regex, names(param))]
        MeanLoad(intercept = sub.param[[grep("int", x = names(sub.param))]],
                 growth = sub.param[[grep("growth", x = names(sub.param))]],
                 alpha = sub.param[["alpha"]],
                 HI = seq(0, 1, 0.2))
    })
    split.L
}

wtest <- myfn2(data = simdata,
      group.name = "group1:group2",
      param = set.of.param[[1]])

test <- lapply(1:2, function(i) myfn2(data = simdata, group.name = "group1:group2",
      param = set.of.param[[i]]))


max(test[[1]][[1]][1], test[[2]][[1]][1])

sapply(1:length(test[[1]][[1]]), function(i) max(test[[1]][[1]][i],test[[2]][[1]][i]))

data.frame(
    do.call("rbind", by(test, column, mean))
)

data.frame(test[[1]])

names(test[[1]])

simdata[1,]

head(simdata)
data.frame(test)
    

    
  
      ## select from our ugly paramter collection
sub.param <-
    param[grepl(par.regex, names(param))]
  

MeanLoad(alpha=sub.param[names(sub.param) %in% "alpha"],
         intercept=sub.param[grepl("int\\.", names(sub.param))],
         growth=sub.param[grepl("growth\\.", names(sub.param))],
         HI=x$HI)

simpara



## Plot:

DF <- data.frame(est = MeanLoad(intersect= extract.fun(hypo$ML.estimated.parameters, "intersect", gpe),
                                growth = extract.fun(hypo$ML.estimated.parameters, "growth", gpe),
                                alpha = hypo$ML.estimated.parameters[names(hypo$ML.estimated.parameters) == "alpha"],                                    HI = HIplot),
                 LB=MeanLoad(intersect= extract.fun(hypo$ML.lower.parameters, "intersect", gpe),
                             growth = extract.fun(hypo$ML.lower.parameters, "growth", gpe),
                             alpha = hypo$ML.lower.parameters[names(hypo$ML.lower.parameters) == "alpha"],                                    HI = HIplot),
                 UB=MeanLoad(intersect= extract.fun(hypo$ML.upper.parameters, "intersect", gpe),
                             growth = extract.fun(hypo$ML.upper.parameters, "growth", gpe),
                             alpha = hypo$ML.upper.parameters[names(hypo$ML.upper.parameters) == "alpha"],                                    HI = HIplot))
names(DF) <- paste0(names(DF), ".", gpe)
DF
}

M <- cbind(dfplot.func(H3, "male"), dfplot.func(H3, "female"), HIplot)

ggplot(M)+
    geom_line(aes(x=HIplot, y=est.male), color="blue")+
    geom_line(aes(x=HIplot, y=est.female), color="red")+
    geom_ribbon(aes(x=HIplot, ymax=UB.male, ymin=LB.male), fill= "blue", alpha=.3)+
    geom_ribbon(aes(x=HIplot, ymax=UB.female, ymin=LB.female), fill= "red", alpha=.3)+
    geom_point(data = alicedata[alicedata$group == "male",], aes(x=HI, y=loads), color="blue", size=0.5)+
    geom_point(data = alicedata[alicedata$group == "female",], aes(x=HI, y=loads), color="red", size=0.5)+
    theme_bw()+ theme(legend.position="none")

