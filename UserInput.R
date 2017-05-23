


## Get it somehow from the formula
fo <- loads~ HI + group1*group2

## the holy grail of how to set this up
MMat <- model.matrix(object=fo, data=fakedata)

## This will get difficult
full.MMat <- model.matrix(fo, data=fakedata,
                          contrasts.arg=
                              list(group1=diag(nlevels(fakedata$group1)),
                                   group2=diag(nlevels(fakedata$group2))))


## or we could extract something from here
terms(fo, data=fakedata)

## this does not seem helpful
get_all_vars(fo, fakedata)

## why does this not help?
with(fakedata, levels(group1:group2))
