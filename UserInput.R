


## Get it somehow from the formula
fo <- loads~ HI*group1*group2

## the holy grail of how to set this up
MMat <- model.matrix(object=fo, data=fakedata)

## or we could extract something from here
terms(fo, data=fakedata)

## this does not seem helpful
get_all_vars(fo, fakedata)

## why does this not help?
with(fakedata, levels(group1:group2))
