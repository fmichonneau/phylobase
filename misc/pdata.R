setClass("pdata", representation(x="vector", y="vector"))
setMethod("[","pdata",function(x,i, j,...,drop=TRUE)new("pdata",x=x@x[i],y=x@y[i]))

# x <- new("pdata", x=c("a","b", "c", "d", "3"), y=c(1:5))
#>x[c(2,4)]
#An object of class “pdata”
#Slot "x":
#[1] "b" "d"
#
#Slot "y":
#[1] 2 4



# doesn't work
#setClass("track", representation("list", comment="character", metadata="vector"), contains="list", prototype(list(), comment="", metadata=NA))
#setMethod("[","track",function(x,i, j,...,drop=TRUE)new("track", list(lapply(x, function(x, i, j, ..., drop=TRUE) x@.Data[i]))))

# this works, how to incorporate into method above?
#> lapply(x, function(x, i=2, j, ..., drop=TRUE) x@.Data[i])
#$x
#[1] "b"

#$y
#[1] 2

# this works, but list structure is destroyed
#> mapply(function(x, i, j, ..., drop=TRUE) x@.Data[i], x, 2)
#  x   y 
#"b" "2" 

