
library(mvQuad)
library(mvtnorm)

Quadrature = function( myGrid, f, ... ){
  Nodes = getNodes(myGrid)
  Weights = getWeights(myGrid)
  # Compute integral
  Integral = rep(NA, nrow(Nodes))
  for( i in seq_along(Integral) ){
    Integral[i] = f( Nodes[i,], ... ) * Weights[i]
  }
  return( sum(Integral) )
}

Sampling = function( dim, min=-Inf, max=Inf, SampVar, f, n=1000, ... ){
  Nodes = matrix( runif(n=n*dim, min=min, max=max), ncol=dim )
  Weights = rep( (max-min)^dim/n, length=n )
  # Compute integral
  Integral = rep(NA, nrow(Nodes))
  for( i in seq_along(Integral) ){
    Integral[i] = f( Nodes[i,], ... ) * Weights[i]
  }
  return( sum(Integral) )
}

#################
# Package example
#################

myGrid <- createNIGrid(dim=2, type="GLe", level=5)
rescale(myGrid, domain=rbind(c(-1,1),c(-1,1)))
getNodes(myGrid)
getWeights(myGrid)

print(myGrid)
plot(myGrid, col="blue")
myFun <- function(x){
   1 - x[,1]^2 * x[,2]^2
}
quadrature(myFun, myGrid)

#################
# Rescale example
#################

C = matrix(c(2,0.9,0.9,2),2)
m = c(-.5, .3)
par(mfrow=c(1,3), mar=c(2,2,0,0) )

myGrid <- createNIGrid(dim=2, type="GHe", level=5,  ndConstruction="sparse")

rescale(myGrid, m=m, C=C, dec.type=0)
plot(myGrid, col="red")

rescale(myGrid, m=m, C=C, dec.type=1)
plot(myGrid, col="green")

rescale(myGrid, m=m, C=C, dec.type=2)
plot(myGrid, col="blue")


#################
# Custom example
#################
Mean = c(1,1)
Cov = matrix(c(1,0.5,0.5,1),ncol=2)
f = dmvnorm

myGrid <- createNIGrid(dim=2, type="GHe", level=5,  ndConstruction="sparse")
rescale(myGrid, m=Mean, C=Cov, dec.type=1)

plot(myGrid)
quadrature(f=f, mean=Mean, sigma=Cov, myGrid)
Quadrature(f=f, mean=Mean, sigma=Cov, myGrid)
Sampling(f=f, mean=Mean, sigma=Cov, dim=2, min=-2, max=4, n=10000)


################
# Calculate exepctation for a 2D function
################
Mean = c(1,1)
Cov = matrix(c(1,0.5,0.5,1),ncol=2)
f = function(x, mean, sigma){
  Like = dmvnorm(x, mean=mean, sigma=sigma)
  Val = sum(x^2) + x[1]
  return( Like*Val )
}

myGrid <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
rescale(myGrid, m=Mean, C=Cov, dec.type=1)

# Explore quadrature grid
plot(myGrid)
print(myGrid)

# Compare methods
Quadrature(f=f, mean=Mean, sigma=Cov, myGrid)
Sampling(f=f, dim=2, mean=Mean, sigma=Cov, min=-2, max=4, n=10000)


################
# Calculate exepctation for a 3D function
################
Mean = c(1,1,1)
Cov = matrix(c(1,0.5,0.2,0.5,1,0.3,0.2,0.3,1),ncol=3)
f = function(x, mean, sigma){
  Like = dmvnorm(x, mean=mean, sigma=sigma)
  Val = sum(x^2) + x[1]
  return( Like*Val )
}

myGrid <- createNIGrid(dim=3, type="GHe", level=4,  ndConstruction="sparse")
rescale(myGrid, m=Mean, C=Cov, dec.type=1)

# Explore quadrature grid
plot(myGrid)
print(myGrid)

# Compare methods
Quadrature(f=f, mean=Mean, sigma=Cov, myGrid)
Sampling(f=f, dim=3, mean=Mean, sigma=Cov, min=-2, max=4, n=10000)








