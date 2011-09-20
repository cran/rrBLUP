distance <- function(G) {
#assumes genotypes on [-1,1] scale
#can handle fractional values
#missing alleles omitted

n <- nrow(G) #number of lines

D <- 0*diag(n)
for (i in 1:(n-1)) {
for (j in (i+1):n) {
D[i,j] <- sqrt(mean((G[i,]-G[j,])^2,na.rm=TRUE))/2
D[j,i] <- D[i,j]
}
}
D
}
