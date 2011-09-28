impute <- function(G, n.neighbor = "ALL", D = NULL) {

mean.rmNA <- function(x) {
	mean(x,na.rm=TRUE)	
}

n <- nrow(G) #number of lines

if (toupper(n.neighbor) == "ALL") {
	#impute with mean of population
	pop.mean <- apply(G,2,mean.rmNA)
	for (i in 1:n) {
  		missing <- which(is.na(G[i,]))
		G[i,missing] <- pop.mean[missing]
	}
} else {
	if (is.null(D)) {D <- as.matrix(dist(G))}
	
	for (i in 1:n) {
		neighbors <- sort(D[i,-i],index.return=TRUE)$ix[1:n.neighbor]
  		missing <- which(is.na(G[i,]))
		G[i,missing] <- apply(matrix(G[neighbors,missing],n.neighbor,length(missing)),2,mean.rmNA)
	}
}
if (length(which(is.nan(G))) > 0) {stop("Use more neighbors.")}
G
}
