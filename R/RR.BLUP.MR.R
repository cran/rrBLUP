RR.BLUP.MR <-
function(y,G,Z,X = NULL, K = NULL, Nfold = 3, NumReduced = 8, MinMarkers = 10, SE="FALSE") {

n <- length(y)
if (is.null(X)) {
  p <- 1
  X <- matrix(rep(1,n),n,1)
}
p <- dim(X)[2]  
if (is.null(p)) {
  p <- 1
  X <- matrix(X,length(X),1)
}

rX <- qr(X)$rank
stopifnot(rX==p)  #must be full rank design matrix
stopifnot(dim(X)[1]==n)

m <- dim(G)[2]  # number of markers
t <- dim(G)[1]

stopifnot(dim(Z)[1]==n)
stopifnot(dim(Z)[2]==t)

if (is.null(K)) {
  K <- G%*%t(G)/m
}

stopifnot(nrow(K)==ncol(K))
stopifnot(nrow(K)==t)

NumMarkers <- round(exp(seq(log(MinMarkers),log(m),length.out=NumReduced)))

#Partition lines for cross-validation
sets <- rep(0,t)
ix <- sort(runif(t),index.return=TRUE)$ix
SetSize <- floor(t/Nfold)
for (k in 1:(Nfold-1)) {sets[ix[((k-1)*SetSize+1):(k*SetSize)]] <- k}
sets[ix[((Nfold-1)*SetSize+1):t]] <- Nfold

CorrResults <- matrix(rep(0,Nfold*NumReduced),NumReduced,Nfold)

for (j in 1:Nfold)  {
  test.lines <- which(sets==j)
  train.lines <- setdiff(1:t,test.lines)
  test <- NULL
  for (line in test.lines) {
     test <- union(test,which(Z[,line]==1))
  }
  train <- setdiff(1:n,test)
  n.train <- length(train)
  n.test <- length(test)

  Scores <- GWA(y[train],G=G[train.lines,],Z=Z[train,train.lines],X=X[train,],K=K[train.lines,train.lines])
  ScoreIdx <- sort(Scores,decreasing=TRUE,index.return=TRUE)$ix

  for (i in 1:NumReduced) {
   MarkerList <- ScoreIdx[1:NumMarkers[i]]
   if (NumMarkers[i] < t) {
     #markers as random effects
     soln <- RR.BLUP(y[train],X=X[train,],Z=Z[train,train.lines]%*%G[train.lines,MarkerList],K=diag(NumMarkers[i]))
     Ypred <- X[test,]%*%soln$beta + Z[test,test.lines]%*%G[test.lines,MarkerList]%*%soln$u
   } else {
     #lines as random effects
     Zaug <- cbind(Z[train,train.lines],matrix(rep(0,n.train*length(test.lines)),n.train,length(test.lines)))
     Gaug <- rbind(G[train.lines,MarkerList],G[test.lines,MarkerList])
     soln <- RR.BLUP(y[train],X=X[train,],Z=Zaug,K=Gaug%*%t(Gaug)/length(MarkerList))
     Ypred <- X[test,]%*%soln$beta + Z[test,test.lines]%*%soln$u[(length(train.lines)+1):t]
   }
   CorrResults[i,j] <- cor.test(y[test],Ypred,method=c("pearson"))$estimate
  } #end i

} #end j 
  
 MeanCorr <- rowMeans(CorrResults)
 #find maximum 
 maxIdx <- which.max(MeanCorr)
 Scores <- GWA(y,G=G,Z=Z,X=X,K=K)
 ScoreIdx <- sort(Scores,decreasing=TRUE,index.return=TRUE)$ix
 MarkerList <- ScoreIdx[1:NumMarkers[maxIdx]]
 soln <- RR.BLUP(y,X=X,Z=Z%*%G[,MarkerList],K=diag(NumMarkers[maxIdx]),SE=SE)
 u <- rep(0,m)
 u[MarkerList] <- soln$u
 
 if (SE == "FALSE") {
   list(CV=cbind(NumMarkers,MeanCorr),beta=soln$beta,u=u)
 } else {
   u.SE <- rep(0,m)
   u.SE[MarkerList] <- soln$u.SE
   list(CV=cbind(NumMarkers,MeanCorr),beta=soln$beta,beta.SE=soln$beta.SE,u=u,u.SE=u.SE)
 }
} #end function

