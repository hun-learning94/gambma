W=diag(rnorm(1000)^2)
 
X=matrix(rnorm(500000),1000,500)
 
A=(as.vector(t(rep(1,1000))%*%W%*%rep(1,1000)))^(-1)*rep(1,1000)%*%t(rep(1,1000))%*%W%*%X
B=colSums(diag(W)/sum(diag(W))*X)
 
max(abs(A[1,]-B))
 
library(microbenchmark)
microbenchmark(colSums(diag(W)/sum(diag(W))*X),colMeans(X))
