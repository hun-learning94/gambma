
xi=c(0.05,0.1,0.15,0.2,0.5,0.55,0,0.6,0.9,0.95)
K=length(xi)

x=ppoints(1000)

# b1=function(x,j){
#   (pmax(0,(x-xi[j])^3)-pmax(0,(x-xi[K])^3))/(xi[K]-xi[j]) -
#     (pmax(0,(x-xi[K-1])^3)-pmax(0,(x-xi[K])^3))/(xi[K]-xi[K-1])
# } 
# b2=function(x,j){
#   (pmax(0,(x-xi[j])^3)-pmax(0,(x-xi[K])^3))/(xi[K]-xi[j]) -
#     (pmax(0,(x-xi[1])^3)-pmax(0,(x-xi[K])^3))/(xi[K]-xi[1])
# }

b1=function(x,j){
  (abs(x-xi[j])^3-abs(x-xi[K])^3)/(xi[K]-xi[j]) - 
    (abs(x-xi[K-1])^3-abs(x-xi[K])^3)/(xi[K]-xi[K-1])
}
b2=function(x,j){
  (abs(x-xi[j])^3-abs(x-xi[K])^3)/(xi[K]-xi[j]) - 
    (abs(x-xi[1])^3-abs(x-xi[K])^3)/(xi[K]-xi[1])
}


W1=matrix(NA_real_,length(x),K)
W1[,1]=1
W1[,2]=x
for(j in 1:(K-2)) W1[,j+2]=b1(x,j)

W2=matrix(NA_real_,length(x),K)
W2[,1]=1
W2[,2]=x
for(j in 1:(K-2)) W2[,j+2]=b2(x,j+1)

W1
W2

P1=W1%*%solve(crossprod(W1))%*%t(W1)
P2=W2%*%solve(crossprod(W2))%*%t(W2)

P1;P2

sum(abs(P1-P2))
all.equal(P1, P2)

y=rnorm(1000,0,100)

max(abs(P1%*%y-P2%*%y))
