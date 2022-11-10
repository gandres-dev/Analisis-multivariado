# ---------------------------------------------------------------------------
# ----------------------------------------------------------------------------




xdat  = read.table("cars.dat")
#delete first column
x = xdat[, -1]

#define variable names
colnames(x) = c("economy", "service", "value", "price", "design", "sporty",
  "safety", "easy")

#x  = as.matrix(x1)
#y  = t((t(x)-mean(x))/sd(x))

r  = cor(x)



###########################
#performs factor analysis using principal factors method
#returns the un-rotated loadings

factpf = function(r){
  n         = ncol(r)
  redcormat = r
  diag(redcormat) = apply(abs(r-diag(1,nrow=n,ncol=n)),2,max)
  xeig      = eigen(redcormat)
  xval      = xeig$values
  xvec      = xeig$vectors
  for (i in 1:length(xval[xval>0])){
    xvec[,i]= xvec[,i]*sqrt(xval[i])}
  rownames(xvec) = colnames(r)
  return(xvec[,xval>0])
}

########################################

#performs factor analysis using iterated principal factors method
#returns the communalties after each iteration and the final un-rotated loadings.

factiter = function(r,niter=5,maxfactors=ncol(r)){
  n         = ncol(r)
  temp      = matrix(0,nrow=n,ncol=n)
  comm      = matrix(0,nrow=niter+1,ncol=n)
  y         = factpf(r)
  m         = ncol(y)
  temp[1:n,1:m] = y
  comm[1,]  = apply(temp^2,1,sum)                 
  for (i in 2:(niter+1)){
    redcormat = r
    diag(redcormat) = comm[i-1,]
    xeig    = eigen(redcormat)
    m = min(maxfactors,length(xeig$values[xeig$values>0]))
    for (j in 1:m){
      xeig$vectors[,j] = xeig$vectors[,j]*sqrt(xeig$values[j])} 
    temp[1:n,1:m] = xeig$vectors[1:n,1:m]
    comm[i,]= apply(temp[1:n,1:m]^2,1,sum)
  if(i ==niter+1){
  	pp = sum(xeig$values[1:m])/ncol(r)
  	cat("proportion  variance explained by",m,"factors is=",pp,"\n")
  }
  }
  f.loadings= temp[1:n,1:m]
  rownames(f.loadings) = colnames(r)
  return(list(comm=comm,f.loadings=f.loadings))
}

 
##############################################
PFA         = factiter(r,niter=5,maxfactor=2)
#factor analysis - iterated principal factors method (PFM)
#number of iterations: 10, maximal number of factors to be extracted: 2
loadings    = PFA$f.loadings
colnames(loadings) = c("Q1","Q2")
communality = PFA$comm
specific_variance = 1-PFA$comm
#b           = cbind(loadings,communality,specific_variance)
#round(b,4)  #iterating the PFM 10 times
# xtable(b,digits=4)
#varim = varimax(loadings) #rotating loadings using function varimax
#print(b)

#windows(width=8,height=8)
opar        = par(mfrow=c(1,2))
lab         = colnames(x) #adding labels

#cat("original loadings (columns of Q)")
#print(loadings)
#commu= apply(loadings^2,1,sum)
#print(commu)
vl   = cbind(loadings[,1],loadings[,2])
com  = diag(vl%*%t(vl))  
psi  = diag(r)-diag(vl%*%t(vl))
tbl  = cbind(loadings[,1],loadings[,2],com,psi)
cat("\n")
cat("Factor Analysis","\n")
print(tbl)


plot(c(-1.1,1.1),c(-1.1,1.1),type="n",main="Car means (loadings)",xlab=expression(q[1]),ylab=expression(q[2])) #plotting... [KONECNE!]
ucircle     = cbind(cos((0:360)/180*3.14159),sin((0:360)/180*3.14159))
points(ucircle,type="l",lty="dotted")
abline(h = 0)
abline(v = 0)
text(loadings[,1:2],labels=lab,col="black",xpd=NA)

#manual rotation NOT VARIMAX
#theta       = 3.14159*7/12
theta = 5.0
rot         = rbind(c(cos(theta),sin(theta)),c(-sin(theta),cos(theta)))
qr          = loadings%*%rot
#cat("\n")
#cat("loadings (columns of Q) after rotation")
#print(qr)
#commu= apply(qr^2,1,sum)
#print(commu)
vl   = cbind(qr[,1],qr[,2])
com  = diag(vl%*%t(vl))  
psi  = diag(r)-diag(vl%*%t(vl))
tbl  = cbind(qr[,1],qr[,2],com,psi)
cat("\n")
cat("Factor Analysis after manual rotation","\n")
print(tbl)

# Now Apply varimax rotation 
var  = varimax(loadings)                #rotates the factor loadings matrix
load = var$loadings              #estimated factor loadings after varimax
vl   = cbind(load[,1],load[,2])
com  = diag(vl%*%t(vl))           #communalities are calculated
psi  = diag(r)-diag(vl%*%t(vl)) #specific variances are calculated
tbl  = cbind(load[,1],load[,2],com,psi)
cat("\n")
cat("Factor Analysis after varimax rotation","\n")
#cat("loadings","\n")
#print(load)
#cat("commun","\n")
#print(com)
#cat("psi","\n")
#print(psi)
print(tbl)



#windows(width=8,height=8)

plot(c(-1.1,1.1),c(-1.1,1.1),type="n",main="Car means (rotation)",xlab=expression(q[1]),ylab=expression(q[2])) 
points(ucircle,type="l",lty="dotted")
abline(h = 0)
abline(v = 0)
text(qr[,1:2],labels=lab,col="black",xpd=NA)

par(opar)
