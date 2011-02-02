if (FALSE) {

# Some useful circular functions
as.cir<-function(x){as.circular(x,type = "angles",units ="radians",template = "none", modulo ="asis", zero = 0, rotation = "counter")}
dvon<-Vectorize(function(x,mu,kappa){dvonmises(as.cir(x),as.cir(mu),kappa)},"mu")
pvon<-Vectorize(function(x,mu,kappa){pvonmises(as.cir(x),as.cir(mu),kappa)},"mu")
qvon<-Vectorize(function(p,mu,kappa){uniroot(function(x){pvon(x,mu,kappa)-p},interval=c(0+1e-6,2*pi-1e-6),tol=1e-6,maxiter=1000)$root},c("p","mu"))

# Kernel circular density
dkernel_circular<-function(x,sample,h){
	density.circular(z=circular(x),x=circular(sample),bw=h,control.circular=list(type="angles",units="radians",template="none",modulo="asis",zero=0,rotation="counter"))$y
}

# Kernel circular distribution
# ttheta: angles where dk=dkernel_circular is computed to gain speed
pkernel_circular<-function(x,sample,h,ttheta=NULL,dk=NULL){
 	# sapply(x,function(x){sum(pvon(x,mu=sample,kappa=h))/length(sample)}) #not efficient at all
	 if(is.null(dk) || is.null(ttheta)){
		f<-splinefun(seq(0,2*pi,by=1e-3),dkernel_circular(seq(0,2*pi,by=1e-3),sample,h))
	 } else {
		f<-splinefun(ttheta,dk)
	 }
	 F<-Vectorize(function(x){integrate(f,0,x)$value},"x")
	 return(F(x))
}

# Cross validation by MSE
mse=function(h,samp){
	t1<-integrate(function(x){dkernel_circular(x,sample=samp,h)^2},lower=0,upper=2*pi,abs.tol=1e-4)$value
	r<-sapply(1:length(samp),function(i){dkernel_circular(x=samp[i],sample=samp[-i],h)})
	t2<-2*sum(r)/length(samp)
	return(t1-t2)
}

cross_mse=function(sample,interval){optimize(function(h){mse(h,sample)},interval=interval,tol=1e-4)}

# Cross validation by ML
ml=function(h,samp){
	s<-sapply(1:length(samp),function(i){log(dkernel_circular(x=samp[i],sample=samp[-i],h))})
	return(1/length(samp)*sum(s))
}

cross_ml=function(sample,interval){optimize(function(h){-ml(h,sample)},interval=interval,tol=1e-4)}


## Some examples ##

# Example 1

sample<-rvonmises(200,circular(pi),10)

bw_mse<-cross_mse(sample,interval=c(0,100))

bw_ml<-cross_ml(sample,interval=c(0,100))

theta<-seq(0,2*pi,by=0.01)

hist(as.numeric(sample)%%(2*pi),freq=F,breaks=seq(0,2*pi,len=30),main="Histogram",xlab="theta")
lines(theta,dkernel_circular(theta,sample,bw_mse$minimum),col="red",lwd=2)
lines(theta,dkernel_circular(theta,sample,bw_ml$minimum),col="blue",lwd=2)
legend("topright",lwd=2,legend=c("MSE","ML"),col=c("red","blue"))
rug(as.numeric(sample)%%(2*pi))


# Example 2

sample<-numeric(200)
sample[1:100]<-rvonmises(100,circular(pi/2),10)
sample[101:200]<-rvonmises(100,circular(3*pi/2),5)

bw_mse<-cross_mse(sample,interval=c(0,100))

bw_ml<-cross_ml(sample,interval=c(0,100))

theta<-seq(0,2*pi,by=0.01)

hist(as.numeric(sample)%%(2*pi),freq=F,breaks=seq(0,2*pi,len=30),main="Histogram",xlab="theta")
lines(theta,dkernel_circular(theta,sample,bw_mse$minimum),col="red",lwd=2)
lines(theta,dkernel_circular(theta,sample,bw_ml$minimum),col="blue",lwd=2)
legend("top",lwd=2,legend=c("MSE","ML"),col=c("red","blue"))
rug(as.numeric(sample)%%(2*pi))

}
