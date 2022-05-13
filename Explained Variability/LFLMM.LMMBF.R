emn <-
function(X){
	N=ncol(X)
	geta = apply(X,1,max)-10
	X = X-geta
	print(range(-geta))
	#x[x==Inf&!is.na(x)]=max(x[x<Inf&!is.na(x)],na.rm=T)
	p=c(0.9,rep(0.1/N,N))
	lkhd=0
	lkhd0=-10000000000
	for(itr in 1:1000){
		d=c(exp(-geta)*p[1]+exp(X)%*%p[-1])
		z=cbind(exp(-geta)*p[1],t(t(exp(X))*p[-1]))/d
		p=c(apply(z,2,mean,na.rm=T))
		p=p/sum(p)
		lkhd=mean(log(d),na.rm=T)
		if(lkhd>lkhd0 && abs(lkhd-lkhd0)<1e-7){cat("Converged.\n");return(z[,c(2:(N+1),1)])}else{print(c(lkhd0,lkhd,p)); lkhd0=lkhd}
	}
}
LMMBF <-
function(Y, res.nlm, colname, ncont=2, Resid=NULL, geta=100, AllClus=F){
	
	logDet = function(x){eval = eigen(x,T)[[1]]; sum(log(eval))}
	
	H = res.nlm$H
	K=ncol(H)
	M=nrow(H)
	
	k = seq(K)[names(res.nlm$nh)==colname]
	Mk= sum(H[,k])
	
	Z = res.nlm$Z
	Zk= Z[,H[,k]==1]
	psi = res.nlm$psi
	psi[k] = psi[k]*geta
	u = c(H%*%psi)
	s2 = res.nlm$sigma2
	yhat = c(res.nlm$X%*%res.nlm$beta)
	# Ytilde = t(t(Y)-yhat) - (res.nlm$a)%*%t(res.nlm$delta)

	ZtZ = t(Z)%*%Z
	A = ZtZ + diag(u^(-2))
	YZ = as.matrix(Y%*%Z) - rep(1,nrow(Y))%*%(t(yhat)%*%Z)
	
	CON = NULL
	DB=NULL
	DBSE=NULL
	if(Mk==1){
		BF = array(0,c(nrow(Y),1))
		DB = array(0,c(nrow(Y),1))
		DBSE = array(0,c(nrow(Y),1))
	}else{
		for(i in 2:ncont){ CON=c(CON, getContrasts2(Mk,i,NULL,NULL)) };   #CON=CON[c(15,2^c(3:0))]
		BF = array(0,c(nrow(Y),length(CON)+1))
		DB = array(0,c(nrow(Y),length(CON)))
		DBSE = array(0,c(nrow(Y),length(CON)))
	}
	Z0 = cbind(Z[,H[,k]==0])
	H0 = H[H[,k]==0,-k]
	ZtZ0 = t(Z0)%*%Z0
	u0 = c(H0%*%psi[-k])
	A0 = ZtZ0 + diag(u0^(-2))
print(dim(A0))
	YZ0 = as.matrix(Y%*%Z0) - rep(1,nrow(Y))%*%(t(yhat)%*%Z0) # Ytilde%*%Z0
	if(!is.null(Resid)){
		resicoef = t(solve(A,t(YZ)))[match(Resid,rownames(Y)),]
		resicoef[,H[,k]==1]=0
        resicoef[,H[,1]==1] = 0
		Re = as.matrix(t(Y[match(Resid,rownames(Y)),])-yhat) - Z%*%t(resicoef)
	}
	M0 = M-Mk
	if(Mk==1){
		BF[,1] = (apply((YZ%*%solve(A))*YZ,1,sum) - apply((YZ0%*%solve(A0))*YZ0,1,sum))/2/s2 - logDet(diag(M)+t(t(ZtZ)*u)*u)/2 + logDet(diag(M0)+(t(t(ZtZ0)*u0)*u0))/2
		DB[,1] = (YZ%*%solve(A))[,H[,k]==1]
		DBSE[,1] = sqrt(solve(A)[H[,k]==1,H[,k]==1]*s2)
	}else{
		BF[,length(CON)+1] = (apply((YZ%*%solve(A))*YZ,1,sum) - apply((YZ0%*%solve(A0))*YZ0,1,sum))/2/s2 - logDet(diag(M)+t(t(ZtZ)*u)*u)/2 + logDet(diag(M0)+(t(t(ZtZ0)*u0)*u0))/2
		for(l in 1:length(CON)){
			print(l)
			H1 = rbind(H[H[,k]==0,], H[H[,k]==1,][seq(ncol(CON[[l]])),])
			u1 = c(H1%*%psi)
			Z1 = cbind(Z[,H[,k]==0],Zk%*%CON[[l]])
			ZtZ1 = t(Z1)%*%Z1
			A1 = ZtZ1 + diag(u1^(-2))
			YZ1 = as.matrix(Y%*%Z1) - rep(1,nrow(Y))%*%(t(yhat)%*%Z1) # Ytilde%*%Z1
			M1 = M-(Mk-ncol(CON[[l]]))
			BF[,l] = (apply((YZ1%*%solve(A1))*YZ1,1,sum) - apply((YZ0%*%solve(A0))*YZ0,1,sum))/2/s2 - logDet(diag(M1)+(t(t(ZtZ1)*u1)*u1))/2 + logDet(diag(M0)+t(t(ZtZ0)*u0)*u0)/2
			B1 = (YZ1%*%solve(A1))[,H1[,k]==1]
			Vg1 = solve(A1)[H1[,k]==1,H1[,k]==1]
			if(ncol(B1)==2){
				DBSE[,l] = sqrt(sum(c(Vg1)*c(1,-1,-1,1))*s2)
				DB[,l]=B1[,2]-B1[,1]
			}
		}
	}
	cat("Posterior mean calculation")
	u = c(H%*%res.nlm$psi)
	A = diag(1/u/u) + ZtZ
	B = t(solve(A,t(YZ)))[,H[,k]==1,drop=F]
	Vb = solve(A)[H[,k]==1,H[,k]==1]
	SE = sqrt(s2)%*%t(sqrt(diag(Vb)))
	
	ltsr = pnorm(0, DB, DBSE)
	ltsr[ltsr<0.5] = 1-ltsr[ltsr<0.5]
	if(Mk>2){if(AllClus){z = emn(BF)}else{z = emn(BF[,1:length(CON)])}}else{z=emn(BF[,1,drop=F])}
	ltsr = ltsr*z[,1:ncol(ltsr)]

	if(ncol(ltsr)>1){
		ltsrname =  unlist(lapply(CON,function(x1,lev){paste(paste(lev[x1[,1]==1],collapse=","),paste(lev[x1[,2]==1],collapse=","),sep="_vs_")},lev=colnames(B)))
		if(ncol(ltsr)>length(CON)){
			ltsrname=c(ltsrname, "All")
		}
		colnames(ltsr) = ltsrname
	}	
	list(bf=BF, beta=B, se=SE, deltabeta=DB, deltabetase=DBSE, contrasts=CON, ltsr=ltsr, z=z, residual=Re)
}
getContrasts2 <-
function(n,k,v,res){
if(length(v)==n){
	if(length(unique(v))==k){if(sum(unique(v)==seq(k))==k){
		tmp = matrix(as.numeric(model.matrix(~0+factor(v))),n)
		res=c(res,list(tmp))
	}}
}else{
	for(i in 1:k){
		res = getContrasts2(n,k,c(v,i),res)
	}
}
return(res)
}
Gost <-
function(res,id=NULL,contid,pos=T,gid=F,th=0.5,thb=1){
	if(is.null(res$contrasts)){
		if(is.na(pos)){
			lab="both"
			flag=res$ltsr[,1]>th
		}else{
			if(pos){
				lab="pos"
				flag=res$ltsr[,1]>th&res$deltabeta[,1]>log(thb)
			}else if(!pos){
				lab="neg"
				flag=res$ltsr[,1]>th&res$deltabeta[,1]<log(1/thb)
			}
		}
	}else{
		if(is.null(id)){id=seq(length(res$contrasts))[unlist(lapply(res$contrasts,function(xx){max(apply(xx==contid,2,sum))}))==length(contid)]}
		if(is.na(pos)){
			flag=res$ltsr[,id]>th
			lab=paste(paste(res$levels[res$contrasts[[id]][,2]>0],collapse="-"),"(both)")
		}else{
			if(pos){
				flag=res$ltsr[,id]>th&res$deltabeta[,id]>log(thb)
				print(res$levels[res$contrasts[[id]][,2]>0])
				lab=paste(res$levels[res$contrasts[[id]][,2]>0],collapse="-")
			}else{
				flag=res$ltsr[,id]>th&res$deltabeta[,id]<log(1/thb)
				print(res$levels[res$contrasts[[id]][,1]>0])
				lab=paste(res$levels[res$contrasts[[id]][,1]>0],collapse="-")
			}
		}
	}
	print(table(is.na(flag))); 
	if(gid){
		rownames(res$beta)[flag]
	}else{
		query=list(rownames(res$beta)[flag])
		names(query)=lab
		gost(query)
		#gost(query,exclude_iea=T)
		#gost(query,custom_bg=rownames(res$beta),exclude_iea=T)
	}
}
LFLMM <-
function(Y, X, a=rep(0,nrow(Y)), beta=NULL, psi=NULL, delta=NULL, Psi=NULL, theta=NULL, omega=NULL, subset=NULL, nLatent=0, forced=F, ITRMAX=200, PLOT=F, nexpcells=NULL, xgamma, s0=NULL){

	library(Matrix)

	if(nrow(X)!=ncol(Y)){warning("Input TPM matrix Y is not compatible with covariate matrix X."); return()}
	if(!is.data.frame(X)){warning("Covariate matrix X is not a data frame."); return()}
	
	if(is.null(subset)){subset=rep(T, nrow(X))}
	nh = 1
	Z = rep(1,sum(subset))
	isnum = 0
	for(i in seq(ncol(X))){
		if(is.character(X[,i]) || is.factor(X[,i])){
			isnum = c(isnum, 0)
			Z1 = array(0,c(sum(subset),length(table(X[subset,i]))))
			Z1[!is.na(X[subset,i]),] = model.matrix(~0+X[subset,i])
			colnames(Z1) = gsub("X\\[subset, i\\]","",colnames(model.matrix(~0+X[subset,i])))
		}else{
			isnum = c(isnum, 1)
			Z1 = matrix(scale(as.numeric(X[subset,i])),sum(subset))
			Z1[is.na(Z1)]=0
			colnames(Z1) = colnames(X)[i]
		}
		Z  = cbind(Z, Z1)
		nh = c(nh, ncol(Z1))
	}
	names(nh)=c("Intercept", names(X))
	
	isnum = cumsum(1-isnum)*isnum
	Collapse=function(x,flag){
		x1=x[cumsum(flag)==0]
		x2=x[rev(cumsum(rev(flag)))==0]
		res = sum(x[flag])
		names(res) = paste(names(x[flag]),collapse="/")
		return(c(x1, res, x2))
	}
	if(sum(isnum)>0){for(i in seq(max(isnum))){
		if(sum(isnum==i)>2){
			if(!forced){
				mflag = "n" #readline(paste("Do you want to merge ", paste(names(nh)[isnum==i],collapse=","), " factors? [N/y]: ",sep=""))
			}else{
				mflag="n"
			}
			if(sum(mflag%in%c("y","Y"))){
				nh    = Collapse(nh, isnum==i)
				isnum = Collapse(isnum, isnum==i)
			}
		}
	}}
	if(nLatent>0){nh=c(nh, LatentFactor=nLatent);}
print(nh)	
	K = length(nh)
	H = diag(K)[rep(seq(K),nh),]
	X = Z

	# Latent factors
	#Z = cbind(Z, array(0,c(nrow(Z),nLatent)))
	if(!is.null(Psi)){
		Z = cbind(Z, Psi)
	}else{
		Z = cbind(Z, array(0,c(nrow(Z),nLatent)))
	}

	# Likelihood
	logDet = function(x){eval = eigen(x,T)[[1]]; sum(log(eval))}
	Kurt=function(V){diag(V)%*%t(diag(V))+2*V^2}
	Solve = function(X){if(length(X)==1){c(1/X)}else{X=eigen(X,T); r=X[[1]]; r[r<1e-7]=1; X[[2]]%*%diag(1/r)%*%t(X[[2]])}}
	lkhd = function(beta, psi, delta, omega, theta, s0, Y, X, Z, H, Yt, YYt, YX, YZ, y2){
		J = nrow(Y)
		N = ncol(Y)
		M = ncol(Z)
		z = c(apply(Z, 2, sum))
		ZtZ  = t(Z)%*%Z # for W
		
		# A = (ZtZ+U^-2) MAT
		u = c(H%*%psi)
		A = ZtZ+diag(1/u^2)

		#print("sigma2")
		# signa^2 | beta delta Psi(Z) psi(A)
		if(nLatent==0){
			ZtRYt   = t(YZ) # t(as.matrix(Y%*%(Z*R)))
		}else{
			ZtRYt   = t(as.matrix(Y%*%Z))
		}
		ZtRX    = t(Z)%*%X
		ZtRYtmb = (ZtRYt - c(ZtRX%*%beta))
		AinvZtRYtmb = solve(A, ZtRYtmb)
		yR2y = y2 # colSums((Yt*R)^2)
		Xb   = c(X%*%beta)
		YR2Xb = YX%*%beta # as.matrix(Y%*%(Xb*R^2))
		if(is.null(s0)){
			s2 = c( ( yR2y - 2*YR2Xb + sum((Xb)^2) ) - colSums(AinvZtRYtmb*ZtRYtmb) )
			if(!is.null(nexpcells)){
				tau    = (N/2+theta) / (s2/2+theta*c(nexpcells%*%omega))
				logtau = log(tau/(N/2+theta)) + digamma((N/2+theta))
				resglm = GammaGlm(tau[xgamma>0.05],logtau[xgamma>0.05],nexpcells[xgamma>0.05,,drop=F],PLOT=1,x=xgamma[xgamma>0.05]); omega=resglm$omega; theta=resglm$theta
				cat("omega theta=");print(c(omega,theta))
				s2 = (s2/2+theta*c(nexpcells%*%omega)) / (N/2+theta)
			}else{
				s2 = (s2+1)/(N+1)
			}
		}else{s2=s0}

		#print("beta")
		# beta | delta sigma2 Psi(Z) psi(A)
		XtR2X = t(X)%*%(X)
		y = colSums(Y/s2) / sum(1/s2)
		beta0=0.01
		beta = c(Solve(diag(ncol(X))*beta0 + XtR2X-t(ZtRX)%*%solve(A)%*%ZtRX) %*% (t(X)%*%y-t(ZtRX)%*%solve(A,t(Z)%*%y)))
		Xb   = c(X%*%beta)

		#print("Psi")
		# Psi | beta delta psi(A) sigma2
		PL=NULL
		if(nLatent>0){
			cat("Latent factor")
			M2 = M-nLatent
			A2 = t(Z[,1:M2])%*%Z[,1:M2] + diag(1/u[1:M2]^2)
				RXb   = c(X%*%beta)
				YR2Xb = c(YX%*%beta)
				YRZ   = YZ[,1:M2]
				ZtRXb = c(t(Z[,1:M2])%*%RXb)
				CAind = c(YRZ%*%solve(A2, ZtRXb))
				YVYt  = YYt

			YVYt = t(YVYt-YR2Xb)-YR2Xb + sum(RXb^2) - YRZ%*%solve(A2)%*%t(YRZ)
			YVYt = t(YVYt+CAind)+CAind - sum(solve(A2,ZtRXb)*ZtRXb)
			YVYt = t(YVYt/sqrt(s2))/sqrt(s2)/J
			cat("Eigen decomp")
			PL = Eigen(YVYt,nLatent+5)
			
			flage=(apply(PL$evec^2,2,max)/apply(PL$evec^2,2,sum))<0.3
			cat("good pcs=");print(sum(flage[1:nLatent]))
			lambda = PL$eval[1:nLatent]; lambda[lambda<1]=1
			Z[,(M2+1):M] = ( as.matrix((Yt)%*%(PL$evec[,1:nLatent]/sqrt(s2))) - RXb%*%t(colSums(PL$evec[,1:nLatent]/sqrt(s2))) )     %*%diag(sqrt(1-1/lambda))
			ZtZ = t(Z)%*%Z
			A = ZtZ+diag(1/u^2)
			ZtRX  = t(Z)%*%(X) # for W
			ZtRYt = t(as.matrix(Y%*%(Z))) # for W
		}
		
		# Var(xi)
		Vx = ZtZ - ZtZ %*% solve(A, ZtZ)
		hess = diag(psi)%*%t(H)%*%Kurt(Vx)%*%H%*%diag(psi)*2/N + t(H)%*%(Vx * t(solve(A)/u)/u)%*%H*2 - diag(c(t(H)%*%diag(Vx)))*2

		# likelihood
		res = sum(log(s2))*N/2 + J * logDet(diag(M) + t(ZtZ*u)*u)/2
		
		print("psi")
		# psi | beta delta Psi(Z) sigma2
		ZtRYtmb = (ZtRYt - c(ZtRX%*%beta)) 
		AinvZtRYtmb = solve(A, ZtRYtmb)
		attr(res, "gradient") = - (c( t(H)%*%apply(t(t(AinvZtRYtmb^2)/s2)/u^3,1,sum) ) - J * c( t(H)%*%(diag(solve(A,ZtZ))/u) ))
		attr(res, "hessian") = - (hess+t(hess))/2 * J
		attr(res, "beta") = beta
		attr(res, "omega") = omega
		attr(res, "theta") = theta
		attr(res, "psi") = psi
		attr(res, "sigma2") = s2
		attr(res, "Z") = Z
		attr(res, "X") = X
		attr(res, "H") = H
		attr(res, "delta") = delta
		attr(res, "PL") = PL
		res
	}
	
	# Matrix prep
	lkhd.all=NULL
	res.min=NULL
	if(is.null(psi)){psi=rep(0.01,length(nh)); if(nLatent>0){psi[length(psi)]=1}}
	if(is.null(beta)){beta=rep(0.1,ncol(X))}
	if(is.null(delta)){delta=rep(0,nrow(Z))}
	if(is.null(omega) && !is.null(nexpcells)){
		s2 = apply(Y,1,var)
		omega = coef(lm(s2~0+nexpcells))
		lambda = 1/c(nexpcells%*%omega)
		theta = mean(lambda^2/(1/s2-lambda)^2)
		print(c(omega,theta))
	}
	
	Y  = Y[,subset]
	Yt = NULL; YYt=NULL; YX=NULL; YZ=NULL;
	if(nLatent>0){print("matrix prep 1")
		Yt = t(Y)
		YYt=as.matrix(Y%*%Yt)
		YX=as.matrix(Y%*%X)
		YZ=as.matrix(Y%*%Z)
		y2 = rowSums(Y^2)
	}else{print("matrix prep 2")
		print("Y")
		print(dim(Y))
		print(Y[1:5, 1:5])

		print("Z")
                print(dim(Z))
                print(Z[1:5, 1:5])
		
		Y = as.matrix(Y)
		YZ=as.matrix(Y%*%Z)
		YX=as.matrix(Y%*%X)
		y2 = rowSums(Y^2)
	}
	print("matrix prep end")

	# Iteration start
	convstat=F
	for(i in 1:ITRMAX){
		cat(paste("[",i,"]"))
		tmp = lkhd(beta=beta, psi=psi, delta=delta, omega=omega, theta=theta, s0=s0, Y=Y, Z=Z, X=X, H=H, Yt=Yt, YYt=YYt, YZ=YZ, YX=YX, y2=y2)
		#return(tmp)
		if(is.null(lkhd.all) || tmp<min(lkhd.all)){res.min=tmp}
		lkhd.all=c(lkhd.all,as.numeric(tmp))
		if(PLOT)plot(lkhd.all)
		if(nLatent>0){
			K2 = ncol(H)-1
			psi[1:K2] = psi[1:K2] - solve(attr(tmp,"hessian")[1:K2,1:K2]+2*diag(0.00001/psi[1:K2]^2-1),attr(tmp,"gradient")[1:K2]+2*(-0.00001/psi[1:K2]+psi[1:K2]))
		}else{
			psi = psi - solve(attr(tmp,"hessian")+2*diag(0.00001/psi^2-1),attr(tmp,"gradient")+2*(-0.00001/psi+psi))
		}
		#cat("psi=");
		print(psi)
		beta =attr(tmp, "beta")
		delta =attr(tmp, "delta")
		omega =attr(tmp, "omega")
		theta =attr(tmp, "theta")
		#cat("omega=");print(c(omega,theta))	
		Z = attr(tmp, "Z") # Psi
		if(length(lkhd.all)>50){lkhd.all=rev(rev(lkhd.all)[1:50])}
		if(i>1 & abs(diff(rev(lkhd.all)[1:2])/rev(lkhd.all)[1])<1e-7){convstat=T; cat("Converged\n");break}
	}
	c(list(lkhd=-as.numeric(res.min)), attributes(res.min), list(nh=nh), Converged=convstat)
}
