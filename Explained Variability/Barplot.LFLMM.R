Barplot <-
function(x, plotfile, title="", col=1){
	
	getUpper95=function(psi, V){
		denom = sum(psi^2)+1
		psivar=NULL
		for(k in 1:length(psi)){
			grad    =         - (2*psi)*psi[k]^2/denom^2
			grad[k] = grad[k] + (2*psi[k])/denom
			psivar = c(psivar, sum(t(V*grad)*grad))
		}
		
		psi^2/denom + 1.96 * sqrt(psivar)
	}

	png(plotfile)

	remids=c(1)
	V = solve(x$hessian)[-remids,-remids]
	labs = names(x$nh)[-remids]
	psi  = x$psi[-remids]
	ord  = rev(order(psi^2))
	par(mar=c(8, 4.1, 4.1, 2.1), mgp=c(2,0.5,0),family="Liberation Sans")
	
	varexp = psi[ord]^2/(1+sum(psi^2))*100
	varse = getUpper95(psi,V)[ord]*100
	labs_plot <- gsub(".", " ", labs, fixed=T)
	labs_plot <- gsub("X ", "%", labs_plot, fixed=T)	

	xcoords = barplot(varexp, name=labs_plot[ord],las=2, ylim=c(0,max(varse)+2), col=col,border=col, main=title)
	mtext("Variance Explained (%)",2,line=2)
		
	segments(xcoords, varexp, xcoords, varse)
	segments(xcoords-0.3, varse, xcoords+0.3, varse)
	
	dev.off()

	res = cbind(varexp, varse)
	rownames(res)=labs_plot[ord]
	res
}
