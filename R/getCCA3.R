#input: XYZ: biological data sets X,Y; design data set Z
#       numCV number of canonical variables
#
#output:    cc3.weights.{xyz}: matrices with weight vectors for each canonical variable (columnwise)
#           cc3.CV.{xyz}: matrices with canonical variables (columnwise)
#           corr: matrix with correlation coefficient for pairwise (X-Y,X-Z,YZ) and correlation of can. variable: sum(corr(X-Z),corr(Y-Z))/2 (columnwise)
#           lambda: best lambda.{xyz} for each can. var. (columnwise)
getCCA3<-function(X,Y,Z,end=c(0.3,0.5,3),step=c(0.01,0.01,0.2),numCV=10,n.r=10, max.counter.test=10){
	

	if (dim(X)[1]!=dim(Y)[1] || dim(Y)[1]!=dim(Z)[1])
		stop("Data matrix has wrong sample size.")

	cc3.weight.x<-matrix(ncol=0,nrow=dim(X)[2])
	cc3.weight.y<-matrix(ncol=0,nrow=dim(Y)[2])
	cc3.weight.z<-matrix(ncol=0,nrow=dim(Z)[2])
	cc3.CV.x<-matrix(ncol=0,nrow=dim(X)[1])
	cc3.CV.y<-matrix(ncol=0,nrow=dim(X)[1])
	cc3.CV.z<-matrix(ncol=0,nrow=dim(X)[1])

    rownames(cc3.weight.x)<-colnames(X)
    rownames(cc3.weight.y)<-colnames(Y)
    rownames(cc3.weight.z)<-colnames(Z)
    rownames(cc3.CV.x)<-rownames(X)
    rownames(cc3.CV.y)<-rownames(X)
    rownames(cc3.CV.z)<-rownames(X)


	all.lambdas<-matrix(ncol=0,nrow=3)
	all.corr<-c()

    Z<-apply(Z,2,function(y){if (var(y)==0)  y-mean(y) else (y-mean(y))/var(y)})

    dims<-c(dim(X)[2],dim(Y)[2],dim(Z)[2])
	
	#determine pseudo matrix once - Z does not change	
    if (abs(det(var(Z)))>10^-20){ #if var(Z) is invertible
        Zp<-solve(var(Z))%*%t(Z)    
    }else{
	    Zp <- (diag(1/sqrt(diag(var(Z)))))%*%t(Z)#otherwise: regularize
    }
	

    

	canVar=1
    while (canVar <= numCV){#calculate the canonical variables
		    #normalize data matrices
		    X<-apply(X,2,function(y) {if (var(y)==0) y-mean(y) else (y-mean(y))/var(y)})
            Y<-apply(Y,2,function(y) {if (var(y)==0) y-mean(y) else (y-mean(y))/var(y)})

		    #determine Pseudomatrices for every can. corr. step, including ridge regression
		    Xp <- (diag(1/sqrt(diag(var(X)))))%*%t(X)#strong regularization
		    Yp <- (diag(1/sqrt(diag(var(Y)))))%*%t(Y)

            #Pseudomatrices * ...
		    XpZ<-Xp%*%Z
		    YpZ<-Yp%*%Z
            ZpX<-Zp%*%X
            ZpY<-Zp%*%Y


    		results <- get.best.lambdas(X,Y,Z,end=end,n.r=n.r,step=step,max.counter.test=max.counter.test)#return best combination of sparsity parameters

            if (is.null(results)) {break}
            lambdas <- c(results$best.lambda.x,results$best.lambda.y,results$best.lambda.z)
            corr.scca <- results$corr
            xj<-results$bestVector[1:dims[1]]
            yj<-results$bestVector[(dims[1]+1):(dims[1]+dims[2])]
            zj<-results$bestVector[(dims[1]+dims[2]+1):(dims[1]+dims[2]+dims[3])]

		    
		    #update data matrices by removing the latent variable - only for X and Y
            zi<-Z%*%zj
		    xi<-X%*%xj
		    reg<-apply(X,2, function(x) {lm(x~xi)})
		    X<-sapply(reg,function(x){x[[2]]})	
		
		    yi<-Y%*%yj
		    reg<-apply(Y,2, function(x) {lm(x~yi)})
		    Y<-sapply(reg,function(x){x[[2]]})	
		
            
		    all.lambdas<-cbind(all.lambdas,lambdas)
		    all.corr<-c(all.corr,corr.scca)
		    cc3.weight.x<-cbind(cc3.weight.x,xj)
		    cc3.weight.y<-cbind(cc3.weight.y,yj)
		    cc3.weight.z<-cbind(cc3.weight.z,zj)
            cc3.CV.x<-cbind(cc3.CV.x,xi)
            cc3.CV.y<-cbind(cc3.CV.y,yi)
            cc3.CV.z<-cbind(cc3.CV.z,zi)

            canVar<-canVar+1
    }#while canVar

 list(cc3.weight.x=cc3.weight.x,cc3.weight.y=cc3.weight.y,cc3.weight.z=cc3.weight.z, cc3.CV.x=cc3.CV.x,cc3.CV.y=cc3.CV.y,cc3.CV.z=cc3.CV.z, corr=all.corr,lambda=all.lambdas, num.CV=numCV)
}



