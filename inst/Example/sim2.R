
#2 Treatments, 5 Types, 4 Replicates

create.Sim.Dataset<-function() {

    n=12;#number of experiments
    #Design Data Set
    #3 Types (e.g. mutants), each 2 replicates, 2 treatments 
    Type1<-c(rep(1,4),rep(0,2*4))
    Type2<-c(rep(0,4),rep(1,4),rep(0,1*4))
    Rep<-c(rep(c(1,0),3*2))
    Treat<-c(rep(c(rep(1,2),rep(0,2)),3))


    Z<-cbind(Type1,Type2,Rep,Treat)
    
    #Anzahl Gene/Metabolite pro Design Vektor
    genAnz<-c(3,4,  2,  2)
    metAnz<-c(3,3,  3,  5)

    #Erzeuge "sinnvolle" Daten
    X<-matrix(ncol=0,nrow=n)
    Y<-matrix(ncol=0,nrow=n)
    coln=c()
    for (j in 1:length(genAnz)){
        for (k in 1:genAnz[j]){
	        m<-runif(1,3,12)
	        sig<-abs(rnorm(1,m*0.01,m*0.01*5))
            gf<-c()
	        for (i in 1:n)
		        gf<-c(gf,abs(rnorm(1,m,sig))*Z[i,j]+rnorm(1))
            X<-cbind(X,gf)
        }
        if (j<3)
            coln<-c(coln,paste(rep(paste("Type",j,"_",sep=""),genAnz[j]),c(1:genAnz[j]),sep=""))
        if (j>2 & j<4)
            coln<-c(coln,paste(rep(paste("Rep_",sep=""),genAnz[j]),c(1:genAnz[j]),sep=""))
        if (j==4)
            coln<-c(coln,paste(paste("Treat_",c(1:genAnz[j]),sep="")))
        

    }
    colnames(X)<-coln
    coln<-c()
    for (j in 1:length(metAnz)){
        for (k in 1:metAnz[j]){
	        m<-runif(1,0.1,12)
	        sig<-abs(rnorm(1,m*0.01,m*0.01*3))
            gf<-c()
	        for (i in 1:n)
		        gf<-c(gf,abs(rnorm(1,m,sig))*Z[i,j]+rnorm(1))
            Y<-cbind(Y,gf)
        }
        if (j<3)
            coln<-c(coln,paste(rep(paste("Type",j,"_",sep=""),metAnz[j]),c(1:metAnz[j]),sep=""))
        if (j>2 & j<4)
            coln<-c(coln,paste(rep(paste("Rep_",sep=""),metAnz[j]),c(1:metAnz[j]),sep=""))
        if (j==4)
            coln<-c(coln,paste(paste("Treat_",c(1:metAnz[j]),sep="")))
        

    }
    colnames(Y)<-coln

    rownames(X)<-c(1:n)
    rownames(Y)<-c(1:n)
    #Rest: Zufallsdaten
    for (j in 1:5){
        m<-runif(1,0.1,12)
	    sig<-abs(rnorm(1,m*0.01,2))
        random<-c()
	    for (i in 1:n)
		    random<-c(random,abs(rnorm(1,m,sig))+rnorm(1))
        X<-cbind(X,random)
    }
    for (j in 1:5){
        m<-runif(1,0.1,12)
	    sig<-abs(rnorm(1,m*0.01,2))
        random<-c()
	    for (i in 1:n)
		    random<-c(random,abs(rnorm(1,m,sig))+rnorm(1))
        Y<-cbind(Y,random)
    }


    write.table(X,"simX.txt",sep="\t")
    write.table(Y,"simY.txt",sep="\t")
    write.table(Z,"simZ.txt",sep="\t")
}

test.CCA.forSim<-function(){
    # source("get.best.lambdas.R")
    # source("getCCA3.R")
    # source("plotCCA.R")
    # source("save.CCA.R")
    # source("scca.function3.R")
  
    # library(spCCA)  

    X2<-read.table("simX.txt",sep="\t")
    Y2<-read.table("simY.txt",sep="\t")
    Z2<-read.table("simZ.txt",sep="\t")

    CCA3 <- getCCA3(X=X2, Y=Y2, Z=Z2, numCV=4)

}





