
#compare patterns of two eigenvectors: return manhattan distance by default; other distances possible
comparePattern<-function(x,y, dm="manhattan"){#, diff
    v<-x-y
    if (dm == "manhattan")
        len<-sum(abs(v))
    if (dm == "euklid")
        len<-sum(v*v)
    return(len)
}

#cluster weight vectors by assigning them to one pattern in list if mean cluster vector is similar to new z.vector;
#make new cluster if not similar to any existing cluster
#count number of vectors associated to one cluster, save correlation for cluster
#cluster items:
#item [[1]]: correlation
#item [[2]]: combined vector of x, y, and z-weight
#item [[3]]: counter for number of vectors in cluster
addPattern<-function(List.Pattern, corr, xyz.vector, z.vector){
        
        if (length(List.Pattern)==0) {
            List.Pattern<-list(list(corr,matrix(xyz.vector,ncol=1,nrow=length(xyz.vector)),1));  
            #print(muster) 
        } else 
        {
            found<-F
            for (i in 1:length(List.Pattern)){
                items<-List.Pattern[[i]]
                s<-dim(items[[2]])[1]-length(z.vector)+1
                e<-dim(items[[2]])[1]
                if (comparePattern(z.vector,rowMeans(as.matrix(items[[2]][s:e,])))<=0.5){#gleiches Muster
                    found<-T
                    List.Pattern[[i]][[3]]<-List.Pattern[[i]][[3]]+1
                    List.Pattern[[i]][[1]]=List.Pattern[[i]][[1]]+corr
                    List.Pattern[[i]][[2]]<-cbind(List.Pattern[[i]][[2]],xyz.vector)
                    break
                }
            }
            if (!found) {
                List.Pattern<-append(List.Pattern,list(list(corr,matrix(xyz.vector,ncol=1,nrow=length(xyz.vector)),1)))
            }
        }
        return(List.Pattern)

}


