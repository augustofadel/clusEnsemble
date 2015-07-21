indices<-function(dat,arq,dissim,indic=c('within.cluster.ss','avg.silwidth','pearsongamma','dunn','dunn2','entropy','wb.ratio','ch','widestgap','sindex' ,  "C_index","Calinski_Harabasz","Davies_Bouldin","Gamma","PBM","SD_Scat","SD_Dis","S_Dbw","dbs","GDI11","GDI12","GDI13","GDI21","GDI22","GDI23","GDI31","GDI32","GDI33","GDI41","GDI42","GDI43","GDI51","GDI52","GDI53")) {
  require(fpc,quietly=T)          #funcoes: cluster.stats
  #require(cluster,quietly=T)      #funcoes: silhouette
  #require(stringr,quietly=T)      #funcoes: str_detect
  require(clusterCrit,quietly=T)  #funcoes: intCriteria, getCriteriaNames
  require(pdfCluster,quietly=T)   #funcoes: dbs
  
  ###PACOTE clusterCrit
  #apply(sol[[1]]$sol[-dim(sol[[1]]$sol)[1],],2,function(x) {
  #  intCriteria(traj=as.matrix(dados[[1]]),
  #              part=as.integer(x),
  #              crit=c("Silhouette"))
  #})
  ###

  indicTF<-indic %in% getCriteriaNames(isInternal=T)
  if ("dbs" %in% indic) {
    indval<-matrix(NA,nrow=dim(arq)[2],ncol=length(indic)+1,dimnames=list(dimnames(arq)[[2]],c(indic[indic!="dbs"],"dbs.mean","dbs.median")))
    indicTF<-c(indicTF,F)
  } else {
    indval<-matrix(NA,nrow=dim(arq)[2],ncol=length(indic),dimnames=list(dimnames(arq)[[2]],indic))
  }
  pb<-txtProgressBar(min=0,max=dim(arq)[[2]],style=3)
  for (j in 1:dim(arq)[[2]]) {
    #cat(j,"...",sep="")
    if (length(unique(arq[,j]))>1) {
      arq[arq[,j]==0,j]<-as.integer(max(arq[,j])+1)   #substitui zeros - funcao intCriteria nao suporta objetos com rotulo=0
      if ("dbs" %in% indic) {
        aux<-dbs(x=dat,cluster=arq[,j])@dbs
        indval[j,"dbs.mean"]<-mean(aux)
        indval[j,"dbs.median"]<-median(aux)
      }
      if (any(!indicTF)) {
        aux<-cluster.stats(dissim,arq[,j],silhouette=T)
        indval[j,!indicTF & is.na(indval[j,])]<-unlist(aux)[indic[!indicTF & is.na(indval[j,])]]
      }
      if (any(indicTF)) {
        aux<-intCriteria(traj=dat,part=arq[,j],crit=indic[indicTF])
        indval[j,indicTF]<-unlist(aux)[tolower(indic[indicTF])]
      }
    }
    #cat(" Concluido.\n")
    setTxtProgressBar(pb,j)
  }
  
  #for (j in dimnames(arq)[[2]]) {
  #  cat(j,"...",sep="")
  #  if (length(unique(arq[,j]))>1) {
  #    aux<-cluster.stats(dissim,arq[,j],silhouette=T)
  #    for (i in indic) {
  #      indval[j,i]<-aux[[i]]
  #    }
  #  }
  #  cat(" Concluido.\n")
  #}
  
  return(indval)
}