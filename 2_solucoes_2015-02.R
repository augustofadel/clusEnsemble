solucoes<-function(arq,p=10,inic=100,k.dens=c(3,4,5,10,15,20,50),tam_pop=100,tot_geracoes=500,out=1) {
  #
  #CARREGAR PACOTES
  require(stats,quietly=T)                          #funcoes: dist, kmeans
  require(fpc,quietly=T)                            #funcoes: kmeansruns, pamk, cluster.stats, dbscan
  #require(mclust,quietly=T)                         #funcoes: 
  #require(methods,quietly=T)                        #funcoes: 
  require(cluster,quietly=T)                        #funcoes: pam, silhouette
  #require(pdfCluster,quietly=T,warn.conflicts=F)    #funcoes: 
  require(birch,quietly=T)                          #funcoes: birch
  #require(stringr,quietly=T)                        #funcoes:   
  require(utils,quietly=T)                          #funcoes: txtProgressBar, setTxtProgressBar
  require(NbClust,quietly=T)                        #funcoes: NbClust
  #
  #SELECIONAR VARIAVEIS NUMERICAS
  #objeto<-arq[,sapply(arq,class)!="numeric"]
  #arq<-arq[,sapply(arq,class)=="numeric"]
  #
  #OBTER MATRIZ DE DISSIMILARIDADES
  dissim<-as.matrix(dist(arq,method="euclidean",upper=T))
  #
  #DEFINIR NUMERO DE GRUPOS
  #cat("\nDefinindo numero de grupos...\n")
  #sol<-kmeansruns(arq,krange=2:ceiling(sqrt(dim(arq)[1])),criterion="asw",runs=inic,scaledata=F,critout=F,plot=F)
  #k<-sol$bestk
  #k.sil<-max(sol$crit)
  #sol<-pamk(dissim,krange=2:ceiling(sqrt(dim(arq)[1])),criterion="ch",usepam=T,scaling=F,diss=T,critout=F)
  #k<-c(k,sol$nc)
  #k.sil<-c(k.sil,max(sol$crit))
  #sol<-NbClust(arq,diss=NULL,distance="euclidean",min.nc=2,max.nc=ceiling(sqrt(dim(arq)[1])),method="ward.D",index="db")
  #k<-c(k,sol$Best.nc[[1]])
  #k.sil<-c(k.sil,sol$Best.nc[[2]])
  #aux<-NULL
  #for (i in k) {aux<-c(aux,i+amp)}
  #k<-unique(aux[aux>1])
  #cat("Concluido.\n")
  #
  #DEFINIR RAIOS DE ABRANGENCIA
  cat("\nDefinindo raios de abrangencia... ")
  distk<-apply(t(apply(dissim,2,function(x) sort(x,partial=(k.dens+1))[(k.dens+1)])),2,function(y) sort(y,method="quick"))
  dimnames(distk)<-list(NULL,as.character(k.dens))
  param<-c(
    apply(distk,2,median),
    apply(distk,2,max),
    apply(apply(distk,2,function(x) {pico_d(10,x)}),2,max),
    apply(apply(distk,2,function(x) {pico_d(20,x)}),2,max)
  )
  cat("Concluido.\n")
  #
  #SOLUCOES
  cat("\nExecutando algoritmos de agrupamento...\n")
  #sol<-sil<-tempo<-NULL
  sol<-tempo<-nome<-NULL
  n<-dim(arq)[1]*.20
  k=2:ceiling(sqrt(dim(arq)[1]))
  for (i in k) {
    cat("\nPara k =",i,"\n")
    #
    #K-MEANS
    cat("K-Means... ")
    for (j in 1:p) {
      t<-proc.time()
      sol.temp<-kmeans(arq,centers=i,nstart=inic,algorithm="Hartigan-Wong")
      tempo<-c(tempo,round((proc.time()-t)[[3]],5))    
      sol<-cbind(sol,sol.temp$cluster)
      #sil<-c(sil,cluster.stats(dissim,sol.temp$cluster,silhouette=T)$avg.silwidth)
      nome<-c(nome,paste("k-means k=",i,"_",j,sep=""))
    }
    cat("Concluido.\n")
    #
    #PAM
    #cat("PAM... ")
    #t<-proc.time()
    #sol.temp<-pam(dissim,k=i,diss=T,cluster.only=F,keep.diss=F,keep.data=F)
    #tempo<-c(tempo,round((proc.time()-t)[[3]],5))
    #sol<-cbind(sol,sol.temp$clustering)
    #sil<-c(sil,sol.temp$silinfo[[3]])
    #cat("Concluido.\n")
    #
    #CLARA
    cat("CLARA... ")
    for (j in 1:p) {
      t<-proc.time()
      sol.temp<-clara(arq,k=i,samples=10,sampsize=n,metric="euclidean",stand=F,pamLike=T,medoids.x=F,keep.data=F)
      tempo<-c(tempo,round((proc.time()-t)[[3]],5))
      sol<-cbind(sol,sol.temp$clustering)
      #sil<-c(sil,sol.temp$silinfo[[3]])
      nome<-c(nome,paste("clara k=",i,"_",j,sep=""))
    }
    cat("Concluido.\n")
    #
    #BRKGA MEDOIDS
    cat("BRKGA Medoids...\n")
    #for (j in 1:p) {
      sol.temp<-BKGA_medoids_af(diretorio=NULL,arquivo=arq,clusters=i,tam_pop,tot_geracoes,pelite=0.20,pmutant=0.15,pr=0.7,distk=dissim)
      tempo<-c(tempo,sol.temp$tempo)
      sol<-cbind(sol,sol.temp$clustering[,2])
      #sil<-c(sil,cluster.stats(dissim,sol.temp$clustering[,2],silhouette=T)$avg.silwidth)
      nome<-c(nome,paste("brkga k=",i,sep=""))
    #}
    cat("\nConcluido.\n")
  }
  cat("\nExecutando algoritmos MRDBSCAN e BIRCH...\n")
  pb<-txtProgressBar(min=0,max=length(param),style=3)
  for (i in 1:length(param)) {
    #
    #MRDBSCAN
    t<-proc.time()
    sol.temp<-dbscan(arq,eps=param[i],MinPts=as.numeric(names(param)[1]),scale=F,method="hybrid",showplot=0)
    tempo<-c(tempo,round((proc.time()-t)[[3]],5))
    sol<-cbind(sol,sol.temp$cluster)
    #if (length(unique(sol.temp$cluster))>1) {
    #  if (out==0) {
    #    sil<-c(sil,summary(silhouette(sol.temp$cluster[sol.temp$cluster!=0],dissim[sol.temp$cluster!=0,sol.temp$cluster!=0]))$avg.width)
    #  } else {
    #    sil<-c(sil,summary(silhouette(sol.temp$cluster,dissim))$avg.width)
    #  }
    #} else { sil<-c(sil,NA) }
    setTxtProgressBar(pb,i-.5)
    #
    #BIRCH
    t<-proc.time()
    sol.temp<-birch(as.matrix(arq),radius=param[i],compact=as.numeric(names(param)[1]))
    tempo<-c(tempo,round((proc.time()-t)[[3]],5))
    aux<-cbind(unlist(sol.temp$members),rep(1:length(sol.temp$members),sol.temp[[1]]))
    aux<-aux[order(aux[,1]),2]
    sol<-cbind(sol,aux)
    #sil<-c(sil,cluster.stats(dissim,aux,silhouette=T)$avg.silwidth)
    setTxtProgressBar(pb,i)
  }
  cat("\nConcluido.\n\n")
  #
  #RESULTADOS
  #sol<-rbind(sol,sil,tempo)
  sol<-rbind(sol,tempo)
  dimnames(sol)[[2]]<-c(nome,paste(rep(c("dbscan r=","birch r="),length(param)),rep(round(param,2),each=2)," d=",rep(names(param),each=2),sep=""))
  #
  #SAIDA
  #assign(x="dissim",value=dissim,pos=.GlobalEnv)
  return(list(sol=sol,dissim=dissim))
}

#sol<-lapply(dados,comb.completo)
