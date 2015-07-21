combinacao<-function(arq,indval,maximo) {
  #
  #CARREGAR PACOTES E FUNCOES
  require(clue,quietly=T)
  #
  #CARREGAR SOLUCOES (elimina solucoes sem indice de validacao)
  #ind<-indval[!is.na(indval[,criterio]),criterio]   #vetor criterios
  ind<-indval[apply(indval,1,function(x) {all(!is.na(x))}),]   #criterios
  #rownames(arq)<-make.unique(rownames(arq))
  sol<-as.data.frame(arq[!(rownames(arq) %in% c('tempo')),apply(indval,1,function(x) {all(!is.na(x))})])
  #sol<-arq[!(dimnames(arq)[[1]] %in% c('tempo')),apply(indval,1,function(x) {all(!is.na(x))})]
  #tempo<-arq[(dimnames(arq)[[1]] %in% c('tempo')),apply(indval,1,function(x) {all(!is.na(x))})]
  #
  #DEFINICAO DOS CONJUNTOS DE PARTICOES BASE
  #if (criterio=='avg.silwidth') { ifelse (sum(ind>0.7)>1,ind.lim<-c(.25,.5,.7),ind.lim<-c(min(summary(ind[ind>=.25])[[1]],.25),summary(ind[ind>=.25])[[3]],summary(ind[ind>=.25])[[5]])) }   #descarta particoes fora do criterio de validacao (criterio combinado: RETZER E SHAN ou quartis sils>.25)
  #obj<-matrix(NA,dim(ind)[1],dim(ind)[2],dimnames=list(NULL,dimnames(ind)[[2]]))
  #ifelse(sum(maximo)>1,obj[,maximo]<-apply(ind[,maximo],2,function(x) { dimnames(ind)[[1]][order(-x,tempo)]}),obj[,maximo]<-dimnames(ind)[[1]][order(-ind[,maximo],tempo)])
  #ifelse(sum(!maximo)>1,obj[,!maximo]<-apply(ind[,!maximo],2,function(x) { dimnames(ind)[[1]][order(x,tempo)]}),obj[,!maximo]<-dimnames(ind)[[1]][order(ind[,!maximo],tempo)])
  #rank<-rep(NA,dim(ind)[1])
  #names(rank)<-dimnames(ind)[[1]]
  #for (i in names(rank)) {rank[i]<-sum(apply(obj,2,function(x) {which(x==i)}))}
  #q<-summary(rank)[c(2,3,5)]
  #conj<-sapply(q,function(x){rank<=x})
  algs<-c('k-means','clara','brkga','dbscan','birch')
  brp<-matrix(NA,nrow=length(algs),ncol=dim(ind)[2],dimnames=list(algs,dimnames(ind)[[2]]))
  for (alg in algs) {
    aux<-ind[str_detect(dimnames(ind)[[1]],alg),]
    ifelse(sum(maximo)>1,brp[alg,maximo]<-dimnames(aux)[[1]][apply(aux[,maximo],2,which.max)],brp[alg,maximo]<-dimnames(aux)[[1]][which.max(aux[,maximo])])
    ifelse(sum(!maximo)>1,brp[alg,!maximo]<-dimnames(aux)[[1]][apply(aux[,!maximo],2,which.min)],brp[alg,!maximo]<-dimnames(aux)[[1]][which.min(aux[,!maximo])])
  }
  #sem solucoes repetidas
  conj<-list(brp=unique(as.vector(brp)),brp.B=NULL,brp.C=NULL)
  for (alg in algs[1:3]) {conj$brp.B<-c(conj$brp.B,conj$brp[str_detect(conj$brp,alg)])}
  for (alg in algs[4:5]) {conj$brp.C<-c(conj$brp.C,conj$brp[str_detect(conj$brp,alg)])}
  #com solucoes repetidas
  #conj<-list(brp=as.vector(brp),brp.B=NULL,brp.C=NULL)
  #for (alg in algs[1:3]) {conj$brp.B<-c(conj$brp.B,conj$brp[str_detect(conj$brp,alg)])}
  #for (alg in algs[4:5]) {conj$brp.C<-c(conj$brp.C,conj$brp[str_detect(conj$brp,alg)])}
  #
  #DECLARACAO DAS MATRIZES DE RESULTADOS
  sol.qtd<-sapply(conj,length)
  consenso<-matrix(NA,nrow=dim(sol)[1],ncol=length(conj),dimnames=list(NULL,names(conj)))   #matriz de particoes consenso  
  result<-matrix(NA,max(sol.qtd)+1,4*(length(conj)),dimnames=list(NULL,c(paste(c("metodo","NMI consenso","ANMI","SNMI"),rep(paste("(",names(conj),")",sep=""),each=4),sep=" "))))   #matriz de resultados (avaliar utilizar 'metodo' como rotulos das linhas ou utilizar data frame)
  #result[1:dim(sol)[[2]],1:2]<-cbind(names(rank),rank)
  #
  #EXECUTAR COMBINACOES DE AGRUPAMENTO
  cat("Executando combinacao de agrupamentos...\n\n")
  tempo<-NULL
  for (j in 1:length(conj)) {
    sol.v<-as.list(sol[,conj[[j]]])
    sol.v<-lapply(sol.v,as.cl_hard_partition)   #converte "vetores-particao" em objeto particao do pacote clue
    comb<-cl_ensemble(list=sol.v)    #particoes iniciais consideradas na combinacao
    #cat("\nConjunto",j,"\nValor maximo do criterio combinado:",q[j],"\nNumero total de solucoes:",sol.qtd[j],sep=" ")
    cat("Conjunto",names(conj[j]),sep=" ","\n")
    t<-proc.time()
    cons<-cl_consensus(comb,method="GV3",control=list(nruns=100,method="CG"))  #particao concenso (metodo: maximiza co-associacao) - ERRO GERENCIAMENTO MEMORIA
    tempo<-c(tempo,round((proc.time()-t)[[3]],5))
    cons<-as.cl_hard_partition(cons)
    consenso[,j]<-cons[[1]]
    sol.v<-as.list(cbind(sol[,conj[[j]]],as.vector(cons[[1]])))    #matriz de informacoes mutuas normalizadas
    names(sol.v)[length(sol.v)]<-"consenso"
    sol.v<-lapply(sol.v,as.cl_hard_partition)
    NMI<-matrix(NA,nrow=length(sol.v),ncol=length(sol.v)+1,dimnames=list(names(sol.v),c(names(sol.v)[-length(sol.v)],"ANMI","SNMI")))
    for (l in 1:(dim(NMI)[1]-1)) {
      NMI[,l]<-cl_agreement(x=sol.v,y=sol.v[[l]],method="NMI")
    }
    for (l in 1:dim(NMI)[1]) { NMI[l,(dim(NMI)[2]-1)]<-mean(NMI[l,-l],na.rm=T) ; NMI[l,(dim(NMI)[2])]<-sum(NMI[l,-l],na.rm=T) }
    result[1:(sol.qtd[j]+1),(j+(j-1)*3)]<-c(conj[[j]],"consenso")
    result[1:sol.qtd[j],(j+(j-1)*3+1)]<-head(NMI[dim(NMI)[1],],dim(NMI)[2]-2)
    result[1:(sol.qtd[j]+1),(j+(j-1)*3+2)]<-NMI[,(dim(NMI)[2]-1)]
    result[1:(sol.qtd[j]+1),(j+(j-1)*3+3)]<-NMI[,dim(NMI)[2]]
  }
  #
  #SAIDA
  return(list(sol.qtd=sol.qtd,consenso=consenso,result=result,tempo=tempo))
}