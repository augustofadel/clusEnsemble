#### 0 - INICIALIZACAO ####
#diretorios
path.fun="~/Desktop/SOBRAPO"                    #localizacao das funcoes
path.dat="~/Desktop/SOBRAPO/dados2"             #localizacao dos arquivos de dados
path.res="~/Desktop/SOBRAPO/resultados_SE"        #localizacao dos arquivos de saida
#pacotes
library(stringr)
#funcoes
source(file.path(path.fun,"1_carregar_2015-02.r",fsep=.Platform$file.sep))
source(file.path(path.fun,"2_funcoesaux.r",fsep=.Platform$file.sep))
source(file.path(path.fun,"2_BKGA_medoids_af_2015-01.r",fsep=.Platform$file.sep))
source(file.path(path.fun,"2_solucoes_2015-02.r",fsep=.Platform$file.sep))
source(file.path(path.fun,"3_indices_2015-02.r",fsep=.Platform$file.sep))
source(file.path(path.fun,"4_combinacao_2015-03.r",fsep=.Platform$file.sep))


#### 1 - CARREGAR INSTANCIAS ####
#Executa funcao 'carregar': carrega todos os arquivos com extensao '.csv' existentes em 'path. Cria lista 'dados', na qual cada objeto corresponde a uma instancia.
dados<-carregar(path=path.dat,                           #localizacao dos arquivos de dados
                std=c(T))#std=c(rep(F,4),rep(T,11),rep(F,10)))     #Padronizar variaveis (atributos)?
                


#### 2 - OBTER SOLUCOES DE AGRUPAMENTO ####
#EXECUCAO ALTERNATIVA (salva os resultados a cada iteracao - evita problemas gerenciamento memoria):
for (i in names(dados)) {
  cat("\n\nARQUIVO:",i,"\n")
  sol<-solucoes(
    arq=dados[[i]],
    p=10,                                   #k-means e CLARA: numero de execucoes para cada valor de k
    inic=100,                               #k-means: inicializacoes aleatorias de centroides
    k.dens=c(3,4,5,10,15,20,50),            #DBSCAN e BIRCH: criterios de densidade
    tam_pop=100,                            #BRKGA-M: tamanho da populacao inicial de vetores de chaves aleatorias
    tot_geracoes=500,                       #BRKGA-M: total de geracoes
    out=1                                   #parametro utilizado no calculo da silhueta media das solucoes obtidas a partir do dbscan: 0=desconsiderar objetos rotulados como outliers
  )
  save(sol,file=file.path(path.res,paste(str_sub(i,1,nchar(i)-2),"_sol.RData",sep=""),fsep=.Platform$file.sep))
}
#Recompor lista de solucoes e matrizes de similaridade:
arquivos<-list.files(path=path.res,pattern="_sol.RData")
sol.aux<-sapply(str_replace(arquivos,"_sol.RData",".i"),function(x) NULL)
for (i in arquivos) {
  load(file.path(path.res,i,fsep=.Platform$file.sep))
  sol.aux[[str_replace(i,"_sol.RData",".i")]]<-sol
}
sol<-sol.aux
rm(sol.aux)
nome="solucoes.RData"
save(sol,file=file.path(path.res,nome,fsep=.Platform$file.sep))


#### 3 - CALCULAR INDICES DE VALIDACAO ####
nome="solucoes.RData"                                             #nome do arquivo contendo as solucoes
load(file.path(path.res,nome,fsep=.Platform$file.sep))
#dados<-carregar(path=path.dat,std=c(rep(F,5),rep(T,11),rep(F,10)))                  #carregar dados


for (i in names(sol)) {
  cat("\n\nARQUIVO:",i,"\n\nCalculando indices de validacao:\n\n")
  indval<-indices(dat=as.matrix(dados[[i]]),
                  arq=apply(sol[[i]]$sol[!(dimnames(sol[[i]]$sol)[[1]] %in% c('tempo')),],2,as.integer),
                  dissim=sol[[i]]$dissim,
                  #indic=c("avg.silwidth","pearsongamma","C_index","Calinski_Harabasz","Davies_Bouldin","GDI22","Gamma","PBM","SD_Scat","SD_Dis","dbs"))
                  indic=c("avg.silwidth","S_Dbw","dbs"))
  save(indval,file=file.path(path.res,paste(str_sub(i,1,nchar(i)-2),"_indice.RData",sep=""),fsep=.Platform$file.sep))
  #rank<-as.data.frame(cbind(sil=ind,tempo=sol[[i]]$sol[(dimnames(sol[[i]]$sol)[[1]] %in% c('tempo')),]))
  #write.table(rank[order(rank$sil,-rank$tempo,decreasing=T),],file=file.path(path.res,str_replace(i,".csv","_rank.csv"),fsep=.Platform$file.sep),na="",sep=";",dec=",",quote=T,row.names=T,col.names=T,fileEncoding="MS-ANSI")
}
#Recompor lista de indices de validacao:
arquivos<-list.files(path=path.res,pattern="_indice.RData")
indval.aux<-sapply(str_replace(arquivos,"_indice.RData",".i"),function(x) NULL)
for (i in arquivos) {
  load(file.path(path.res,i,fsep=.Platform$file.sep))
  indval.aux[[str_replace(i,"_indice.RData",".i")]]<-indval
}
indval<-indval.aux
rm(indval.aux)
nome="indices.RData"
save(indval,file=file.path(path.res,nome,fsep=.Platform$file.sep))


#### 4 - OBTER COMBINACOES DE AGRUPAMENTO ####
nome="solucoes.RData"                                             #nome do arquivo contendo as solucoes
load(file.path(path.res,nome,fsep=.Platform$file.sep))
nome="indices.RData"
load(file.path(path.res,nome,fsep=.Platform$file.sep))

#Executar funcao 'combinacao' utilizando estrutura de repeticao:
#EXECUCAO ALTERNATIVA (salva os resultados a cada iteracao - problemas gerenciamento memoria):
for (i in names(sol)) {
  cat("\n\nARQUIVO:",i,"\n\n")
  comb<-combinacao(arq=sol[[i]]$sol,
                   indval=indval[[i]],
                   maximo=c(T,F,T,T))
  
  #Exportar resultados:
  write.table(comb$result,file=file.path(path.res,paste(str_sub(i,1,nchar(i)-2),"_NMI.csv"),fsep=.Platform$file.sep),na="",sep=";",dec=",",quote=F,row.names=F,col.names=T,fileEncoding="MS-ANSI")
  cat("Gerado arquivo de NMIs.\n")
  write.table(comb$consenso,file=file.path(path.res,paste(str_sub(i,1,nchar(i)-2),"_consenso.csv"),fsep=.Platform$file.sep),sep=";",dec=",",quote=F,row.names=F,col.names=T,fileEncoding="MS-ANSI")
  cat("Gerado arquivo de particoes consenso.\n")
  
  #Carregar resultados
  #comb<-list(
  #  result=read.table(file.path(path.res,paste(str_sub(i,1,nchar(i)-2),"_NMI.csv"),fsep=.Platform$file.sep),sep=";",dec=",",fileEncoding="MS-ANSI",header=F,skip=1,stringsAsFactors=F),
  #  consenso=read.table(file.path(path.res,paste(str_sub(i,1,nchar(i)-2),"_consenso.csv"),fsep=.Platform$file.sep),sep=";",dec=",",fileEncoding="MS-ANSI",header=F,skip=1,stringsAsFactors=F),
  #  sol.qtd=1:3)
  
  #Criar e exportar diagrama comparativo criterios de composicao dos conjuntos de particoes base:
  m<-seq(3,dim(comb$result)[2],4)
  #lim<-range(as.numeric(comb$result[,m]),na.rm=T)
  lim<-range(as.numeric(unlist(comb$result[,m])),na.rm=T)
  #grafico ANMIs particoes consenso
  ANMI.cons<-NULL
  for (j in m) {
    ANMI<-as.numeric(comb$result[,j])[!is.na(comb$result[,j])]
    ANMI.cons<-c(ANMI.cons,ANMI[length(ANMI)])
  }
  jpeg(file=file.path(path.res,paste(str_sub(i,1,nchar(i)-2),"_NMI.jpeg"),fsep=.Platform$file.sep),width=7.9,height=6,units="cm",res=300,quality=100,bg="transparent")
  par(cex=.5,mar=c(4,4,4,1))
  plot(1:length(comb$sol.qtd),ANMI.cons,ylim=lim*c(1,1.01),xlim=c(0,length(m)+1),xaxt='n',pch=19,lwd=1,col=4,type="b",main=paste("Diagrama comparativo\nInstância:",str_sub(i,1,nchar(i)-2),sep=" "),xlab="Conjuntos de particoes base",ylab="ANMI")
  text(1:length(comb$sol.qtd),y=ANMI.cons,labels=comb$sol.qtd,pos=3,offset=.3,cex=.5)
  text(1:length(comb$sol.qtd),par("usr")[3],labels=c('A','B','C'),adj=c(1.1,1.1),xpd=TRUE,cex=.9) #labels=names(comb$sol.qtd),srt=45
  #boxplots ANMIs particoes base X conjunto de particoes base
  for (j in 1:length(m)) {
    ANMI<-as.numeric(comb$result[,m[j]])[!is.na(comb$result[,m[j]])]
    boxplot(ANMI[-length(ANMI)],add=T,at=j)
  }
  dev.off()
  cat("Gerado diagrama comparativo.\n")
  
  cat("\nConcluido.\n")
}