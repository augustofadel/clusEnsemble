carregar<-function(path,std) {
  aux<-getwd()
  setwd(path)
  arquivos<-list.files(path)
  if (length(std)!=length(arquivos)) {cat("\nERRO: o atributo std deve ter um elemento para cada arquivo a ser carregado.\n\n")}
  cat("\nCarregando arquivos...\n\n")
  dados<-lapply(arquivos,function(x) {read.table(x,header=F,skip=0,sep=";",dec=",",nrows=2000,comment.char="",fileEncoding="UTF-16")}) #nrows=2000 para economizaer memoria. Alterar em caso de arquivos maiores ou excluir em caso de arquivos de tamanho desconhecido | comment.char="" p/ economizar memoria.
  names(dados)<-arquivos
  for (i in arquivos) {
    cat(i,": identificadas",dim(dados[[i]])[2],"variaveis e",dim(dados[[i]])[1],"observacoes.\n")
    if(std[which(arquivos==i)]) {dados[[i]]<-scale(dados[[i]],scale=T,center=T)} #padroniza variaveis numericas
  }
  setwd(aux)
  cat("\nTotal de arquivos carregados: ",length(dados),".\n",sep="")
  return(dados)
}