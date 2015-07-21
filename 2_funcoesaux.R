## funcoes "combina", "produto_cartesiano" e "prod_constroi"
## calcula diferenca maxima entre todos os elementos de 2 segmentos consecutivos do vetor de distancias Vdist
dist_max<-function (i,intervalos,V) {
  x<-as.matrix(intervalos[i,1]:intervalos[i,2]) #indices dos elementos que compoem o Vdist[i]
  y<-as.matrix(intervalos[i+1,1]:intervalos[i+1,2]) #indices dos elementos que compoem o Vdist[i+1]
  a<-matrix(apply(x,1,function(x) apply(y,1,function(y) rbind(x,y))),nrow=2) #indices dos elementos 2 a 2 (produto cartesiano)
  max(apply(a,2,function(a) abs(V[a[1]]-V[a[2]]))) #maior diferenca entre todas as diferencas 2 a 2
}
## funcao pico_d
pico_d<-function(nintervalos,V) {
  resto<-length(V)%%nintervalos
  if (resto>0) { V<-V[-(1:resto)] }
  n<-length(V)
  incremento<-n%/%nintervalos
  intervalos<-cbind(seq(1,n,incremento),seq(incremento,n,incremento))
  i<-matrix(1:(nintervalos-1))
  b<-NULL
  b=apply(i,1,function(i) dist_max(i,intervalos,V))
  return(b)
}