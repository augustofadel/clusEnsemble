BKGA_medoids_af<-function(diretorio,arquivo,clusters,tam_pop=100,tot_geracoes=500,pelite=0.20,pmutant=0.15,pr=0.7,distk=dissim)
{
#clusters = n?mero de clusters
#tam_pop = n?mero de solu??es do algoritmo gen?tico
#######################################################################################
leitura_dados<-function(arquivo)
{dados=as.matrix(arquivo)
 #dados=as.matrix(read.table(paste(diretorio,"~/Desktop/SOBRAPO/dados/",arquivo,sep=""),sep="\t",dec="."))
 xmed=apply(dados,2,mean)
 sdmed=apply(dados,2,sd)
 numero_objetos=dim(dados)[1]
 numero_atributos=dim(dados)[2]
 for(i in 1:numero_atributos) {dados[,i]=(dados[,i]-xmed[i])/sdmed[i]}
 #distk=as.matrix(dist(dados))
 return(list(dados=dados,distk=distk,numero_objetos=numero_objetos))
}
#####################################################################################
path_relinking<-function(gelite,felite,clusters,numero_objetos,distk,te)
{
fobj<-unique(felite)
qmelhores<-min(10,length(fobj))
idf<-sort(match(fobj,felite))
#Escolhe duas solu??es do conjunto elite
escolhidas<-1:qmelhores
combinacoes<-combn(qmelhores,2)

solucao<<-NULL
combina_path<-function(combinacoes)
{
 v1<-combinacoes[1]
 v2<-combinacoes[2]
 A<-rbind(gelite[idf[escolhidas[v1]],],gelite[idf[escolhidas[v2]],])
 candidatos<-setdiff(1:numero_objetos,as.vector(A))
 m1<-NULL
 m2<-NULL
 for(i in 1:clusters)
  {medoid_s1<-A[1,i]
   medoid_s2<-A[2,i]
   w1<-which.min(distk[candidatos,medoid_s1])
   m1<-c(m1,candidatos[w1])
   candidatos<-candidatos[-w1]
   w2<-which.min(distk[candidatos,medoid_s2])
   m2<-c(m2,candidatos[w2])
   candidatos<-candidatos[-w2]
  }
 m3<-NULL
 m4<-NULL
 for(i in 1:clusters)
  {medoid_s1<-A[1,i]
   medoid_s2<-A[2,i]
   w3<-which.min(distk[candidatos,medoid_s1])
   m3<-c(m3,candidatos[w3])
   candidatos<-candidatos[-w3]
   w4<-which.min(distk[candidatos,medoid_s2])
   m4<-c(m4,candidatos[w4])
  }
 m5<-NULL
 m6<-NULL
 for(i in 1:clusters)
  {medoid_s1<-A[1,i]
   medoid_s2<-A[2,i]
   w5<-which.min(distk[candidatos,medoid_s1])
   m5<-c(m5,candidatos[w5])
   candidatos<-candidatos[-w5]
   w6<-which.min(distk[candidatos,medoid_s2])
   m6<-c(m6,candidatos[w6])
  }
 m7<-NULL
 m8<-NULL
 for(i in 1:clusters)
  {medoid_s1<-A[1,i]
   medoid_s2<-A[2,i]
   w7<-which.min(distk[candidatos,medoid_s1])
   m7<-c(m7,candidatos[w7])
   candidatos<-candidatos[-w7]
   w8<-which.min(distk[candidatos,medoid_s2])
   m8<-c(m8,candidatos[w8])
  }
 A<-rbind(A[1,],m1,m2,m3,m4,m5,m6,m7,m8,A[2,])
 row.names(A)<-NULL
 for(k in 1:clusters)
  {for(z in 2:10){solucao<<-rbind(solucao,c(A[z,k],A[1,-k]))}}
  
 v1<-combinacoes[2]
 v2<-combinacoes[1]
 A<-rbind(gelite[idf[escolhidas[v1]],],gelite[idf[escolhidas[v2]],])
 candidatos<-setdiff(1:numero_objetos,as.vector(A))
 m1<-NULL
 m2<-NULL
 for(i in 1:clusters)
  {medoid_s1<-A[1,i]
   medoid_s2<-A[2,i]
   w1<-which.min(distk[candidatos,medoid_s1])
   m1<-c(m1,candidatos[w1])
   candidatos<-candidatos[-w1]
   w2<-which.min(distk[candidatos,medoid_s2])
   m2<-c(m2,candidatos[w2])
   candidatos<-candidatos[-w2]
  }
 m3<-NULL
 m4<-NULL
 for(i in 1:clusters)
  {medoid_s1<-A[1,i]
   medoid_s2<-A[2,i]
   w3<-which.min(distk[candidatos,medoid_s1])
   m3<-c(m3,candidatos[w3])
   candidatos<-candidatos[-w3]
   w4<-which.min(distk[candidatos,medoid_s2])
   m4<-c(m4,candidatos[w4])
  }
 m5<-NULL
 m6<-NULL
 for(i in 1:clusters)
  {medoid_s1<-A[1,i]
   medoid_s2<-A[2,i]
   w5<-which.min(distk[candidatos,medoid_s1])
   m5<-c(m5,candidatos[w5])
   candidatos<-candidatos[-w5]
   w6<-which.min(distk[candidatos,medoid_s2])
   m6<-c(m6,candidatos[w6])
  }
 m7<-NULL
 m8<-NULL
 for(i in 1:clusters)
  {medoid_s1<-A[1,i]
   medoid_s2<-A[2,i]
   w7<-which.min(distk[candidatos,medoid_s1])
   m7<-c(m7,candidatos[w7])
   candidatos<-candidatos[-w7]
   w8<-which.min(distk[candidatos,medoid_s2])
   m8<-c(m8,candidatos[w8])
  }
 A<-rbind(A[1,],m1,m2,m3,m4,m5,m6,m7,m8,A[2,])
 row.names(A)<-NULL
 for(k in 1:clusters)
  {for(z in 2:10){solucao<<-rbind(solucao,c(A[z,k],A[1,-k]))}}
  
  

 return(0)
}
soluca<-apply(combinacoes,2,function(combinacoes) combina_path(combinacoes))
vq<-apply(solucao,1,function(x) length(table(x)))
vq<-which(vq==clusters)
solucao<-solucao[vq,]
uz<-apply(solucao,1,function(x) x%*%c(10^(1:clusters)))
uw<-sort(unique(uz))
uq<-match(uw,uz)
solucao<-solucao[uq,]
si<-t(apply(solucao,1,function(x) fitness(x,distk,clusters,numero_objetos)))
pv<-which.min(si[,1])
return(list(a=solucao[pv,],b=si[pv,]))
}


#####################################################################################
#####################################################################################
############################ Opera??es nos Medoids###################################
#Constru??o da Popula??o
pop_gen<-function(psize,N)
{popa<-matrix(runif(psize*N,0,1),nrow=psize,ncol=N,byrow=TRUE)
return(popa)
}

####################################################################################
minimare<-function(vet,k)
{kvet<-NULL
for(i in 1:k) {p=which.min(vet);kvet<-c(kvet,p);vet[p]<-Inf}
return(kvet)
}
#########################################################################################
decoder<-function(pop,clusters,linesp)
{decod<-matrix(apply(pop,1,function(x) minimare(x,clusters)),dim(pop)[1],clusters,byrow=TRUE)
return(decod)
}



#######################################################################################
fitness<-function(g,distk,clusters,numero_objetos)
{
fk<-NULL
nk<-NULL
objx<-as.matrix(c(1:numero_objetos))
fobj<-sum(apply(objx,1,function(x) min(distk[g,x])))
return(c(fobj/numero_objetos,g))
}


########################################################################################
########################################################################################

crossover_uniforme<-function(gelite,gnonelite,pe,pelite,pnelite,psort,N,sq)
{
uniform<-function(ge,gn,ab,pe,N,sq)
{pw<-runif(N,0,1)
 gp<-which(pw<pe)
 gk<-gn[ab[2],]
 gk[gp]<-ge[ab[1],gp]
 return(gk)
}
ab=cbind(sample(pelite,psort,replace=TRUE),sample(pnelite,psort,replace=TRUE))
gnew=t(apply(ab,1,function(ab) uniform(gelite,gnonelite,ab,pe,N,sq)))
return(gnew)
}
######################################################################################
############### Programa Principal ######################
#########################################################
dados<-leitura_dados(arquivo)
distk=dados$distk
numero_objetos=dados$numero_objetos
#########################################################
tempo=proc.time()
fbest_mult<-Inf
#tam_pop<-min(max(round(numero_objetos*0.25),50),150)
size_population<-tam_pop
pelite<-round(pelite*size_population)
pmutant<-round(pmutant*size_population)
f<-pop_gen(size_population*10,numero_objetos)
f<-f[(order(f[,1]))[1:size_population],]
#Decodifica??o para a solu??o do problema de otimiza??o combinat?ria
g<-decoder(f,clusters,size_population)
objetos=c(1:numero_objetos)
#C?lculo da Fun??o Objetivo
ft<-t(apply(g,1,function(x) fitness(x,distk,clusters,numero_objetos)))
fbest<-Inf
i<-0
sq<-1:numero_objetos
#library(snowfall)
#sfInit(parallel=TRUE,cpus=2)
pb<-txtProgressBar(min=0,max=tot_geracoes,style=3)
while(i<tot_geracoes)
    {i<-i+1
     #cat("Gera??o ",i,"\n")
     flush.console()
     pq<-order(ft[,1])
     f<-f[pq,] #Ordena pela Fitness
     ft<-ft[pq,]
     g<-g[pq,]
     fmin<-ft[1,1]
     if (fmin<fbest)
       {fbest<-fmin
        gbest<-g[1,]
        solution<-ft[1,]
        loop<-i
        #cat("Solu??o na Gera??o ",i," = ",round(fbest,8),"\n")
        flush.console()
       }
     felite<-f[1:pelite,]
     fnonelite<-f[(pelite+1):size_population,] #N?o Elite
     fmutant<-pop_gen(psize=pmutant,numero_objetos)#Mutantes
     fnovos<-crossover_uniforme(felite,fnonelite,pr,pelite,(size_population-pelite),(size_population-pelite-pmutant),numero_objetos,sq)
     f<-rbind(fmutant,fnovos)
     g<-rbind(g[1:pelite,],decoder(f,clusters,size_population))
     f<-rbind(felite,f)
     glk<-g[(pelite+1):size_population,]
     ftk<-t(apply(glk,1,function(x) fitness(x,distk,clusters,numero_objetos)))
     ft<-rbind(ft[1:pelite,],ftk)
     #####################################################################################
     u<-unique(ft[,1])
     rx<-length(u)
     if (rx<length(ft[,1]))
       {cx<-length(ft[,1])-rx
        fn<-pop_gen(psize=cx,numero_objetos)
        px<-match(u,ft[,1])
        pz<-setdiff(c(1:size_population),px)
        gn<-decoder(fn,clusters,cx)
        ftn<-t(apply(gn,1,function(x) fitness(x,distk,clusters,numero_objetos)))
        ft[pz,]<-ftn
        g[pz,]<-gn
        f[pz,]<-fn
       }
     setTxtProgressBar(pb,i)
     
    }
tempo=round((proc.time()-tempo)[[3]],5)
tempo1<-proc.time()
x<-path_relinking(g,ft[,1],clusters,numero_objetos,distk,round(pelite*size_population))
medoids_finais<-x$a
oi<-as.matrix(1:numero_objetos)
clustering<-cbind(oi,apply(oi,1,function(x) which.min(distk[x,medoids_finais])))
tempo1=round((proc.time()-tempo1)[[3]],5)
solution<-list(fobj=solution[1],meds=solution[(2:(clusters+1))],tempo=tempo,loop=loop,distk=distk,
f=f,ft=ft,fpath=x,tempo1=tempo1,clustering=clustering)
#cat("Solu??o BRKGA           ",solution$fobj,"Tempo = ",tempo,"\n")
#cat("Solu??o Path Relinking  ",solution$fpath$b[1],"Tempo = ",tempo1,"\n")
flush.console()

return(solution)
}



#a<-BKGA_medoids_af("J:","DS1-200DATA.i",3)