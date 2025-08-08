coll = c(  "#E69F00",  "#56B4E9","#CC79A7", "#009E73")

A = read.csv("isotopedata.csv")
INFO = read.csv("MesoNeoData.tsv", sep = "\t")
b = which(!is.na(A$Sample.ID) &  
            !is.na(A$Material) &  
            !is.na(A[[8]]) & !is.na(A[[9]]) & A[[1]] %in% INFO[[1]])
A = A[b,c(1,2,8,9)]
Bone = A[A[[2]] == "bone",]

AGES = c()
for (i in 1:nrow(Bone)) {
  a = INFO$ageAverage[which(INFO[[1]] == Bone[i,1])]
  AGES = c(AGES, a)
}
Bone = cbind(Bone, AGES)
matrixx = matrix(0, nrow = 0, ncol = 6)
for (i in 1:nrow(Bone)) {
  if (file.exists(paste0("Summaries/", Bone[[1]][i], ".txt" ))) {
    temp = scan(paste0("Summaries/", Bone[[1]][i], ".txt" ))
  } else {
    temp = c(-1,-1,-1,-1,-1,-1)
  }
  matrixx = rbind(matrixx, temp)
}
Bone = cbind(Bone, matrixx[,1:4])
Bone = Bone[which(Bone[[8]] >= 0),]


Bone1 = Bone[which(Bone[[6]] > Bone[[7]] & Bone[[6]] > Bone[[8]] & Bone[[6]] > Bone[[9]]),]
Bone2 = Bone[which(Bone[[7]] > Bone[[6]] & Bone[[7]] > Bone[[8]] & Bone[[7]] > Bone[[9]]),]
Bone3 = Bone[which(Bone[[8]] > Bone[[7]] & Bone[[8]] > Bone[[6]] & Bone[[8]] > Bone[[9]]),]
Bone4 = Bone[which(Bone[[9]] > Bone[[7]] & Bone[[9]] > Bone[[8]] & Bone[[9]] > Bone[[6]]),]
Bone1 = Bone1[which(Bone1[[5]] > 5000 &  Bone1[[5]] < 8000),]
Bone3 = Bone3[which(Bone3[[5]] > 5000 &  Bone3[[5]] < 8000),]
Bone4 = Bone4[which(Bone4[[5]] > 5000 &  Bone4[[5]] < 8000),]

Bone1 = Bone1[which(Bone1[[3]] < -18),]
Bone3 = Bone3[which(Bone3[[3]] < -18),]
Bone4 = Bone4[which(Bone4[[3]] < -18),]



plot(Bone1[[5]] , Bone1[[3]] ,  xlim = c(8000,4800), ylim = c(6,20),
    main = NULL,  pch=19,
     col = coll[1], xlab =  "Years Before Present",  
     ylab = expression(delta^15 * "N% AIR"  ))
points(Bone3[[4]], Bone3[[3]], col =  coll[3],  pch=19)
points(Bone4[[4]], Bone4[[3]], col =  coll[4],  pch=19)


########################################################
library(ggplot2)
library(plotrix)

LIST = list()
LIST[[1]]=Bone1
LIST[[2]]=Bone2
LIST[[3]]=Bone3
LIST[[4]]=Bone4
indexx =1
for (data in LIST) {
  if (indexx != 2){
for (i in 1:nrow(data)) {
  temp = c(data[i,6], data[i,7], data[i,8], data[i,9])
  print(temp)
  draw.circle(as.numeric(data[[5]][i]),as.numeric(data[[4]][i]),radius=300/4, col=coll[indexx])
  
floating.pie(as.numeric(data[[5]][i]), as.numeric(data[[4]][i]), temp, 
             radius = 200/4, col = coll)

}}
  indexx = indexx + 1
}


NEO = sort(c(Bone1[[4]], Bone3[[4]]))
YAM = sort(Bone4[[4]])

pval = ks.test(NEO, YAM)$p.value
legend("bottomright", legend = c("ANA", "CHG", "EHG", "WHG") ,fill = c("#E69F00",  "#56B4E9", "#009E73","#CC79A7"))

par(fig = c(0.492,1, 0.492, 1), new = T)  
 

boxplot(YAM, NEO, ylab = expression(delta^15 * "N% AIR"  ), col = c("white","white"), range= 0,ylim = c(8,17))

temp  = quantile(NEO)

axis(1, at=1:2, labels=c("CHG + EHG", "WHG + ANA"))
rect(1.6+0.005,temp[2]- 0.1, (2.395  + (1.6+0.005))/2 , temp[4] , col = "#CC79A7", border = NA)
rect((2.395  + (1.6+0.005))/2-0.0001,temp[2]- 0.1, 2.395   , temp[4] , col = "#E69F00", border = NA)

temp  = quantile(YAM)

rect(1.6+0.005-1,temp[2]- 0.1, (2.395  + (1.6+0.005))/2 -1, temp[4] , col ="#56B4E9", border = NA)
rect((2.395  + (1.6+0.005))/2-0.0001 -1,temp[2]- 0.1, 2.395 -1  , temp[4] , col = "#009E73", border = NA)



boxplot(YAM, NEO, ylab = expression(delta^15 * "N% AIR"  ), col = c(rgb(0,0,0,0),rgb(0,0,0,0)), add =T, range= 0)


# Add significance bar and p-value
segments(1, max(c(YAM, NEO)) + 1, 2, max(c(YAM, NEO)) + 1)
segments(1, max(c(YAM, NEO)) + 1, 2, max(c(YAM, NEO)) + 1)
segments(1, max(c(YAM, NEO)) + 1, 1, max(c(YAM, NEO)) +0.5)
segments(2, max(c(YAM, NEO)) + 1, 2, max(c(YAM, NEO)) +0.5)

text(1.5, max(c(YAM, NEO)) + 1.5, paste("p =", signif(pval, digits = 3)), cex = .7)

 


