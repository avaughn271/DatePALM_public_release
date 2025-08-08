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



plot(Bone1[[5]] , Bone1[[3]] ,  xlim = c(11000,0), ylim = c(6,20),
    main =  NULL,  pch=19,
     col = coll[1], xlab =  "Years Before Present",  
     ylab =  expression(delta^15 * "N% AIR"  ))
points(Bone2[[4]], Bone2[[3]], col = coll[2],  pch=19)
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
for (i in 1:nrow(data)) {
  temp = c(data[i,6], data[i,7], data[i,8], data[i,9])
  print(temp)
  draw.circle(as.numeric(data[[5]][i]),as.numeric(data[[4]][i]),radius=300/1.2, col=coll[indexx])
  
floating.pie(as.numeric(data[[5]][i]), as.numeric(data[[4]][i]), temp, 
             radius = 200/1.2, col = coll)

}
  indexx = indexx + 1
}
legend("bottomright", legend = c("ANA", "CHG", "EHG", "WHG") ,fill = c("#E69F00",  "#56B4E9", "#009E73","#CC79A7"))

