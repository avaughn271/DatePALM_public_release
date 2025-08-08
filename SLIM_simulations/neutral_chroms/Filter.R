A = read.table("temp.txt", header = T)
B = A[[2]]
for (i in 1:length(B)) {
  B[i] = strsplit(B[i], ":")[[1]][2]
}
B = as.numeric(B)
B = round(B-5000, -4)

 
indices = c()
for (bloc in unique(B)) {
  temp = (which(B == bloc))
  if (length(temp) == 1 ) {
    indices = c(indices, temp)
  } else {
    temp2 = temp[1]
    val = abs(A$beta[temp2] /  A$se[temp2]  )
    for (kk in 2:length(temp)) {
      if (abs(A$beta[temp[kk]] /  A$se[temp[kk]]  )  > val) {
        val =  abs(A$beta[temp[kk]] /  A$se[temp[kk]]  ) 
        temp2 = temp[kk]
      } 
    }
    indices = c(indices, temp2)
  }
}
A = A[indices,]
A = A[A$pval <= 5 * 10**(-8),]
write.table(A, "temp2.txt" , row.names = F, quote = F, sep ="\t")