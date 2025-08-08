A = readLines("directional.txt")
separate = c(grep("#OUT", A), length(A) + 1)

#individual i contains genomes 2i and 2i + 1


temppos = read.table("temp2.txt", header =T)[[2]]

for (i in 1:length(temppos)) {
  temppos[i] = strsplit(temppos[i], ":")[[1]][2]
}
temppos = as.numeric(temppos)

B = list()
for (j in 2:length(separate)) {
  B[[j - 1]] = A[separate[j - 1]:(separate[j] - 1)]
}
B[[length(B) - 1]] = B[[length(B) - 1]][1:(length(B[[length(B) - 1]]) - 1)]

FULLPOPULATION = B[[length(B)]]
mut = which(FULLPOPULATION == "Mutations:")
indiv = which(FULLPOPULATION == "Individuals:")
genom = which(FULLPOPULATION == "Genomes:")
mutationss =  FULLPOPULATION[(mut+1):(indiv-1)]
genomess = FULLPOPULATION[(genom+1):length(FULLPOPULATION)]
mutationmatrix = matrix("", nrow = length(mutationss), ncol  = 9)
for (i in 1:length(mutationss)) {
  mutationmatrix[i,] = strsplit(mutationss[i], " ")[[1]]
}
colnames(mutationmatrix) = c("tempid", "permid", "mtype", "pos", "effect", 
                             "dom", "pop", "origintime", "copies")
mutationmatrix = mutationmatrix[,-c(3,6,7, 8)]
mode(mutationmatrix) = "numeric"
Largemutationmatrix = mutationmatrix[order(mutationmatrix[,1]), ]

checkk = Largemutationmatrix[,1] + 1
if (0 != sum(checkk - 1:length(checkk))) {print("PROBLEM")}

print(which(table(Largemutationmatrix[,3]) > 1 ))

for (mutationn in 1:nrow(Largemutationmatrix)) {
  
  if (Largemutationmatrix[mutationn,3] %in% temppos) {
  
  
  
  
  
  countss = Largemutationmatrix[mutationn,5]
  tempp = countss/length(genomess)
  if (tempp > 0.5) {
    tempp = 1 - tempp
  }
  if (tempp > 0.005) {
  permid = Largemutationmatrix[mutationn, 2]
  
  A = c()
  toprint = ""
  for (timeindex in 1:(length(B) - 1)) {
    mutationss =  B[[timeindex]][3:(length(B[[timeindex]]) - 3) ]
    mutationmatrix = matrix("", nrow = length(mutationss), ncol  = 9)
    for (i in 1:length(mutationss)) {
      mutationmatrix[i,] = strsplit(mutationss[i], " ")[[1]]
    }
    colnames(mutationmatrix) = c("tempid", "permid", "mtype", "pos", "effect", 
                                 "dom", "pop", "origintime", "copies")
    mutationmatrix = mutationmatrix[,-c(3,6,7, 8)]
    mode(mutationmatrix) = "numeric"
    mutationmatrix = mutationmatrix[order(mutationmatrix[,1]), ]
    
    checkk = mutationmatrix[,1] + 1
    if (0 != sum(checkk - 1:length(checkk))) {print("PROBLEM")}
    
    genome1 = as.numeric(strsplit(B[[timeindex]][length(B[[timeindex]])], " ")[[1]][-c(1,2)]) + 1
    genome2 = as.numeric(strsplit(B[[timeindex]][length(B[[timeindex]]) - 1], " ")[[1]][-c(1,2)]) + 1
    
    genome1permids = c()
    copies = 0
    if (permid %in% mutationmatrix[genome1, 2]) {
      copies = copies + 1
    }
    if (permid %in% mutationmatrix[genome2, 2]) {
      copies = copies + 1
    }
    if (copies == 2) {
      toprint = paste(toprint,  toString(length(B) - timeindex + 0.001) ,  "-Inf", "0.00", sep = " ")
      toprint = paste0(toprint, "\n")
      toprint = paste(toprint,  toString(length(B) - timeindex + 0.001) ,  "-Inf", "0.00", sep = " ")
      toprint = paste0(toprint, "\n")
    } else if (copies == 0) {
      toprint = paste(toprint,  toString(length(B) - timeindex + 0.001) ,   "0.00","-Inf", sep = " ")
      toprint = paste0(toprint, "\n")
      toprint = paste(toprint,  toString(length(B) - timeindex + 0.001) ,  "0.00", "-Inf", sep = " ")
      toprint = paste0(toprint, "\n")
    } else {
      toprint = paste(toprint,  toString(length(B) - timeindex + 0.001) ,  "-Inf",  "0.00", sep = " ")
      toprint = paste0(toprint, "\n")
      toprint = paste(toprint,  toString(length(B) - timeindex + 0.001) ,  "0.00", "-Inf", sep = " ")
      toprint = paste0(toprint, "\n")
    }
  }
 writeLines(toprint, paste0("SNPs/" ,  toString(permid), "_" ,  toString(Largemutationmatrix[mutationn, 3]),  "_A_G.txt")   )
}}
}
