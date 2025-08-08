A = readLines("directional.txt")
separate = c(grep("#OUT", A), length(A) + 1)




#individual i contains genomes 2i and 2i + 1

B = list()
for (j in 2:length(separate)) {
  B[[j - 1]] = A[separate[j - 1]:(separate[j] - 1)]
}

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
mutationmatrix = mutationmatrix[,-c( 3,6,7, 8)]
mode(mutationmatrix) = "numeric"
mutationmatrix = mutationmatrix[order(mutationmatrix[,1]), ]
TraitVector = rep(0, length(genomess))

checkk = mutationmatrix[,1] + 1
if (0 != sum(checkk - 1:length(checkk))) {print("PROBLEM")}

for (j in 1:length(genomess)) { # can think about adding environmental noise here as well.
  a = as.numeric(strsplit(genomess[j], " ")[[1]][-c(1,2)]) + 1
  TraitVector[j] = TraitVector[j] + sum(mutationmatrix[a,4])
}
TraitVectorIndiv = rep(0, length(genomess) / 2)
for (j in 1:(length(  genomess) / 2)) {
  TraitVectorIndiv[j] = TraitVector[j + j] + TraitVector[j + j - 1] 
}



print("a")
MasterList = list()
  for (j in 1:length(genomess)) { 
 MasterList[[j]] = as.numeric(strsplit(genomess[j], " ")[[1]][-c(1,2)]) + 1
 MasterList[[j]] = sort( MasterList[[j]])
}
print("b")

binary_search <- function(x, vec) {
  low <- 1
  high <- length(vec)
  while (low <= high) {
    mid <- floor((low + high) / 2)
    if (vec[mid] == x) return(TRUE)
    else if (vec[mid] < x) low <- mid + 1
    else high <- mid - 1
  }
  return(FALSE)
}


#GenotypeMatrix = matrix(0, nrow = length(genomess), ncol = nrow(mutationmatrix))
#for (j in 1:length(genomess)) {
#  a = as.numeric(strsplit(genomess[j], " ")[[1]][-c(1,2)]) + 1
#  GenotypeMatrix[j, a] = 1
#}
GenotypeMatrixIndiv = rep(0, length(genomess) / 2 )
GenotypeMatrix = rep(0, length(genomess))

towrite = "LD_block\tvariant\tderived_allele\tminor_allele\tminor_AF\tbeta\tse\tpval"
print(nrow(mutationmatrix))
for (j in 1:nrow(mutationmatrix)) {
  countss = mutationmatrix[j,5]
  tempp = countss/length(genomess)
  if (tempp > 0.5) {
    tempp = 1 - tempp
  }
  if (tempp > 0.005) {   #MAF!!!!
  print(j)
  if (mutationmatrix[j,5] != 0 & mutationmatrix[j,5]  != length(genomess)) {
    
    for (genometemp in 1:length(GenotypeMatrix)) {
      GenotypeMatrix[genometemp] = binary_search(j, MasterList[[genometemp]]) + 0
    }
    
    for (indivtemp in 1:length(GenotypeMatrixIndiv)) {
      GenotypeMatrixIndiv[indivtemp] = GenotypeMatrix[indivtemp+indivtemp] + GenotypeMatrix[indivtemp+indivtemp-1] 
    }
    
 model <- lm(TraitVectorIndiv ~ GenotypeMatrixIndiv)
  coeff = (summary(model)$coefficients[2,])
  freqq = (mutationmatrix[j,5]/length(genomess))
  if (freqq < 0.5) {
    temp = paste("G", freqq, sep = "\t")
  } else {
    temp = paste("A" , 1 - freqq, sep = "\t")
  }
  
  towrite = paste(towrite,  paste(
    toString(j) , paste0(mutationmatrix[j,2] , ":" ,mutationmatrix[j,3] , ":A:G"), "G",
  temp, coeff[1], coeff[2],  coeff[4], sep = "\t")  , sep ="\n")

  }
}}

writeLines(towrite, "temp.txt")
