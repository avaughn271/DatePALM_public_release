jointfile =  read.table("temp2.txt", header = T)
 
for (row in 1:nrow(jointfile)) {
  if (jointfile$derived_allele[row] == jointfile$minor_allele[row]) {
    allelefreq = jointfile$minor_AF[row]
  }
  else {
    allelefreq = 1.0 - jointfile$minor_AF[row]
  }
  name = jointfile$variant[row]
  name = gsub(":","_", name)
  writeLines(toString(allelefreq), paste0("ModernFreq/", name, ".txt"))
}
