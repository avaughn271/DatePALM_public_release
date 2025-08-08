Ntrials=30
rm -f -r  Results
mkdir Results
#LOOKS GOOD BUT NEED TO CHANGE MAF IN PROCESSDIRECT
for (( i=1; i<=$Ntrials; i++ ))
do

rm -f  directional.txt
rm -f  temp.txt
rm -f  temp2.txt

slim 13_6_directional.slim > directional.txt

Rscript processdirect.R
Rscript Filter.R

rm -r -f SNPs 
mkdir SNPs
Rscript processtimeseries.R

rm -r -f ModernFreq 
mkdir ModernFreq 
Rscript   getModernDerived.R  

./run_inference.sh  
./run_palm.sh


mv output.txt Results/output${i}.txt

done
