
rm -r -f SNPFOLDER
mkdir SNPFOLDER

yourfilenames=`ls SNPs/*.txt`

for eachfile in $yourfilenames
do

fullfilename=$(echo "$eachfile" | cut -d "/" -f 2)

modernfreq=$(head -n 1  ModernFreq/${fullfilename})
echo $modernfreq
 echo $fullfilename


v2=${fullfilename/_A_G.txt/}

python3.10 ~/desktop/palm/palm-master/snp_lik.py  --popFreq  $modernfreq --ancientHaps SNPs/${fullfilename} --out SNPFOLDER/${v2}  --N 100000 --df 600 --tCutoff 531

done