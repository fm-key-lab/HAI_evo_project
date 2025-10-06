# script to extract krakekn results across different taxonomix levels
# repuires lauch in working directory
# results in: 4-bracken/sumBrackTest_*.txt


## Input in Snakefile: expand("2-kraken2/{sampleID}_krakenRep.txt",sampleID=SAMPLE_ls),

## create own input:
#input_samples=`awk -F',' '{if (NR>1) {print $2}}' samples.csv | sort | uniq`
#input=''
#echo ${input_samples}
#
#for subject in ${input_samples}; do
# input+='2-kraken2/'${subject}'_krakenRep_bracken.txt ' 
#done
#echo ${input} 



# Read input from CL
man_switch="off"
while [ $# != "0" ] ; do
 if [[ ${man_switch} == 'off' ]]; then
  input=$1
  man_switch="on"
 else
  input="$input $1"
 fi
shift
done


for f in ${input}; do
 sid=$(basename ${f} '_krakenRep_bracken.txt')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="S"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 4-bracken/sumBrackTest_species.txt

for f in ${input}; do
 sid=$(basename ${f} '_krakenRep_bracken.txt')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="G"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 4-bracken/sumBrackTest_genus.txt

for f in ${input}; do
 sid=$(basename ${f} '_krakenRep_bracken.txt')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="F"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 4-bracken/sumBrackTest_family.txt

for f in ${input}; do
 sid=$(basename ${f} '_krakenRep_bracken.txt')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="O"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 4-bracken/sumBrackTest_order.txt

for f in ${input}; do
 sid=$(basename ${f} '_krakenRep_bracken.txt')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="C"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 4-bracken/sumBrackTest_class.txt

for f in ${input}; do
 sid=$(basename ${f} '_krakenRep_bracken.txt')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="P"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 4-bracken/sumBrackTest_phylum.txt

for f in ${input}; do
 sid=$(basename ${f} '_krakenRep_bracken.txt')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="D"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 4-bracken/sumBrackTest_domain.txt

