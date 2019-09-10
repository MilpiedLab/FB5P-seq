#!/bin/sh

# Prepare file containing the list of plates and samples to be analyzed
# Author: DONG Chuang CIML/PMlab 05.16.2019

if [ $# -lt 1 ] ; then
	echo  "\nAuthor: DONG Chuang CIML/PMlab 05.16.2019"
	echo  "USAGE: $0 <command>"
	echo  "Use for Facs based 5' prime End Sc RNAseq Project"
	echo  "Prepare file containing the list of plates and samples to be analyzed"
	echo  "For Example	$0 custom_180416_h_HuPhysioB_1_metadata.csv\n"
	exit 1
fi

numPlate=0
nameSample=""

while IFS=, read -r col1 col2 col3
do	
	if [ $col1 == "Sample" ]
	then
		nameSample=$col3
	
	elif [ $col1 == "Plate" ] 
	then
		numPlate=$col3
	fi
done < $1

#echo -e "Sample name of Metadata is $nameSample\n"
#echo -e "Number plates of Metadata is $numPlate\n"

mkdir Input
cd Input

here=$PWD
cpt=0
for i in `ls ../Rawdata/*.gz | cut -d"/" -f3`
	do
		#echo $i
		if [[ $i == *R1_001.fastq.gz* ]];
		then
  			fold_Name=`echo $i|cut -d"_" -f 2`
			fold_Name="plate"`echo ${fold_Name:1:2}`
			mkdir $fold_Name
			cpt=$(($cpt + 1))
			cd $fold_Name
			#echo "partial match"
			ln -s ../../Rawdata/$i ./
		else
			#echo "no match"
			ln -s ../../Rawdata/$i ./
			cd $here
		fi
	done

if [ $cpt == $numPlate ]
then
	echo "successsssssssssssssssssssful"
else
	echo "nomber of plate is difference with metadata"
	echo -e "Number plates of Metadata is $numPlate\n"
	echo -e "Number plates of built is $cpt\n"
fi

wait
sync
sleep 1
cd ..

ls Input/*/*fastq.gz | cut -d"/" -f2,3 | cut -d"_" -f1-4 | uniq |\
awk '{split($0, res,"/" ); print res[1],"\t",res[2]}' > config/plates_samples_list.txt

echo "successsssssssssssssssssssful" 
