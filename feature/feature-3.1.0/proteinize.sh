#!/bin/sh

i=0
for filename in ./data/pdb/*; do
	basename=$(basename -- "$filename")
	file="${basename%%.*}"
	./bin/featurize -l proteins.properties $file > "feature_files/$file.ff"
	#echo $file
	#i=$((i + 1))
	#if [ $i -eq 3 ];
#		then echo good;
#		break;
#	fi;
done

echo
echo All files featurized
echo

