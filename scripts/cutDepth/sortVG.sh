for regionDir in /hive/users/sblum/ga4gh/kmers/graphs/vg/*
do
	region=${regionDir##*/}
	if [ $region!='CENX' ]
		then
		echo $region
		for inFile in ${regionDir}/*
		do
			noDirInFile=${inFile##*/}
			noDirnoExtInFile=${noDirInFile%.vg}
			outDir=/cluster/home/sblum/ga4gh/cutDepth/sortedVG/${region}
			outFile=${outDir}/${noDirnoExtInFile}-sorted.json
			if [ ! -f "${outFile}" ]
				then
				echo Sorting ${outFile}...
				vg ids -s $inFile | vg view -j - > "$outFile"
			fi
		done
	fi
done