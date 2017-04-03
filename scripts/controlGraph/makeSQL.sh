set -ex

for region in BRCA1 BRCA2 MHC SMA LRC_KIR
do

	echo "Processing ${region}"

	#Shift vcf
	./shiftVCF.py --ref "${region}.fa" --shift 200 --inFile "${region}/${region}.vcf"  --outFile "${region}/${region}_shifted.vcf"

	#Sort (output is vcf)
	vcf-sort ${region}/${region}_shifted.vcf > ${region}/${region}_sorted.vcf

	#Bgzip
	bgzip -c ${region}/${region}_sorted.vcf > ${region}/${region}_sorted.vcf.gz

	#Index
	tabix -f -p vcf ${region}/${region}_sorted.vcf.gz

	#Get contig, start, and end
	contig=`cat ${region}/${region}.bed | awk '{print $1}'`
	start=`cat ${region}/${region}.bed | awk '{print $2+1}'`
	end=`cat ${region}/${region}.bed | awk '{print $3}'`

	#Construct vg
	vg construct -v ${region}/${region}_sorted.vcf.gz -r /cluster/home/hickey/ga4gh/vg/bakeoff/fa_GRCh38/${contig}.fa -R ${contig}:${start}-${end} > ${region}/${region}.vg
	
	#Convert to sql
	vg2sg -s ${region}/${region}.vg ${region}/graph.fa ${region}/graph.sql


done