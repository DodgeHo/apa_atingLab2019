#!/bin/bash -l
echo "Make sure that fastq-dump in $PATH and configure it appropriately"

wkd=01_wkd
fastqd=${wkd}/fastq

echo "Check if fastq directory exists ..."
if [ ! -d $fastqd ]; then
	mkdir -p $fastqd
	echo "created ${fastqd}"
fi

rm -rf ${fastqd}/*.fastq.gz

echo "Start to download PolyA-seq FASTQ files ..."

#	- SRP083252 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86178)
#	- SRP083254 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86180)

srrs=(SRR4091084)

for srr in "${srrs[@]}"
do
  echo "fastq-dump -F -O ${fastqd} $srr"
  fastq-dump -F -O ${fastqd} $srr | pv -p -t -e -r -b > /dev/null
done

echo "Done."
echo "Check if FASTQ files are located in ${fastqd}"
