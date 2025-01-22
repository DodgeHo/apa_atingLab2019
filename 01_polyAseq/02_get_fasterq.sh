#!/bin/bash -l
echo "Make sure that fasterq-dump is in \$PATH and configure it appropriately"

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

#srrs=(SRR4091084 SRR4091085 SRR4091086 SRR4091087 SRR4091088 SRR4091089 SRR4091090 SRR4091091 SRR4091104 SRR4091105 SRR4091106 SRR4091107 SRR4091108 SRR4091109 SRR4091110 SRR4091111 SRR4091113 SRR4091115 SRR4091117 SRR4091119)
srrs=(SRR4091084)
for srr in "${srrs[@]}"
do
  echo "Downloading $srr ..."
  fasterq-dump -O ${fastqd} $srr | pv -p -t -e -r -b > /dev/null
done

echo "Compressing FASTQ files ..."
for srr in "${srrs[@]}"
do
  gzip ${fastqd}/${srr}.fastq
done

echo "Done."
echo "Check if FASTQ files are located in ${fastqd}"
