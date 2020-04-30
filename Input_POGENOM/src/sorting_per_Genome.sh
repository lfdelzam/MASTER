#!/bin/bash -l

err_report() {
    echo "Error on line $1 - script sorting_per_Genome.sh"
    exit 1
}

trap 'err_report $LINENO' ERR

#---arguments

i=$1
wdir=$2
MAGs_path=$3
threads=$4
dts=$5
pdir=$6

# --- Main ---

if [ "$(ls -A $wdir/$MAGs_path/$i)" ] # if directory is not empty
  then
    files=($(ls $wdir/$MAGs_path/$i/*.bam))  # list of bam files
    if (( $(echo "${#files[@]} > 1" | bc -l) )) # if there are more than 1 bam file
       then
          snakemake -s snakefiles/step2 --config mag_name="$i" -j $threads -F --quiet  # It merges BAM files and generates VCF files

    elif (( $(echo "${#files[@]} == 1" | bc -l) ))
        then
          snakemake -s snakefiles/step2_B --config mag_name="$i" -j $threads -F --quiet  # It generates VCF files
    fi
else
  mkdir -p $wdir/06_VCF/$dts/$pdir
  echo "The genome $i has not BAM file that passes the filter breadth and coverage. A vcf file cannot be created" > $wdir/06_VCF/$dts/$pdir/$i"_samples.txt"
  echo "The genome $i has not BAM file that passes the filter breadth and coverage. Directory $MAGs_path/$i is empty, and a vcf file cannot be created"
fi
