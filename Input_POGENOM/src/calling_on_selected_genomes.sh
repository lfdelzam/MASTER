#!/bin/bash -l

err_report() {
    echo "Error on line $1 - script sorting_vf.sh"
    exit 1
}

trap 'err_report $LINENO' ERR

#---arguments

wdir=$1
MAGs_path=$2
threads=$3
dts=$4
pdir=$5

#---- Main
declare -a mags
mags=()
cd $wdir/$MAGs_path/
mags=($(ls -d *)) # list of Genomes
cd $wdir

for i in "${mags[@]}"
do
if [ "$(ls -A $wdir/$MAGs_path/$i)" ]  # if directory is not empty
  then
    files=($(ls $wdir/$MAGs_path/$i/*.bam)) # list of BAM files for this genome
    if (( $(echo "${#files[@]} > 1" | bc -l) )) # If there are more than 1 bam file (several samples)
       then
            snakemake -s snakefiles/step2_calling_several_BAMs --config mag_name="$i" -j $threads -F  # It merges BAM files and generates VCF files

    elif (( $(echo "${#files[@]} == 1" | bc -l) ))
        then
          snakemake -s snakefiles/step2_B_calling_one_BAM --config mag_name="$i" -j $threads -F  # It generates VCF files
    fi
else
  mkdir -p $wdir/06_VCF/$dts/$pdir
  echo "The genome $i has not BAM file that passes the filter breadth and coverage. A vcf file cannot be created" > $wdir/06_VCF/$dts/$pdir/$i"_samples.txt"
  echo "The genome $i has not BAM file that passes the filter breadth and coverage. Directory $MAGs_path/$i is empty, and a vcf file cannot be created"
fi
done
