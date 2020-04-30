#!/bin/bash -l

err_report() {
  echo "Error on line $1 - script create_prefilt_Reads_subdir.sh"
  exit 2
}

trap 'err_report $LINENO' ERR

#---arguments---

wd=$(pwd)
fract=$1
mags_ext=$2
reads_ext=$3
dataset=$4
Dts=$5

#--- Main ----

mkdir -p $dataset/Reads/fraction_$fract

Rds=($(ls $wd/RAW_DATA/Reads/$Dts/*$reads_ext))

for r in "${Rds[@]}"
   do
      read_file=$(basename $r)
       if ! test -s $wd/$dataset/Reads/fraction_$fract/$read_file  #If file doesn't exit or if it exists but it is empty
          then seqtk sample -s100 $r $fract > $wd/$dataset/Reads/fraction_$fract/$read_file
               lines_in=$(cat $wd/$dataset/Reads/fraction_$fract/$read_file | wc -l)
  	           number_reads=$( echo $lines_in/4 | bc )
               gzip -q $wd/$dataset/Reads/fraction_$fract/$read_file
               echo "      Subset $read_file created - Number of reads in subset: $number_reads"
       else
            echo "      Subset $read_file already created"
       fi
   done
