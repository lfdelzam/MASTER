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

mkdir -p $dataset/Reads/$Dts/fraction_$fract

Rds=($(ls $wd/RAW_DATA/Reads/$Dts/*$reads_ext))

for r in "${Rds[@]}"
   do
      read_file=$(basename $r)
       if ! test -s $wd/$dataset/Reads/$Dts/fraction_$fract/$read_file  #If file doesn't exit or if it exists but it is empty
          then
          all_lines_in=$(gzip -cd $r | wc -l)
          all_number_reads=$( echo $all_lines_in/4 | bc )
          subsample=$( echo $all_number_reads*$fract | bc ) # to make sure that paired reads will have equal number of reads. Sometimes, when using directly fraction in seqtk, the final number of reads may be different
          seqtk sample -s100 $r $subsample > $wd/$dataset/Reads/$Dts/fraction_$fract/temponame
          lines_in=$(cat $wd/$dataset/Reads/$Dts/fraction_$fract/temponame | wc -l)
          number_reads=$( echo $lines_in/4 | bc )
          gzip -c $wd/$dataset/Reads/$Dts/fraction_$fract/temponame > $wd/$dataset/Reads/$Dts/fraction_$fract/$read_file
          rm $wd/$dataset/Reads/$Dts/fraction_$fract/temponame
          echo "      Subset $read_file created - Number of reads in subset: $number_reads"
       else
            echo "      Subset $read_file already created"
       fi
   done
