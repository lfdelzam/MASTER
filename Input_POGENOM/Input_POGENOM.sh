#!/bin/bash -l
wd=$(pwd)
err_report() {
    echo "Error on line $1 - script Input_POGENOM.sh"
    cd $wd
    if [ -z "$mag" ]; then mag=$(fullname=$(basename $(ls RAW_DATA/Genomes/$dataset/*$genomes_ext | head -1)) ; echo ${fullname%.*}); fi
    if [ -z "$samples" ]; then samples=$(fulln=$(basename $(ls RAW_DATA/Reads/$dataset/*$fwd_index$reads_ext | head -1)) ; echo ${fulln%$fwd_index$reads_ext}); fi
    mess1="Check if the 'ip_env' has been activated - command:\n    'conda activate ip_env'"
    mess2="Use the command:\n    'snakemake -s snakefiles/step1_pogenom_input --unlock'"
    mess3="Use the command:\n    'snakemake -s snakefiles/step_pogenom_input --config my_mag='"$mag"' my_samples='"$samples"' --unlock'"
    mess="or\n    'snakemake -s snakefiles/step2 --config mag_name='"$mag"' --unlock' or\n    'snakemake -s snakefiles/step2_B --config mag_name='"$mag"' --unlock'"
    mess_inst="If you haven't installed the required software, please $mess1\n    Then install software with the command:\n     conda install -c bioconda -c conda-forge prodigal hmmer"
    if [ "$1" == 86 ]; then echo -e "if you are using conda, check if the 'ip_env' has been activated - command:\n    'conda activate ip_env' "; fi
    if [ "$1" == 88 ]; then echo -e "TIP 1 - Look at $wd/log_files/samples_filter_$dataset.log\nTIP 2 - $mess1\nTIP 3 - Use the command:\n    'snakemake -s snakefiles/step_filter --unlock'\n     and run the pipeline again"; fi
    if [ "$1" == 104 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset.$mag.coverage_breadth.log\nTIP 2 - $mess1\nTIP 3 - $mess3\n     and run the pipeline again"; fi
    if [ "$1" == 106 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset.$mag"_vcf_files.log"\nTIP 2 - $mess1\nTIP 3 - $mess3 $mess\n     and run the pipeline again"; fi
    if [ "$1" == 109 ] || [ "$1" == 123 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset.$mag"_gff_file.log"\nTIP 2 - $mess_inst\nTIP 3 - Use the command:\n    'snakemake -s snakefiles/Prokaryotic_gene_calling_annotations --config my_mag='"$mag"' --unlock'\n     and run the pipeline again"; fi
    if [ "$1" == 123 ]; then echo -e "TIP 1 - Look at\n    $wd/log_files/$dataset.$mag.coverage_breadth.log or\n    $wd/log_files/$dataset.$mag"_vcf_files.log"\nTIP 2 - $mess1\nTIP 3 - $mess3\n     and run the pipeline again"; fi
    if [ "$1" == 129 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset"_Genomes_coverage_breadth.log"\nTIP 2 - $mess1\nTIP 3 - $mess2\n     and run the pipeline again"; fi
    if [ "$1" == 131 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset"_Genomes_vcf_files.log"\nTIP 2 - $mess1\nTIP 3 - $mess2 $mess\n     and run the pipeline again"; fi
    if [ "$1" == 134 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset"_gff_files.log"\nTIP 2 - $mess_inst\nTIP 3 -Use the command:\n    'snakemake -s snakefiles/Prokaryotic_gene_calling_annotations --config my_mag="all_genomes" --unlock'\n     and run the pipeline again"; fi
    if test -f "temporal"; then rm temporal; fi
    exit 1
}
trap 'err_report $LINENO' ERR

#Default options
configFile=$wd/config_files/Input_POGENOM_config.json
#----Argument parse----------------------
for a in "$@"
do
case $a in
  -d=*|--path_to_config_file=*)
  if [ -z "${a#*=}" ];  then
  echo "value to argument -d No supplied"
  exit 0
  else configFile="${a#*=}"
  fi
  shift # past argument
  ;;

  *)
  echo -e "\nUsage: bash Input_POGENOM.sh [options]\n -d=<absolute path to configFile. Default=$configFile>\n"
  echo -e 'Description:\nThis program executes a pipeline that generates the required input files for POGENOM.\nThe aim of this pipeline is to increase the reproducibility of the data analysis, and to simplify the use of POGENOM.\nPOGENOM is a computer program that calculates several population genetic parameters for a genome in relation to a set of samples (https://github.com/EnvGen/POGENOM).'
  exit 0
  ;;

esac
done
if [[ "$configFile" != /* ]]; then
    echo "Please provide an absoltute path to configfile e.g., bash Input_POGENOM.sh '/absolute/path/to/configfile' "; exit 1
fi

cat $configFile | sed s/"[{|}]"//g | sed s/":"/"="/g | sed s/",$"//g | sed s/'"'//g | sed s/" "//g > temporal
. temporal

#Checking key parameters setting
if [[ "$workdir" != /* ]]; then
    echo "Please provide an absoltute path to the working directory in configfile e.g., 'workdir': '/absolute/path/to/working_directory/' "
    echo "Please double-check the absolute path to working directory $workdir and to configfile $configFile"
    rm temporal
    exit 1
fi

options=("$dataset" "$min_coverage" "$min_breadth" "$min_bsq_for_cov_median_calculation" "$threads" "$genomes_ext" "$reads_ext" "$fwd_index" "$rev_index" "$bowtie2_params" "$mapqual" "$freebayes_parameters" "$vcffilter_qual")
for o in "${options[@]}"; do if [ -z "$o" ]; then echo "A key parameter is undefined, please check in the config_files/Input_POGENOM_config.json file the parameters used"; rm temporal; exit 1; fi; done

if  [[ "$mode" == prefilt ]] && ( [ -z "$fraction" ] || [ -z "$temp_sub_Reads_dir" ] ); then
    echo 'A key parameter in "mode" : "prefilt" is undefined, please check in the config_files/Input_POGENOM_config.json file the parameters used'; rm temporal; exit 1
fi

if [ "$annotation" == yes ]; then
    if [[ "$pfam_db_path" != /* ]]; then echo "Please provide an absoltute path to the Pfam database in configfile"; rm temporal; exit 1; fi
    if [ -z "$evalue_pfam" ]; then echo "'evalue_pfam' parameter for annotation is undefined, please check in the config_files/Input_POGENOM_config.json file the parameters used"; rm temporal; exit 1; fi
fi

echo "INFO: Starting Input_POGENOM pipeline - Working directory: $workdir"
if [[ $snakemake_extra_params == *","* ]]; then extra_params=$( echo $snakemake_extra_params | sed s/","/" "/g); else extra_params=$snakemake_extra_params; fi
mkdir -p $workdir/log_files

# ----Using prefilt mode - full workflow - With Sample pre-screening
# --- main - mode prefilt
if  [[ "$mode" == prefilt ]]; then
         cd $workdir # Sample pre-screening
         echo "INFO: Generating Reads subsets - Fraction used $fraction"
         bash src/create_prefilt_Reads_subdir.sh $fraction $genomes_ext $reads_ext $temp_sub_Reads_dir $dataset
         echo "INFO: Calculating Genome Median coverage - sub-samples - Median coverage threshold $min_coverage"
         snakemake -s snakefiles/step_filter -j $threads $extra_params 2>log_files/samples_filter_$dataset.log

         if [[ "$remove_subreads" == yes ]] && test -d "$temp_sub_Reads_dir/Reads"; then
              echo "WARNING: You have chosen to remove $temp_sub_Reads_dir/Reads/"
              rm -rf $temp_sub_Reads_dir/Reads/
         fi
         result_dir="PREFILT/"$dataset"/params_cov_"$min_coverage"_mpq_"$mapqual"_bq_"$min_bsq_for_cov_median_calculation"_fr_"$fraction
         file_empty=$(grep -v "#" $result_dir/Selected_samples_Genomes.txt | wc -l)  # File generated by snakefiles/step_filter (line 88) and contains the list of Genomes and selected samples
             if [ "$file_empty" -eq 0 ]; then
                echo -e "INFO: With the current parameter setting: Dataset $dataset - Fraction $fraction - Median coverage threshold $min_coverage - Min-base quality $min_bsq_for_cov_median_calculation - Mapping quality $mapqual\n      There is no Genome - sample with Estimated Median Coverage higher than threshold.\n      A vcf file cannot be created\n"
             else
                echo "INFO: Calculating Genome Median coverage and breadth - Dataset: $dataset - Median coverage threshold: $min_coverage - Breadth threshold: $min_breadth %"
                grep -v "#" $result_dir/Selected_samples_Genomes.txt | while read line  # Each line contains the Genome name and it corresponding list of selected samples
                do # Mapping and coverage/breadth calculations
                  mag=$(echo $line | cut -d " " -f1) # Genome name
                  samples=$(echo $line | cut -d " " -f2)  #list of selected samples
                  snakemake -s snakefiles/step_pogenom_input step1_all --config my_mag="$mag" my_samples="$samples" -j $threads $extra_params 2> log_files/$dataset.$mag.coverage_breadth.log
                  echo "INFO: Generating VCF files - Genome $mag"
                  snakemake -s snakefiles/step_pogenom_input vcf --config my_mag="$mag" my_samples="$samples" -j $threads $extra_params 2> log_files/$dataset.$mag"_vcf_files.log"
                  if [ "$annotation" == yes ]; then
                        echo "INFO: Generating gff file - Genome $mag"
                        snakemake -s snakefiles/Prokaryotic_gene_calling_annotations --config my_mag="$mag" -j $threads $extra_params 2> log_files/$dataset.$mag"_gff_file.log"
                  fi
                done
             fi
             no_genome=$(grep "#" $result_dir/Selected_samples_Genomes.txt | wc -l) # List of Genomes without samples ( coverage and breadth below thresholds)
             if [ "$no_genome" -ne 0 ]; then
                 echo "**********************************************"
                 echo "The following Genome(s) has(have) not been analysed"
                 grep "#" $result_dir/Selected_samples_Genomes.txt
                 echo -e "**********************************************\n"
             fi
         echo "INFO: Input_POGENOM pipeline is done !!!"
rm temporal
exit 0
fi
# --- End of prefilt mode

# --- Option when analysing a dataset without prefilt (without Sample pre-screening)
cd $workdir
    echo "INFO: Calculating Genome Median coverage and breadth - Dataset: $dataset - Median coverage threshold: $min_coverage - Breadth threshold: $min_breadth %"
    snakemake -s snakefiles/step1_pogenom_input step1_all -j $threads $extra_params 2> log_files/$dataset"_Genomes_coverage_breadth.log"
    echo "INFO: Generating VCF files - Dataset $dataset"
    snakemake -s snakefiles/step1_pogenom_input vcf -j $threads $extra_params 2> log_files/$dataset"_Genomes_vcf_files.log"
    if [ "$annotation" == yes ]; then
       echo "INFO: Generating gff files - Dataset $dataset"
      snakemake -s snakefiles/Prokaryotic_gene_calling_annotations --config my_mag="all_genomes" -j $threads $extra_params 2> log_files/$dataset"_gff_files.log"
    fi

rm temporal

echo 'INFO: Input_POGENOM pipeline is done !!!'
