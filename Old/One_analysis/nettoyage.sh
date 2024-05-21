#!/bin/bash

# Define parameters
fastq_dir=$1
type_data=$2
primer_removed=$3
output_path=$4
type_database=$5
forward_primer=$6
reverse_primer=$7
# Just for 16S or 18S
amplicon_length_variable=$8 # taille amplicon variable ? "Yes" or "No" expected
# Just for 16S or 18S non variable
amplicon_length=$9

# Parameters for cutadapt mean
total_sum_length=0
total_number_seq=0

# Parameters to remove primer
# Modify them if you use different primer
#forward_primer_16S="CCTACGGGNGGCWGCAG"
#reverse_primer_16S="GACTACHVGGGTATCTAATCC"
#forward_primer_18S="TGCGGCTTAATTYGACTCAAC"
#reverse_primer_18S="GCATCACAGAYCTGTT"
#forward_primer_ITS="CTTGGTCATTTAGAGGAAGTAA"
#reverse_primer_ITS="GCTGCGTTCTTCATCGATGC"

# Parameters for figaro
figaro="/Users/ambre/figaro/figaro/"
#amplicon_length_16S=440
#amplicon_length_18S=443
forward_primer_length=${#forward_primer}
reverse_primer_length=${#reverse_primer}
figaro_update="TRUE"

# Parameters for R script to clean data with DADA2
rscript="/usr/local/bin"
script_dada2_dir="/Users/ambre/Desktop"


# IPG
#rscript="/usr/bin"
#script_dada2_dir="/media/linuxipg/1.42.6-8754"
#figaro="/media/linuxipg/1.42.6-8754/Isaure/figaro/figaro"


#### CUTADAPT ####
# The aim is to cut all sequences to the same length

if [[ ($type_data == "16S" || $type_data == "18S") && $amplicon_length_variable == "No" ]]
then
    echo "---------Type of gene studied : $type_data"
    
    cutadapt_dir=$(readlink -f "$output_path/cutadapt")
    figaro_dir=$(readlink -f "$output_path/figaro")
    
    # Check if directory exists, if not create it
    if [ ! -d "$cutadapt_dir" ]; then
        mkdir -p "$cutadapt_dir"
    fi

    if [ ! -d "$figaro_dir" ]; then
        mkdir -p "$figaro_dir"
    fi
    
    for file in $fastq_dir/*.fastq.gz;
    do
        # Check if there are matching files
        [ -e "$file" ] || continue
        length=$(gzcat "$file" | awk 'NR%4 == 2 {print length($0)}') #length of sequences for the file
        sum_length=$(gzcat "$file" | awk 'NR%4 == 2 {sum+=length($0)} END {print sum}') #sum of sequence lengths for the file
        number_seq=$(gzcat "$file" | awk '!(NR % 4) {k++}END{print k}') #number of sequences for the file
        mean=$((sum_length/number_seq)) #mean of length for the file
        
        total_sum_length=$((total_sum_length + sum_length))
        total_number_seq=$((total_number_seq + number_seq))
    
        echo "Number of sequences in the file $file : $number_seq" >> "$output_path/output_parameters.txt"
        echo "Average length of sequences in the file $file : $mean"
        echo "Average length of sequences in the file $file : $mean" >> "$output_path/output_parameters.txt"
    done
    
    total_mean=$((total_sum_length / total_number_seq))

    echo "Total average length of sequences in all files : $total_mean"

    echo "Length at which the sequences was cut : $total_mean" >> "$output_path/output_parameters.txt"

    for file in $fastq_dir/*.fastq.gz
    do
        output_file="$cutadapt_dir/$(basename "$file")"
        # Cutadapt is used to cut sequences at the same length, defined previously
        cutadapt -l "$total_mean" -m "$total_mean" -o "$output_file" "$file"
        echo "File $file processed and new file saved in $output_file."
    done


#### FIGARO ####
# The aim is to have the parameters form filterAndTrim() function

    python3 $figaro/figaro.py -i $cutadapt_dir -o $figaro_dir -a $amplicon_length -f $forward_primer_length -r $reverse_primer_length
    
    if [ $? -ne 0 ]
    then
        figaro_update="FALSE"
    fi
fi


#### CUTADAPT ####
# The aim is to cut primer from sequences
## Condition : ask if primer already been removed, if yes don't do this step

if [[ $primer_removed == "No" && ($type_data == "16S" || $type_data == "18S") ]]
then
    fastq_no_primer=$(readlink -f "$output_path/fastq_no_primer")
    # Check if there are matching files
    if [ ! -d "$fastq_no_primer" ]; then
        mkdir -p "$fastq_no_primer"
    fi
    
    for file_R1 in $fastq_dir/*_R1*.fastq.gz
    do
        #Find file R2 associated with R1 (with file name)
        file_R2=$(echo "$file_R1" | sed 's/_R1/_R2/g')
        echo "---------File $file_R1 associated with file $file_R2 being processed..."
    
        # Define output files
        output_R1="$fastq_no_primer/$(basename "$file_R1")"
        output_R2="$fastq_no_primer/$(basename "$file_R2")"

        if [[ $type_data == "18S" ]]
        then
            cutadapt -g $forward_primer -G $reverse_primer -o "$output_R1" -p "$output_R2" "$file_R1" "$file_R2"
        elif [[ $type_data == "16S" ]]
        then
            cutadapt -g $forward_primer -G $reverse_primer -o "$output_R1" -p "$output_R2" "$file_R1" "$file_R2"
        fi
    done
    echo "Primers were removed with cutadapt. " >> "$output_path/output_parameters.txt"
fi


#### R SCRIPT DADA2 ####

if [[ ($type_data == "16S" || $type_data == "18S") && $amplicon_length_variable == "No" && $primer_removed == "Yes" ]]
then
    echo "Command line used : $rscript/Rscript $script_dada2_dir/nettoyage_metagenomique_ITS.R -f $fastq_dir -t $type_data -o $output_path --type_database $type_database -j $figaro_dir -p $primer_removed --figaro $figaro_update -l $amplicon_length_variable -a $amplicon_length --primer_forward $forward_primer --primer_reverse $reverse_primer" >> "$output_path/output_parameters.txt"
    $rscript/Rscript $script_dada2_dir/nettoyage_metagenomique_ITS.R -f $fastq_dir -t $type_data -o $output_path --type_database $type_database -j $figaro_dir -p $primer_removed --figaro $figaro_update -l $amplicon_length_variable -a $amplicon_length --primer_forward $forward_primer --primer_reverse $reverse_primer
    
elif [[ ($type_data == "16S" || $type_data == "18S") && $amplicon_length_variable == "No" && $primer_removed == "No" ]]
then
    echo "Command line used : $rscript/Rscript $script_dada2_dir/nettoyage_metagenomique_ITS.R -f $fastq_no_primer -t $type_data -o $output_path --type_database $type_database -j $figaro_dir -p $primer_removed --figaro $figaro_update -l $amplicon_length_variable -a $amplicon_length --primer_forward $forward_primer --primer_reverse $reverse_primer" >> "$output_path/output_parameters.txt"
    $rscript/Rscript $script_dada2_dir/nettoyage_metagenomique_ITS.R -f $fastq_no_primer -t $type_data -o $output_path --type_database $type_database -j $figaro_dir -p $primer_removed --figaro $figaro_update -l $amplicon_length_variable -a $amplicon_length --primer_forward $forward_primer --primer_reverse $reverse_primer
    
elif [[ ($type_data == "16S" || $type_data == "18S") && $amplicon_length_variable == "Yes" && $primer_removed == "Yes"  ]]
then
    echo " Command line used : $rscript/Rscript $script_dada2_dir/nettoyage_metagenomique_ITS.R -f $fastq_dir -t $type_data -o $output_path --type_database $type_database -p $primer_removed -l $amplicon_length_variable --primer_forward $forward_primer --primer_reverse $reverse_primer" >> "$output_path/output_parameters.txt"
    $rscript/Rscript $script_dada2_dir/nettoyage_metagenomique_ITS.R -f $fastq_dir -t $type_data -o $output_path --type_database $type_database -p $primer_removed -l $amplicon_length_variable --primer_forward $forward_primer --primer_reverse $reverse_primer
elif [[ ($type_data == "16S" || $type_data == "18S") && $amplicon_length_variable == "Yes" && $primer_removed == "No"  ]]
then
    echo " Command line used : $rscript/Rscript $script_dada2_dir/nettoyage_metagenomique_ITS.R -f $fastq_no_primer -t $type_data -o $output_path --type_database $type_database -p $primer_removed -l $amplicon_length_variable --primer_forward $forward_primer --primer_reverse $reverse_primer" >> "$output_path/output_parameters.txt"
    $rscript/Rscript $script_dada2_dir/nettoyage_metagenomique_ITS.R -f $fastq_no_primer -t $type_data -o $output_path --type_database $type_database -p $primer_removed -l $amplicon_length_variable --primer_forward $forward_primer --primer_reverse $reverse_primer
elif [[ $type_data == "ITS" ]]
then
    echo " Command line used : $rscript/Rscript $script_dada2_dir/nettoyage_metagenomique_ITS.R -f $fastq_dir -t $type_data -o $output_path --type_database $type_database -p $primer_removed --primer_forward $forward_primer --primer_reverse $reverse_primer " >> "$output_path/output_parameters.txt"
    $rscript/Rscript $script_dada2_dir/nettoyage_metagenomique_ITS.R -f $fastq_dir -t $type_data -o $output_path --type_database $type_database -p $primer_removed --primer_forward $forward_primer --primer_reverse $reverse_primer
fi
