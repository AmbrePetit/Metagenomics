#!/bin/bash

# Define parameters
fastq_dir=$1
cutadapt_dir=$2
figaro_dir=$3
fastq_no_primer=$4
type_data=$5
primer_removed=$6
output_path=$7 #must be the same as the path given to the script nettoyage_metagenomique.R

# Parameters for cutadapt mean
total_sum_length=0
total_number_seq=0

# Parameters for figaro
figaro="/Users/ambre/figaro/figaro/"
amplicon_length=440
forward_primer_length=17
reverse_primer_length=21

# Parameters for cutadapt to remove primer
# Modify them if you use different primer
forward_primer_16S="CCTACGGGNGGCWGCAG"
reverse_primer_16S="GACTACHVGGGTATCTAATCC"
forward_primer_18S="TGCGGCTTAATTYGACTCAAC"
reverse_primer_18S="GCATCACAGAYCTGTT"

####CUTADAPT####
# The aim is to cut all sequences to the same length

if [[ $type_data == "16S" ]]
then
    echo "---------Type of gene studied : $type_data"
    
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
        file_sortie="$cutadapt_dir/$(basename "$file")"
        # Cutadapt is used to cut sequences at the same length, defined previously
        cutadapt -l "$total_mean" -m "$total_mean" -o "$file_sortie" "$file"
        echo "File $file processed and new file saved in $file_sortie."
    done


####FIGARO####
# The aim is to have the parameters form filterAndTrim() function

    python3 $figaro/figaro.py -i $cutadapt_dir -o $figaro_dir -a $amplicon_length -f $forward_primer_length -r $reverse_primer_length
fi


####CUTADAPT####
# The aim is to cut primer from sequences
## Condition : ask if primer already been removed, if yes don't do this step

if [[ $primer_removed == "false" ]]
then
    
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
            cutadapt -g $forward_primer_18S -G $reverse_primer_18S -o "$output_R1" -p "$output_R2" "$file_R1" "$file_R2"
        else
            cutadapt -g $forward_primer_16S -G $reverse_primer_16S -o "$output_R1" -p "$output_R2" "$file_R1" "$file_R2"
        fi
    done
    echo "Primers were removed with cutadapt. " >> "$output_path/output_parameters.txt"
fi
