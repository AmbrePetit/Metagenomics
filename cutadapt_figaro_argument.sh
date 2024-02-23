#!/bin/bash

# Define parameters
fastq_dir=$1
cutadapt_dir=$2
figaro_dir=$3
fastq_no_primer=$4
type_data=$5
primer_removed=$6

# Parameters for cutadapt mean
min_mean=999999

# Parameters for figaro
amplicon_length=440
forward_primer_length=17
reverse_primer_length=21

# Parameters for cutadapt to remove primer
# Modify them if you use different primer
forward_primer_16S="CCTACGGGNGGCWGCAG"
reverse_primer_16S="GACTACHVGGGTATCTAATCC"
forward_primer_18S="TGCGGCTTAATTYGACTCAAC"
reverse_primer_18S="GCATCACAGAYCTGTT"

# Extract the path of the directory given as parameter
#output_path=$(dirname "$(dirname "$fastq_dir")")
output_path=$(dirname "$fastq_dir")
echo output_path

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
    
        echo "Number of sequences in the file $file : $number_seq" >> "$output_path/output_parameters.txt"
        echo "Average length of sequences in the file $file : $mean"
    
        # Check if the current mean is lower than the current minimum mean
        if [[ $mean < $min_mean ]]
        then
            min_mean="$mean"
        fi
    done

    echo "The smallest average of the files processed is : $min_mean"

    min_mean=$(expr "$min_mean" - 5)

    echo "Length at which the sequences was cut : $min_mean" >> "$output_path/output_parameters.txt"

    for file in $fastq_dir/*.fastq.gz
    do
        file_sortie="$cutadapt_dir/$(basename "$file")"
        # Cutadapt is used to cut sequences at the same length, defined previously
        cutadapt -l "$min_mean" -m "$min_mean" -o "$file_sortie" "$file"
        echo "File $file processed and new file saved in $file_sortie."
    done


####FIGARO####
# The aim is to have the parameters form filterAndTrim() function

    python3 /Users/ambre/figaro/figaro/figaro.py -i $cutadapt_dir -o $figaro_dir -a $amplicon_length -f $forward_primer_length -r $reverse_primer_length
    #python3 figaro.py -i $cutadapt_dir -o $figaro_dir -a $amplicon_length -f $forward_primer_length -r $reverse_primer_length

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
