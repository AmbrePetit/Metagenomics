#!/bin/bash

echo "-------------------- Run of data cleaning and taxonomic assignment programs --------------------"
# Define parameters
multiple_analysis=$1 #different sequencing output to analyse together ? "Yes" or "No" expected
base_dir=$2 #path to the main fold containing fastq files
output_dir=$3 #path to output fold
primer_removed=$4 #primer already removed? "Yes" or "No" expected
type_data=$5 #16S, 18S or ITS
type_database=$6 #silva, pr2 or unite
forward_primer=$7 #forward primer sequence
reverse_primer=$8 #reverse primer sequence
# Just for 16S or 18S
amplicon_length_variable=$9 #are amplicon lengths variable ? "Yes" or "No" expected
# Just for 16S or 18S non variable
amplicon_length=${10} #expected amplicon length

# Parameters for figaro
figaro="/media/linuxipg/1.42.6-8754/Isaure/figaro/figaro"
forward_primer_length=${#forward_primer}
reverse_primer_length=${#reverse_primer}

# Parameters for R script to clean data with DADA2
rscript="/usr/bin"
script_dada2_dir="/media/linuxipg/1.42.6-8754/Ambre/Script"

# Add informations about the run in the file
echo "Different sequencing outputs to analyze together ? $multiple_analysis " > "$output_dir/output_parameters.txt"
echo "Type of data for analysis : $type_data" >> "$output_dir/output_parameters.txt"
echo "Database used : $type_database" >> "$output_dir/output_parameters.txt"
echo "Primers sequence used : " >> "$output_dir/output_parameters.txt"
echo "   Forward primer = $forward_primer" >> "$output_dir/output_parameters.txt"
echo "   Reverse primer = $reverse_primer" >> "$output_dir/output_parameters.txt"
echo "Primer were already been removed ? $primer_removed" >> "$output_dir/output_parameters.txt"
echo "Version of tools used during this analysis : " >> "$output_dir/output_parameters.txt"
echo "   Cutadapt : $(cutadapt --version)" >> "$output_dir/output_parameters.txt"
echo "   R : $($rscript/Rscript --version)" >> "$output_dir/output_parameters.txt"
echo "----------- Start Bash program -----------" 



#### CUTADAPT ####
# The aim is to cut all sequences to the same length
for sub_dir in "$base_dir"/*/
do
    fastq_dir=$(readlink -f "$sub_dir/fastq")
    
    echo $sub_dir
    echo $fastq_dir

    ## Quality control ##
    for file in $fastq_dir/*.fastq.gz;
    do
        [ -e "$file" ] || continue
            
        # Check if files contains more than 10000 reads, if not delete file
        number_reads=$(gzcat "$file" | awk '!(NR % 4) {k++}END{print k}') #number of reads for the file
        echo "Number of reads in the file $file : $number_reads " >> "$output_dir/output_parameters.txt"
            
        if [[ "$number_reads" -lt 10000 ]]; then
            echo "Number of reads is less than 10,000. Deleting file: $file" >> "$output_dir/output_parameters.txt"
            rm "$file"
            continue
        fi
        
        length=$(gzcat "$file" | awk 'NR%4 == 2 {print length($0)}') #length of sequences for the file
        sum_length=$(gzcat "$file" | awk 'NR%4 == 2 {sum+=length($0)} END {print sum}') #sum of sequence lengths for the file
        number_seq=$(gzcat "$file" | awk '!(NR % 4) {k++}END{print k}') #number of sequences for the file
        mean=$((sum_length/number_seq)) #mean of length for the file
            
        total_sum_length=$((total_sum_length + sum_length))
        total_number_seq=$((total_number_seq + number_seq))
        
    done
    
    total_mean=$((total_sum_length / total_number_seq))
    echo "Total average length of sequences in all files : $total_mean" >> "$output_dir/output_parameters.txt"
    
    if [[ ($type_data == "16S" || $type_data == "18S") && $amplicon_length_variable == "No" ]]; then
        min_length_amplicon_expected=$((amplicon_length / 2 - forward_primer_length - reverse_primer_length))
        if [[ total_mean -lt $min_length_amplicon_expected ]]; then
            echo "Error : Total average length of sequences in all files is less than $min_length_amplicon_expected bp (minimum length amplicon expected), exit program."
            exit 1
        fi
    elif [[ ($type_data == "ITS") || (($type_data == "16S" || $type_data == "18S") && $amplicon_length_variable == "Yes") ]]; then
        if [[ total_mean -lt 100 ]]; then
            echo "Error : Total average length of sequences in all files is less than 100 bp, exit program."
            exit 1
        fi
    fi

    if [[ ($type_data == "16S" || $type_data == "18S") && $amplicon_length_variable == "No" ]]
    then
        echo "---------Type of gene studied : $type_data"
        
        cutadapt_dir="$sub_dir/cutadapt"
        figaro_dir="$sub_dir/figaro"
        
        # Parameters for cutadapt mean
        total_sum_length=0
        total_number_seq=0
        
        figaro_update="TRUE"
        
        # Check if directory exists, if not create it
        if [ ! -d "$cutadapt_dir" ]; then
            mkdir -p "$cutadapt_dir"
        fi

        if [ ! -d "$figaro_dir" ]; then
            mkdir -p "$figaro_dir"
        fi
        

        echo "Length at which the sequences was cut to run figaro : $total_mean" >> "$output_dir/output_parameters.txt"

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
        fastq_no_primer="$sub_dir/fastq_no_primer"
        
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

            cutadapt -g $forward_primer -G $reverse_primer -o "$output_R1" -p "$output_R2" "$file_R1" "$file_R2"
        done
        echo "Primers were removed with cutadapt for files in $fastq_dir. " >> "$output_dir/output_parameters.txt"
    fi
done

#### R SCRIPT DADA2 ####

if [[ $multiple_analysis == "No" ]]
then
    if [[ $type_data == "16S" || $type_data == "18S" ]]
    then
        if [[ $amplicon_length_variable == "No" ]]
        then
            echo "Command line used to run R script : $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -f $base_dir -t $type_data -o $output_dir --type_database $type_database -j $figaro_dir -p $primer_removed --figaro $figaro_update -a $amplicon_length --primer_forward $forward_primer --primer_reverse $reverse_primer" >> "$output_dir/output_parameters.txt"
            $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -f $base_dir -t $type_data -o $output_dir --type_database $type_database -j $figaro_dir -p $primer_removed --figaro $figaro_update -a $amplicon_length --primer_forward $forward_primer --primer_reverse $reverse_primer
        
        elif [[ $amplicon_length_variable == "Yes" ]]
        then
            echo "Command line used to run R script : $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -f $base_dir -t $type_data -o $output_dir --type_database $type_database -p $primer_removed -l $amplicon_length_variable --primer_forward $forward_primer --primer_reverse $reverse_primer" >> "$output_dir/output_parameters.txt"
            $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -f $base_dir -t $type_data -o $output_dir --type_database $type_database -p $primer_removed -l $amplicon_length_variable --primer_forward $forward_primer --primer_reverse $reverse_primer
        fi
    elif [[ $type_data == "ITS" ]]
    then
        echo "Command line used to run R script : $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -f $base_dir -t $type_data -o $output_dir --type_database $type_database -p $primer_removed --primer_forward $forward_primer --primer_reverse $reverse_primer " >> "$output_dir/output_parameters.txt"
        $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -f $base_dir -t $type_data -o $output_dir --type_database $type_database -p $primer_removed --primer_forward $forward_primer --primer_reverse $reverse_primer
    fi
        
elif [[ $multiple_analysis == "Yes" ]]
then
    if [[ $type_data == "16S" || $type_data == "18S" ]]
    then
        if [[ $amplicon_length_variable == "No" ]]
        then
            echo "Command line used to run R script : $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -m $multiple_analysis -f $base_dir -t $type_data -o $output_dir --type_database $type_database -p $primer_removed -a $amplicon_length --primer_forward $forward_primer --primer_reverse $reverse_primer" >> "$output_dir/output_parameters.txt"
            $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -m $multiple_analysis -f $base_dir -t $type_data -o $output_dir --type_database $type_database -p $primer_removed -a $amplicon_length --primer_forward $forward_primer --primer_reverse $reverse_primer
            
        elif [[ $amplicon_length_variable == "Yes" ]]
        then
            echo "Command line used to run R script : $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -m $multiple_analysis -f $base_dir -t $type_data -o $output_dir --type_database $type_database -p $primer_removed -l $amplicon_length_variable --primer_forward $forward_primer --primer_reverse $reverse_primer" >> "$output_dir/output_parameters.txt"
            $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -m $multiple_analysis -f $base_dir -t $type_data -o $output_dir --type_database $type_database -p $primer_removed -l $amplicon_length_variable --primer_forward $forward_primer --primer_reverse $reverse_primer
        fi
        
    elif [[ $type_data == "ITS" ]]
    then
        echo "Command line used to run R script : $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -m $multiple_analysis -f $base_dir -t $type_data -o $output_dir --type_database $type_database -p $primer_removed --primer_forward $forward_primer --primer_reverse $reverse_primer" >> "$output_dir/output_parameters.txt"
        $rscript/Rscript $script_dada2_dir/cleaning_multiple_seq.R -m $multiple_analysis -f $base_dir -t $type_data -o $output_dir --type_database $type_database -p $primer_removed --primer_forward $forward_primer --primer_reverse $reverse_primer
    fi
fi
