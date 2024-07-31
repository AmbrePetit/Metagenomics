#!/bin/bash


fasta_file=$1 #path to asv.fasta file, output from dada2
output_dir=$2 #path to save alignment and newick files

FastTree_dir="/Users/ambre/Downloads"

mafft --auto $fasta_file > $output_dir/alignment.fasta
$FastTree_dir/FastTree -nt $output_dir/alignment.fasta > $output_dir/phyloTree.newick
