This project takes place during my internship at Institut Pasteur Guadeloupe for my master's degree in Bioinformatics. 

The aim of this project is to perform a metabarcoding analysis of Illumina sequencing fastQ data for 16S, 18S and ITS barcodes. It enables taxonomic assignment of sample data and compositional analysis. 

The bash and R programs are used to clean fastQ data (primer removal, end truncation, etc.) and perform taxonomic assignment using tools such as cutadapt, figaro and the R package DADA2.

The python program allows analysis of the abundance tables provided by the dada2 analysis, and dynamically displays various statistics (alpha and beta diversity) and graphs. 
