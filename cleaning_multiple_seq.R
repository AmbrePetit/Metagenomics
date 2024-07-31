#!/usr/bin/env Rscript
# https://urldefense.com/v3/__https://rdrr.io/bioc/dada2/man__;!!JFdNOqOXpB6UZW0!sfdiZIJNYKiEsUmod6uhw5UthyLRYaqI7Z0UfdX9aALgEuLSAWCO01yYmdEZlvygqVD392vWB0PrC1NKHgHE4EK8BQz77AZlLpg$ 

print("----------- Start R program -----------")

##################
### libraries
##################

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(jsonlite))

##################
### parse arguments
##################

# list of options

option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false",
        dest="verbose", help="Print little output"),
    make_option(c("-m", "--multiple_ana"), dest="mult_ana", default="No",
        help="several different sequencing outputs to analyze together? 'Yes' or 'No' expected"),
    make_option(c("-f","--fastq"), dest="fastq_dir", default="None",
        help="path to FASTQ files [MANDATORY]"),
    make_option(c("-t","--type"), dest="type_data",
        help="type of data that you analyse (16S, 18S or ITS) [MANDATORY]"),
    make_option(c("-o", "--output"), dest="output", default="None",
        help="path to output directory [MANDATORY]"),
    make_option("--type_database", dest="type_database", default="Silva",
        help="type of database used for taxonomic assignment"),
    make_option("--pattern_R1", dest="opt_fns_R1", default="_R1",
        help="Pattern in fastq files name to identify R1"),
    make_option("--pattern_R2", dest="opt_fns_R2", default="_R2",
        help="Pattern in fastq files name to identify R2"),
    make_option("--pattern_samples", dest="pattern_samples", default="_R"),
    make_option(c("-p", "--primer"), dest="primer",
        help="have the primers already been removed ? 'Yes' or 'No' expected"),
    make_option("--figaro", dest="figaro_update", default="TRUE",
        help="Have figaro done an output and succesfully worked ? 'Yes' or 'No' expected"),
    make_option(c("-l", "--amplicon_variable"), dest="amplicon_length_variable", default="No",
        help="Are amplicon lengths variable ? 'Yes' or 'No' expected"),
    make_option(c("-a", "--amplicon_length"), dest="amplicon_length",
        help="Expected length of amplicons"),
    make_option("--primer_forward", dest="primer_fwd",
        help="Forward primer used for amplicon analysis"),
    make_option("--primer_reverse", dest="primer_rev",
        help="Reverse primer used for amplicon analysis")
 )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

print("Arguments:")
print(opt)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# argument verification
if(opt$output == "None"){print("--output is a mandatory option; exiting now");q("no", 1,FALSE)}
if(opt$fastq_dir == "None"){print("--fastq is a mandatory option; exiting now");q("no", 1,FALSE)}
if(opt$type_data!="16S" && opt$type_data!="18S" && opt$type_data!="ITS"){print("--type is 16S, 18S or ITS only; exiting now");q("no",1,FALSE)}

# assign all output dir variables
# 1st check if path ends with a "/"
if(substrRight(opt$output, 1) != '/'){opt$output <- paste(opt$output, "/", sep="")}

# empty list to store sequence tables
seqtab_list <- list()

# databases
database_dir="/Users/ambre/Desktop/MASTER/M2/Stage/database/"
if (opt$type_data=="16S"){
  Silva_file <- paste0(database_dir, "silva_nr99_v138.1_train_set.fa")
  db <- "Silva"
}else if (opt$type_data=="ITS"){
  Unite_file <- paste0(database_dir, "sh_general_release_dynamic_25.07.2023.fasta")
  db <- "Unite"
}else if (opt$type_data=="18S"){
    if(opt$type_database == "PR2"){
      PR2_file <- paste0(database_dir, "pr2_version_5.0.0_SSU_dada2.fasta.gz")
      db <- "PR2"
    }else{
      Silva_file <- paste0(database_dir, "silva_nr_v132_train_set.fa")
      db <- "Silva"
    }
}

# creation results files
dada2_dir <- paste(opt$output, "dada2/", sep="")
dir.create(dada2_dir)
# Créer un data frame vide avec les noms des colonnes souhaitées
column_names <- c("sample", "input", "filtered", "denoised", "merged", "tabled", "nonchim")
empty_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(empty_df) <- column_names
write.table(empty_df, str_c(dada2_dir, "number_reads.csv"), sep = ",")


##################
### main
##################

# Get the list of subdirectories within the specified input directory
subdirs <- list.dirs(opt$fastq_dir, recursive = FALSE)
print(subdirs)
k <- 1
var_list <- list()

# Loop through each input directory
for (input_dir in subdirs) {
    print("input_dir")
    print(input_dir)
    subfiles <- list.dirs(input_dir, recursive = FALSE)
    print("subfiles")
    print(subfiles)
    # make a list of all the fastq files in a directory and separate forward and reverse
    # if primer not be removed before analysis, fastq files to take are the outputs of cutadapt
    if (opt$primer=="No" && (opt$type_data=="16S" || opt$type_data=="18S")){
        fastqfile <- paste0(input_dir, "/fastq_no_primer")
        fns <- sort(list.files(fastqfile, full.names = TRUE))
    }else{
        fns <- sort(list.files(subfiles, full.names = TRUE))
    }
    
    fns <- fns[str_detect(basename(fns), ".fastq")]
    fns_R1 <- fns[str_detect(basename(fns), opt$opt_fns_R1)]
    fns_R2 <- fns[str_detect(basename(fns), opt$opt_fns_R2)]

    # examination of the names
    print("fns_R1")
    print(fns_R1)
    print("fns_R2")
    print(fns_R2)

    # extract samples names (format : NAMESAMPLE.fastq)
    # pattern in files name to separate samples
    sample.names <- str_split(basename(fns_R1), pattern = opt$pattern_samples, simplify = TRUE)
    sample.names <- sample.names[,1]
    
    # creation results files
    filtered_dir <- paste(input_dir, "/filtered/", sep="")
    dir.create(filtered_dir)

    # remove primers for ITS gene
    if (opt$type_data=="ITS" && opt$primer=="No") {
        print("--------------Cutting primer for ITS gene--------")
        allOrients <- function(primer) {
            # Create all orientations of the input sequence
            require(Biostrings)
            dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
            orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna), RevComp = Biostrings::reverseComplement(dna))
            return(sapply(orients, toString))  # Convert back to character vector
        }
        FWD.orients <- allOrients(opt$primer_fwd)
        REV.orients <- allOrients(opt$primer_rev)

        fns_R1.filtN <- file.path(input_dir, "filtN", basename(fns_R1)) # Put N-filtered files in filtN/ subdirectory
        fns_R2.filtN <- file.path(input_dir, "filtN", basename(fns_R2))
        filterAndTrim(fns_R1, fns_R1.filtN, fns_R2, fns_R2.filtN, maxN = 0, multithread = TRUE)

        primerHits <- function(primer, fn) {
            # Counts number of reads in which the primer is found
            nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
            return(sum(nhits > 0))
        }
        
        print("Checking the presence of primer")
        FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fns_R1.filtN)
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fns_R2.filtN)
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fns_R1.filtN)
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fns_R2.filtN)
        print(FWD.ForwardReads)
        print(FWD.ReverseReads)
        print(REV.ForwardReads)
        print(REV.ReverseReads)
        
        print("cutadapt")
        cutadapt <- "/usr/local/bin/cutadapt"
        system2(cutadapt, args = "--version")
        path.cut <- file.path(input_dir, "cutadapt")
        if(!dir.exists(path.cut)) dir.create(path.cut)
        fns_R1.cut <- file.path(path.cut, basename(fns_R1))
        fns_R2.cut <- file.path(path.cut, basename(fns_R2))
        FWD.RC <- dada2:::rc(opt$primer_fwd)
        REV.RC <- dada2:::rc(opt$primer_rev)
        # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
        R1.flags <- paste("-g", opt$primer_fwd, "-a", REV.RC)
        # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
        R2.flags <- paste("-G", opt$primer_rev, "-A", FWD.RC)
        # Run Cutadapt
        for(i in seq_along(fns_R1)) {
          system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-o", fns_R1.cut[i], "-p", fns_R2.cut[i], fns_R1.filtN[i], fns_R2.filtN[i]))
        }
        
        print("Checking presence of primer after cutadapt")
        FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fns_R1.cut)
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fns_R2.cut)
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fns_R1.cut)
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fns_R2.cut)
        
        print(FWD.ForwardReads)
        print(FWD.ReverseReads)
        print(REV.ForwardReads)
        print(REV.ReverseReads)
        
    }

    # add samples names
    print("--------------Add sample names--------")
    filt_R1 <- str_c(filtered_dir, sample.names, "_R1_filt.fastq")
    filt_R2 <- str_c(filtered_dir, sample.names, "_R2_filt.fastq")
    names(filt_R1) <- sample.names
    names(filt_R2) <- sample.names
    print("filt_R1")
    print(filt_R1)
    print("filt_R2")
    print(filt_R2)

    print("------------Filter and trim------------")

    truncQ <- 2  #trunc sequences if phred quality below threshold

    if (opt$figaro_update == "TRUE" && (opt$amplicon_length_variable == "No" && (opt$type_data=="16S" || opt$type_data=="18S"))) {
      # Load JSON file and extract the first parameter of trimPosition
      json_file <- sort(list.files(subfiles, full.names = TRUE))
      json_file <- json_file[str_detect(basename(json_file), "trimParameters.json")]
      json_file <- jsonlite::fromJSON(json_file)
      
      found_params <- FALSE
      i <- 1

      while (!found_params && i <= length(json_file$trimPosition)) {
        trim_position <- json_file$trimPosition[[i]]
        maxEE <- json_file$maxExpectedError[[i]]
        
        trim_first_parameters <- trim_position[1]
        trim_second_parameters <- trim_position[2]
        
        maxEE_first_parameters <- maxEE[1]
        maxEE_second_parameters <- maxEE[2]
        
        if (maxEE_first_parameters <= 8 && maxEE_second_parameters <= 8) {
          found_params <- TRUE
        } else {
          i <- i + 1
        }
      }

      if (found_params) {
        print("trim_parameters : ")
        print(trim_first_parameters)
        print(trim_second_parameters)

        print("maxEE : ")
        print(maxEE_first_parameters)
        print(maxEE_second_parameters)
      } else {  # if all maxEE parameters are greater than 8, take the first parameters for trim position and (5,5) for maxEE
          trim_position <- json_file$trimPosition[1]
          trim_first_parameters <- trim_position[[1]][1]
          trim_second_parameters <- trim_position[[1]][2]
          maxEE_first_parameters <- 5
          maxEE_second_parameters <- 5
      }

      #remove primer length to figaro parameters
      print("figaro parameters with primers")
      print(trim_first_parameters)
      print(trim_second_parameters)
      trim_first_parameters <- trim_first_parameters - nchar(opt$primer_fwd)
      trim_second_parameters <- trim_second_parameters - nchar(opt$primer_rev)
      print("figaro parameters with no primers")
      print(trim_first_parameters)
      print(trim_second_parameters)


      out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen=c(trim_first_parameters,trim_second_parameters), maxN=0, maxEE=c(maxEE_first_parameters, maxEE_second_parameters), rm.phix=TRUE, compress=FALSE)
      print("------------filterAndTrim function------------")
      cat(paste0("out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen=c(",trim_first_parameters,",",trim_second_parameters,"), maxN=0, maxEE=c(",maxEE_first_parameters,",",maxEE_second_parameters,"), rm.phix=TRUE, compress=FALSE)"))
    }else if (opt$type_data=="ITS"){
      out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, maxN=0, maxEE=c(2,5), minLen = 50, rm.phix=TRUE, compress=FALSE)
      print("------------filterAndTrim function------------")
      print("out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, maxN=0, maxEE=c(2,5), minLen = 50, rm.phix=TRUE, compress=FALSE)")
    }else if (opt$figaro_update == "FALSE" || (opt$amplicon_length_variable == "Yes" && (opt$type_data=="16S" || opt$type_data=="18S"))){
      truncLens <- c(290, 260, 240, 220)
      maxEE_first_parameters <- 2
      maxEE_second_parameters <- 5
      maxEE <- c(maxEE_first_parameters, maxEE_second_parameters) #sequences removed if score below threshold (mean of quality score per base)
      best_truncLen <- NULL
      best_nb_reads <- 0
      best_result <- NULL

      getN <- function(x) sum(getUniques(x))

      for (tl in truncLens) {
          print(tl)
          for (offset in c(-20, -40, -60)) {
              print(tl + offset)
              tryCatch({
                  out_test <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen=c(tl, tl + offset), maxN=0, maxEE=maxEE, rm.phix=TRUE, compress=FALSE)
                  print(out_test)
                  
                  err_R1 <- learnErrors(filt_R1, multithread=TRUE)
                  err_R2 <- learnErrors(filt_R2, multithread=TRUE)
              
                  # statistic analyse
                  dada_R1 <- dada(filt_R1, err = err_R1, multithread = FALSE, pool = FALSE)
                  dada_R2 <- dada(filt_R2, err = err_R2, multithread = FALSE, pool = FALSE)
                  filtered <- sapply(dada_R1, getN)
                  print("filtered")
                  print(filtered)
                  
                  # merge forward and reverse
                  mergers <- mergePairs(dada_R1, filt_R1, dada_R2, filt_R2, verbose = TRUE)
                  merged <- sapply(mergers, getN)
                  print("merged")
                  print(merged)
                  
                  # check if this truncLen combination has more merged reads than the previous best one
                  if (sum(merged) > best_nb_reads) {
                      best_nb_reads <- sum(merged)
                      best_truncLen <- c(tl, tl + offset)
                      trim_first_parameters <- tl
                      trim_second_parameters <- tl + offset
                      out <- out_test
                      
                  }
              }, error = function(e) { #allows to not stop the program if the length of the reads is smaller than the truncLen threshold
                  if (grepl("Not all provided files exist", conditionMessage(e))) {
                      cat("Some input samples had no reads pass the filter.", tl, "+", offset, ". Next combination.\n")
                  }else{
                      cat("Error encountered while processing the combination", tl, "+", offset, ":", conditionMessage(e), "\n")}
                })
          }
    }
      print("Best truncLen combination :")
      print(best_truncLen)
      print("------------filterAndTrim function------------")
      print("out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen=c(",trim_first_parameters,",",trim_second_parameters,"), maxN=0, maxEE=c(",maxEE_first_parameters,",",maxEE_second_parameters,"), rm.phix=TRUE, compress=FALSE)")
    }

    if (opt$type_data=="16S" || opt$type_data=="18S"){

        allOrients <- function(primer) {
        # Create all orientations of the input sequence
        require(Biostrings)
        dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
        orients <- c(Forward = dna, Reverse = Biostrings::reverse(dna))
        return(sapply(orients, toString))  # Convert back to character vector
        }
        
        FWD.orients <- allOrients(opt$primer_fwd)
        REV.orients <- allOrients(opt$primer_rev)

        primerHits <- function(primer, fn) {
            # Counts number of reads in which the primer is found
            nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
            return(sum(nhits > 0))
        }

        print("Checking the presence of primer after filtering")
        FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filt_R1)
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = filt_R2)
        print(FWD.ForwardReads)
        print(REV.ReverseReads)
        
    }

    print("------------Custering------------")

    err_R1 <- learnErrors(filt_R1, multithread=TRUE)
    err_R2 <- learnErrors(filt_R2, multithread=TRUE)

    # statistic analyse
    dada_R1 <- dada(filt_R1, err = err_R1, multithread = FALSE, pool = FALSE)
    dada_R2 <- dada(filt_R2, err = err_R2, multithread = FALSE, pool = FALSE)

    # merge forward and reverse
    mergers <- mergePairs(dada_R1, filt_R1, dada_R2, filt_R2, verbose = TRUE)
    seqtab <- makeSequenceTable(mergers)
    
    print("Number of sequences by lengths in seqtab")
    print(table(nchar(getSequences(seqtab))))

    # remove sequences that are longer or shorter than expected amplicons lengths with 10% error
    if (opt$amplicon_length_variable == "No" && (opt$type_data=="16S" || opt$type_data=="18S")){
        lower_bound <- as.numeric(opt$amplicon_length) * 0.9
        upper_bound <- as.numeric(opt$amplicon_length) * 1.1
        seqtab <- seqtab[, nchar(colnames(seqtab)) %in% lower_bound:upper_bound]
        print("Number of sequences by lengths in seqtab after removing longer or shorter sequences than expected")
        print(table(nchar(getSequences(seqtab))))
    }
    
    # suppression of chimera
    print("suppression of chimera")
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
    print("% of non chimeras : ")
    print(sum(seqtab.nochim)/sum(seqtab) * 100)
    print("total number of sequences : ")
    print(sum(seqtab.nochim))
    
    # creation of csv with reads number in each sample after each step
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dada_R1, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
    rownames(track) <- sample.names
    print(track)
    write.table(track, str_c(dada2_dir, "number_reads.csv"), append = TRUE, sep = ",", col.names = FALSE)
    
    # creation output file with parameters
    print("output file with parameters")
    fichier_temp <- paste(opt$output, "output_parameters.txt", sep="")
    contenu <- readLines(fichier_temp)
    if(opt$type_data == "ITS"){
        contenu <- c(contenu, paste("\nInformations for sequencing output ", k, " : \n"), paste("maxEE forward : 2 \n maxEE reverse : 5 \n "), paste("minLen : 50 \n"))
    }else{
        contenu <- c(contenu, paste("\nInformations for sequencing output ", k, " : \n"), paste(" truncLen forward : ", trim_first_parameters, "\n truncLen reverse : ", trim_second_parameters), paste("\n maxEE forward : ", maxEE_first_parameters, "\n maxEE reverse : ", maxEE_second_parameters))
    }
    writeLines(contenu, paste0(opt$output, "output_parameters.txt"))
    
    
    if (opt$mult_ana == "Yes"){
        # Save seqtab in a format which keep all informations
        saveRDS(seqtab.nochim, file = paste0(k, "_seqtab.rds"))
        k <- k+1
    }
    
}

if (opt$mult_ana == "Yes"){
    for (l in 1:(k-1)){
        var_list[[l]] <- readRDS(file = paste0(l, "_seqtab.rds"))
    }
    
    print("var_list")
    for (i in seq_along(var_list)) {
        # Extract current seqtab
        seqtab_varlist <- var_list[[i]]
    
        # Check if seqtab contains lines' names
        sample_names <- rownames(seqtab_varlist)
            
        # Display name of the first sample
        print(paste("Seqtab", i, ":", sample_names))
        }
  
    # Once the analysis for each input directory is complete, merge the sequence tables together
    seqtab.nochim <- do.call(mergeSequenceTables, var_list)
    print("merged seqtab")
    sample_names <- rownames(seqtab.nochim) # Display name of samples
    print(sample_names) # Display only the first elements
}

# Change the name of the sequences to store them in the taxonomy table
seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% rownames_to_column(var = "sequence") %>%
  rowid_to_column(var = "ASVNumber") %>% mutate(ASVNumber = sprintf("asv%04d", ASVNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))

df <- seqtab.nochim_trans
seq_out <- Biostrings::DNAStringSet(df$sequence)

names(seq_out) <- df$ASVNumber

# Exportation of fasta file with all uniques sequences found in samples
Biostrings::writeXStringSet(seq_out, str_c(dada2_dir, "asv_no_taxo.fasta"), compress = FALSE, width = 20000)

# Assignation of the taxonomy with the right database according to --type
if(opt$type_data == "16S"){
        print("-------Assignation 16S----------")
        taxa <- assignTaxonomy(seqtab.nochim, refFasta = Silva_file, minBoot = 50, outputBootstraps = TRUE, verbose = TRUE, tryRC = TRUE)
}else if(opt$type_data == "ITS"){
        print("-------Assignation ITS----------")
        taxa <- assignTaxonomy(seqtab.nochim, refFasta = Unite_file, minBoot = 50, outputBootstraps = TRUE, verbose = TRUE, tryRC = TRUE)
}else if(opt$type_data == "18S"){
        print("--------Assignation 18S---------")
        if(opt$type_database == "PR2"){
          PR2_tax_levels <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
          taxa <- assignTaxonomy(seqtab.nochim, refFasta = PR2_file, taxLevels = PR2_tax_levels, minBoot = 50, outputBootstraps = TRUE, verbose = TRUE, tryRC = TRUE)}
        else{
          taxa <- assignTaxonomy(seqtab.nochim, refFasta = Silva_file, minBoot = 50, outputBootstraps = TRUE, verbose = TRUE, tryRC = TRUE)}
        }
print("------End assignation-----")

# Creation of database file
taxa_tax <- as.data.frame(taxa$tax)
taxa_boot <- taxa$boot
taxa_boot <- as.data.frame(taxa$boot) %>% rename_all(~str_c(., "_boot"))
seqtab.nochim_trans <- taxa_tax %>% bind_cols(taxa_boot) %>% bind_cols(seqtab.nochim_trans)

# Exportation
write_csv(seqtab.nochim_trans, str_c(dada2_dir, "database.csv"))

# Creation fasta
print("df <- seqtab.nochim_trans")
df <- seqtab.nochim_trans
print("df[is.na(df)] <- NA")
df[is.na(df)] <- "NA"
print("seq_out <- Biostrings::DNAStringSet(df$sequence)")
seq_out <- Biostrings::DNAStringSet(df$sequence)

print("names(seq_out) ........")
names(seq_out) <- str_c(df$ASVNumber, df$Domain, df$Supergroup, df$Division, df$Subdivision, df$Class, df$Order, df$Family, df$Genus, df$Species, sep = "|" )

print("Biostrings::writeXStringSet ..........")
Biostrings::writeXStringSet(seq_out, str_c(dada2_dir, "asv.fasta"), compress = FALSE, width = 20000)

# Creation csv for analysis
print("Creation csv")
if(opt$type_database == "PR2"){
    i<-9 ; j<-21     #if you're using an old version of PR2 where there are only 8 levels, remove 1 from i and j
}else if(opt$type_database == "Silva"){
    i<-6 ; j<-15
}else if(opt$type_database == "Unite"){
    i<-7 ; j<-17
}

# Creation taxo file
print("creation taxo file")
tax <- seqtab.nochim_trans[,1:i]
tax <- cbind(seqtab.nochim_trans$ASVNumber, tax)
colnames(tax)[1] <- "ASVNumber"
write_csv(as_tibble(tax), file = str_c(dada2_dir, "taxo.csv"))

# Creation asv file
print("creation asv file")
asv <- seqtab.nochim_trans[,j:length(seqtab.nochim_trans)]
asv <- cbind(seqtab.nochim_trans$ASVNumber, asv)
colnames(asv)[1] <- "ASVNumber"
write_csv(as_tibble(asv), file = str_c(dada2_dir, "asv.csv"))



print("-----------End--------")

