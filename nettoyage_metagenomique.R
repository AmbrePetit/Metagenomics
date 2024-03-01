#!/usr/bin/env Rscript
# https://rdrr.io/bioc/dada2/man

##################
### libraries
##################

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(ggplot2))
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
    make_option(c("-f","--fastq"), dest="fastq_dir", default="None",
        help="path to fastq files directory with primer removed [MANDATORY]"),
    make_option(c("-t","--type"), dest="type_data", default="18S",
                help="type of data that you analyse (16S or 18S)"),
    make_option(c("-o", "--output"), dest="output", default="None",
        help="path to output directory [MANDATORY]"),
    make_option("--type_database", dest="type_database",
                help="type of database used for taxonomic assignment"),
    make_option(c("-j", "--json"), dest="json_file", default="None",
                help="path to json file directory, output from Figaro"),
    make_option("--pattern_R1", dest="opt_fns_R1", default="_R1",
        help="Pattern in fastq files name to identify R1"),
    make_option("--pattern_R2", dest="opt_fns_R2", default="_R2",
        help="Pattern in fastq files name to identify R2"),
    make_option("--pattern_samples", dest="pattern_samples", default="_R"),
    make_option(c("-p", "--primer"), dest="primer", help="have the primers already been removed ? 'true' or 'false' expected")
 )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# argument verification
if(opt$output == "None"){print("--output is a mandatory option; exiting now");q("no", 1,FALSE)}
if(opt$fastq_dir == "None"){print("--fastq is a mandatory option; exiting now");q("no", 1,FALSE)}
if(opt$type_data!="16S" && opt$type_data!="18S"){print("--type is 16S or 18S only; exiting now");q("no",1,FALSE)}
if (opt$type_data=="16S") {
  if(opt$json_file == "None"){print("--json is a mandatory option; exiting now");q("no", 1,FALSE)}}

# assign all output dir variables
# 1st check if path ends with a "/"
if(substrRight(opt$output, 1) != '/'){opt$output <- paste(opt$output, "/", sep="")}
if(substrRight(opt$json_file, 1) != '/'){opt$json_file <- paste(opt$json_file, "/", sep="")}

# check output
print("opt$output")
print(opt$output)

# creation results files
filtered_dir <- paste(opt$output, "filtered/", sep="")
dada2_dir <- paste(opt$output, "dada2/", sep="")
dir.create(filtered_dir)
dir.create(dada2_dir)

# databases
database_dir="/Users/ambre/Desktop/Stage/database/"
if (opt$type_data=="16S") {
Silva_file <- paste0(database_dir, "silva_nr99_v138.1_train_set.fa")
}else{
  if (opt$type_data=="18S") { 
    if(opt$type_database == "PR2"){
      PR2_file <- paste0(database_dir, "pr2_version_5.0.0_SSU_dada2.fasta.gz")
    }else{
      Silva_file <- paste0(database_dir, "silva_nr_v138_train_set.fa")}
  }
}

##################
### main
##################

print("-----------Start--------")

# fastq examination

# make a list of all the fastq files in a directory and separate forward and reverse
# if primer not be removed before analysis, fastq files to take are the outputs of cutadapt
fns <- sort(list.files(opt$fastq_dir, full.names = TRUE))
fns <- fns[str_detect(basename(fns), ".fastq")]
print("fns")
print(fns)
fns_R1 <- fns[str_detect(basename(fns), opt$opt_fns_R1)]
fns_R2 <- fns[str_detect(basename(fns), opt$opt_fns_R2)]

# examination of the names
print("fns_R1")
print(fns_R1)
print("fns_R2")
print(fns_R2)

# extract samples names (format : NAMESAMPLE_XXX.fastq)
# pattern in files name to separate samples
sample.names <- str_split(basename(fns_R1), pattern = opt$pattern_samples, simplify = TRUE)
print("sample.names")
print(sample.names)
sample.names <- sample.names[,1]

# add samples names
print("--------------Add sample names--------")
filt_R1 <- str_c(filtered_dir, sample.names, "_R1_filt.fastq")
filt_R2 <- str_c(filtered_dir, sample.names, "_R2_filt.fastq")
names(filt_R1) <- sample.names
names(filt_R2) <- sample.names

print("------------Filter and trim------------")

truncQ <- 2  #trunc sequences if phred quality below threshold
maxEE <- c(2,5) #sequences removed if score below threshold (mean of quality score per base)

if (opt$type_data=="16S") {
  # Load JSON file and extract the first parameter of trimPosition
  json_file <- paste0(opt$json_file, "trimParameters.json")
  json_file <- jsonlite::fromJSON(json_file)
  trim_position <- json_file$trimPosition[1]
  trim_first_parameters <- trim_position[[1]][1]
  trim_second_parameters <- trim_position[[1]][2]
  print("trim_parameters : ")
  print(trim_first_parameters)
  print(trim_second_parameters)
  out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen=c(trim_first_parameters,trim_second_parameters), maxN=0, maxEE=maxEE, truncQ=2, rm.phix=TRUE, compress=FALSE)
  best_truncLen <- c(trim_first_parameters,trim_second_parameters)
  print("------------filterAndTrim function for 16S gene------------")
  print("out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen=trim_parameters, maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE, compress=FALSE)")
}else{
  truncLens <- c(290, 260, 240, 220)
  best_truncLen <- NULL
  best_nb_reads <- 0
  best_result <- NULL

  getN <- function(x) sum(getUniques(x))

  for (tl in truncLens){
      print(tl)
      for (offset in c(-20, -40, -60)) {
          print(tl + offset)
          tryCatch({
              out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen=c(tl, tl + offset), maxN=0, maxEE=maxEE, truncQ=2, rm.phix=TRUE, compress=FALSE)
              print(out)
          
              err_R1 <- learnErrors(filt_R1, multithread=TRUE)
              err_R2 <- learnErrors(filt_R2, multithread=TRUE)
          
              # statistic analyse
              dada_R1 <- dada(filt_R1, err = err_R1, multithread = FALSE, pool = FALSE)
              dada_R2 <- dada(filt_R2, err = err_R2, multithread = FALSE, pool = FALSE)
              print("sapply(dada_R1, getN)")
              filtered <- sapply(dada_R1, getN)
              print("filtered")
              print(filtered)
          
              # merge forward and reverse
              mergers <- mergePairs(dada_R1, filt_R1, dada_R2, filt_R2, verbose = TRUE)
              print("sapply(mergers, getN)")
              merged <- sapply(mergers, getN)
              print("merged")
              print(merged)
          
              # check if this truncLen combination has more merged reads than the previous best one
              if (sum(merged) > best_nb_reads) {
                  best_nb_reads <- sum(merged)
                  best_truncLen <- c(tl, tl + offset)
              }
          }, error = function(e) {  #allows to not stop the program if the length of the reads is smaller than the truncLen threshold
              if (grepl("Not all provided files exist", conditionMessage(e))) {
                  cat("Some input samples had no reads pass the filter.", tl, "+", offset, ". Next combination.\n")
              }else{
                  cat("Error encountered while processing the combination", tl, "+", offset, ":", conditionMessage(e), "\n")}
            })
      }
    }
  print("Best truncLen combination :")
  print(best_truncLen)
  out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen=best_truncLen, maxN=0, maxEE=maxEE, truncQ=2, rm.phix=TRUE, compress=FALSE)
}

print("------------Custering------------")

err_R1 <- learnErrors(filt_R1, multithread=TRUE)
err_R2 <- learnErrors(filt_R2, multithread=TRUE)

# statistic analyse
dada_R1 <- dada(filt_R1, err = err_R1, multithread = FALSE, pool = FALSE)
dada_R2 <- dada(filt_R2, err = err_R2, multithread = FALSE, pool = FALSE)

# merge forward and reverse
mergers <- mergePairs(dada_R1, filt_R1, dada_R2, filt_R2, verbose = TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

t_seqtab <- t(seqtab)
table(nchar(getSequences(seqtab)))

# suppression of chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
paste0("% of non chimeras : ", sum(seqtab.nochim)/sum(seqtab) * 100)
paste0("total number of sequences : ", sum(seqtab.nochim))

getN <- function(x) sum(getUniques(x))

# variable inspection
print("class(mergers)")
class(mergers)
print("length(mergers)")
length(mergers)

print("class(dada_R1)")
class(dada_R1)
print("length(dada_R1)")
length(dada_R1)

# creation of csv with reads number in each sample after each step
print("mergers")
print(head(mergers))
print("seqtab")
print(head(seqtab))
print("seqtab.nochim")
print(head(seqtab.nochim))
track <- cbind(out, sapply(dada_R1, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
write_csv(data.frame(track), str_c(dada2_dir, "number_reads.csv"))

# change the name of the sequences to store them in the taxonomy table
seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% rownames_to_column(var = "sequence") %>% 
  rowid_to_column(var = "ASVNumber") %>% mutate(ASVNumber = sprintf("asv%04d", 
                                                                    ASVNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))

df <- seqtab.nochim_trans
seq_out <- Biostrings::DNAStringSet(df$sequence)

names(seq_out) <- df$ASVNumber

# exportation of fasta file with all uniques sequences found in samples
Biostrings::writeXStringSet(seq_out, str_c(dada2_dir, "asv_no_taxo.fasta"), compress = FALSE, 
                            width = 20000)

# assignation of the taxonomy with the right database according to --type
if(opt$type_data == "16S"){
        print("-------Assignation 16S----------")
        db <- "Silva"
        taxa <- assignTaxonomy(seqtab.nochim, refFasta = Silva_file, minBoot = 50, verbose = TRUE, tryRC = TRUE)
}else{
        print("--------Assignation 18S---------")
        if(opt$type_database == "PR2"){
          PR2_tax_levels <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
          db <- "PR2"
          taxa <- assignTaxonomy(seqtab.nochim, refFasta = PR2_file, taxLevels = PR2_tax_levels, minBoot = 50, outputBootstraps = TRUE, verbose = TRUE, tryRC = FALSE)}
        else{
          db <- "Silva"
          taxa <- assignTaxonomy(seqtab.nochim, refFasta = Silva_file, minBoot = 50, verbose = TRUE, tryRC = TRUE)}
        }
print("------End assignation-----")

# creation of database file
if(opt$type_data == "16S"){
  print("----------16S----------")
  taxa <- as.data.frame(taxa)
  seqtab.nochim <- t(seqtab.nochim)
  seqtab.nochim_trans <- cbind(taxa, seqtab.nochim_trans)
  print("-----------End 16S--------")
}else{
  print("----------18S----------")
  if(opt$type_database == "PR2"){
    print("taxa$tax")
    print(head(taxa$tax))
    print("taxa$boot")
    print(head(taxa$boot))
    taxa_tax <- as.data.frame(taxa$tax)
    print("class(taxa$boot)")
    str(taxa$boot)
    class(taxa$boot)
    taxa_boot <- taxa$boot
    taxa_boot <- as.data.frame(taxa$boot) %>% rename_all(~str_c(., "_boot"))
    seqtab.nochim_trans <- taxa_tax %>% bind_cols(taxa_boot) %>% bind_cols(seqtab.nochim_trans)
  }else{
    taxa <- as.data.frame(taxa)
    seqtab.nochim <- t(seqtab.nochim)
    seqtab.nochim_trans <- cbind(taxa, seqtab.nochim_trans)
  }
  print("-----------End 18S--------")
}

# exportation
write_csv(seqtab.nochim_trans, str_c(dada2_dir, "database.csv"))

# creation fasta
print("df <- seqtab.nochim_trans")
df <- seqtab.nochim_trans
print("df")
print(head(df))
print("df[is.na(df)] <- NA")
df[is.na(df)] <- "NA"
print("seq_out <- Biostrings::DNAStringSet(df$sequence)")
seq_out <- Biostrings::DNAStringSet(df$sequence)

print("names(seq_out) ........")
names(seq_out) <- str_c(df$ASVNumber, df$Domain, df$Supergroup, df$Division, df$Subdivision, df$Class, df$Order, df$Family,
                        df$Genus, df$Species, sep = "|")

print("Biostrings::writeXStringSet ..........")
Biostrings::writeXStringSet(seq_out, str_c(dada2_dir, "asv.fasta"), compress = FALSE, width = 20000)

# creation csv for phyloseq
print("Creation csv for phyloseq")
if(opt$type_database == "PR2"){
    i<-9 ; j<-20     #if you're using an old version of PR2 where there are only 8 levels, remove 1 from i and j
}else{i<-6 ; j<-14}  #if there are species, change to : i <- 7 and j<- 15

# creation taxo file
print("creation taxo file")
tax <- seqtab.nochim_trans[,1:i]
tax <- cbind(seqtab.nochim_trans$ASVNumber, tax)
colnames(tax)[1] <- "ASVNumber"
write_csv(as_tibble(tax), file = str_c(dada2_dir, "taxo.csv"))

# creation asv file
print("creation asv file")
asv <- seqtab.nochim_trans[,j:length(seqtab.nochim_trans)]
asv <- cbind(seqtab.nochim_trans$ASVNumber, asv)
colnames(asv)[1] <- "ASVNumber"
write_csv(as_tibble(asv), file = str_c(dada2_dir, "asv.csv"))

# creation output file with parameters
print("creation output file with parameters")
fichier_temp <- paste0(opt$output, "output_parameters.txt")
contenu <- readLines(fichier_temp)
contenu <- c(contenu, "truncLen : ", best_truncLen, paste("maxEE : 2,5", "\n"), paste("Type of data for analysis : ", opt$type_data, "\n"), paste("Database used : ", db, "\n"), paste("Primer were already been removed ? ", opt$primer))
writeLines(contenu, paste0(opt$output, "output_parameters.txt"))


print("-----------End--------")

