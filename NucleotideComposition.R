# ###################################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2023
# Email:  nitishimtech@gmail.com
# Date: 14 March, 2023
# Script Name: code_name.R
#####################################################################
# ##################  Script Description: ###########################
# R code for Nucleotide copmposition
##################### Notes: ########################################
# This function will work for both DNA and RNA
#####################################################################
################## SET WORKING DIRECTORY ############################
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "/home/nmishra/"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
#####################################################################
#####################################################################
#suppressMessages(suppressWarnings(library(rDNAse)))

# Function to convert DNA in RNA
DNA_to_RNA <- function(x){
  x <- gsub("T", "U", x)
  return(x)
}

# Function to convert RNA in DNA

RNA_to_DNA <- function(x){
  x <- gsub("U", "T", x)
  return(x)
}

# Function for making k-mers index
index_kmer <- function (k, alphabet = "ACGT") 
{
  dict = unlist(strsplit(alphabet, ""))
  make_index = list()
  temp = rep(dict, k)
  make_index = sort(unique(apply(expand.grid(split(temp, rep(1:k, 
                                                             each = 4))), 1, paste, collapse = "")))
  return(make_index)
}

NucleotideComposition <- function (x, k = 2, upto = FALSE, normalize = FALSE, reverse = FALSE, Nucleotide="DNA")
{
  if(Nucleotide=="DNA")
  {
    dict = c("A", "C", "G", "T")
    make_index = list()
    kmer = c()
    if (upto) 
      l = 1
    else l = k
    for (j in l:k) {
      #make_index = make_kmer_index(j, alphabet = "ACGT")
      make_index <- index_kmer(k, alphabet = "ACGT") 
      xPaste = c()
      n = nchar(x)
      xSplit = strsplit(x, split = "")[[1]]
      for (i in 1:j) {
        temp = xSplit[i:(n - j + i)]
        xPaste = paste0(xPaste, temp)
      }
      temp_kmer = summary(factor(xPaste, levels = make_index), 
                          maxsum = 4^j)
      if (reverse) {
        reverse_index = sapply(make_index, revchars)
        reverse_index = chartr("ACGT", "TGCA", reverse_index)
        temp_kmer_reverse = summary(factor(xPaste, levels = reverse_index), 
                                    maxsum = 4^j)
        cmp = which((make_index > reverse_index) - (make_index < 
                                                      reverse_index) < 0)
        temp_kmer[cmp] = temp_kmer[cmp] + temp_kmer_reverse[cmp]
        cmp = which((make_index > reverse_index) - (make_index < 
                                                      reverse_index) <= 0)
        temp_kmer = temp_kmer[cmp]
      }
      if (normalize) 
        temp_kmer = temp_kmer/sum(temp_kmer)
      kmer = append(kmer, temp_kmer)
    }
    return(kmer)
  }
  if(Nucleotide=="RNA")
  {
    dict = c("A", "C", "G", "U")
    make_index = list()
    kmer = c()
    if (upto) 
      l = 1
    else l = k
    for (j in l:k) {
      make_index = index_kmer(j, alphabet = "ACGU")
      xPaste = c()
      n = nchar(x)
      xSplit = strsplit(x, split = "")[[1]]
      for (i in 1:j) {
        temp = xSplit[i:(n - j + i)]
        xPaste = paste0(xPaste, temp)
      }
      temp_kmer = summary(factor(xPaste, levels = make_index), 
                          maxsum = 4^j)
      if (reverse) {
        reverse_index = sapply(make_index, revchars)
        reverse_index = chartr("ACGU", "UGCA", reverse_index)
        temp_kmer_reverse = summary(factor(xPaste, levels = reverse_index), 
                                    maxsum = 4^j)
        cmp = which((make_index > reverse_index) - (make_index < 
                                                      reverse_index) < 0)
        temp_kmer[cmp] = temp_kmer[cmp] + temp_kmer_reverse[cmp]
        cmp = which((make_index > reverse_index) - (make_index < 
                                                      reverse_index) <= 0)
        temp_kmer = temp_kmer[cmp]
      }
      if (normalize) 
        temp_kmer = temp_kmer/sum(temp_kmer)
      kmer = append(kmer, temp_kmer)
    }
    return(kmer)
  }
}




################################################################################
################################################################################
## ******** Binary composition ******** ##

BinaryComposition <- function (x, reverse = FALSE, Nucleotide="DNA"){
  if (reverse==TRUE) {
    x <- Biostrings::reverse(x)
  }
  if(Nucleotide=="DNA")
  {
    Y <- mgsub::mgsub(x, pattern = c("A", "C", "G", "T"), replacement = c("1000", "0100", "0010","0001"))
    Y <- as.numeric(unlist(strsplit(Y, "")))
    return(Y)
  }
  if(Nucleotide=="RNA"){ 
    Y <- mgsub::mgsub(x, pattern = c("A", "C", "G", "U"), replacement = c("1000", "0100", "0010","0001"))
    Y <- as.numeric(unlist(strsplit(Y, "")))
    return(Y)
  }
}

################################################################################
################################################################################
# **************** Running example ****************** #
x <- "AATTCATGCGTCCGGACTTCTGCCTCGAGCCGCCGTACACTGGGCCCTGCAAAGCTC"
NucleotideComposition(x, reverse = FALSE, Nucleotide = "DNA")
NucleotideComposition(x, reverse = FALSE, k = 1, Nucleotide = "DNA")

X <- DNA_to_RNA(x)
NucleotideComposition(X, reverse = FALSE, Nucleotide = "RNA")
NucleotideComposition(X, reverse = FALSE, Nucleotide = "RNA", k = 1)

BinaryComposition(x, Nucleotide = "DNA", reverse = FALSE)
BinaryComposition(x, Nucleotide = "DNA", reverse = TRUE)
BinaryComposition(DNA_to_RNA(X), Nucleotide = "RNA", reverse = FALSE)

