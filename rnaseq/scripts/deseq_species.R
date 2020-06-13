#deseq_species.R 

#SET UP THE DATA TO RUN DESEQ
library(DESeq2)
library(cowplot)
library(tidyverse)
rm(list=ls())

source('rnaseq/scripts/rnaseq_functions.R')

# LOAD THE SUBSET DO ANALYZE  -----------------------------------------------------------
#(see rnaseq/sex_compare_groups for these)
source('rnaseq/species_compare_groups/Ppicta_v_Gholbrooki.R')


# RUN DESEQ ---------------------------------------------------------------


run_deseq_species = function(counts, coldata, PREFIX, CONTRAST){
  coldata$Organism = factor(coldata$Organism)
  coldata$sex = factor(coldata$sex)
  #get variance stabilized counts for timepoint 2
  dds<-DESeqDataSetFromMatrix(counts,
                              colData = coldata, 
                              design = formula(~ Organism)
  )
  
  INDEPDENDENT_FILTERING = FALSE
  dds <- DESeq(dds)
  resultsNames(dds)
  res = results(dds,
                contrast = CONTRAST,
                independentFiltering=INDEPDENDENT_FILTERING)
  rld = rlog(dds)
  
  outname = paste('rnaseq/results_files/', PREFIX, '_deseqResults.Rdata', sep='')
  save(res,coldata, file=outname)
  rldOut = paste('rnaseq/results_files/', PREFIX, '_rld.Rdata', sep='')
  save(rld,coldata, file=rldOut)
  print('--------------')
  print(paste('Run for', PREFIX, '...', sep=''))
  print(paste('rld saves as', rldOut))
  print(paste('deseq results saved as', outname))
}

#RUN FOR FEMALES
run_deseq_species(fcounts,
                  fcoldata,
                  PREFIX="Ppicta_v_Gholbrooki_females",
                  CONTRAST = c("Organism", "Micropoecilia_picta", "Gambusia_holbrooki"))

#RUN FOR MALES
run_deseq_species(mcounts,
                  mcoldata,
                  PREFIX="Ppicta_v_Gholbrooki_males",
                  CONTRAST = c("Organism", "Micropoecilia_picta", "Gambusia_holbrooki"))



