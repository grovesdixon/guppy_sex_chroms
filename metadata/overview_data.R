#overview_data.R
#get overall impression of the data
#comes in 2 SRA projects

library(tidyverse)
rm(list=ls())

# UPLOAD THE TWO SRA RUN TABLES -------------------------------------------

#PRJNA528814
sra1 = read_tsv('metadata/PRJNA528814_SraRunTable.txt') %>% 
  dplyr::select(BioSample,
         Assay_Type,
         Organism,
         Run,
         Sample_Name,
         individual,
         sex,
         tissue,
         BioProject)

#PRJNA353986
sra2 = read_tsv('metadata/PRJNA353986_SraRunTable.txt') %>% 
  mutate(individual = 'none') %>% 
  dplyr::rename(Assay_Type=`Assay Type`,
         Sample_Name = `Sample Name`) %>% 
  dplyr::select(BioSample,
         Assay_Type,
         Organism,
         Run,
         Sample_Name,
         individual,
         sex,
         tissue,
         BioProject)

#bind them
sra = rbind(sra1, sra2)
sra

#write out
sra %>% 
  write_csv(path='metadata/mySraRunTable.csv')


# GET SOME STATS ----------------------------------------------------------


#overall counts
sum1 = sra %>% 
  group_by(BioProject,
           Assay_Type,
           Organism,
           individual,
           Sample_Name,
           sex) %>% 
  summarize(N=n())
sum1
sum1 %>% 
  write_csv(path='~/Desktop/sampleCounts.csv')


#overall counts another way
sum2 = sra %>% 
  group_by(BioProject,
           Assay_Type,
           Organism,
           sex) %>% 
  summarize(Nrun=n(),
            Nsample = length(unique(BioSample)))
sum2

sum2 %>% 
  mutate(xNchr = Nsample*24)

#where Nrun gives us the number of Runs, and Nsample gives us the number of actual biological samples
#will build commands to concatenate all the runs from the same samples

sum2 %>% 
  group_by(Assay_Type) %>% 
  summarize(tot_Nrun = sum(Nrun),
            tot_Nsample = sum(Nsample))


#check by species
sum2 %>% 
  group_by(Assay_Type, Organism, sex) %>% 
  summarize(tot_Nrun = sum(Nrun),
            tot_Nsample = sum(Nsample))


# PRINT OUT COMMANDS FOR CONCATENATING MULTIPLE FASTQS --------------------
#based on sum1, we have many duplicate runs for individual samples
#looking at sum2, the number of samples for all species except reticulata is 3 male and 3 female,
#which matches with what is given in table S1 of the manuscript
#will use the sraRun tables to print out catting commands to combine the trim files from the multi-run samples

#get samples and individuals
get_samples = function(sra, organism){
  sub = sra %>% 
    filter(Organism==organism)
  sub.samples = sub$Sample_Name
  sub.inds = sub$individual
  length(unique(sub.inds)) #(so sample splits by assay type, individual does not)
  print(paste('Total unique individuals =', length(unique(sub.inds))))
  print(paste('Total unique samples =', length(unique(sub.samples))))
  print('returning samples')
  return(unique(sub.samples))
}

#print catting commands to concatenate the bisulfite *.trim files
build_cat_commands = function(sra, reps){
  spp.commands = c()
  totRuns = sra %>% 
    filter(Sample_Name %in% reps) %>% 
    nrow()
  for (r in reps){
    sub=sra[sra$Sample_Name==r,]
    assayType=as.character(unique(sub$Assay_Type))
    if (length(assayType) > 1){
      print('ERROR!')
      break
    }
    at=assayType[1]
    runs = sub$Run
    fruns=paste(runs, "1.trim", sep="_")
    rruns=paste(runs, "2.trim", sep="_")
    forString = paste(fruns, collapse=" ")
    revString = paste(rruns, collapse=" ")
    forOut = paste(paste(">", r), "_1_", at, "_catted.fq", sep='')
    revOut = paste(paste(">", r), "_2_", at, "_catted.fq", sep='')
    commandFor = paste( paste("cat", forString), forOut)
    commandRev = paste( paste("cat", revString), revOut)
    print("---------")
    print(r)
    print(paste(r, "runs:"))
    print(runs)
    print("commands:")
    print(commandFor)
    print(commandRev)
    spp.commands = append(spp.commands, commandFor)
    spp.commands = append(spp.commands, commandRev)
  }
  print('========')
  print(paste('Total runs =', totRuns))
  print(paste('Total samples =', length(reps)))
  print(paste('Total cat commands made (including paired ends) =', length(spp.commands)))
  return(spp.commands)
}

#----- Poecilia wingei
#DNA = samples 291-296; expect 7 runs per sample
#RNA = samples 201-203 + 265-267; expect 4 runs per sample
#12 total samples
#6*7+6*4 == 66 total runs
p.w.samples = get_samples(sra, 'Poecilia wingei')
pw.commands = build_cat_commands(sra, reps=unique(p.w.samples))
length(pw.commands)

#----- Poecilia picta
#DNA = samples 247,248,265-268; expect 22 runs per sample
#RNA = samples 013-019; expect 4 runs per sample
#12 total samples
#6*22+6*4 == 156 total runs
p.p.samples = get_samples(sra, 'Micropoecilia picta')#note genus change compared to table S1
pp.commands = build_cat_commands(sra, reps=unique(p.p.samples))
length(pp.commands)

#----- Poecilia latipinna
#DNA = samples 269-272 + 289-290; expect 22 runs per sample
#RNA = samples 002,004,005,006,007,012; expect 4 runs per sample
#12 total samples
#6*22+6*4 == 156 total runs
p.l.samples = get_samples(sra, 'Poecilia latipinna')
pl.commands = build_cat_commands(sra, reps=unique(p.l.samples))
length(pl.commands)

#----- Gambusia holbrooki
#DNA = samples 241-246; expect 22 runs per sample
#RNA = samples 204-208 + 281; expect 4 runs per sample
#12 total samples
#6*22+6*4 == 156 total runs
g.h.samples = get_samples(sra, 'Gambusia holbrooki')
gh.commands = build_cat_commands(sra, reps=unique(g.h.samples))
length(gh.commands)


#----- Poecilia reticulata
#!!NOTE, THESE ARE ON SRA AS SINGLE RUNS PER SAMPLE, SO DON'T NEED TO DO CATTING
#THIS IS JUST LEFT HERE FOR REFERENCE
#DNA = expect 2 female and 26 male = 28
#RNA = expect 4 female and 11 male = 15
#43 total samples
#2+26+4+11 == 156 total runs
p.r.samples = get_samples(sra, 'Poecilia reticulata')
pr.commands = build_cat_commands(sra, reps=unique(p.r.samples))
length(pr.commands)


#write everything out
all.commands = c(pw.commands,
                 pp.commands,
                 pl.commands,
                 gh.commands,
                 pr.commands)
length(all.commands)
2*length(unique(sra[sra$Organism != 'Poecilia reticulata',]$Sample_Name))
2*length(unique(sra[sra$Organism == 'Poecilia reticulata',]$Sample_Name))
#match up at 182 cat commands
#expect 182 total *_catted.fq files

#write out
fileConn<-file('metadata/cat_guppy_trim_files.txt')
writeLines(all.commands, fileConn)
close(fileConn)



# CREATE ANOTHER SET OF COMMANDS FOR RENAMING INDIVIDUAL RUNS -------------

unique(sra$Organism)
nrow(sra)*2 #matches total number of .trim files I have
species = sra$Organism
species = sub("Gambusia holbrooki", "Gambusia", species)
species = sub("Micropoecilia picta", "picta", species)
species = sub("Poecilia latipinna", "latipinna", species)
species = sub("Poecilia wingei", "wingei", species)
species = sub("Poecilia reticulata", "reticulata", species)
assay = sra$Assay_Type
assay = sub("RNA-Seq", "RNA", assay)
sex = if_else(sra$sex == 'female',
              'F',
              'M')
ind = sra$individual
ret_num = sapply(sra$Sample_Name, function(x) strsplit(x, '_')[[1]][2])
ret_source= sapply(sra$Sample_Name, function(x) strsplit(x, '_')[[1]][3])
ret_ind = sub('Male', '', sra$Sample_Name)
ret_ind = sub('Female', '', ret_ind)
ret_ind = sub('_', '', ret_ind)
ret_ind = sub('_', '', ret_ind)
ret_ind = sub('_', '', ret_ind)

ind = if_else(sra$Organism != "Poecilia reticulata",
              sra$individual,
              ret_ind)
my_tag0 = paste(species, sex, sep='_')
my_tag1 = paste(my_tag0, ind, sep='_')
my_tag = paste(my_tag1, assay, sep='_')
final_tag = paste(my_tag, sra$Run, sep='_')
length(unique(final_tag))


raw_path = '/scratch/02260/grovesd/mank_rebuttal/raw_fastqs/'
r1_files = paste(raw_path, sra$Run, '_1.trim', sep='')
r2_files = paste(raw_path, sra$Run, '_2.trim', sep='')
r1_new_names = paste(final_tag, '_1.trim', sep='')
r2_new_names = paste(final_tag, '_2.trim', sep='')
commandFor = paste('ln -s ', r1_files, ' ./', r1_new_names, sep='')
commandRev = paste('ln -s ', r1_files, ' ./', r2_new_names, sep='')
#write out
fileConn<-file('metadata/rename_single_run_commands.txt')
writeLines(append(commandFor, commandRev), fileConn)
close(fileConn)
#then move to where you want to make the re-named symbolic links, and execute these
#once this is done, check that numbers match expectations
table(sra$Organism)
sra %>% 
  group_by(Organism, Assay_Type) %>% 
  summarize(N=n())

sra %>% 
  group_by(Organism, Assay_Type, sex) %>% 
  summarize(N=n())

