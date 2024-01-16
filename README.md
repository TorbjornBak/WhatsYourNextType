# WhatsYourNextType
WhatsYourNextType is a pipeline for HLA typing Oxford Nanopore reads, developed in a special
course at DTU healthtech. 

The pipelines HLA typing performance is not perfect yet, with correct predictions ranging from 55.56\% for some loci,
to 94.44\% for other loci.

## To Install

First step is to clone the repository into a folder on you pc by using the command:
$ git clone https://github.com/TorbjornBak/WhatsYourNextType

The pipeline uses conda so if you do not have Miniconda installed, please go this conda website: 
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html


In the new terminal run this command to install the blastdatabase in the correct location:
```
bash Setup/blastinstaller.sh
```
or manually:
```
conda env create -f Setup/WYNenvironment.yml
conda activate WYNT
wget https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_gen.fasta -O Data/hla_gen.fasta
makeblastdb -in Data/hla_gen.fasta -out Data/HLAdatabase/HLA_db -dbtype nucl -title HLA_db
```
For macOS, use curl instead of wget:
```
curl https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_gen.fasta -O Data/hla_gen.fasta
```

## Basic usage
Running the pipeline for the first time may take longer than expected as all of the conda environments need to be created. Subsequent runs will be faster. 

For unix, use the `local` or `gridion` profile tag.

For macOS, use the `mac` profile tag. The difference here is that singularity does not work on Mac, so docker is used instead. This means, that docker has to be installed and setup properly. 

```
nextflow run WYNT.nf --fastqfile (LongReadFastqFile) --primerlist (PathToListOfPrimers*)
    \\ -profile (local/gridion/mac) --blastdb (pathToBlastDB)
    \\ (-resume) --coverage (coverage INT)
```

Command line arguments:
```
--fastqfile: The fastqfile should be long read ONT. It is possible to parse multiple files
     to the program by using a glob pattern (i.e. "*").
    To make this work, the path to the fasqfiles need to written in quotation marks ""

--primerlist: The primerlist should be a tab seperated file with two columns. The first is
    the names of the HLA types (i.e. HLAA, HLAB, HLAC, DQB1, DPB1, DRB1). The aminoacid
    sequence of the corresponding primer for each type should be listed in the second
    column (see Data/HLATypesList.txt as an example)

-profile: If run with GridION profile, each process is allocated more processing power
    and more RAM

--blastdb: During the install the of the blast database the folder
     is located at ./Data/HLAdatabase/.
    Should you wish to use another database, then you can specify it here

-resume: Should the pipeline crash, it is possible to resume where it left off,
    when bugs have been fixed. 

--coverage: The max coverage of each type of HLA, meaning after the reads have been split
     into their respective bins. They can be downsampled further to tweak coverage. 
```

## Troubleshooting
Switching between using mamba and conda for handling virtual environments is quite easy. 
As standard, conda is used, but if you want to switch to the faster mamba, 
change the boolean parameter "conda.useMamba" to true in the nextflow.config file.

The assembler (Flye) has a randomizing element to it. Therefore, results may vary slightly and sometimes the assembly may even fail. 
There is an option to make the assembler deterministic, but this heavily increases the run time. 

