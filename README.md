# WhatsYourNextType

## To Install

First step is to clone the repository into a folder on you pc by using the command:
$ git clone https://github.com/TorbjornBak/WhatsYourNextType

The pipeline uses conda so if you do not have Miniconda installed. It can be installed by running these commands
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
sh ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
```
For macOS, use curl instead of wget
Restart the terminal

In the new terminal run this comman to install the blastdatabase in the correct location:
```
sh Scripts/blastinstaller.sh
```
or manually:
```
conda install mamba -c conda-forge

```
Close and reopen the terminal, then continue with:
```
conda env create -f Setup/WYNenvironment.yml
conda activate WYNT
wget http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta -O Data/hla_gen.fasta
makeblastdb -in Data/hla_gen.fasta -out Data/HLAdatabase/HLA_db -dbtype nucl -title HLA_db
```
For macOS, use curl instead of wget:
```
curl http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta -O Data/hla_gen.fasta
```

## To Run 
Running the pipeline for the first time may take longer than expected as all of the conda environments used need to be created. Subsequent runs will be faster. 

```
nextflow run WYNT.nf --fastqfile (LongReadFastqFile) --primerlist (PathToListOfPrimers*) -profile (local/gridion) --blastdb (pathToBlastDB) (-resume) --coverage (coverage INT)

--fastqfile: The fastqfile should be long read ONT. It is possible to parse multiple files to the program by using a glob pattern (i.e. "*"). To make this work, the path to the fasqfiles need to written in quotation marks ""

--primerlist: The primerlist should be a tab seperated file with two columns. The first is the names of the HLA types (i.e. HLAA, HLAB, HLAC, DQB1, DPB1, DRB1). The aminoacid sequence of the corresponding primer for each type should be listed in the second column (see Data/HLATypesList.txt as an example)

-profile: If run with GridION profile, each process is allocated more processing power and more RAM

--blastdb: During the install the of the blast database the folder is located at ./Data/HLAdatabase/
Should you wish to use another database, then you can specify it here

-resume: Should the pipeline crash, it is possible to resume where it left off, if bugs have been fixed. 

--coverage: The max coverage of each type of HLA. So after the reads have been split into their respective bins. They can be downsampled further to tweak coverage. 
```
## Troubleshooting
Sometimes you may be unable to run mamba and should switch to conda instead. Do this by setting the parameter "conda.useMamba" to false in the nextflow.config file

The assembler (Flye) has a randomizing element to it. Therefore, results may vary slightly and sometimes the assembly may even fail. There is an option to make the assembler deterministic, but this heavily increases the run time. 

--
